#include "bp_dec/bp_dec.h"
#include "OSD/OSD.h"
#include "lib_rand/lib_rand.h"
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int DS_MAX_ITER = MAX_ITER;
int DS_DEBUG_OSD = 0;
static int DS_TRACE_SINGLE = 0;

typedef struct {
    int n, mx, mz, k, checks, vars, g_rows;
    const char *name;
} Code;

static Code code;
static uint8_t H[M][N];
static uint8_t G[M_G][N];
static double prior[N][Q];
static double prob_storage[N][Q];
static double *prob[N];
static uint8_t last_run[N];
static uint64_t ds_rng_state=1;
#pragma omp threadprivate(DS_MAX_ITER, DS_DEBUG_OSD)
#pragma omp threadprivate(DS_TRACE_SINGLE, prior, prob_storage, prob, last_run, ds_rng_state)
static double ds_uniform(void) {
    uint64_t z=(ds_rng_state+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    z^=z>>31;
    return (double)(z>>11)*(1.0/9007199254740992.0);
}

static void usage(const char *exe) {
    printf("DS-CSS MBP4+ADOSD\n"
           "Usage: %s --code 18-original|18-ds3|36|36-4-6 [options]\n"
           "  --p-data P      physical depolarizing probability (default 0.01)\n"
           "  --p-meas P      initial probability for each virtual check qubit (default 0.01)\n"
           "  --p-meas-x P    X-syndrome virtual Z-error probability\n"
           "  --p-meas-z P    Z-syndrome virtual X-error probability\n"
           "  --iterations N  maximum MBP iterations (default 100)\n"
           "  --order N       -2 adaptive APOSD/RSR, -1 disabled, 0..12 fixed order (default 4 for 18-ds3, otherwise 2)\n"
           "  --shots N       Monte Carlo shots (default 1000)\n"
           "  --threads N     Monte Carlo worker threads (default 1)\n"
           "  --target-failures N stop after N logical failures (optional; --shots is the maximum)\n"
           "  --always-osd    invoke ADOSD even when BP already matches the syndrome\n"
           "  --seed N        random seed (default 1)\n"
           "  --self-test     run matrix, fault, and logical checks\n", exe);
}

static int parse_int(const char *s, int *out) {
    char *end = NULL; long v;
    errno = 0; v = strtol(s, &end, 10);
    if (errno || !s[0] || *end || v < INT_MIN || v > INT_MAX) return 0;
    *out = (int)v; return 1;
}
static int parse_double(const char *s, double *out) {
    char *end = NULL; errno = 0; *out = strtod(s, &end);
    return !errno && s[0] && !*end && isfinite(*out);
}
static const char *arg_value(int argc, char **argv, int *i, const char *arg) {
    size_t n = strlen(arg);
    if (!strncmp(argv[*i], arg, n) && argv[*i][n] == '=') return argv[*i] + n + 1;
    if (!strcmp(argv[*i], arg) && *i + 1 < argc) return argv[++*i];
    return NULL;
}
static void fail(const char *s) { fprintf(stderr, "Error: %s\n", s); exit(2); }

static FILE *open_code_file(const char *dir, const char *name) {
    char path[1024]; FILE *fp;
    snprintf(path, sizeof(path), "%s\\%s", dir, name);
    fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); exit(2); }
    return fp;
}

static void append_result(const char *which, double p_data, uint64_t failures,
                          int shots, double logical_error_rate) {
    char filename[128];
    FILE *fp;
    /* Approximate 95% binomial interval: two standard errors. */
    double error_bar = 2.0 * sqrt(logical_error_rate * (1.0 - logical_error_rate) / shots);
    long file_size;

    if (snprintf(filename, sizeof(filename), "%s.txt", which) >= (int)sizeof(filename))
        fail("code name is too long for result filename");
    fp = fopen(filename, "a+");
    if (!fp) {
        fprintf(stderr, "Cannot append results to %s: %s\n", filename, strerror(errno));
        exit(2);
    }
    if (fseek(fp, 0, SEEK_END) != 0 || (file_size = ftell(fp)) < 0) {
        fprintf(stderr, "Cannot inspect result file %s\n", filename);
        fclose(fp);
        exit(2);
    }
    if (file_size == 0) fprintf(fp, "//p_data\tlogical_error_rate\terror_bar\n");
    fprintf(fp, "%.12g\t%.12g\t%.12g\n", p_data, logical_error_rate, error_bar);
    if (fclose(fp) != 0) {
        fprintf(stderr, "Cannot finish writing result file %s\n", filename);
        exit(2);
    }
    printf("results_appended_to=%s\n", filename);
}

/* Parse the upstream GF(4) alist and expand it into bounded runtime storage. */
static void read_alist(FILE *fp, int expected_n, int expected_m, uint8_t *dst,
                       int stride, int cap_n, int cap_m) {
    int n, m, max_col, max_row, i, j, col, value;
    int *cols, *rows;
    if (fscanf(fp, "%d%d%d%d", &n, &m, &max_col, &max_row) != 4 ||
        n != expected_n || m != expected_m || n > cap_n || m > cap_m)
        fail("matrix dimensions or alist header are invalid");
    cols = calloc((size_t)n, sizeof(*cols)); rows = calloc((size_t)m, sizeof(*rows));
    if (!cols || !rows) fail("out of memory reading matrix");
    for (i=0;i<n;i++) if (fscanf(fp,"%d",&cols[i])!=1 || cols[i]<0 || cols[i]>max_col) fail("bad column degrees");
    for (i=0;i<m;i++) if (fscanf(fp,"%d",&rows[i])!=1 || rows[i]<0 || rows[i]>max_row) fail("bad row degrees");
    /* Skip column adjacency: row adjacency is the authoritative matrix. */
    for (i=0;i<n;i++) for (j=0;j<max_col;j++) if (fscanf(fp,"%d%d",&col,&value)!=2) fail("truncated column adjacency");
    memset(dst, 0, (size_t)cap_m * (size_t)stride);
    for (i=0;i<m;i++) {
        for (j=0;j<max_row;j++) {
            if (fscanf(fp,"%d%d",&col,&value)!=2) fail("truncated row adjacency");
            if (j<rows[i]) {
                if (col<1 || col>n || value<1 || value>3 || dst[(size_t)i*stride+col-1]) fail("invalid or duplicate GF(4) matrix entry");
                dst[(size_t)i*stride+col-1]=(uint8_t)value;
            } else if (col || value) fail("bad alist padding");
        }
    }
    free(cols); free(rows); fclose(fp);
}

static void configure_code(const char *which) {
    const char *dir;
    if (!strcmp(which,"18-original")) { code=(Code){18,6,6,6,12,30,24,"[[18,6,3]] original"}; dir="..\\codes\\18-original"; }
    else if (!strcmp(which,"18-ds3")) { code=(Code){18,6,6,6,12,30,24,"[[18,6,3]] DS-distance-3"}; dir="..\\codes\\18-ds3"; }
    else if (!strcmp(which,"36")) { code=(Code){36,16,16,12,32,68,48,"[[36,12,5]]"}; dir="..\\codes\\36"; }
    else if (!strcmp(which,"36-4-6")) { code=(Code){36,23,23,4,46,82,40,"[[36,4,6]] sparse redundant DS-distance-5"}; dir="..\\codes\\36-4-6"; }
    else fail("--code must be 18-original, 18-ds3, 36, or 36-4-6");
    read_alist(open_code_file(dir,"DS_A4.txt"),code.vars,code.checks,&H[0][0],N,N,M);
    read_alist(open_code_file(dir,"DS_G4s.txt"),code.vars,code.g_rows,&G[0][0],N,N,M_G);
}

static void build_sparse(a_matrix_GFQ *A, g_matrix_GFQ *normalizer) {
    int n,m;
    A->NN=code.vars; A->MM=code.checks; A->max_num_n=0; A->max_num_m=0;
    for(n=0;n<code.vars;n++) {
        int count=0; for(m=0;m<code.checks;m++) if(H[m][n]) count++;
        A->num_n[n]=count; if(count>A->max_num_n) A->max_num_n=count;
        A->nlist[n]=calloc((size_t)(count?count:1),sizeof(int32_t));
        A->nlist_val[n]=calloc((size_t)(count?count:1),sizeof(GFQ_t));
        A->ni2j[n]=calloc((size_t)(count?count:1),sizeof(int32_t));
    }
    for(m=0;m<code.checks;m++) {
        int count=0; for(n=0;n<code.vars;n++) if(H[m][n]) count++;
        A->num_m[m]=count; if(count>A->max_num_m) A->max_num_m=count;
        A->mlist[m]=calloc((size_t)(count?count:1),sizeof(int32_t));
        A->mlist_val[m]=calloc((size_t)(count?count:1),sizeof(GFQ_t));
        A->mj2i[m]=calloc((size_t)(count?count:1),sizeof(int32_t));
    }
    {
        int ci[N]={0}, ri[M]={0};
        for(m=0;m<code.checks;m++) for(n=0;n<code.vars;n++) if(H[m][n]) {
            int a=ci[n]++, b=ri[m]++;
            A->nlist[n][a]=m; A->nlist_val[n][a]=H[m][n]; A->ni2j[n][a]=b;
            A->mlist[m][b]=n; A->mlist_val[m][b]=H[m][n]; A->mj2i[m][b]=a;
        }
    }
    normalizer->NN=code.vars; normalizer->MM=M_G; normalizer->max_num_n=0; normalizer->max_num_m=0;
    for(n=0;n<code.vars;n++) {
        int count=0; for(m=0;m<code.g_rows;m++) if(G[m][n]) count++;
        normalizer->num_n[n]=count; if(count>normalizer->max_num_n) normalizer->max_num_n=count;
        normalizer->nlist[n]=calloc((size_t)(count?count:1),sizeof(int32_t));
        normalizer->nlist_val[n]=calloc((size_t)(count?count:1),sizeof(GFQ_t));
        normalizer->ni2j[n]=calloc((size_t)(count?count:1),sizeof(int32_t));
    }
    for(m=0;m<M_G;m++) {
        int count=0; for(n=0;n<code.vars;n++) if(m<code.g_rows && G[m][n]) count++;
        normalizer->num_m[m]=count; if(count>normalizer->max_num_m) normalizer->max_num_m=count;
        normalizer->mlist[m]=calloc((size_t)(count?count:1),sizeof(int32_t));
        normalizer->mlist_val[m]=calloc((size_t)(count?count:1),sizeof(GFQ_t));
    }
    {
        int ci[N]={0}, ri[M_G]={0};
        for(m=0;m<code.g_rows;m++) for(n=0;n<code.vars;n++) if(G[m][n]) {
            int a=ci[n]++, b=ri[m]++;
            normalizer->nlist[n][a]=m; normalizer->nlist_val[n][a]=G[m][n];
            normalizer->mlist[m][b]=n; normalizer->mlist_val[m][b]=G[m][n];
        }
    }
}

static void setup_osd(OSD *osd) {
    int m,n; RV_t mask;
    initOSD(osd); osd->osdw=-2; osd->RankH=(uint32_t)code.checks;
    memset(osd->parity_intact[0],0,0); /* allocations are row-separated */
    for(m=0;m<M;m++) memset(osd->parity_intact[m],0,(size_t)osd->P_col*sizeof(RV_t));
    for(m=0;m<code.checks;m++) for(n=0;n<code.vars;n++) if(H[m][n]) {
        unsigned pos=(unsigned)n>>QShift;
        unsigned shift=(QMask-((unsigned)n&QMask))<<1;
        mask=(H[m][n]==1)?2:1;
        osd->parity_intact[m][pos] |= mask<<shift;
    }
}

static void set_priors(QBP_Ctl *bp, a_matrix_GFQ *A, double p_data, double p_meas_x, double p_meas_z) {
    int n,m,j,q;
    for(n=0;n<N;n++) {
        for(q=0;q<Q;q++) prior[n][q]=0.0;
        if(n<code.n) { prior[n][0]=1-p_data; prior[n][1]=prior[n][2]=prior[n][3]=p_data/3.0; }
        else if(n<code.vars) {
            int xsector=(n<code.n+code.mx);
            double p=xsector?p_meas_x:p_meas_z;
            prior[n][0]=1-p; prior[n][xsector?2:1]=p;
        } else prior[n][0]=1.0;
        memcpy(bp->pn[n],prior[n],sizeof(prior[n]));
    }
    for(m=0;m<A->MM;m++) for(j=0;j<A->num_m[m];j++) {
        n=A->mlist[m][j]; q=A->mlist_val[m][j];
        bp->dqml[m][j]=2.0*(prior[n][0]+prior[n][q])-1.0;
    }
}

static void parity_syndrome(const uint8_t *e, uint8_t *s) {
    int m,n; memset(s,0,(size_t)code.checks);
    for(m=0;m<code.checks;m++) for(n=0;n<code.vars;n++) if(H[m][n] && e[n])
        s[m] ^= (uint8_t)(((e[n]>>1)&H[m][n]) ^ (e[n]&(H[m][n]>>1)));
}

static int correction_matches_syndrome(const uint8_t *e, const uint8_t *syndrome) {
    uint8_t check[M];
    parity_syndrome(e,check);
    return memcmp(check,syndrome,(size_t)code.checks)==0;
}

static int valid_virtual(const uint8_t *e) {
    int n; for(n=code.n;n<code.vars;n++) {
        uint8_t v=e[n];
        if(n<code.n+code.mx) { if(v!=0 && v!=2) return 0; }
        else if(v!=0 && v!=1) return 0;
    }
    for(n=code.vars;n<N;n++) if(e[n]) return 0;
    return 1;
}

static void repair_virtual(uint8_t *e, const uint8_t *syndrome) {
    int m,n; uint8_t *pred=calloc((size_t)code.checks,1); if(!pred) fail("out of memory repairing virtual components");
    for(m=0;m<code.checks;m++) {
        for(n=0;n<code.n;n++) if(H[m][n] && e[n]) pred[m]^=(uint8_t)(((e[n]>>1)&H[m][n]) ^ (e[n]&(H[m][n]>>1)));
        if(m<code.mx) e[code.n+m]=((pred[m]^syndrome[m])?2:0);
        else e[code.n+code.mx+(m-code.mx)]=((pred[m]^syndrome[m])?1:0);
    }
    free(pred);
}

static int logical_failure(const uint8_t *truth, const uint8_t *guess) {
    int row,n; uint8_t syndrome;
    for(row=0;row<code.g_rows;row++) {
        syndrome=0;
        for(n=0;n<code.n;n++) if(G[row][n]) {
            uint8_t d=truth[n]^guess[n];
            syndrome^=(uint8_t)(((d>>1)&G[row][n]) ^ (d&(G[row][n]>>1)));
        }
        if(syndrome) return 1;
    }
    return 0;
}

static uint8_t *post_decode_ds(OSD *osd, const uint8_t *syndrome,
                               const double **probabilities, const uint8_t *last_for,
                               int order) {
    if(order==-2) return post_decodeOSDfull(osd,syndrome,probabilities,last_for,(uint32_t)DS_MAX_ITER);
    return post_decodeOSDw(osd,syndrome,probabilities,last_for);
}

static int decode_fault(QBP_Ctl *bp, a_matrix_GFQ *A, OSD *osd, uint8_t fault_pos, uint8_t pauli,
                        double p_data, double p_meas_x, double p_meas_z, int always_osd) {
    uint8_t truth[N]={0}, syndrome[M]={0}, decoded[N]={0};
    unsigned iter=0; int ok, bp_consistent, osd_called=0, trace=(DS_TRACE_SINGLE && fault_pos==4 && pauli==3);
    truth[fault_pos]=pauli; parity_syndrome(truth,syndrome);
    if(trace) {
        printf("TRACE fault: data[%u]=Y, no measurement faults; syndrome=",fault_pos);
        for(int m=0;m<code.checks;m++) printf("%u",syndrome[m]);
        printf("\n");
    }
    ok=Qbp_init20(bp,syndrome,A,p_data,0.0); set_priors(bp,A,p_data,p_meas_x,p_meas_z);
    if(!ok) {
        while(iter<(unsigned)DS_MAX_ITER) { iter++; if(Qbp_dec20(bp,A)) {ok=1;break;} }
    }
    bp_consistent=correction_matches_syndrome(bp->tt,syndrome);
    if(trace) {
        printf("TRACE BP hard decision:");
        for(int n=0;n<code.vars;n++) printf(" %u",bp->tt[n]);
        printf("\n");
        {
            int n=4;
            double sum=bp->qn[n][0]+bp->qn[n][1]+bp->qn[n][2]+bp->qn[n][3];
            printf("TRACE posterior data[%d]: I=%.17g X=%.17g Z=%.17g Y=%.17g\n",n,bp->qn[n][0]/sum,bp->qn[n][1]/sum,bp->qn[n][2]/sum,bp->qn[n][3]/sum);
        }
        for(int n=code.n;n<code.vars;n++) {
            double sum=bp->qn[n][0]+bp->qn[n][1]+bp->qn[n][2]+bp->qn[n][3];
            printf("TRACE posterior virtual[%d]: I=%.17g X=%.17g Z=%.17g Y=%.17g\n",n-code.n,bp->qn[n][0]/sum,bp->qn[n][1]/sum,bp->qn[n][2]/sum,bp->qn[n][3]/sum);
        }
    }
    if((bp_consistent && !always_osd) || osd->osdw==-1) {
        memcpy(decoded,bp->tt,sizeof(decoded));
    } else {
        for(int n=0;n<N;n++) {
            double sum=bp->qn[n][0]+bp->qn[n][1]+bp->qn[n][2]+bp->qn[n][3];
            if(sum<=0 || !isfinite(sum)) { prob[n][0]=1; prob[n][1]=prob[n][2]=prob[n][3]=0; }
            else for(int q=0;q<Q;q++) prob[n][q]=bp->qn[n][q]/sum;
            last_run[n]=1;
        }
        int old_trace=DS_DEBUG_OSD; DS_DEBUG_OSD=trace;
        osd_called=1;
        memcpy(decoded,post_decode_ds(osd,syndrome,(const double **)prob,last_run,osd->osdw),sizeof(decoded));
        DS_DEBUG_OSD=old_trace;
    }
    if(trace) {
        printf("TRACE control result: bp_converged=%s bp_iterations=%u osd_called=%s decoder_source=%s\n",
               ok?"true":"false",iter,osd_called?"true":"false",osd_called?"ADOSD":"BP");
        printf("TRACE decoder correction before virtual repair:");
        for(int n=0;n<code.vars;n++) printf(" %u",decoded[n]);
        printf("\n");
    }
    repair_virtual(decoded,syndrome);
    if(trace) { printf("TRACE final correction after virtual repair:"); for(int n=0;n<code.vars;n++) printf(" %u",decoded[n]); printf("\n"); }
    if(!valid_virtual(decoded)) return 0;
    if(trace) printf("TRACE logical_failure=%d\n",logical_failure(truth,decoded));
    if(logical_failure(truth,decoded)) return 0;
    return 1;
}

static int matrix_tests(a_matrix_GFQ *A) {
    int m,n,failures=0; uint8_t e[N]={0}, s[M], reference[M];
    /* Every allowed single fault must have a reproducible augmented syndrome. */
    for(n=0;n<code.n;n++) for(uint8_t p=1;p<=3;p++) {
        memset(e,0,sizeof(e)); e[n]=p; parity_syndrome(e,s);
        Quan_GenSyndrome(A,e,reference);
        if(memcmp(s,reference,(size_t)code.checks)) failures++;
    }
    for(m=0;m<code.checks;m++) {
        memset(e,0,sizeof(e)); e[code.n+m]=(m<code.mx)?2:1; parity_syndrome(e,s);
        if(s[m]!=1 || !valid_virtual(e)) failures++;
        for(n=0;n<code.checks;n++) if(n!=m && s[n]) failures++;
    }
    return failures;
}

static int self_test(a_matrix_GFQ *A, g_matrix_GFQ *g, OSD *osd, QBP_Ctl *bp, const char *which, int order,
                     double p_data, double p_meas_x, double p_meas_z, int always_osd) {
    int i, errors=matrix_tests(A);
    uint8_t zero[N]={0}, syndrome[M]={0}, diff[N]={0};
    if(errors) { fprintf(stderr,"matrix/single-fault consistency: FAIL (%d)\n",errors); return 1; }
    printf("matrix and allowed single-fault syndrome consistency: PASS\n");
    if(DS_TRACE_SINGLE) {
        int correct;
        osd->osdw=(int8_t)order;
        correct=decode_fault(bp,A,osd,4,3,p_data,p_meas_x,p_meas_z,always_osd);
        printf("TRACE focused q4-Y result=%s\n",correct?"PASS":"FAIL");
        return correct?0:1;
    }
    if(strcmp(which,"18-original")!=0) {
        int correct=1; osd->osdw=(int8_t)order;
        for(i=0;i<code.n;i++) for(int p=1;p<=3;p++) if(!decode_fault(bp,A,osd,(uint8_t)i,(uint8_t)p,p_data,p_meas_x,p_meas_z,always_osd)) { correct=0; fprintf(stderr,"single-fault decode failed at physical qubit %d Pauli %d\n",i,p); }
        for(i=code.n;i<code.vars;i++) if(!decode_fault(bp,A,osd,(uint8_t)i,(uint8_t)(i<code.n+code.mx?2:1),p_data,p_meas_x,p_meas_z,always_osd)) { correct=0; fprintf(stderr,"single-fault decode failed at virtual qubit %d\n",i-code.n); }
        osd->osdw=order;
        if(!correct) { fprintf(stderr,"single-fault end-to-end decode: FAIL\n"); return 1; }
        printf("single-fault end-to-end physical and measurement decoding: PASS\n");
    }
    if(!is_zero_vec(zero,N)) errors++;
    /* A physical operator in the G span commutes with all rows, while
       a one-qubit Pauli outside the code normalizer must be detected by G. */
    for(i=0;i<code.n && errors==0;i++) for(int p=1;p<=3;p++) {
        memset(diff,0,sizeof(diff)); diff[i]=(uint8_t)p;
        if(!logical_failure(diff,zero)) errors++;
    }
    memset(diff,0,sizeof(diff)); diff[code.n]=(code.mx?2:1);
    if(logical_failure(diff,zero)) errors++;
    if(errors) { fprintf(stderr,"physical logical-test setup: FAIL\n"); return 1; }
    printf("physical normalizer logical-success check: PASS (measurement-only differences ignored)\n");
    Qbp_init20(bp,syndrome,A,0.01,0.0);
    set_priors(bp,A,0.01,0.01,0.01);
    if(bp->pn[code.n][1]!=0.0 || bp->pn[code.n][3]!=0.0 || bp->pn[code.n][0]!=1.0-0.01 || bp->pn[code.n][2]!=0.01 ||
       bp->pn[code.n+code.mx][2]!=0.0 || bp->pn[code.n+code.mx][3]!=0.0 || bp->pn[code.n+code.mx][1]!=0.01) {
        fprintf(stderr,"virtual prior X-sector: %.17g %.17g %.17g %.17g; Z-sector: %.17g %.17g %.17g %.17g\n",
          bp->pn[code.n][0],bp->pn[code.n][1],bp->pn[code.n][2],bp->pn[code.n][3],
          bp->pn[code.n+code.mx][0],bp->pn[code.n+code.mx][1],bp->pn[code.n+code.mx][2],bp->pn[code.n+code.mx][3]);
        fail("zero-probability virtual Pauli components are not preserved");
    }
    printf("zero-probability virtual Pauli components: PASS\n");
    (void)g; (void)osd;
    return 0;
}

int main(int argc, char **argv) {
    const char *which=NULL; double p_data=0.01,p_meas=0.01,p_meas_x=0.01,p_meas_z=0.01;
    int have_meas_x=0,have_meas_z=0,have_order=0,iterations=100,order=2,shots=1000,target_failures=0,seed=1,self=0,always_osd=0,threads=1,i;
    int shots_run=0;
    uint64_t next_shot=0;
    uint64_t failures=0; double total_iterations=0.0;
    setvbuf(stdout,NULL,_IONBF,0);
    DS_TRACE_SINGLE = getenv("DS_TRACE_SINGLE") != NULL;
    for(i=1;i<argc;i++) {
        const char *v;
        if(!strcmp(argv[i],"--help")||!strcmp(argv[i],"-h")){usage(argv[0]);return 0;}
        if(!strcmp(argv[i],"--self-test")){self=1;continue;}
        if(!strcmp(argv[i],"--always-osd")){always_osd=1;continue;}
        if((v=arg_value(argc,argv,&i,"--code"))) {which=v;continue;}
        if((v=arg_value(argc,argv,&i,"--p-data"))) {if(!parse_double(v,&p_data))fail("invalid --p-data");continue;}
        if((v=arg_value(argc,argv,&i,"--p-meas"))) {if(!parse_double(v,&p_meas))fail("invalid --p-meas");continue;}
        if((v=arg_value(argc,argv,&i,"--p-meas-x"))) {if(!parse_double(v,&p_meas_x))fail("invalid --p-meas-x");have_meas_x=1;continue;}
        if((v=arg_value(argc,argv,&i,"--p-meas-z"))) {if(!parse_double(v,&p_meas_z))fail("invalid --p-meas-z");have_meas_z=1;continue;}
        if((v=arg_value(argc,argv,&i,"--iterations"))) {if(!parse_int(v,&iterations))fail("invalid --iterations");continue;}
        if((v=arg_value(argc,argv,&i,"--order"))) {if(!parse_int(v,&order))fail("invalid --order");have_order=1;continue;}
        if((v=arg_value(argc,argv,&i,"--shots"))) {if(!parse_int(v,&shots))fail("invalid --shots");continue;}
        if((v=arg_value(argc,argv,&i,"--threads"))) {if(!parse_int(v,&threads))fail("invalid --threads");continue;}
        if((v=arg_value(argc,argv,&i,"--target-failures"))) {if(!parse_int(v,&target_failures))fail("invalid --target-failures");continue;}
        if((v=arg_value(argc,argv,&i,"--seed"))) {if(!parse_int(v,&seed))fail("invalid --seed");continue;}
        fprintf(stderr,"Unknown or incomplete option: %s\n",argv[i]); usage(argv[0]); return 2;
    }
    if(!which) fail("--code is required");
    if(!have_meas_x) p_meas_x=p_meas;
    if(!have_meas_z) p_meas_z=p_meas;
    if(p_data<0 || p_data>=0.75 || p_meas<0 || p_meas>=0.5 || p_meas_x<0 || p_meas_x>=0.5 || p_meas_z<0 || p_meas_z>=0.5 || iterations<1 || iterations>MAX_ITER || order< -2 || order>12 || shots<1 || threads<1 || threads>256 || target_failures<0 || (target_failures>0 && p_data==0 && p_meas_x==0 && p_meas_z==0)) fail("parameter out of range or target failures cannot be reached");
    if(always_osd && order==-1) fail("--always-osd cannot be combined with --order -1");
    if(threads>shots) threads=shots;
    omp_set_dynamic(0);
    configure_code(which);
    if(!have_order && !strcmp(which,"18-ds3")) order=4;
    a_matrix_GFQ Aobj; g_matrix_GFQ Gobj;
    a_matrix_GFQ *A=&Aobj; g_matrix_GFQ *g=&Gobj;
    memset(A,0,sizeof(*A)); memset(g,0,sizeof(*g));
    build_sparse(A,g);
    printf("Code: %s; data=%d virtual=%d checks=%d; p_data=%.8g p_meas_x=%.8g p_meas_z=%.8g; MBP iterations=%d ADOSD order=%d always_osd=%s shots=%d seed=%d threads=%d\n",
        code.name,code.n,code.mx+code.mz,code.checks,p_data,p_meas_x,p_meas_z,iterations,order,always_osd?"true":"false",shots,seed,threads);
    if(self) {
        QBP_Ctl bpobj; OSD osdobj;
        QBP_Ctl *bp=&bpobj; OSD *osd=&osdobj;
        memset(bp,0,sizeof(*bp));
        alloc_QBPC(A,bp); setup_osd(osd);
        OSD_set_runtime_shape(osd,(uint32_t)code.vars,(uint32_t)code.checks,(uint32_t)code.checks);
        osd->osdw=(int8_t)order;
        for(i=0;i<N;i++) prob[i]=prob_storage[i];
        DS_MAX_ITER=iterations; DS_DEBUG_OSD=0; ds_rng_state=(uint64_t)(uint32_t)seed; seed64((uint64_t)(uint32_t)seed); rnd256_init();
        return self_test(A,g,osd,bp,which,order,p_data,p_meas_x,p_meas_z,always_osd);
    }
#pragma omp parallel num_threads(threads) shared(A, p_data, p_meas_x, p_meas_z, iterations, order, shots, target_failures, always_osd, seed, failures, total_iterations, shots_run, next_shot)
    {
        QBP_Ctl bpobj; OSD osdobj; QBP_Ctl *worker_bp=&bpobj; OSD *worker_osd=&osdobj;
        int worker_id=omp_get_thread_num();
        uint64_t worker_seed=(uint64_t)(uint32_t)seed + UINT64_C(0x9e3779b97f4a7c15)*(uint64_t)(worker_id+1);
        memset(worker_bp,0,sizeof(*worker_bp));
        alloc_QBPC(A,worker_bp); setup_osd(worker_osd);
        OSD_set_runtime_shape(worker_osd,(uint32_t)code.vars,(uint32_t)code.checks,(uint32_t)code.checks);
        worker_osd->osdw=(int8_t)order;
        for(int pi=0;pi<N;pi++) prob[pi]=prob_storage[pi];
        DS_MAX_ITER=iterations; DS_DEBUG_OSD=0; DS_TRACE_SINGLE=0;
        ds_rng_state=worker_seed; seed64(worker_seed); rnd256_init();
        for(;;) {
        uint64_t shot_index, failure_snapshot;
        #pragma omp atomic capture
        { shot_index=next_shot; next_shot++; }
        if(shot_index>=(uint64_t)shots) break;
        #pragma omp atomic read
        failure_snapshot=failures;
        if(target_failures>0 && failure_snapshot>=(uint64_t)target_failures) break;
        uint8_t truth[N]={0}, syndrome[M]={0}, decoded[N]={0}, diff[N]={0}, logical[M_G]={0};
        int n, bp_consistent; unsigned iter=0; int ok;
        for(n=0;n<code.n;n++) {
            double r=ds_uniform();
            if(r<p_data) truth[n]=(uint8_t)(1+(int)(ds_uniform()*3.0));
        }
        for(n=code.n;n<code.vars;n++) if(ds_uniform()<(n<code.n+code.mx?p_meas_x:p_meas_z)) truth[n]=(n<code.n+code.mx)?2:1;
        parity_syndrome(truth,syndrome);
        ok=Qbp_init20(worker_bp,syndrome,A,p_data,0.0);
        set_priors(worker_bp,A,p_data,p_meas_x,p_meas_z);
        if(!ok) {
            while(iter<(unsigned)iterations) { iter++; if(Qbp_dec20(worker_bp,A)) {ok=1;break;} }
        }
        bp_consistent=correction_matches_syndrome(worker_bp->tt,syndrome);
        if((bp_consistent && !always_osd) || order==-1) memcpy(decoded,worker_bp->tt,sizeof(decoded));
        else {
            for(n=0;n<N;n++) {
                double sum=0; for(int q=0;q<Q;q++) sum+=worker_bp->qn[n][q];
                if(sum<=0 || !isfinite(sum)) { prob[n][0]=1; prob[n][1]=prob[n][2]=prob[n][3]=0; }
                else for(int q=0;q<Q;q++) prob[n][q]=worker_bp->qn[n][q]/sum;
                last_run[n]=1;
            }
            memcpy(decoded,post_decode_ds(worker_osd,syndrome,(const double **)prob,last_run,order),sizeof(decoded));
        }
        repair_virtual(decoded,syndrome);
        if(!valid_virtual(decoded)) fail("decoder produced a forbidden virtual Pauli component");
        uint8_t check[M]; parity_syndrome(decoded,check);
        if(memcmp(check,syndrome,(size_t)code.checks)) fail("decoded correction is inconsistent with observed syndrome");
        for(n=0;n<code.vars;n++) diff[n]=truth[n]^decoded[n];
        Quan_DegSyndrome(g,diff,logical);
        if(!is_zero_vec(logical,M_G)) {
            #pragma omp atomic update
            failures++;
        }
        #pragma omp atomic update
        total_iterations+=iter;
        #pragma omp atomic update
        shots_run++;
        }
        free_OSD(worker_osd);
    }
    double logical_error_rate=(shots_run>0)?(double)failures/shots_run:0.0;
    printf("logical_error_rate=%.10g failures=%llu shots=%d average_iterations=%.6f",
        logical_error_rate,(unsigned long long)failures,shots_run,shots_run?total_iterations/shots_run:0.0);
    if(target_failures>0) printf(" target_failures=%d reached=%s",target_failures,failures>=(uint64_t)target_failures?"true":"false");
    printf("\n");
    append_result(which,p_data,failures,shots_run,logical_error_rate);
    return 0;
}
