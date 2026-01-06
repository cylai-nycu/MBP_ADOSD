#include "../ldpc_parm.h"

/* (By convention (systematic form), the righthand M*M matrix is an invertible matrix)
a_matrix *A to
[1 1 0 1;
 0 0 1 1;
 1 1 1 1]
will have N=4; M=3,
  *num_n = [2 2 2 3], max_num_n = 3,
  *num_m = [3 2 4], max_num_m = 4,
  *nlist = [c1 c2 c3 c4] with *c1=[1 3 0], *c2=[1 3 0], *c3=[2 3 0], *c4=[1 2 3], all non-zero -1 as index;
  *mlist = [r1 r2 r3] with *r1=[1 2 4 0], *r2=[3 4 0 0], *r3=[1 2 3 4], all non-zero -1 as index;
  *ni2j = [j1 j2 j3 j4] with *j1=[0 0], *j2=[1 1], *j3=[0 2], *j4=[2 1 3];
          j1 means in col_1 the two ones are both the first ones in row order
*/
/*typedef struct {
  int32_t NN , MM ;
  int32_t max_num_n ;
  int32_t max_num_m ;
  int32_t num_n[N];
  int32_t num_m[M];
  int32_t *nlist[N];
  int32_t *mlist[M];
  int32_t *ni2j[N];
} a_matrix ;*/
typedef struct {    // note: use 2N to support stabilizer check matrix in GF2 form
  int32_t NN , MM ;
  int32_t max_num_n ;
  int32_t max_num_m ;
  int32_t num_n[2*N];
  int32_t num_m[M];
  int32_t *nlist[2*N];
  int32_t *mlist[M];
  int32_t *ni2j[2*N];
} a2_matrix ;


/* (By convention (systematic form), the righthand M*M matrix is an invertible matrix FOR LINEAR CODE)
a_matrix_GFQ *A to
[1 2 2 1 0;
 0 1 2 2 1;
 1 0 1 2 2;
 2 1 0 1 2]
will have N=5; M=4,
  *num_n = [3 3 3 4 3], max_num_n = 4,
  *num_m = [4 4 4 4], max_num_m = 4,
  *nlist = [c1 c2 c3 c4 c5]     with *c1=[1 3 4 0], *c2=[1 2 4 0], *c3=[1 2 3 0], *c4=[1 2 3 4], *c5=[2 3 4 0], all non-zero -1 as index;
  *nlist_val = [v1 v2 v3 v4 v5] with *v1=[1 1 2 0], *v2=[2 1 1 0], *v3=[2 2 1 0], *v4=[1 2 2 1], *v5=[2 3 4 0];
  *mlist = [r1 r2 r3 r4]        with *r1=[1 2 3 4], *r2=[2 3 4 5], *r3=[1 3 4 5], *r4=[1 2 4 5], all non-zero -1 as index;;
  *mlist_val = [u1 u2 u3 u4 u5] with *u1=[1 2 2 1], *u2=[1 2 2 1], *u3=[1 1 2 2], *u4=[2 1 1 2];
  *ni2j = [j1 j2 j3 j4] with *j1=[0 0 0], *j2=[1 0 1], *j3=[2 1 1], *j4=[3 2 2 2], *j5=[3 3 3];
          j1 means in col_1 the three non-zeros are all the first non-zeros in row order
*/
typedef struct {
  int32_t NN , MM ;
  int32_t max_num_n ;
  int32_t max_num_m ;
  int32_t num_n[N];
  int32_t num_m[M];
  int32_t *nlist[N];  GFQ_t *nlist_val[N];
  int32_t *mlist[M];  GFQ_t *mlist_val[M];
  int32_t *ni2j[N];
  int32_t *mj2i[M];
} a_matrix_GFQ ;

typedef struct {
  int32_t NN , MM ;
  int32_t max_num_n ;
  int32_t max_num_m ;
  int32_t num_n[N];
  int32_t num_m[M_G];
  int32_t *nlist[N];    GFQ_t *nlist_val[N];
  int32_t *mlist[M_G];  GFQ_t *mlist_val[M_G];
  int32_t *ni2j[N];
} g_matrix_GFQ ;


typedef struct {  // used with g_matrix_GFQ
  uint8_t reset;    // set to 1 when init
  //double  p_ch;     // channel error probability (BSC generalized for "Quantum")
  GFQ_t target_z[M];    // target syndrome vector to converge

  double pn[N][Q];  // For each of the N variable nodes, probability vector "RELATED" to {I,X,Z,Y}; (vertical step init)
  //-- COMMENTS IGNORE idx [q], which implicitly runs from q=0 to q<[Q-1]=3 for {X,Z,Y}
  // HORIZONTAL:
  double dm[M];     // dm[m] = rm[m]^0 - rm[m]^1          // (MacKay '99 IT paper's delta_r)
  double *dml[M];       // to alloc as d[m][l] with 0 <= l < num_m[m]  (like MacKay's (49))
  double (*rml[M])[2];  // to alloc as rml[m][l] with 0 <= l < num_m[m], each rml[m][l] buffers r_ml^0 and r_ml^1
  //double *rml[M];       // to alloc as rml[m][l] with 0 <= l < num_m[m], each rml[m][l] buffers r_ml^0 (and so r_ml^1 = 1 - r_ml^0)
  double *dqml[M];      // dqml[m][l] = q_ml^0 - q_ml^1 (see MacKay's (47)-(51))
  // VERTICAL:      // use qlm to denote q_ml (see MacKay's (47)-(51))
  //double *qlm[N];   // to alloc as qlm[l][m] with 0 <= m < num_n[l], each qlm[m][l] buffers q_ml^0 (and q_ml^1 = 1 - q_ml^0)
  double qn[N][Q];  // updated probability vectors of variable nodes (vertical step)

#if EN_AFP_V
  double *dqml_ori[M];    // alternative dqml without applying alpha_v, to support alpha_v prime
#endif

  GFQ_t tt[N];    // updated decoded vector (estimated TX vector after BP)
  GFQ_t zz[M];    // updated syndrome vector

#if DBL_CHK
  GFQ_t tt_last[N];
#endif

  int32_t iter;         // iter cnt
  double qn_avg[N][Q];  // for POS_AVG

#if STOP_BY_C   // These variables do not need to be delivered between iterations,
  double dm_c[M];     // put them here due to:
  double *dqml_c[M];  // can create the data structure of dqml_c with dqml
  GFQ_t zz_c[M];
#endif
  uint8_t is_hit_c;     // indicator for if hit the required condition for STOP_BY_C (CheckNode syndrome OK)


  double ap;      // to pass the parameter alpha prime (AFP)

  int32_t seq[N]; // should contain 0 to N-1 in some order to run serial decoder (dec24)

#if MAJ_SCHE
  int32_t maj_cnt[N];  // add one vote if a stabilizer associated to this qubit has syndrome=1
  int32_t maj_max;     // maximum of maj_cnt
#endif

#if CYC_CHK
  int32_t *cyc_cnt[N];  // cycle count [n][m] (0: no cycles (1 overlap), 1: has one 4-cycle (2 overlap), k: has k+1 overlaps with other rows)

  int32_t **cycs[N];    // [N][sM(n)][m] so that we can for m = 1 to sM(n)
  double (*rml2[M])[2]; // to buffer original (or suppressed) rml
#endif

} QBP_Ctl ;  // "Quantum" Belief Propagation (Sum-Product, Message-Passing) Control block





/*//-- Matrix vector multiply Ax = y
void MatVecMultiply(a_matrix *A, uint8_t *x, uint8_t *y);
//-- generate syndrome Hr = z
#define GenSyndrome(H, r, z)  MatVecMultiply(H, r, z);*/
//-- generate syndrome for "Quantum"
void Quan_GenSyndrome(a_matrix_GFQ *A, GFQ_t *x, GFQ_t *y);   // input check matrix, and will check M = N-K rows
void Quan_DegSyndrome(g_matrix_GFQ *G, GFQ_t *x, GFQ_t *y);   // input generator matrix, and will check M+2K = N+K rows
//void Quan_LogSyndrome(g_matrix_GFQ *G, GFQ_t *x, GFQ_t *y);   // input generator matrix, and will check ONLY 2K rows


//-- return 1 if v is a zero vector; otherwise return 0
uint8_t is_zero_vec(GFQ_t *v, int32_t len);
//-- check Hamming weight of vector v
uint32_t HamWt(GFQ_t *v, int32_t len);
//-- check Hamming distance between x and y
uint32_t HamDist(GFQ_t *x, GFQ_t *y, int32_t len);
//-- get v = y - x over GF4
void VecDiff(GFQ_t *v, GFQ_t *y, GFQ_t *x, int32_t len);



//-- load sparse from file to a_matrix_GFQ *A
void load_A_GFQ(FILE *fpA, a_matrix_GFQ *A);
void load_G_GFQ(FILE *fpA, g_matrix_GFQ *A);

//-- alloc for "Quantum" BP Ctrl Blk
void alloc_QBPC(a_matrix_GFQ *A, QBP_Ctl *bp);



//-- "Quantum" BP init: set bp->tt as all zeros. return 1 if zz pass; otherwise INIT *bp and return 0
uint8_t Qbp_init20(QBP_Ctl *bp , GFQ_t *zz , a_matrix_GFQ *A , double p_ch, double afp);  // zz: syndrome, p_ch: channel err rate, afp: see AFP

//-- "Quantum" BP init WITH "Random Perturbation" (RAND_PERTURB) , based on Qbp_init20
uint8_t Qbp_init31(QBP_Ctl *bp , GFQ_t *zz , a_matrix_GFQ *A , double p_ch, double afp); // zz: syndrome, p_ch: channel err rate, afp: see AFP

//-- "Quantum" BP, parallel update
uint8_t Qbp_dec24(QBP_Ctl *bp, a_matrix_GFQ *A);  // copy from bp_dec24 and modify
//-- "Quantum" BP, serial update (along variable nodes)
uint8_t Qbp_dec20(QBP_Ctl *bp, a_matrix_GFQ *A);  // copy form Qbp_dec24 and modify
//-- "Quantum" BP, serial update (along check nodes)
uint8_t Qbp_dec44(QBP_Ctl *bp, a_matrix_GFQ *A);  // copy form Qbp_dec20 and modify








#if USE_GF2_DEC
typedef struct {  // used with a_matrix
  double  p_ch;     // channel error probability (BSC)
  //uint8_t z[M];     // syndrome vector generated before decoding
  uint8_t target_z[M];    // target syndrome vector to converge

  double pn[2*N][2];  // probability (for =0 and =1) of variable nodes (vertical step init)
  // HORIZONTAL:
  //double rm[M][2];  // row check nodes (horizontal step)  // keep for comments
  double dm[M];     // dm[m] = rm[m][0] - rm[m][1]          // (MacKay '99 IT paper's delta_r)
  double *dml[M];   // to alloc as d[m][l] with 0 <= l < num_m[m]  (MacKay's (49))
  double (*rml[M])[2];  // to alloc as rml[m][l] with 0 <= l < num_m[m], each rml[m][l] buffers r_ml^0 and r_ml^1
  double *dqml[M];      // dqml[m][l] = q_ml^0 - q_ml^1 (see MacKay's (47)-(51))
  // VERTICAL:          // and we will qlm to denote q_ml (see MacKay's (47)-(51))
  double (*qlm[2*N])[2];  // to alloc as qlm[l][m] with 0 <= m < num_n[l], each qlm[m][l] buffers q_ml^0 and q_ml^1
  //double *dlm[2*N];       // dlm[l][m] = qlm[l][m][0] - qlm[l][m][1]
  double qn[2*N][2];  // updated probability (for =0 and =1) of variable nodes (vertical step)
  //double dn[2*N];     // updated d[n] = qn[n][0] - qn[n][1]

  uint8_t tt[2*N];    // updated decoded vector (estimated TX vector after BP)
  uint8_t zz[M];    // updated syndrome vector
} BP2_Ctl ;  // Belief Propagation (Sum-Product, Message-Passing) Control block


//-- Matrix vector multiply Ax = y
void Mat2VecMultiply(a2_matrix *A, uint8_t *x, uint8_t *y);
//-- generate syndrome Hr = z
#define GenSyndrome2(H, r, z)  Mat2VecMultiply(H, r, z);

//-- load sparse from file to a_matrix *A
void load_A2(FILE *fpA, a2_matrix *A);

//-- alloc for BP Ctrl Blk
void alloc_BPC2(a2_matrix *A, BP2_Ctl *bp);


//-- BP init and gen syndrome: bias => A->tt => syndrome chk => return 0 if pass; otherwise init bp and return 1
//uint8_t bp2_init20(BP2_Ctl *bp, double bias[][2], a2_matrix *A, uint8_t mode);  // mode = 'b'ias vec or 's'yndrome decoding
uint8_t bp2_init22(BP2_Ctl *bp, double bias[][2], a2_matrix *A, uint8_t *target_zz);

//-- BP dec (one call one iter): return 1 syndrome OK, 0 still NG
uint8_t bp2_dec20(BP2_Ctl *bp, a2_matrix *A);
//-- base on bp_dec4 but decode to target_z (support SYNDROME_DEC)
uint8_t bp2_dec24(BP2_Ctl *bp, a2_matrix *A);

#endif  // USE_GF2_DEC
