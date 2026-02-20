#include "./ldpc_parm.h"
#include "./lib_rand/lib_rand.h"
#include "./lib_math/fast_math.h"
#include "./bp_dec/bp_dec.h"
#include "./bp_dec/bp_llr.h"
#include "./OSD/OSD.h"




#include <math.h>
#include <time.h>
#include <inttypes.h>
#define min(a, b) ((a<b)?a:b)



/*
        MAIN
                     */
int main(void)
{


 printf("Step 1\n");

clock_t t_start = clock();




    // ###########################################################################################
    int  osdw=OSDW; // the OSD order defined by the user
    int max_iter = MAX_ITER;  // maximum number of iterations

    double alpha_a = 0, alpha_b = AFP, alpha= AFP;
    double afp=0;
    afp = (alpha==0 || alpha==100) ? 0 : 100.0/alpha;



    uint8_t *finalResult;
    uint8_t hardD, q;
    double sum;


    double **prob;   // should send to OSD
    uint8_t *LastRun; // should send to OSD


    prob = calloc(N, sizeof(*prob));
    LastRun = calloc(N, sizeof(*LastRun));
    for(int i = 0; i<N; i++){prob[i]= calloc(4, sizeof(*prob[i])); }

    OSD osdDecoder, *osdDec = &osdDecoder;



  int32_t iter, b, m, target;
  uint8_t *usedTo = calloc(N, sizeof(*usedTo));


  GFQ_t *nn   = calloc(N, sizeof(*nn)); // noise vector
  GFQ_t *diff = calloc(N, sizeof(*diff));
  GFQ_t *zz   = calloc(M, sizeof(*zz)); //syndrome vector
  GFQ_t *zz_G = calloc(M_G, sizeof(*zz_G));
  uint64_t tx_seed[4];

  //- RX
  double p_err, rnd_val=0, p_bias;                  // if(rnd_val) { } // just to prevent compile warning when CH_TYPE not 0



  uint8_t syndrome_ok, init_syndro_ok;

    rnd256_init();    rnd256_init_priv(tx_seed);


//===================================================
//Load the parity check matrix and logical matrix
//===================================================

  FILE *fpA = fopen ( PATH_A , "r" );
  FILE *fpGs = fopen( PATH_Gs, "r" );
  a_matrix_GFQ  Amtx, *A=&Amtx;
  g_matrix_GFQ  Gmtx, *G=&Gmtx;
  load_A_GFQ(fpA, A);
  load_G_GFQ(fpGs, G);
  fclose(fpA);
  fclose(fpGs);



  #if LLR_BP == 0
  QBP_Ctl       BPC , *bp=&BPC ;
  alloc_QBPC(A, bp);   // alloc "Quantum" BP decoder needed resource
  #else
  LBP_Ctl       BPC , *bp=&BPC ;
  alloc_LBPC(A, bp);   // alloc LLR-based Quantum BP decoder resource
  #endif





 printf("Step 2\n");






  uint32_t nErrBit,  nErrQbit;
  double   bler, ber, far,  LER, qber, qfar,  chk_bler, chk_ber,  aiter;
  double rStep1=0, rStep2=0, tStep1=0, tStep2=0;
                    //far=false rate                              //aiter=avg iter

 uint32_t n_ill=0; //#of illegal RSR

  time_t  tStr, tEnd, sec, sec_last=0;  // care that: gcc32 has year 2038 problem (gcc64 no this problem)


  tStr = time(NULL);

    initOSD(osdDec);
    osdDec->cnt_enter_osd = 0;
    osdDec->cnt_stage1_fail = 0;
    osdDec->cnt_stage2_fail = 0;
    osdDec->cnt_osd0 = 0;


    FILE *fp;
    fp = fopen( PATH_A , "r" );
    load_A_OSD(osdDec, fp);
    fclose(fp);
    osdDec -> RankH  = M;
    osdDec -> osdw = osdw;





//*********************************************************************
// Load DEM check matrix, logical matrix, and dem error probabilities
//*********************************************************************
  printf("Step 3\n");

p_err = P_ERR;  //P_ERR: predefined error rate



double *p_dem = NULL;
int    N_dem = 0;
int    cap = 1024;

p_dem = (double*)malloc(cap * sizeof(double));


char filename[256];
if (rev ==0 ){
    snprintf(filename, sizeof(filename),
         "codes/p_dem_d%d_r%d_n%.3gz.txt", distance, distance, p_err);
}else{
    snprintf(filename, sizeof(filename),
         "codes/p_dem_d%d_r%d_n%.3gz_rev.txt", distance, distance, p_err);
}


FILE *fpp = fopen(filename, "r");
if(fpp == NULL){
    printf("Cannot open %s\n", filename);
    exit(1);
}

while(1){
    double tmp;
    if(fscanf(fpp, "%lf", &tmp) != 1) break;

    if(N_dem >= cap){
        cap *= 2;
        p_dem = (double*)realloc(p_dem, cap * sizeof(double));
    }
    p_dem[N_dem++] = tmp;
}
fclose(fpp);

printf("Loaded %d DEM probabilities\n", N_dem);

// ================================
// Monte Carlo Simulation for Quantum Error Correction
// BP + OSDw decoding + alpha optimization
// ================================



printf("Step 4\n");

printf("\n=====decode =====\n");

printf("\n p_err=%e\n", p_err);


//int OSD_call=0;
//int BP_unsat=0;
uint32_t LE_cnt=0;
uint32_t LE_Max=LE_tar;
uint32_t shot_max=Shot_max;
uint32_t sample=0;
uint32_t total_iteration=0;
uint32_t total_iter_BPnoOSD=0;
uint32_t total_iter_OSD=0;

uint32_t cnt_0syndrome=0;

uint32_t totalSneaky = 0;

double time_OSD = 0.0;
double time_BP = 0.0;


for (; LE_cnt< LE_Max &&sample<shot_max;sample++){
clock_t last_print_time = clock();

if(sample % ShowTime == 0){
    clock_t t_now = clock();
    double elapsed = (double)(t_now - t_start) / CLOCKS_PER_SEC;

    printf("Samples = %d,  LE = %d,  LER = %.3e,  Time = %.1f s\n",
           sample, LE_cnt, (sample>0)? (1.0*LE_cnt/sample) : 0.0, elapsed);
    fflush(stdout);
}





/* 1. Generate ONE error */
for (b = 0; b < N; b++) {
    double r = rv_UnifOne();
    if(r < p_dem[b]) nn[b] = 1;
    else             nn[b] = 0;
}

/* 2. Syndrome calculation */
Quan_GenSyndrome(A, nn, zz);

/* 3. BP init */
//p_bias = p_err;

#if LLR_BP==0
init_syndro_ok = Qbp_init20(bp, zz, A, p_dem, afp);
#else
init_syndro_ok = Lbp_init(bp, zz, A, p_bias, afp);
#endif

syndrome_ok = init_syndro_ok;
iter = 0;





/* 4. BP iteration */

if(syndrome_ok){
    cnt_0syndrome++;
}

    clock_t t3 = clock();

while(!syndrome_ok && iter < max_iter){
    iter++;
    syndrome_ok = Qbp_dec20(bp, A);
//#if BY_DEC==20
//    syndrome_ok = Qbp_dec20(bp, A);
//#elif BY_DEC==24 || BY_DEC==25
//    syndrome_ok = Qbp_dec24(bp, A);
//#endif
}
    clock_t t4 = clock();
    time_BP += (double)(t4 - t3) / CLOCKS_PER_SEC;

    total_iteration+=iter;
//    total_iter_BPnoOSD+=iter;

/* 5. OSD if needed */


if(osdw == -1 || syndrome_ok){
    finalResult = bp->tt;
}else{
    //OSD_call++;

    for(b=0;b<N;b++){
        sum=0; for(q=0;q<Q;q++) sum+=bp->qn[b][q];
        for(q=0;q<Q;q++) prob[b][q]=bp->qn[b][q]/sum;
    }

    //========================================
    //Reliable Subset Reduction
    //========================================




    clock_t t0 = clock();
    if(osdw == -2){
        finalResult = post_decodeOSDfull(osdDec, zz, (const double**)prob, LastRun, max_iter);
    }else{
        finalResult = post_decodeOSDw(osdDec, zz, (const double**)prob, LastRun);
    }
//    printf("sneaky=%d\n",osdDec->numOfSneaky);
    totalSneaky += osdDec->numOfSneaky;

    clock_t t1 = clock();
    time_OSD += (double)(t1 - t0) / CLOCKS_PER_SEC;

//    total_iter_BPnoOSD-=iter;
    total_iter_OSD+=iter;
}



// �p�� diff = finalResult XOR nn
VecDiff(diff, finalResult, nn, N);

// �ˬd diff �b�޿� qubit �W�� syndrome
Quan_DegSyndrome(G, diff, zz_G);






if(is_zero_vec(zz_G, M_G)) {
    //printf("Decoding SUCCESS (no logical error)\n");
} else {
    LE_cnt++;
//    printf(" %d",i);
    //printf("X LOGICAL FAILURE\n");
}

}


total_iter_BPnoOSD=total_iteration-total_iter_OSD;

double total_time_sec = 0.0;
clock_t t_end = clock();
total_time_sec = (double)(t_end - t_start) / CLOCKS_PER_SEC;
double avg_time = total_time_sec / sample;   // �C�� decoding �������ɶ� (��)
double avg_iter = 1.0*total_iteration/ sample;
double avg_iter_noOSD = 1.0*total_iter_BPnoOSD/(sample-osdDec->cnt_enter_osd-cnt_0syndrome);


double pL = 1.0 * LE_cnt / sample;
double sigma_k = sqrt((double)LE_cnt);
double err_plus  = 2.0*sigma_k / sample;
double err_minus = 2.0*sigma_k / sample;



FILE *fps;
char outname[256];

if (rev==0 ){
snprintf(outname, sizeof(outname),
         "Results/surface_d%d.txt", distance);
}else{
snprintf(outname, sizeof(outname),
         "Results/surface_d%d_rev.txt", distance);
}

fps = fopen(outname, "a+");
if(fps == NULL){
    printf("Cannot open %s\n", outname);
    exit(1);
}



printf("\np_err=%e LER +err  -err = %e  +%e  -%e\n",
       p_err, pL, err_plus, err_minus);

//printf(" Avg decoding time = %.3e sec    Avg iterations= %.3e \n", avg_time, avg_iter);
//
//printf(" total iterations = %d    \n", total_iteration);
//printf(" total iterations without OSD= %d    \n", total_iter_BPnoOSD);
//printf(" total iterations with OSD= %d    \n", total_iter_OSD);


printf("\nAvg iterations w/o OSD= %.3e  \n", avg_iter_noOSD );
//printf(" Avg iterations no OSD= %.3e     avg iterations OSD= %.3e \n", avg_iter_noOSD, 1.0*(total_iteration-total_iter_BPnoOSD)/osdDec->cnt_enter_osd );
//printf("total time = %.3f sec\n", total_time_sec);
//printf("OSD total time = %.3f sec\n", time_OSD);
//printf("OSD calls = %lld\n", osdDec->cnt_enter_osd);
printf("Time per OSD = %.6f sec\n", time_OSD / osdDec->cnt_enter_osd );
printf("Time per BP iteration = %.6f sec\n",
       time_BP / total_iteration);


        double avgSneaky = (double)totalSneaky / osdDec->cnt_enter_osd;
    printf("Average remaining bits after RSR: %.2f\n", avgSneaky);
    printf("relative code length= %.2f\n", avgSneaky/N );


//printf("\nOSD call=%d  osdDEc>cnt_enter_osd=%d\n", OSD_call,osdDec->cnt_enter_osd);

fprintf(fps,
    "\n%e  %e  %e  %e  \n",
    p_err,    pL,    err_plus,    err_minus
);
fprintf(fps,
    "Time per BP iteration =%.3e\nTime per OSD = %.3e\nAvg iterations w/o OSD=  %.3e\n",
time_BP / total_iteration,   time_OSD / osdDec->cnt_enter_osd,    avg_iter_noOSD
);
fprintf(fps,
    "Sample = %d   LE= %d  OSDw=%d OSD calls= %d\n",
    sample,   LE_cnt,  OSDW, osdDec->cnt_enter_osd
);

fprintf(fps,"Average remaining bits after RSR: %.2f\n", avgSneaky);
fprintf(fps, "Rel_Thr=%.6f  relative code length= %.2f\n", RelThr, avgSneaky/N );
fprintf(fps,"Stage1Fail: %llu times\n", osdDec->cnt_stage1_fail);
fprintf(fps,"Stage2Fail: %llu times\n", osdDec->cnt_stage2_fail);
fprintf(fps,"OSD-0 sufficient: %llu times\n", osdDec->cnt_osd0);
fclose(fps);



printf("\n===== OSD Statistics =====\n");
printf("OSDw=%d\n",OSDW);
printf("Total samples: %llu\n", sample);
printf("Total LER failures: %llu\n", LE_cnt);
printf("Enter OSD: %llu times\n", osdDec->cnt_enter_osd);
printf("Stage1Fail: %llu times\n", osdDec->cnt_stage1_fail);
printf("Stage2Fail: %llu times\n", osdDec->cnt_stage2_fail);
printf("OSD-0 sufficient: %llu times\n", osdDec->cnt_osd0);
//printf("LER = %.6e\n", pL);



    free(usedTo);
    free(nn);
    free(diff);
    free(zz);
    free(zz_G);
	printf("\nMain function done!\n");


    if(osdDec != NULL)
        free_OSD(osdDec);


	getchar();
	return 0;
}
