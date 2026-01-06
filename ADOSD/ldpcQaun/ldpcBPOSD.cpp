#include <math.h>
#include "./ldpc_parm.h"
//#include "./lib_mat/lib_mat.h"
#include "./lib_rand/lib_rand.h"
#include "./bp_dec/bp_dec.h"
#include "./bp_dec/bp_llr.h"
#include "./OSD/Decoder.h"
#include "./OSD/OSDMatrix.h"
#include "math.h"

//################ local functions ################//
#if ITER_LOG
FILE *fpLogE,  *fpLogC ;    FILE *_iter_log,  *_iter_chk ;
static void GenIterLog(a_matrix_GFQ *A, QBP_Ctl *bp)
{
  int32_t  m, b, q; //, n;
  double  sum; //tmp_sum;

  //-- log qn
  /*for(n=0; n < A->NN; n++) {  tmp_sum = bp->pn[n][0]+bp->pn[n][1]+bp->pn[n][2]+bp->pn[n][3];
    fprintf(_iter_log, "%g  %g  %g  %g\n", bp->pn[n][0]/tmp_sum, bp->pn[n][1]/tmp_sum, bp->pn[n][2]/tmp_sum, bp->pn[n][3]/tmp_sum);
  }*/
  for ( b = 0 ; b < N ; b ++ ) {
    sum=0;
    for(q=0;q<Q;q++) { sum += bp->qn[b][q]; }
    for ( q = 0 ; q < Q ; q ++ ) { fprintf(_iter_log, "%g  ", bp->qn[b][q]/sum); }
    fprintf(_iter_log, "\n");
  }

  //-- log dm_c
  /*for(m=0; m < A->MM; m++) { fprintf(_iter_chk, "%g  ", bp->dm_c[m]); }
  fprintf(_iter_chk, "\n");*/
  //-- log dm_d (= 1^z_m * dm_c)
  for(m=0; m < A->MM; m++) { fprintf(_iter_chk, "%g  ", bp->dm_c[m]*((double)bp->target_z[m]*(-2) + 1)); }
  fprintf(_iter_chk, "\n");
}
#endif


/*
        MAIN
                     */
int main1(void)
{
    double invProbTemp1[N];
    double invProbTemp2[N];
    double bEntropyTemp1[N];
    double bEntropyTemp2[N];
    bool staticPlace[N];
    bool checkNonSATvar[N];
    double **prob = new double* [N];    // should send to OSD
    uint8_t *lastFor = new uint8_t [N]; // should send to OSD
    uint8_t usedTo[N];
    uint8_t hardD;
    for(int i = 0; i<N; i++) {usedTo[i] = 0; lastFor[i] = 0; staticPlace[i] = true;}
    for(int i = 0; i<N; i++) prob[i] = new double[4];
    Decoder decoder;
    decoder.addParityMatrix(PATH_A);
    decoder.setRank(K);
    decoder.setOSDw(OSDW);
    FILE *fpRecord = fopen("../record_score", "w");

    double invProb, errInvProb, nonSTCInvProb, bEntropy, errBEntropy, nonSTCBEntropy, checkNonSATInvProb, checkNonSATBEntropy, x, x_1, x_2, x_3;
    uint32_t hammingW, nonSTCplaceW = 0, checkNonSATW, checkNonSATvarW;
    GFQ_t diff[N], zz_G[M_G], zz_temp[M];
//double p_nos = pow(10.0, (-10.0  -  10.0/8 *10)/10);   //#define P_START (-10.0  -  10.0/8 *10)    // 0.0056
//double p_nos = (double)0.013;
//double p_nos = (double)0.01;

//double p_nos = (double)0.1;
//double p_nos = (double)0.2;
//double p_nos = (double)0.1;     // BSC error probability
//double p_nos = pow(10.0, -14.9/10);                   // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = pow(10.0, -14.99/10);    // [[5,1,3]] code, +1 redundant check row (hit Degenerate Err)
//double p_nos = pow(10.0, -15.0/10);     // [[5,1,3]] code, +1 redundant check row (hit Degenerate Err)
//double p_nos = pow(10.0, -15.01/10);                   // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
double p_nos = (double)0.031;                         // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = pow(10.0, -15.1/10);                   // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = (double)0.0178;                        // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = (double)0.0177828;                     // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = pow(10.0, -17.499998559187063/10);     // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = pow(10.0, -17.5/10);     // [[5,1,3]] code, +1 redundant check row (hit Degenerate Err)
//double p_nos = pow(10.0, -17.52/10);    // [[5,1,3]] code, +1 redundant check row (hit Degenerate Err)
//double p_nos = pow(10.0, -17.6/10);                   // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = (double)0.01;
//double p_nos = pow(10.0, -2.0 -1.0/8);
//double p_nos = pow(10.0, -2.0 -3.0/8);
//double p_nos = (double)0.003;
//double p_nos = (double)0.0031;
//double p_nos = (double)2.0e-4;
//double p_nos = (double)0.00017778;                    // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = (double)5.62341e-005;    // [[5,1,3]] code, +1 redundant check row (hit Degenerate Err)
//double p_nos = pow(10.0, -42.5/10);     // [[5,1,3]] code, +1 redundant check row (hit Degenerate Err)
//double p_nos = pow(10.0, -42.6/10);     // [[5,1,3]] code, +1 redundant check row (hit Degenerate Err)
//double p_nos = pow(10.0, -44.0/10);                   // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = (double)3.16228e-005;                  // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = pow(10.0, -44.999996786570648/10);     // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = pow(10.0, -45.0/10);     // [[5,1,3]] code, +1 redundant check row (hit Degenerate Err)
//double p_nos = pow(10.0, -46.0/10);     // [[5,1,3]] code, +1 redundant check row (hit Degenerate Err)
//double p_nos = pow(10.0, -46.7/10);     // [[5,1,3]] code, +1 redundant check row (hit Degenerate Err)
//double p_nos = pow(10.0, -47.0/10);                  // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = pow(10.0, -50.0/10);                  // [[5,1,3]] code, +1 redundant check row (hit Exact Err)

  int test_prt = 1 , iter_prt = 1;

  int32_t iter;
  int32_t b , q, m; // dummy variables

  //- Read file only
  uint8_t rr[N], zz[M], syndrome_ok;
  FILE *fpA = fopen ( PATH_A , "r" ) ;
  FILE *fpGs = fopen ( PATH_Gs , "r" ) ;
  a_matrix_GFQ  Amtx, *A=&Amtx;
  g_matrix_GFQ  Gmtx, *G=&Gmtx;

  double afp;

  #if ITER_LOG
  char W[4]={'I','X','Z','Y'};  // for log error patterns during iterations
  fpLogE = fopen ("_Log_E1", "w");  fpLogC = fopen ("_Log_C1", "w");
  _iter_log = fopen ("_iter_log", "w");  _iter_chk = fopen("_iter_chk_d", "w"); // _chk: use dm ; _chk_c: use dm_c ; _chk_s: use dm_d
  #endif

  //- init
  rnd256_init();
  //rnd_max_show();

  //- load sparse from file to a_matrix_GFQ *A
  load_A_GFQ(fpA, A);
  load_G_GFQ(fpGs, G);

  #if LLR_BP == 0
  double sum;
  QBP_Ctl    BPC , *bp=&BPC ;
  alloc_QBPC(A, bp);   // alloc "Quantum" BP decoder needed resource
  #else
  LBP_Ctl    BPC , *bp=&BPC ;
  alloc_LBPC(A, bp);   // alloc LLR-based Quantum BP decoder resource
  #endif


  //- force TX=all0 and then add noise as vector rr
  memset(rr, 0, sizeof(rr));  // force all0
#if 1
  //rr[0] = 1;  // set first qubit error be X     // make [[5,1,3]] Qbp_dec20 FAIL when p=0.0031 , a=1.6
  //rr[0] = 2;  // set first qubit error be Z
  //rr[0] = 3;  // set first qubit error be Y
  //rr[1] = 1;  // set 2nd qubit error be X       // Qbp_dec20 FAIL when p=0.0031 , a=1.6;    make [[5,1,3]] Qbp_dec24 FAIL when p=0.0178 (+ap=13 PASS, +stopC FAIL again); and dec_44 FAIL
  //rr[1] = 2;  // set 2nd qubit error be Z
  //rr[1] = 3;  // set 2nd qubit error be Y
  //rr[2] = 1;  // set 3rd qubit error be X
  //rr[2] = 2;  // set 3rd qubit error be Z       // make [[5,1,3]] Qbp_dec20 FAIL when p=0.0031 , a=1.6; and dec_44 FAIL
  //rr[2] = 3;  // set 3rd qubit error be Y
  //rr[3] = 1;  // set 4th qubit error be X
  //rr[3] = 2;  // set 4th qubit error be Z
  //rr[3] = 3;  // set 4th qubit error be Y       // make [[5,1,3]] Qbp_dec20 FAIL
  //rr[4] = 1;  // set 5th qubit error be X
  //rr[4] = 2;  // set 5th qubit error be Z       // make [[5,1,3]] Qbp_dec20 FAIL when p=0.0031 , a=1.6; and dec_44 FAIL
  //rr[4] = 3;  // set 5th qubit error be Y       // make [[5,1,3]] Qbp_dec20 FAIL


  //rr[24] = 2; rr[29] = 2; rr[33] = 2; rr[41] = 2; rr[42] = 2;
  //rr[0] = 3; rr[12] = 2; rr[15] = 3; rr[23] = 2; rr[24] = 2; rr[60] = 2; rr[80] = 3; rr[97] = 3; rr[98] = 2; // false alarm
  //rr[10] = 1; rr[28] = 3; rr[30] = 1; rr[49] = 1; rr[70] = 1; rr[89] = 1;
  //rr[6] = 3; rr[21] = 1; rr[29] = 1; rr[30] = 3; rr[33] = 1; rr[49] = 3; rr[55] = 1;
  //rr[57] = 2; rr[79] = 2; rr[81] = 2; rr[93] = 2; rr[95] = 2; rr[97] = 3;
  //rr[9] = 1; rr[11] = 1; rr[12] = 1; rr[86] = 1; rr[97] = 1;    // false alarm, but alpha = 95 solve it
//rr[5] = 1; rr[25] = 1; rr[34] = 1; rr[45] = 1; rr[54] = 3; rr[88] = 3;
rr[10] = 1; rr[28] = 3; rr[30] = 1; rr[49] = 1; rr[70] = 1; rr[89] = 1;

#endif

  for ( b = 0 ; b < N ; b ++ )  printf("%u ", rr[b]);   //error pattern
  printf(" = added error pattern \n");

  //-- generate syndrome Ar = z
  Quan_GenSyndrome(A, rr, zz);

  if(test_prt) {
    for ( b = 0 ; b < M ; b ++ ) { printf("%d ", zz[b]); }  //syndrome
    printf(" = RX syndrome zz \n");
  }

  printf("CW dec start ... \n");

  afp = (AFP==0)?  0 : (double)100.0/AFP;
  if(AFP2) { printf("\n AFP2 != 0 , it is not supported in this test mode! exit now!!!!\n"); exit(1);}

  //-- BP init
#if LLR_BP == 0
  syndrome_ok = Qbp_init20(bp, zz, A, p_nos, afp);     /////////////////////////////////////
  //syndrome_ok = Qbp_init31(bp, zz, A, p_nos, afp);   /////////////////////////////////////
  for ( b = 0 ; b < N ; b ++ ) {
    for ( q = 0 ; q < Q ; q ++ ) { printf("%.2f ", bp->pn[b][q]); }   printf(", ");
  } printf(" = init bp->pn \n");
#else
  syndrome_ok = Lbp_init(bp, zz, A, p_nos, afp);
  for ( b = 0 ; b < N ; b ++ ) {
    for ( q = 0 ; q < Q-1 ; q ++ ) { printf("%g ", bp->LA[b][q]); }   printf(" = init bp->LA[%d] \n", b);
  }
#endif // LLR_BP


  printf("BP init done, checked syndrome_ok = %u \n", syndrome_ok);
  printf("press enter to continue...\n"); getchar();

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start
  #if ITER_LOG    // initial log here
  for(b=0; b<N; b++) { if(rr[b]!=0) fprintf(fpLogE, "%c%d ", W[rr[b]],b+1); }     fprintf(fpLogE, "\n");
  for(b=0; b<M; b++) { if(zz[b]!=0) fprintf(fpLogC, "%d ", b+1); }                fprintf(fpLogC, "\n");
  GenIterLog(A, bp);
  #endif

  //-- BP dec
  iter = 0;
  while ( !syndrome_ok ) {
  //while (iter < MAX_ITER) {
    iter ++;
  #if BY_DEC==20
    syndrome_ok = Qbp_dec20(bp, A);
  #elif BY_DEC==24 || BY_DEC==25 //|| BY_DEC==26
    syndrome_ok = Qbp_dec24(bp, A);
  #elif BY_DEC==44
    syndrome_ok = Qbp_dec44(bp, A);
  #elif BY_DEC==60
    syndrome_ok = Lbp_dec60(bp, A);
  #elif BY_DEC==64
    syndrome_ok = Lbp_dec64(bp, A);
  #elif BY_DEC==84
    syndrome_ok = Lbp_dec84(bp, A);
  #else
    XXX_NG_CFG_XXX
  #endif

  #if ITER_LOG  // per iteration log here
    for(b=0; b<N; b++) { if(bp->tt[b]!=0) fprintf(fpLogE, "%c%d ", W[bp->tt[b]],b+1); }     fprintf(fpLogE, "\n");
    for(b=0; b<M; b++) { if(bp->zz[b]!=0) fprintf(fpLogC, "%d ", b+1); }                    fprintf(fpLogC, "\n");
    GenIterLog(A, bp);
  #endif



    if(iter_prt) {
      printf("  Run Iter %d done: syndrome_ok = %u    ", iter, syndrome_ok);
      printf("  Curr EstErr to Actual Err HamDist %d  is_hit_c %u\n", HamDist(rr, bp->tt, N), bp->is_hit_c);
      for ( b = 0 ; b < N ; b ++ ) { printf("%d ", bp->tt[b]); }  printf(" = estimated bp->tt \n");
      for ( b = 0 ; b < M ; b ++ ) { printf("%d ", bp->zz[b]); }  printf(" = syndrome bp->zz \n");

    #if LLR_BP == 0
      for ( b = 0 ; b < N ; b ++ ) {
        sum=0; for(q=0;q<Q;q++) { sum += bp->qn[b][q]; }
        for ( q = 0 ; q < Q ; q ++ ) {
            prob[b][q] = bp->qn[b][q]/sum;
            printf("%.2f ", prob[b][q]);
        }   printf(", ");
      } printf(" = bp->qn (normalized) \n");
      //getchar();  // STOP per ITERATION
    #else
      for ( b = 0 ; b < N ; b ++ ) {
        for ( q = 0 ; q < Q-1 ; q ++ ) { printf("%g ", bp->GA[b][q]); }   printf(" = bp->La[%d] \n", b);
      }
      //getchar();  // STOP per ITERATION
    #endif
    }


    invProb=0; errInvProb=0; nonSTCInvProb=0; bEntropy=0; errBEntropy=0; nonSTCBEntropy=0; checkNonSATInvProb = 0; checkNonSATBEntropy = 0;
    hammingW=0, checkNonSATW = 0, checkNonSATvarW = 0;

    Quan_GenSyndrome(A, bp->tt, zz_temp);
    VecDiff(diff, zz_temp, zz, M);
    memset(checkNonSATvar, 0, sizeof(checkNonSATvar));
    for(m =0; m<M; m++){
        if(diff[m]!=0){
            checkNonSATW++;
            b = A->num_m[m];
            for(q = 0; q<b; q++ )
                    checkNonSATvar[A->mlist[m][q]] = true;
        }

    }

    for(b = 0; b<N; b++){
            hardD = 0;
            for(q = 1; q<=3; q++){
                if(prob[b][q] > prob[b][hardD])
                    hardD = q;
            }

            if(hardD == usedTo[b])
                lastFor[b]++;
            else{
                lastFor[b] = 1;
                usedTo[b] = hardD;
            }

            //================calculate the score =================
            x = log2(prob[b][hardD]);
            if(prob[b][hardD]!=1)
                x_1 = log2(1-prob[b][hardD]);

            invProb -= x;
            if(prob[b][hardD]!=0 && (prob[b][hardD]!=1) )
                bEntropy -= (prob[b][hardD]* x)+(1-prob[b][hardD])* x_1;
            if(hardD != 0){
                hammingW+=1;
                errInvProb -= x;
                if(prob[b][hardD]!=0 && (prob[b][hardD]!=1) )
                    errBEntropy -= (prob[b][hardD]* x)+(1-prob[b][hardD])* x_1;
            }
            if(checkNonSATvar[b]){
                checkNonSATvarW++;
                checkNonSATInvProb -= x;
                if(prob[b][hardD]!=0 && (prob[b][hardD]!=1) )
                    checkNonSATBEntropy -= (prob[b][hardD]* x)+(1-prob[b][hardD])* x_1;
            }
            //=============================================


    }
    invProb /= N;
    bEntropy /= N;
    errInvProb /= hammingW;
    errBEntropy /= hammingW;
    if(checkNonSATvarW == 0){checkNonSATBEntropy = 0; checkNonSATInvProb= 0;}
    else{
        checkNonSATBEntropy /= checkNonSATvarW;
        checkNonSATInvProb /= checkNonSATvarW;
    }



    //============== calculate non static part =================
    if(iter == 1){
            for(b = 0; b<N; b++) {
                    x = log2(prob[b][hardD]);
                    if(prob[b][hardD]!=1)
                            x_1 = log2(1-prob[b][hardD]);

                     invProbTemp1[b] = x;
                    if(prob[b][hardD]!=0 && (prob[b][hardD]!=1) )
                        bEntropyTemp1[b] = (prob[b][hardD]* x)+(1-prob[b][hardD])* x_1;
            }

    }
    else if(iter == 2){
            for(b = 0; b<N; b++) {
                    x = log2(prob[b][hardD]);
                    if(prob[b][hardD]!=1)
                            x_1 = log2(1-prob[b][hardD]);

                     invProbTemp2[b] = x;
                    if(prob[b][hardD]!=0 && (prob[b][hardD]!=1) )
                        bEntropyTemp2[b] = (prob[b][hardD]* x)+(1-prob[b][hardD])* x_1;
            }

    }
    else if(iter == 3){
            x = 0;  x_1 = 0;  x_2 = 0; x_3 = 0;
            for(b=0; b<N; b++){
                if(lastFor[b]<3){
                    nonSTCplaceW ++;
                    staticPlace[b] = false;
                    x -= invProbTemp1[b];
                    x_1 -= bEntropyTemp1[b];
                    x_2 -= invProbTemp2[b];
                    x_3 -= bEntropyTemp2[b];
                }
            }
           fprintf(fpRecord, "ITER_1_nonStaticInverseOfProbability: %lf \n", (nonSTCplaceW == 0)?0:x/nonSTCplaceW);
            fprintf(fpRecord, "ITER_1_nonStaticBinaryEntropy: %lf \n", (nonSTCplaceW == 0)?0:x_1/nonSTCplaceW);
            fprintf(fpRecord, "ITER_2_nonStaticInverseOfProbability: %lf \n", (nonSTCplaceW == 0)?0:x_2/nonSTCplaceW);
            fprintf(fpRecord, "ITER_2_nonStaticBinaryEntropy: %lf \n", (nonSTCplaceW == 0)?0:x_3/nonSTCplaceW);
            fprintf(fpRecord, "nonStaticPlaceWeight: %d \n", nonSTCplaceW);


            printf("ITER_1_nonStaticInverseOfProbability: %lf \n", (nonSTCplaceW == 0)?0:x/nonSTCplaceW);
            printf("ITER_1_nonStaticBinaryEntropy: %lf \n", (nonSTCplaceW == 0)?0:x_1/nonSTCplaceW);
            printf("ITER_2_nonStaticInverseOfProbability: %lf \n", (nonSTCplaceW == 0)?0:x_2/nonSTCplaceW);
            printf("ITER_2_nonStaticBinaryEntropy: %lf \n", (nonSTCplaceW == 0)?0:x_3/nonSTCplaceW);
            printf("nonStaticPlaceWeight: %d\n", nonSTCplaceW);
    }

    if (iter>=3){
            for( b = 0; b<N; b++){
                if(!staticPlace[b]){
                    x = log2(prob[b][usedTo[b]]);
                    if(prob[b][usedTo[b]]!=1)
                            x_1 = log2(1-prob[b][usedTo[b]]);

                    nonSTCInvProb -= x;
                    if(prob[b][usedTo[b]]!=0 && (prob[b][usedTo[b]]!=1) )
                         nonSTCBEntropy -= (prob[b][usedTo[b]]* x)+(1-prob[b][usedTo[b]])* x_1;
                }
            }
            nonSTCInvProb/=nonSTCplaceW;
            nonSTCBEntropy/=nonSTCplaceW;
    }
    // ===============================================

    fprintf(fpRecord, "\t%d - %u\n", iter, syndrome_ok);
    fprintf(fpRecord, "\t\tInverseOfProbability: %lf \n", invProb);
    fprintf(fpRecord, "\t\tBinaryEntropy: %lf \n", bEntropy);
    fprintf(fpRecord, "\t\tErrorInverseOfProbability: %lf \n", errInvProb);
    fprintf(fpRecord, "\t\tErrorBinaryEntropy: %lf \n", errBEntropy);
    fprintf(fpRecord, "\t\tnonStaticInverseOfProbability: %lf \n", nonSTCInvProb);
    fprintf(fpRecord, "\t\tnonStaticBinaryEntropy: %lf \n", nonSTCBEntropy);
    fprintf(fpRecord, "\t\tcheckNonSATPlaceInverseOfProbability: %lf \n", checkNonSATInvProb);
    fprintf(fpRecord, "\t\tcheckNonSATPlaceBinaryEntropy: %lf \n", checkNonSATBEntropy);
    fprintf(fpRecord, "\t\tHammingWeight: %d\n", hammingW);
    fprintf(fpRecord, "\t\tcheckNonSatisfyWeight: %d\n", checkNonSATW);
    fprintf(fpRecord, "\t\tcheckNonSatisfyVariablePlaceWeight: %d\n\t\t", checkNonSATvarW);

    printf("InverseOfProbability: %lf \n", invProb);
    printf("BinaryEntropy: %lf \n", bEntropy);
    printf("ErrorInverseOfProbability: %lf \n", errInvProb);
    printf("ErrorBinaryEntropy: %lf \n", errBEntropy);
    printf("nonStaticInverseOfProbability: %lf \n", nonSTCInvProb);
    printf("nonStaticBinaryEntropy: %lf \n", nonSTCBEntropy);
    printf("checkNonSATPlaceInverseOfProbability: %lf \n", checkNonSATInvProb);
    printf("checkNonSATPlaceBinaryEntropy: %lf \n", checkNonSATBEntropy);
    printf("HammingWeight: %d\n", hammingW);
    printf("checkNonSatisfyWeight: %d\n", checkNonSATW);
    printf("checkNonSatisfyVariablePlaceWeight: %d\n", checkNonSATvarW);

    decoder.setSyndromeCopy(zz);
    decoder.resetParity();
    uint8_t *finalResult = decoder.post_decode(prob, lastFor);

    hammingW = 0;
    for(int i = 0; i<N; i++){
        cout<<(uint32_t)finalResult[i]<<" ";
        fprintf(fpRecord, "%d ", (uint32_t)finalResult[i]);
        if(finalResult[i] != 0 )    hammingW++;
    }
    fprintf(fpRecord, "\n\t\tFinalHammingWeight: %d\n", hammingW);

    cout<<endl;


//-- "Quantum" degeneracy check
VecDiff(diff, finalResult, rr, N);     // if bp->tt is a degenerate err of nn,
Quan_DegSyndrome(G, diff, zz_G);  // then diff will gen zz_G = all0
if(is_zero_vec(zz_G, M_G)) {
    if(is_zero_vec(diff, N))  {printf("\nDEC OK !! and Hit Exact Error \n"); fprintf(fpRecord, "\n\t\tDEC OK !! and Hit Exact Error \n2\n");}
    else                      {printf("\nDEC OK !! and Hit Degenerate Error !!!! \n"); fprintf(fpRecord ,"\n\t\tDEC OK !! and Hit Degenerate Error !!!! \n1\n");}
}else{
//printf("\nDEC NG !! \n");
Quan_GenSyndrome(A, finalResult, zz_temp);
VecDiff(diff, zz_temp, zz, M);
if(is_zero_vec(diff, M))   {printf("\nDEC NG !! Syndrome match FALSE ALARM !! \n"); fprintf(fpRecord, "\n\t\tDEC NG !! Syndrome match FALSE ALARM !! \n0\n");}
else              {printf("\nDEC NG !! Syndrome NOT match !! \n"); fprintf(fpRecord, "\n\t\tDEC NG !! Syndrome NOT match !! \n\n");}
}


  printf("\n");


    //if ( syndrome_ok || iter==MAX_ITER )   break;     // CODE DEFAULT
    if(iter==MAX_ITER){ break; } else{ syndrome_ok=0; } // TEST FORCE MAX_ITER


  }

  //GFQ_t diff[N], zz_G[M_G];
  //-- "Quantum" degeneracy check
  VecDiff(diff, bp->tt, rr, N);     // if bp->tt is a degenerate err of nn,
  Quan_DegSyndrome(G, diff, zz_G);  // then diff will gen zz_G = all0
  if(is_zero_vec(zz_G, M_G)) {
    if(is_zero_vec(diff, N))  printf("\nDEC OK !! and Hit Exact Error \n");
    else                      printf("\nDEC OK !! and Hit Degenerate Error !!!! \n");
  }else{
    //printf("\nDEC NG !! \n");
    if(syndrome_ok)   printf("\nDEC NG !! Syndrome match FALSE ALARM !! \n");
    else              printf("\nDEC NG !! Syndrome NOT match !! \n");
  }

  fclose(fpA);  fclose(fpGs);
    fclose(fpRecord);
  #if ITER_LOG
  fclose(fpLogE);  fclose(fpLogC);      fclose(_iter_log);  fclose(_iter_chk);
  #endif
	return 0;
}
