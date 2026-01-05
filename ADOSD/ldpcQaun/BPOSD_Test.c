#include "./ldpc_parm.h"
//#include "./lib_mat/lib_mat.h"
#include "./lib_rand/lib_rand.h"
#include "./lib_math/fast_math.h"
#include "./bp_dec/bp_dec.h"
#include "./bp_dec/bp_llr.h"
#include "./OSD/OSD.h"

#include <math.h>
#include <time.h>
#include <inttypes.h>
#define min(a, b) ((a<b)?a:b)

#define MODE (4)
/*  MODE 1 only generate syndrome : OSDw: all
    MODE 2 BP                       OSDw: -1
    MODE 3 (BP+Cal_LF)              OSDw: -1
    MODE 4 (BP+Cal_LF)+OSD          OSDw: -2, 0, 1, 2 */

#define NUM_DATA 500// 5e7
#define NUM_RECORD 10 // 1e6
#define NUM_STORE_TIME 50

int main0(void)
{

// Surface BLER~ 10^-6
double p_nos = 0.033; // d = 15
//double p_nos = 0.025; // d = 13
//double p_nos = 0.017; // d = 11

// ###########################################################################################
    int  osdw=OSDW, max_iter = MAX_ITER, max_iter_1 = MAX_ITER+1;
    double alpha_a = 0, alpha_b = AFP, alpha= AFP;
    scanf("%d", &osdw);
    //osdw = -2;
    //scanf("%lf", &p_nos);

    uint8_t *finalResult;
    uint8_t hardD, q;
    double sum;

    uint8_t usedTo[N];
    double **prob;   // should send to OSD
    uint8_t *lastFor; // should send to OSD
    prob = calloc(N, sizeof(*prob));
    lastFor = calloc(N, sizeof(*lastFor));
    for(int i = 0; i<N; i++){prob[i]= calloc(4, sizeof(*prob[i])); }

    OSD osdDecoder, *osdDec = &osdDecoder;
    FILE *fp;
    clock_t store_time[NUM_STORE_TIME];
    #if OSD_RECORD
    OsdRecord osdRecord, *osdRec = &osdRecord;
    #endif // OSD_RECORD
// ##########################################################################################

  int32_t iter, b, m, target;
  GFQ_t nn[N], diff[N], zz[M], zz_G[M_G];
  uint64_t tx_seed[4];

  //- RX
  double rnd_val=0, p_bias, afp=0;                   if(rnd_val) { } // just to prevent compile warning when CH_TYPE not 0
  FILE *fpA = fopen ( PATH_A , "r" );
  FILE *fpGs = fopen( PATH_Gs, "r" );

  a_matrix_GFQ  Amtx, *A=&Amtx;
  g_matrix_GFQ  Gmtx, *G=&Gmtx;
  uint8_t syndrome_ok, init_syndro_ok;

  //- init proc
  //fast_math_test();
  rnd256_init();    rnd256_init_priv(tx_seed);

  load_A_GFQ(fpA, A);
  load_G_GFQ(fpGs, G);

  #if LLR_BP == 0
  QBP_Ctl       BPC , *bp=&BPC ;
  alloc_QBPC(A, bp);   // alloc "Quantum" BP decoder needed resource
  #else
  LBP_Ctl       BPC , *bp=&BPC ;
  alloc_LBPC(A, bp);   // alloc LLR-based Quantum BP decoder resource
  #endif

  int32_t n_list_sol = 0;  if(n_list_sol){ }    // number of list solutions (and prevent compile warning)

  char st_name[1000] = "";
  char path_St[1000] = "";  //= "St";     //sprintf(path_St,  "St_%u_%u_it%u_dec%u_1e-11%s%s%s", N,K,MAX_ITER,BY_DEC,(FIX_P? "_fixP":""),(ALFA? "_sq2":""),(BETA? "_ofs":""));

  if(CH_TYPE==0)      { } // dep ch, no special mark
  else if(CH_TYPE==1) { strcpy(st_name, "XZ_"); }  // indep X-Z ch
  else                { printf("CH_TYPE %d undefined !!  now exit();", CH_TYPE);  getchar();  exit(1); }

  sprintf(st_name, "%s%s_it%u_dec%u", st_name, PRE_A, max_iter, BY_DEC);  // care: st_name in both I/O , behavior not guaranteed
  if(USE_GF2_DEC) sprintf(st_name, "%s_%s", st_name, "GF2");        // care: st_name in both I/O , behavior not guaranteed
    // also cfg Qbp_init20(): rnd means random order (Fisher-Yates approach) , rndS means Sattolo algorithm

  if(ALFA)    sprintf(st_name, "%s_%s%u", st_name, "ac", ALFA);     // care: st_name in both I/O , behavior not guaranteed

  if(alpha) sprintf(st_name, "%s_%s%.0lf", st_name, "ap", alpha);      // care: st_name in both I/O , behavior not guaranteed

  // --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if(osdw >= 0){sprintf(st_name, "%s_OSDW%d", st_name, osdw);}
  else if(osdw == -1){sprintf(st_name, "%s_NoOSD", st_name);}
  else if(osdw == -2){sprintf(st_name, "%s_OSDfull", st_name);}


  // --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  sprintf(path_St,  "St_Mode%d_%s_pnos_%d", MODE, st_name, (int)(p_nos*1000));

  uint64_t nCW, nIterAcc;
  double  aiter;

    clock_t  tStr, tEnd, sec, currSec, tLast;  // care that: gcc32 has year 2038 problem (gcc64 no this problem)

    nCW=0; nIterAcc=0;
    printf("//==== init p_nos = %g ====// N = %u K = %u   It = %u  BY_DEC %u \n", p_nos,N,K,MAX_ITER,BY_DEC);

    initOSD(osdDec);
    #if OSD_RECORD
    initOsdRecord(osdRec);
    #endif // OSD_RECORD
    fp = fopen( PATH_A , "r" );
    load_A_OSD(osdDec, fp);
    fclose(fp);
    osdDec -> RankH  = N- K;
    osdDec -> osdw = osdw;
    alpha = alpha_a * log10(p_nos) * 10 + alpha_b;
    if(alpha <= 50)    alpha = 0;
    printf("\nNew alpha: %lf\n\n", alpha);

    tLast = tStr = clock();

  do {
    //- add noise nn_(Nx1), where (Nx1) means a column vector

    for ( b = 0 ; b < N ; b ++ ) {
        #if CH_TYPE==0    // dep ch
          rnd_val = rv_UnifOne();
          if(rnd_val < p_nos/3)        { nn[b] = 3; }
          else if(rnd_val < p_nos/3*2) { nn[b] = 2; }
          else if(rnd_val < p_nos)     { nn[b] = 1; }
          else                         { nn[b] = 0; }
        #elif CH_TYPE==1  // indep X-Z ch
          nn[b] = 0;
          if(rv_UnifOne() < p_nos)  { nn[b] |= 0b01; }  // set X error
          if(rv_UnifOne() < p_nos)  { nn[b] |= 0b10; }  // set Z error
        #endif
    }


    //-- symplectic syndrome measurement
    Quan_GenSyndrome(A, nn, zz);

    #if MODE >= 2
    p_bias = p_nos; // do necessary scaling if needed

    #if SORTCOMPARE == LASTFORSORT && MODE >= 3
        if(osdw != -1){
            for ( b = 0 ; b < N ; b ++ ) {
                usedTo[b] = 0;
                lastFor[b] = 0;
            }
        }

    #endif // SORTCOMPARE


    afp =  (alpha==0 || alpha == 100) ?  0 :(double) 100.0/alpha;

    init_syndro_ok = Qbp_init20(bp, zz, A, p_bias, afp);  // contains a subroutine GenInitBiasVec

    syndrome_ok = init_syndro_ok;   iter = 0;

    while ( !init_syndro_ok ) {     iter ++;
        syndrome_ok = Qbp_dec20(bp, A);     ///////////////////////////////////////////////////

        #if MODE >= 3
        if(osdw != -1 && !syndrome_ok){
            for(b = 0; b<N; b++){
                // ==== calculate lastFor =====
                if(bp->tt[b] != usedTo[b]){
                    usedTo[b] = bp->tt[b];
                    lastFor[b] = iter;
                }
            }
        }
        #endif // SORTCOMPARE
          if ( syndrome_ok || iter==max_iter )   break;
    }
    //###################################################################################################

    # if MODE == 4
    if(osdw == -1 || syndrome_ok){  // no OSD or syndrome is OK
        finalResult = bp->tt;
    }
    else{  // Do OSDw or OSDfull
        // ==== calculate probability ===
        //printf("%d ", nCW);
        for(b = 0; b<N; b++){
            lastFor[b] = max_iter_1 - lastFor[b];
            sum=0; for(q=0;q<Q;q++) { sum += bp->qn[b][q]; }
            for ( q = 0 ; q < Q ; q ++ ) {
                prob[b][q] = bp->qn[b][q]/sum;
            }
        }
        if(osdw == -2){
            finalResult = post_decodeOSDfull(osdDec, zz, (const double **) prob, lastFor, max_iter);
        }
        else{
            finalResult = post_decodeOSDw(osdDec, zz, (const double **) prob, lastFor);
        }

        # if OSD_RECORD
            if (osdw == -2){
                osdRec -> nStep1 ++;
                #if STORE_DISTR
                if (osdDec -> status == UsePosd){
                    osdRec -> newMatSizeDistr[osdDec -> numOfUnreliCheck][osdDec -> numOfUnreliBit] ++;
                    if (osdDec -> osd0orNot){
                        osdRec -> newMatSizeDistr_osd0[osdDec -> numOfUnreliCheck][osdDec -> numOfUnreliBit] ++;
                    }
                }else{
                    osdRec -> nStep2 ++;
                    osdRec -> newMatSizeDistr[M][N<<1] ++;
                }
                #else
                if (osdDec -> status != UsePosd) osdRec -> nStep2 ++;
                #endif // STORE_DISTR
            }
        #endif // OSD_RECORD

        # if SORTCOMPARE == LASTFORSORT && STORE_LF && OSDW !=-1
        for(int i = 0; i<N; i++){
            storeLF[lastFor[i]] ++;
        }
        #endif
    }
    #endif // MODE

    nIterAcc += iter;
    #endif // MODE

    nCW++;

    if(nCW % (NUM_RECORD) == 0) {
        tEnd = clock();  currSec = tEnd-tLast; sec = tEnd-tStr;
        printf("nCW: %d, current sec: %d, total sec: %d\n", nCW, currSec, sec);
        tLast = clock();

        aiter=(double)nIterAcc/nCW;

        FILE *fpSt = fopen ( path_St , "w" );   //= fopen ( path_Last , "a+w" );
        # if OSD_RECORD
        printf("P_nos: %g\t AvgIter: %g\t nCW: %llu\tsec: %u\tnStep1: %llu\tnStep2: %llu\t\n",  p_nos, aiter, nCW, sec, osdRec -> nStep1, osdRec -> nStep2);

        fprintf(fpSt,"%g\t%g\t%llu\t%u\t%llu\t%llu\t\n",  p_nos, aiter, nCW, sec, osdRec -> nStep1, osdRec -> nStep2);
        #else
        printf("P_nos: %g\t AvgIter: %g\t nCW: %llu\tsec: %u\t\n",  p_nos, aiter, nCW, sec);
        fprintf(fpSt,"%g\t%g\t%llu\t%u\t\n",  p_nos, aiter, nCW, sec);
        #endif // OSD_RECORD

        #if CAL_TIME
        store_time[(nCW/NUM_RECORD)-1] = currSec;
        for(int i = 0; i<nCW/NUM_RECORD; i++){
            fprintf(fpSt, "%llu\n", store_time[i]);
        }
        #endif // CAL_TIME

        #if STORE_DISTR
        uint32_t NNbit = N<<1;
        uint32_t sum_check = 0;
        for(int i = 1; i<=M; i++){
            for(int j = 1; j<=NNbit; j++){
                fprintf(fpSt, "%10llu", osdRec -> newMatSizeDistr[i][j]);
                sum_check += osdRec -> newMatSizeDistr[i][j];
            }
            fprintf(fpSt, "\n");
        }
        printf("Check: nStep1 = %d, Correct :%d\n", sum_check, (sum_check == osdRec -> nStep1));
        sum_check = 0;
        for(int i = 1; i<=M; i++){
            for(int j = 1; j<=NNbit; j++){
                fprintf(fpSt, "%10llu", osdRec -> newMatSizeDistr_osd0[i][j]);
                sum_check += osdRec -> newMatSizeDistr_osd0[i][j];
            }
            fprintf(fpSt, "\n");
        }
        printf("Number of osd0: %d\n\n", sum_check);
        #endif // STORE_DISTR

        fclose(fpSt);
    }

  } while (nCW < NUM_DATA);

  fclose(fpA);
  fclose(fpGs);
	printf("main func done!\n");
  //free_OSD(osdDec);
//#if LOG_ZE_ai
//  fclose(fpZE);
//#endif
	getchar();
	return 0;
}
