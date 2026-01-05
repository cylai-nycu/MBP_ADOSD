#include "./ldpc_parm.h"
//#include "./lib_mat/lib_mat.h"
#include "./lib_rand/lib_rand.h"
#include "./lib_math/fast_math.h"
#include "./bp_dec/bp_dec.h"
#include "./bp_dec/bp_llr.h"

#include <math.h>
#include <time.h>
#include <inttypes.h>

//FILE *fpErr;

//-- these are only for status, not affecting dec logic
/*#define REC_N 10000
int32_t rec_cnt = 0;   // -2: OFF, -1:ON (and it will then self enable and add), 0:ENABLE
//uint8_t rec_tbl[REC_N][M];  // syndrome recording table (for those that cannot be decoded)
int32_t hit_cnt = 0;*/

/*
        MAIN
                     */
int main4(void)
{
#if 0 // NORMAL_TEST
  //#define N_ERR_STOP (20)  // needed error CWs (at beginning)
  #define N_ERR_STOP (100)  // needed error CWs (at beginning)
  #define N_ERR_ACCURATE  (100) // needed error CWs (when BLER curve starts to drop)
  #define N_ERR_STOP_1em1 (100)  // needed error CWs (after BLER <= 1e-1)
  #define N_ERR_STOP_1em2 (100)  // needed error CWs (after BLER <= 1e-2)
  #define N_ERR_STOP_1em3 (100)   // needed error CWs (after BLER <= 1e-3)

  #define N_ERR_STOP_1em4 (100)   // needed error CWs (after BLER <= 1e-4)
  //#define N_ERR_STOP_1em4 (20)   // needed error CWs (after BLER <= 1e-4)
  //#define N_ERR_STOP_1em4 (5)   // needed error CWs (after BLER <= 1e-4)

  #define N_ERR_STOP_1em5 (100)   // needed error CWs (after BLER <= 1e-5)
  //#define N_ERR_STOP_1em5 (20)   // needed error CWs (after BLER <= 1e-5)
  //#define N_ERR_STOP_1em5 (5)   // needed error CWs (after BLER <= 1e-5)

  #define N_ERR_STOP_1em6 (100)   // needed error CWs (after BLER <= 1e-6)
  //#define N_ERR_STOP_1em6 (20)   // needed error CWs (after BLER <= 1e-6)
  //#define N_ERR_STOP_1em6 (5)   // needed error CWs (after BLER <= 1e-6)

  #define N_ERR_STOP_1em7 (100)   // needed error CWs (after BLER <= 1e-7)
  //#define N_ERR_STOP_1em7 (20)   // needed error CWs (after BLER <= 1e-7)
  //#define N_ERR_STOP_1em7 (5)   // needed error CWs (after BLER <= 1e-7)
  //#define N_ERR_STOP_1em7 (1)   // needed error CWs (after BLER <= 1e-7)

  //#define BLER_STOP (8.0e-2)  // sim stop BLER
  //#define BER_STOP  (8.0e-2)  // sim stop BER
  ////#define BLER_STOP (1.0e-2)  // sim stop BLER
  ////#define BER_STOP  (1.0e-2)  // sim stop BER
  //#define BLER_STOP (4.0e-5)  // sim stop BLER
  //#define BER_STOP  (4.0e-5)  // sim stop BER
  #define BLER_STOP (1.0e-5)  // sim stop BLER
  #define BER_STOP  (1.0e-5)  // sim stop BER
  //#define BLER_STOP (4.0e-6)  // sim stop BLER
  //#define BER_STOP  (4.0e-6)  // sim stop BER
  //#define BLER_STOP (1.0e-6)  // sim stop BLER
  //#define BER_STOP  (1.0e-6)  // sim stop BER
  ////#define BLER_STOP (4.0e-7)  // sim stop BLER
  ////#define BER_STOP  (4.0e-7)  // sim stop BER
  //#define BLER_STOP (1.0e-7)  // sim stop BLER
  //#define BER_STOP  (1.0e-7)  // sim stop BER
  ////#define BLER_STOP (1.2e-10)  // sim stop BLER
  ////#define BER_STOP  (1.2e-10)  // sim stop BER

#elif 0 // QUICK_TEST
  /*//#define N_ERR_STOP (20)  // needed error CWs (at beginning)
  #define N_ERR_STOP (5)  // needed error CWs (at beginning)
  #define N_ERR_ACCURATE  (20) // needed error CWs (when BLER curve starts to drop)
  #define N_ERR_STOP_1em1 (20)  // needed error CWs (after BLER <= 1e-1)
  #define N_ERR_STOP_1em2 (20)  // needed error CWs (after BLER <= 1e-2)
  #define N_ERR_STOP_1em3 (20)   // needed error CWs (after BLER <= 1e-3)
  #define N_ERR_STOP_1em4 (5)   // needed error CWs (after BLER <= 1e-4)
  #define N_ERR_STOP_1em5 (5)   // needed error CWs (after BLER <= 1e-5)
  #define N_ERR_STOP_1em6 (3)   // needed error CWs (after BLER <= 1e-6)
  #define N_ERR_STOP_1em7 (1)   // needed error CWs (after BLER <= 1e-7)*/
  #define N_ERR_STOP (5)    // needed error CWs (at beginning)
  #define N_ERR_ACCURATE  (5) // needed error CWs (when BLER curve starts to drop)
  #define N_ERR_STOP_1em1 (5)  // needed error CWs (after BLER <= 1e-1)
  #define N_ERR_STOP_1em2 (5)  // needed error CWs (after BLER <= 1e-2)
  #define N_ERR_STOP_1em3 (5)   // needed error CWs (after BLER <= 1e-3)
  #define N_ERR_STOP_1em4 (3)   // needed error CWs (after BLER <= 1e-4)
  #define N_ERR_STOP_1em5 (3)   // needed error CWs (after BLER <= 1e-5)
  #define N_ERR_STOP_1em6 (3)   // needed error CWs (after BLER <= 1e-6)
  //#define N_ERR_STOP_1em6 (1)   // needed error CWs (after BLER <= 1e-6)
  #define N_ERR_STOP_1em7 (1)   // needed error CWs (after BLER <= 1e-7)
  //#define BLER_STOP (4.0e-5)  // sim stop BLER
  //#define BER_STOP  (4.0e-5)  // sim stop BER
  #define BLER_STOP (4.0e-5)  // sim stop BLER
  #define BER_STOP  (4.0e-5)  // sim stop BER
  //#define BLER_STOP (4.0e-6)  // sim stop BLER
  //#define BER_STOP  (4.0e-6)  // sim stop BER
  ////#define BLER_STOP (2.0e-6)  // sim stop BLER
  ////#define BER_STOP  (2.0e-6)  // sim stop BER
  //#define BLER_STOP (1.0e-6)  // sim stop BLER
  //#define BER_STOP  (1.0e-6)  // sim stop BER
  ////#define BLER_STOP (1.0e-7)  // sim stop BLER
  ////#define BER_STOP  (1.0e-7)  // sim stop BER
#else
  #define N_ERR_STOP (500)    // needed error CWs (at beginning)
  #define N_ERR_ACCURATE  (500) // needed error CWs (when BLER curve starts to drop)
  #define N_ERR_STOP_1em1 (500)  // needed error CWs (after BLER <= 1e-1)
  #define N_ERR_STOP_1em2 (500)  // needed error CWs (after BLER <= 1e-2)
  #define N_ERR_STOP_1em3 (500)   // needed error CWs (after BLER <= 1e-3)
  #define N_ERR_STOP_1em4 (500)   // needed error CWs (after BLER <= 1e-4)
  #define N_ERR_STOP_1em5 (500)   // needed error CWs (after BLER <= 1e-5)
  #define N_ERR_STOP_1em6 (500)   // needed error CWs (after BLER <= 1e-6)
  #define N_ERR_STOP_1em7 (500)   // needed error CWs (after BLER <= 1e-7)
  //#define BLER_STOP (4.0e-5)  // sim stop BLER
  //#define BER_STOP  (4.0e-5)  // sim stop BER
  #define BLER_STOP (4.0e-2)  // sim stop BLER
  #define BER_STOP  (4.0e-2)  // sim stop BER
#endif


//#define P_START  P_STEP
//#define P_START (-0.0)   // start sim err probability (in dB, e.g. -0.0 means from 10^(0/10) = 1)
//#define P_START (-1.25)   // start sim err probability (in dB, e.g. -1.25 means from 10^(-1.25/10) = 0.75)
//#define P_START (-5.0)    // start sim err probability (in dB, e.g. -5 means from 10^(-5/10) = 0.3162)
//#define P_START (-5.625)  // start sim err probability (in dB, e.g. -5 means from 10^(-5.625/10) = 0.27384)
//#define P_START (-10.0)   // start sim err probability (in dB, e.g. -10 means from 10^(-10/10) = 0.1)
#define P_START (-10.0  -  10.0/8 *4)     // 0.0316
////#define P_START (-10.0  -  10.0/8 *10)    // 0.0056
//#define P_START (-20.0)    // start sim err probability (in dB, e.g. -20 means from 10^(-20/10) = 0.01)
//#define P_START (-25.0)
//#define P_START (-45.0)
/*#define P_STEP  (-10.0/8)  // err probability step (in dB, every times 10^(STEP/10))  // (moved to ldpc_parm.h)*/



  int32_t iter, b;
  GFQ_t nn[N], diff[N], zz[M], zz_G[M_G];
  uint64_t tx_seed[4];

  //- RX
  double p_nos, rnd_val=0, p_bias, afp=0;                   if(rnd_val) { } // just to prevent compile warning when CH_TYPE not 0
  FILE *fpA = fopen ( PATH_A , "r" );
  FILE *fpGs = fopen( PATH_Gs, "r" );


  a_matrix_GFQ  Amtx, *A=&Amtx;
  g_matrix_GFQ  Gmtx, *G=&Gmtx;
  uint8_t syndrome_ok, init_syndro_ok;

  //- init proc
  //fast_math_test();
  rnd256_init();    rnd256_init_priv(tx_seed);
  //- load sparse parity-check matrix from file to a_matrix_GFQ *A
  load_A_GFQ(fpA, A);
  load_G_GFQ(fpGs, G);

  #if LLR_BP == 0
  QBP_Ctl       BPC , *bp=&BPC ;
  alloc_QBPC(A, bp);   // alloc "Quantum" BP decoder needed resource
  #else
  LBP_Ctl       BPC , *bp=&BPC ;
  alloc_LBPC(A, bp);   // alloc LLR-based Quantum BP decoder resource
  #endif

#if USE_GF2_DEC == 0
#elif USE_GF2_DEC == 1
  double bias2[2*N][2];   //uint8_t nn2[2*N];   //uint8_t target_zz2[M];
  FILE *fpA_GF2 = fopen ( PATH_A_GF2 , "r" );
  a2_matrix Amtx_GF2, *A_GF2=&Amtx_GF2;     load_A2(fpA_GF2, A_GF2);
  BP2_Ctl   BPC_GF2 , *bp_GF2=&BPC_GF2;     alloc_BPC2(A_GF2, bp_GF2);
  //-- check afp to prevent compiler warning
  if(afp==0) { } //{ printf("sizeof(bias2) = %d", sizeof(bias2));   getchar(); }
#else
  XXX_NG_GFQ_DEC_SEL_XXX
#endif


  int32_t n_list_sol = 0;  if(n_list_sol){ }    // number of list solutions (and prevent compile warning)
#if LIST_DEC == 0
#elif LIST_DEC == 1
  #if USE_GF2_DEC == 0
  GFQ_t ListSols[LIST_DEC][N];
  #elif USE_GF2_DEC == 1
  GFQ_t ListSols_GF2[LIST_DEC][2*N];
  #endif
#else
  XXX_NG_LIST_DEC_XXX
#endif


  char st_name[1000] = "";
  char path_St[1000] ;  //= "St";     //sprintf(path_St,  "St_%u_%u_it%u_dec%u_1e-11%s%s%s", N,K,MAX_ITER,BY_DEC,(FIX_P? "_fixP":""),(ALFA? "_sq2":""),(BETA? "_ofs":""));
  char path_Last[1000]; //="StLast"; //sprintf(path_Last,"St_%u_%u_it%u_dec%u_1e-11%s%s%s_Last", N,K,MAX_ITER,BY_DEC,(FIX_P? "_fixP":""),(ALFA? "_sq2":""),(BETA? "_ofs":""));
  char path_Scr[1000] ; //= "Scr";   //sprintf(path_Scr,"Scr_%u_%u_it%u_dec%u_1e-11%s%s%s", N,K,MAX_ITER,BY_DEC,(FIX_P? "_fixP":""),(ALFA? "_sq2":""),(BETA? "_ofs":""));

  if(CH_TYPE==0)      { } // dep ch, no special mark
  else if(CH_TYPE==1) { strcpy(st_name, "XZ_"); }  // indep X-Z ch
  else                { printf("CH_TYPE %d undefined !!  now exit();", CH_TYPE);  getchar();  exit(1); }

  sprintf(st_name, "%s%s_it%u_dec%u", st_name, PRE_A, MAX_ITER, BY_DEC);  // care: st_name in both I/O , behavior not guaranteed
  if(USE_GF2_DEC) sprintf(st_name, "%s_%s", st_name, "GF2");        // care: st_name in both I/O , behavior not guaranteed
  if(RND_SCHE)    sprintf(st_name, "%s_%s%u", st_name, "rndS", RND_SCHE); // care: st_name in both I/O , behavior not guaranteed
    // also cfg Qbp_init20(): rnd means random order (Fisher-Yates approach) , rndS means Sattolo algorithm

  if(LIST_DEC)    sprintf(st_name, "%s_%s%u", st_name, "list", LIST_DEC); // care: st_name in both I/O , behavior not guaranteed
  if(ZE)      sprintf(st_name, "%s_%s%u", st_name, "ze", ZE);       // care: st_name in both I/O , behavior not guaranteed
  if(FIX_P)   sprintf(st_name, "%s_%s%u", st_name, "fixP", FIX_P);  // care: st_name in both I/O , behavior not guaranteed
  if(FIX_e)   sprintf(st_name, "%s%s%d",  st_name, "e",    FIX_e);  // care: st_name in both I/O , behavior not guaranteed
  if(AFN)     sprintf(st_name, "%s_%s%u", st_name, "an", AFN);      // care: st_name in both I/O , behavior not guaranteed
  if(ALFA)    sprintf(st_name, "%s_%s%u", st_name, "ac", ALFA);     // care: st_name in both I/O , behavior not guaranteed

  if(AFP2){  //ap: DEC2=DEC1, app: force DEC2=20, aps: force DEC2=24, apc: force DEC2=44
    if(DEC2=='0') sprintf(st_name, "%s_%s",   st_name, "ap");       // care: st_name in both I/O , behavior not guaranteed
    else          sprintf(st_name, "%s_%s%c", st_name, "ap",DEC2);  // care: st_name in both I/O , behavior not guaranteed

    sprintf(st_name, "%s%u", st_name, AFP);     // care: st_name in both I/O , behavior not guaranteed
    if(AFPINC)  sprintf(st_name, "%s,%d", st_name, AFPINC);       // care: st_name in both I/O , behavior not guaranteed
    sprintf(st_name, "%s,%u", st_name, AFP2);   // care: st_name in both I/O , behavior not guaranteed

    if(DBL_CHK) sprintf(st_name, "%s_%s%u", st_name, "dbl", DBL_CHK); // care: st_name in both I/O , behavior not guaranteed

  }else if(AFP) sprintf(st_name, "%s_%s%u", st_name, "ap", AFP);      // care: st_name in both I/O , behavior not guaranteed

  if(GMN)     sprintf(st_name, "%s_%s%u-%u", st_name, "gmn", GMN,LMN);     // care: st_name in both I/O , behavior not guaranteed
  if(CYC_CHK) sprintf(st_name, "%s_%s%u", st_name, "cyc", CYC_CHK); // care: st_name in both I/O , behavior not guaranteed


  if(AFPADD)  sprintf(st_name, "%s_%s%d", st_name, "ad", AFPADD);   // care: st_name in both I/O , behavior not guaranteed
  if(BETA)    sprintf(st_name, "%s_%s%u", st_name, "b", BETA);      // care: st_name in both I/O , behavior not guaranteed
#if EN_AFP_V==0
  if(ALFA_V)  sprintf(st_name, "%s_%s%u", st_name, "av", ALFA_V);   // care: st_name in both I/O , behavior not guaranteed
#else
  if(ALFA_V)  sprintf(st_name, "%s_%s%u", st_name, "avp", ALFA_V);   // care: st_name in both I/O , behavior not guaranteed
#endif
  if(BETA_pn) sprintf(st_name, "%s_%s%u", st_name, "beta", BETA_pn); // care: st_name in both I/O , behavior not guaranteed
  if(POS_AVG) sprintf(st_name, "%s_%s%u", st_name, "pa", POS_AVG);  // care: st_name in both I/O , behavior not guaranteed
  if(STOP_BY_C) sprintf(st_name, "%s_%s", st_name, "stopC");  // care: st_name in both I/O , behavior not guaranteed
  if(MAJ_VOTE)  sprintf(st_name, "%s_%s%u", st_name, "maj", MAJ_VOTE);  // care: st_name in both I/O , behavior not guaranteed
  if(MAJ_SCHE)  sprintf(st_name, "%s_%s%u", st_name, "msch", MAJ_SCHE); // care: st_name in both I/O , behavior not guaranteed

  sprintf(path_St,  "St_%s", st_name);      sprintf(path_Scr,  "St_%s.Scrn", st_name);      sprintf(path_Last,  "St_%s.Last", st_name);
  FILE *fpSt = fopen ( path_St , "a+w" );   FILE *fpScr = fopen ( path_Scr , "a+w" );   FILE *fpLast; //= fopen ( path_Last , "a+w" );

//======================================
FILE *fpRecordcorrect = fopen("../Data/record_correct", "w");
FILE *fpRecordWrong = fopen("../Data/record_wrong", "w");
FILE *fpRecordLong = fopen("../Data/record_Long", "w");
int correctCount = 0;

                // Fa=False, Eb=err bits                    // Q=quantum, Eq=err qubits
  uint64_t nCW, nErr, nFa, nEb, nIterAcc, nStop=N_ERR_STOP,  nErrQ, nFaQ, nEq;
  uint32_t nErrBit,  nErrQbit;
  double   bler, ber, far,  qbler, qber, qfar,  chk_bler, chk_ber,  aiter;
                    //far=false rate                              //aiter=avg iter

  time_t  tStr, tEnd, sec, sec_last=0;  // care that: gcc32 has year 2038 problem (gcc64 no this problem)

  nCW=0; nErr=0; nFa=0; nEb=0; nIterAcc=0; nErrBit=0;  nErrQ=0; nFaQ=0; nEq=0; nErrQbit=0;
  p_nos = pow(10, P_START/10);  chk_bler=1.0; chk_ber=1.0;
  printf("//==== init p_nos = %g ====// N = %u K = %u   It = %u  BY_DEC %u \n", p_nos,N,K,MAX_ITER,BY_DEC);
  fprintf(fpScr, "//==== init p_nos = %g ====// N = %u K = %u   It = %u  BY_DEC %u \n", p_nos,N,K,MAX_ITER,BY_DEC);
  fclose(fpScr);
  tStr = time(NULL);

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


  #if FIX_P == 0
    p_bias = p_nos*1.0; // do necessary scaling if needed
    //p_bias = p_nos*40.0;  if(p_bias>4.6e-2) p_bias=4.6e-2;
    //p_bias = p_nos*10.0; // do necessary scaling if needed
  #else   // fix p_bias
    //p_bias = (double)FIX_P/10 * 1e-2;
    p_bias = (double)FIX_P/10 * 1e-2 * (FIX_e? pow(10.0, FIX_e+3) : 1.0);  // when FIX_e=-2, 1/10*1e-2*10^(FIX_e+3) = 10^-2
    //p_bias = 4.6e-2;   // fix p_bias
  #endif

  #if ZE
    double p_est = (double) HamWt(zz,M)*ZE/100 / N;
    if(p_est > p_nos) p_bias = p_est;
    else              p_bias = p_nos;
  #endif


  #if USE_GF2_DEC //########################################################################################
   #if CH_TYPE==0    // dep ch
    double p_bsc = p_bias*2/3;
   #elif CH_TYPE==1  // indep X-Z ch
    double p_bsc = p_bias;
   #endif

    /*for(b=0; b<N ; b++) { nn2[b] = nn[b]&0b01;   nn2[b+N] = (nn[b]&0b10)>>1; }
    for(b=0; b<2*N;b++) if(nn2[b]){ bias2[b][0] = p_bsc;   bias2[b][1] = (1-p_bsc); }
                        else      { bias2[b][0] = (1-p_bsc);   bias2[b][1] = p_bsc; }
    bp_GF2->p_ch = p_bsc;
    init_syndro_ok = bp2_init20(bp_GF2, bias2, A_GF2, 's');   // before return, will set all bp->tt[n]=0 if mode=='s'  */

    for(b=0; b<2*N; b++) { bias2[b][0] = (1-p_bsc);   bias2[b][1] = p_bsc; }   // like if all nn2[b]=0
    //GenSyndrome2(A_GF2, nn2, target_zz2);  // generate target syndrome from nn2  (SKIP, use zz as the same)
    init_syndro_ok = bp2_init22(bp_GF2, bias2, A_GF2, zz);

    n_list_sol = 0;
    syndrome_ok = init_syndro_ok;   iter = 0;
    //-- iterative BP decoding
    while ( !init_syndro_ok ) {     iter ++;
    #if BY_DEC == 20
      syndrome_ok = bp2_dec20(bp_GF2, A_GF2);
    #elif BY_DEC == 24
      syndrome_ok = bp2_dec24(bp_GF2, A_GF2);
    #endif

    #if LIST_DEC == 0   // DEFAULT_MODE
      if ( syndrome_ok || iter==MAX_ITER )   break;
    #elif LIST_DEC == 1 // BUF_LAST
      if(syndrome_ok) { n_list_sol = 1;   memcpy(ListSols_GF2[0], bp_GF2->tt, sizeof(ListSols_GF2[0])); }  // buffer the last solution
      if(iter==MAX_ITER)  break;
      //if(n_list_sol > 0)  break;  // enable this row if do TEST to match the original behavior (when LIST_DEC=0)
    #endif

    }
    nIterAcc += iter;
    //printf("BP_dec syndrome %s  iter %d  TX-RX vec HamDist %d \n", syndrome_ok? "OK!":"FAIL!!!" , iter , HamDist(nn2, bp_GF2->tt, 2*N));

    #if LIST_DEC == 1   // BUF_LAST
    if(n_list_sol > 0) { memcpy(bp_GF2->tt, ListSols_GF2[0], sizeof(bp_GF2->tt)); }
    #endif

    for(b=0; b<N ; b++) { bp->tt[b] = bp_GF2->tt[b] + bp_GF2->tt[b+N]*2; }    // pass result to bp->tt for statistics

  #else // USE_GF2_DEC==0   // (default mode, decoding by GF4) //###########################################################

    //%%%%%%%%%%%% DEFAULT_RUN that could be with AFP (without AFP2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    afp = (AFP==0)?  0 : (double)100.0/AFP;

    //-- "Quantum" BP init: set bp->tt as all zeros. return 1 if zz pass; otherwise INIT *bp and return 0
    #if LLR_BP == 0
    init_syndro_ok = Qbp_init20(bp, zz, A, p_bias, afp);  // contains a subroutine GenInitBiasVec
    //init_syndro_ok = Qbp_init31(bp, zz, A, p_bias, afp);   /////////////////////////////////////
    #else
    init_syndro_ok = Lbp_init(bp, zz, A, p_bias, afp);  // contains a subroutine GenInitLLRxyz
    #endif // LLR_BP

    n_list_sol = 0;
    syndrome_ok = init_syndro_ok;   iter = 0;
    //-- iterative BP decoding
    while ( !init_syndro_ok ) {     iter ++;
    #if BY_DEC==20
      syndrome_ok = Qbp_dec20(bp, A);     ///////////////////////////////////////////////////
    #elif BY_DEC==24 || BY_DEC==25 //|| BY_DEC==26
      syndrome_ok = Qbp_dec24(bp, A);     ///////////////////////////////////////////////////
    #elif BY_DEC==44
      syndrome_ok = Qbp_dec44(bp, A);     ///////////////////////////////////////////////////
    #elif BY_DEC==60
      syndrome_ok = Lbp_dec60(bp, A);     ///////////////////////////////////////////////////
    #elif BY_DEC==64
      syndrome_ok = Lbp_dec64(bp, A);     ///////////////////////////////////////////////////
    #elif BY_DEC==84
      syndrome_ok = Lbp_dec84(bp, A);     ///////////////////////////////////////////////////
    #endif

      /*for ( b = 0 ; b < N ; b ++ ) { printf("%u ", nn[b]); }      printf(" = nn \n\n");
      for ( b = 0 ; b < N ; b ++ ) { printf("%u ", bp->tt[b]); }  printf(" = bp->tt \n\n");
      nErrBit = HamDist(nn, bp->tt, N);
      printf("nCW %"PRIu64"  iter %u  syndrome_ok %u   nErrBit %u \n", nCW,iter,syndrome_ok,nErrBit);
      getchar();*/

    #if LIST_DEC == 0   // DEFAULT_MODE
      if ( syndrome_ok || iter==MAX_ITER )   break;
    #elif LIST_DEC == 1 // BUF_LAST
      if(syndrome_ok) { n_list_sol = 1;   memcpy(ListSols, bp->tt, sizeof(ListSols[0])); }  // buffer the last solution
      if(iter==MAX_ITER)  break;
      //if(n_list_sol > 0)  break;  // enable this row if do TEST to match the original behavior (when LIST_DEC=0)
    #endif
    }

    nIterAcc += iter;

    #if LIST_DEC == 1   // BUF_LAST
    if(n_list_sol > 0) { memcpy(bp->tt, ListSols[0], sizeof(bp->tt)); }
    #endif


    //%%%%%%%%%%%% AFP2_RUN that is additional for AFP2!=0. Run this only when the previous run's syndrome is NG  %%%%%%%%%%%%%%%//
    #if AFP2
    int32_t  trial = 0;
      #if DBL_CHK
      memcpy(bp->tt_last, bp->tt, sizeof(bp->tt_last));
      init_syndro_ok = 0;  // force to find at least 2 times
      #endif
    while(!init_syndro_ok) // AFP2_main_while
    {
      if(AFPINC==0)   afp = (AFP2==0)?  0 : (double)100.0/AFP2;
      else { trial++; afp = (double)100.0 / (AFP + AFPINC*trial); }

      init_syndro_ok = Qbp_init20(bp, zz, A, p_bias, afp);  // contains a subroutine GenInitBiasVec

      n_list_sol = 0;
      syndrome_ok = init_syndro_ok;   iter = 0;
      while ( !init_syndro_ok ) {     iter ++;
      //-- ap
      #if DEC2 == '0'
        #if BY_DEC==20
        syndrome_ok = Qbp_dec20(bp, A);
        #elif BY_DEC==24 || BY_DEC==25 //|| BY_DEC==26
        syndrome_ok = Qbp_dec24(bp, A);
        #elif BY_DEC==44
        syndrome_ok = Qbp_dec44(bp, A);
        #else
        XXX_BY_QBP_DEC_XXX
        #endif
      //-- ap + force schedule
      #elif DEC2 == 'p'
        syndrome_ok = Qbp_dec20(bp, A);   //app
      #elif DEC2 == 's'
        syndrome_ok = Qbp_dec24(bp, A);   //aps
      #elif DEC2 == 'c'
        syndrome_ok = Qbp_dec44(bp, A);   //apc
      #else
        XXX_NG_DEC2_XXX
      #endif


      #if LIST_DEC == 0   // DEFAULT_MODE
        if ( syndrome_ok || iter==MAX_ITER )   break;
      #elif LIST_DEC == 1 // BUF_LAST
        if(syndrome_ok) { n_list_sol = 1;   memcpy(ListSols, bp->tt, sizeof(ListSols[0])); }  // buffer the last solution
        if(iter==MAX_ITER)  break;
        //if(n_list_sol > 0)  break;  // enable this row if do TEST to match the original behavior (when LIST_DEC=0)
      #endif
      }

      nIterAcc += iter;

      #if LIST_DEC == 1   // BUF_LAST
      if(n_list_sol > 0) { memcpy(bp->tt, ListSols[0], sizeof(bp->tt)); }
      #endif

      if(afp==(double)100.0/AFP2)   break;

      init_syndro_ok = syndrome_ok; // if this decoding instance syndrome_ok, we say next instance init_syndro_ok (will not run next instance)
      #if DBL_CHK == 1
      //if(init_syndro_ok)
        if(HamDist(bp->tt, bp->tt_last, N) != 0)
          init_syndro_ok = 0;    // force syndrome NG to continue AFP2_main_while
      memcpy(bp->tt_last, bp->tt, sizeof(bp->tt_last));
      #elif DBL_CHK == 2
      //if(afp==(double)100.0/AFP)                      init_syndro_ok = 0;  // find at least 2 times
      //else if(HamDist(bp->tt, bp->tt_last, N) == 0)   init_syndro_ok = 0;  // find until outputs different (for 2 times)
      if(HamDist(bp->tt, bp->tt_last, N) == 0)   init_syndro_ok = 0;   // force find until outputs different (for 2 times)
      memcpy(bp->tt_last, bp->tt, sizeof(bp->tt_last));
      #elif DBL_CHK == 3
      //if(afp==(double)100.0/AFP)  init_syndro_ok = 0;  // find at least 2 times
      //else if(iter < MAX_ITER)    init_syndro_ok = 0;  // find next if iter max not hit
      if(iter < MAX_ITER)    init_syndro_ok = 0;   // force find next if iter max not hit
      else{ // iter max is hit: solution is by the last time
        memcpy(bp->tt, bp->tt_last, sizeof(bp->tt));   break;
      }
      memcpy(bp->tt_last, bp->tt, sizeof(bp->tt_last));
      #endif
    } // END of AFP2_main_while: while(!init_syndro_ok)...
    #endif  // END of "#if AFP2"
    //%%%%%%%%%%%% END OF additional run for AFP2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

  #endif  // END of USE_GF2_DEC==0   //########################################################################################


    //==== statistics
    nCW++;
    nErrBit = HamDist(nn, bp->tt, N);
    if(syndrome_ok) {
      if(nErrBit==0) { }  // fine only when both syndrome and HamDist pass
      else { nErr++; nEb+=nErrBit; nFa++; }
    } else { nErr++; nEb+=nErrBit;}

    //-- "Quantum" degeneracy check
    nErrQbit = nErrBit;   // start with classical error
    if(syndrome_ok && nErrQbit!=0) {
      VecDiff(diff, bp->tt, nn, N);     // if bp->tt is a degenerate err of nn,
      Quan_DegSyndrome(G, diff, zz_G);  // then diff will gen zz_G = all_zero
      //Quan_LogSyndrome(G, diff, zz_G);  // then diff will gen zz_G = all_zero (here only check last 2K logical rows)
      if(is_zero_vec(zz_G, M_G)) { nErrQbit = 0; } // force 0 if it's degenerate error
    }
    if(nErrQbit!=0) {
      nErrQ++;  nEq+=nErrQbit;
      if(syndrome_ok) nFaQ++;
//=========================================
for(int i = 0; i<N; i++)
       fprintf(fpRecordWrong, "%d ", nn[i]);
fprintf(fpRecordWrong, "\n");
    }
    else{
//=========================================
if(iter>=80){
    for(int i = 0; i<N; i++)
        fprintf(fpRecordLong, "%d ", nn[i]);
fprintf(fpRecordLong, "\n");
}
else if(correctCount<200){
     for(int i = 0; i<N; i++)
        fprintf(fpRecordcorrect, "%d ", nn[i]);
fprintf(fpRecordcorrect, "\n");
correctCount++;
}

    }
    //-- MISMATCH happen
    if(0){  // CODE DEFAULT
    //if(nErrQbit!=0) {   // FOR TEST
    //if(nErrQbit!=0 && HamWt(nn,N)>=7) {   // FOR TEST
    //if(nErrQbit!=0 && HamWt(nn,N)>=6 && HamWt(bp->tt,N)<=4) {   // FOR TEST
    //if(nErrQbit!=0 && iter==MAX_ITER) {   // FOR TEST
    //if(nErrBit!=nErrQbit) {   // FOR TEST
    //if(!syndrome_ok && (nErrBit==0 || nErrQbit==0 || nErrBit!=nErrQbit) ) { // FOR TEST2
      for ( b = 0 ; b < N ; b ++ ) { printf("%u ", b%10); }          printf(" = b \n\n");
      for ( b = 0 ; b < M ; b ++ ) { printf("%u ", zz[b]); }      printf(" = zz \n\n");
      for ( b = 0 ; b < N ; b ++ ) { printf("%u ", nn[b]); }      printf(" = nn \n\n");
      for ( b = 0 ; b < N ; b ++ ) { printf("%u ", bp->tt[b]); }  printf(" = bp->tt \n\n");
      for ( b = 0 ; b < N ; b ++ ) { printf("%u ", diff[b]); }    printf(" = diff \n\n");
      for ( b = 0 ; b < M_G ; b ++ ) { printf("%u ", zz_G[b]); }    printf(" = zz_G \n\n");
      for ( b = 0 ; b < N ; b ++ ) { printf("%d ", bp->seq[b]); }    printf(" = bp->seq \n\n");
      printf("syndrome_ok %d  nErrBit %u  nErrQbit %u  nCW %llu  iter %d  wt(nn) %u  is_hit_c %u\n", syndrome_ok, nErrBit, nErrQbit, nCW, iter, HamWt(nn,N), bp->is_hit_c);
      printf("If hit MISMATCH like syndrome NG but nErrQbit==0, then suggest CTRL-C to exit and debug; \n");   getchar(); // FOR TEST2
    }
    /*if(rec_cnt>=0) if(nErrQbit!=0 && HamWt(nn,N)<=8) {
      hit_cnt ++;
      uint8_t add_new = 1;
      for(b=0; b<rec_cnt; b++) { if(HamDist(zz, rec_tbl[b], M) == 0) add_new = 0; }
      if(add_new) { memcpy(rec_tbl[rec_cnt], zz, sizeof(zz));   rec_cnt++; }
      //printf("\n rec_cnt %d \n", rec_cnt);
    }*/

  #define LOG_ZE_ai 0   // only used with AFP2
  #if LOG_ZE_ai
    static FILE *fpZE;  static uint8_t ZE_1st_time = 1;  if(ZE_1st_time){ fpZE = fopen("_log_ZE", "w");  ZE_1st_time = 0; }
    if(nErrQbit==0) fprintf(fpZE, "%d\t%d\t%d\n", HamWt(zz,M), HamWt(nn,N), (AFP + AFPINC*trial));
  #endif // LOG_ZE_ai


    if(nCW % (1000*1) == 0) { tEnd = time(NULL);  sec = tEnd-tStr; }
    if(sec >= sec_last+10) {
      fpScr= fopen ( path_Scr , "a+w" );    fpLast = fopen ( path_Last , "a+w" );   //fpSt = fopen ( path_St , "a+w" );
      tEnd = time(NULL);  sec = tEnd-tStr;  sec_last=sec;

      bler=(double)nErr/nCW;    ber=(double)nEb/(nCW*N);  far=(double)nFa/nCW;      aiter=(double)nIterAcc/nCW;
      qbler=(double)nErrQ/nCW;  qber=(double)nEq/(nCW*N); qfar=(double)nFaQ/nCW;

      //printf(" now p %g , nCW %"PRIu64" nErr %"PRIu64" nEb %"PRIu64"  nErrQ %"PRIu64" nEq %"PRIu64" => qBLER %g qBER %g  BLER %g BER %g  AvgIter %g  @ %ld sec \n",
        //p_nos , nCW,nErr,nEb,nErrQ,nEq, qbler,qber,bler,ber,aiter , sec);
      //- skip BLER and BER when show
      printf(" now p %g , nCW %llu nErr %llu nEb %llu  nErrQ %llu nEq %llu => qBLER %g qBER %g  AvgIter %g  @ %ld sec\n",
        p_nos , nCW,nErr,nEb,nErrQ,nEq, qbler,qber, aiter , sec );
      //printf(" now p %g , nCW %"PRIu64" nErr %"PRIu64" nEb %"PRIu64"  nErrQ %"PRIu64" nEq %"PRIu64" => qBLER %g qBER %g  AvgIter %g  @ %ld sec  rec_cnt %d %d\n",
        //p_nos , nCW,nErr,nEb,nErrQ,nEq, qbler,qber, aiter , sec , rec_cnt,hit_cnt);
      fprintf(fpScr, " now p %g , nCW %llu nErr %llu nEb %llu  nErrQ %llu nEq %llu => qBLER %g qBER %g  BLER %g BER %g  AvgIter %g  @ %ld sec \n",
        p_nos , nCW,nErr,nEb,nErrQ,nEq, qbler,qber,bler,ber,aiter , sec);

      //- row data in fpSt: p_nos, bler, far, ber,  aiter,  nCW, nErr, nEb,  sec,  nEq, nErrQ, qber, qbler, qfar
      fprintf(fpLast,"%g\t%g\t%g\t%g\t%g\t%llu\t%llu\t%llu\t%ld\t%llu\t%llu\t%g\t%g\t%g\n",
        p_nos,bler,far,ber,aiter, nCW,nErr,nEb, sec, nEq,nErrQ,qber,qbler,qfar);

      if(nErrQ >= N_ERR_STOP_1em7){
        if(1 > qbler && qbler > 1e-1) nStop = N_ERR_ACCURATE;
        if(qbler <= 1e-1)  nStop = N_ERR_STOP_1em1;  // start to accelerate
        if(qbler <= 1e-2)  nStop = N_ERR_STOP_1em2;  // start to accelerate more ...
        if(qbler <= 1e-3)  nStop = N_ERR_STOP_1em3;
        if(qbler <= 1e-4)  nStop = N_ERR_STOP_1em4;
        if(qbler <= 1e-5)  nStop = N_ERR_STOP_1em5;
        if(qbler <= 1e-6)  nStop = N_ERR_STOP_1em6;
        if(qbler <= 1e-7)  nStop = N_ERR_STOP_1em7;
      }

      fclose(fpScr);    fclose(fpLast);   //fclose(fpSt);
    }

    if (nErrQ >= nStop ) {
//NEXT:
      fpScr= fopen ( path_Scr , "a+w" );    fpLast = fopen ( path_Last , "a+w" );   //fpSt = fopen ( path_St , "a+w" );
      tEnd = time(NULL);  sec = tEnd-tStr;

      bler=(double)nErr/nCW;    ber=(double)nEb/(nCW*N);  far=(double)nFa/nCW;      aiter=(double)nIterAcc/nCW;
      qbler=(double)nErrQ/nCW;  qber=(double)nEq/(nCW*N); qfar=(double)nFaQ/nCW;

      chk_bler = qbler;   chk_ber = qber;

      //printf("p_nos %g  qBLER %g qBER %g  BLER %g BER %g  AvgIter %g  for nCW %"PRIu64" nErr %"PRIu64" nEb %"PRIu64"  nErrQ %"PRIu64" nEq %"PRIu64"  @ %ld sec \n",
        //p_nos,qbler,qber,bler,ber,aiter , nCW,nErr,nEb,nErrQ,nEq , sec);
      //- skip BLER and BER when show
      printf("p_nos %g  qBLER %g qBER %g  AvgIter %g  for nCW %llu nErr %llu nEb %llu  nErrQ %llu nEq %llu  @ %ld sec \n",
        p_nos,qbler,qber,aiter , nCW,nErr,nEb,nErrQ,nEq , sec);
      fprintf(fpScr, "p_nos %g  qBLER %g qBER %g  BLER %g BER %g  AvgIter %g  for nCW %llu nErr %llu nEb %llu  nErrQ %llu nEq %llu  @ %ld sec \n",
        p_nos,qbler,qber,bler,ber,aiter , nCW,nErr,nEb,nErrQ,nEq , sec);

      //- row data in fpSt: p_nos, bler, far, ber,  aiter,  nCW, nErr, nEb,  sec,  nEq, nErrQ, qber, qbler, qfar
      fprintf(fpSt,  "%g\t%g\t%g\t%g\t%g\t%llu\t%llu\t%llu\t%ld\t%llu\t%llu\t%g\t%g\t%g\n",
        p_nos,bler,far,ber,aiter, nCW,nErr,nEb, sec, nEq,nErrQ,qber,qbler,qfar);
      fprintf(fpLast,"%g\t%g\t%g\t%g\t%g\t%llu\t%llu\t%llu\t%ld\t%llu\t%llu\t%g\t%g\t%g\n",
        p_nos,bler,far,ber,aiter, nCW,nErr,nEb, sec, nEq,nErrQ,qber,qbler,qfar);

      nCW=0; nErr=0; nFa=0; nEb=0; nIterAcc=0; nErrBit=0;  nErrQ=0; nFaQ=0; nEq=0; nErrQbit=0;
      p_nos *= pow(10, P_STEP/10);

      if(1 > qbler && qbler > 1e-1) nStop = N_ERR_ACCURATE;
      if(qbler <= 1e-1)  nStop = N_ERR_STOP_1em1;  // start to accelerate
      if(qbler <= 1e-2)  nStop = N_ERR_STOP_1em2;  // start to accelerate more ...
      if(qbler <= 1e-3)  nStop = N_ERR_STOP_1em3;
      if(qbler <= 1e-4)  nStop = N_ERR_STOP_1em4;
      if(qbler <= 1e-5)  nStop = N_ERR_STOP_1em5;
      if(qbler <= 1e-6)  nStop = N_ERR_STOP_1em6;
      if(qbler <= 1e-7)  nStop = N_ERR_STOP_1em7;

      fclose(fpScr);    fclose(fpLast);   //fclose(fpSt);

      //if(qbler <= 1e-3) { if(rec_cnt==-1)  rec_cnt=0; }  // enable record cnt
    }

  } while (chk_bler>BLER_STOP || chk_ber>BER_STOP);

  fclose(fpA);
  fclose(fpGs);
	printf("main func done!\n");
  fpScr= fopen ( path_Scr , "a+w" );
  fprintf(fpScr, "main func done!\n");
  fclose(fpSt); fclose(fpLast); fclose(fpScr);
//#if LOG_ZE_ai
//  fclose(fpZE);
//#endif
	getchar();
	return 0;
}
