#include "./ldpc_parm.h"
#include "./DecConfig.h"
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

//FILE *fpErr;

//-- these are only for status, not affecting dec logic
/*#define REC_N 10000
int32_t rec_cnt = 0;   // -2: OFF, -1:ON (and it will then self enable and add), 0:ENABLE
//uint8_t rec_tbl[REC_N][M];  // syndrome recording table (for those that cannot be decoded)
int32_t hit_cnt = 0;*/

/*
        MAIN
                     */
int main(void)
{



#if 1 // NORMAL_TEST
  #define HLER_tot (1000)  //# of logical error to be collected at High error rates
  #define MLER_tot (400)  //# of logical error to be collected at High error rates
  #define LLER_tot (100)  //# of logical error to be collected at Low error rates

  //#define N_ERR_STOP (20)  // needed error CWs (at beginning)
  //#define N_ERR_STOP (5000)  // needed error CWs (at beginning)
  //#define N_ERR_ACCURATE  (5000) // needed error CWs (when BLER curve starts to drop)
  #define N_ERR_STOP (200)  // needed error CWs (at beginning)
  #define N_ERR_ACCURATE  (200) // needed error CWs (when BLER curve starts to drop)

  #define N_ERR_STOP_1em1 (HLER_tot)  // needed error CWs (after BLER <= 1e-1)
  #define N_ERR_STOP_1em2 (HLER_tot)  // needed error CWs (after BLER <= 1e-2)
  #define N_ERR_STOP_1em3 (HLER_tot)   // needed error CWs (after BLER <= 1e-3)

  #define N_ERR_STOP_1em4 (LLER_tot)   // needed error CWs (after BLER <= 1e-4)
  //#define N_ERR_STOP_1em4 (20)   // needed error CWs (after BLER <= 1e-4)
  //#define N_ERR_STOP_1em4 (5)   // needed error CWs (after BLER <= 1e-4)

  #define N_ERR_STOP_1em5 (LLER_tot)   // needed error CWs (after BLER <= 1e-5)
  //#define N_ERR_STOP_1em5 (20)   // needed error CWs (after BLER <= 1e-5)
  //#define N_ERR_STOP_1em5 (5)   // needed error CWs (after BLER <= 1e-5)

  #define N_ERR_STOP_1em6 (LLER_tot)   // needed error CWs (after BLER <= 1e-6)
  //#define N_ERR_STOP_1em6 (20)   // needed error CWs (after BLER <= 1e-6)
  //#define N_ERR_STOP_1em6 (5)   // needed error CWs (after BLER <= 1e-6)

  #define N_ERR_STOP_1em7 (LLER_tot)   // needed error CWs (after BLER <= 1e-7)
  //#define N_ERR_STOP_1em7 (20)   // needed error CWs (after BLER <= 1e-7)
  //#define N_ERR_STOP_1em7 (5)   // needed error CWs (after BLER <= 1e-7)
  //#define N_ERR_STOP_1em7 (1)   // needed error CWs (after BLER <= 1e-7)

  //#define BLER_STOP (0.1)  // sim stop BLER
  //#define BER_STOP  (0.1)  // sim stop BER
  //#define BLER_STOP (8.0e-2)  // sim stop BLER
  //#define BER_STOP  (8.0e-2)  // sim stop BER
  //#define BLER_STOP (1.0e-3)  // sim stop BLER
  //#define BER_STOP  (1.0e-3)  // sim stop BER
  //#define BLER_STOP (4.0e-5)  // sim stop BLER
  //#define BER_STOP  (4.0e-5)  // sim stop BER
  #define BLER_STOP (8.0e-6)  // sim stop BLER
  #define BER_STOP  (8.0e-6)  // sim stop BER
  //#define BLER_STOP (2.0e-6)  // sim stop BLER
  //#define BER_STOP  (2.0e-6)  // sim stop BER
  //#define BLER_STOP (1.0e-6)  // sim stop BLER
  //#define BER_STOP  (1.0e-6)  // sim stop BER
  //#define BLER_STOP (8.0e-7)  // sim stop BLER
  //#define BER_STOP  (8.0e-7)  // sim stop BER
  //#define BLER_STOP (1.0e-7)  // sim stop BLER
  //#define BER_STOP  (1.0e-7)  // sim stop BER
  //#define BLER_STOP (8.0e-8)  // sim stop BLER
  //#define BER_STOP  (8.0e-8)  // sim stop BER
  //#define BLER_STOP (1.2e-8)  // sim stop BLER
  //#define BER_STOP  (1.2e-8)  // sim stop BER

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
  #define N_ERR_STOP (10000)  // needed error CWs (at beginning)
  #define N_ERR_ACCURATE  (10000) // needed error CWs (when BLER curve starts to drop)
  #define N_ERR_STOP_1em1 (10000)  // needed error CWs (after BLER <= 1e-1)
  #define N_ERR_STOP_1em2 (10000)  // needed error CWs (after BLER <= 1e-2)
  #define N_ERR_STOP_1em3 (10000)   // needed error CWs (after BLER <= 1e-3)
  #define N_ERR_STOP_1em4 (10000)   // needed error CWs (after BLER <= 1e-4)
  #define N_ERR_STOP_1em5 (10000)   // needed error CWs (after BLER <= 1e-5)
  #define N_ERR_STOP_1em6 (10000)   // needed error CWs (after BLER <= 1e-6)
  #define N_ERR_STOP_1em7 (10000)   // needed error CWs (after BLER <= 1e-7)
  #define BLER_STOP (2.0e-2)  // sim stop BLER
  #define BER_STOP  (2.0e-2)  // sim stop BER
#endif

//#define P_START P_STEP*8// for normal set p_step 8
//#define P_START P_STEP*16// for zoom, set p_step 32
//#define P_START  P_STEP
//#define P_START (-0.0)   // start sim err probability (in dB, e.g. -0.0 means from 10^(0/10) = 1)
//#define P_START (-1.25)   // start sim err probability (in dB, e.g. -1.25 means from 10^(-1.25/10) = 0.75)
//#define P_START (-5.0)    // start sim err probability (in dB, e.g. -5 means from 10^(-5/10) = 0.3162)
//#define P_START (-7.85)  // start sim err probability (in dB, e.g. -5 means from 10^(-6.7/10)=0.21379620895)
//#define P_START (-6.98970004)  // start sim err probability (in dB, e.g. 10^(-6.98970004/10)=0.2000000)


// P_start= 10*log10(target initial error rate)
//#define P_START (-10.0)   // start sim err probability (in dB, e.g. -10 means from 10^(-10/10) = 0.1)
//#define P_START (-10.0  -  10.0/8 *2)     // 0.0562
////#define P_START (-10.0  -  10.0/8 *10)    // 0.0056
//#define P_START (-14.202)     // 0.038
//#define P_START (-10.0  -  10.0/8 *4) //-15.0    // 0.0316
//#define P_START (-15.739)     // 0.0266747
//#define P_START (-16.25)     // 0.0237137
#define P_START (-16.9897)     // 0.02
//#define P_START (-17.5)     // 0.0177828
//#define P_START (-18.0)     //
//#define P_START (-20.0)    // start sim err probability (in dB, e.g. -20 means from 10^(-20/10) = 0.01)
//#define P_START (-30.0)
//#define P_START (-45.0)
/*#define P_STEP  (-10.0/8)  // err probability step (in dB, every times 10^(STEP/10))  // (moved to ldpc_parm.h)*/

    // ###########################################################################################
    int  osdw=OSDW, max_iter = MAX_ITER, max_iter_1;
    double alpha_a = 0, alpha_b = AFP, alpha= AFP;
    //scanf("%d", &max_iter);
    //scanf("%d%d", &max_iter, &osdw);
    //scanf("%d%d%d", &max_iter, &osdw, &alpha_b);
    //scanf("%d%d", &max_iter, &alpha_b);
    //scanf("%d", &osdw);
    max_iter_1 = max_iter + 1;

    uint8_t *finalResult;
    uint8_t hardD, q;
    double sum;
    # if USE_GF2_DEC == 0
        uint8_t usedTo[N];
        double **prob;   // should send to OSD
        uint8_t *LastRun; // should send to OSD
        prob = calloc(N, sizeof(*prob));
        LastRun = calloc(N, sizeof(*LastRun));
        for(int i = 0; i<N; i++){prob[i]= calloc(4, sizeof(*prob[i])); }
    # else
        int N2 = 2*N;
        uint8_t usedToGF2[N2];
        double **probGF2;
        uint8_t *LastRunGF2;
        LastRunGF2 = calloc(N2, sizeof(*LastRunGF2));
        probGF2 = calloc(N2, sizeof(*probGF2));
        for(int i = 0; i<N2; i++){probGF2[i]= calloc(2, sizeof(*probGF2[i])); }
    # endif
    OSD osdDecoder, *osdDec = &osdDecoder;
    FILE *fp;

    # if STORE_WRONG
        char path_saveWrong[1000];
    #endif
    # if SORTCOMPARE == LastRunSORT  && STORE_LF
        uint64_t *storeLF = calloc(max_iter+2, sizeof(*storeLF));
        char path_saveLF[1000];
    # endif // SORTCOMPARE
    #if OSD_RECORD
    OsdRecord osdRecord, *osdRec = &osdRecord;
    #endif // OSD_RECORD

    FILE *fpErr = fopen("err", "r");
    // ##########################################################################################

  int32_t iter, b, m, target;
  GFQ_t nn[N], diff[N], zz[M], zz_G[M_G];
  uint64_t tx_seed[4];

  //- RX
  double p_err, rnd_val=0, p_bias, afp=0;                   if(rnd_val) { } // just to prevent compile warning when CH_TYPE not 0
  FILE *fpA = fopen ( PATH_A , "r" );
  FILE *fpGs = fopen( PATH_Gs, "r" );

  a_matrix_GFQ  Amtx, *A=&Amtx;
  g_matrix_GFQ  Gmtx, *G=&Gmtx;
  uint8_t syndrome_ok, init_syndro_ok;

  //- init proc
  //fast_math_test();
  rnd256_init();    rnd256_init_priv(tx_seed);
  //for(int i = 0; i<N/2; i++)
  //      rnd_val = rv_UnifOne();

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


    char path_BER[1000] ;  //= "BER";     //sprintf(path_BER,  "St_%u_%u_it%u_dec%u_1e-11%s%s%s", N,K,MAX_ITER,BY_DEC,(FIX_P? "_fixP":""),(ALFA? "_sq2":""),(BETA? "_ofs":""));


  char path_St[1000] ;  //= "St";     //sprintf(path_St,  "St_%u_%u_it%u_dec%u_1e-11%s%s%s", N,K,MAX_ITER,BY_DEC,(FIX_P? "_fixP":""),(ALFA? "_sq2":""),(BETA? "_ofs":""));
  char path_Last[1000]; //="StLast"; //sprintf(path_Last,"St_%u_%u_it%u_dec%u_1e-11%s%s%s_Last", N,K,MAX_ITER,BY_DEC,(FIX_P? "_fixP":""),(ALFA? "_sq2":""),(BETA? "_ofs":""));
  char path_Scr[1000] ; //= "Scr";   //sprintf(path_Scr,"Scr_%u_%u_it%u_dec%u_1e-11%s%s%s", N,K,MAX_ITER,BY_DEC,(FIX_P? "_fixP":""),(ALFA? "_sq2":""),(BETA? "_ofs":""));

  if(CH_TYPE==0)      { } // dep ch, no special mark
  else if(CH_TYPE==1) { strcpy(st_name, "XZ_"); }  // indep X-Z ch
  else                { printf("CH_TYPE %d undefined !!  now exit();", CH_TYPE);  getchar();  exit(1); }

  sprintf(st_name, "%s%s_it%u_dec%u", st_name, PRE_A, max_iter, BY_DEC);  // care: st_name in both I/O , behavior not guaranteed
  //sprintf(st_name, "%s_r%f", st_name,RelThr); //display the value of reliability threshold in the output file

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

  }else if(alpha) sprintf(st_name, "%s_%s%.0lf", st_name, "ap", alpha);      // care: st_name in both I/O , behavior not guaranteed

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
  if(POSD_FLIPPING_ALL_BIT)  sprintf(st_name, "%s_FlipAllBits", st_name);

  // --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# if SORTCOMPARE == LastRunSORT
    sprintf(st_name, "%s_LastRun", st_name);
#elif SORTCOMPARE == ENTROPYSORT
    sprintf(st_name, "%s_Entropy", st_name);
#elif SORTCOMPARE == MAXSORT
    sprintf(st_name, "%s_max", st_name);
#elif SORTCOMPARE == BINARYSORT
    sprintf(st_name, "%s_Binary", st_name);
# endif // SORTCOMPARE

  if(osdw >= 0){sprintf(st_name, "%s_OSDW%d", st_name, osdw);}
  else if(osdw == -1){sprintf(st_name, "%s_NoOSD", st_name);}
  else if(osdw == -2){sprintf(st_name, "%s_ADOSD", st_name);}



  # if OPTalpha && USE_GF2_DEC==0
  char path_record[1000];
  sprintf(st_name, "%s_AfpOPT", st_name);
  sprintf(path_record, "RC_%s", st_name);
  FILE *fpRecord;
  fpRecord = fopen(path_record, "w");
  fclose(fpRecord);

  const int aTotal = 20;
  const int maxAlpha = 185;
  int alpha_try1[] = {90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145,  150, 155, 160, 165, 170, 175, 180, 185};
  int alpha_try2[9];
  int *alpha_try;
  int alpha_check[aTotal], accWrong[aTotal];
  int phaseFinish;
  uint32_t record_nCW[aTotal];
  int alpha_length, n_survive, target_wrong, alpha1 = 1000, alpha2 = 1000, alpha_min, nnCW, ii, minError, maxCW, bestalpha;
  double gap;

  time_t TStr = time(NULL), TEnd, Sec = 0, Sec_last = 0;

  double p_try[] = {-11.875, -11.25, -10.625}; // 0.0649, 0.0749, 0.0866 alpha will getting smaller

    initOSD(osdDec);
    #if OSD_RECORD
    initOsdRecord(osdRec);
    #endif // OSD_RECORD

    fp = fopen( PATH_A , "r" );
    load_A_OSD(osdDec, fp);
    fclose(fp);
    osdDec -> RankH  = N- K;
    osdDec -> osdw = 0;

  for(int i = 0; i<6; i++){
        phaseFinish = 0;
        printf("Phase %d alpha opt\n", i);
        nnCW = 0;
        p_err = pow(10, p_try[i/2]/10);
        if( (i&1) == 0){
            target_wrong = 100;
            alpha_try = alpha_try1;
            alpha_min = min(alpha1, alpha2);
            n_survive = 0;
            for(int j = 0; j<aTotal; j++){
                accWrong[j] = 0;
                if(alpha_try1[j] <= alpha_min){
                    alpha_check[j] = 0;
                    n_survive ++;
                }
                else{
                    alpha_check[j] = 1;
                }
            }
            alpha_length = n_survive;
        }
        else{
            target_wrong = 60;
            alpha_length = 9;
            alpha_try = alpha_try2;
            n_survive = 9;
            for(int j = 0; j<9; j++){
                accWrong[j] = 0;
                alpha_check[j] = 0;
            }
        }
        // finish initialize alpha_length, alpha_try, alpha_check, acc_wrong, target_wrong, n_survive

        do {
            for ( b = 0 ; b < N ; b ++ ) {
                  rnd_val = rv_UnifOne();
                  if(rnd_val < p_err/3)        { nn[b] = 3; }
                  else if(rnd_val < p_err/3*2) { nn[b] = 2; }
                  else if(rnd_val < p_err)     { nn[b] = 1; }
                  else                         { nn[b] = 0; }
            }
            Quan_GenSyndrome(A, nn, zz);
            p_bias = p_err;
          for(ii = 0; ii<alpha_length; ii++){

                if(alpha_check[ii])
                    continue;
                for ( b = 0 ; b < N ; b ++ ) {
                    usedTo[b] = 0;
                    LastRun[b] = 1;
                }
                afp =  (alpha_try[ii]==100) ?  0 :(double) 100.0/alpha_try[ii];
                init_syndro_ok = Qbp_init20(bp, zz, A, p_bias, afp);  // contains a subroutine GenInitBiasVec
                syndrome_ok = init_syndro_ok;   iter = 0;

                while ( !init_syndro_ok ) {
                    iter ++;
                    syndrome_ok = Qbp_dec20(bp, A);
                    if(osdw == -1 || !syndrome_ok){
                        for(b = 0; b<N; b++){
                            // ==== calculate LastRun =====
                            hardD = bp->tt[b];
                            if(hardD == usedTo[b]){
                                LastRun[b] = LastRun[b]+1;
                            }
                            else{
                                LastRun[b] = 1;
                                usedTo[b] = hardD;
                            }
                        }
                    }
                    // MAX_ITER
                    if ( syndrome_ok || iter==max_iter )   break;
                }

                if(!syndrome_ok){
                         for(b = 0; b<N; b++){
                            sum=0; for(q=0;q<Q;q++) { sum += bp->qn[b][q]; }
                            for ( q = 0 ; q < Q ; q ++ ) {
                                prob[b][q] = bp->qn[b][q]/sum;
                            }
                        }
                        //finalResult = post_decodeOSDfull(osdDec, zz, prob, LastRun, max_iter);
                        finalResult = post_decodeOSDw(osdDec, zz, prob, LastRun);
                }
                else{     finalResult = bp->tt;     }

                if((HamDist(nn, finalResult, N))!=0) {
                  VecDiff(diff, finalResult, nn, N);     // if bp->tt is a degenerate err of nn,
                  Quan_DegSyndrome(G, diff, zz_G);  // then diff will gen zz_G = all_zero
                  if((is_zero_vec(zz_G, M_G))==0) {
                        accWrong[ii] ++;
                        if(accWrong[ii] == target_wrong){
                            n_survive--;
                            alpha_check[ii] = 1;
                            record_nCW[ii] = nnCW+1;
                        }
                    }
                }
            }

            nnCW++;
            // delete some non necessary part
            if((nnCW % 100)==0){
                minError = target_wrong;
                for(ii = 0; ii<alpha_length; ii++){
                    if(alpha_check[ii] == 0)
                        minError = min(minError, accWrong[ii]);
                }
                minError = (minError+1)*5;
                for(ii = 0; ii<alpha_length; ii++){
                    if((alpha_check[ii] == 0) && (minError <= accWrong[ii])){
                        alpha_check[ii] = 1;
                        n_survive --;
                        record_nCW[ii] = nnCW+1;
                    }
                }
            }

            if(nnCW % (1000) == 0) { TEnd = time(NULL);  Sec = TEnd-TStr;
                if(Sec >= Sec_last+10) {
                        Sec_last=Sec;
                        minError = target_wrong;
                        printf(" now p %g , #CW %llu @ %ld sec \n     ", p_err , nnCW , Sec);
                        for(ii = 0; ii<alpha_length; ii++){
                            if(alpha_check[ii] == 0){
                                printf("%3d %2d    ", alpha_try[ii], accWrong[ii]);
                                minError = min(minError, accWrong[ii]);
                            }
                        }
                        printf("\n Now LER: %lf\n", (double)minError/nnCW);
                }
            }
            if (n_survive <= 1 || nnCW >= 2e6) {
                    phaseFinish = 1;
                    TEnd = time(NULL);  Sec = TEnd-TStr;  Sec_last=Sec;

                    if(n_survive == 0){
                        maxCW = 1;
                        for(ii = 0; ii<alpha_length; ii++){
                            if(maxCW < record_nCW[ii]){
                                bestalpha = alpha_try[ii];
                                maxCW = record_nCW[ii];
                            }
                        }
                    }
                    else{
                        minError = target_wrong<<1;
                        for(ii = 0; ii<alpha_length; ii++){
                            if(alpha_check[ii] == 0){
                                record_nCW[ii] = nnCW;
                                if(minError > accWrong[ii]){
                                    bestalpha = alpha_try[ii];
                                    minError = accWrong[ii];
                                }
                            }
                        }
                    }

                    fpRecord = fopen(path_record, "a");
                    fprintf(fpRecord, "%d\t%d\t%lf\n", i, alpha_length, p_err);
                    for(ii = 0; ii<alpha_length; ii++){
                        fprintf(fpRecord, "%d\t%d\t%d\n", alpha_try[ii], accWrong[ii], record_nCW[ii]);
                    }
                    fclose(fpRecord);

                    printf("Finish %d phase opt best alpha = %d ===============\n", i, bestalpha);
                    if(bestalpha >= maxAlpha){
                        printf("alpha is too small");
                        system("pause");
                        return 1;
                    }

                    if((i&1) == 0){
                        bestalpha -= 4;
                        for(ii = 0; ii<9; ii++)
                            alpha_try2[ii] = bestalpha+ii;
                    }
                    else{
                        switch(i){
                        case 1:
                            alpha1 = bestalpha;
                            break;
                        case 3:
                            alpha2 = bestalpha;
                            break;
                        case 5:
                            alpha_a = (alpha1-alpha2)/(p_try[0]-p_try[1]);
                            alpha_b = alpha1 - alpha_a*p_try[0];
                            if(alpha_a >0) {
                                alpha_a = 0;
                                alpha_b = alpha1;
                            }
                            gap = alpha_a*p_try[2]+alpha_b-bestalpha;
                            printf("\n\n alpha = %lf * log10(p_err) * 10+ %lf, actual diff = %lf\n\n", alpha_a, alpha_b, gap);
                            break;
                        }
                    }
            }
          } while (!phaseFinish);
    }

    fpRecord = fopen(path_record, "a");
    fprintf(fpRecord, "%lf\t%lf\t%ld\n", alpha_a, alpha_b, Sec);
    fclose(fpRecord);

  # endif // OPTalpha

  # if STORE_WRONG
      sprintf(path_saveWrong, "St_%s.Wrong_error", st_name);
      FILE *fpWrong = fopen(path_saveWrong, "w");
      fclose(fpWrong);
  # endif // STORE_WRONG
  # if SORTCOMPARE == LastRunSORT && STORE_LF
      sprintf(path_saveLF, "St_%s.LF", st_name);
      FILE *fpLF = fopen(path_saveLF, "w");
      fclose(fpLF);
  # endif // SORTCOMPARE




    sprintf(path_St,  "Results/St/St_%s", st_name);      sprintf(path_Scr,  "Results/Scrn/St_%s.Scrn", st_name);      sprintf(path_Last,  "Results/Last/St_%s.Last", st_name);
  FILE *fpSt = fopen ( path_St , "a+w" );
  FILE *fpScr = fopen ( path_Scr , "a+w" );
  FILE *fpLast; //= fopen ( path_Last , "a+w" );


  sprintf(path_BER,  "Results/BER_%s.txt", st_name);
  FILE *fpBER = fopen ( path_BER , "w" );

                // Fa=False, Eb=err bits                    // Q=quantum, Eq=err qubits
  uint64_t nCW, nErr, nFa, nEb, nIterAcc, nStop=N_ERR_STOP,  nErrQ, nFaQ, nEq;
  uint32_t nErrBit,  nErrQbit;
  double   bler, ber, far,  LER, qber, qfar,  chk_bler, chk_ber,  aiter;
  double rStep1=0, rStep2=0, tStep1=0, tStep2=0;
                    //far=false rate                              //aiter=avg iter

  time_t  tStr, tEnd, sec, sec_last=0;  // care that: gcc32 has year 2038 problem (gcc64 no this problem)

  nCW=0; nErr=0; nFa=0; nEb=0; nIterAcc=0; nErrBit=0;  nErrQ=0; nFaQ=0; nEq=0; nErrQbit=0;
  p_err = pow(10, P_START/10);  chk_bler=1.0; chk_ber=1.0;
  printf("//==== init p_err = %g ====// N = %u K = %u   It = %u  BY_DEC %u \n", p_err,N,K,MAX_ITER,BY_DEC);
  fprintf(fpScr, "//==== init p_err = %g ====// N = %u K = %u   It = %u  BY_DEC %u \n", p_err,N,K,MAX_ITER,BY_DEC);
  fclose(fpScr);
  tStr = time(NULL);

  # if SORTCOMPARE == LastRunSORT && STORE_LF
        fpLF = fopen(path_saveLF, "a");
        fprintf(fpLF, "%lf ", p_err);
      for(int i = 0; i<max_iter; i++)   fprintf(fpLF, "0 "); fprintf(fpLF, "\n");
      fclose(fpLF);
      memset(storeLF, 0, (max_iter+1)*sizeof(*storeLF));
  # endif // SORTCOMPARE

    initOSD(osdDec);
    fp = fopen( PATH_A , "r" );
    load_A_OSD(osdDec, fp);
    fclose(fp);
    osdDec -> RankH  = N- K;
    osdDec -> osdw = osdw;
    alpha = alpha_a * log10(p_err) * 10 + alpha_b;
    if(alpha <= 50)    alpha = 0;
    //printf("\nNew alpha: %lf\n\n", alpha);
  do {
    //- add noise nn_(Nx1), where (Nx1) means a column vector
    # if 1
//printf("True error\n");
        for ( b = 0 ; b < N ; b ++ ) {
            #if CH_TYPE==0    // dep ch
              rnd_val = rv_UnifOne();
              if(rnd_val < p_err/3)        { nn[b] = 3; }
              else if(rnd_val < p_err/3*2) { nn[b] = 2; }
              else if(rnd_val < p_err)     { nn[b] = 1; }
              else                         { nn[b] = 0; }

            #elif CH_TYPE==1  // indep X-Z ch
              nn[b] = 0;
              if(rv_UnifOne() < p_err)  { nn[b] |= 0b01; }  // set X error
              if(rv_UnifOne() < p_err)  { nn[b] |= 0b10; }  // set Z error
            #endif
        }

    #else
        //printf("True error\n");
        for(b = 0; b<N; b++){
            fscanf(fpErr, "%d", &nn[b]);
        //    printf("%d ", nn[b]);
        }
        //printf("\n");
    #endif // 1

    //-- symplectic syndrome measurement
    Quan_GenSyndrome(A, nn, zz);

//for(b = 0; b<M; b++){  if(zz[b]) printf("%d ", b); }                 printf("\n");
//for(b = 0; b<N; b++){  if(nn[b]) printf("%d %d  ", b, nn[b]); } printf("\n");
//for(b = 0; b<N; b++){  printf("%d ", b%10); }                 printf("\n");
//for(b = 0; b<N; b++){  printf("%d ", nn[b]); }                 printf("\n\n");


  #if FIX_P == 0
    p_bias = p_err*1.0; // do necessary scaling if needed
    //p_bias = p_err*40.0;  if(p_bias>4.6e-2) p_bias=4.6e-2;
    //p_bias = p_err*10.0; // do necessary scaling if needed
  #else   // fix p_bias
    //p_bias = (double)FIX_P/10 * 1e-2;
    p_bias = (double)FIX_P/10 * 1e-2 * (FIX_e? pow(10.0, FIX_e+3) : 1.0);  // when FIX_e=-2, 1/10*1e-2*10^(FIX_e+3) = 10^-2
    //p_bias = 4.6e-2;   // fix p_bias
  #endif

  #if ZE
    double p_est = (double) HamWt(zz,M)*ZE/100 / N;
    if(p_est > p_err) p_bias = p_est;
    else              p_bias = p_err;
  #endif

    // _______________________________________________________________________________________________________________
    #if SORTCOMPARE == LastRunSORT && OSDW  != -1
        #if USE_GF2_DEC == 0
        for ( b = 0 ; b < N ; b ++ ) {
            usedTo[b] = 0;
            LastRun[b] = 0;
        }
        #else
        for ( b = 0 ; b < N2 ; b ++ ) {
            usedToGF2[b] = 0;
            LastRunGF2[b] = 0;
        }
        #endif
    #endif // SORTCOMPARE
    // _______________________________________________________________________________________________________________


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

    // ==================================================================================================
    #if SORTCOMPARE == LastRunSORT
    if(osdw != -1 && !syndrome_ok){
        for(b = 0; b<N2; b++){
            // ==== calculate LastRun =====
            if(bp_GF2->tt[b] != usedToGF2[b]){
                usedToGF2[b] = bp->tt[b];
                LastRunGF2[b] = iter;
            }
        }
    }

    #endif // SORTCOMPARE
      // ==================================================================================================

    #if LIST_DEC == 0   // DEFAULT_MODE
      if ( syndrome_ok || iter==max_iter )   break;
    #elif LIST_DEC == 1 // BUF_LAST
      if(syndrome_ok) { n_list_sol = 1;   memcpy(ListSols_GF2[0], bp_GF2->tt, sizeof(ListSols_GF2[0])); }  // buffer the last solution
      if(iter==max_iter)  break;
      //if(n_list_sol > 0)  break;  // enable this row if do TEST to match the original behavior (when LIST_DEC=0)
    #endif

    }
    nIterAcc += iter;
    //printf("BP_dec syndrome %s  iter %d  TX-RX vec HamDist %d \n", syndrome_ok? "OK!":"FAIL!!!" , iter , HamDist(nn2, bp_GF2->tt, 2*N));

    #if LIST_DEC == 1   // BUF_LAST
    if(n_list_sol > 0) { memcpy(bp_GF2->tt, ListSols_GF2[0], sizeof(bp_GF2->tt)); }
    #endif

    // ==================================================================================================
    if(osdw == -1 || syndrome_ok){ // no OSD or syndrome is OK
        for(b=0; b<N ; b++) { bp->tt[b] = bp_GF2->tt[b] + bp_GF2->tt[b+N]*2; }    // pass result to bp->tt for statistics
        finalResult = bp->tt;
    }

    else{                 // Do OSDw or OSDfull
        for(b = 0; b<N2; b++){
            LastRun[b] = max_iter_1 - LastRun[b];
            sum = bp_GF2->qn[b][0]+bp_GF2->qn[b][1];
            probGF2[b][0] = bp_GF2->qn[b][0]/sum;
            probGF2[b][1] = bp_GF2->qn[b][1]/sum;
        }
        if(osdw == -2){
             finalResult = post_decodeOSDfull(osdDec, zz, probGF2, LastRunGF2, max_iter);
        }
        else{
            finalResult = post_decodeOSDw(osdDec, zz, probGF2, LastRunGF2);
        }
        # if SORTCOMPARE == LastRunSORT && STORE_LF
        for(int i = 0; i<N2; i++)
            storeLF[LastRunGF2[i]] ++;
        #endif
    }
    // ==================================================================================================

  #else // USE_GF2_DEC==0   // (default mode, decoding by GF4) //###########################################################
    //%%%%%%%%%%%% DEFAULT_RUN that could be with AFP (without AFP2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // afp = (AFP==0)?  0 : (double)100.0/AFP;
    afp =  (alpha==0 || alpha == 100) ?  0 :(double) 100.0/alpha;
    // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

    // ==================================================================================================
    #if SORTCOMPARE == LastRunSORT && OSDW !=-1
    if(osdw == -1 || !syndrome_ok){
        if(osdw != -1 && !syndrome_ok){
            for(b = 0; b<N; b++){
                // ==== calculate LastRun =====
                if(bp->tt[b] != usedTo[b]){
                    usedTo[b] = bp->tt[b];
                    LastRun[b] = iter;
                }
            }
        }
         //for(b = 0; b<N; b++){ printf("%d ", bp->tt[b]);}            printf("\n");
    }
    #endif // SORTCOMPARE
      // ==================================================================================================
    // MAX_ITER
    #if LIST_DEC == 0   // DEFAULT_MODE
      if ( syndrome_ok || iter==max_iter )   break;
    #elif LIST_DEC == 1 // BUF_LAST
      if(syndrome_ok) { n_list_sol = 1;   memcpy(ListSols, bp->tt, sizeof(ListSols[0])); }  // buffer the last solution
      if(iter==MAX_ITER)  break;
      //if(n_list_sol > 0)  break;  // enable this row if do TEST to match the original behavior (when LIST_DEC=0)
    #endif
    }
    //###################################################################################################
    if(osdw == -1 || syndrome_ok){  // no OSD or syndrome is OK
        finalResult = bp->tt;
    }
    else{                                               // Do OSDw or OSDfull
        // ==== calculate probability ===
        //printf("%d ", nCW);
        for(b = 0; b<N; b++){
            LastRun[b] = max_iter_1 - LastRun[b];
            sum=0; for(q=0;q<Q;q++) { sum += bp->qn[b][q]; }
            for ( q = 0 ; q < Q ; q ++ ) {
                prob[b][q] = bp->qn[b][q]/sum;
            }
        }
        if(osdw == -2){
            finalResult = post_decodeOSDfull(osdDec, zz, (const double **) prob, LastRun, max_iter);
        }
        else{
            finalResult = post_decodeOSDw(osdDec, zz, (const double **) prob, LastRun);
        }

        # if SORTCOMPARE == LastRunSORT && STORE_LF && OSDW !=-1
        for(int i = 0; i<N; i++){
            storeLF[LastRun[i]] ++;
        }
        #endif

{
    //for(b = 0; b<N; b++){  printf("%d ", b%10); }                 printf("\n");
    //for(b = 0; b<N; b++){ printf("%d ", bp->tt[b]);}            printf("\n");
    //for(b = 0; b<N; b++){ printf("%d ", finalResult[b]);}       printf("\n");
    //for(b = 0; b<N; b++){ printf("%d ", nn[b]);}                    printf("\n");

    //if(!weirdFlag)
    //    system("pause");
}
    }
    //###################################################################################################

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
    #if 0
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
    }
    #elif 1
    nCW++;
    nErrBit = HamDist(nn, finalResult, N);

      if(nErrBit==0) { }  // fine only when both syndrome and HamDist pass
      else { nErr++; nEb+=nErrBit; nFa++; }

    //-- "Quantum" degeneracy check
    nErrQbit = nErrBit;   // start with classical error
    if(nErrQbit!=0) {
      VecDiff(diff, finalResult, nn, N);     // if bp->tt is a degenerate err of nn,
      Quan_DegSyndrome(G, diff, zz_G);  // then diff will gen zz_G = all_zero
      //Quan_LogSyndrome(G, diff, zz_G);  // then diff will gen zz_G = all_zero (here only check last 2K logical rows)
      if(is_zero_vec(zz_G, M_G)) { nErrQbit = 0; } // force 0 if it's degenerate error
    }
    if(nErrQbit!=0) {
        nErrQ++;  nEq+=nErrQbit;
        nFaQ++;
        //=========================================
        # if STORE_WRONG
        fpWrong = fopen(path_saveWrong, "a");
            for(int i = 0; i<N; i++){
                fprintf(fpWrong, "%d ", nn[i]);
//printf("%d ", nn[i]);
            }
            fprintf(fpWrong, "\n");
//printf("\n");
//printf("Wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
      fclose(fpWrong);
        #endif // STORE_WRONG
        // ========================================
        //printf("Damn...\n");
        //system("pause");
    }
    else{
        //printf("Yes\n");
    }
    #endif // 0

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


    if(nCW % (1000*1) == 0) { tEnd = time(NULL);  sec = tEnd-tStr;
        # if SORTCOMPARE == LastRunSORT && STORE_LF && OSDW !=-1
        if(nCW % 1000000 == 0){
            fpLF = fopen(path_saveLF, "a");
            for(int i = 1; i<=max_iter+1; i++)    fprintf(fpLF, "%llu ", storeLF[i]);
            fprintf(fpLF, "\n");
            fclose(fpLF);
        }
        #endif
        if(sec >= sec_last+10) {
          fpScr= fopen ( path_Scr , "a+w" );    fpLast = fopen ( path_Last , "a+w" );   //fpSt = fopen ( path_St , "a+w" );
          tEnd = time(NULL);  sec = tEnd-tStr;  sec_last=sec;

        bler=(double)nErr/nCW;    ber=(double)nEb/(nCW*N);  far=(double)nFa/nCW;      aiter=(double)nIterAcc/nCW;
        LER=(double)nErrQ/nCW;  qber=(double)nEq/(nCW*N); qfar=(double)nFaQ/nCW;





                printf("    now p %-10g  #CW %-10"PRIu64"  #LErr %-6"PRIu64"  => LER=%-10.2e  AvgIter=%.2f  @ %10ld sec\n",
        p_err, nCW, nErrQ, LER, aiter, sec);



          fprintf(fpScr, " now p %g , #CW %llu     #ErrQ %llu => LER %g    BLER %g    AvgIter %g  @ %ld sec\n",
                               p_err , nCW,            nErrQ,      LER,    bler, aiter, sec);

           fprintf(fpLast,"%g\t%g\t%g\t%g\t%g\t%llu\t%llu\t%llu\t%ld\t%llu\t%llu\t%g\t%g\t%g\n",
        p_err,bler,far,ber,aiter, nCW,nErr,nEb, sec, nEq,nErrQ,qber,LER,qfar);

          if(nErrQ >= N_ERR_STOP_1em7){
            if(1 > LER && LER > 1e-1) nStop = N_ERR_ACCURATE;
            if(LER <= 1e-1)  nStop = N_ERR_STOP_1em1;  // start to accelerate
            if(LER <= 1e-2)  nStop = N_ERR_STOP_1em2;  // start to accelerate more ...
            if(LER <= 1e-3)  nStop = N_ERR_STOP_1em3;
            if(LER <= 1e-4)  nStop = N_ERR_STOP_1em4;
            if(LER <= 1e-5)  nStop = N_ERR_STOP_1em5;
            if(LER <= 1e-6)  nStop = N_ERR_STOP_1em6;
            if(LER <= 1e-7)  nStop = N_ERR_STOP_1em7;
          }

          fclose(fpScr);    fclose(fpLast);   //fclose(fpSt);
        }
    }

    if (nErrQ >= nStop) {
//NEXT:
      fpScr= fopen ( path_Scr , "a+w" );    fpLast = fopen ( path_Last , "a+w" );   //fpSt = fopen ( path_St , "a+w" );
      tEnd = time(NULL);  sec = tEnd-tStr;

      bler=(double)nErr/nCW;    ber=(double)nEb/(nCW*N);  far=(double)nFa/nCW;      aiter=(double)nIterAcc/nCW;
      LER=(double)nErrQ/nCW;  qber=(double)nEq/(nCW*N); qfar=(double)nFaQ/nCW;


      chk_bler = LER;   chk_ber = qber;



       printf("\np_err %-10g  LER %-10.2e   AvgIter %.2f  #CW %-10"PRIu64"  #LErr %-6"PRIu64"    @ %10ld sec \n\n",
        p_err,LER, aiter , nCW, nErrQ, sec);


      fprintf(fpScr, "\np_err %g  LER %g qBER %g  BLER %g BER %g  AvgIter %.2f  for nCW %llu nErr %llu nEb %llu  nErrQ %llu nEq %llu  @ %ld sec\n",
        p_err,LER,qber,bler,ber,aiter , nCW,nErr,nEb,nErrQ,nEq , sec);

       fprintf(fpSt,  "%g\t%g\t%g\t%g\t%g\t%llu\t%llu\t%llu\t%ld\t%llu\t%llu\t%g\t%g\t%g\n",
        p_err,bler,far,ber,aiter, nCW,nErr,nEb, sec, nEq,nErrQ,qber,LER,qfar);
      fprintf(fpLast,"%g\t%g\t%g\t%g\t%g\t%llu\t%llu\t%llu\t%ld\t%llu\t%llu\t%g\t%g\t%g\n",
        p_err,bler,far,ber,aiter, nCW,nErr,nEb, sec, nEq,nErrQ,qber,LER,qfar);

     fprintf(fpBER,  "%g\t %10e\n", p_err,LER);



      nCW=0; nErr=0; nFa=0; nEb=0; nIterAcc=0; nErrBit=0;  nErrQ=0; nFaQ=0; nEq=0; nErrQbit=0;

        // ==========================
      double step;
      (LER>1e-5)? (step = P_STEP):(step = P_STEP/2);
      p_err *= pow(10, step/10);
    alpha = alpha_a * log10(p_err) * 10 + alpha_b;
    if(alpha <= 50)    alpha = 0;
    //printf("\nNew alpha: %lf\n\n", alpha);

    # if SORTCOMPARE == LastRunSORT && STORE_LF && OSDW !=-1
    fpLF = fopen(path_saveLF, "a");
    for(int i = 1; i<=max_iter+1; i++)    fprintf(fpLF, "%llu ", storeLF[i]);
        fprintf(fpLF, "\n");
        fprintf(fpLF, "%lf ", p_err);
      for(int i = 0; i<max_iter; i++)   fprintf(fpLF, "0 "); fprintf(fpLF, "\n");
      memset(storeLF, 0, (max_iter+1)*sizeof(*storeLF));
      fclose(fpLF);
  # endif // SORTCOMPARE
      //p_err *= pow(10, P_STEP/10);
        // ==========================

      if(1 > LER && LER > 1e-1) nStop = N_ERR_ACCURATE;
      if(LER <= 1e-1)  nStop = N_ERR_STOP_1em1;  // start to accelerate
      if(LER <= 1e-2)  nStop = N_ERR_STOP_1em2;  // start to accelerate more ...
      if(LER <= 1e-3)  nStop = N_ERR_STOP_1em3;
      if(LER <= 1e-4)  nStop = N_ERR_STOP_1em4;
      if(LER <= 1e-5)  nStop = N_ERR_STOP_1em5;
      if(LER <= 1e-6)  nStop = N_ERR_STOP_1em6;
      if(LER <= 1e-7)  nStop = N_ERR_STOP_1em7;

      fclose(fpScr);    fclose(fpLast);   //fclose(fpSt);

      //if(LER <= 1e-3) { if(rec_cnt==-1)  rec_cnt=0; }  // enable record cnt
    }

  } while (chk_bler>BLER_STOP || chk_ber>BER_STOP);


    fprintf(fpBER,  "%s\n", path_BER);
  fclose(fpBER);

  fclose(fpA);
  fclose(fpGs);
	printf("\nMain function done!\n");
  fpScr= fopen ( path_Scr , "a+w" );
  fprintf(fpScr, "main func done!\n");
  fclose(fpSt); fclose(fpLast); fclose(fpScr);
  free_OSD(osdDec);
//#if LOG_ZE_ai
//  fclose(fpZE);
//#endif
	getchar();
	return 0;
}
