#include <math.h>
#include "./bp_dec.h"
#include "./bp_llr.h"
#include "../lib_rand/lib_rand.h"
#include "../lib_math/fast_math.h"

//#include <time.h>   // test rand runtime


//-- alloc for LLR-based Quantum BP Ctrl Blk
void alloc_LBPC(a_matrix_GFQ *A, LBP_Ctl *bp)
{
  int32_t  m; //, n;

  for(m=0; m < A->MM; m++) {
    //bp->GAmj[m]= calloc(A->num_m[m], sizeof( *(bp->GAmj[m]) ) );  //printf("sizeof( *(bp->GAmj[m]) ) = %u \n", sizeof( *(bp->GAmj[m]) ));
    bp->Lmj[m] = calloc(A->num_m[m], sizeof( *(bp->Lmj[m] ) ) );  //printf("sizeof( *(bp->Lmj[m]) )  = %u \n", sizeof( *(bp->Lmj[m]) ));
    bp->Dmj[m] = calloc(A->num_m[m], sizeof( *(bp->Dmj[m] ) ) );  //printf("sizeof( *(bp->Dmj[m]) )  = %u \n", sizeof( *(bp->Dmj[m]) )); getchar();
  #if 0 //EN_AFP_V
    bp->dqml_ori[m] = calloc(A->num_m[m], sizeof( *(bp->dqml_ori[m]) ) );
  #endif
  #if 0 //STOP_BY_C
    bp->dqml_c[m] = calloc(A->num_m[m], sizeof( *(bp->dqml_c[m]) ) );
  #endif
  }

  //====================================================
//#if 0 //CYC_CHK
//  int32_t i, j, i2, m2, j2, n2, ii,ii2;
//
//  for(n=0; n < A->NN; n++){
//    bp->cyc_cnt[n] = calloc( A->num_n[n] , sizeof( *(bp->cyc_cnt[n]) ) );
//    bp->cycs[n] = calloc( A->num_n[n] , sizeof( *(bp->cycs[n]) ) );
//    for(i=0 ; i < A->num_n[n]; i++)  bp->cycs[n][i] = calloc( A->num_n[n], sizeof(*(bp->cycs[n][i])) );
//  }
//  for(m=0; m < A->MM; m++) { 
//    bp->rml2[m] = calloc(A->num_m[m], sizeof( *(bp->rml2[m]) ) );
//  }
//
//  // ======== to init bp->cyc_cnt ======== //
//  //
//  // n        n2
//  // m ------ m  ... (run all n2)
//  // |        |
//  // m2------ 1?
//  // :
//  // and (run all m2)
//
//  for(n=0 ; n < A->NN; n++)
//  for(i=0 ; i < A->num_n[n]-1; i++) {
//    m = A->nlist[n][i];
//    //----------- (run all m2)
//    for(i2=i+1; i2 < A->num_n[n]; i2++) {
//      m2 = A->nlist[n][i2];
//      j = A->ni2j[n][i] + 1;    // +1 to get first n2
//      j2 = A->ni2j[n][i2] + 1;  // +1 to (if match will) get n2
//      for( ; j < A->num_m[m]; j++) {  // (run all n2)
//        n2 = A->mlist[m][j];
//        for( ; j2 < A->num_m[m2]; j2++) // check if hit 1? (as above fig): if hit => a four-cycle at bit n's i,i2 
//          if(A->mlist[m2][j2]==n2) { ii = A->mj2i[m][j];  ii2 = A->mj2i[m2][j2];   // need to convert to bit n2's ii,ii2 
//            bp->cyc_cnt[n][i]++; bp->cyc_cnt[n][i2]++; bp->cyc_cnt[n2][ii]++; bp->cyc_cnt[n2][ii2]++;}
//          else if(A->mlist[m2][j2]>n2) break; // to move to next j (i.e. next n2)
//      }
//    }
//    //-----------------------
//  }
//
//  #if 0 // set to 0 to disable
//  for(n=0; n < A->NN; n++) { // debug print
//    for(i=0 ; i < A->num_n[n]; i++) { printf("%d ", bp->cyc_cnt[n][i]); }
//    printf(" = bp->cyc_cnt[%d]\n", n);  //getchar();
//  } getchar();
//  #endif
//
//  // ======== to init bp->cycs (a little bit similar to above) ======== //
//  for(n=0 ; n < A->NN; n++)
//  for(i=0 ; i < A->num_n[n]-1; i++) {
//    m = A->nlist[n][i];
//    //----------- (run all m2)
//    for(i2=i+1; i2 < A->num_n[n]; i2++) {
//      m2 = A->nlist[n][i2];
//      j = A->ni2j[n][i] + 1;    // +1 to get first n2
//      j2 = A->ni2j[n][i2] + 1;  // +1 to (if match will) get n2
//      for( ; j < A->num_m[m]; j++) {  // (run all n2)
//        n2 = A->mlist[m][j];
//        for( ; j2 < A->num_m[m2]; j2++) // check if hit 1? (as above fig): if hit => a four-cycle at bit n's i,i2 
//          if(A->mlist[m2][j2]==n2) { ii = A->mj2i[m][j];  ii2 = A->mj2i[m2][j2];   // need to convert to bit n2's ii,ii2 
//            bp->cycs[n][i][i2]++; bp->cycs[n][i2][i]++; bp->cycs[n2][ii][ii2]++; bp->cycs[n2][ii2][ii]++;}
//          else if(A->mlist[m2][j2]>n2) break; // to move to next j (i.e. next n2)
//      }
//    }
//    //-----------------------
//  }
//
//  #if 0 // set to 0 to disable
//  for(n=0; n < A->NN; n++) { // debug print
//    for(i=0; i < A->num_n[n]; i++) {
//      for(i2=0; i2 < A->num_n[n]; i2++) { printf("%d ", bp->cycs[n][i][i2]); }
//      printf("   ");
//    }
//    printf(" = bp->cycs[%d]\n", n);  //getchar();
//  } getchar();
//  #endif
//
//#endif  // END of CYC_CHK
//
//
//#if 0 //BY_DEC != 25    // default serial order
//  for(n=0; n<N; n++) { bp->seq[n] = n; }
//
//  #if RND_SCHE == 1
//  //time_t  tStr;  tStr = time(NULL);   for(n=0; n<0xfffffff; n++) { m = rv_Nminus1(N); } // a discrete r.v. over {0, 1, 2, ..., N-1}
//  //printf("rv_Nminus1 runtime : %ld sec \n", time(NULL) - tStart); getchar();
//
//  int32_t n1, tmp_n;   uint64_t seq_seed[4];   rnd256_init_priv(seq_seed);
//  for(n=0; n<N-1; n++) {                  //-- cite https://danluu.com/sattolo/
//    //n1 = priv_Nminus1(N, seq_seed);       // a discrete r.v. over {0, 1, 2, ..., N-1}   // Fisher-Yates approach (may multi-cycles)
//    n1 = priv_IntAB(n+1, N-1, seq_seed);  // a discrete r.v. over {n+1,n+2, ..., N-1}   // Sattolo approach  (guaranteed one cycle)
//    tmp_n = bp->seq[n];   bp->seq[n] = bp->seq[n1];   bp->seq[n1] = tmp_n;  // swap idxs at [n] and [n1]
//  }
//  #endif
//
//#elif 0 //BY_DEC == 25
//  #if ORI==4115 || ORI==13134   // [[13,1,3]] code, Surface Code, X (M=N-1), extended 4 redundant rows
//          // [[41,1,5]] code, Surface Code, X (M=N-1)
//    #if Surf_L == 2
//  int32_t seq[N] = {0,1,2,7,12,11,10,5,3,4,9,8,6};   // circ
//    #elif Surf_L == 4
//  int32_t seq[N] = {0,1,2,3,4, 13,22,31,40, 39,38,37,36, 27,18,9, 5,6,7,8, 17,26,35, 34,33,32, 23,14, 10,11,12, 21,30, 29,28, 19, 15,16, 25,24,20};   // circ
//  //int32_t seq[N] = {0,4,40,36, 1,2,3, 13,22,31, 39,38,37, 27,18,9, 5,6,7,8, 17,26,35, 34,33,32, 23,14, 10,11,12, 21,30, 29,28, 19, 15,16, 25,24,20};  // corner first
//  //int32_t seq[N] = {0,1,2,3,4, 13,22,31,40, 39,38,37,36, 27,18,9, 5,6,7,8, 10,11,12, 14,15,16,17, 19,20,21, 23,24,25,26, 28,29,30, 32,33,34,35};      // cseq: circ then seq
//  //int32_t seq[N] = {20,24,25,16,15,19,28,29,30,21,12,11,10,14,23,32,33,34,35,26,17,8,7,6,5,9,18,27,36,37,38,39,40,31,22,13,4,3,2,1,0};   // circ_inv (tested, not good)
//    #elif Surf_L == 6
//  int32_t seq[N] = {0,1,2,3,4,5,6,19,32,45,58,71,84,83,82,81,80,79,78,65,52,39,26,13,7,8,9,10,11,12,25,38,51,64,77,76,75,74,73,72,59,46,33,20,14,15,16,17,18,31,44,57,70,69,68,67,66,53,40,27,21,22,23,24,37,50,63,62,61,60,47,34,28,29,30,43,56,55,54,41,35,36,49,48,42};   // circ
//    #elif Surf_L == 8
//  int32_t seq[N] = {0,1,2,3,4,5,6,7,8,25,42,59,76,93,110,127,144,143,142,141,140,139,138,137,136,119,102,85,68,51,34,17,9,10,11,12,13,14,15,16,33,50,67,84,101,118,135,134,133,132,131,130,129,128,111,94,77,60,43,26,18,19,20,21,22,23,24,41,58,75,92,109,126,125,124,123,122,121,120,103,86,69,52,35,27,28,29,30,31,32,49,66,83,100,117,116,115,114,113,112,95,78,61,44,36,37,38,39,40,57,74,91,108,107,106,105,104,87,70,53,45,46,47,48,65,82,99,98,97,96,79,62,54,55,56,73,90,89,88,71,63,64,81,80,72};   // circ
//    #elif Surf_L == 10
//  int32_t seq[N] = {0,1,2,3,4,5,6,7,8,9,10,31,52,73,94,115,136,157,178,199,220,219,218,217,216,215,214,213,212,211,210,189,168,147,126,105,84,63,42,21,11,12,13,14,15,16,17,18,19,20,41,62,83,104,125,146,167,188,209,208,207,206,205,204,203,202,201,200,179,158,137,116,95,74,53,32,22,23,24,25,26,27,28,29,30,51,72,93,114,135,156,177,198,197,196,195,194,193,192,191,190,169,148,127,106,85,64,43,33,34,35,36,37,38,39,40,61,82,103,124,145,166,187,186,185,184,183,182,181,180,159,138,117,96,75,54,44,45,46,47,48,49,50,71,92,113,134,155,176,175,174,173,172,171,170,149,128,107,86,65,55,56,57,58,59,60,81,102,123,144,165,164,163,162,161,160,139,118,97,76,66,67,68,69,70,91,112,133,154,153,152,151,150,129,108,87,77,78,79,80,101,122,143,142,141,140,119,98,88,89,90,111,132,131,130,109,99,100,121,120,110};   // circ
//    #elif Surf_L == 12
//  int32_t seq[N] = {0,1,2,3,4,5,6,7,8,9,10,11,12,37,62,87,112,137,162,187,212,237,262,287,312,311,310,309,308,307,306,305,304,303,302,301,300,275,250,225,200,175,150,125,100,75,50,25,13,14,15,16,17,18,19,20,21,22,23,24,49,74,99,124,149,174,199,224,249,274,299,298,297,296,295,294,293,292,291,290,289,288,263,238,213,188,163,138,113,88,63,38,26,27,28,29,30,31,32,33,34,35,36,61,86,111,136,161,186,211,236,261,286,285,284,283,282,281,280,279,278,277,276,251,226,201,176,151,126,101,76,51,39,40,41,42,43,44,45,46,47,48,73,98,123,148,173,198,223,248,273,272,271,270,269,268,267,266,265,264,239,214,189,164,139,114,89,64,52,53,54,55,56,57,58,59,60,85,110,135,160,185,210,235,260,259,258,257,256,255,254,253,252,227,202,177,152,127,102,77,65,66,67,68,69,70,71,72,97,122,147,172,197,222,247,246,245,244,243,242,241,240,215,190,165,140,115,90,78,79,80,81,82,83,84,109,134,159,184,209,234,233,232,231,230,229,228,203,178,153,128,103,91,92,93,94,95,96,121,146,171,196,221,220,219,218,217,216,191,166,141,116,104,105,106,107,108,133,158,183,208,207,206,205,204,179,154,129,117,118,119,120,145,170,195,194,193,192,167,142,130,131,132,157,182,181,180,155,143,144,169,168,156};   // circ
//    #elif Surf_L == 14
//  int32_t seq[N] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,43,72,101,130,159,188,217,246,275,304,333,362,391,420,419,418,417,416,415,414,413,412,411,410,409,408,407,406,377,348,319,290,261,232,203,174,145,116,87,58,29,15,16,17,18,19,20,21,22,23,24,25,26,27,28,57,86,115,144,173,202,231,260,289,318,347,376,405,404,403,402,401,400,399,398,397,396,395,394,393,392,363,334,305,276,247,218,189,160,131,102,73,44,30,31,32,33,34,35,36,37,38,39,40,41,42,71,100,129,158,187,216,245,274,303,332,361,390,389,388,387,386,385,384,383,382,381,380,379,378,349,320,291,262,233,204,175,146,117,88,59,45,46,47,48,49,50,51,52,53,54,55,56,85,114,143,172,201,230,259,288,317,346,375,374,373,372,371,370,369,368,367,366,365,364,335,306,277,248,219,190,161,132,103,74,60,61,62,63,64,65,66,67,68,69,70,99,128,157,186,215,244,273,302,331,360,359,358,357,356,355,354,353,352,351,350,321,292,263,234,205,176,147,118,89,75,76,77,78,79,80,81,82,83,84,113,142,171,200,229,258,287,316,345,344,343,342,341,340,339,338,337,336,307,278,249,220,191,162,133,104,90,91,92,93,94,95,96,97,98,127,156,185,214,243,272,301,330,329,328,327,326,325,324,323,322,293,264,235,206,177,148,119,105,106,107,108,109,110,111,112,141,170,199,228,257,286,315,314,313,312,311,310,309,308,279,250,221,192,163,134,120,121,122,123,124,125,126,155,184,213,242,271,300,299,298,297,296,295,294,265,236,207,178,149,135,136,137,138,139,140,169,198,227,256,285,284,283,282,281,280,251,222,193,164,150,151,152,153,154,183,212,241,270,269,268,267,266,237,208,179,165,166,167,168,197,226,255,254,253,252,223,194,180,181,182,211,240,239,238,209,195,196,225,224,210};   // circ
//    #else
//  XXX_NOT_SUPPORT_XXX
//    #endif
//  memcpy(bp->seq , seq , sizeof(bp->seq));
//
//  #elif ORI==2515 // [[25,1,5]] code, Surface Code ROTATED (also support XZZX)
//    /*#if Surf_d == 3
//  int32_t seq[N] = {1,3,7,9, 2,8, 4,6, 5};  // 4 groups
//    #elif Surf_d == 5
//  int32_t seq[N] = {1,3,5,11,13,15,21,23,25, 2,4,12,14,22,24, 6,8,10,16,18,20, 7,9,17,19};  // 4 groups
//    #elif Surf_d == 7
//  int32_t seq[N] = {1,3,5,7,15,17,19,21,29,31,33,35,43,45,47,49, 2,4,6,16,18,20,30,32,34,44,46,48, 8,10,12,14,22,24,26,28,36,38,40,42, 9,11,13,23,25,27,37,39,41};  // 4 groups
//    #else
//  XXX_NOT_SUPPORT_XXX
//    #endif */
//  //-- now auto-done as follows
//  int32_t seq[N];   int32_t  gp,row,col;  // partition into 4 groups
//  n=0;
//  for(gp=0; gp < 4; gp++) {
//    for(row=gp/2; row < Surf_d; row+=2) 
//    for(col=gp%2; col < Surf_d; col+=2) {
//      seq[n] = row*Surf_d + col + 1;    n++;      //printf("%d ", n);
//    }
//  }                                               //printf("\n");
//  //for(n=0;n<N;n++) { printf("%d ", seq[n]);}      printf("\n"); getchar();
//  for(n=0;n<N;n++) { seq[n] = seq[n] - 1; }
//  memcpy( bp->seq , seq , sizeof(bp->seq) );
//
//  #else
//  XXX_NOT_SUPPORT_XXX
//  #endif
//#endif
}



static void HardDecisionLLR(LBP_Ctl *bp, a_matrix_GFQ *A)
{
  int32_t  n;
  for(n=0; n < A->NN; n++) {
    if(bp->GA[n][0] > 0)
      if(bp->GA[n][1] > 0)
        if(bp->GA[n][2] > 0)  bp->tt[n] = 0;
        else                  bp->tt[n] = 3;
      else  // bp->GA[n][1] < 0
        if(bp->GA[n][2] < bp->GA[n][1]) bp->tt[n] = 3;
        else                            bp->tt[n] = 2;
    else// bp->GA[n][0] < 0
      if(bp->GA[n][1] < bp->GA[n][0])
        if(bp->GA[n][2] < bp->GA[n][1]) bp->tt[n] = 3;
        else                            bp->tt[n] = 2;
      else//bp->GA[n][0]<bp->GA[n][1]
        if(bp->GA[n][2] < bp->GA[n][0]) bp->tt[n] = 3;
        else                            bp->tt[n] = 1;
  }
}





//-- restrict LLR too-large
#define LLR_MAX (35)  // IEEE754 double tanh(35.9/2) = 0.999999999999999 , but tanh(36/2) = 1
#define prevent_llr_ov(llr) { \
  if(llr > LLR_MAX)     { llr =  LLR_MAX; } \
  else if(llr<-LLR_MAX) { llr = -LLR_MAX; } \
}
//-- prevent LLR to be too close to zero
#define LLR_MIN (1.0e-10)  // support up to row weight 30, since tanh(1e-10/2)^30 = 9.313225746154796e-310 ~= |DBL_MIN|
#define prevent_llr_0(llr) { \
  if(llr < LLR_MIN) { \
    if(llr >= 0)          { llr =  LLR_MIN; } \
    else if(llr>-LLR_MIN) { llr = -LLR_MIN; } \
  } \
}
//-- restrict LLR in reasonable range
#define restrict_llr(llr) { \
  if(llr > LLR_MAX) { llr = LLR_MAX; } \
  else if(llr < LLR_MIN) { \
    if(llr >= 0)          { llr =  LLR_MIN; } \
    else if(llr>-LLR_MIN) { llr = -LLR_MIN; } \
    else if(llr<-LLR_MAX) { llr = -LLR_MAX; } \
  } \
}


#if 0 //STOP_BY_C   //: syndrome generated by Check Nodes (refer to ChkNode_GenSyndrome();)
static void ChkNode_GenSyndrByLLR(a_matrix_GFQ *A, LBP_Ctl *bp)
{
  int32_t  m, n, j, q;

  for(m=0; m < A->MM; m++) {
    bp->dm_c[m] = 1.0;                  // reset dm_c
    for(j=0; j < A->num_m[m]; j++) {
      n = A->mlist[m][j];
      q = A->mlist_val[m][j]; // {X,Z,Y} val {1,2,3}

      switch (q) {
      case 1: qlm[0]=(bp->qn[n][0]+bp->qn[n][1]); qlm[1]=(bp->qn[n][2]+bp->qn[n][3]); break;
      case 2: qlm[0]=(bp->qn[n][0]+bp->qn[n][2]); qlm[1]=(bp->qn[n][1]+bp->qn[n][3]); break;
      case 3: qlm[0]=(bp->qn[n][0]+bp->qn[n][3]); qlm[1]=(bp->qn[n][1]+bp->qn[n][2]); break;
      }

    #if 1 // SOFT dqml_c: more flexible, for if need to do more computation by dqml_c
      qlm[0] = qlm[0]/(qlm[0] + qlm[1]);
      bp->dqml_c[m][j] = qlm[0]*2 - 1;
    #else // HARD dqml_c: more simpler, for if only need to do judgment by sign of dqml_c
      if(qlm[0]>qlm[1]) bp->dqml_c[m][j] =  1;
      else              bp->dqml_c[m][j] = -1;
    #endif

      bp->dm_c[m] *= bp->dqml_c[m][j];  // update dm_c
    }

    if(bp->dm_c[m]>0) bp->zz_c[m] = 0;
    else              bp->zz_c[m] = 1;
  }

  if( 0 == HamDist(bp->zz_c, bp->target_z, A->MM) )   bp->is_hit_c = 1;
  else                                                bp->is_hit_c = 0;
}
#endif


#define h_exp exp                   // ideal
#define h_log log                   // ideal
#define h_tanh2(a)  (tanh(a/2.0))   // ideal
#define h_atanh2(x) (atanh(x)*2)    // ideal
//-- fast lib: seems only 10% fast, 
/*#define h_exp exp_fast
#define h_log log_fast
#define h_tanh2(a)  ((h_exp(a)-1.0)/(h_exp(a)+1.0))
#define h_atanh2(x) (h_log((1.0+x)/(1.0-x)))*/


#define h_x(LA) ( h_log( ( 1.0+h_exp(-LA[0]) ) / (h_exp(-LA[1])+h_exp(-LA[2])) ) )
#define h_y(LA) ( h_log( ( 1.0+h_exp(-LA[1]) ) / (h_exp(-LA[0])+h_exp(-LA[2])) ) )
#define h_z(LA) ( h_log( ( 1.0+h_exp(-LA[2]) ) / (h_exp(-LA[0])+h_exp(-LA[1])) ) )
#define h_xyz(LA,type) ( \
  ((type)==0? ( h_log( ( 1.0+h_exp(-LA[0]) ) / (h_exp(-LA[1])+h_exp(-LA[2])) ) ) : \
              ((type)==1? ( h_log( ( 1.0+h_exp(-LA[1]) ) / (h_exp(-LA[0])+h_exp(-LA[2])) ) ) : \
                          ( h_log( ( 1.0+h_exp(-LA[2]) ) / (h_exp(-LA[0])+h_exp(-LA[1])) ) ) \
              ) \
  ) \
)
//-- the function (originally called h) lambda_type(Lambda_n) for q type (0,1,2) <-> (X,Y,Z)
/*static double h_lambda(double *LA, int32_t type)
{
  switch(type){
  case 0:   return( h_log( ( 1.0+h_exp(-LA[0]) ) / (h_exp(-LA[1])+h_exp(-LA[2])) ) );   break;
  case 1:   return( h_log( ( 1.0+h_exp(-LA[1]) ) / (h_exp(-LA[0])+h_exp(-LA[2])) ) );   break;
  case 2:   return( h_log( ( 1.0+h_exp(-LA[2]) ) / (h_exp(-LA[0])+h_exp(-LA[1])) ) );   break;
  //default:  printf("NG q type = %u , must debug, suggest CTRL-C to exit;\n", type); getchar();
  }
  printf("NG q type = %u , must debug, suggest CTRL-C to exit;\n", type); getchar();
  return 0.0;
}*/

//#define bplus(a,b) (    h_log( ( 1.0 + (h_tanh2(a) * h_tanh2(b)) ) / ( 1.0 - (h_tanh2(a) * h_tanh2(b)) ) )    )
#define bplus(a,b)  (  h_atanh2( h_tanh2(a) * h_tanh2(b) )  )
#define bminus(a,b) (  h_atanh2( h_tanh2(a) / h_tanh2(b) )  )
/*static double boxplus(double a, double b)
{
  double x;
  //x = tanh(a/2) * tanh(b/2);
  x = h_tanh2(a) * h_tanh2(b);

  //y = 0.5 * h_log((1.0+x)/(1.0-x));   // y = atanh(x)

  //y = h_log((1.0+x)/(1.0-x));         // y = 2.0*atanh(x)
  //return( h_log((1.0+x)/(1.0-x)) );
  return( h_atanh2(x) );
}*/




// initLLR: init LLR vector over {X,Z,Y}
// p_ch: channel err rate (dep ch: eps = p_ch/3;  inde X-Z ch: pX = pZ = p_ch)
static void GenInitLLR(double *initLLR, double p_ch)  // also apply CH_TYPE and ALFA_V insider
{
  if(Q!=4)  { printf("Q=%u , but must =4 ! exit() !!", Q);  getchar();  exit(1); }

#if 0 //AFN
#endif

  //-- init_p
//#if CH_TYPE==0    // dep ch

//#elif CH_TYPE==1  // indep X-Z ch
  //initLLR[0] = h_log( (1.0-p_ch) / (p_ch/3) );
  initLLR[0] = h_log( (1.0-p_ch) *3/p_ch );
  if(fabs(initLLR[0]) > LLR_MAX) { printf("seems too-large magnitude with init Lambda %g , continue ??? \n", initLLR[0]); getchar(); }
  initLLR[1] = initLLR[0];
  initLLR[2] = initLLR[0];
//#else
  //XXX_NG_CH_TYPE_XXX
//#endif

  /*{ int32_t q;
    printf("ideal_LA[0 to Q-1] = ");
    for(q=0;q<Q-1;q++) { printf(" %g ", log( (1.0-p_ch) / (p_ch/3) ) ); }
    printf("\n");
    printf("init_LLR[0 to Q-1] = ");
    for(q=0;q<Q-1;q++) { printf(" %g ", initLLR[q]); }
    printf("\n"); getchar();
  }*/

  //#if ALFA_V == 0   // DEFAULT MODE
  //#else   // enable ALFA_V
  //#endif
}



//-- LLR-based Quantum BP init: set bp->tt as all zeros. return 1 if zz pass; otherwise INIT *bp and return 0
// copy and modified from Qbp_init20()
uint8_t Lbp_init(LBP_Ctl *bp , GFQ_t *zz , a_matrix_GFQ *A , double p_ch, double afp) // zz: syndrome, p_ch: channel err rate, afp: see AFP
{
  int32_t  m , n , j, q;
  //double LA[N][Q-1];

  //-- always reset output bp->tt as all zeros
  memset(bp->tt , 0 , sizeof(bp->tt));

  //-- return 1 if zz pass
  if( is_zero_vec(zz, M) )  return 1; 
  //-- AFP init
  bp->ap = afp;

  //-- generate initial bias vector 
  double initLLR[Q-1];  // init Lambda vec over {X,Z,Y}
  GenInitLLR(initLLR, p_ch);  // also apply CH_TYPE and ALFA_V insider

  //-- initial LA
  for(n=0; n < A->NN; n++)
  for(q=0;q<Q-1;q++) { bp->LA[n][q] = initLLR[q]; }

  /*for(n=0; n < A->NN; n++) {
    printf("bp->LA[%d][0 to Q-1] = ", n);
    for(q=0;q<Q-1;q++) { printf(" %g ", bp->LA[n][q]); }
    printf("\n");
  } getchar();*/


  //-- initial GA   // can skip, as every iteration will init again
#if 0 //ITER_LOG  // may also enable if PARTIAL_UPD_BY_qn or DYNAMIC_pn
  memcpy ( bp->GA , bp->LA, sizeof(bp->GA) );
  ChkNode_GenSyndrByLLR(A, bp);
#endif

  //-- init Lmj[m][j] , according to the edge = X, Z, or Y
  for(m=0; m < A->MM; m++)
  for(j=0; j < A->num_m[m]; j++) {
    n = A->mlist[m][j];
    q = A->mlist_val[m][j];       // {X,Z,Y} val {1,2,3}
    bp->Lmj[m][j] = h_xyz(bp->LA[n], q-1); // q type = q-1
  #if 1 // LLR_MAX and LLR_MIN
    restrict_llr(bp->Lmj[m][j]);
  #endif
  }


  /*double idealLLR, idealGA;
  idealLLR = log( (1.0-p_ch) / (p_ch/3) );    idealGA = log( (1.0+exp(-idealLLR)) / (exp(-idealLLR)+exp(-idealLLR)) );
  printf("\n idealLLR = %g , \n idealGA  = %g \n\n", idealLLR, idealGA);
  for(m=0; m < A->MM; m++) {  // print along m
    printf("bp->Lmj[%d][1 to M(n)] pairs (q, GA) =", m);
    for(j=0; j < A->num_m[m]; j++) {
      //n = A->mlist[m][j];
      q = A->mlist_val[m][j];       // {X,Z,Y} val {1,2,3}
      printf(" ( %d , %g ) ,", q, bp->Lmj[m][j]); // q type = q-1
    }
    printf("\n");
  } getchar();
  int32_t  i;
  for(n=0; n < A->NN; n++) {  // print along n
    printf("bp->Lmj[1 to M(n)][%d] pairs (q, GA) =", n);
    for(i=0; i < A->num_n[n]; i++) {
      m = A->nlist[n][i];
      j = A->ni2j[n][i];
      q = A->nlist_val[n][i];       // {X,Z,Y} val {1,2,3}
      printf(" ( %d , %g ) ,", q, bp->Lmj[m][j]); // q type = q-1
    }
    printf("\n");
  } getchar();*/

  //-- test function boxplus (MACRO bplus) OK!!
  /*double a = -3.45377 , b = 5.13937;
  //printf("boxplus(%g, %g) = %g \n", a, b, boxplus(a, b));
  printf("bplus(%g, %g) = %g \n", a, b, bplus(a, b));
  double ta = (exp(a) - 1.0) / (exp(a) + 1.0);   double tb = (exp(b) - 1.0) / (exp(b) + 1.0);
  double x = ta*tb ,  y = log((1.0+x)/(1.0-x));
  //double y = 2.0 * atanh(x);
  printf("ta = %g , tb = %g ,  x= %g  ideal_bplus(%g, %g) = %g \n", ta, tb, x, a, b, y);
  printf("bminus(%g, %g) = %g \n", y, a, bminus(y, a));
  getchar();*/


#if 0 //EN_AFP_V  // init dqml_ori without applying ALFA_V
#endif


  #if 0 //RND_SCHE == 2
  int32_t n1, tmp_n;   uint64_t seq_seed[4];   rnd256_init_priv(seq_seed);
  for(n=0; n<N-1; n++) {                  //-- cite https://danluu.com/sattolo/
    //n1 = priv_Nminus1(N, seq_seed);       // a discrete r.v. over {0, 1, 2, ..., N-1}   // Fisher-Yates approach (may multi-cycles)
    n1 = priv_IntAB(n+1, N-1, seq_seed);  // a discrete r.v. over {n+1,n+2, ..., N-1}   // Sattolo approach  (guaranteed one cycle)
    tmp_n = bp->seq[n];   bp->seq[n] = bp->seq[n1];   bp->seq[n1] = tmp_n;  // swap idxs at [n] and [n1]
  }
  #endif

  //-- special serial schedule
  //init_votes(bp , zz , A);
  #if 0 //MAJ_VOTE
  #endif
  #if 0 //MAJ_SCHE
  #endif

  //-- set target syndrome
  memcpy ( bp->target_z , zz, sizeof(bp->target_z) );
  //-- INIT *bp done and return 0
  bp->reset = 1;  bp->iter = 0;   //bp->is_hit_c = 0;
  return 0;
}


#define TRK_INF_NAN 0   // also track if the resultant Dm is 0

//-- Parallel LLR-based BP4
uint8_t Lbp_dec60(LBP_Ctl *bp, a_matrix_GFQ *A)   // copy form Qbp_dec20 and modify
{
  int32_t  m, n, i, j, q;
  double Dmj;

  bp->reset = 0;

  //======== Horizontal Step ========//
  for(m=0; m < A->MM; m++) {

    /*bp->Dm[m] = bp->Lmj[m][0];  //printf("Dm %g \n", bp->Dm[m]);  // ASSUME no zero-weight row
    for(j=1; j < A->num_m[m]; j++) { bp->Dm[m] = bplus(bp->Dm[m], bp->Lmj[m][j]); }  // (Dm ready after here)
    //for(j=1; j < A->num_m[m]; j++) { printf("Dm %g \n", bp->Dm[m]); }
    if(bp->target_z[m]) bp->Dm[m] = -bp->Dm[m];  // zm 0 -> 1,  zm 1 -> -1;    (now Dm ready)*/

    bp->Dm[m] = (bp->target_z[m]==0)? bp->Lmj[m][0] : -bp->Lmj[m][0];   // ASSUME no zero-weight row; also handle z_m here
    for(j=1; j < A->num_m[m]; j++) { bp->Dm[m] = bplus(bp->Dm[m], bp->Lmj[m][j]); }  // (Dm ready after here)

    #if TRK_INF_NAN
    if(bp->Dm[m]==INFINITY || bp->Dm[m]==-INFINITY || isnan(bp->Dm[m]) || bp->Dm[m]==0)  printf("\n\n hit bp->Dm[m] = %g  \n\n", bp->Dm[m]);
    #endif

    /*if(isnan(atanh(-1.1)))        printf("\n  find atanh(-1.1) = NAN (%g) \n", atanh(-1.1));
    if(atanh(-1.0) == -INFINITY)  printf("\n  find atanh(-1.0) = -INFINITY (%g) \n", atanh(-1.0));
    if(atanh(1.0) == INFINITY)    printf("\n  find atanh(1.0) = INFINITY (%g)\n", atanh(1.0));
    if(isnan(atanh(1.1)))         printf("\n  find atanh(-1.1) = NAN (%g)\n", atanh(1.1));
    printf("\n\n  atanh(-1.1) = %g , atanh(-1.0) = %g , atanh(-0.9999) = %g , atanh(0.0) = %g , atanh(-0.0) = %g , atanh(0.9999) = %g , atanh(1.0) = %g , atanh(1.1) = %g \n\n", 
      atanh(-1.1), atanh(-1.0), atanh(-0.9999), atanh(0.0), atanh(-0.0), atanh(0.9999), atanh(1.0), atanh(1.1));  //getchar();
    printf("\n\n  atanh(NAN) = %g , atanh(-INFINITY) = %g , atanh(INFINITY) = %g , atanh(NAN) = %g \n\n", 
      atanh(atanh(-1.1)), atanh(atanh(-1.0)), atanh(atanh(1.0)), atanh(atanh(1.1)));  getchar();*/

    for(j=0; j < A->num_m[m]; j++) {  // update Dmj
    #if 1 //EN_AFP_V==0
      bp->Dmj[m][j] = bminus(bp->Dm[m], bp->Lmj[m][j]);
      //bp->Dmj[m][j] *= ((double)bp->target_z[m]*(-2) + 1);  // zm 0 -> 1,  zm 1 -> -1;
    #else
    #endif
      //restrict_d2(...); 
    }
  }

  /*printf("iter = %d \n", bp->iter);
  for(m=0; m < A->MM; m++){
    for(j=0; j < A->num_m[m]; j++) {
      //n = A->mlist[m][j];
      q = A->mlist_val[m][j];       // {X,Z,Y} val {1,2,3}
      printf(" ( %d , %g ) ,", q, bp->Lmj[m][j]); // q type = q-1
    }
    printf(" = bp->Lmj[%d][1 to M(n)] pairs (q, GA)\n", m);
  }
  for(m=0; m < A->MM; m++) { printf(" %g ", bp->Dm[m]); }
  printf(" = bp->Dm \n");
  for(m=0; m < A->MM; m++){
    printf("bp->Dmj[%d][1 to M(n)] pairs (q, GA) =", m);
    for(j=0; j < A->num_m[m]; j++) {
      //n = A->mlist[m][j];
      q = A->mlist_val[m][j];       // {X,Z,Y} val {1,2,3}
      printf(" ( %d , %g ) ,", q, bp->Dmj[m][j]); // q type = q-1
    }
    printf("\n");
  } getchar();*/


  //======== Vertical Step (update GA[n][q] first , update Lmj[m][j] later) ========//
  memcpy ( bp->GA , bp->LA, sizeof(bp->GA) );  // init GA as LA

  //-- loop all variable nodes 
  for(n=0; n < A->NN; n++) {
    #if 0 //MAJ_SCHE
    #endif

    for(i=0; i < A->num_n[n]; i++) {  // update GA  
      m = A->nlist[n][i];
      j = A->ni2j[n][i];

      #if 0 //CYC_CHK
      // buffer Dmj before change
      #endif
      #if ALFA
      bp->Dmj[m][j] *= (100.0/ALFA);
      #endif

      Dmj = bp->Dmj[m][j];   // shorthand notation

      //-- apply AFP
      #if 1 //GMN==0
      if(bp->ap!=0){ Dmj *= bp->ap; }
      #else
      #endif

      #if 0 //BETA
      #endif
      #if 0 //CYC_CHK
      //buffer suppressed Dmj
      #endif

      q = A->mlist_val[m][j]; // {X,Z,Y} val {1,2,3}
      switch (q) {
      case 1: bp->GA[n][1]+=Dmj;  bp->GA[n][2]+=Dmj;  break;  // edge is X: q=1, type=0
      case 2: bp->GA[n][0]+=Dmj;  bp->GA[n][2]+=Dmj;  break;  // edge is Z: q=2, type=1
      case 3: bp->GA[n][0]+=Dmj;  bp->GA[n][1]+=Dmj;  break;  // edge is Y: q=3, type=2
      default:  printf("NG q = %u , must debug, suggest CTRL-C to exit;\n", q); getchar();
      }
    }
  }
  //======== Vertical Step (update GA done) ========//
  /*int32_t b;
      printf("dec iter %d \n", bp->iter);
      for ( b = 0 ; b < N ; b ++ ) {
        for ( q = 0 ; q < Q-1 ; q ++ ) { printf("%g ", bp->GA[b][q]); }   printf(" = bp->La[%d] \n", b);
      }
      getchar();  // STOP per ITERATION*/

  //-------- hard decision --------//
  bp->iter += 1;  // always do this first
  HardDecisionLLR(bp, A);

  //-------- Syndrome Checking --------//
  Quan_GenSyndrome(A, bp->tt, bp->zz);
  #if 0 //STOP_BY_C
  #endif

  //-- check if meet target_z 
  if(0 == HamDist(bp->zz, bp->target_z, A->MM)) { return 1; }
  else if(bp->iter == MAX_ITER)                 { return 0; }
  else                                          { } // DELAY_RETURN_0 since need to update Lmj[m][j]

  //======== Vertical Step (update Lmj[m][j] here) ========//
  for(n=0; n < A->NN; n++) {
    for(i=0; i < A->num_n[n]; i++) {  // update Lmj
      m = A->nlist[n][i];
      j = A->ni2j[n][i];
      q = A->mlist_val[m][j]; // {X,Z,Y} val {1,2,3}
      #if 0 //AFPADD
      #endif

      #if CYC_CHK==0  // CODE_DEFAULT
      #else // CYC_CHK
      #endif

      #if 0 //EN_AFP_V
      #endif
      #if 0 //ALFA_V
      #endif

      //- update Lmj // recall that // q type = q-1
      bp->Lmj[m][j] = h_xyz(bp->GA[n], q-1) - bp->Dmj[m][j];

      #if 1 // LLR MAX and LLR_MIN
      restrict_llr(bp->Lmj[m][j]);
      #endif
    }
    #if 0 //MAJ_SCHE
    #endif
  }
  #if 0 //MAJ_SCHE
  #endif

  return 0;   // do the DELAY_RETURN_0 now
}



//-- Serial LLR-based BP4
uint8_t Lbp_dec64(LBP_Ctl *bp, a_matrix_GFQ *A)   // copy from Qbp_dec24 and modify
{
  int32_t  m, n, i, j, q;//, n2;
  double  Dmj;

  #if RND_SCHE == 3
  #endif

 if(bp->reset == 1) {   //======== make Dm ready first ========//
  for(m=0; m < A->MM; m++) {
    bp->Dm[m] = (bp->target_z[m]==0)? bp->Lmj[m][0] : -bp->Lmj[m][0];   // ASSUME no zero-weight row; also handle z_m here
    for(j=1; j < A->num_m[m]; j++) { bp->Dm[m] = bplus(bp->Dm[m], bp->Lmj[m][j]); }  // (Dm ready after here)
    #if TRK_INF_NAN
    if(bp->Dm[m]==INFINITY || bp->Dm[m]==-INFINITY || isnan(bp->Dm[m]) || bp->Dm[m]==0)  printf("\n\n hit bp->Dm[m] = %g  \n\n", bp->Dm[m]);
    #endif
  }
 }
  bp->reset = 0;

  memcpy ( bp->GA , bp->LA, sizeof(bp->GA) );  // init GA as LA
  //memset(bp->GA, 0, sizeof(bp->GA));  // init GA as 0 (add LA later)  

  //======== main process of this iteration ========//
  for(n=0; n < A->NN; n++) {    // ORI dec64
  //for(n2=0; n2 < A->NN; n2++) { // FOR dec65 (also support dec64)
  //  n = bp->seq[n2];            // FOR dec65 (also support dec64)

    #if 0 //MAJ_SCHE
    #endif

    //-------- Combine from m \in \sM(n) --------//
    for(i=0; i < A->num_n[n]; i++) {
      m = A->nlist[n][i];
      j = A->ni2j[n][i];

      #if 1 //EN_AFP_V==0
      bp->Dmj[m][j] = bminus(bp->Dm[m], bp->Lmj[m][j]);
      #else
      #endif

      #if 0 //CYC_CHK
      // buffer rml before change
      #endif
      #if ALFA
      bp->Dmj[m][j] *= (100.0/ALFA);                    // used for "init GA as LA"
      #endif

      Dmj = bp->Dmj[m][j];   // shorthand notation

      //-- apply AFP
      #if 1 //GMN==0
      if(bp->ap!=0){ Dmj *= bp->ap; }                    // used for "init GA as LA"
      #else
      #endif

      #if 0 //BETA
      #endif
      #if 0 //CYC_CHK
      //buffer suppressed Dmj
      #endif

      q = A->mlist_val[m][j]; // {X,Z,Y} val {1,2,3}
      switch (q) {
      case 1: bp->GA[n][1]+=Dmj;  bp->GA[n][2]+=Dmj;  break;  // edge is X: q=1, type=0
      case 2: bp->GA[n][0]+=Dmj;  bp->GA[n][2]+=Dmj;  break;  // edge is Z: q=2, type=1
      case 3: bp->GA[n][0]+=Dmj;  bp->GA[n][1]+=Dmj;  break;  // edge is Y: q=3, type=2
      default:  printf("NG q = %u , must debug, suggest CTRL-C to exit;\n", q); getchar();
      }
    }
    /*for(q=0; q<Q-1; q++) {  // GA has been (sum of Dmj) ; now to let GA = LA + (sum of Dmj)  // used for "init GA as 0"
      #if ALFA
      bp->GA[n][q] = bp->LA[n][q] + bp->GA[n][q]*(100.0/ALFA);
      #elif AFP
      bp->GA[n][q] = bp->LA[n][q] + bp->GA[n][q]*(bp->ap);
      #else
      bp->GA[n][q] = bp->LA[n][q] + bp->GA[n][q];
      #endif
    }*/


    //-------- Dispatch to m \in \sM(n) --------//
    for(i=0; i < A->num_n[n]; i++) {
      m = A->nlist[n][i];
      j = A->ni2j[n][i];
      q = A->mlist_val[m][j]; // {X,Z,Y} val {1,2,3}
      #if 0 //AFPADD
      #endif

      #if CYC_CHK==0  // CODE_DEFAULT
      #else // CYC_CHK
      #endif

      #if 0 //EN_AFP_V
      #endif
      #if 0 //ALFA_V
      #endif

      //- update Lmj // recall that // q type = q-1
      bp->Lmj[m][j] = h_xyz(bp->GA[n], q-1) - bp->Dmj[m][j];

      #if 1 // LLR MAX and LLR_MIN
      restrict_llr(bp->Lmj[m][j]);
      #endif


      //- make dm more updated (KEY STEP for Serial BP)
      bp->Dm[m] = bplus(bp->Dmj[m][j], bp->Lmj[m][j]);
      #if TRK_INF_NAN
      if(bp->Dm[m]==INFINITY || bp->Dm[m]==-INFINITY || isnan(bp->Dm[m]) || bp->Dm[m]==0)  printf("\n\n hit bp->Dm[m] = %g  \n\n", bp->Dm[m]);
      #endif
    }

    #if 0 //MAJ_SCHE
    #endif
  }
  #if 0 //MAJ_SCHE
  #endif


  //======== hard decision ========//
  bp->iter += 1;  // always do this first
  HardDecisionLLR(bp, A);

  //======== Syndrome Checking ========//
  Quan_GenSyndrome(A, bp->tt, bp->zz);
  #if 0 //STOP_BY_C
  #endif

  //-- check if meet target_z 
  if( 0 == HamDist(bp->zz, bp->target_z, A->MM) )  return 1;
  else                                             return 0;
}



//-- Serial LLR-based BP4 (along checks)
uint8_t Lbp_dec84(LBP_Ctl *bp, a_matrix_GFQ *A)   // copy form Qbp_dec44 and modify
{
  int32_t  m, n, j, q; //i;
  double  Dmj;

  if(bp->reset == 1) { // init GA and default Dmj
    memcpy ( bp->GA , bp->LA, sizeof(bp->GA) );  // init GA as LA
    for(m=0; m < A->MM; m++) for(j=0; j < A->num_m[m]; j++) { bp->Dmj[m][j] = 0.0; }
  }
  bp->reset = 0;

  //-- refresh GA
  // skip first, as need more memories to handle AFP


  //================ serial update along check nodes ================//
  for(m=0; m < A->MM; m++) {
    //-------- pull Lmj and update Dm --------// (Dmj already ready)
    for(j=0; j < A->num_m[m]; j++) {
      Dmj = bp->Dmj[m][j];   // shorthand notation
    #if 0 //AFPADD
    #endif
      n = A->mlist[m][j];
      q = A->mlist_val[m][j]; // {X,Z,Y} val {1,2,3}
      switch (q) {
      case 1: bp->GA[n][1]-=Dmj;  bp->GA[n][2]-=Dmj;  break;  // edge is X: q=1, type=0
      case 2: bp->GA[n][0]-=Dmj;  bp->GA[n][2]-=Dmj;  break;  // edge is Z: q=2, type=1
      case 3: bp->GA[n][0]-=Dmj;  bp->GA[n][1]-=Dmj;  break;  // edge is Y: q=3, type=2
      default:  printf("NG q = %u , must debug, suggest CTRL-C to exit;\n", q); getchar();
      }

    #if 0 //EN_AFP_V
    #endif
    #if 0 //ALFA_V
    #endif

      //- update Lmj // recall that // q type = q-1
      bp->Lmj[m][j] = h_xyz(bp->GA[n], q-1); //- bp->Dmj[m][j];  // -Dmj has been done above

      #if 1 // LLR MAX and LLR_MIN
      restrict_llr(bp->Lmj[m][j]);
      #endif

      #if 0 //STOP_BY_C
      #endif

      //////// UPDATING dm HERE !!!! ////////////////
      if(j==0)  bp->Dm[m] = (bp->target_z[m]==0)? bp->Lmj[m][0] : -bp->Lmj[m][0];   // ASSUME no zero-weight row; also handle z_m here
      else      bp->Dm[m] = bplus(bp->Dm[m], bp->Lmj[m][j]);
    }
    #if TRK_INF_NAN
    if(bp->Dm[m]==INFINITY || bp->Dm[m]==-INFINITY || isnan(bp->Dm[m]) || bp->Dm[m]==0)  printf("\n\n hit bp->Dm[m] = %g  \n\n", bp->Dm[m]);
    #endif

    //-------- dispatch Dmj --------//  //// to make GA more updated
    for(j=0; j < A->num_m[m]; j++) {
      #if 1 //EN_AFP_V==0
      bp->Dmj[m][j] = bminus(bp->Dm[m], bp->Lmj[m][j]);
      #else
      #endif

      #if 0 //CYC_CHK
      // buffer Dmj before change
      #endif
      #if ALFA
      bp->Dmj[m][j] *= (100.0/ALFA);
      #endif

      Dmj = bp->Dmj[m][j];   // shorthand notation

      //-- apply AFP
      #if 1 //GMN==0
      if(bp->ap!=0){ Dmj *= bp->ap; }
      #else
      #endif

      #if 0 //BETA
      #endif
      #if 0 //CYC_CHK
      //buffer suppressed Dmj
      #endif

      n = A->mlist[m][j];
      q = A->mlist_val[m][j]; // {X,Z,Y} val {1,2,3}
      switch (q) {
      case 1: bp->GA[n][1]+=Dmj;  bp->GA[n][2]+=Dmj;  break;  // edge is X: q=1, type=0
      case 2: bp->GA[n][0]+=Dmj;  bp->GA[n][2]+=Dmj;  break;  // edge is Z: q=2, type=1
      case 3: bp->GA[n][0]+=Dmj;  bp->GA[n][1]+=Dmj;  break;  // edge is Y: q=3, type=2
      default:  printf("NG q = %u , must debug, suggest CTRL-C to exit;\n", q); getchar();
      }
    } //-------- dispatch to each dml --------// DONE
  } //======== END OF "serial update along check nodes" ========//


  //======== hard decision ========//
  bp->iter += 1;  // always do this first
  HardDecisionLLR(bp, A);

  //======== Syndrome Checking ========//
  Quan_GenSyndrome(A, bp->tt, bp->zz);
  #if 0 //STOP_BY_C
  #endif

  //-- check if meet target_z 
  if( 0 == HamDist(bp->zz, bp->target_z, A->MM) )  return 1;
  else                                             return 0;

}


