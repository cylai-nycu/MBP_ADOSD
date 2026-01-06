#include <math.h>
#include "./ldpc_parm.h"
//#include "./lib_mat/lib_mat.h"
#include "./lib_rand/lib_rand.h"
#include "./bp_dec/bp_dec.h"
#include "./bp_dec/bp_llr.h"


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
    sum=0; for(q=0;q<Q;q++) { sum += bp->qn[b][q]; }
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
int main0(void)
{
//double p_nos = pow(10.0, (-10.0  -  10.0/8 *10)/10);   //#define P_START (-10.0  -  10.0/8 *10)    // 0.0056
//double p_nos = (double)0.013;
//double p_nos = (double)0.01;

double p_nos = (double)0.1;
//double p_nos = (double)0.2;
//double p_nos = (double)0.1;     // BSC error probability
//double p_nos = pow(10.0, -14.9/10);                   // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = pow(10.0, -14.99/10);    // [[5,1,3]] code, +1 redundant check row (hit Degenerate Err)
//double p_nos = pow(10.0, -15.0/10);     // [[5,1,3]] code, +1 redundant check row (hit Degenerate Err)
//double p_nos = pow(10.0, -15.01/10);                   // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
//double p_nos = (double)0.031;                         // [[5,1,3]] code, +1 redundant check row (hit Exact Err)
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
  int32_t b , q; // dummy variables

  //- RX
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
  rr[3] = 3;  // set 4th qubit error be Y       // make [[5,1,3]] Qbp_dec20 FAIL
  //rr[4] = 1;  // set 5th qubit error be X
  //rr[4] = 2;  // set 5th qubit error be Z       // make [[5,1,3]] Qbp_dec20 FAIL when p=0.0031 , a=1.6; and dec_44 FAIL
  //rr[4] = 3;  // set 5th qubit error be Y       // make [[5,1,3]] Qbp_dec20 FAIL
#elif ORI==202
  rr[0] = 1;  // set first qubit error be X
#elif ORI==2022
  rr[0] = 1;  // set first qubit error be X
  rr[0+N/2] = 1;  // set first qubit error be X   // to TEST [[2,0,2]] L=2
#elif ORI==513 || ORI==5131
  rr[3] = 3;  // set 4th qubit error be Y       // make [[5,1,3]] Qbp_dec20 FAIL
#elif ORI==5132
  rr[3] = 3;  // set 4th qubit error be Y       // make [[5,1,3]] Qbp_dec20 FAIL
  rr[3+N/2] = 3;  // set 4th qubit error be Y       // to TEST [[5,1,3]] L=2
#elif ORI==1115
  //rr[0]=2; rr[4]=2; rr[6]=2;  // s=81: 	1 0 0 0 1 0 1 0 0 0
  //rr[0]=2; rr[2]=2; rr[8]=2;  // s=261: 1 0 1 0 0 0 0 0 1 0
  //rr[0]=2; rr[1]=2; rr[3]=2;  // s=11:  1 0 1 1 0 0 0 0 0 0
  rr[0]=2; rr[1]=2; rr[2]=2;  // s=7:  1 1 1 0 0 0 0 0 0 0
#elif ORI==4115   // Surface Codes
  //-- for L=4 @ p_nos = pow(10.0, -2.0 -1.0/8);
  //rr[0+36]=1;  rr[1+36]=1;   //rr[0+36]=3;  rr[1+36]=3;   // +36 NG
  //rr[0+4]=2;  rr[9+4]=2;   // +3,+4 NG
  //rr[0+31]=1;  rr[9+31]=3;   // +31 NG  (hit max iter)
  //rr[0+35]=3;  rr[5+35]=2;   // +35 NG  (hit max iter)
  //rr[0+35]=1;  rr[5+35]=2;   // +35 NG  (hit max iter)
  //rr[35]=1;  rr[39]=2;   // NG if rr[35,39] = [1,2] or [3,2];
  //-- for L=4 @ p_nos = pow(10.0, -2.0 -3.0/8);
  rr[5]=3;  rr[10]=3;

#elif ORI==2515   // Surface Code ROTATED
  #if Surf_d == 9
  //rr[1]=3;  rr[2]=1;                      // make L=9 NG if not fixP 0.013  @ #define P_START (-10.0  -  10.0/8 *10)    // 0.0056
  //rr[0]=3;  rr[1]=3;  rr[2]=3;  rr[3]=3;  // make L=9 NG if use dec20       @ #define P_START (-10.0  -  10.0/8 *10)    // 0.0056
  //rr[5]=1;  rr[14]=3;  rr[55]=3;  rr[71]=2;  rr[75]=3;    //autoP: will flip bits 4 5 14 24 45 55 62 71 75 (fixP flips 6 bits in these); dec20 autoP only flips: 14 55 65 75 (and fixed at 14 55 75)
  //rr[25]=1; rr[33]=1; rr[61]=1; rr[62]=1; rr[64]=2; rr[72]=2; rr[73]=3;   // make L=9 NG if use dec20   @ #define P_START (-20.0)   // 0.01
  //rr[22]=2; rr[29]=2; rr[30]=2; rr[51]=1; rr[58]=3; rr[70]=2; rr[71]=2;   // make L=9 NG if use dec20   @ #define P_START (-20.0)   // 0.01
  //rr[3]=1; rr[23]=3; rr[34]=1; rr[60]=2; rr[69]=2; rr[70]=1;  // make L=9 NG if use dec20   @ #define P_START (-20.0)   // 0.01
  //rr[9]=2; rr[22]=3; rr[32]=2; rr[60]=1; rr[72]=3; rr[77]=1;  // make L=9 NG if use dec20   @ #define P_START (-20.0)   // 0.01
  //rr[16]=3; rr[27]=2; rr[43]=1; rr[50]=3; rr[69]=2; rr[74]=1; // make L=9 NG if use dec20   @ #define P_START (-20.0)   // 0.01
  //rr[24]=2; rr[28]=1; rr[32]=3; rr[43]=1; rr[62]=2; rr[76]=1;  // make L=9 NG if use dec20   @ #define P_START (-20.0)   // 0.01

  //-- L=9 manual cases @ P_START 0.01                                    // behaviors @ //#define P_START (-10.0  -  10.0/8 *10)    // 0.0056
  //rr[63-1]=1; rr[62-1]=2; rr[61-1]=3; rr[60-1]=2; rr[68-1]=3; rr[77-1]=1; // CASE 1-1: dec20 NG (trap in 6 errors), dec24 ap65 ok at iter 23
  //rr[63-1]=2; rr[62-1]=2; rr[61-1]=3; rr[60-1]=2; rr[68-1]=3; rr[77-1]=1; // CASE 1-2: dec20 NG (trap in 6 errors), dec24 ap65 ok at iter 38

  //rr[19-1]=2; rr[20-1]=1; rr[29-1]=3; rr[25-1]=3; rr[33-1]=3; rr[5-1]=1;  // CASE 2-1: dec20 NG (trap in 6 errors), dec24 ap65 ok at iter 4
  //rr[19-1]=2; rr[20-1]=1; rr[29-1]=3; rr[25-1]=2; rr[33-1]=3; rr[5-1]=1;  // CASE 2-2: dec20 NG (trap in 4 errors), dec24 ap65 ok at iter 2

  //rr[19-1]=2; rr[20-1]=2; rr[29-1]=3; rr[25-1]=2; rr[33-1]=3; rr[5-1]=1;  // CASE 2-3: dec20 NG (trap in 3 errors), dec24 ap65 ok at iter 2
  //rr[19-1]=2; rr[20-1]=2; rr[29-1]=3; rr[25-1]=3; rr[33-1]=2; rr[5-1]=1;  // CASE 2-4: dec20 NG (trap in 4 errors), dec24 ap65 ok at iter 3

  //rr[19-1]=2; rr[20-1]=2; rr[21-1]=3; rr[14-1]=2; rr[22-1]=3; rr[5-1]=1;  // CASE 3-1: dec20 NG (trap in 10 errors) dec24 ap65 ok at iter 12
  //rr[19-1]=2; rr[20-1]=2; rr[21-1]=3; rr[14-1]=3; rr[22-1]=2; rr[5-1]=1;  // CASE 3-2: dec20 NG (trap in 6 errors), dec24 ap65 ok at iter 20

  //rr[18-1]=3; rr[36-1]=3;

  #elif Surf_d == 7
  //-- L=7 manual cases @ P_START 0.013
  //rr[15-1]=2; rr[16-1]=2; rr[23-1]=3; rr[25-1]=2; rr[19-1]=3; rr[3-1]=1;  // CASE1: dec20 NG , dec24 ap65 ok at iter 5
  //rr[15-1]=2; rr[16-1]=2; rr[23-1]=3; rr[25-1]=3; rr[19-1]=2; rr[3-1]=1;  // CASE2: dec20 NG , dec24 ap65 at iter 6
  //rr[15-1]=2; rr[16-1]=2; rr[23-1]=3; rr[33-1]=2; rr[39-1]=3; rr[3-1]=1;  // CASE3: dec20 NG , dec24 ap65 ok at iter 4
  //rr[15-1]=2; rr[16-1]=2; rr[23-1]=3; rr[25-1]=3; rr[19-1]=2; rr[4-1]=1;  // CASE4: dec20 NG , dec24 ap65 ok at iter 6
  //rr[15-1]=2; rr[16-1]=2; rr[23-1]=3; rr[33-1]=2; rr[39-1]=3; rr[4-1]=1;  // CASE5: dec20 NG , dec24 ap65 ok at iter 4
  //rr[15-1]=2; rr[16-1]=2; rr[23-1]=3; rr[33-1]=2; rr[39-1]=3; rr[40-1]=1; rr[4-1]=1;  // CASE6: dec20 NG , dec24   ap65 ok at iter 4
  rr[15-1]=2; rr[16-1]=2; rr[23-1]=3; rr[33-1]=2; rr[39-1]=3; rr[40-1]=3; rr[4-1]=1;  // Y CASE7: dec20 NG , dec24 ap65 ok at iter 10
  //rr[15-1]=2; rr[16-1]=2; rr[23-1]=3; rr[33-1]=2; rr[39-1]=3; rr[40-1]=3; rr[4-1]=1; rr[19-1]=2; rr[13-1]=3;  // CASE8: dec20 NG , dec24 ap65 ok at iter 13
  //rr[15-1]=2; rr[16-1]=2; rr[23-1]=3; rr[33-1]=2; rr[39-1]=3; rr[40-1]=3; rr[4-1]=1; rr[13-1]=2; rr[14-1]=2;  // Y CASE9:  dec20 NG (final 3 errs) , dec24 ap65 ok at iter  4 
  ////rr[15-1]=2; rr[16-1]=2; rr[23-1]=3; rr[33-1]=2; rr[39-1]=3; rr[40-1]=3; rr[4-1]=1; rr[6-1]=1;  rr[7-1]=1;   // Y CASE10: dec20 NG (final 5 errs) , dec24 ap65 ok at iter 13 
  //rr[15-1]=2; rr[16-1]=2; rr[23-1]=3; rr[33-1]=2; rr[39-1]=3; rr[40-1]=3; rr[4-1]=1; rr[7-1]=1;  rr[14-1]=1;  // Y CASE11: dec20 NG (final 4 errs) , dec24 ap65 ok at iter  8  

  #elif Surf_d == 5
  //-- L=5 manual cases @ P_START 0.013
  //rr[11-1]=2; rr[12-1]=2; rr[17-1]=3; rr[13-1]=2; rr[9-1]=3; rr[5-1]=1;   // dec24 ap65 NG at iter 2 (False Alarm)
  //rr[11-1]=2; rr[12-1]=2; rr[17-1]=3; rr[13-1]=2; rr[9-1]=3; rr[3-1]=1;   // dec24 ap65 NG at iter 2 (False Alarm)

  //rr[11-1]=2; rr[12-1]=2; rr[17-1]=3; rr[23-1]=2; rr[19-1]=3; rr[5-1]=1;  // dec24 ap65 OK at iter 3
  //rr[11-1]=2; rr[12-1]=2; rr[17-1]=3; rr[23-1]=2; rr[19-1]=3; rr[3-1]=1;  // dec24 ap65 NG at iter 15 (False Alarm)
  //rr[11-1]=2; rr[12-1]=2; rr[7-1]=3; rr[23-1]=2; rr[19-1]=3; rr[3-1]=1;  // dec24 ap65 NG at iter 15 (False Alarm)

  //rr[11-1]=2; rr[12-1]=2; rr[17-1]=3; rr[23-1]=1; rr[19-1]=3; rr[5-1]=1;  // dec24 ap65 OK at iter 7
  //rr[11-1]=2; rr[12-1]=2; rr[17-1]=3; rr[23-1]=1; rr[19-1]=3; rr[3-1]=1;  // dec24 ap65 NG at iter 6 (False Alarm)

  //rr[11-1]=2; rr[12-1]=2; rr[17-1]=3; rr[23-1]=1; rr[19-1]=1; rr[5-1]=1;  // dec24 ap65 OK at iter 4
  //rr[11-1]=2; rr[12-1]=2; rr[17-1]=3; rr[23-1]=1; rr[19-1]=1; rr[20-1]=3;  // dec24 ap65 OK at iter 4
  //rr[11-1]=2; rr[12-1]=2; rr[17-1]=3; rr[23-1]=1; rr[19-1]=1; rr[24-1]=3;  // dec24 ap65 NG at iter 3  (False Alarm)

  //rr[7-1]=3; rr[8-1]=3; rr[13-1]=3; rr[18-1]=3; rr[19-1]=3; // NG
  //rr[6-1]=3; rr[8-1]=3; rr[13-1]=3; rr[18-1]=3; rr[20-1]=3; // NG
  //rr[3-1]=2; rr[8-1]=1; rr[18-1]=1; rr[23-1]=2; // NG
  rr[2-1]=3; rr[12-1]=3; rr[13-1]=3; rr[14-1]=3; rr[24-1]=3; // NG (damping) / dec24_ai65 OK (exact) / dec24_ac65 OK (exact) 
  //rr[2-1]=3; rr[12-1]=3; rr[13-1]=1; rr[14-1]=3; rr[24-1]=3; rr[8-1]=1; rr[18-1]=1; // NG
  //rr[12-1]=3; rr[13-1]=3; rr[14-1]=3; rr[8-1]=1; rr[18-1]=1; // NG (damping, then converge wt-2) / dec24_ai65 OK (deg) / dec24_ac65 OK (deg) 

  #else
  //rr[5-1]=3;  // OK
  //rr[1-1]=1; rr[9-1]=1; // NG (false alarm)
  //rr[1-1]=2; rr[9-1]=2; // OK
  //rr[1-1]=3; rr[9-1]=3; // OK
  //rr[1-1]=1; rr[9-1]=2; // NG (converged low-wt)
  //rr[1-1]=1; rr[9-1]=3; // NG (converged low-wt)
  //rr[1-1]=2; rr[9-1]=3; // OK
  rr[2-1]=3; rr[8-1]=3; // NG (damping) / dec20_ai120-870 OK / dec20_ac150-190 OK
  //rr[1-1]=2; rr[9-1]=2; rr[3-1]=1; rr[7-1]=1; // OK
  //rr[1-1]=1; rr[9-1]=1; rr[3-1]=2; rr[7-1]=2; // NG (false alarm)
  //rr[1-1]=3; rr[9-1]=3; rr[3-1]=3; rr[7-1]=3; // NG (damping; then converge to I) / dec24_ai90 OK (deg) / dec24_ac90 OK (deg) 
  //rr[2-1]=2; rr[8-1]=2; rr[4-1]=1; rr[6-1]=1; // the same as previous
  //rr[2-1]=1; rr[8-1]=1; rr[4-1]=2; rr[6-1]=2; // NG (false alarm)
  //rr[2-1]=3; rr[8-1]=3; rr[4-1]=3; rr[6-1]=3; // OK (deg)
  //rr[2-1]=1; rr[8-1]=1; rr[4-1]=1; rr[6-1]=1; // OK (deg)
  //rr[2-1]=2; rr[8-1]=2; rr[4-1]=2; rr[6-1]=2; // OK (deg)
  //rr[1-1]=3; rr[9-1]=3; rr[3-1]=3; rr[7-1]=1; // NG (converge to I)
  //rr[1-1]=3; rr[9-1]=3; rr[3-1]=3; rr[7-1]=2; // NG (damping; then converge to wt-2)
  //rr[1-1]=3; rr[9-1]=3; rr[3-1]=3; rr[7-1]=3; rr[5-1]=3; // NG (false alarm)
  #endif

#endif

  for ( b = 0 ; b < N ; b ++ )  printf("%u ", rr[b]);
  printf(" = added error pattern \n");

  //-- generate syndrome Ar = z
  Quan_GenSyndrome(A, rr, zz);
  if(test_prt) {
    for ( b = 0 ; b < M ; b ++ ) { printf("%d ", zz[b]); }
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
        for ( q = 0 ; q < Q ; q ++ ) { printf("%.2f ", bp->qn[b][q]/sum); }   printf(", ");
      } printf(" = bp->qn (normalized) \n");
      //getchar();  // STOP per ITERATION
    #else
      for ( b = 0 ; b < N ; b ++ ) {
        for ( q = 0 ; q < Q-1 ; q ++ ) { printf("%g ", bp->GA[b][q]); }   printf(" = bp->La[%d] \n", b);
      }
      //getchar();  // STOP per ITERATION
    #endif
    }

    //if ( syndrome_ok || iter==MAX_ITER )   break;     // CODE DEFAULT
    if(iter==MAX_ITER){ break; } else{ syndrome_ok=0; } // TEST FORCE MAX_ITER
  }

  GFQ_t diff[N], zz_G[M_G];
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

  #if ITER_LOG
  fclose(fpLogE);  fclose(fpLogC);      fclose(_iter_log);  fclose(_iter_chk);
  #endif

	return 0;
}
