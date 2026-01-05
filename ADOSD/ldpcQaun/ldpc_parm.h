#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#define USE_GF2_DEC 0


/*
Define the code to be simulated.
*/

//#define ORI 202// [[2,0,2]] code
//#define ORI 2022// [[2,0,2]] code, L=2
//#define ORI 412// [[4,1,2]] code
//#define ORI 422// [[4,2,2]] code
//#define ORI 513// [[5,1,3]] code
//#define ORI 5131// [[5,1,3]] code, +1 redundant check row
//#define ORI 5132// [[5,1,3]] code, L=2

//#define ORI 713// [[7,1,3]] code
//#define ORI 1115// [[11,1,5]] code
//#define ORI 1023// [[10,2,3]] code, Toric Code TWISTED (cyclic based on XZZX)
//#define ORI 1624// [[16,2,4]] code, Toric Code ROTATED  (also support XZZX)

//#define ORI 2515 // [[25,1,5]] code, Surface Code ROTATED (also support XZZX)
//#define Surf_d  13  // d_min (only ODD)


//#define ORI 488   // (4.8.8)  color codes
//#define hex_d 5 //distance


//#define ORI 251580// [[25,1,5]] code, Surface Code ROTATED  (ALL Redundant ROWs wt<=4, e.g., d=5 has M=80)
//#define ORI 13134// [[13,1,3]] code, Surface Code, X (M=N-1), extended 4 redundant rows    (NOT USED)
//#define ORI 4115// [[41,1,5]] code, Surface Code, X (M=N-1)                               (WE CHOOSE THIS)
//#define ORI 1823// [[18,2,3]] code, Toric Code                                            (NOT USED)
//#define ORI 18232// [[18,2,3]] code, Toric Code, 2 redundant rows                          (WE CHOOSE THIS)
//#define ORI 1915// [[19,1,5]] code, Color Code
//#define ORI 7134// [[7,1,3]] (color) code, cyclic shift all rows (4 more redundant rows, with M=7+7 in H)
//#define ORI 11357// [[113,57,12]] code from Grassl's table
//#define ORI 12928// [[129,28]] hypergraph-product code by [7,4,3] and [15,7,5] BCH codes
//#define ORI 25632// [[256,32]] bicycle code
//#define ORI 25634// [[256,34]] bicycle code (check matrix A4 has 2 redundant rows)

//#define ORI 80040016// [[800,400]] GB16 code
//#define ORI 800400// [[800,400]] bicycle code
//#define ORI 800404// [[800,404]] bicycle code (check matrix A4 has 4 redundant rows)
//#define ORI 800418// [[800,418]] bicycle code (check matrix A4 has 18 redundant rows)
//#define ORI 3786// [[3786,946]] bicycle code
//#define ORI 12628// [[126,28]] GB code

//#define ORI 88224// [[882,24]] GHP code
//#define ORI 88248// [[882,48]] GHP code
//#define ORI 882482// [[882,48]] GB code

//#define ORI 1844// color on torus




//#define ORI 36999   // [[369,9,9]] lifted-connected surface code
//#define ORI 17577   // [[175,7,7]] lifted-connected surface code
//#define ORI 6555   // [[65,5,5]] lifted-connected surface code
//#define ORI 1533   // [[15,3,3]] lifted-connected surface code

//#define ORI 625258   // [[625,25,8]] HP code


#define ORI 7212   // [[72,12,6]] BB cod
//#define ORI 908   // [[90,8,10]] BB codee
//#define ORI 1088   // [[108,8,10]] BB code
//#define ORI 14412   // [[144,12,12]] BB code
//#define ORI 28812   // [[288,12,18]] BB code


//#define ORI 1441237   // [[144,12,12]] GB code rho=8  HYL

//#define ORI 500194   // [[500,194]] code
//#define ORI 500188   // [[500,188]] code
//#define ORI 50019465   // [[500,194]] code  65 GenTannerSch
//#define ORI 500196   // [[500,196]] code  GT




//#define ORI 392322  // [[392,32]] QGB code, rw=12, de-shih
//#define ORI 3923272  // [[392,32,7]] code, QGB code, HYL

//#define ORI 3923214  // [[392,32,14]] QGB code,  HYL
//#define ORI 3804610  // [[380,46,10]] QGB code, rw=13, HYL
//#define ORI 3924814  // [[392,48,14]] QGB code


//#define ORI 39232  // [[392,32]] QT code  De-shih
//#define ORI 3923212// [[392,32,12]] code, QT code, HYL  1523
//#define ORI 39248  // [[392,48]] QT code
//#define ORI 39236  // [[392,36]] QT code

//#define ORI 68628  // [[686,28,14]] QGT code 68628_2017
//#define ORI 68634  // [[686,34,13]] QGT code 68634_1306


//#define ORI 68830 // [[688,30,19]] GB code  688_30_16



//#define ORI 4804012  // [[480,40,12]] QGB code, #6 rw= HYL

//#define ORI 144128  // [[144,128,8]] QGB code, 144 checks
//#define ORI 1441282  // [[144,128,8]] QGB code, rw=7
//#define ORI 1441283  // [[144,128]] QGB code, rw=7  de-shih

//#define ORI 5121748  // [[512,174,8]] SPC code
//#define ORI 5121744  // [[512,174]] QGB code #4  d=5
//#define ORI 5121746  // [[512,174]] QGB code #15  d=6
//#define ORI 5121747  // [[512,174]] QGB code  rw 9    QGB512174_21_G4



//#define ORI 51217484  // [[512,174]] QGB code weight 10  de-shih
//#define ORI 51217482  // [[512,174]] QGB code weight 8  de-shih

//#define ORI 51217483  // [[512,174]] QGB code #4  d=6




/*
Define the decoder parameters
*/

#define LLR_BP  0   // 0: off (i.e., use linear BP), 1: use LLR-BP




#define OSDW  (-2)       // -1 no OSD, -2 APOSD, and others
#define RelThr  (0.999995)  // reliability threshold
#define MAX_ITER  (100)
#define AFP   160     // alpha prime, set 0 to off, set to 120 means 120/100  (\alpha in the MBP paper)

/*adjust alpha adaptively in AMBP */
// MBP:set AFP to some number, while set AFP2=0. AFPINC=0
// Larger AFP means smaller step size
#define AFP2  0     // set to non-zero to enable:  if AFPINC!=0 further, then must set a non-zero AFP, e.g. AFP=100
#define AFPINC  (0) // (what AFP2 used:) 0: off; otherwise, afp refers to AFP+AFPINC*n for the additional n-th trial until reach AFP2
                    // (-1)


#define BY_DEC  20  // 20 is parallel dec, 24 is serial dec, 25 is circular schedule (Surface code only, 26 is syndrome-based schedule (not enabled now)
                    //                     44 is serial dec along check nodes
                    // 60 is parallel LLR, 64 is serial LLR
                    //                     84 is serial LLR along check nodes
//random schedule
#define RND_SCHE 0  // 0: off, random schedule 1: init , 2: per CW, 3: per iter



//#define P_STEP  (-0.1)  // err probability step (in dB, every time 10^(STEP/10))
//#define P_STEP  (-0.02)  // err probability step (in dB, every time 10^(STEP/10))
//#define P_STEP  (-10.0/4)  // err probability step (in dB), every time p=p*10^(STEP/10))
// -10.0/n
#define P_STEP  (-10.0/8)  // err probability step (in dB, every time 10^(STEP/10))
#define FIX_P (0)       // set 0 to off, set to 46 means 4.6e-2 = 0.046 if FIX_e = 0 or -2
#define FIX_e (0)     // set 0 to off, set to -2 means 4.6e-2 if FIX_P = 46

#define ALFA  0     // set 0 to off, set to 120 means 100/120  (\alpha_c in refined BP)





  #define DEC2    '0' // (what AFP2 used:) '0': DEC2 follows BY_DEC,  'p': DEC2 by dec20,  's': DEC2 by dec24,  'c': DEC2 by dec44
  #define DBL_CHK 0   // if this afp has the same output as last afp (or this afp hits max iter), then stop


#define AFPADD  0   // alpha prime post added, set 0 to off, set to 120 means to add 100/120
                    // When message normalized by 100/AFP, the memory weight, originally (1 - 100/AFP), can be adjusted here

#define BETA  0     // set 0 to off, this is in linear domain (for efficiency, try to set to something like 2, 4, 8, 16, ...)

#define ALFA_V  0   // alpha_v, set 0 to off, set to 12 means 10/12
#define EN_AFP_V 0  // 0: disable alpha_v prime, 1: enable (and alpha_v controlled by ALFA_V)



#define MAJ_VOTE 0  // 0: off, 1: vote by syndrome (test maj_voting(); threshold j=1,2,3, not improved), 2: vote outer bits less reliable (need ORI=4115) (tested not good)
#define MAJ_SCHE 0  // majority voting schedule. 0: off, 1: on (soft strategy), 2: on (hard strategy), 4 on (hard, 2 times)

#define CH_TYPE 0   // 0: dep ch,  1: indep X-Z

#define POS_AVG 0   // 0 is off, or set POS_AVG>=2: (old*(POS_AVG-1) + new)/POS_AVG


#define AFN   0     // alpha for init_pn, set 0 to off, set to 12 means 10/12

#define CYC_CHK 0   // 0: off, 1: check short cycles (and adjust edge weights)

#define ZE 0        // 0: off; e.g. if ZE=35, then init_p = wt_E/N with an estimated wt_E = wt(z)*ZE/100


#define GMN 0       // 0: off; e.g., if gain GMN=140, then (1/AFP)*=140/100 if |\sN(m)|=2.
#define LMN 0       // default 200; when GMN is enabled, the (1/AFP)*=140/100 will be upper bounded by LMN


#define LIST_DEC 0  // 0: off, 1: buffer the last syndrome-matched solution


#define STOP_BY_C 0   // to simulate the effect: stop at check node when syndrome match (0 Off, 1 ON, 2 TEST force OK)
#define ITER_LOG (1 && STOP_BY_C)   // 0: off, 1: log iteration data (must enable STOP_BY_C)


#define BETA_pn   0   // set 0 to off, set to 120 means 100/120

#define OPTalpha (0)

#define RV_LEN  (64)  // row variable, when doing OSD, we can store RV_LEN bits in a variable
#define POSD_FLIPPING_ALL_BIT (1)
#define OSD_RECORD (1)

# if RV_LEN == 64
    typedef uint64_t RV_t;
    #define QShift 5
    #define BShift 6
    #define QMask 31
    #define BMask  63
# else
    typedef uint32_t RV_t;
    #define QShift 4
    #define BShift 5
    #define QMask 15
    #define BMask  31
# endif // RV_LEN

#define STORE_WRONG 0// store the error that fail solve
#define STORE_LF 0
#define CAL_TIME 1
#define STORE_DISTR 0






/*
Define the parameters for codes and
the source files for the  parity check matrix and normalizer matrix

N: number of qubits
K: number of logical qubits
M: number of parity checks
PATH_A: parity check matrix
PATH_Gs: normalizer matrix
*/

#if ORI==202   // [[2,0,2]] code
  #define N (2)
  #define K (0)
  #define M   (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("Q202_A4")
  #define PATH_Gs ("Q202_G4s")
  #define PRE_A    ("202")



#elif ORI==1441237  // [[144,12,12]] QGB code r=8
  #define N (144)
  #define K (12)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_Gs ("codes/QGB144_12_37_G4.txt")
    #define PATH_A  ("codes/QGB144_12_37_H4.txt")
  #define PRE_A   ("QGB144_12_37")



#elif ORI==144128  // [[144,128,8]] QGB code nonfull rank matrix
  #define N (144)
  #define K (12)
  #define M (144)
  #define M_G   (N+K)
  #define PATH_Gs ("codes/QGB144128_G4.txt")
    #define PATH_A  ("codes/QGB144128_A4.txt")
  #define PRE_A   ("QGB144128")


#elif ORI==1441282  // [[144,128,8]] QGB code rw=7 #24
  #define N (144)
  #define K (12)
  #define M (N-K)
  #define M_G  (N+K)
  #define PATH_Gs ("codes/QGB14412_24_G4.txt")
    #define PATH_A  ("codes/QGB14412_24_A4.txt")
  #define PRE_A   ("QGB144128_24")

#elif ORI==1441283  // [[144,128,8]] QGB code rw=7 de-shih
  #define N (144)
  #define K (12)
  #define M (N-K)
  #define M_G  (N+K)
  #define PATH_Gs ("codes/QGB14412_w7_G4s")
    #define PATH_A  ("codes/QGB14412_w7_A4")
  #define PRE_A   ("QGB144128_w7ds")

#elif ORI==5121744  //  [[512,174]] QGB code #4  d=5
  #define N (512)
  #define K (174)
  #define M (N-K)
  #define M_G  (N+K)
  #define PATH_Gs ("codes/QGB5121744_G4.txt")
    #define PATH_A  ("codes/QGB5121744_A4.txt")
  #define PRE_A   ("QGB5121744")

#elif ORI==5121746  //  [[512,174]] QGB code #15  d=6
  #define N (512)
  #define K (174)
  #define M (N-K)
  #define M_G  (N+K)
  #define PATH_Gs ("codes/QGB5121746_G4.txt")
    #define PATH_A  ("codes/QGB5121746_A4.txt")
  #define PRE_A   ("QGB5121746")


#elif ORI==5121747  //  [[512,174]] QGB code  rw 9    QGB512174_21_G4
  #define N (512)
  #define K (174)
  #define M (N-K)
  #define M_G  (N+K)
  #define PATH_Gs ("codes/QGB512174_21_G4.txt")
    #define PATH_A  ("codes/QGB512174_21_A4.txt")
  #define PRE_A   ("QGB512174_21")




#elif ORI==5121748  //  [[512,174,8]] SPC code
  #define N (512)
  #define K (174)
  #define M (N-K)
  #define M_G  (N+K)
  #define PATH_Gs ("codes/QSPC512174_G4s.txt")
  #define PATH_A  ("codes/QSPC512174_A4.txt")
  #define PRE_A   ("QSPC5121748")

#elif ORI==51217484  //  [[512,174,8]] QGB code weight 10  de-shih
  #define N (512)
  #define K (174)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_Gs ("codes/QGB512174_w10_G4s.txt")
  #define PATH_A  ("codes/QGB512174_w10_A4.txt")
  #define PRE_A   ("QGB512174_w10")


#elif ORI==51217482  //  [[512,174,8]] QGB code weight 8  de-shih
  #define N (512)
  #define K (174)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_Gs ("codes/QGB512174_w8_G4s.txt")
  #define PATH_A  ("codes/QGB512174_w8_A4.txt")
  #define PRE_A   ("QGB512174_w8")







#elif ORI==2022   // [[2,0,2]], L=2
  #define N (4)
  #define K (0)
  #define PATH_A  ("Q202L2_A4")
  #define PATH_Gs ("Q202_G4s")

#elif ORI==412   // [[4,1,2]] code
  #define N (4)
  #define K (1)
  #define PATH_A  ("Q412_A4")
  #define PATH_Gs ("Q412_G4s")

#elif ORI==422   // [[4,2,2]] code
  #define N (4)
  #define K (2)
  #define PATH_A  ("Q422_A4")
  #define PATH_Gs ("Q422_G4s")

#elif ORI==513     // [[5,1,3]] code
  #define N (5)
  #define K (1)
  #define M   (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/Q513_A4.txt")
  #define PATH_Gs ("codes/Q513_G4.txt")
  #define PRE_A    ("513")
  //#define PRE_A    ("513_DYNqn")    // DYNAMIC_qn (only for test, not used)
  //#define PRE_A    ("513_DYNpn")    // DYNAMIC_pn
 #if USE_GF2_DEC
  #define PATH_A_GF2  ("Q_A2_La")
 #endif

#elif ORI==5131     // [[5,1,3]] code, +1 redundant check row
  #define N (5)
  #define K (1)
  #define M (N-K+1)
  #define PATH_A  ("codes/Q5131_A4.txt")
  #define PATH_Gs ("codes/Q513_G4.txt")
  #define PRE_A    ("5131")

#elif ORI==5132     // [[5,1,3]], L=2
  #define N (10)
  #define K (2)
  #define PATH_A  ("Q_L2v6_A4")
  #define PATH_Gs ("Q_G4s")

#elif ORI==713     // [[7,1,3]] code
  #define N (7)
  #define K (1)
  #define PATH_A  ("codes/Q713_A4")
  #define PATH_Gs ("codes/Q713_G4s")
  #define PRE_A    ("713")
 #if USE_GF2_DEC
  #define PATH_A_GF2  ("Q713_A2_La")
 #endif


#elif ORI==1115   // [[11,1,5]] code
  #define N (11)
  #define K (1)
  #define PATH_A  ("Q11_1_5_A4")
  #define PATH_Gs ("Q11_1_5_G4s")
  #define PRE_A   ("11_1_5")
  /*#define PATH_A  ("Q11_1_5_std_A4")
  #define PATH_Gs ("Q11_1_5_std_G4s")
  #define PRE_A   ("11_1_5_std")*/
  /*#define PATH_A  ("Q11_1_5_sparse_rand_A4")
  #define PATH_Gs ("Q11_1_5_sparse_rand_G4s")
  #define PRE_A   ("11_1_5_sparse_rand")*/


/*#elif ORI==13134  // [[13,1,3]] code, Surface Code, X (M=N-1), extended 4 redundant rows    (NOT USED)
  #define Surf_L  2  // Length of surface lattice
  #define N (13)
  #define K (1)
  #define M (N-K+4)
  #define M_G   (N+K)
  #define PATH_A  ("Surfxe_2x2_A4")
  #define PATH_Gs ("Surfx_2x2_G4s")
  #define PRE_A   ("surfXe_2x2")
*/

#elif ORI==4115   // [[41,1,5]] code, Surface Code, X (M=N-1)                               (WE CHOOSE THIS)
  #define Surf_L  4  // Length of surface lattice  (only EVEN) ///////////////////////////////////////
  #define K (1)
    #if Surf_L==2
  #define N (13)  // N = L^2 + (L+1)^2
  #define PATH_A  ("Surfx_2x2_A4")
  #define PATH_Gs ("Surfx_2x2_G4s")
  #define PRE_A   ("surfX_2x2")
    #elif Surf_L==4
  #define N (41)  // N = L^2 + (L+1)^2
  #define PATH_A  ("Surfx_4x4_A4")
  #define PATH_Gs ("Surfx_4x4_G4s")
  #define PRE_A   ("surfX_4x4")
    #elif Surf_L>=6
  #define PATH_A  ("Surfx_" STR(Surf_L) "x" STR(Surf_L) "_A4")
  #define PATH_Gs ("Surfx_" STR(Surf_L) "x" STR(Surf_L) "_G4s")
  #define PRE_A   ("surfX_" STR(Surf_L) "x" STR(Surf_L))
    #endif

/*#elif ORI==1823  // [[18,2,3]] code, Toric Code                                            (NOT USED)
  #define Surf_L  17  // Length of surface lattice ///////////////////////////////////////////////////
  #define N (Surf_L*Surf_L*2)
  #define K (2)
  #define PATH_A  ("Surfo2_" STR(Surf_L) "x" STR(Surf_L) "_A4")
  #define PATH_Gs ("Surfo2_" STR(Surf_L) "x" STR(Surf_L) "_G4s")
  #define PRE_A   ("surfo2_" STR(Surf_L) "x" STR(Surf_L))
*/


#elif ORI==625258  // [[625,25,8]] HP code
  #define N (625)
  #define K (25)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QHP_625_25_8_A4.txt")
  #define PATH_Gs ("codes/QHP_625_25_8_G4s.txt")
  #define PRE_A   ("QHP625_25")


#elif ORI==6555  // [[65,5,5]] lifted-connected surface code
  #define N (65)
  #define K (5)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QLC655_A4.txt")
  #define PATH_Gs ("codes/QLC655_G4s.txt")
  #define PRE_A   ("QLC655")

#elif ORI==1533  // [[15,3,3]] lifted-connected surface code
  #define N (15)
  #define K (3)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QLC153_A4.txt")
  #define PATH_Gs ("codes/QLC153_G4s.txt")
  #define PRE_A   ("QLC153")

#elif ORI==17577  // [[175,7,7]] lifted-connected surface code
  #define N (175)
  #define K (7)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QLC1757_A4.txt")
  #define PATH_Gs ("codes/QLC1757_G4s.txt")
  #define PRE_A   ("QLC1757")

#elif ORI==36999  // [[369,9,9]] lifted-connected surface code
  #define N (369)
  #define K (9)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QLC3699_A4.txt")
  #define PATH_Gs ("codes/QLC3699_G4s.txt")
  #define PRE_A   ("QLC3699")



#elif ORI==39232  // [[392,32]] code, QT code
  #define N (392)
  #define K (32)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QT39232_A4")
  #define PATH_Gs ("codes/QT39232_G4s")
  #define PRE_A   ("QT39232")

#elif ORI==3923212  // [[392,32,12]] code, QT code 1523
  #define N (392)
  #define K (32)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QGT392_32_1523_H4.txt")
  #define PATH_Gs ("codes/QGT392_32_1523_G4.txt")
  #define PRE_A   ("QGT392_32_1523")




#elif ORI==39236  // [[392,36]] code, QGT code
  #define N (392)
  #define K (36)
  #define M (356)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QGT39236_1463_H4.txt")
  #define PATH_Gs ("codes/QGT39236_1463_G4.txt")
  #define PRE_A   ("QGT39236")


#elif ORI==392322  // [[392,32]] code, QGB code, de-shih
  #define N (392)
  #define K (32)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QGB39232_w12_A4")
  #define PATH_Gs ("codes/QGB39232_w12_G4s")
  #define PRE_A   ("QGB39232")


#elif ORI==3923214  // [[392,32,14]]  31 code, QGB code, HYL
  #define N (392)
  #define K (32)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QGB392_32_31_H4.txt")
  #define PATH_Gs ("codes/QGB392_32_31_G4.txt")
  #define PRE_A   ("QGB3923214")


  #elif ORI==3924814  // [[392,48,14]]  rw=13,   7 code, QGB code, HYL
  #define N (392)
  #define K (48)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QGB392_48_7_H4.txt")
  #define PATH_Gs ("codes/QGB392_48_7_G4.txt")
  #define PRE_A   ("QGB3924814")

#elif ORI==3804610  // [[380,46,10]] QGB code, rw=13, HYL
  #define N (380)
  #define K (46)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QGB380_46_6_H4.txt")
  #define PATH_Gs ("codes/QGB380_46_6_G4.txt")
  #define PRE_A   ("QGB3804610")


#elif ORI==3923272  // [[392,32,7]] code, QGB code, HYL
  #define N (392)
  #define K (32)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QGB392_32_31_A4.txt")
  #define PATH_Gs ("codes/QGB392_32_31_G4.txt")
  #define PRE_A   ("QGB392327_2")

#elif ORI==3923217  // [[392,32,7]] code, QGB code, HYL rw=12 #17
  #define N (392)
  #define K (32)
  #define M (N)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QGB3923272_A4.txt")
  #define PATH_Gs ("codes/QGB3923272_G4.txt")
  #define PRE_A   ("QGB392327_2")

#elif ORI==39248  // [[392,48]] code, QT code
  #define N (392)
  #define K (48)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QT392481462_A4.txt")
  #define PATH_Gs ("codes/QT392481462_G4.txt")
  #define PRE_A   ("QT39248")


#elif ORI==4804012  // [[480,40,12]] code, QGB code, HYL  #6
  #define N (480)
  #define K (40)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QGB48040_6_A4.txt")
  #define PATH_Gs ("codes/QGB48040_6_G4.txt")
  #define PRE_A   ("QGB4804012")


#elif ORI==500194  // [[500,194,4]] code,  HYL code
  #define N (500)
  #define K (194)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/Q500194_A4.txt")
  #define PATH_Gs ("codes/Q500194_G4.txt")
  #define PRE_A   ("500194")

#elif ORI==50019465  // [[500,194,4]] code,  HYL code  65 GTS code
  #define N (500)
  #define K (194)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/GTS50019465_A4.txt")
  #define PATH_Gs ("codes/GTS50019465_G4.txt")
  #define PRE_A   ("GTS500194_65")

#elif ORI==500196  // [[500,196]] code,  GT
  #define N (500)
  #define K (196)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/GT5001961280_A4.txt")
  #define PATH_Gs ("codes/GT5001961280_G4.txt")
  #define PRE_A   ("GT500196")


#elif ORI==500188  // [[500,194,4]] code,  HYL code
  #define N (500)
  #define K (188)
    #define M (N-K)
    #define M_G   (N+K)
  #define PATH_A  ("codes/Q500188_A4.txt")
  #define PATH_Gs ("codes/Q500188_G4.txt")
  #define PRE_A   ("500188")



  #elif ORI==68628  // [[686,28,14]] code, QGTS code 686_28_2017, HYL
  #define N (686)
  #define K (28)
  #define M (658)
  #define M_G   (N+K)
  #define PATH_A  ("codes/GTS68628_2017_H4.txt")
  #define PATH_Gs ("codes/GTS68628_2017_G4.txt")
  #define PRE_A   ("QGT68628_2017")

#elif ORI==68634  // [[686,34,13]] code, QGT code 686_34_1306, HYL
  #define N (686)
  #define K (34)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QGT68634_1306_H4.txt")
  #define PATH_Gs ("codes/QGT68634_1306_G4.txt")
  #define PRE_A   ("QGT68634_1306")

#elif ORI==68830  // [[688,30,19]] code, QGB code 688_30_16, HYL
  #define N (688)
  #define K (30)
  #define M (N)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QGB68830_16_H4.txt")
  #define PATH_Gs ("codes/QGB68830_16_G4.txt")
  #define PRE_A   ("QGB68830_16_688check")

#elif ORI==488      // (4,8,8) color codes arXiv:1108.5738
  #define N (hex_d*hex_d +2*hex_d -1) /2
  #define K 1
  #define M (N-K)
  #define PATH_Gs ("codes/Color_488_d" STR(hex_d) "_G4s")
  #define PATH_A  "codes/Color_488_d" STR(hex_d) "_A4"
  #define PRE_A   ("color_488_d" STR(hex_d))



#elif ORI==28812   // [[288,12,18]] BB code
  #define N (288)
  #define K (12)
  #define M (N-K)
  #define M_G  (N+K)
  #define PATH_A  ("codes/QBB28812_A4.txt")
  #define PATH_Gs ("codes/QBB28812_G4s.txt")
  #define PRE_A   ("QBB28812")

#elif ORI==14412   // [[144,12,12]] BB code
  #define N (144)
  #define K (12)
  #define M (N-K)
  #define M_G  (N+K)
  #define PATH_A  ("codes/QBB14412_A4.txt")
  #define PATH_Gs ("codes/QBB14412_G4s.txt")
  #define PRE_A   ("QBB14412")


#elif ORI==1088   // [[108,8,10]] BB codee
  #define N (108)
  #define K (8)
  #define M (N-K)
  #define M_G  (N+K)
  #define PATH_A  ("codes/QBB1088_A4.txt")
  #define PATH_Gs ("codes/QBB1088_G4s.txt")
  #define PRE_A   ("QBB1088")


#elif ORI==908   // [[90,8,10]] BB codee
  #define N (90)
  #define K (8)
  #define M (N-K)
  #define M_G  (N+K)
  #define PATH_A  ("codes/QBB908_A4.txt")
  #define PATH_Gs ("codes/QBB908_G4s.txt")
  #define PRE_A   ("QBB908")

#elif ORI==7212   // [[72,12,6]] BB cod
  #define N (72)
  #define K (12)
  #define M (N-K)
  #define M_G  (N+K)
  #define PATH_A  ("codes/QBB7212_A4.txt")
  #define PATH_Gs ("codes/QBB7212_G4s.txt")
  #define PRE_A   ("QBB7212")


#elif ORI==11357   // [[113,57,12]] code from Grassl's table
  #define N (113)
  #define K (57)
  #define PATH_A  ("Q113_57_12_poly2_A4")
  #define PATH_Gs ("Q113_57_12_poly2_G4s")
  #define PRE_A    ("113_57_12_poly2")


#elif ORI==7134    // [[7,1,3]] (color) code, cyclic shift all rows (4 more redundant rows, with M=7+7 in H)
  #define hex_d  5  // must be 3 or 5 per the following #if
  #define N ((3*hex_d*hex_d +1) /4)
  #define K (1)

    #if 0 //-- only for d=3
  #define PATH_Gs ("Color_hex_cyc_" STR(hex_d) "_G4s")

  #define M (N+N)
  #define PATH_A  ("Color_hex_cyc_all_" STR(hex_d) "_A4")
  #define PRE_A   ("color_hex_cyc_all_" STR(hex_d))
  //#define PATH_A  ("Color_hex_cycintlv_all_" STR(hex_d) "_A4")
  //#define PRE_A   ("color_hex_cycintlv_all_" STR(hex_d))

  /*#define M (4+4) //-- test less rows                   // set to 4+4, 5+5, or 6+6
  #define PATH_A  ("Color_hex_cyc_r4_" STR(hex_d) "_A4")  // set to r4,  r5,  or r6
  #define PRE_A   ("color_hex_cyc_r4_" STR(hex_d))        // set to r4,  r5,  or r6*/

    #else //-- only for d=5
  #define PATH_Gs ("Color_hex_" STR(hex_d) "_G4s")

  /*#define M (N-K+6)   // HX has 3 redundant rows (H has 3+3=6 more redundant rows)
  #define PATH_A  ("Color_hexR3_" STR(hex_d) "_A4")
  #define PRE_A   ("color_hexR3_" STR(hex_d))*/

  /*#define M (N-K+8)   // HX has 4 redundant rows (H has 4+4=8 more redundant rows)
  #define PATH_A  ("Color_hexR4_" STR(hex_d) "_A4")
  #define PRE_A   ("color_hexR4_" STR(hex_d))*/

  #define M (N-K)   // HX has 0 redundant rows (that is full-rank based on hexR4)
  #define PATH_A  ("Color_hexR0_" STR(hex_d) "_A4")
  #define PRE_A   ("color_hexR0_" STR(hex_d))
    #endif


#elif ORI==1023    // [[10,2,3]] code, Toric Code TWISTED  (cyclic based on XZZX)
  /* OLD USELESS
  //#define N 5
  //#define K (1)
  #define N 50    // 10,20,50
  #define K (2)
  #define PATH_Gs ("Surf_XZ_twist_" STR(N) "_G4s")
    #if 1
  #define M (N)
  #define PATH_A  ("Surf_XZ_twistR_" STR(N) "_A4")
  #define PRE_A   ("surf_XZ_twistR_" STR(N))
    #else
  #define M (N-K)
  #define PATH_A  ("Surf_XZ_twist_" STR(N) "_A4")
  #define PRE_A   ("surf_XZ_twist_" STR(N))
    #endif*/

  /* old, only for Twist=1
  #define Surf_L  10  // d=L+1 if L even, d=L if L odd  (L=2: [[5,1,3]], L=3: [[10,2,3]], L=4: [[17,1,5]],  L=5: [[26,2,5]], ...)
  #define N (Surf_L*Surf_L+1)
  #define K (2-(N%2))
  #define PATH_Gs ("Surf_XZ_cyc_" STR(Surf_L) "_G4s")
    #if 1
  #define M (N)
  #define PATH_A  ("Surf_XZ_cycR_" STR(Surf_L) "_A4")
  #define PRE_A   ("surf_XZ_cycR_" STR(Surf_L))
    #else
  #define M (N-K)
  #define PATH_A  ("Surf_XZ_cyc_" STR(Surf_L) "_A4")
  #define PRE_A   ("surf_XZ_cyc_" STR(Surf_L))
    #endif*/

  #define Surf_L  17   // L
  #define Twist_R 16 // R, needs to coprime with L
  #define N (Surf_L*Surf_L + Twist_R*Twist_R)
  #define K (2-(N%2))
  #define D (Surf_L + Twist_R)
  #define PATH_Gs ("Surf_XZ_tw" STR(Twist_R) "_L" STR(Surf_L) "_G4s")
  #define M (N)
  #define PATH_A  ("Surf_XZ_tw" STR(Twist_R) "_L" STR(Surf_L) "_a_A4")
  #define PRE_A   ("surf_XZ_tw" STR(Twist_R) "_L" STR(Surf_L))


#elif ORI==18232  // [[18,2,3]] code, Toric Code, 2 redundant rows                          (WE CHOOSE THIS)
  #define Surf_L  11  // Length of surface lattice  (only ODD) ////////////////////////////////////
  #define N (Surf_L*Surf_L*2)
  #define K (2)
  #define M (N-K+2)                                               // (typically surfo is full-rank and surfo2 is redundant)
  #define PATH_A  ("Surfo_" STR(Surf_L) "x" STR(Surf_L) "_A4")    // but this case surfo is redundant
  #define PATH_Gs ("Surfo2_" STR(Surf_L) "x" STR(Surf_L) "_G4s")  // and surfo2 is full-rank for this case
  #define PRE_A   ("surfo_" STR(Surf_L) "x" STR(Surf_L))


#elif ORI==800400   // [[800,400]] bicycle code
  #define N (800)
  #define K (400)
  #define M (N-K)
  #define PATH_A  ("codes/Q800_400_rng7_A4")
  #define PATH_Gs ("codes/Q800_400_rng7_G4s")
  #define PRE_A    ("QB800_400_rng7")
//  #define PATH_A  ("Q800_400_rng7Hjmin_A4")
//  #define PATH_Gs ("Q800_400_rng7Hjmin_G4s")
//  #define PRE_A    ("800_400_rng7Hjmin")
//#if USE_GF2_DEC
//  #define PATH_A_GF2 ("Q800_400_rng7Hjmin_A2_La")
// #endif

#elif ORI==800404  // [[800,404]] bicycle code (check matrix A4 has 4 redundant rows)
  #define N (800)
  #define K (404)
  #define M (N-K+4)
  #define PATH_A  ("Q800_404_wv16_rng7PE_A4")
  #define PATH_Gs ("Q800_404_wv16_rng7PE_G4s")
  #define PRE_A    ("800_404_wv16_rng7PE")

#elif ORI==800418  // [[800,418]] bicycle code (check matrix A4 has 18 redundant rows)
  #define N (800)
  #define K (418)
  #define M (N-K+18)
  #define PATH_A  ("Q800_418_wv20_QE_A4")
  #define PATH_Gs ("Q800_418_wv20_QE_G4s")
  #define PRE_A    ("800_418_wv20_QE")


#elif ORI==80040016  // [[800,400]]  GB rw15   800400_16
  #define N (800)
  #define K (400)
  #define M (N-K)
  #define PATH_A  ("codes/QGB800400_16_H4.txt")
  #define PATH_Gs ("codes/QGB800400_16_G4.txt")
  #define PRE_A    ("QGB800400_16")


#elif ORI==1915   // [[19,1,5]] Color code
  #define K (1)
  #define hex_d  17  // d_min (only ODD) ////////////////////////////////////
  #define D (hex_d)
    #if 1   // (6,6,6) style: triangle plane with hex-cell
  #define N ((3*hex_d*hex_d +1) /4)

  //-- for d=3,5,7,9,...
  #define PATH_Gs ("Color_hex_" STR(hex_d) "_G4s")
  #define PATH_A  ("Color_hex_" STR(hex_d) "_A4")
  #define PRE_A   ("color_hex_" STR(hex_d))
  //#define PATH_A  ("Color_hex_intlv_" STR(hex_d) "_A4")
  //#define PRE_A   ("color_hex_intlv_" STR(hex_d))
 #if USE_GF2_DEC
  #define PATH_A_GF2 ("Color_hex_" STR(hex_d) "_A2_La")
 #endif

  //-- only for d=3
  /*#define PATH_Gs ("Color_hex_cyc_" STR(hex_d) "_G4s")
  //#define PATH_A  ("Color_hex_cyc_" STR(hex_d) "_A4")
  //#define PRE_A   ("color_hex_cyc_" STR(hex_d))
  #define PATH_A  ("Color_hex_cycintlv_" STR(hex_d) "_A4")
  #define PRE_A   ("color_hex_cycintlv_" STR(hex_d))*/

    #elif 0    // (4,8,8) style
  #define N ((hex_d*hex_d +2*hex_d -1) /2)
  #define PATH_Gs ("Color_488_d" STR(hex_d) "_G4s")
  #define PATH_A  ("Color_488_d" STR(hex_d) "_A4")
  #define PRE_A   ("color_488_d" STR(hex_d))
    #endif







#elif ORI==3786      // [[3786,946]] bicycle code
  #define N (3786)
  #define K (946)
  #define PATH_A  ("Q3786_946_rng7Hk24_A4")
  #define PATH_Gs ("Q3786_946_rng7Hk24_G4s")
  #define PRE_A    ("3786_946_rng7Hk24")
 #if USE_GF2_DEC
  #define PATH_A_GF2 ("Q3786_946_rng7Hk24_A2_La")
 #endif

  /*#define K (1894)
  #define PATH_A  ("Q3786_1894_rng7k24_A4")
  #define PATH_Gs ("Q3786_1894_rng7k24_G4s")*/


#elif ORI==25632   // [[256,32]] bicycle code
  #define N (256)
  #define K (32)
  #define PATH_A  ("Q256_32_rng3_A4")
  #define PATH_Gs ("Q256_32_rng3_G4s")
  #define PRE_A    ("256_32_rng3")
  //#define PRE_A    ("256_32_rng3_PARqn")  // PARTIAL_UPD_BY_qn
  //#define PRE_A    ("256_32_rng3_DYNpn")    // DYNAMIC_pn
  //#define PRE_A    ("256_32_rng3_DYNpnALL") // DYNAMIC_pnALL
 #if USE_GF2_DEC
  #define PATH_A_GF2 ("Q256_32_rng3_A2_La")
 #endif

#elif ORI==25634   // [[256,34]] bicycle code (check matrix A4 has 2 redundant rows)
  #define N (256)
  #define K (34)
  #define M (N-K+2)
  #define PATH_A  ("Q256_34_rng3E_A4")
  #define PATH_Gs ("Q256_34_rng3E_G4s")
  #define PRE_A    ("256_34_rng3E")



#elif ORI==251580    // [[25,1,5]] code, Surface Code ROTATED (ALL Redundant ROWs wt<=4, e.g., d=5 has M=80)
  #define Surf_d  3  // d_min (only ODD)
  #define N (Surf_d*Surf_d)
  #define K (1)
  //-- CSS style (XXXX or ZZZZ)
  #define PATH_Gs ("Surf_rot_" STR(Surf_d) "x" STR(Surf_d) "_G4s")
  #define PRE_A   ("surf_rot_" STR(Surf_d) "x" STR(Surf_d) "_redun")
  #define PATH_A  ("Surf_rot_" STR(Surf_d) "x" STR(Surf_d) "_redun_A4")
  #define M ( (Surf_d==3)? (26) : ( (Surf_d==5)? (80) : ( (Surf_d==7)? (158) : ( (Surf_d==9)? (260) : (0) ) ) ) )

  //-- XZZX  style
  // (NOT CONSIDERED YET)


#elif ORI==2515    // [[25,1,5]] code, Surface Code ROTATED (also support XZZX)
  #define N (Surf_d*Surf_d)
  #define K (1)
  //#define D (Surf_d)
  #define M (N-K)
  //-- CSS style (XXXX or ZZZZ)
  #define PATH_A  ("codes/Surf_rot_" STR(Surf_d) "x" STR(Surf_d) "_A4")
  #define PATH_Gs ("codes/Surf_rot_" STR(Surf_d) "x" STR(Surf_d) "_G4s")
  #define PRE_A   ("rot_surf_" STR(Surf_d) "x" STR(Surf_d))
  //#define PRE_A   ("surf_rot_" STR(Surf_d) "x" STR(Surf_d) "_DYNpn")
  //-- XZZX  style
  /*#define PATH_A  ("Surf_XZ_rot_" STR(Surf_d) "x" STR(Surf_d) "_A4")
  #define PATH_Gs ("Surf_XZ_rot_" STR(Surf_d) "x" STR(Surf_d) "_G4s")
  #define PRE_A   ("surf_XZ_rot_" STR(Surf_d) "x" STR(Surf_d))*/


#elif ORI==1624    // [[16,2,4]] code, Toric Code ROTATED
  //-- CSS style (XXXX or ZZZZ)
  #define Surf_L  10      // L=d_min in this case (only EVEN)
  #define N (Surf_L*Surf_L)
  #define K (2)
  #define D (Surf_L)
  #define M (N-K+2)
  #define PATH_A  ("Surf_rotO2_" STR(Surf_L) "x" STR(Surf_L) "_A4")   // rotO2 is redundant
  #define PATH_Gs ("Surf_rotO_" STR(Surf_L) "x" STR(Surf_L) "_G4s")   // rotO is full-rank
  #define PRE_A   ("surf_rotO_" STR(Surf_L) "x" STR(Surf_L))
  //#define PRE_A   ("surf_rotO_" STR(Surf_L) "x" STR(Surf_L) "_dlt-16")
 #if USE_GF2_DEC
  #define PATH_A_GF2 ("Surf_rotO2_" STR(Surf_L) "x" STR(Surf_L) "_A2_La")
 #endif
  //-- XZZX style
  /*#define Surf_L  5  // L=d_min in this case (support both EVEN and ODD)
  #define N (Surf_L*Surf_L)
  #define K (2-(N%2))
  #define M (N)
  #define PATH_A  ("Surf_XZ_rotOr_" STR(Surf_L) "x" STR(Surf_L) "_A4")   // rotOr has 2 (or 1) redundant rows if d is EVEN (or ODD)
  #define PATH_Gs ("Surf_XZ_rotO_" STR(Surf_L) "x" STR(Surf_L) "_G4s")   // rotO is full-rank
  #define PRE_A   ("surf_XZ_rotO_" STR(Surf_L) "x" STR(Surf_L))*/


#elif ORI==12928   // [[129,28]] hypergraph-product code by [7,4,3] and [15,7,5] BCH codes
  #define N (129)
  #define K (28)
  #define PATH_Gs ("bch_hyper_G4s")

  //-- normal H
  #define M (N-K)
  #define PATH_A  ("bch_hyper_A4")
  #define PRE_A   ("bch_hpyper_129_28")
    #if USE_GF2_DEC
  #define PATH_A_GF2  ("bch_hyper_A2_La")
    #endif
  //-- flip H as column wt >= 2
  /*#define M (N-K)
  #define PATH_A  ("bch_hyper_c2more_AX_AZ_v2_A4")
  #define PRE_A   ("bch_hpyper_129_28_c2more_AX_AZ_v2")
  //#define PATH_A  ("bch_hyper_c2more_v2_A4")
  //#define PRE_A   ("bch_hpyper_129_28_c2more_v2")*/

  //-- Cyclic-Extended H
  /*#define M (210)
  #define PATH_A  ("bch_hyper_Ext_A4")
  #define PRE_A   ("bch_hpyper_Ext_129_28")*/

#elif ORI==12628   // [[126,28]] GB code
  #define N (126)
  #define K (28)
  #define M_G (N+K)
  #define PATH_Gs ("Q126_28c_G4s")
  #define M (N)
  #define PATH_A  ("Q126_28c_A4")
  #define PRE_A    ("126_28c")
  /*#define M (N-K)
  #define PATH_A  ("Q126_28c_N-K_A4")
  #define PRE_A    ("126_28c_N-K")*/
    //#if USE_GF2_DEC
    //#define PATH_A_GF2 ("Q126_28_A2_La")  // to add H2 matrix by col_vec if needed
    //#endif



#elif ORI==88224   // [[882,24]] GHP code
  #define N (882)
  #define K (24)
 // #define D(16)
  #define PATH_Gs ("codes/GHP882_24_16_G4.txt")
  #define M (882)
  #define M_G   (N+2*K)
  #define PATH_A  ("codes/GHP882_24_16_H4.txt")
  #define PRE_A    ("GHP882_24")


#elif ORI==88248   // [[882,48]] GHP code
  #define N (882)
  #define K (48)
 // #define D(16)
  #define PATH_Gs ("codes/Q88248_G4s")
  #define M (834)
  #define M_G   (N+K)
  #define PATH_A  ("codes/Q88248_A4")
  #define PRE_A    ("882_48")


#elif ORI==882482   // [[882,48]] GB code, de-shih
  #define N (882)
  #define K (48)
 // #define D(16)
  #define PATH_Gs ("codes/QGB88248_w8_G4s")
  #define M (834)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QGB88248_w8_A4")
  #define PRE_A    ("GB882_48_w8")

 #elif ORI==1844    // [[18,4,4]] code, Color code (Toric & ROTATED)
  #define d_min  24  // d_min (only be a multiple of 4)
  #define N (d_min*d_min*9/8)
  #define K (4)
  #define D (d_min)
  #define PATH_Gs ("Color_toric_rot_d" STR(d_min) "_fulrank_G4s")
  #define PRE_A   ("Color_toric_rot_d" STR(d_min) )
  #define PATH_A  ("Color_toric_rot_d" STR(d_min) "_A4")
  #define M       (N)
  #define M_G     (N+K)

#else
  XXXXXX
#endif


//#if ORI==5131 || ORI==13134 || ORI==1023 || ORI==1624 || ORI==251580 || ORI==18232 || ORI==12928 || ORI==25634 || ORI==800404 || ORI==800418 || ORI==7134 || ORI == 12628|| ORI==88248 || ORI == 1844
//  // M already defined
//#else
//  #define M   (N-K)
//#endif
//
//# if ORI == 1624 || ORI == 1844 || ORI == 2515 || ORI==88248 || ORI == 1023 || ORI== 1915
//    // D already defined
//#else
//    # define D 0
//# endif // ORI

# define D 0

//#define M_G   (N+K)


#define Q 4       // GF(Q=4)
typedef uint8_t   GFQ_t;   // current only support up to GF(256)

//extern FILE *fpErr;




