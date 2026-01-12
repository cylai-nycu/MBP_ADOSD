#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#define USE_GF2_DEC 0


/*
Define the quantum code to be simulated.
*/

#define OSDW  (-2)       // -1 no OSD, -2 APOSD, and others
#define ADOSDw (2)      // 0 for OSD-0, 2 for budget complexity approximately equal to OSD-2

#define RelThr  (0.99)  // reliability threshold
#define MAX_ITER  (10)

#define P_ERR 0.01
#define LE_tar   100000
#define Shot_max 100000

#define ShowTime 1000



#define rev 1   // 0: with  additional X detectors    1:  without  additional X detectors


//rotated surface codes without  additional X detectors

//#define ORI 55116       // [[9,1,3]] rotated surface code, circuit level  with only Z detectors
//#define ORI 301172      // [[25,1,5]] rotated surface code, circuit level with only Z detectors
//#define ORI 8831192     // [[49,1,7]] rotated surface code, circuit level  with only Z detectors
#define ORI 19451400    // [[81,1,9]] rotated surface code, circuit level  with only Z detectors


//rotated surface codes with  additional X detectors
//#define ORI 219124   // [[9,1,3]] rotated surface code, circuit level
//#define ORI 16771120   // [[25,1,5]] rotated surface code, circuit level
//#define ORI 54711336   // [[49,1,7]] rotated surface code, circuit level//
//#define ORI 127051720   // [[81,1,9]] rotated surface code, circuit level
//#define ORI 2448311320 // [[121,1,11]] rotated surface code, circuit level






/*

Define the decoder parameters

*/

#define LLR_BP  0   // 0: off (i.e., use linear BP), 1: use LLR-BP





#define AFP   150     // alpha prime, set 0 to off, set to 120 means 120/100  (\alpha in the MBP paper)

/* Adjust alpha adaptively in AMBP
 MBP:   Set AFP to a fixed value, e.g., 150, and set AFP2 = 0, AFPINC = 0.
         A larger AFP corresponds to a smaller step size.
 AMBP:  Set AFP as the starting value, AFP2 as the ending value, and AFPINC as the increment.
         Example: AFP = 160, AFP2 = 30, AFPINC = -1

*/

#define AFP2  0     // set to non-zero to enable:  if AFPINC!=0 further, then must set a non-zero AFP, e.g. AFP=100
#define AFPINC  (0) // (what AFP2 used:) 0: off; otherwise, afp refers to AFP+AFPINC*n for the additional n-th trial until reach AFP2
                    // (-1)
//random schedule
//AMBP with a random schedule usually performs better
#define RND_SCHE 0  // 0: off, random schedule 1: init , 2: per CW, 3: per iter



#define BY_DEC  20  // 20 is parallel dec, 24 is serial dec, 25 is circular schedule (Surface code only, 26 is syndrome-based schedule (not enabled now)
                    //                     44 is serial dec along check nodes
                    // 60 is parallel LLR, 64 is serial LLR
                    //                     84 is serial LLR along check nodes





//#define P_STEP  (-0.1)  // err probability step (in dB, every time 10^(STEP/10))
//#define P_STEP  (-0.02)  // err probability step (in dB, every time 10^(STEP/10))
//#define P_STEP  (-10.0/4)  // err probability step (in dB), every time p=p*10^(STEP/10))
// -10.0/n
#define P_STEP  (-10.0/8)  // err probability step (in dB, every time 10^(STEP/10))
#define FIX_P (0)       // set 0 to off, set to 46 means 4.6e-2 = 0.046 if FIX_e = 0 or -2
#define FIX_e (0)     // set 0 to off, set to -2 means 4.6e-2 if FIX_P = 46

#define ALFA  0     // set 0 to off, set to 120 means 100/120  (\alpha_c in refined BP)





  //#define DEC2    '0' // (what AFP2 used:) '0': DEC2 follows BY_DEC,  'p': DEC2 by dec20,  's': DEC2 by dec24,  'c': DEC2 by dec44
 // #define DBL_CHK 0   // if this afp has the same output as last afp (or this afp hits max iter), then stop


#define AFPADD  0   // alpha prime post added, set 0 to off, set to 120 means to add 100/120
                    // When message normalized by 100/AFP, the memory weight, originally (1 - 100/AFP), can be adjusted here

#define BETA  0     // set 0 to off, this is in linear domain (for efficiency, try to set to something like 2, 4, 8, 16, ...)

#define ALFA_V  0   // alpha_v, set 0 to off, set to 12 means 10/12
#define EN_AFP_V 0  // 0: disable alpha_v prime, 1: enable (and alpha_v controlled by ALFA_V)





#define CH_TYPE 0   // 0: dep ch,  1: indep X-Z

#define POS_AVG 0   // 0 is off, or set POS_AVG>=2: (old*(POS_AVG-1) + new)/POS_AVG


#define AFN   0     // alpha for init_pn, set 0 to off, set to 12 means 10/12

#define CYC_CHK 0   // 0: off, 1: check short cycles (and adjust edge weights)



#define GMN 0       // 0: off; e.g., if gain GMN=140, then (1/AFP)*=140/100 if |\sN(m)|=2.
#define LMN 0       // default 200; when GMN is enabled, the (1/AFP)*=140/100 will be upper bounded by LMN


#define LIST_DEC 0  // 0: off, 1: buffer the last syndrome-matched solution


#define STOP_BY_C 0   // to simulate the effect: stop at check node when syndrome match (0 Off, 1 ON, 2 TEST force OK)
#define ITER_LOG (1 && STOP_BY_C)   // 0: off, 1: log iteration data (must enable STOP_BY_C)


#define BETA_pn   0   // set 0 to off, set to 120 means 100/120




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
M_G: number of normalizers, usually N+K
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


#elif ORI==55116     // [[9,1,3]] surface code, circuit level    with only Z detectors
  #define N (55)
  #define K (1)
  #define M   (16)
  #define M_G   (1)
  #define PATH_A  ("codes/surface55_1_16_H4rev.txt")
  #define PATH_Gs ("codes/surface55_1_16_G4rev.txt")
  #define PATH_P  ("codes/surface55_1_16")
  #define PRE_A    ("55116")
  #define distance 3


#elif ORI==301172     // [[25,1,5]] surface code, circuit level with only Z detectors
  #define N (301)
  #define K (1)
  #define M   (72)
  #define M_G   (1)
  #define PATH_A  ("codes/surface301_1_72_H4rev.txt")
  #define PATH_Gs ("codes/surface301_1_72_G4rev.txt")
  #define PRE_A    ("301172")
  #define distance 5

#elif ORI==8831192     // [[49,1,7]] surface code, circuit level  with only Z detectors
  #define N (883)
  #define K (1)
  #define M   (192)
  #define M_G   (1)
  #define PATH_A  ("codes/surface883_1_192_H4rev.txt")
  #define PATH_Gs ("codes/surface883_1_192_G4rev.txt")
  #define PRE_A    ("8831192")
  #define distance 7


#elif ORI==19451400     // [[81,1,9]] surface code, circuit level  with only Z detectors
  #define N (1945)
  #define K (1)
  #define M   (400)
  #define M_G   (1)
  #define PATH_A  ("codes/surface1945_1_400_H4rev.txt")
  #define PATH_Gs ("codes/surface1945_1_400_G4rev.txt")
  #define PRE_A    ("19451400")
  #define distance 9





#elif ORI==219124     // [[9,1,3]] surface code, circuit level  with  additional X detectors
  #define N (219)
  #define K (1)
  #define M   (24)
  #define M_G   (1)
  #define PATH_A  ("codes/surface219_1_24_H4.txt")
  #define PATH_Gs ("codes/surface219_1_24_G4.txt")
  #define PATH_P  ("codes/surface219_1_24")
  #define PRE_A    ("219124")
  #define distance 3


#elif ORI==16771120     // [[25,1,5]] surface code, circuit level  with  additional X detectors
  #define N (1677)
  #define K (1)
  #define M   (120)
  #define M_G   (1)
  #define PATH_A  ("codes/surface1677_1_120_H4.txt")
  #define PATH_Gs ("codes/surface1677_1_120_G4.txt")
  #define PRE_A    ("16771120")
  #define distance 5


#elif ORI==127051720     // [[81,1,9]] surface code, circuit level  with  additional X detectors
  #define N (12705)
  #define K (1)
  #define M   (720)
  #define M_G   (1)
  #define PATH_A  ("codes/surface12705_1_720_H4.txt")
  #define PATH_Gs ("codes/surface12705_1_720_G4.txt")
  #define PRE_A    ("127051720")
  #define distance 9


#elif ORI==54711336     // [[49,1,7]] surface code, circuit level  with  additional X detectors
  #define N (5471)
  #define K (1)
  #define M   (336)
  #define M_G   (1)
  #define PATH_A  ("codes/surface5471_1_336_H4.txt")
  #define PATH_Gs ("codes/surface5471_1_336_G4.txt")
  #define PRE_A    ("54711336")
  #define distance 7

#elif ORI==2448311320     // [[121,1,11]] surface code, circuit level   with  additional X detectors
  #define N (24483)
  #define K (1)
  #define M   (1320)
  #define M_G   (1)
  #define PATH_A  ("codes/surface24483_1_1320_H4.txt")
  #define PATH_Gs ("codes/surface24483_1_1320_G4.txt")
  #define PRE_A    ("2448311320")
  #define distance 11





#elif ORI==713     // [[7,1,3]] code
  #define N (7)
  #define K (1)
  #define PATH_A  ("codes/Q713_A4")
  #define PATH_Gs ("codes/Q713_G4s")
  #define PRE_A    ("713")
 #if USE_GF2_DEC
  #define PATH_A_GF2  ("Q713_A2_La")
 #endif




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
  #define M_G   (N+K)
  #define PATH_A  ("codes/QBB28812_A4.txt")
  #define PATH_Gs ("codes/QBB28812_G4s.txt")
  #define PRE_A   ("QBB28812")

#elif ORI==14412   // [[144,12,12]] BB code
  #define N (144)
  #define K (12)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QBB14412_A4.txt")
  #define PATH_Gs ("codes/QBB14412_G4s.txt")
  #define PRE_A   ("QBB14412")


#elif ORI==1088   // [[108,8,10]] BB codee
  #define N (108)
  #define K (8)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QBB1088_A4.txt")
  #define PATH_Gs ("codes/QBB1088_G4s.txt")
  #define PRE_A   ("QBB1088")


#elif ORI==908   // [[90,8,10]] BB codee
  #define N (90)
  #define K (8)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QBB908_A4.txt")
  #define PATH_Gs ("codes/QBB908_G4s.txt")
  #define PRE_A   ("QBB908")

#elif ORI==7212   // [[72,12,6]] BB cod
  #define N (72)
  #define K (12)
  #define M (N-K)
  #define M_G   (N+K)
  #define PATH_A  ("codes/QBB7212_A4.txt")
  #define PATH_Gs ("codes/QBB7212_G4s.txt")
  #define PRE_A   ("QBB7212")





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



# define D 0

//#define M_G   (N+K)


#define Q 4       // GF(Q=4)  type of errors

typedef uint8_t   GFQ_t;   // current only support up to GF(256)

//extern FILE *fpErr;




