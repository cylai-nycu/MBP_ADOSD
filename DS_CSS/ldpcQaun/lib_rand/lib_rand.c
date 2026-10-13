#include <inttypes.h>
#include <stdio.h>
#include <math.h>
#include "lib_rand.h"

//- some test show
void rnd_max_show(void){
  uint64_t rndmax_val; //= (uint64_t)( (double)RND_MAX );
  double rndmax_dbl; //= RND_MAX;
  uint64_t tmprand = next64();
  printf("tmprand = %"PRIx64"\n", tmprand );
  printf("tmprand: LSB32 %x , >>28 %x , >>32 %x\n", (uint32_t)tmprand, (uint32_t)(tmprand>>28), (uint32_t)(tmprand>>32) );

  printf("size double %d \n", sizeof(rndmax_dbl));
  rndmax_dbl = (double)(RND_MAX);
  rndmax_val = (uint64_t)(rndmax_dbl);
  printf("rndmax_dbl val = %f  cast back as uint64_t = %016"PRIx64" \n", rndmax_dbl, rndmax_val);
  rndmax_dbl = (double)(RND_MAX);
  rndmax_val = (uint64_t)(rndmax_dbl-1);
  printf("rndmax_dbl val = %f , (uint64_t)(rndmax_dbl-1) = %016"PRIx64" \n", rndmax_dbl, rndmax_val);

  double rndrmp_dbl;
  rndrmp_dbl = (double)(RND_MAX - 0x3fe);
  rndmax_val = (uint64_t)(rndrmp_dbl);
  printf("(double)(RND_MAX - 0x3fe) = %f , cast back as uint64_t = %016"PRIx64" \n", rndrmp_dbl, rndmax_val);
  rndrmp_dbl = (double)(RND_MAX - 0x3ff);
  rndmax_val = (uint64_t)(rndrmp_dbl);
  printf("(double)(RND_MAX - 0x3ff) = %f , cast back as uint64_t = %016"PRIx64" \n", rndrmp_dbl, rndmax_val);
  rndrmp_dbl = (double)(RND_MAX - 0x400);
  rndmax_val = (uint64_t)(rndrmp_dbl);
  printf("(double)(RND_MAX - 0x400) = %f , cast back as uint64_t = %016"PRIx64" \n", rndrmp_dbl, rndmax_val);
  rndrmp_dbl = (double)(RND_MAX - 0x401);
  rndmax_val = (uint64_t)(rndrmp_dbl);
  printf("(double)(RND_MAX - 0x401) = %f , cast back as uint64_t = %016"PRIx64" \n", rndrmp_dbl, rndmax_val);
  printf("When close to RND_AMX, the rounding step is 11 bits , but the rounding resolution is 10 bits by\n");
  printf("check 0xffffffffffffffff - 0x3ff = %016"PRIx64" round up\n", 0xffffffffffffffff - 0x3ff);
  printf("check 0xffffffffffffffff - 0x400 = %016"PRIx64" round down\n", 0xffffffffffffffff - 0x400);
  printf("\n");

  printf("Though there are 54 =64-10 effective bits \n");
  printf("(int)( (double)(0x400)/RND_MAX * pow(2,54) )  is  %d \n", (int)((double)(0x400)/RND_MAX * pow(2,54)) );
  printf("(int)( (double)(0x3ff)/RND_MAX * pow(2,54) )  is  %d \n", (int)((double)(0x3ff)/RND_MAX * pow(2,54)) );
  printf("when amplitude is small, the floating-point resolution is good: cf rv_UnifOne(): \n");
  printf("(double)(RND_MAX-1)/RND_MAX == 1.0  is TRUE or FALSE: %d \n", (double)(RND_MAX-1)/RND_MAX == 1.0);
  printf("(double)(RND_MAX-0x3ff)/RND_MAX == 1.0  is TRUE or FALSE: %d \n", (double)(RND_MAX-0x3ff)/RND_MAX == 1.0);
  printf("(double)(RND_MAX-0x400)/RND_MAX == 1.0  is TRUE or FALSE: %d \n", (double)(RND_MAX-0x400)/RND_MAX == 1.0);
  printf("(double)(0x400)/RND_MAX == 0.0  is TRUE or FALSE: %d \n", (double)(0x400)/RND_MAX == 0.0);
  printf("(double)(0x3ff)/RND_MAX == 0.0  is TRUE or FALSE: %d \n", (double)(0x3ff)/RND_MAX == 0.0);
  printf("(double)(0x1)/RND_MAX == 0.0  is TRUE or FALSE: %d \n", (double)(0x1)/RND_MAX == 0.0);
  printf("(double)(0x0)/RND_MAX == 0.0  is TRUE or FALSE: %d \n", (double)(0x0)/RND_MAX == 0.0);
  printf("\n");


  /*uint64_t tmprnd256, tmpnxt256, tmp_cnt=0;
  tmprnd256 = next();
  while(1){
    tmp_cnt++;
    tmpnxt256 = next();

    //-- to show there could be count more than 16 bits when there are 16 bits the same
    if( (tmprnd256 & 0xffff000000000000) == (tmpnxt256 & 0xffff000000000000) ){
      printf("rand256 hit ori when cnt = %08x %08x \n", (uint32_t)(tmp_cnt>>32) , (uint32_t)tmp_cnt );
      //getchar();
      break;
    }

    if(tmp_cnt==0)
      printf("cnt wrap\n");
  }*/

}



//-- uniform r.v. over [0.0, 1.0]
inline double rv_UnifOne(void)
{
	return ((double)next() / RND_MAX);
}

inline double priv_UnifOne(uint64_t s_priv[])
{
	return ((double)next_priv(s_priv) / RND_MAX);
}


//-- uniform discrete r.v. over {0, 1, 2, ..., N-1}
inline uint64_t rv_Nminus1(uint64_t given_N)
{
  uint64_t rt;  // result
#if 0 //-- THIS IS SLOWER (NB runs 13 sec) --//
  do{ rt = next(); } while( rt > RND_MAX - ((RND_MAX % given_N)+1)%given_N);
    // e.g. say that if RND_MAX=15, then for the given N = [1  2  3  4  5  6  7  8]
    // redo if rt > 15 - mod((mod(15, [1:8]))+1, [1:8]) = [15 15 14 15 14 11 13 15]
    //                                which case may redo:  x  x  o  x  o  o  o  x
  return(rt%given_N);
#else  //-- THIS IS FASTER (NB runs 10 sec)-- //  rm < N-1: x  x  o  x  o  o  o  x
  uint64_t rm = RND_MAX % given_N;  // remainder,     rm = [0  1  0  3  0  3  1  7]
  uint64_t bd = RND_MAX - rm;       // bound,         bd = [15 14 15 12 15 12 14 8]
  if(rm < given_N-1) {      // rm is in {0,1,2,...,N-2}: needs to cut the tail
    do{ rt = next(); } while( rt >= bd);
    return(rt % given_N);
  } else return(next() % given_N);  // rm = L-1: can output directly
#endif
}

inline uint64_t priv_Nminus1(uint64_t given_N, uint64_t s_priv[])
{
  uint64_t rt;  // result
  uint64_t rm = RND_MAX % given_N;  // remainder,     rm = [0  1  0  3  0  3  1  7]
  uint64_t bd = RND_MAX - rm;       // bound,         bd = [15 14 15 12 15 12 14 8]
  if(rm < given_N-1) {      // rm is in {0,1,2,...,N-2}: needs to cut the tail
    do{ rt = next_priv(s_priv); } while( rt >= bd);
    return(rt % given_N);
  } else return(next_priv(s_priv) % given_N);  // rm = L-1: can output directly
}

//-- uniform discrete r.v. over {A, A+1, A+2, ..., B}, must with A<=B
inline int64_t rv_IntAB(int64_t A, int64_t B)
{
  return(A + rv_Nminus1(B-A+1));
}

inline int64_t priv_IntAB(int64_t A, int64_t B, uint64_t s_priv[])
{
  return(A + priv_Nminus1(B-A+1, s_priv));
}
