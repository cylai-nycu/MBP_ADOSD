#include <stdint.h>

#define RND_MAX	  0xffffffffffffffff              // = UINT64_MAX = 2^64 -1
#define RND_DMAX	((double)0xffffffffffffffff)    // (most CPU will get) 2^64 if cast with "rounding"

//- some ref:
// https://nullprogram.com/blog/2017/09/21/
// http://xoshiro.di.unimi.it/

//- implemented in splitmix64.c
uint64_t next64();

//- implemented in xoshiro256starstar.c
void rnd256_init();   //- must call first
void rnd256_init_priv(uint64_t s[]);

uint64_t next(void);
uint64_t next_priv(uint64_t s_priv[]);


//==== implemented in lib_rand.c ====// (START)
void rnd_max_show(void);  //- some test show

//-- uniform r.v. over [0.0, 1.0]
double rv_UnifOne(void);
double priv_UnifOne(uint64_t s_priv[]);

//-- uniform discrete r.v. over {0, 1, 2, ..., N-1}
uint64_t rv_Nminus1(uint64_t given_N);
uint64_t priv_Nminus1(uint64_t given_N, uint64_t s_priv[]);

//-- uniform discrete r.v. over {A, A+1, A+2, ..., B}, must with A<=B
int64_t rv_IntAB(int64_t A, int64_t B);
int64_t priv_IntAB(int64_t A, int64_t B, uint64_t s_priv[]);


//==== implemented in lib_rand.c ====// (END)
