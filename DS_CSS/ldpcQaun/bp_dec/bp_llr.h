#include "../ldpc_parm.h"


typedef struct {  // used with g_matrix_GFQ
  uint8_t reset;      // set to 1 when init
  GFQ_t target_z[M];  // target syndrome vector to converge

  // INIT and VERTICAL:
  double LA[N][Q-1];  // \Lambda_n^W in the algorithm
  double GA[N][Q-1];  // \Gamma_n^W in the algorithm
  //double (*GAmj[M])[Q-1];    // \Gamma_mn^W in the algorithm (n\to m)
  double *Lmj[M];     // \lambda_mn in the algorithm  (n\to m)

  // HORIZONTAL:
  double Dm[M];       // \boxplus of Dmn[m][n] for all n of some m
  double *Dmj[M];     // \Delta_mn in the algorithm (n\to m)


  GFQ_t tt[N];    // updated decoded vector (estimated TX vector after BP)
  GFQ_t zz[M];    // updated syndrome vector


  int32_t iter;         // iter cnt

#if 0 //STOP_BY_C   // These variables do not need to be delivered between iterations,
  double dm_c[M];     // put them here due to:
  double *dqml_c[M];  // can create the data structure of dqml_c with dqml
  GFQ_t zz_c[M];
#endif
  uint8_t is_hit_c;     // indicator for if hit the required condition for STOP_BY_C (CheckNode syndrome OK)

  double ap;      // to pass the parameter alpha prime (AFP)

  int32_t seq[N]; // should contain 0 to N-1 in some order to run serial decoder (dec64)

#if 0 //MAJ_SCHE
  int32_t maj_cnt[N];  // add one vote if a stabilizer associated to this qubit has syndrome=1
  int32_t maj_max;     // maximum of maj_cnt
#endif

#if 0 //CYC_CHK
  int32_t *cyc_cnt[N];  // cycle count [n][m] (0: no cycles (1 overlap), 1: has one 4-cycle (2 overlap), k: has k+1 overlaps with other rows)

  int32_t **cycs[N];    // [N][sM(n)][m] so that we can for m = 1 to sM(n)
  double (*rml2[M])[2]; // to buffer original rml
#endif

} LBP_Ctl ;  // "LLR-based" Quantum Belief Propagation (Sum-Product, Message-Passing) Control block



//-- alloc for LLR-based Quantum BP Ctrl Blk
void alloc_LBPC(a_matrix_GFQ *A, LBP_Ctl *bp);

//-- LLR-based Quantum BP init: set bp->tt as all zeros. return 1 if zz pass; otherwise INIT *bp and return 0
uint8_t Lbp_init(LBP_Ctl *bp , GFQ_t *zz , a_matrix_GFQ *A , double p_ch, double afp); // zz: syndrome, p_ch: channel err rate, afp: see AFP

//-- Parallel LLR-based BP4
uint8_t Lbp_dec60(LBP_Ctl *bp, a_matrix_GFQ *A);  // copy form Qbp_dec20 and modify
//-- Serial LLR-based BP4
uint8_t Lbp_dec64(LBP_Ctl *bp, a_matrix_GFQ *A);  // copy from Qbp_dec24 and modify
//-- Serial LLR-based BP4 (along checks)
uint8_t Lbp_dec84(LBP_Ctl *bp, a_matrix_GFQ *A);  // copy form Qbp_dec44 and modify
