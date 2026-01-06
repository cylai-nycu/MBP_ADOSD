#ifndef OSD_H
#define OSD_H

#include "../ldpc_parm.h"
#include <stdint.h>
#define LASTFORSORT     0
#define ENTROPYSORT     1
#define MAXSORT         2
#define BINARYSORT      3
#define SORTCOMPARE LASTFORSORT

typedef struct{
    uint64_t nStep1, nStep2;
    #if STORE_DISTR
    uint64_t newMatSizeDistr[M+1][(N<<1)+1];    // Distribution of new matrix of POSD.
    uint64_t newMatSizeDistr_osd0[M+1][(N<<1)+1];
    #endif // STORE_DISTR
}OsdRecord;

enum OsdStatus{
    UsePosd = 0,
    Stage1Fail = 1,
    Stage2Fail = 2,
    NormalOsd = 3
};

typedef struct{
    uint32_t index;
    double reliBin;
    double reliQua;
}Index_Rely ;


typedef struct
{
    int8_t osdw;                    // If osdw = -1, do not use osdw, and won't generate probability
    enum OsdStatus status;
    uint64_t usedTime;
    uint8_t osd0orNot;
	uint32_t NN, NNbit, MM, P_col, RankH;
	int32_t Dist;
    uint32_t n, m, p_col;
    RV_t mask;
	RV_t **parity;                  // MM*P_col
	RV_t **parity_intact;           // MM*P_col
	uint8_t *syndrome;              // MM
	uint8_t *tempSyndrome;          // MM

	uint32_t numOfUnreliCheck;
    uint32_t numOfUnreliBit;
    uint32_t numOfUnreliP_col;
    uint32_t numOfSneaky;
	Index_Rely *permu;              // NNbit    use this to permute the index according to reliability
	RV_t *unreliMask;               // P_col,   mask the bit that is unreliable.
    uint32_t *unreliCheck;          // MM       store the unreliable check
    uint32_t *unreliP_col;          // P_col    store the unreliable P_col
	RV_t *pivotMask;                // P_col,   mask the pivot when generate codeword
    int32_t *rowPivotPlace;         // MM,      if the row has pivot, record where the pivot is in P (-1: not indep row, -2: reliable check)
    int8_t *rowPivotStatus;         // MM       use to record the status of the bit in terms of the row
    uint32_t *rowPivotBrother;      // MM       use to record the where the brother is in terms of the row
    int32_t *pivotBrother;          // NNbit    use to know what row belongs to brother
    int32_t *matchBrother;          // NNbit
    uint32_t *sneakyBit;            // NNbit,   Store the number of pivot has found before that bit, we can use it to see it is reliable or not
    uint32_t **bit2RowPos;          // NNbit*MM to let syndrome update more convenient, since we know what position should change
    int32_t *bitColWeight;          // NNbit    use to know the column weight of each bit in parity, if >sqrt(N) may be normalizer

    uint32_t bestErrorWeight;
    RV_t *error;                    // P_col
    RV_t *mostPossibleError;        // P_col
    uint8_t* returnError;           // NN
}OSD;


// private function
static void initParity(RV_t **parity,  RV_t **parity_intact);
static void getSyndrome(const RV_t *error, const RV_t **parity, uint8_t *syndrome, uint32_t row, uint32_t col); // getSyndrome
static void printParity(const RV_t **parity, uint32_t row, uint32_t col);  // printParity(this->parity, this->MM, this->P_col);
static void printcError(const RV_t *error, uint32_t len);                            // printError(error, P_col);
static void printError(const uint8_t *error, uint32_t len);                        // printError(error, NN);
static void printSyndrome(const uint8_t *syndrome, uint32_t len);    // printSyndrome(syndrome, MM)

// Calculation
static uint32_t decideOSDW(uint32_t x);
static void decideRelyGF2(Index_Rely *permu, const double **prob, const uint8_t *lastFor);
static void decideRelyGF2_None(Index_Rely *permu, const double **prob);
static void decideRelyGF4(Index_Rely *permu, const double **prob, const uint8_t *lastFor);
static void decideRelyGF4_Entropy(Index_Rely *permu, const double **prob);
static void decideRelyGF4_Max(Index_Rely *permu, const double **prob);
static void decideRelyGF4_None(Index_Rely *permu, const double **prob);
static uint32_t decideUnreliGF2(Index_Rely *permu, const double **prob, const uint8_t *lastFor, RV_t *unreliMask, uint32_t lastForLimit, double reliLimit);
static uint32_t decideUnreliGF2_None(Index_Rely *permu, const double **prob, RV_t *unreliMask, double reliLimit);
static uint32_t decideUnreliGF4(Index_Rely *permu, const double **prob, const uint8_t *lastFor, RV_t *unreliMask, uint32_t lastForLimit, double reliLimit);
static uint32_t decideUnreliGF4_Entropy(Index_Rely *permu, const double **prob, RV_t *unreliMask, uint32_t ETPLimit, double reliLimit);
static uint32_t decideUnreliGF4_Max(Index_Rely *permu, const double **prob, RV_t *unreliMask, uint32_t maxLimit, double reliLimit);
static uint32_t decideUnreliGF4_None(Index_Rely *permu, const double **prob, RV_t *unreliMask, double reliLimit);
static void hardDecisionGF2(const double **prob, RV_t *err, uint32_t NN, uint32_t P_col);
static void hardDecisionGF4(const double **prob, RV_t *err, uint32_t NN, uint32_t P_col);
static uint32_t decideRelyCheck( RV_t **parity,  RV_t *unreliMask, int32_t *rowPivotPlace, uint32_t *unreliCheck);   // return num of unreliable check
static uint32_t decideUnreliP_col(const RV_t *unreliMask, uint32_t *unreliP_col);
static int32_t reliCheckSAT(RV_t **parity, RV_t *error, const uint8_t *syndrome, const uint32_t *rowPivotPlace);
static void updateUnreliSynd( RV_t **parity, RV_t *error, uint8_t *syndrome,  const RV_t *unreliMask, const uint32_t *unreliCheck, uint32_t numOfUnreliCheck);
static void GJonParity(OSD *osdDec, int32_t all);
static uint32_t parityHasSol(const uint8_t *syndrome, const uint32_t *unreliCheck, const int32_t *rowPivotPlace, uint32_t numOfUnreliCheck);
static int32_t genBit2RowPos(RV_t **parity, uint32_t **bit2RowPos, int32_t *bitColWeight, uint32_t *sneakyBit, uint32_t numOfSneaky, uint32_t *unreliCheck, uint32_t numOfUnreliCheck, int32_t dist);
static void updateMostPossibleError(RV_t *mostPossibleError, RV_t *error, uint8_t *syndrome, uint32_t *unreliCheck, int32_t *rowPivotPlace);
static uint32_t OSD_prepare2Flip(OSD *osdDec);
static void OSD_w_recursive(OSD *osdDec, int32_t iter, int32_t start, int32_t weight);
static void OSD_lambda_recursive(OSD *osdDec, int32_t start, int32_t weight);
static void syndromeXOR1Col(const RV_t **parity, uint32_t pos, uint32_t shift, uint8_t *syndrome, const uint32_t *unreliCheck);
static int32_t syndromeXOR1Col_wPos(const int *bit2RowPos, uint8_t *syndrome, uint32_t bitColWeight, int8_t *rowPivotStatus, uint32_t *rowPivotBrother, int32_t weight, int32_t updWeight);
static int compareQua(const void *a, const void *b);
static int compareBin(const void *a, const void *b);

// Public
OSD* OSD_new();
void initOSD(OSD *osdDec);
void free_OSD(OSD *osdDec);
void load_A_OSD(OSD *osdDec, FILE *fp);
//void getSyndromeFromIntactH(RV_t *error, uint8_t *syndrome, uint32_t MM, uint32_t P_col) {getSyndrome(error, parity_intact, syndrome, MM, P_col);}
uint8_t* post_decodeOSDw(OSD *osdDec, const uint8_t* synd, const double **prob, const uint8_t *lastFor);
uint8_t* post_decodeOSDfull(OSD *osdDec, const uint8_t* synd, const double **prob, const uint8_t *lastFor, uint32_t max_iter);


uint8_t numOf1(RV_t num);
void compressError(const uint8_t *error, RV_t *cError, int NN, int P_col);
void restoreError(const RV_t *cError, uint8_t *error, int BB);

int compareT(const void *a, const void *b);

// Public OsdRecord
void initOsdRecord(OsdRecord *osdRec);

#endif
