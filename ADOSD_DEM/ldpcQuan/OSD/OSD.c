#include "OSD.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#define getQPos(x) (x>>QShift)
#define getBPos(x) (x>>BShift)
#define getQShift(x) ((QMask-(x&QMask))<<1)
#define getBShift(x)  (BMask-(x&BMask))
#define max(a, b) ((a>b)?a:b)
#define min(a, b) ((a<b)?a:b)
#define abs(a) ((a>=0)?a:-a)

#define DBG_PRT_OSD 0
static const int dbg_prtOSD = DBG_PRT_OSD;
static const int dbg_prt=0;
static uint32_t  acRow, acCol, acBit, maxTrial, maxLayer;
static uint32_t NN, MM, P_col, NNbit;
static RV_t *tempError;
#if CAL_TIME
static uint64_t tStart, tEnd, diff;
#endif // CAL_TIME
static uint32_t numOfTrialComb;
static struct timeval currentTime;


uint8_t numOf1Arr[256];

int compareT(const void *a, const void *b)
{
      int c = *(int *)a;
      int d = *(int *)b;
      if(c < d) {return -1;}                //return  -1 means a < b
      else if (c == d) {return 0;}          //return   0 means a = b
      else return 1;                         //return  1 means a>b
}

OSD* OSD_new(){
    OSD temp, *osdDec = &temp;
    initOSD(osdDec);
    return osdDec;
}

void initOSD(OSD *this){
    this->NN = N;
    this->NNbit = this->NN<<1;
    this->MM = M;
    this->P_col = this->NN >> QShift;    // the number of col (in compress mode) do parity should have
    if(this->NN & QMask)
        this->P_col ++;
    this->RankH = N-K;
    this->osdw = OSDW;
    this->Dist = D;


    // ==== initialize counters ====
    this->cnt_enter_osd = 0;
    this->cnt_stage1_fail = 0;
    this->cnt_stage2_fail = 0;
    this->cnt_osd0 = 0;


    NN = this->NN, MM = this->MM, P_col = this->P_col, NNbit = this->NNbit;

    this->parity = calloc(this->MM, sizeof( *this->parity));
    this->parity_intact = calloc(this->MM, sizeof( *this->parity_intact));
    for(uint32_t m = 0; m<this->MM; m++){
        this->parity[m] = calloc(this->P_col, sizeof(*this->parity[m]));
        this->parity_intact[m] = calloc(this->P_col, sizeof(*this->parity_intact[m]));
    }
    this->syndrome = calloc(this->MM, sizeof(*this->syndrome));
    this->tempSyndrome = calloc(this->MM, sizeof(*this->tempSyndrome));

    this->permu = calloc(this->NNbit, sizeof(*this->permu));
    this->unreliMask = calloc(this->P_col, sizeof(*this->unreliMask));
    this->unreliCheck = calloc(this->MM, sizeof(*this->unreliCheck));
    this->unreliP_col = calloc(this->P_col, sizeof(*this->unreliP_col));
    this->pivotMask = calloc(this->P_col, sizeof(*this->pivotMask));
    this->rowPivotPlace = calloc(this->MM, sizeof(*this->rowPivotPlace));
    this->rowPivotStatus = calloc(this->MM, sizeof(*this->rowPivotStatus));
    this->rowPivotBrother = calloc(this->MM, sizeof(*this->rowPivotBrother));
    this->pivotBrother = calloc(this->NNbit, sizeof(*this->pivotBrother));
    this->matchBrother = calloc(this->NNbit, sizeof(*this->matchBrother));
    this->sneakyBit = calloc(this->NNbit, sizeof(*this->sneakyBit));
    this->bit2RowPos = calloc(this->NNbit, sizeof( *this->bit2RowPos));
    for(uint32_t n = 0; n<this->NNbit; n++){
        this->bit2RowPos[n] = calloc(this->MM, sizeof(*this->bit2RowPos[n]));
    }
    this->bitColWeight = calloc(this->NNbit, sizeof(*this->bitColWeight));


    this->error = calloc(this->P_col, sizeof(*this->error));
    this->mostPossibleError = calloc(this->P_col, sizeof(*this->mostPossibleError));
    this->returnError = calloc(this->NN, sizeof(*this->returnError));

    for(int i = 0; i<this->MM; i++)
        this->unreliCheck[i] = i;
    for(int i =0; i<this->P_col; i++)
        this->unreliP_col[i] = i;
    memset(this->unreliMask, 255, P_col*sizeof(*this->unreliMask));


    // We want to calculate the number of error we find and calculate the weight
    maxTrial = 1+(N+K)+(N+K)*(N+K-1)/2;
    numOfTrialComb =maxTrial; // the number of trial we use in osd2

    uint32_t squ = maxTrial;
    maxLayer = 0;
    while(squ >>=1)
        maxLayer++;
        if(OSDW ==-2){
        if(ADOSDw == 0)
            printf("ADOSD-0 \n");
        else{
            printf("Maximum Trial = %d\n", maxTrial);
            printf("Maximum layer = %d (osdw)\n",maxLayer);
        }
        }

    // Global thing
    int num;
    for(int i = 0; i<256; i++){
        num = i;
        num ^= num>>4;
        num ^= num>>2;
        num ^= num>>1;
        numOf1Arr[i] = (uint8_t)(num&1);
    }
    acRow = this->MM, acCol = this->P_col, acBit = this->NNbit;
    tempError = calloc(this->P_col, sizeof(*tempError));
}

void initOsdRecord(OsdRecord *this){
    this->nStep1 = 0, this->nStep2 = 0;
    uint32_t NNbit = N<<1;
    # if STORE_DISTR
    for(int i = 1; i<=M; i++){
        for(int j = 1; j<=NNbit; j++){
            this -> newMatSizeDistr[i][j] = 0;
            this -> newMatSizeDistr_osd0[i][j] = 0;
        }
    }
    #endif // STORE_DISTR
}

void free_OSD(OSD *this){
    static uint32_t m;

    for(m = 0; m<this->MM; m++){
        free(this->parity[m]);
        free(this->parity_intact[m]);
    }

    for(m = 0; m<this->NNbit; m++)
        free(this->bit2RowPos[m]);
    free(this->parity);
    free(this->parity_intact);
    free(this->bit2RowPos);
    free(this->syndrome);
    free(this->tempSyndrome);
    free(this->permu);
    free(this->unreliMask);
    free(this->unreliCheck);
    free(this->unreliP_col);
    free(this->pivotMask);
    free(this->rowPivotPlace);
    free(this->rowPivotStatus);
    free(this->rowPivotBrother);
    free(this->pivotBrother);
    free(this->matchBrother);
    free(this->sneakyBit);
    free(this->bitColWeight);
    free(this->error);
    free(this->mostPossibleError);
    free(this->returnError);
    free(tempError);
}




/*

ADOSD with full decoding history

*/
uint8_t* post_decodeOSDfull(OSD *this, const uint8_t* synd, const double **probability, const uint8_t *lastFor, uint32_t max_iter){
    static uint32_t fail, weight;
    static uint32_t n, m, p_col, temp;

    fail = 0;
    this -> osd0orNot = 0;
    // record time and osd rate
    this -> status = UsePosd;
//    # if CAL_TIME
//    gettimeofday(&currentTime, NULL);
//    tStart = currentTime.tv_sec * (uint64_t)1e6 + currentTime.tv_usec;
//    #endif // CAL_TIME

    // ==== copy the target syndrome + parity =======
    memcpy(this->syndrome, synd, MM*sizeof(*this->syndrome));
    initParity(this->parity, this->parity_intact);
//    #if DBG_PRT_OSD
//        //printf("Initial parity - "); printParity(this->parity, MM, P_col);
//        //if(dbg_prtOSD){printf("Target syndrome - "); printSyndrome(this->syndrome, MM);}
//    #endif // DBG_PRT_OSD

    // === Step1: decide unreliable part + hard decision ==

        hardDecisionGF4(probability, this->error, NN, P_col);
        #if SORTCOMPARE == LASTFORSORT
        this->numOfUnreliBit = decideUnreliGF4(this->permu, probability, lastFor, this->unreliMask, max_iter, RelThr);    // gen Permu
//        #elif SORTCOMPARE == ENTROPYSORT
//        this->numOfUnreliBit = decideUnreliGF4_Entropy(this->permu, probability, this->unreliMask, 1-RelThr, RelThr);    // gen Permu
//        #elif SORTCOMPARE == MAXSORT
//        this->numOfUnreliBit = decideUnreliGF4_Max(this->permu, probability, this->unreliMask, RelThr, RelThr);    // gen Permu
//        #elif SORTCOMPARE == BINARYSORT
//        this->numOfUnreliBit = decideUnreliGF4_None(this->permu, probability,  this->unreliMask, RelThr);    // gen Permu
        #endif // SORTCOMPARE


	//#if DBG_PRT_OSD
//        printf("Error from pre-decoder - NN=%d P_col=%d Qshift=%d   ",NN,P_col,QShift);
//        printcError(this->error, P_col);
//        for(int i = 0; i<NN; i++){ printf("%d ", i%10);}printf("\n");
//        printf("Unreliable mask:\n");
//	    for(p_col = 0; p_col<P_col; p_col++){
//            for(int j = BMask; j>=0; j--)
//                printf("%d ", ((this->unreliMask[p_col]>>j)&1));
//            printf("\n");
//	    }
	//#endif // DBG_PRT_OSD

//	printf("number of sneaky=%d\n", this->numOfUnreliBit);

	//  === Step2: find the check with only reli bit =======
    this->numOfUnreliP_col = decideUnreliP_col(this->unreliMask, this->unreliP_col);
    this->numOfUnreliCheck =  decideRelyCheck(this->parity_intact, this->unreliMask, this->rowPivotPlace, this->unreliCheck);
    acRow = this->numOfUnreliCheck;
    acBit = this->numOfUnreliBit;
    acCol = this->numOfUnreliP_col;

    #if DBG_PRT_OSD && 1
        printf("Row pivot place:\n");
        for(int i=0; i<M; i++){
            printf("%d ", this->rowPivotPlace[i]);
        }
        printf("\n");
        printf("Unreliable check: \n");
        for(int i = 0; i<this->numOfUnreliCheck; i++){
            printf("%d ", this->unreliCheck[i]);
        }
        printf("\n");
        printf("Unreliable P_col: %d\n", this->numOfUnreliP_col);
        for(int i = 0; i<this->numOfUnreliP_col; i++){
            printf("%d ", this->unreliP_col[i]);
        }
        printf("\n");
    #endif // DBG_PRT_OSD

    // === Step3: Check reliable check is satisfy ========
    if(!reliCheckSAT(this->parity_intact, this->error, synd, this->rowPivotPlace)){
        fail = 1;
        this -> status = Stage1Fail;
        this->cnt_stage1_fail++;   //count
//        printf("stage1 fail \n");
//
    }
    else if (dbg_prtOSD){
        printf("success at first step!!!!\n");
    }

    // === Step4: permu + Get new syndrome =========
    if(fail == 0){
        #if SORTCOMPARE == BINARYSORT
                qsort(this->permu, this->numOfUnreliBit, sizeof(*this->permu), compareBin);
        #else
                qsort(this->permu, this->numOfUnreliBit, sizeof(*this->permu), compareQua);
        #endif // SORTCOMPARE
        updateUnreliSynd(this->parity_intact, this->error, this->syndrome, this->unreliMask, this->unreliCheck, this->numOfUnreliCheck);
        GJonParity(this, POSD_FLIPPING_ALL_BIT);


        if(!parityHasSol(this->syndrome, this->unreliCheck, this->rowPivotPlace, this->numOfUnreliCheck)){
            this-> status = Stage2Fail;
            fail = 1;
            this->cnt_stage2_fail++;   // count
//            printf("stage2 fail \n");

        }
        else if (dbg_prtOSD)
            printf("success at second step!!!!\n");

        # if DBG_PRT_OSD
            printf("Unreliable bit:\n");
            //for(n = 0; n<this->numOfUnreliBit; n++){
            //    printf("%4d %c ",(this->permu[n].index>>1), ((this->permu[n].index &1) ? 'X': 'Z'));
            //}

            for(n = 0; n<this->numOfUnreliBit; n++){
                printf("bit: %4d, qibit: %4d, %c error, with lastFor: %lf, reliability: %lf\n",
                       this->permu[n].index, (this->permu[n].index>>1), ((this->permu[n].index &1) ? 'X': 'Z'), (this->permu[n].reliQua), (this->permu[n].reliBin));
            }
            printf("\n");
            # if 0
            printf("After GJ parity - \n"); for(int i = 0; i<NN; i++){ printf("%d ", i%10);}printf("\n"); printParity(this->parity, MM, P_col);
            printSyndrome(this->syndrome, MM);

            printf("Row Pivot Place: \n");
            for(m= 0; m<MM; m++)
                printf("%d ", this->rowPivotPlace[m]);
            printf("\n");
            printf("Pivot Place: \n");
            for(p_col = 0; p_col<P_col; p_col++){
                for(int j = BMask; j>=0; j--)
                    printf("%d ", ((this->pivotMask[p_col]>>j)&1));
                printf("\n");
            }
            # endif
            printf("Sneaky bit:\n");
            for(int i = 0; i<this->numOfSneaky; i++)
                printf("%d ", this->sneakyBit[i]);
            printf("\n");
        # endif // DBG_PRT_OSD
    }

    // === Step 5: run through all possible =============
    /*
    determine the order of ADOSD
    */
    if(fail == 0){
        this -> osd0orNot = genBit2RowPos(this->parity, this->bit2RowPos, this->bitColWeight, this->sneakyBit, this->numOfSneaky, this->unreliCheck, this->numOfUnreliCheck, this->Dist);
        if(this -> osd0orNot){
            this->osdw = 0;
            this->cnt_osd0++;    // count where OSD-0 suffices
//            printf("OSD-0 suffices \n");
        }
        else{// determine the order of ADOSD
            if (ADOSDw == 0)
                this->osdw = 0;                 // force to order-0
            else
                this->osdw = decideOSDW(this->numOfSneaky);  // decide order based on sneaky bits

        }
        weight = OSD_prepare2Flip(this);

        this->cnt_enter_osd++;   // count


        OSD_w_recursive(this, 0, 0, weight);
        //this->osdw = min(maxLayer, this->numOfSneaky);
        //OSD_lambda_recursive(this, 0, weight);
    }
    else{
        acRow = MM;
        acBit = NNbit;
        acCol = P_col;
        for(int i = 0; i<MM; i++)
            this->unreliCheck[i] = i;
        for(int i =0; i<P_col; i++)
            this->unreliP_col[i] = i;
        memset(this->unreliMask, 255, P_col*sizeof(*this->unreliMask));
        this->osdw = 2;
        post_decodeOSDw(this, synd, probability, lastFor);
    }

    restoreError(this->mostPossibleError, this->returnError, NN);

    #if DBG_PRT_OSD
        printf("Final: Most possible error with weight: %d\n", this->bestErrorWeight);
        printf("Final error - "); printcError(this->mostPossibleError, P_col);

        getSyndrome(this->mostPossibleError, this->parity_intact, this->tempSyndrome, MM, P_col);
        (memcmp(this->tempSyndrome, synd, MM*sizeof(*synd))==0)? printf("Syndrome match\n"):printf("Syndrome fail\n");
        printf("Restore - ");
        printError(this->returnError, NN);
        printf("\n\n");
    #endif // DBG_PRT_OSD

    // record the value in this calculation
    #if CAL_TIME
    gettimeofday(&currentTime, NULL);
    tEnd = currentTime.tv_sec * (uint64_t)1e6 + currentTime.tv_usec;
    (tEnd<tStart)?(diff = 0):(diff = tEnd-tStart);
    this -> usedTime = diff;
    #endif // CAL_TIME

    return this->returnError;
}



/*

Regular OSD-w

*/
uint8_t* post_decodeOSDw(OSD *this, const uint8_t* synd, const double **probability, const uint8_t *lastFor){
    static uint32_t weight;
    static uint32_t n, m, p_col;

    this-> status = NormalOsd;
    this->cnt_enter_osd++;   // count

    // ====== reset parity and target syndrome ==========
    memcpy(this->syndrome, synd, MM*sizeof(*this->syndrome));
    initParity(this->parity, this->parity_intact);
    #if DBG_PRT_OSD
        //printf("Initial parity - "); printParity(this->parity, MM, P_col);
        //if(dbg_prtOSD){printf("Target syndrome - "); printSyndrome(this->syndrome, MM);}
    #endif // DBG_PRT_OSD

    // ==== doing hard decision and get the reliability ====
    #if USE_GF2_DEC
        hardDecisionGF2(probability, this->error, NN, P_col);
        # if SORTCOMPARE == LASTFOR
            decideRelyGF2(this->permu, probability, lastFor);    // gen Permu
        #else
            decideRelyGF2_None(this->permu, probability);    // gen Permu
        #endif // SORTCOMPARE
    #else
        hardDecisionGF4(probability, this->error, NN, P_col);
        #if SORTCOMPARE == LASTFORSORT
            decideRelyGF4(this->permu, probability, lastFor);    // gen Permu
        #elif SORTCOMPARE == ENTROPYSORT
            decideRelyGF4_Entropy(this->permu, probability);
        #elif SORTCOMPARE == MAXSORT
            decideRelyGF4_Max(this->permu, probability);
        #elif SORTCOMPARE == BINARYSORT
            decideRelyGF4_None(this->permu, probability);
        #endif // SORTCOMPARE
	#endif // USE_GF2_DEC

	// ============== sort ================
	# if SORTCOMPARE == BINARYSORT
        qsort(this->permu, NNbit, sizeof(*this->permu), compareBin);
	# else
        qsort(this->permu, NNbit, sizeof(*this->permu), compareQua);
    #endif

    #if DBG_PRT_OSD
        printf("Error from pre-decoder - ");
        printcError(this->error, P_col);
for(int i = 0; i<NN; i++){ printf("%d ", i%10);}printf("\n");

        printf("After decide Rely, the rearrangment: \n");
        for(n = 0; n<NNbit; n++){
            printf("bit: %4d, qibit: %4d, %c error, with lastFor: %lf, reliability: %lf\n",
                   this->permu[n].index, (this->permu[n].index>>1), ((this->permu[n].index &1) ? 'X': 'Z'), (this->permu[n].reliQua), (this->permu[n].reliBin));
        }
    #endif


    // ========== doing Gaussian elimination ======
    /* just show where the pivot (bit form) in each row.
       It will apply to all -1 first, and then record the pivot when finding it. */
    for(m = 0; m<MM; m++)
        this->rowPivotPlace[m] = -1;
    GJonParity(this, 1);
    #if DBG_PRT_OSD
        //printf("After GJ parity - "); printParity(this->parity, MM, P_col);
        //printf("After GJ syndrome - "); printSyndrome(this->syndrome, MM);
        printf("Row Pivot Place: \n");
        for(m= 0; m<MM; m++)
            printf("%d ", this->rowPivotPlace[m]);
        printf("\n");
        printf("Pivot Place: \n");
        for(p_col = 0; p_col<P_col; p_col++){
            for(int j = BMask; j>=0; j--)
                printf("%d ", ((this->pivotMask[p_col]>>j)&1));
            printf("\n");
        }
        printf("Sneaky bit: %d\n", this->numOfSneaky);
        for(n = 0; n<this->numOfSneaky; n++){
            printf("%u ", this->sneakyBit[n]);
        }
        printf("\n");
    #endif // DBG_PRT_OSD

    // ========= doing OSDw =================
    genBit2RowPos(this->parity, this->bit2RowPos, this->bitColWeight, this->sneakyBit, this->numOfSneaky, this->unreliCheck, MM, this->Dist);
    weight = OSD_prepare2Flip(this);
    OSD_w_recursive(this, 0, 0, weight);
    restoreError(this->mostPossibleError, this->returnError, NN);

    #if DBG_PRT_OSD
        printf("Final: Most possible error with weight: %d\n", this->bestErrorWeight);
        printf("Final error - "); printcError(this->mostPossibleError, P_col);

        getSyndrome(this->mostPossibleError, this->parity_intact, this->tempSyndrome, MM, P_col);
        (memcmp(this->tempSyndrome, synd, MM*sizeof(*synd))==0)? printf("Syndrome match\n"):printf("Syndrome fail\n");
        printf("Restore - ");
        printError(this->returnError, NN);
        printf("\n\n");
    #endif
    return this->returnError;
}

int compareQua(const void *a, const void *b){
    const Index_Rely *A = (Index_Rely*)a, *B = (Index_Rely*)b;

    if((A->reliQua) > (B->reliQua))
        return 1;
    else if((A->reliQua) < (B->reliQua))
        return -1;
    else{
        if ( (A->reliBin) > (B->reliBin))
            return 1;
        else if ((A->reliBin) < (B->reliBin))
            return -1;
    }
    return 0;

}

int compareBin(const void *a, const void *b){
    const Index_Rely *A = (Index_Rely*)a, *B = (Index_Rely*)b;

    if ( (A->reliBin) > (B->reliBin))
        return 1;
    else if ((A->reliBin) < (B->reliBin))
        return -1;
    else
        return 0;
}



uint32_t decideOSDW(uint32_t x){
    /*
     * Determine the effective OSD order (osdw) based on a complexity budget.
     *
     * x         : number of unreliable (sneaky) bits
     * maxLayer  : maximum allowed OSD order (upper bound)
     * maxTrial  : maximum allowed number of tested error patterns
     *
     * The function finds the largest OSD order w such that
     *
     *      sum_{k=0}^{w} C(x, k) <= maxTrial ,
     *
     * i.e., the total number of tested error patterns does not exceed
     * the prescribed complexity budget.
     *
     * Return value:
     *      Effective OSD order (osdw).
     */

    static uint64_t result, i, sum, comb;

    // If the number of unreliable bits is small,
    // we can afford full enumeration up to order x.
    if(x <= maxLayer){
        return x;
    }

    // sum  : cumulative number of tested patterns
    // comb : C(x, i), number of patterns at OSD order i
    sum = 1;     // C(x,0)
    comb = 1;    // C(x,0)

    for(i = 1; i <= x; i++){
        // Compute binomial coefficient C(x, i) iteratively
        comb = comb * (x + 1 - i) / i;

        // If adding OSD-i patterns exceeds the complexity budget,
        // stop and return the previous order (i-1)
        if(sum + comb > maxTrial){
            return i - 1;
        }

        // Accumulate total number of tested patterns
        sum += comb;
    }

    // In theory, this case is reached only if all orders up to x are affordable
    return x;
}

void decideRelyGF4(Index_Rely *permu, const double **prob, const uint8_t *lastFor){
    /* Since I already change the order of X, Z in parity check matrix
    So, the order of parity and error are all correct
    The nth bit error is correspond to nth column of parity*/
    static uint32_t nbit, n;

    for(nbit = 0, n = 0; n<NN; nbit++, n++){
        // Reliability of nth qubit of Z error
        permu[nbit].index = nbit;
        permu[nbit].reliQua = lastFor[n];
        permu[nbit].reliBin = max(prob[n][0]+prob[n][1], prob[n][2]+prob[n][3]);

        nbit++;
        // Reliability of nth qubit of X error
        permu[nbit].index = nbit;
        permu[nbit].reliQua = lastFor[n];
        permu[nbit].reliBin = max(prob[n][0]+prob[n][2], prob[n][1]+prob[n][3]);
    }
}

//void decideRelyGF4_Entropy(Index_Rely *permu, const double **prob){
//    static uint32_t nbit, n;
//    static double reliQ;
//
//    for(nbit = 0, n = 0; n<NN; nbit++, n++){
//        reliQ = prob[n][0]*log2(prob[n][0]) + prob[n][1]*log2(prob[n][1]) + prob[n][2]*log2(prob[n][2]) + prob[n][3]*log2(prob[n][3]);
//        // Reliability of nth qubit of Z error
//        permu[nbit].index = nbit;
//        permu[nbit].reliQua = reliQ;
//        permu[nbit].reliBin = max(prob[n][0]+prob[n][1], prob[n][2]+prob[n][3]);
//
//        nbit++;
//        // Reliability of nth qubit of X error
//        permu[nbit].index = nbit;
//        permu[nbit].reliQua = reliQ;
//        permu[nbit].reliBin = max(prob[n][0]+prob[n][2], prob[n][1]+prob[n][3]);
//    }
//}

//void decideRelyGF4_Max(Index_Rely *permu, const double **prob){
//    static uint32_t nbit, n;
//    static double reliQ;
//
//    for(nbit = 0, n = 0; n<NN; nbit++, n++){
//        reliQ = max(max(max(prob[n][0], prob[n][1]), prob[n][2]), prob[n][3]);
//        // Reliability of nth qubit of Z error
//        permu[nbit].index = nbit;
//        permu[nbit].reliQua = reliQ;
//        permu[nbit].reliBin = max(prob[n][0]+prob[n][1], prob[n][2]+prob[n][3]);
//
//        nbit++;
//        // Reliability of nth qubit of X error
//        permu[nbit].index = nbit;
//        permu[nbit].reliQua = reliQ;
//        permu[nbit].reliBin = max(prob[n][0]+prob[n][2], prob[n][1]+prob[n][3]);
//    }
//}
//
//void decideRelyGF4_None(Index_Rely *permu, const double **prob){
//    static uint32_t nbit, n;
//
//    for(nbit = 0, n = 0; n<NN; nbit++, n++){
//        // Reliability of nth qubit of Z error
//        permu[nbit].index = nbit;
//        permu[nbit].reliBin = max(prob[n][0]+prob[n][1], prob[n][2]+prob[n][3]);
//
//        nbit++;
//        // Reliability of nth qubit of X error
//        permu[nbit].index = nbit;
//        permu[nbit].reliBin = max(prob[n][0]+prob[n][2], prob[n][1]+prob[n][3]);
//    }
//}

//void decideRelyGF2(Index_Rely *permu, const double **prob, const uint8_t *lastFor){
//    /* Since I already change the order of X, Z in parity check matrix
//    So, the order of parity and error are all correct
//    The nth bit error is correspond to nth column of parity*/
//    static uint32_t prob_pos, n, nbit;
//    for(n = 0, nbit=0; n<NN; nbit++, n++){
//         Reliability of nth qubit of Z error
//        prob_pos = n+NN;
//        permu[nbit].index = nbit;
//        permu[nbit].reliQua = lastFor[prob_pos];
//        permu[nbit].reliBin = max(prob[prob_pos][0], prob[prob_pos][1]);
//
//        nbit++;
//         Reliability of nth qubit of X error
//        prob_pos = n;
//        permu[nbit].index = nbit;
//        permu[nbit].reliQua = lastFor[prob_pos];
//        permu[nbit].reliBin = max(prob[prob_pos][0], prob[prob_pos][1]);
//    }
//}
//
//void decideRelyGF2_None(Index_Rely *permu, const double **prob){
//    static uint32_t prob_pos, n, nbit;
//    for(n = 0, nbit=0; n<NN; nbit++, n++){
//         Reliability of nth qubit of Z error
//        prob_pos = n+NN;
//        permu[nbit].index = nbit;
//        permu[nbit].reliBin = max(prob[prob_pos][0], prob[prob_pos][1]);
//
//        nbit++;
//         Reliability of nth qubit of X error
//        prob_pos = n;
//        permu[nbit].index = nbit;
//        permu[nbit].reliBin = max(prob[prob_pos][0], prob[prob_pos][1]);
//    }
//}

uint32_t decideUnreliGF4(Index_Rely *permu, const double **prob, const uint8_t *lastFor, RV_t *unreliMask, uint32_t lastForLimit, double reliLimit){
    /* Since I already change the order of X, Z in parity check matrix
    So, the order of parity and error are all correct
    The nth bit error is correspond to nth column of parity*/
    static uint32_t nbit, n, tempLast, count, pos, shift;
    static double tempProb;

    memset(unreliMask, 0, P_col*sizeof(*unreliMask));
    count = 0;
    for(nbit = 0, n = 0; n<NN; nbit++, n++){
//        // Reliability of nth qubit of Z error =====================
//        tempProb = max(prob[n][0]+prob[n][1], prob[n][2]+prob[n][3]);
//        if(lastFor[n] < lastForLimit || tempProb <= reliLimit){
//            permu[count].index = nbit;
//            permu[count].reliQua = lastFor[n];
//            permu[count].reliBin = tempProb;
//            pos = getBPos(nbit);
//            shift = getBShift(nbit);
//            unreliMask[pos] |= (((RV_t)1)<<shift);
//            count++;
//        }

        nbit++;
        // Reliability of nth qubit of X error =====================
        tempProb = max(prob[n][0]+prob[n][2], prob[n][1]+prob[n][3]);

//        printf("lastfor[n]=%d tempProb=%.6e reliLimit=%.6e\n",lastFor[n], tempProb, reliLimit);

//        if(lastFor[n] < lastForLimit || tempProb <= reliLimit){
//  LastRun  is not used now
            if( tempProb <= reliLimit){
            permu[count].index = nbit;
            permu[count].reliQua = lastFor[n];
            permu[count].reliBin = tempProb;
            pos = getBPos(nbit);
            shift = getBShift(nbit);
            unreliMask[pos] |= (((RV_t)1)<<shift);
            count++;
        }
    }

    return count;
}

// ETPLimit could be 0.00005
uint32_t decideUnreliGF4_Entropy(Index_Rely *permu, const double **prob, RV_t *unreliMask, uint32_t ETPLimit, double reliLimit){
    static uint32_t nbit, n, tempLast, count, pos, shift;
    static double tempProb, tempETP;

    ETPLimit = -ETPLimit;
    memset(unreliMask, 0, P_col*sizeof(*unreliMask));
    count = 0;
    for(nbit = 0, n = 0; n<NN; nbit++, n++){
        tempETP = prob[n][0]*log2(prob[n][0]) + prob[n][1]*log2(prob[n][1]) + prob[n][2]*log2(prob[n][2]) + prob[n][3]*log2(prob[n][3]);
        // Reliability of nth qubit of Z error =====================
        tempProb = max(prob[n][0]+prob[n][1], prob[n][2]+prob[n][3]);
        if(tempETP < ETPLimit || tempProb <= reliLimit){
            permu[count].index = nbit;
            permu[count].reliQua = tempETP;
            permu[count].reliBin = tempProb;
            pos = getBPos(nbit);
            shift = getBShift(nbit);
            unreliMask[pos] |= (((RV_t)1)<<shift);
            count++;
        }

        nbit++;
        // Reliability of nth qubit of X error =====================
        tempProb = max(prob[n][0]+prob[n][2], prob[n][1]+prob[n][3]);
        if(tempETP < ETPLimit  || tempProb <= reliLimit){
            permu[count].index = nbit;
            permu[count].reliQua = tempETP;
            permu[count].reliBin = tempProb;
            pos = getBPos(nbit);
            shift = getBShift(nbit);
            unreliMask[pos] |= (((RV_t)1)<<shift);
            count++;
        }
    }

    return count;
}

//uint32_t decideUnreliGF4_Max(Index_Rely *permu, const double **prob, RV_t *unreliMask, uint32_t maxLimit, double reliLimit){
//    static uint32_t nbit, n, tempLast, count, pos, shift;
//    static double tempProb, tempMax;
//
//    memset(unreliMask, 0, P_col*sizeof(*unreliMask));
//    count = 0;
//    for(nbit = 0, n = 0; n<NN; nbit++, n++){
//        tempMax =  max(max(max(prob[n][0], prob[n][1]), prob[n][2]), prob[n][3]);
//        // Reliability of nth qubit of Z error =====================
//        tempProb = max(prob[n][0]+prob[n][1], prob[n][2]+prob[n][3]);
//        if(tempMax < maxLimit || tempProb <= reliLimit){
//            permu[count].index = nbit;
//            permu[count].reliQua = tempMax;
//            permu[count].reliBin = tempProb;
//            pos = getBPos(nbit);
//            shift = getBShift(nbit);
//            unreliMask[pos] |= (((RV_t)1)<<shift);
//            count++;
//        }
//
//        nbit++;
//        // Reliability of nth qubit of X error =====================
//        tempProb = max(prob[n][0]+prob[n][2], prob[n][1]+prob[n][3]);
//        if(tempMax < maxLimit || tempProb <= reliLimit){
//            permu[count].index = nbit;
//            permu[count].reliQua =tempMax;
//            permu[count].reliBin = tempProb;
//            pos = getBPos(nbit);
//            shift = getBShift(nbit);
//            unreliMask[pos] |= (((RV_t)1)<<shift);
//            count++;
//        }
//    }
//
//    return count;
//}
//
//uint32_t decideUnreliGF4_None(Index_Rely *permu, const double **prob, RV_t *unreliMask, double reliLimit){
//    static uint32_t nbit, n, tempLast, count, pos, shift;
//    static double tempProb;
//
//    memset(unreliMask, 0, P_col*sizeof(*unreliMask));
//    count = 0;
//    for(nbit = 0, n = 0; n<NN; nbit++, n++){
//        // Reliability of nth qubit of Z error =====================
//        tempProb = max(prob[n][0]+prob[n][1], prob[n][2]+prob[n][3]);
//        if( tempProb <= reliLimit){
//            permu[count].index = nbit;
//            permu[count].reliBin = tempProb;
//            pos = getBPos(nbit);
//            shift = getBShift(nbit);
//            unreliMask[pos] |= (((RV_t)1)<<shift);
//            count++;
//        }
//
//        nbit++;
//        // Reliability of nth qubit of X error =====================
//        tempProb = max(prob[n][0]+prob[n][2], prob[n][1]+prob[n][3]);
//        if(tempProb <= reliLimit){
//            permu[count].index = nbit;
//            permu[count].reliBin = tempProb;
//            pos = getBPos(nbit);
//            shift = getBShift(nbit);
//            unreliMask[pos] |= (((RV_t)1)<<shift);
//            count++;
//        }
//    }
//
//    return count;
//}

//uint32_t decideUnreliGF2(Index_Rely *permu, const double **prob, const uint8_t *lastFor, RV_t *unreliMask, uint32_t lastForLimit, double reliLimit){
//    /* Since I already change the order of X, Z in parity check matrix
//    So, the order of parity and error are all correct
//    The nth bit error is correspond to nth column of parity*/
//    static uint32_t prob_pos, n, nbit, tempLast, count, pos, shift;
//    static double tempProb;
//
//    memset(unreliMask, 0, P_col*sizeof(*unreliMask));
//    count = 0;
//    for(n = 0, nbit=0; n<NN; nbit++, n++){
//        // Reliability of nth qubit of Z error
//        prob_pos = n+NN;
//        tempProb = max(prob[prob_pos][0], prob[prob_pos][1]);
//        if(lastFor[prob_pos] <= lastForLimit || tempProb <= reliLimit){
//            permu[count].index = nbit;
//            permu[count].reliQua = lastFor[prob_pos];
//            permu[count].reliBin = tempProb;
//            pos = getBPos(nbit);
//            shift = getBShift(nbit);
//            unreliMask[pos] |= (((RV_t)1)<<shift);
//            count++;
//        }
//
//        nbit++;
//        // Reliability of nth qubit of X error
//        prob_pos = n;
//        tempProb = max(prob[prob_pos][0], prob[prob_pos][1]);
//        if(lastFor[prob_pos] <= lastForLimit || tempProb <= reliLimit){
//            permu[count].index = nbit;
//            permu[count].reliQua = lastFor[prob_pos];
//            permu[count].reliBin = tempProb;
//            pos = getBPos(nbit);
//            shift = getBShift(nbit);
//            unreliMask[pos] |= (((RV_t)1)<<shift);
//            count++;
//        }
//    }
//    return count;
//}

//uint32_t decideUnreliGF2_None(Index_Rely *permu, const double **prob, RV_t *unreliMask, double reliLimit){
//    /* Since I already change the order of X, Z in parity check matrix
//    So, the order of parity and error are all correct
//    The nth bit error is correspond to nth column of parity*/
//    static uint32_t prob_pos, n, nbit, tempLast, count, pos, shift;
//    static double tempProb;
//
//    memset(unreliMask, 0, P_col*sizeof(*unreliMask));
//    count = 0;
//    for(n = 0, nbit=0; n<NN; nbit++, n++){
//        // Reliability of nth qubit of Z error
//        prob_pos = n+NN;
//        tempProb = max(prob[prob_pos][0], prob[prob_pos][1]);
//        if( tempProb <= reliLimit){
//            permu[count].index = nbit;
//            permu[count].reliBin = tempProb;
//            pos = getBPos(nbit);
//            shift = getBShift(nbit);
//            unreliMask[pos] |= (((RV_t)1)<<shift);
//            count++;
//        }
//
//        nbit++;
//        // Reliability of nth qubit of X error
//        prob_pos = n;
//        tempProb = max(prob[prob_pos][0], prob[prob_pos][1]);
//        if(tempProb <= reliLimit){
//            permu[count].index = nbit;
//            permu[count].reliBin = tempProb;
//            pos = getBPos(nbit);
//            shift = getBShift(nbit);
//            unreliMask[pos] |= (((RV_t)1)<<shift);
//            count++;
//        }
//    }
//    return count;
//}

//void hardDecisionGF2(const double **prob, RV_t *error, uint32_t NN, uint32_t P_col){ // len is code length NN
//    static uint32_t pos, shift, n;
//    memset(error, 0, P_col*sizeof(*error));
//
//    for(n = 0; n<NN; n++){
//        if(prob[n][1] > 0.5){
//            pos = getQPos(n);
//            shift = getQShift(n);
//            error[pos] |= ((RV_t)1)<<shift;
//        }
//        if(prob[n+NN][1] > 0.5){
//            pos = getQPos(n);
//            shift = getQShift(n);
//            error[pos] |= ((RV_t)2)<<shift;
//        }
//    }
//}

void hardDecisionGF4(const double **prob, RV_t *error, uint32_t NN, uint32_t P_col){
    static uint32_t n, pos, shift;
    static RV_t index;
    memset(error, 0, P_col*sizeof(*error));

    for(n = 0; n<NN; n++){
        index = 0;
        for(int i = 1; i < 4; i++){
            if(prob[n][i] > prob[n][index])
                index = i;
        }
        pos = getQPos(n);
        shift = getQShift(n);
        error[pos] |= index<<shift;
    }
}

uint32_t decideUnreliP_col(const RV_t *unreliMask, uint32_t *unreliP_col){
    static uint32_t count, n;
    count = 0;

    for(n = 0; n<P_col; n++){
        if((unreliMask[n])!=(RV_t)0){
            unreliP_col[count++] = n;
        }
    }
    return count;
}

uint32_t decideRelyCheck( RV_t **parity, RV_t *unreliMask, int32_t *rowPivotPlace, uint32_t *unreliCheck){
    static uint32_t m, p_col, count;
    count = 0;
    for(m = 0; m<MM; m++){
        rowPivotPlace[m] = -2;
        for(p_col = 0; p_col<P_col; p_col++){
            if((parity[m][p_col] & unreliMask[p_col]) != 0){
                rowPivotPlace[m] = -1;
                unreliCheck[count++] = m;
                break;
            }
        }
    }
    return count;
}

int32_t reliCheckSAT( RV_t **parity,  RV_t *error, const uint8_t *syndrome, const uint32_t *rowPivotPlace){
    static uint32_t m, p_col;
    static RV_t result;
    for(m = 0; m<MM; m++){
        if(rowPivotPlace[m] == -2){
            result = 0;
            for(p_col = 0; p_col<P_col; p_col++){
                result ^= error[p_col] & parity[m][p_col];
            }
            if(numOf1(result) != syndrome[m]){
                return 0;
            }
        }
    }
    return 1;
}

void updateUnreliSynd( RV_t **parity,  RV_t *error, uint8_t *syndrome,  const RV_t *unreliMask, const uint32_t *unreliCheck, uint32_t numOfUnreliCheck){
    static uint32_t m, p_col, i;
    static RV_t result;
    for(p_col = 0; p_col<P_col; p_col++){
        tempError[p_col] = (~unreliMask[p_col]) & error[p_col];
    }

    for(i = 0; i<numOfUnreliCheck; i++){
        m = unreliCheck[i];
        result = 0;
        for(p_col = 0; p_col<P_col; p_col++){
            result ^= tempError[p_col] & parity[m][p_col];
        }
        syndrome[m] ^= (uint8_t)numOf1(result);
    }
}

// use pemu to do Gaussian on parity
void  GJonParity(OSD *this, int32_t all){
    static uint32_t count , pivot, upbound, pos, shift, n, m, i, j, jj, k, updCheck, setPivot, numOfSneaky, bufferCount;
    static uint32_t buffer[N];
    static RV_t pivotMaskInit , mask;
    pivotMaskInit = 0;
    pivotMaskInit = ~pivotMaskInit;

    /* pivotMask is to establish the mask when we generating the codeword at the end.
       some of the column will not join to calculate the syndrome.
       It will set to all 1, and when finding the pivot, set it to 0. */
    memset(this->pivotMask, 255, P_col*sizeof(*this->pivotMask));

    numOfSneaky = 0;

    count = 0, upbound = 0, bufferCount = 0;
    for(n = 0; n<acBit; n++){
        setPivot = 0;
        pivot = this->permu[n].index;
        pos = getBPos(pivot);
        shift = getBShift(pivot);
        mask = ((RV_t)1)<< shift;

        for(k = upbound; k<acRow; k++){
            m = this->unreliCheck[k];
            if((this->parity[m][pos] & mask) && (this->rowPivotPlace[m] == -1)){        // find pivot
                setPivot = 1;
                count++;
                this->rowPivotPlace[m] = pivot;

                this->pivotMask[pos] ^= mask;
                if(k == upbound){
                    upbound++;
                    while((upbound < acRow) && (this->rowPivotPlace[this->unreliCheck[upbound]] != -1))
                        upbound++;
                }

                // to elliminate entire col with only one pivot (GaussianGF2_1col)
                for(i = 0; i<acRow; i++){
                    updCheck = this->unreliCheck[i];
                    if((this->parity[updCheck][pos] & mask)!=0  && updCheck != m){
                        for(jj = 0; jj < acCol; jj++){
                            j = this->unreliP_col[jj];
                            this->parity[updCheck][j] ^= (this->parity[m][j]);
                        }
                        this->syndrome[updCheck] ^= this->syndrome[m];
                    }
                }
//printf("%d %d\n", m, pivot);
//for(int i = 0; i<NN; i++){ printf("%d ", i%10);}printf("\n");
//printParity(this->parity, this->MM, this->P_col);
                break;
            }
        }

        if(setPivot == 0){
            if ((this->unreliMask[pos]>>shift)&1){
                buffer[bufferCount] = pivot;
                bufferCount ++;
            }
        }else{
            for(i = 0; i< bufferCount; i++){
                this->sneakyBit[numOfSneaky] = buffer[i];
                numOfSneaky ++;
            }
            bufferCount = 0;
        }


        // if all == 1, we want to add all bits into the flipping list
        if(count == this->RankH){
            break;
        }
    }

    if (all){
        for(i = 0; i< bufferCount; i++){
            this->sneakyBit[numOfSneaky] = buffer[i];
            numOfSneaky ++;
        }

        if(count == this->RankH){
            n++;
            for(; n<acBit; n++){
                this->sneakyBit[numOfSneaky] = this->permu[n].index;
                numOfSneaky ++;
            }
        }
    }

    this->numOfSneaky = numOfSneaky;
    /*if(count != this->RankH){
        printf("The rank of this parity check matrix is %d not %d\n", count, this->RankH);
    }*/
    return;
}

uint32_t parityHasSol(const uint8_t *syndrome, const uint32_t *unreliCheck, const int32_t *rowPivotPlace, uint32_t numOfUnreliCheck){
    static uint32_t i, m;
    for(i = 0; i<numOfUnreliCheck; i++){
        m = unreliCheck[i];
        if(rowPivotPlace[m] == -1 && syndrome[m] != 0){
            return 0;
        }
    }
    return 1;
}

int32_t genBit2RowPos(RV_t **parity, uint32_t **bit2RowPos, int32_t *bitColWeight, uint32_t *sneakyBit, uint32_t numOfSneaky, uint32_t *unreliCheck, uint32_t numOfUnreliCheck, int32_t dist){
    static uint32_t n, nn, pos, shift, m, mm;
    static int32_t osd0orNot;
//printParity(parity, acRow, acCol);
    osd0orNot = 1;
    memset(bitColWeight, 0, NNbit*sizeof(*bitColWeight));
    for(n = 0; n<numOfSneaky; n++){
        nn = sneakyBit[n];
        pos = getBPos(nn);
        shift  = getBShift(nn);
        for(m = 0; m<numOfUnreliCheck; m++){
            mm = unreliCheck[m];
            if((parity[mm][pos]>>shift)&1){
                bit2RowPos[nn][bitColWeight[nn]++] = mm;
            }
        }
        if(bitColWeight[nn] >= (dist-1)){
            osd0orNot = 0;
        }
    }
    return osd0orNot;
}

uint32_t OSD_prepare2Flip(OSD *this){
    static uint32_t p_col, m, n, pos, shift, qShift, i, brother;
    static int32_t pivot, brotherM;
    static RV_t mask, result;
    this->bestErrorWeight = NN;

    int weight = 0;
    for(p_col = 0; p_col<P_col; p_col++){
        this->error[p_col] &= this->pivotMask[p_col];
        for(mask = ((RV_t)3)<<(QMask<<1); mask>0; mask>>=2){
            if(mask & this->error[p_col])
                weight++;
        }
    }
    // regen the syndrome
    for(i = 0; i<acRow; i++){
        m = this->unreliCheck[i];
        if(this->rowPivotPlace[m] >= 0){
            result = 0;
            for(p_col = 0; p_col<P_col; p_col++)
                result ^= (this->error[p_col] & this->unreliMask[p_col]) & this->parity[m][p_col];
            this->syndrome[m] ^= numOf1(result);
        }
    }

    // establish rowPivotStatus
    memset(this->rowPivotStatus, 255, MM*sizeof(*this->rowPivotStatus));
    memset(this->matchBrother, 255, NNbit*sizeof(*this->matchBrother));
    memset(this->pivotBrother, 255, NNbit*sizeof(*this->pivotBrother));
    for(m = 0; m<MM; m++){
        if((pivot = this->rowPivotPlace[m])>=0){
            pos = getBPos(pivot);
            shift  = getBShift(pivot);
            qShift = 126&shift;
            if((this->pivotMask[pos] >> qShift)&3){         // his brother is not a pivot, consider 1 or 0
                brother = pivot^1;
                this->pivotBrother[brother] = m;
                if((this->error[pos]>>qShift)&3){           // but his brother is already 1, set to 4, we don't have to change it
                    this->rowPivotStatus[m] = 1;
                }
                else{                                       // his brother is 0, weight is depend on me
                    this->rowPivotStatus[m] = 0;
                    if(this->syndrome[m]) weight++;
                }
            }
            else{                                           // his brother is also a pivot, consider 0, 1, 2, 3, 4
                if((brotherM = this->matchBrother[pivot]) >= 0){   // he find his brother
                    this->rowPivotBrother[m] = brotherM;
                    this->rowPivotBrother[brotherM] = m;
                    this->rowPivotStatus[m] = 2;
                    this->rowPivotStatus[brotherM] = 2;
                    if((this->syndrome[m] !=0) || (this->syndrome[brotherM]!=0)){ weight++;}
                }
                else{                                       // he is the first brother
                    brother = pivot^1;
                    this->matchBrother[brother] = m;
                }
            }
        }
    }

    #if DBG_PRT_OSD
        printf("New Syndrome - "); printSyndrome(this->syndrome, MM);
        printf("Error that mask out - "); printcError(this->error, P_col);
        printf("The weight before generate valid error: %d\n", weight);
        printf("Pivot status\n");
        printf("Pivot Brother:\n");
        for(int i = 0; i<NNbit; i++)
            printf("%d %d  ", i, this->pivotBrother[i]);
        printf("\n");
        printf("Row pivot Place:\n");
        for(int i = 0; i<MM; i++)
            printf("%d %d  ", i, this->rowPivotPlace[i]);
        printf("\n");
        printf("Row pivot Status:\n");
        for(int i = 0; i<MM; i++)
            printf("%d %d  ", i, this->rowPivotStatus[i]);
        printf("\n");
    #endif // DBG_PRT_OSD

    return weight;
}

void OSD_w_recursive(OSD *this, int32_t iter, int32_t start, int32_t weight){
    /*  error, syndrome are already prepared
        weight: the error weight when we get into here
        iter: how many bit have we flipped
        start: which bit should we start to flip */

    static RV_t mask;
    int32_t nowWeight, brotherM;
    uint32_t pivot, pos, shift;
//printf("%d ", weight);
    if(weight < this->bestErrorWeight){
        this->bestErrorWeight = weight;

        updateMostPossibleError(this->mostPossibleError, this->error, this->syndrome, this->unreliCheck, this->rowPivotPlace);
        #if DBG_PRT_OSD
        printf("\nBest weight: %d\n", this->bestErrorWeight);
        printcError(this->mostPossibleError, P_col);
        #endif // DBG_PRT_OSD
/*restoreError(this->mostPossibleError, this->returnError, NN);
printcError(this->mostPossibleError, P_col);
int count1= 0;
for(int i = 0; i<NN; i++){
    if(this->returnError[i]) count1++;
}
if(count1!=weight)
    system("pause");*/
    }

    if(iter == this->osdw){
        return;
    }
    // doing recursive part
    RV_t nowMask, brotherMask;
    uint32_t brother, change;
    int nn;
    for(int i = start; i<this->numOfSneaky; i++){
        nn = this->sneakyBit[i];
        pos = getBPos(nn);
        shift = getBShift(nn);
        if((this->pivotMask[pos]>>shift)&1){
            nowMask = 1;
            nowMask <<= shift;
            this->error[pos] ^= nowMask;
            nowWeight = weight;
            // deal with added weight + change the status of pivot
            brotherMask = ((shift&1)? (nowMask>>1):(nowMask<<1));
            if((brotherM = this->pivotBrother[nn]) >=0){     // his brother is a pivot
                change = 1;
                if(this->error[pos]&nowMask){                       // I'm now an error
                    this->rowPivotStatus[brotherM] = 1;
                    if(this->syndrome[brotherM] == 0) nowWeight = weight+1;
                }
                else{                                                               // I'm now not an error
                    this->rowPivotStatus[brotherM] = 0;
                    if(this->syndrome[brotherM] == 0) nowWeight = weight -1;
                }
            }
            else{
                if((this->error[pos] & brotherMask) == 0)
                    (this->error[pos]&nowMask)?(nowWeight = weight+1):(nowWeight = weight -1); // I'm now an error or not
            }
             // deal with syndrome, syndrome XOR1Col
            nowWeight = syndromeXOR1Col_wPos(this->bit2RowPos[nn], this->syndrome, this->bitColWeight[nn], this->rowPivotStatus, this->rowPivotBrother, nowWeight, 1);
            //syndromeXOR1Col(this->parity, pos, shift, this->syndrome, this->unreliCheck);
            // osd resursive
            OSD_w_recursive(this, iter+1, i+1, nowWeight);   // 1
            // all things recover
            this->error[pos] ^= nowMask;
            if(change == 1){
               this->rowPivotStatus[brotherM] = 1-this->rowPivotStatus[brotherM];
            }
            syndromeXOR1Col_wPos(this->bit2RowPos[nn], this->syndrome, this->bitColWeight[nn], this->rowPivotStatus, this->rowPivotBrother, weight, 0);
            //syndromeXOR1Col(this->parity, pos, shift, this->syndrome, this->unreliCheck);
        }
    }
}

void OSD_lambda_recursive(OSD *this, int32_t start, int32_t weight){
    /*  error, syndrome are already prepared
        weight: the error weight when we get into here
        iter: how many bit have we flipped
        start: which bit should we start to flip */

    static RV_t mask;
    int32_t nowWeight, brotherM;
    uint32_t pivot, pos, shift;

    if(weight < this->bestErrorWeight){
        this->bestErrorWeight = weight;
        //printf("Best weight: %d\n", this->bestErrorWeight);
        updateMostPossibleError(this->mostPossibleError, this->error, this->syndrome, this->unreliCheck, this->rowPivotPlace);
    }

    if(start == this->osdw){
        return;
    }

    // doing recursive part
    RV_t nowMask, brotherMask;
    uint32_t brother, change;
    int nn;
    for(int i = start; i<this->osdw; i++){
        nn = this->sneakyBit[i];
        pos = getBPos(nn);
        shift = getBShift(nn);

        if((this->pivotMask[pos]>>shift)&1){
            nowMask = 1;
            nowMask <<= shift;
            this->error[pos] ^= nowMask;

            // deal with added weight + change the status of pivot
            brotherMask = ((shift&1)? (nowMask>>1):(nowMask<<1));
            if((brotherM = this->pivotBrother[nn]) >=0){     // his brother is a pivot
                change = 1;
                if(this->error[pos]&nowMask){                       // I'm now an error
                    this->rowPivotStatus[brotherM] = 1;
                    if(this->syndrome[brotherM] == 0) nowWeight = weight+1;
                }
                else{                                                               // I'm now not an error
                    this->rowPivotStatus[brotherM] = 0;
                    if(this->syndrome[brotherM] == 0) nowWeight = weight -1;
                }
            }
            else{
                if(this->error[pos] & brotherMask){
                    nowWeight  = weight;
                }
                else{
                    (this->error[pos]&nowMask)?(nowWeight = weight+1):(nowWeight = weight -1); // I'm now an error or not
                }
            }

             // deal with syndrome, syndrome XOR1Col
            nowWeight = syndromeXOR1Col_wPos(this->bit2RowPos[nn], this->syndrome, this->bitColWeight[nn], this->rowPivotStatus, this->rowPivotBrother, nowWeight, 1);
            //syndromeXOR1Col(this->parity, pos, shift, this->syndrome, this->unreliCheck);

            // osd resursive
            OSD_lambda_recursive(this, i+1, nowWeight);   // 1
            // all things recover
            this->error[pos] ^= nowMask;
            if(change == 1){
               this->rowPivotStatus[brotherM] = 1-this->rowPivotStatus[brotherM];
            }
            syndromeXOR1Col_wPos(this->bit2RowPos[nn], this->syndrome, this->bitColWeight[nn], this->rowPivotStatus, this->rowPivotBrother, weight, 0);
            //syndromeXOR1Col(this->parity, pos, shift, this->syndrome, this->unreliCheck);
        }
    }
}

void syndromeXOR1Col(const RV_t **parity, uint32_t pos,uint32_t shift, uint8_t *syndrome, const uint32_t *unreliCheck){
    static uint32_t m, i;
    for(i = 0; i<acRow; i++){
        m = unreliCheck[i];
        syndrome[m] ^= (parity[m][pos]>>shift)&1;
    }
}

int32_t syndromeXOR1Col_wPos(const int *bit2RowPos, uint8_t *syndrome, uint32_t bitColWeight, int8_t *rowPivotStatus, uint32_t *rowPivotBrother, int32_t weight, int32_t updWeight){
    static uint32_t m, i;
    for(i = 0; i<bitColWeight; i++){
        m = bit2RowPos[i];
        syndrome[m] ^= 1;
        if(updWeight){
            if(rowPivotStatus[m] == 2){  // has brother
                if(!syndrome[rowPivotBrother[m]]){  // brother is 0
                    (syndrome[m])?(weight++):(weight--);
                }
            }
            else if(rowPivotStatus[m] == 0){
                (syndrome[m])?(weight++):(weight--);
            }
        }
    }
    return weight;
}


void updateMostPossibleError(RV_t *mostPossibleError, RV_t *error, uint8_t *syndrome, uint32_t *unreliCheck, int32_t *rowPivotPlace){
    static uint32_t m, pivot, pos, shift, qShift, i;

    // choose the error with weight near target weight, if same distance choose minimum one
    //int distNow = abs(nowWeight-this->targetWeight), distBest = abs(this->bestErrorWeight-this->targetWeight);
    //if( distNow < distBest  || (distNow == distBest && nowWeight<this->bestErrorWeight)){

   memcpy(mostPossibleError, error, P_col*sizeof(*error));
    for(i = 0; i<acRow; i++){
        m = unreliCheck[i];
        if(syndrome[m]){
            pivot = rowPivotPlace[m];
            pos = getBPos(pivot);
            shift = getBShift(pivot);
            mostPossibleError[pos] ^= ((RV_t)1)<<shift;
        }
    }
    #if DBG_PRT_OSD
        //printSyndrome(this->syndrome, this->MM);
        printcError(mostPossibleError, P_col);
    #endif // DBG_PRT_OSD

}

void load_A_OSD(OSD *this, FILE *fp){
    uint32_t TempN, TempM, max_num_n, dumb;
    uint32_t *num_data = (uint32_t*)calloc((size_t)this->NN, sizeof(*num_data));
	int32_t err;
	uint32_t i, total, pos, shift, val;
	uint32_t n, m, p_col;
	RV_t mask;
	// read basic parameter ============================================
	err = fscanf(fp,"%d%d", &TempN, &TempM);
	if(TempN != this->NN || TempM != this->MM){printf("Load A basic parameter NG!!\n"); exit(1); }
	err = fscanf(fp,"%d%d", &max_num_n, &dumb);
	if(err==EOF) { printf("scan A basic parm NG !!\n"); getchar(); exit(1); }
	if(dbg_prt) printf("[NN = %d, MM = %d, max_num_n = %d, max_num_m = %d]\n", TempN, TempM, max_num_n, dumb);

	// read how many index are there in the row and column =========================
    if(dbg_prt) printf("\n num_n = ");
    for(n = 0; n < this->NN; n++){
        err = fscanf(fp, "%d", &(num_data[n]));
        if(dbg_prt) printf("%d ", num_data[n]);
    }
    if(err==EOF) {printf("scan A num_n NG !!\n"); getchar(); exit(1); }

    // throw out the row part
    for(m = 0; m < this->MM; m++)
        err = fscanf(fp, "%d", &dumb);
    if(dbg_prt) printf("\n");

	// read where the index is and the weight ==================================
    for(n = 0; n<this->NN; n++){
        total = num_data[n];
        pos = getQPos(n);
        shift = getQShift(n);
        for(i = 0; i < total; i++){
            err = fscanf(fp, "%d%d", &m, &val);
            if(val == 1)
                mask = 2;
            else if(val == 2)
                mask = 1;

            if(dbg_prt) printf("n = %d, i = %d, m = %d, pos = %d, shift = %d, val = %d\n", n, i, m, pos, shift, val);
            mask <<= shift;

            this->parity_intact[--m][pos] |= mask;
        }
        for(; i<max_num_n; i++)
            err = fscanf(fp, "%d%d", &dumb, &dumb);
    }
    free(num_data);

    #if DBG_PRT_OSD
        //printParity(this->parity_intact, this->MM, this->P_col);
    #endif // DBG_PRT_OSD
}

void initParity(RV_t **parity,  RV_t **parity_intact){
    uint32_t copysize = P_col * sizeof(*parity[0]), m;
    for(m = 0; m<MM; m++){
        memcpy(parity[m], parity_intact[m], copysize);
    }
}

void printParity(const RV_t **parity, uint32_t row, uint32_t col){
    static uint32_t dumb, m, p_col;
    static int8_t n;
    static RV_t mask;
    printf("Parity check matrix\n");

    for(m = 0; m<row; m++){
        for(p_col = 0; p_col<col; p_col++){
            mask = ((RV_t)3)<<(QMask<<1);
            for (n = QMask; n>=0; n--){
                dumb = (mask & parity[m][p_col]) >> (n<<1);
                printf("%d ", dumb);
                mask >>= 2;
            }
        }
        printf("\n");
    }
}

void printcError(const RV_t *error, uint32_t len){
    static uint32_t dumb, p_col;
    static int8_t n;
    static RV_t mask;
    printf("Error (compressed)\n");

    for(p_col = 0; p_col<len; p_col++){
        mask = ((RV_t)3)<<(QMask<<1);
        for (n = QMask; n>=0; n--){
            dumb = (mask & error[p_col]) >> (n<<1);
            printf("%d ", dumb);
            mask >>= 2;
        }
    }
    printf("\n");
}

void printError(const uint8_t *error, uint32_t len){
    static uint32_t n;

    printf("Error\n");
    for(n= 0; n<len; n++){
        printf("%d ", (uint32_t)error[n]);
    }
    printf("\n");
}

void getSyndrome(const RV_t *error, const RV_t **parity, uint8_t *syndrome, uint32_t row, uint32_t col){
    static RV_t result;
    static uint32_t m, p_col;

    for(m = 0; m < row; m++){
        result = 0;
        for(p_col = 0; p_col < col; p_col++){
            result ^= ((parity[m][p_col]) & (error[p_col]));
        }
        syndrome[m] = numOf1(result);
    }
}

void printSyndrome(const uint8_t *syndrome, uint32_t len){
    static uint32_t m;
    printf("Syndrome:\n");
    for(m = 0; m<len; m++){
        printf("%d ", (uint32_t) syndrome[m]);
    }
    printf("\n");
}


/* (By convention (systematic form), the righthand M*M matrix is an invertible matrix FOR LINEAR CODE)
a_matrix_GFQ *A to
[1 2 2 1 0;
 0 1 2 2 1;
 1 0 1 2 2;
 2 1 0 1 2]
will have N=5; M=4,
	*num_n = [3 3 3 4 3], max_num_n = 4,
	*num_m = [4 4 4 4], max_num_m = 4,
	*nlist = [c1 c2 c3 c4 c5]     with *c1=[1 3 4 0], *c2=[1 2 4 0], *c3=[1 2 3 0], *c4=[1 2 3 4], *c5=[2 3 4 0], all non-zero -1 as index;
	*nlist_val = [v1 v2 v3 v4 v5] with *v1=[1 1 2 0], *v2=[2 1 1 0], *v3=[2 2 1 0], *v4=[1 2 2 1], *v5=[1 2 2 0];
	*mlist = [r1 r2 r3 r4]        with *r1=[1 2 3 4], *r2=[2 3 4 5], *r3=[1 3 4 5], *r4=[1 2 4 5], all non-zero -1 as index;;
	*mlist_val = [u1 u2 u3 u4 u5] with *u1=[1 2 2 1], *u2=[1 2 2 1], *u3=[1 1 2 2], *u4=[2 1 1 2];
	*ni2j = [j1 j2 j3 j4 j5] with *j1=[0 0 0], *j2=[1 0 1], *j3=[2 1 1], *j4=[3 2 2 2], *j5=[3 3 3];
			j1 means in col_1 the three non-zeros are all the first non-zeros in row order
*/


uint8_t numOf1(RV_t num){
    #if RV_LEN == 64
        num ^= num>>32;
    #endif // RV_LEN

    num ^= num>>16;
    num ^= num>>8;
    return numOf1Arr[(num&255)];
}

void compressError(const uint8_t *error, RV_t *cError, int NN, int P_col){
    static uint32_t n, pos, shift;
    static RV_t mask;
    memset(cError, 0, P_col*sizeof(*cError));

    for(n = 0; n<NN; n++){
        if(mask = error[n]){
            pos = getQPos(n);
            shift = getQShift(n);
            cError[pos] |= mask<<shift;
        }
    }
}

void restoreError(const RV_t *cError, uint8_t *error, int NN){
    static uint32_t n;
    static uint32_t pos, shift;
    for(n = 0; n<NN; n++){
        pos = getQPos(n);
        shift = getQShift(n);
        error[n] = (cError[pos]>>shift)&3;
    }
}
