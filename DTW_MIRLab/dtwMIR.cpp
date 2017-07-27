/*
 * Syntax: dtwMIR(queryMat, frameMat)
 * each column of the input matrices is a vector representing the preprocessed signal at a given time
 * return a matrix of accumulated cost and the optimal accumulated cost
 * */

#include <math.h>
#include "mex.h"

/*
 * Calculate cost (Hamming distance) for two vectors
 * ptrQuery is a pointer to the query matrix
 * queryCol is the column index of query
 * ptrRef is a pointer to the reference matrix
 * refCol is the column index of the reference
 * M is the dimension of each vector
 * */
int getCost(const double *ptrQuery, const int queryCol, const double *ptrRef, const int refCol, const int M) {
    int accum = 0;
    for (int m = 0; m < M; m++) {
        if (ptrQuery[queryCol * M + m] != ptrRef[refCol * M + m]) {
            ++accum;
        }
    }
    return accum;
}

/*
 * Generate a cost matrix
 * Cost[m][n] is the cost between the mth query vector and the nth reference vector
 * */
void getCostMat(double *ptrCost, const double *ptrQuery, const double *ptrRef, const int queryN, const int refN,
    const int M) {
    for (int n = 0; n < refN; n++) {
        for (int m = 0; m < queryN; m++) {
            ptrCost[n * queryN + m] = getCost(ptrQuery, m, ptrRef, n, M);
        }
    }
    return;
}

/*
 * Cauculate the matrix of accumulated cost
 * The base case assumes that no penalty would be given to start at any time
 * Each cell potentially has three options to fill in,
 * representing tempo of {(1, 1), (1, 2), (2, 1)}
 * */
void getSMat(double *ptrS, double* ptrBtMat, const double *ptrCost, const int numRow, const int numCol) {
   // base case
    for (int col = 0; col < numCol; col++) {
        ptrS[col * numRow] = ptrCost[col * numRow]; // no penalty for starting at any place
        ptrBtMat[col * numRow] = col * numRow; // set corresponding offset values
    }
   // DP
    for (int col = 0; col < numCol; col++) {
        for (int row = 1; row < numRow; row++) {
            double bestSoFar = INFINITY;
            int bestOffset = -1;
           // option 1: query and ref has the same tempo
            if (row - 1 >= 0 && col - 1 >= 0) {
                double opt1 = ptrS[(col - 1) * numRow + (row - 1)];
                if (opt1 != INFINITY && opt1 + ptrCost[col * numRow + row] < bestSoFar) {
                    bestSoFar = opt1 + ptrCost[col * numRow + row];
                    bestOffset = (col - 1) * numRow + (row - 1);
                }
            }
           // option 2: query is at twice the tempo as ref
            // penalize the cost weight in this case
            if (row - 2 >= 0 && col - 1 >= 0) {
                double opt2 = ptrS[(col - 1) * numRow + (row - 2)];
                if (opt2 != INFINITY && opt2 + 2 * ptrCost[col * numRow + row] < bestSoFar) {
                    bestSoFar = opt2 + 2 * ptrCost[col * numRow + row];
                    bestOffset = (col - 1) * numRow + (row - 2);
                }
            }
           // option 3: query is at half the tempo as ref
            if (row - 1 >= 0 && col - 2 >= 0) {
                double opt3 = ptrS[(col - 2) * numRow + (row - 1)];
                if (opt3 != INFINITY && opt3 + ptrCost[col * numRow + row] < bestSoFar) {
                    bestSoFar = opt3 + ptrCost[col * numRow + row];
                    bestOffset = (col - 2) * numRow + (row - 1);
                }
            }
            // consider the best option
            ptrS[col * numRow + row] = bestSoFar;
            ptrBtMat[col * numRow + row] = bestOffset;
        }
    }
    return;
}

/*
 * Return the optimal accumulated cost,
 * the smallest value at the last row across all columns
 * */
int getMinSCostCol(const double *ptrS, const int numRow, const int numCol) {
    double optimal = INFINITY;
    int index = -1;
    for (int col = 0; col < numCol; col++) {
        if (ptrS[col * numRow + numRow - 1] < optimal) {
            optimal = ptrS[col * numRow + numRow - 1];
            index = col;
        }
    }
    return index;
}

/*
 * Return the optimal offset,
 * the column index of the first row where the optimal path begins
 * */
int getOffset(const double *ptrBtMat, const int optCol, const int numRow, const int numCol) {
    int row = numRow - 1;
    int col = optCol;
    while (row > 0) {
        int index = row + col * numRow;
        index = ptrBtMat[index];
        row = index % numRow;
        col = index / numRow;
    }
    return col;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
#define S plhs[0]
#define optSCost   plhs[1]
#define btMat plhs[2]
#define optOffset plhs[3]
#define queryMat    prhs[0]
#define refMat    prhs[1]
   // read input
    int queryM = mxGetM(queryMat);
    int queryN = mxGetN(queryMat);
    int refM = mxGetM(refMat);
    int refN = mxGetN(refMat);
    int M;
    if (queryM == refM) {
        M = queryM;
    } else {
        M = 1;
    }
   // set matrices
    mxArray *costMat = mxCreateDoubleMatrix(queryN, refN, mxREAL);
    S = mxCreateDoubleMatrix(queryN, refN, mxREAL);
    optSCost = mxCreateDoubleMatrix(1, 1, mxREAL);
    btMat = mxCreateDoubleMatrix(queryN, refN, mxREAL);
    optOffset = mxCreateDoubleMatrix(1, 1, mxREAL);
   // set pointers
    double *ptrQuery = mxGetPr(queryMat);
    double *ptrRef = mxGetPr(refMat);
    double *ptrCost = mxGetPr(costMat);
    double *ptrS = mxGetPr(S);
    double *ptrOptSCost = mxGetPr(optSCost);
    double *ptrBtMat = mxGetPr(btMat);
    double *ptrOptOffset = mxGetPr(optOffset);
   // generate cost matrix to costMat
    getCostMat(ptrCost, ptrQuery, ptrRef, queryN, refN, M);
    getSMat(ptrS, ptrBtMat, ptrCost, queryN, refN);

    int minSCostCol = getMinSCostCol(ptrS, queryN, refN);
    // find optimal accumulated cost and convert it to score
    ptrOptSCost[0] = 1 - (ptrS[minSCostCol * queryN + queryN - 1] / (queryM * queryN));
    // find optimal offset
    ptrOptOffset[0] = getOffset(ptrBtMat, minSCostCol, queryN, refN);
    
    return;
}
