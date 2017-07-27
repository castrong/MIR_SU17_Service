/*
 * Syntax: segmentalDTW_backtrace(C,queryLengths)
 * C is a cost matrix containing segment-level path scores where rows correspond to different query files and 
 * columns correspond to different reference frames
 * queryLengths is a row vector containing the lengths of the query files in frames
 * returns a row vector containing the optimal endpoints of each query file
 * */

#include <math.h>
#include "mex.h"
#include <vector>
using namespace std;

/*
 * Computes cumulative cost matrix
 */
void computeS(const double *ptrCostMat, const int rows,const int cols, const int32_T *queryLens, double *ptrS, int32_T *ptrBtMat) {
  // initialize first row
  for (int col = 0; col < cols; col++) {

    double bestScore = ptrCostMat[ col * rows ];
    int btIndex = 1;

    if ( col - 1 >= 0 && ptrS[ (col - 1) * rows ] < bestScore) {
      bestScore = ptrS[ (col - 1) * rows];
      btIndex = 0;
    }

    ptrS[ col * rows ] = bestScore;
    ptrBtMat[ col * rows ] = btIndex;
  }
  // DP
  for (int row = 1; row < rows; row++) {
    int hopSize = queryLens[row] / 2; // min match length for this query file
    for (int col = 0; col < cols; col++) {

      double bestScore = INFINITY;
      int btIndex = -1;

      // option 1: (0,-1)
      if (col - 1 >= 0 && ptrS[ (col - 1) * rows + row ] < bestScore) {
	bestScore = ptrS[ (col - 1) * rows + row ];
	btIndex = 0;
      }

      // option 2: (-1,-hopSize)
      if (col - hopSize >= 0 && ptrS[ (col - hopSize) * rows + row - 1 ] + ptrCostMat[ col * rows + row ] < bestScore) {
	bestScore = ptrS[ (col - hopSize) * rows + row - 1 ] + ptrCostMat[ col * rows + row ];
	btIndex = 1;
      }
      
      ptrS[ col * rows + row ] = bestScore;
      ptrBtMat[ col * rows + row ] = btIndex;
    }
  }
}

/*
 * Find optimal path using backtracing
 */
void determinePath(const int32_T *ptrBtMat, const int rows, const int cols, const int32_T *queryLens, int32_T *ptrEndPts) {
  int row = rows - 1;
  int col = cols - 1;
  int hopSize;
  while (row >= 0) {
    int btIndex = ptrBtMat[ col * rows + row ];
    if (btIndex == 0) { // (0,-1)
      col = col - 1;
    }
    else if (btIndex == 1) { // (-1,-hopSize)
      ptrEndPts[row] = col + 1; // return result in Matlab-style indexing
      hopSize = queryLens[row] / 2;
      row = row - 1;
      col = col - hopSize;
    }
    else { // invalid
      mexPrintf("Invalid backtrace index!  Row %d, col %d, btIndex %d", row, col, btIndex);
      return;
    }
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
#define endPts plhs[0]
#define costMat prhs[0]
#define queryLengths prhs[1]
#define S plhs[1]
#define btMat plhs[2]

    // initialize
    int numQueries = mxGetM(costMat);
    int numRefFrames = mxGetN(costMat);
    double *ptrCostMat = mxGetPr(costMat);
    int32_T *ptrQueryLengths = (int32_T *) mxGetData(queryLengths);
    endPts = mxCreateNumericMatrix(1, numQueries, mxINT32_CLASS, mxREAL);
    int32_T *ptrEndPts = (int32_T *) mxGetData(endPts);
    S = mxCreateDoubleMatrix(numQueries, numRefFrames, mxREAL);
    double *ptrS = mxGetPr(S);
    btMat = mxCreateNumericMatrix(numQueries,numRefFrames, mxINT32_CLASS, mxREAL);
    int32_T *ptrBtMat = (int32_T *) mxGetData(btMat);
    
    // calculate cumulative cost matrix
    computeS(ptrCostMat,numQueries,numRefFrames,ptrQueryLengths,ptrS,ptrBtMat);

    // backtrace
    determinePath(ptrBtMat,numQueries,numRefFrames,ptrQueryLengths,ptrEndPts);

    return;
}
