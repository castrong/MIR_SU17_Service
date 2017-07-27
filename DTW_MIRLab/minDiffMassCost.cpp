/*
 * Syntax: minDiffMassCost(queryMat, refMat)
 * queryMat is a matrix specifying the CQT coefficients of the query file, one column per frame
 * refMat is the same matrix for the reference file
 * returns a matrix of the pairwise frame costs, where we use the minimum differential mass cost function
 * */

#include <math.h>
#include "mex.h"

/*
 * Calculate pairwise cost between query and reference files
 * The cost used here is a minimum differential mass
 * */
void getCostMat(double *ptrCost, const double *ptrQuery, const double *ptrRef, const int queryFrms, const int refFrms,
    const int nbins) {
    for (int r = 0; r < refFrms; r++) {
        for (int q = 0; q < queryFrms; q++) {
	  double accum = 0;
	  for (int bin = 0; bin < nbins; bin++) {
	    double diff = ptrQuery[q * nbins + bin] - ptrRef[r * nbins + bin];
	    if (diff > 0) {
	      accum += diff;
	    }
	  }
	  ptrCost[r * queryFrms + q] = accum;
        }
    }
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
#define costMat plhs[0]
#define queryMat    prhs[0]
#define refMat    prhs[1]

   // read input
    int nbins = mxGetM(queryMat);
    int nQueryFrms = mxGetN(queryMat);
    int nRefFrms = mxGetN(refMat);

   // set pointers
    costMat = mxCreateDoubleMatrix(nQueryFrms, nRefFrms, mxREAL);
    double *ptrQuery = mxGetPr(queryMat);
    double *ptrRef = mxGetPr(refMat);
    double *ptrCost = mxGetPr(costMat);

   // generate cost matrix to costMat
    getCostMat(ptrCost, ptrQuery, ptrRef, nQueryFrms, nRefFrms, nbins);

    return;
}
