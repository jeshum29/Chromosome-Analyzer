/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier)
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2008 The MathWorks, Inc.
 *
 *========================================================*/
/* $Revision: 1.1.10.1 $ */

#include "mex.h"
#include "matrix.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<conio.h>
#include "im3d.h"

// sort ascending
int sort(const void *x, const void *y) {
    return (*(double*)x - *(double*)y);
}

/* The computational routine */
void compRoutine(double *regions_in, double *surfs_in, double *im_out_1, double *im_out_2, double *inds, mwSize dims[3]){
    
    double *nnbs;
    
    int nrr, ncc, nzz, rr, cc, zz, ind, ind_kk, ii, jj, kk, iii, jjj, kkk, nn, nnn, NN, NN_b, nn_b;
    int xx, control_param, counter, counter2;
    int tmp;
    int conn;
    int sz[3];
    double maxVal, maxVal_2;
    
    int overwrite_flag;
    double val;
    
    
    nrr = (int)dims[0];
    ncc = (int)dims[1];
    nzz = (int)dims[2];
    sz[0] = nrr; sz[1] = ncc; sz[2] = nzz;
    
    NN = nrr*ncc*nzz;
    
    // neighborhood size
    rr = (int)1;
    cc = (int)1;
    zz = (int)1;
    
    conn = inds[0];
    
    nnbs = (double*)malloc(conn*sizeof(double));
    
    for (nn = 0; nn < NN; nn++){
        
        tmp = (int)fmod(nn, nrr*ncc);
        ii = (int)fmod(tmp, nrr);
        jj = (tmp - ii)/nrr;
        kk = (nn - tmp)/(nrr*ncc);
        
        val = getPixel(surfs_in, ii, jj, kk, sz);
        
        // only looking for borders between regions
        
        if (val==0) {
            continue;
        }
        
         getNeighborhood(regions_in, nnbs, ii, jj, kk, sz, conn);
//         
         qsort(nnbs, conn, sizeof(double), sort);
         maxVal = *(nnbs + conn - 1);
//         
        // if the pixel is surrounded by background
        if (maxVal == 0) {
            continue;
        }
        
         counter = 0;
         maxVal_2 = 0;
         
        // look for a second nonzero value in the nhood, which means the pixel is a border
        for (nnn = 0; nnn < conn; nnn++){
            
            if (*(nnbs+nnn) != maxVal && *(nnbs+nnn) != 0) {
                
                counter++;
                maxVal_2 = maxVal;
                maxVal = *(nnbs+nnn);

            }
        }
         
         if (counter==2) {
            setPixel(im_out_1, ii, jj, kk, sz, maxVal);
            setPixel(im_out_2, ii, jj, kk, sz, maxVal_2);
         }
    
    }
    
    free(nnbs);
    
} //end function


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double nrr, ncc, nzz, val;
    double *regions_in, *surfs_in, *im_out_1, *im_out_2, *inds;
    
    mwSize dims[3];
    
    int *start_of_pr;
    int bytes_to_copy;
    int NN;
    int nn;
    
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=0) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","three output required.");
    }
    
    
    regions_in = mxGetPr(prhs[0]);
    surfs_in = mxGetPr(prhs[1]);
    im_out_1 = mxGetPr(prhs[2]);
    im_out_2 = mxGetPr(prhs[3]);
    
    inds = mxGetPr(prhs[4]);
    
    nrr = (mxGetDimensions(prhs[0]))[0];
    ncc = (mxGetDimensions(prhs[0]))[1];
    nzz = (mxGetDimensions(prhs[0]))[2];
    
    dims[0] = nrr;
    dims[1] = ncc;
    dims[2] = nzz;
            
    compRoutine(regions_in, surfs_in, im_out_1, im_out_2, inds, dims);
    
    
}








