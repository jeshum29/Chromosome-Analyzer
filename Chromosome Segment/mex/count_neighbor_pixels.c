/*==========================================================
 * find_isolated_pixels_3.c
 *
 * this function operates on a 3d labeled binary image
 * it returns all unconnected pixels
 *
 * Keith Cheveralls
 * Dernburg Lab
 * August 2011
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

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
void compRoutine(double *im_in, double *im_out, double *inds, mwSize dims[3]){
    
    double *nnbs;
    
    int nrr, ncc, nzz, rr, cc, zz, ind, ind_kk, ii, jj, kk, iii, jjj, kkk, nn, nnn, NN, NN_b, nn_b;
    int xx, control_param, counter, counter2;
    int tmp;
    int conn;
    int sz[3];
    double maxVal;
    
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
        
        *(im_out+nn) = 0;
        
    }
    
    for (nn = 0; nn < NN; nn++){
        
        tmp = (int)fmod(nn, nrr*ncc);
        ii = (int)fmod(tmp, nrr);
        jj = (tmp - ii)/nrr;
        kk = (nn - tmp)/(nrr*ncc);
        
        val = getPixel(im_in, ii, jj, kk, sz);
        
        // only looking for borders between regions
        
        if (val==0) {
            continue;
        }
        
        getNeighborhood(im_in, nnbs, ii, jj, kk, sz, conn);
        
        counter = 0;
        
        for (nnn = 0; nnn < conn; nnn++){
            
            if (*(nnbs+nnn) != 0) {
                
                counter++;

            }
        }
        
        
    setPixel(im_out, ii, jj, kk, sz, counter);
    
    }
    
    free(nnbs);
    
} //end function


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double nrr, ncc, nzz, val;
    double *im_in, *im_out, *inds;
    
    mwSize dims[3];
    
    int *start_of_pr;
    int bytes_to_copy;
    int NN;
    int nn;
    
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","three output required.");
    }
    
    inds = mxGetPr(prhs[1]);
    im_in = mxGetPr(prhs[0]);
    
    nrr = (mxGetDimensions(prhs[0]))[0];
    ncc = (mxGetDimensions(prhs[0]))[1];
    nzz = (mxGetDimensions(prhs[0]))[2];
    
// printf("nrr = %f\n",nrr);
// printf("ncc = %f\n",ncc);
// printf("nzz = %f\n",nzz);
    
    dims[0] = nrr;
    dims[1] = ncc;
    dims[2] = nzz;
    
    NN = (int)nrr*ncc*nzz;
    
    im_out = (double*)malloc(NN*sizeof(double));
    
//     plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
//     total_distance = mxGetPr(plhs[2]);
    
    compRoutine(im_in, im_out, inds, dims);
    
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    start_of_pr = (int *)mxGetData(plhs[0]);
    bytes_to_copy = NN * mxGetElementSize(plhs[0]);
    memcpy(start_of_pr, im_out, bytes_to_copy);
    
    free(im_out);
    
}


