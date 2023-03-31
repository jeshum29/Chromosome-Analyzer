/*==========================================================
 * find_all_borders_3.c
 *
 * this function operates on two 3d binary images
 * it returns the labels of the regions in the label image that
 * neighbor the region(s) in the first image
 *
 * In practice, we expect the number of neighboring regions to be small (< 10)
 * this function returns only the first ten for simplicity of memory management
 *
 * Keith Cheveralls
 * Dernburg Lab
 * Sept 2011
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
void compRoutine(double *im_in, double *im_out,  double *inds, mwSize dims[3]){
    
    double *nnbs;
    
    int nrr, ncc, nzz, rr, cc, zz, ind, ind_kk, ii, jj, kk, iii, jjj, kkk, nn, nnn, NN, NN_b, nn_b;
    int xx, control_param, counter, counter2;
    int tmp;
    int conn;
    int sz[3];
    double maxVal;
    
    int overwrite_flag;
    double val;
    int numNeighbors;
    
    nrr = (int)dims[0];
    ncc = (int)dims[1];
    nzz = (int)dims[2];
    sz[0] = nrr; sz[1] = ncc; sz[2] = nzz;
    
    NN = nrr*ncc*nzz;
    
    conn = inds[0];
    numNeighbors = 0;
    
    nnbs = (double*)malloc(conn*sizeof(double));
    
    for (nn = 0; nn < NN; nn++){
        
        tmp = (int)fmod(nn, nrr*ncc);
        ii = (int)fmod(tmp, nrr);
        jj = (tmp - ii)/nrr;
        kk = (nn - tmp)/(nrr*ncc);
        
        val = getPixel(im_in, ii, jj, kk, sz);
        
        if (val==0) {
            continue;
        }
        
        getNeighborhood(im_in, nnbs, ii, jj, kk, sz, conn);
        
        // if the pixel is surrounded by background
        if (maxVal == 0) {
            continue;
        }
        
        numNeighbors = 0;
        
        for (nnn = 0; nnn < conn; nnn++){
            if (nnbs[nnn]==1)
                numNeighbors++;
        }
        
        if (numNeighbors==conn)
            setPixel(im_out, ii, jj, kk, sz, 1);
           
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
    int total_num_neighbors;
    
    total_num_neighbors = 10;
    
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=0) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","three output required.");
    }
    
    im_in = mxGetPr(prhs[0]);
    im_out = mxGetPr(prhs[1]);
    inds = mxGetPr(prhs[2]);
    
    nrr = (mxGetDimensions(prhs[0]))[0];
    ncc = (mxGetDimensions(prhs[0]))[1];
    nzz = (mxGetDimensions(prhs[0]))[2];
    
    dims[0] = nrr;
    dims[1] = ncc;
    dims[2] = nzz;
    
    compRoutine(im_in, im_out, inds, dims);
    
    
}


