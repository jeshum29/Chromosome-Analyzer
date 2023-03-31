/*==========================================================
 * skel3.c
 *
 * This function operates on a 3d binary image. 
 *
 * It returns the skeleton of the "shapes" in the input image.
 *
 * The calling syntax is 
 * 
 * skeleton = skel3(binaryImage)
 *
 * The code is adopted from a Java version written for ImageJ, the citation for which is
 *
      *  This work is an implementation by Ignacio Arganda-Carreras of the
      * 3D thinning algorithm from Lee et al. "Building skeleton models via 3-D
      * medial surface/axis thinning algorithms. Computer Vision, Graphics, and
      * Image Processing, 56(6):462â??478, 1994." 
      * Based on the ITK version from Hanno Homann 
      * <a href="http://hdl.handle.net/1926/1292"> http://hdl.handle.net/1926/1292</a>
      *
      *  More information at Skeletonize3D homepage:
      *  http://imagejdocu.tudor.lu/doku.php?id=plugin:morphology:skeletonize3d:start
      *
      * @version 1.0 11/19/2008
      * @author Ignacio Arganda-Carreras <ignacio.arganda@uam.es>
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
#include "skel.h"



 
/* The computational routine */
void compRoutine( int *outputImage, int *outputImage2, mwSize dims[3]){
    
    int eulerLUT[256];
    
    int i,j,k,n;
    int iter = 1;
    int unchangedBorders = 0;
    int currentBorder;
    int isBorderPoint;
    int numberOfNeighbors;
    int neighbor[27];
    int index[3];
    int sz[3];
    int noChange;
    int counter;
    int val, bin;
    int *simpleBorderPointsX, *simpleBorderPointsY, *simpleBorderPointsZ;
    
    int x, y, z, rows, cols, depth;
    
    int test[3];
    
    rows = (int)dims[0];
    cols = (int)dims[1];
    depth = (int)dims[2];
    
    simpleBorderPointsX = (int*)malloc(rows*cols*depth*sizeof(int));
    simpleBorderPointsY = (int*)malloc(rows*cols*depth*sizeof(int));
    simpleBorderPointsZ = (int*)malloc(rows*cols*depth*sizeof(int));

    sz[0] = rows; 
    sz[1] = cols; 
    sz[2] = depth;

    fillEulerLUT( eulerLUT );

    
//     // Loop through the image several times until there is no change.
    while( unchangedBorders < 6 ){  // loop until no change for all the six border types
    
        unchangedBorders = 0;
        counter = 0;
        
        for(  currentBorder = 1; currentBorder <= 6; currentBorder++) {
            
            // Loop through the image.
            for ( z = 0; z < depth; z++){
                for ( y = 0; y < cols; y++){
                    for ( x = 0; x < rows; x++){
                        
                        
                        // check if point is foreground
                        val = getPixel(outputImage, x, y, z, sz);                    

                        if ( val==0 ){
                            continue;         
                        }
                        
                        // check 6-neighbors if point is a border point of type currentBorder
                        isBorderPoint = 0;
                        // North
                        if( currentBorder == 1 && getPixel(outputImage, x, y-1, z, sz) == 0 )
                            isBorderPoint = 1;
                        // South
                        if( currentBorder == 2 && getPixel(outputImage, x, y+1, z, sz) == 0 )
                            isBorderPoint = 1;
                        // East
                        if( currentBorder == 3 && getPixel(outputImage, x+1, y, z, sz) == 0 )
                            isBorderPoint = 1;
                        // West
                        if( currentBorder == 4 && getPixel(outputImage, x-1, y, z, sz) == 0 )
                            isBorderPoint = 1;
                        // Up
                        if( currentBorder == 5 && getPixel(outputImage, x, y, z+1, sz) == 0 )
                            isBorderPoint = 1;
                        // Bottom
                        if( currentBorder == 6 && getPixel(outputImage, x, y, z-1, sz) == 0 )
                            isBorderPoint = 1;
                        
                                                
                        if( isBorderPoint==0 )
                        {
                            continue;         // current point is not deletable
                        }

                        // check if point is the end of an arc
                        numberOfNeighbors = -1;   // -1 and not 0 because the center pixel will be counted as well

                        getNeighborhood(outputImage, neighbor, x, y, z, sz);
                        
                        
                        for(  i = 0; i < 27; i++ ) // i =  0..26
                        {
                            if( neighbor[i] == 1 )
                                numberOfNeighbors++;
                        }
                        
                        if( numberOfNeighbors == 1 )
                        {
                            continue;         // current point is not deletable
                        }
                        
                        // Check if point is Euler invariant
                        
                         bin = isEulerInvariant(neighbor, eulerLUT);
                     
                        if (bin!=0)
                        {
                            continue;         // current point is not deletable
                        }
                       
                        
                        // Check if point is simple (deletion does not change connectivity in the 3x3x3 neighborhood)
                        
                         bin = isSimplePoint(neighbor);
                        
                        if ( bin==0 )
                        {
                            continue;         // current point is not deletable
                        }


                        // add all simple border points to a list for sequential re-checking
                        
                        *(simpleBorderPointsX + counter) = x;
                        *(simpleBorderPointsY + counter) = y;
                        *(simpleBorderPointsZ + counter) = z;
                        
                        counter++;
                                                
                    }
                }
            }
                        
            
           // sequential re-checking to preserve connectivity when deleting in a parallel way
            
            noChange = 1;
            
            for ( i = 0;    i < counter ; i++)
            {
                
                x = *(simpleBorderPointsX + i);
                y = *(simpleBorderPointsY + i);
                z = *(simpleBorderPointsZ + i);
                                
                // 1. Set simple border point to 0
                
                setPixel( outputImage, x, y, z, sz, 0);

                // 2. Check if neighborhood is still connected
                getNeighborhood(outputImage, neighbor, x, y, z, sz);

                bin = isSimplePoint(neighbor);
                
                if ( bin==0 ){
                    
                    // we cannot delete current point, so reset
                    setPixel( outputImage, x, y, z, sz, 1);

                } else {
                    noChange = 0;
                }
                
            }
            
            if ( noChange==1 ) {
                unchangedBorders++;
            }
            
            counter = 0;
                            
        } // end currentBorder for loop
        
        // Progress bar iterations
        iter++;
    //    printf("iter = %d\n", iter);
        
        if (iter > 10){
        break;
    }
        
    } // end while
    
    free(simpleBorderPointsX);
    free(simpleBorderPointsY);
    free(simpleBorderPointsZ);

    
} /* end computeThinImage */


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double nrr, ncc, nzz;
    double *im_in_double, *im_out_double;
    int *im_in, *im_out, *im_in_copy;
    double val;
    int valI;
    mwSize dims[3];
    
    int *start_of_pr;
    int bytes_to_copy;
    int NN;
    int nn;
    
    if(nrhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","three output required.");
    }
    
    im_in_double = mxGetPr(prhs[0]);
    
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
    
    
    im_out_double = (double*)malloc(NN*sizeof(double));
    im_in = (int*)malloc(NN*sizeof(int));
    im_in_copy = (int*)malloc(NN*sizeof(int));
    im_out = (int*)malloc(NN*sizeof(int));
    
    
    for (nn = 0; nn < NN; nn++){
        *(im_in + nn) = (int)*(im_in_double + nn);
        *(im_in_copy + nn) = 0;
    };
    
  compRoutine(im_in, im_in_copy, dims);
  
    
    for (nn = 0; nn < NN; nn++){
        *(im_out_double + nn) = (double) *(im_in + nn);
    };
    
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    start_of_pr = (int *)mxGetData(plhs[0]);
    bytes_to_copy = NN * mxGetElementSize(plhs[0]);
    memcpy(start_of_pr, im_out_double, bytes_to_copy);
        
    free(im_out);
    free(im_in);
    free(im_out_double);
    free(im_in_copy);
    
}


