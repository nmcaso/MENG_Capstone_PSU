#include "mex.h"
#include <stdio.h>

struct DataSizes {
    size_t m_rows;                              // number of rows
    size_t m_cols;                              // number of columns
    size_t m_depth;                             // number of sensor elements
    size_t m_numel;                             // number of elements in index matrix
    size_t image_elements;                      // image elements
    size_t rf_rows;                             // rfdata number of rows
    size_t rf_cols;                             // rfdata number of columns
};

// DAS Function
void dasindex(unsigned int *M, double *rf, double* img_out, struct DataSizes sz)
{
    
    

}

// Gatewy /MATLAB Wrapper function:
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Groom Inputs ------------------------------------------------------------
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs index matrix and rfdata required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Too many output arguments.");
    }
    // make sure the first input argument is an array of unsigned int 32.
    if( mxIsComplex(prhs[0]) || !mxIsUint32(prhs[0]) ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Unsupported index matrix datatype. Function currently supports: UInt32");
    }
    // ensure the second input argument is type double
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
        ("MyToolbox:arrayProduct:notDouble","Input rfdata must be of datatype double or single (not complex)");
    }    

    // preallocate and setup variables ------------------------------------------------
        
    double * img;                               // beamformed image pointer
    unsigned int * M;                           // Index Matrix pointer
    double * rfdata;                            // rfdata matrix pointer
    const mwSize * dims;                        // M dimensions

    M = (unsigned int*) mxGetPr(prhs[0]);       //import Index Matrix
    rfdata = (double *) mxGetPr(prhs[1]);       //import rfdata
    dims = mxGetDimensions(prhs[0]);            // get the dimensions of the index matrix
    rfdims = mxGetDimensions(prhs[1]);

    //some error checking
    if (rfdims[0] != dims[2]) {
        mexErrMsgIdAndTxt("MyToolbox:DAS_Index:InvalidInput","Mismatched dimensions on either rfdata or the index matrix");        
    }

    struct DataSizes var_sz = {
        dims[0],
        dims[1],
        dims[2],
        dims[0] * dims[1] * dims[2],
        dims[0] * dims[1],
        mxGetM(prhs[1]),
        mxGetN(prhs[1])
    };
    
    // create the output image matrix
    plhs[0] = mxCreateDoubleMatrix(var_sz.m_rows, var_sz.m_cols, mxREAL);
    img = (double*) mxGetPr(plhs[0]);

    // Run the DAS function
    dasindex(M, rfdata, img, var_sz);

}
