#include "mex.h"
#include <stdio.h>

struct DataSizes {
    size_t m_rows;                              // number of rows
    size_t m_cols;                              // number of columns
    size_t n_sens;                              // number of sensor elements
    size_t m_numel;                             // number of elements in index matrix
    size_t n_pixels;                            // image elements
    size_t frame_size;                          // rfdata number of rows
    size_t n_frames;                            // rfdata number of frames
};

// DAS Function
void dasindex_single_frame(unsigned int *M, double *rf, double* img_out, struct DataSizes sz)
{

    //memory indexing
    for (int sen = 0; sen < sz.n_sens; ++sen) {
        for(int pixel = 0; pixel < sz.n_pixels; ++pixel) {
            *(img_out + pixel) += *(rf + *(M + pixel + sz.n_pixels * sen));
        }
    }

}

void dasindex_multi_frame(unsigned int* M, double *rf, double* img_out, struct DataSizes sz) {

    int frame_img, frame_rf;

    //outer loop for each frame
    for (int frame = 0; frame < sz.n_frames; ++frame) {
        frame_img = frame*sz.n_pixels;
        frame_rf = frame*sz.n_sens*sz.frame_size;

        // inner loop for memory indexing as before
        for (int sen = 0; sen < sz.n_sens; ++sen) {
            for(int pixel = 0; pixel < sz.n_pixels; ++pixel) {
                *(img_out + pixel + frame_img) += *(rf + *(M + pixel + sz.n_pixels * sen) + frame_rf);
            }
        }
    }

}

// Gatewy /MATLAB Wrapper function. 0.057809 seconds average for a 400x400x500 delay matrix
void mexFunction( int nlhs, mxArray *outputs[], int nrhs, const mxArray *inputs[]) {
    
    // Groom Inputs ------------------------------------------------------------
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs index matrix and rfdata required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Too many output arguments.");
    }
    // make sure the first input argument is an array of unsigned int 32.
    if( mxIsComplex(inputs[0]) || !mxIsUint32(inputs[0]) ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Unsupported index matrix datatype. Function currently supports: UInt32");
    }
    // ensure the second input argument is type double
    if( !mxIsDouble(inputs[1]) || mxIsComplex(inputs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input rfdata must be of datatype double or single (not complex)");
    }    

    // preallocate and setup variables ------------------------------------------------
        
    double * img;                               // beamformed image pointer
    unsigned int * M;                           // Index Matrix pointer
    double * rfdata;                            // rfdata matrix pointer
    const mwSize * dims;                        // M dimensions
    const mwSize * rfdims;

    M = (unsigned int*) mxGetPr(inputs[0]);       //import Index Matrix
    rfdata = (double *) mxGetPr(inputs[1]);       //import rfdata
    dims = mxGetDimensions(inputs[0]);            //get the dimensions of the index matrix
    rfdims = mxGetDimensions(inputs[1]);          //get the dimensions of rfdata

    //some error checking
    if (rfdims[1] != dims[2]) {
        mexErrMsgIdAndTxt("MyToolbox:DAS_Index:InvalidInput","Mismatched dimensions on either rfdata or the index matrix");        
    }

    struct DataSizes var_sz = {
        dims[0],
        dims[1],
        dims[2],
        dims[0] * dims[1] * dims[2],
        dims[0] * dims[1],
        rfdims[0],
        1
    };
    
    if (mxGetNumberOfDimensions(inputs[1]) == 3) {
        
        // create the output image matrix
        var_sz.n_frames = rfdims[2];
        const mwSize output_dimensions[3] ={var_sz.m_rows, var_sz.m_cols, var_sz.n_frames};
        outputs[0] = mxCreateNumericArray(3, &output_dimensions, mxDOUBLE_CLASS, mxREAL);

        img = (double*) mxGetPr(outputs[0]);

        //run the multiframe DAS function
        dasindex_multi_frame(M, rfdata, img, var_sz);
    
    
    } else {

        // create the output image matrix
        outputs[0] = mxCreateDoubleMatrix(var_sz.m_rows, var_sz.m_cols, mxREAL);
        img = (double*) mxGetPr(outputs[0]);

        // Run the DAS function
        dasindex_single_frame(M, rfdata, img, var_sz);
    }

}