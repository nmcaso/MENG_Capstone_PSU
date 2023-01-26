#include "C:/Program Files/MATLAB/R2022b/extern/include/mex.h"
#include "C:/Program Files/MATLAB/R2022b/toolbox/parallel/gpu/extern/include/gpu/mxGPUArray.h"
#include <chrono>

#define timer       std::chrono::high_resolution_clock
#define timertime   std::chrono::high_resolution_clock::time_point
#define timecast    std::chrono::duration_cast<std::chrono::microseconds>
#define timesecs    std::chrono::microseconds

//a convenient structure
typedef struct sizes {
            size_t imax;
            size_t jmax;
            size_t kmax;
            size_t xysize;
            size_t M_numel;
} sizes;

//a device function to reduce the last several warps of the reduction faster.
__device__ void warpreduce(volatile double* s_matrix, int thread_vector) {
    s_matrix[thread_vector] += s_matrix[thread_vector + 32];
    s_matrix[thread_vector] += s_matrix[thread_vector + 16];
    s_matrix[thread_vector] += s_matrix[thread_vector + 8];
    s_matrix[thread_vector] += s_matrix[thread_vector + 4];    
    s_matrix[thread_vector] += s_matrix[thread_vector + 2];
    s_matrix[thread_vector] += s_matrix[thread_vector + 1];
}

//kernel Creates a beamformed image array on the GPU.
__global__ void DAS_Index_GPU(const unsigned int* indmat, const double* rfdata, double* img3d) {    

    //nice easy 1D kernel
    int z = threadIdx.x;
    int xy = blockIdx.x;
    int xyz = xy * blockDim.x + z; 

    img3d[xyz] = rfdata[indmat[xyz]];
}

//going to have to transpose this to be faster. I'll do that tomorrow. take a look at https://github.com/shwina/cuper/blob/master/cuTranspose/transpose3d.cu, dev_transpose_102_in_place

//kernel to take the sum in the coalesced first dimension of the 3D array
__global__ void DAS_3DSUM(double* matrix_3d, double* matrix_2d) {

    extern __shared__ double shared_matrix_data[];

    int thread_vector   = threadIdx.x; //1:256
    int all_threads     = blockIdx.x * blockDim.x*2 + threadIdx.x; //1:256 + 1:160k*256*2

    //first add during global load
    shared_matrix_data[thread_vector] = matrix_3d[all_threads] + matrix_3d[all_threads + blockDim.x];
    __syncthreads();
    
    //interleaved addition for 2 interations (no loop = no overhead)
    if(thread_vector < 128) {
        shared_matrix_data[thread_vector] += shared_matrix_data[thread_vector + 128];} __syncthreads();
    if(thread_vector < 64) {
        shared_matrix_data[thread_vector] += shared_matrix_data[thread_vector + 64];} __syncthreads();

    //run the unwrapped warp reduction function
    if(thread_vector < 32) warpreduce(shared_matrix_data, thread_vector);

    //copy to global memory
    if(thread_vector == 0) matrix_2d[blockIdx.x] = shared_matrix_data[0];
}

//wrapper into MATLAB
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

    mxInitGPU();
    const char*         const errId = "parallel:gpu:gpudasindex:InvalidInput"; //error message

    //Input Error handling
    if (nrhs!=2) {
        mexErrMsgIdAndTxt(errId, "Expected 2 inputs: Index Matrix M (3D), rfdata matrix (2D)");
    } 
    if (nlhs!=1 && nlhs!=0) {
        mexErrMsgIdAndTxt(errId, "Invalid number of output arguments, only 1 allowed.");
    }
    if (!mxIsGPUArray(prhs[0])) {
        mexErrMsgIdAndTxt(errId,"Index Matrix M must be a gpuArray.");
    }
    if (!mxIsGPUArray(prhs[1])) {
        mexErrMsgIdAndTxt(errId,"rfdata must be a gpuArray.");
    }

    // Load mx gpu array objects
    const mxGPUArray*   M           = mxGPUCreateFromMxArray(prhs[0]);
    const mxGPUArray*   rfdata      = mxGPUCreateFromMxArray(prhs[1]);
    const mwSize*       dims_3      = mxGPUGetDimensions(M);
    const mwSize        rfdims      = mxGPUGetNumberOfDimensions(rfdata);
    const mwSize        mdims       = mxGPUGetNumberOfDimensions(M);

    if (mdims != 3) {
        mexErrMsgIdAndTxt(errId, "Index Matrix M must have 3 dimensions");
    }

    const mwSize        xysz[2]     = {dims_3[0], dims_3[1]};

    sizes sz1 = {
    //define block and thread sizes
   
    //Get the number of x/y/z threads
    dims_3[1],
    dims_3[0],
    dims_3[2],
    dims_3[0] * dims_3[1],
    dims_3[2] * dims_3[0] * dims_3[1]
    };

    const mwSize* mnum = &sz1.M_numel;

    //Verify that inputs are correct classes before extracting the pointer.
    if (rfdims != 2) {
        mexErrMsgIdAndTxt(errId, "rfdata input must have 2 dimensions (more than 1 frame per function call is not supported)");
    }
    if (mxGPUGetClassID(M) != mxUINT32_CLASS) {
        mexErrMsgIdAndTxt(errId, "Index Matrix M must have class 'uint32'");
    }
    if (mxGPUGetClassID(rfdata) != mxDOUBLE_CLASS) {
        mexErrMsgIdAndTxt(errId, "rfdata must have class 'double'");
    } 
    if(dims_3[2] != 512) {
        mexErrMsgIdAndTxt(errId, "Incorrect number of sensors detected. The third dimension of the index matrix should be the number of sensors, 512. Contact the developer at caso.nathan@gmail.com if you are trying to use a different number of sensors.");
    }

    // Extract a pointer to the input data on the device.
    const unsigned int*     M_dvc       = (const unsigned int*) (mxGPUGetDataReadOnly(M));
    const double*           rfdata_dvc  = (const double*)       (mxGPUGetDataReadOnly(rfdata));

    // GPU Arrays and device-side pointers for the 3-D image and 2-D image
    mxGPUArray* img_3d      = mxGPUCreateGPUArray(1, mnum, mxDOUBLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
    double*     img3_dvc    = (double*)(mxGPUGetData(img_3d));

    mxGPUArray* img_2d      = mxGPUCreateGPUArray(2, xysz, mxDOUBLE_CLASS, mxREAL, MX_GPU_INITIALIZE_VALUES);
    double*     img2_dvc    = (double*)(mxGPUGetData(img_2d));

    //Call the DAS Indexing kernel
    // DAS_Index_GPU<<<block, threads>>>(M_dvc, rfdata_dvc, img3_dvc, img2_dvc, sz1);
    DAS_Index_GPU <<<sz1.xysize, 512>>> (M_dvc, rfdata_dvc, img3_dvc);
    cudaDeviceSynchronize();

    //Call the flattening kernel
    DAS_3DSUM <<<sz1.xysize, 256, 256*sizeof(double)>>>(img3_dvc, img2_dvc);
    cudaDeviceSynchronize();
    
    //Get the result as a gpuArray
    plhs[0] = mxGPUCreateMxArrayOnGPU(img_2d); //takes about 38 microsecs

    //cleanup (takes about 8 microseconds)
    mxGPUDestroyGPUArray(M);
    mxGPUDestroyGPUArray(rfdata);
    mxGPUDestroyGPUArray(img_3d);
    mxGPUDestroyGPUArray(img_2d);
}
