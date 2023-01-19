#include "C:/Program Files/MATLAB/R2022b/extern/include/mex.h"
#include "C:/Program Files/MATLAB/R2022b/toolbox/parallel/gpu/extern/include/gpu/mxGPUArray.h"
#include <chrono>

#define timer       std::chrono::high_resolution_clock
#define timertime   std::chrono::high_resolution_clock::time_point
#define timecast    std::chrono::duration_cast<std::chrono::microseconds>
#define timesecs    std::chrono::microseconds

//a convenient structure
typedef struct sizes {
            size_t i;
            size_t j;
            size_t k;
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
__global__ void DAS_Index_GPU(const unsigned int* indmat, const double* rfdata, double* img3d, double* img, sizes sz) {    
    // define device variables
    int x = threadIdx.x + blockDim.x * blockIdx.x;
    int y = threadIdx.y + blockDim.y * blockIdx.y;
    int z = threadIdx.z + blockDim.z * blockIdx.z;

    //memory index
    if(z < sz.kmax && y < sz.jmax && x < sz.imax) {
        int M_index             = x + y * sz.imax + z * sz.xysize;
        int Img_Index           = z + x * sz.kmax + y * sz.kmax*sz.imax;
        *(img3d + Img_Index) = *(rfdata + *(indmat + M_index));
    }
}

__global__ void DAS_3DSUM(double* matrix_3d, double* matrix_2d, sizes sz) {

    extern __shared__ double shared_matrix_data[];

    int thread_vector   = threadIdx.x; //1:512
    int all_threads     = blockIdx.x * blockDim.x*2 + threadIdx.x; //1:512 + 1:160k*512

    //first add during global load, we've done the first add from 512 - 256 elements
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
    if (nrhs!=5) {
        mexErrMsgIdAndTxt(errId, "Expected 5 inputs: M, rfdata, thread/block x, thread/block y, thread/block z");
    } 
    if (nlhs!=1 && nlhs!=0) {
        mexErrMsgIdAndTxt(errId, "Invalid number of output arguments, only 1 allowed.");
    }
    if (!mxIsGPUArray(prhs[0])) {
        mexErrMsgIdAndTxt(errId, "Index Matrix M must be a GPU array.");
    } 
    if (!mxIsGPUArray(prhs[1])) {
        mexErrMsgIdAndTxt(errId, "rfdata must be a GPU array.");
    } 
    if (!mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4])) {
        mexErrMsgIdAndTxt(errId, "Block dimensions must be of datatype 'double'");
    }

    // Load mx gpu array objects
    const mxGPUArray*   M           = mxGPUCreateFromMxArray(prhs[0]);
    const mxGPUArray*   rfdata      = mxGPUCreateFromMxArray(prhs[1]);
    const mwSize*       dims_3      = mxGPUGetDimensions(M);
    const mwSize        xysz[2]     = {dims_3[0], dims_3[1]};

    sizes sz1 = {
    //define block and thread sizes
    (unsigned long long)*(double*)mxGetData(prhs[2]),
    (unsigned long long)*(double*)mxGetData(prhs[3]),
    (unsigned long long)*(double*)mxGetData(prhs[4]),
    
    //Get the number of x/y/z threads
    dims_3[1],
    dims_3[0],
    dims_3[2],
    dims_3[0] * dims_3[1],
    dims_3[2] * dims_3[0] * dims_3[1]
    };

    //round up to the nearest integer of the dimension length divided by the thread   dimension
    int blockx = sz1.imax/sz1.i + (sz1.imax % sz1.i != 0);
    int blocky = sz1.jmax/sz1.j + (sz1.jmax % sz1.j != 0);
    int blockz = sz1.kmax/sz1.k + (sz1.kmax % sz1.k != 0);

    //grid and block 3D arrays
    dim3    threads(sz1.i, sz1.j, sz1.k);
    dim3    block(blockx, blocky, blockz);
    size_t  total_threads = sz1.i * sz1.j * sz1.k;

    const mwSize* mnum = &sz1.M_numel;
    //Verify that inputs are correct classes before extracting the pointer.
    if (mxGPUGetClassID(M) != mxUINT32_CLASS) {
        mexErrMsgIdAndTxt(errId, "Index Matrix M must have class 'uint32'");
    }
    if (mxGPUGetClassID(rfdata) != mxDOUBLE_CLASS) {
        mexErrMsgIdAndTxt(errId, "rfdata must have class 'double'");
    } 
    if (total_threads > 1024 || total_threads <= 1) {
        mexErrMsgIdAndTxt(errId, "Block dimensions product must be between 1 and 1024");
    }

    // Extract a pointer to the input data on the device.
    const unsigned int*     M_dvc       = (const unsigned int*) (mxGPUGetDataReadOnly(M));
    const double*           rfdata_dvc  = (const double*)       (mxGPUGetDataReadOnly(rfdata));

    // Create a GPUArray to hold the result and get its underlying pointer.
    mxGPUArray* img_3d      = mxGPUCreateGPUArray(1, mnum, mxDOUBLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
    double*     img3_dvc    = (double*)(mxGPUGetData(img_3d));

    mxGPUArray* img_2d      = mxGPUCreateGPUArray(2, xysz, mxDOUBLE_CLASS, mxREAL, MX_GPU_INITIALIZE_VALUES);
    double*     img2_dvc    = (double*)(mxGPUGetData(img_2d));

    timertime start8 = timer::now();
    // Call the kernel
    DAS_Index_GPU<<<block, threads>>>(M_dvc, rfdata_dvc, img3_dvc, img2_dvc, sz1);
    cudaDeviceSynchronize();

    //Call the other kernel
    DAS_3DSUM <<<sz1.xysize, 256, 256*sizeof(double)>>>(img3_dvc, img2_dvc, sz1);
    cudaDeviceSynchronize();

    timertime stop8 = timer::now();
    timesecs indexgpu = timecast(stop8 - start8);
    printf("Kernel Index: %d\n",indexgpu);
    
    //Get the result as a gpuArray
    plhs[0] = mxGPUCreateMxArrayOnGPU(img_2d); //takes about 38 microsecs

    //cleanup
    mxGPUDestroyGPUArray(M);
    mxGPUDestroyGPUArray(rfdata);
    mxGPUDestroyGPUArray(img_3d);
    mxGPUDestroyGPUArray(img_2d);
    //cleanup takes about 8 microsecs
}
