#include "Header Files\CUDADAS.cuh"

using namespace std;

 struct sizes {
                        int i;
                        int j;
                        int k;
                        int imax;
                        int jmax;
                        int kmax;
                        int xysize;
    };

//Creates a beamformed image array on the GPU.
__global__ void DAS_Index_GPU(float* rfptr, unsigned int* indmat, double* reconptr, sizes sz) {    
    
    //preallocate
    int x               = threadIdx.x + blockDim.x*blockIdx.x;
    int y               = threadIdx.y + blockDim.y*blockIdx.y;
    int z               = threadIdx.z + blockDim.z*blockIdx.z;
    unsigned int        final_index;

    //memory indexing (works correctly)
     if(x < sz.imax && y < sz.jmax && z < sz.kmax) {
        final_index = *(indmat + x + y*sz.imax + z*sz.xysize);
        *(reconptr+x+sz.imax*y + z*sz.xysize) = *(rfptr+final_index);
     }
}

//This function takes the GPU array form DAS_Index_GPU (in 3D) and does an average along the z-dimension, thereby flattening it to 2 dimensions for image viewing
__global__ void DAS_Index_Flatten(double* recon_3d, double* recon_2d, sizes sz) {

    //preallocate
    int x               = threadIdx.x + blockDim.x*blockIdx.x;
    int y               = threadIdx.y + blockDim.y*blockIdx.y;
    int z               = threadIdx.z + blockDim.z*blockIdx.z;

    // does not coalesce properly, see gpudasindex.cu mex function for the correct function)
    if(x < sz.imax && y < sz.jmax && z < sz.kmax) {
        *(recon_2d+x+sz.imax*y + z*sz.xysize) += *(recon_3d + x + sz.imax*y + z*sz.xysize);
    }

}

//Wrapper function to call from C++ Code that will initialize variables on the GPU and call the kernel with proper syntax.
void indexgpuwrapper(Dataset data2, IndexMatrix indmat, double* cpureconptr) {
    
    struct              sizes sz1;

    //define block and thread sizes
    sz1.imax            = indmat.M_cols;
    sz1.jmax            = indmat.M_rows;
    sz1.kmax            = indmat.M_depth;
    sz1.xysize          = sz1.imax*sz1.jmax;
    
    //Get the number of threads in each direction
    sz1.i               = 16;
    sz1.j               = 16;
    sz1.k               = 4;
        
    //round up to the nearest integer of the dimension length divided by the thread dimension
    int blockx          = indmat.M_rows/sz1.i + (indmat.M_rows % sz1.i != 0);
    int blocky          = indmat.M_cols/sz1.j + (indmat.M_cols % sz1.j != 0);
    int blockz          = indmat.M_depth/sz1.k +(indmat.M_depth % sz1.k != 0);

    //grid and block 3D arrays
    dim3                grid(sz1.i, sz1.j, sz1.k);
    dim3                block(blockx, blocky, blockz);

    //preallocate a GPU pointer and make convenient size variables
    int     rfsize      = data2.rfdata_cols*data2.rfdata_rows;
    double*  recontemp  = new double [indmat.M_numel]();
    double*  reconptr   = new double [sz1.xysize]();
    unsigned int* gpu_index   = new unsigned int [indmat.M_numel];
    float*  gpu_rf      = new float [rfsize];

    // preallocate space on the GPU for the big variables
    cudaMallocManaged   (&reconptr,     indmat.M_numel  *sizeof(double)); //**use if outputting index matrix
    cudaMallocManaged   (&recontemp,    indmat.M_numel  *sizeof(double));
    cudaMallocManaged   (&gpu_rf,       rfsize          *sizeof(float));
    cudaMallocManaged   (&gpu_index,    indmat.M_numel  *sizeof(int));

    //load the big variables onto the GPU
    cudaMemcpy          (gpu_index,     indmat.M,       indmat.M_numel*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy          (gpu_rf,        data2.rfptr,    rfsize*sizeof(float),        cudaMemcpyHostToDevice);

    //set up a number of runs to perform
    int                 numavgs = 1;

    auto start3         = chrono::high_resolution_clock::now();

    //run the DAS algorithm on the GPU a number of times equal to numavgs.
    for(int i = 0; i < numavgs; ++i) {
        DAS_Index_GPU       <<< block, grid >>>     (gpu_rf, gpu_index, reconptr, sz1);
        // DAS_Index_Flatten   <<<block,grid>>>        (recontemp,reconptr,sz1); //comment this line if you're outputting the index matrix
        cudaDeviceSynchronize                       ();
    }

    auto stop3          = chrono::high_resolution_clock::now();
    auto duration3      = std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3);
    cout << "GPU Computing Time: " << static_cast<float>(duration3.count())/static_cast<float>(numavgs) << " ms" << endl;

    auto start4         = chrono::high_resolution_clock::now();
    // Copy the GPU array back to the CPU
    cudaMemcpy          (cpureconptr, reconptr, indmat.M_numel*sizeof(double), cudaMemcpyDeviceToHost); //copies the index matrix out
    cudaFree            (reconptr);
    cudaFree            (gpu_rf);
    cudaFree            (gpu_index);
    cudaFree            (recontemp);
    auto stop4          = chrono::high_resolution_clock::now();
    auto duration4      = std::chrono::duration_cast<std::chrono::microseconds>(stop4 - start4);
    cout << "GPU back to CPU and Free Memory time: " << static_cast<float>(duration4.count()) << " ms" << endl;
}