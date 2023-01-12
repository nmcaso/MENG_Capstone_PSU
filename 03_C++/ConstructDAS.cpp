#include "Header Files\ConstructDAS.hpp"

using namespace std;

int main() {
    //set input path and number of averages.
    const char* winpath = "C:/Users/cason/OneDrive/Documents/PSU/Project/03_C++/LeafData";
    // const char* linpath = "/mnt/c/Users/cason/OneDrive/Documents/PSU/Project/03_C++/LeafData/leaf.mat";
    int                 num_avgs = 1;

    const char* usepath = winpath;

    //Initialize the classes ---------------------------------------------------------------------------------
    auto start          = std::chrono::high_resolution_clock::now();
    Dataset             leafdata(usepath);
    ImageArea           leafarea(-0.02,0.02,-0.02,0.02,0.0001);
    SensorArray         array1(usepath);
    DAS                 indmatdas, sparsemmdas;
    IndexMatrix         indexer(leafarea,array1,leafdata);
    COOMatrix           coo_sparse_mat(leafarea,array1,leafdata);

    //stop the clock, print the performance time
    auto end            = std::chrono::high_resolution_clock::now();
    auto duration       = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    cout                << "Overhead Time: " << static_cast<float>(duration.count()) << " ms" << endl;

    //DAS by index matrix method -----------------------------------------------------------------------------
    auto start2         = std::chrono::high_resolution_clock::now();
    //cout << "Enter number of DAS iterations: "; cin >> num_avgs; cout << "Calculating..." << endl;
    
    for(int i = 0; i < num_avgs; i++) {
        indmatdas.DAS_Index(leafdata,leafarea,indexer,array1);
    }
    auto end2           = std::chrono::high_resolution_clock::now();
    auto duration2      = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
    cout                << "Average CPU C++ Computing time for 3D INDEX is: " << static_cast<float>(duration2.count())/static_cast<float>(num_avgs) << " ms" << endl;

    //DAS by Sparse Matrix Multiplication ---------------------------------------------------------------------
    auto start3         = std::chrono::high_resolution_clock::now();
    cout << "Enter number of DAS iterations: "; cin >> num_avgs; cout << "Calculating..." << endl;
    
    for(int i = 0; i < num_avgs; i++) {
        sparsemmdas.DAS_COO_SPMULT(leafdata,leafarea,coo_sparse_mat,array1,true);
    }
    auto end3           = std::chrono::high_resolution_clock::now();
    auto duration3      = std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3);
    cout                << "Average CPU C++ Computing time for MMULT is: " << static_cast<float>(duration3.count())/static_cast<float>(num_avgs) << " ms" << endl;

    //DAS Index on GPU! ---------------------------------------------------------------------------------------
    cout                << "Calculating on GPU...";
    // double* reconarray2 = new double [indexer.M_cols*indexer.M_rows]();
    double* reconarray2 = new double [indexer.M_numel] (); //use this if your're outputting the index matrix from the GPU
    indexgpuwrapper     (leafdata, indexer, reconarray2);

    //Export the reconstructed array to a binary file -------------------------------------------------------------------
    ofstream            ind_cpu_out, ind_gpu_out, spm_cpu_out;
    ind_cpu_out.open    ("cpp_indexCpuOut.bin");
    ind_gpu_out.open    ("cpp_indexGpuOut.bin");
    spm_cpu_out.open    ("cpp_spmltCpuOut.bin");
       
    int xysize          = leafarea.x_arr_sz*leafarea.y_arr_sz;
    ind_cpu_out.write   (reinterpret_cast<char*>(indmatdas.recon_array), 160000*8);
    
    for(int i = 0; i < xysize; ++i) {
    ind_cpu_out.write   (reinterpret_cast<char*>(indmatdas.recon_array+i), 8);
    spm_cpu_out.write   (reinterpret_cast<char*>(sparsemmdas.recon_array+i), 8);
    ind_gpu_out.write   (reinterpret_cast<const char *>(reconarray2+i),4);
    }
    
    for(int j = 0; j < indexer.M_numel; ++j) { //use this loop if you're outputting the index matrix
    ind_gpu_out.write   (reinterpret_cast<char*>(indexer.M+j), 4); 
    }

    ind_cpu_out.close   ();
    ind_gpu_out.close   ();
    spm_cpu_out.close   ();

}

