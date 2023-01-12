#include "Header Files/DAS.hpp"

using namespace std;

DAS::DAS() {}
DAS::~DAS() {}

//function that replicates Matlab's vector addition for 2 vectors, except it outputs a 1D vector of dimension 1 x m*n.
double* matrix_vecadd(double* vec1ptr, double* vec2ptr, int vec1sz, int vec2sz) {
    double* cvec = new double [vec1sz*vec2sz];
        for(int i = 0; i < vec1sz; ++i) {
        for(int j = 0; j < vec2sz; ++j){
            *(cvec+vec1sz*i+j) = *(vec1ptr+i) + *(vec2ptr+j);
        }
    }
    return cvec;
}

IndexMatrix::~IndexMatrix() {}

//Index matrix class generates the entire 3D-ish index matrix
IndexMatrix::IndexMatrix(ImageArea area1, SensorArray array1, Dataset data1) {

    M_rows          = area1.x_arr_sz;
    M_cols          = area1.y_arr_sz;
    M_depth         = array1.x0_sz;
    M_numel         = M_depth*M_cols*M_rows;

    //preallocate 2 vectors and a matrix
    double* vec1ptr = new double [area1.x_arr_sz];
    double* vec2ptr = new double [area1.y_arr_sz];
    M               = new unsigned int [area1.x_arr_sz*area1.y_arr_sz*array1.x0_sz];

    //loop to create a long array of numbers.
    for(int num = 0; num < array1.x0_sz; ++num) {
    
        //Get x0 and z0 for a sensor, the factor of 0.001 converts from m to mm
        double x0val = *(array1.x0ptr + num)*0.001;
        double z0val = *(array1.z0ptr + num)*0.001;

        //find the matrix indices
        for(int i = 0; i < area1.x_arr_sz; ++i) {
            *(vec1ptr+i) = pow(*(area1.xptr+i) - x0val,2);
        }
        for(int i = 0; i < area1.y_arr_sz; ++i) {
            *(vec2ptr+i) = pow(*(area1.yptr+i) - z0val,2);
        }

        //interpolate to every index
        double* vec3ptr = matrix_vecadd(vec1ptr, vec2ptr, area1.x_arr_sz, area1.y_arr_sz);

        //Apply square root, divide by the speed of sound, and phase appropriately for the size of the rfdata array.
        for(int i = 0; i < area1.x_arr_sz*area1.y_arr_sz; ++i) {
            *(M+i+num*area1.x_arr_sz*area1.y_arr_sz) = round(sqrt(*(vec3ptr+i))*1/data1.c*1/data1.dt)+num*data1.rfdata_rows-1;
        }
    }
}

COOMatrix::~COOMatrix() {}

COOMatrix::COOMatrix(ImageArea area1, SensorArray array1, Dataset data1) {

    const int xysize    = area1.x_arr_sz*area1.y_arr_sz;
    M_numel             = xysize*array1.x0_sz;

    //the Rows is given by a repeating sequence of 0 to the number of sensors
    M_R                 = new int[M_numel];
    for(int i = 0; i < M_numel; ++i) {
        *(M_R+i) = i%xysize;
    }

    //preallocate 2 vectors and a matrix
    double* vec1ptr     = new double [area1.x_arr_sz];
    double* vec2ptr     = new double [area1.y_arr_sz];
    M_C                 = new int [xysize*array1.x0_sz];

    //loop to create Column matrix (same as the IndexMatrix!)
    for(int num = 0; num < array1.x0_sz; ++num) {
    
        //Get x0 and z0 for a sensor, the factor of 0.001 converts from m to mm
        float x0val = *(array1.x0ptr + num)*0.001;
        float z0val = *(array1.z0ptr + num)*0.001;

        //find the matrix indices
        for(int i = 0; i < area1.x_arr_sz; ++i) {
            *(vec1ptr+i) = pow(*(area1.xptr+i) - x0val,2);
        }
        for(int i = 0; i < area1.y_arr_sz; ++i) {
            *(vec2ptr+i) = pow(*(area1.yptr+i) - z0val,2);
        }

        //interpolate to every index
        double* vec3ptr = matrix_vecadd(vec1ptr, vec2ptr, area1.x_arr_sz, area1.y_arr_sz);

        //Apply square root, divide by the speed of sound, and phase appropriately for the size of the rfdata array.
        for(int i = 0; i < xysize; ++i) {
            *(M_C+i+num*xysize) = round(sqrt(*(vec3ptr+i))*1/data1.c*1/data1.dt)+num*data1.rfdata_rows-1;
        }

    }

    //values = 1 for every Row and Col option
    M_V             = new int[xysize*array1.x0_sz];
    for(int i = 0; i < xysize*array1.x0_sz; ++i) {
        *(M_V+i) = 1;
    }

     //put everything together into 1 long vector
    M_RCV               = new int[M_numel*2];
    for(int i = 0; i < M_numel; ++i){
        *(M_RCV + 2*i+0) = *(M_R+i);
        *(M_RCV + 2*i+1) = *(M_C+i);
        // *(M_RCV + 3*i+2) = *(M_V+i);
    }

}

//Creates a beamformed image array
void DAS::DAS_Index(Dataset data2, ImageArea area1, IndexMatrix indmat, SensorArray sensarr) {    
    
    //preallocate
    xysize              = area1.y_arr_sz*area1.x_arr_sz;
    recon_array         = new double [xysize]();

    //memory index through each value to beamform the image
    for (int i = 0; i < sensarr.x0_sz; ++i) {
        for(int j = 0; j < xysize; ++j) {
            *(recon_array+j) += *(data2.rfptr+*(indmat.M+j+i*xysize));
        }
    }

    //divide by the number of sensors
    for (int i = 0; i < xysize; ++i) {
        *(recon_array+i) /= sensarr.x0_sz;
    }
}

void DAS::DAS_COO_SPMULT(Dataset data1, ImageArea area1, COOMatrix mmat, SensorArray sensarr, bool vals_all_one) {

    //preallocate
    xysize              = area1.y_arr_sz*area1.x_arr_sz;
    recon_array         = new double [xysize]();

    //memory index
    for (int i = 0; i < mmat.M_numel; ++i) {
            *(recon_array+*(mmat.M_RCV+2*i)) += *(data1.rfptr + *(mmat.M_RCV + 2*i+1)); 
        }

}

void DAS::DAS_ELL_SPMULT(Dataset data1, ImageArea area1, ELLMatrix mmat, SensorArray sensarr, bool vals_all_one) {

    //preallocate
    xysize              = area1.y_arr_sz*area1.x_arr_sz;
    recon_array         = new double[xysize]();

    // mmat.M


}