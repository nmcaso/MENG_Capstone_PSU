#include <iostream>
#include <cmath>
#include <chrono>
#include "Header Files\mat.h"
#include "Header Files\matrix.h"
#include "Header Files\ImageArea.hpp"
#include "Header Files\Dataset.hpp"

using namespace std;

//include guard
#ifndef _IndMatClass_
#define _IndMatClass_

class IndexMatrix {

public:

//attributes
unsigned int* M;
int     M_rows;
int     M_cols;
int     M_depth;
int     M_numel;
bool*   errorvec;

//methods
IndexMatrix(ImageArea area1, SensorArray array1, Dataset data1);
~IndexMatrix();

};

#endif //_IndMatClass_

#ifndef COOMatClass
#define COOMatClass   
class COOMatrix {

public:

//attributes
unsigned int* M;
int*    M_R;
int*    M_C;
int*    M_V;
int*    M_RCV;
int     M_numel;

//methods
COOMatrix(ImageArea area1, SensorArray array1, Dataset data1);
~COOMatrix();
};

#endif //Sparse matrix multiplication version of 

#ifndef ELLMatrix_Class
#define ELLMatrix_Class
class ELLMatrix {

public:

//attributes
int*    indices;
int     col_num;

};

#endif

//include guard
#ifndef _DASClass_
#define _DASClass_

//create a class of mainly functions that will act upon the dataset.
class DAS {

public:

//atributes!
double* recon_array;
int xysize;

//default constructor/destructor:
DAS();
~DAS();

//DAS by indexing in a large index matrix
void DAS_Index(Dataset data2, ImageArea area1, IndexMatrix indmat, SensorArray sensarr);

//DAS by multiplying by a large sparse matrix in coordinate format
void DAS_COO_SPMULT(Dataset data1, ImageArea area1, COOMatrix mmat, SensorArray sensarr, bool vals_all_one);

//DAS by multiplying by sparce matrix in ELLPACK format
void DAS_ELL_SPMULT(Dataset data1, ImageArea area1, ELLMatrix mmat, SensorArray sensarr, bool vals_all_one);

//Static method
double* matrix_vecadd(double* vec1ptr, double* vec2ptr, int vec1sz, int vec2sz);
};

#endif //_DASClass_