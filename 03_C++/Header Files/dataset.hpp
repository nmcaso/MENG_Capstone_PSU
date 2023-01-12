#include <iostream>
#include <cmath>
#include "Header Files\mat.h"
#include "Header Files\matrix.h"

using namespace     std;

//include guards
#ifndef             _DatasetClass_
#define             _DatasetClass_

//class to represent and load a MATLAB dataset that we might use for DAS.
class Dataset {

public:
//attributes
const char*         filepath;
mxArray*            eletemp;
mxArray*            fstemp;
mxArray*            ctemp;
mxArray*            rfdata;
const char*         *var_names;

//array pointers
float*              rfptr;

//copy the constants to c++ variables
double              ele;
double              c;
double              fs;
double              dt;

//array dimensions
int                 rfdata_rows;
int                 rfdata_cols;

//Constructor
Dataset             (const char* path);
~Dataset            ();

//Method Declarations
void                load_dataset();
};

#endif //for dataset class

//include guards again
#ifndef             _SensArrClass_
#define             _SensArrClass_

// Class to respresent the sensor array configuration (which is really x0 and z0 only). This is very similar to the dataset class
class SensorArray {

public:
//attributes:
const char*         filepath;

//matlab API arrays
mxArray*            x0;
mxArray*            z0;

//memory indices
float*              x0ptr;
float*              z0ptr;

//dimensions
int                 x0_sz;
int                 z0_sz;

//methods:
SensorArray         (const char* path2);
~SensorArray        ();

//void laoder
void                load_sensorarray();
};

#endif //sensor array class guard
