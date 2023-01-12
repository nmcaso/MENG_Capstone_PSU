#include                "Header Files\dataset.hpp"

using namespace         std;

Dataset::Dataset        (const char* fpath) {
    filepath            = fpath;
    load_dataset        ();
}

Dataset::~Dataset       () {}

SensorArray::SensorArray(const char* fpath2) {
    filepath            = fpath2;
    load_sensorarray    ();
}

SensorArray::~SensorArray() {}

void Dataset::load_dataset() {
    //open the the matlab file with the DAS structure in it using the MATLAB API:
    MATFile* matlabfile = matOpen(filepath, "r");
    
    int                 ndir;
    char** dir          = matGetDir(matlabfile, &ndir);

    const char* aa      = dir[0];
    eletemp             = matGetVariable(matlabfile, aa);
    const char* ab      = dir[1];
    fstemp              = matGetVariable(matlabfile, ab);
    const char* ac      = dir[2];
    rfdata              = matGetVariable(matlabfile, ac);
    const char* ad      = dir[3];
    ctemp               = matGetVariable(matlabfile, ad);
    
    int num1            = matClose(matlabfile);

    //get pointers to the mxArrays
    rfptr               = (float*) mxGetPr(rfdata);

    double* temptr1     = (double*) mxGetPr(eletemp);
    double* temptr2     = (double*) mxGetPr(fstemp);
    double* temptr3     = (double*) mxGetPr(ctemp);

    ele                 = *temptr1;
    fs                  = *temptr2;
    c                   = *temptr3*1000.00;
    dt                  = 1/(fs*1000000.00);

    rfdata_rows         = mxGetM(rfdata);
    rfdata_cols         = mxGetN(rfdata);
}

void SensorArray::load_sensorarray() {

    //open the the matlab file and import only x0 and z0
    MATFile* matlabfile2= matOpen(filepath, "r");
    
    int                 ndir;
    char** dir          = matGetDir(matlabfile2, &ndir);

    const char* ae      = dir[4];
    x0                  = matGetVariable(matlabfile2, ae);
    const char* af      = dir[5];
    z0                  = matGetVariable(matlabfile2, af);
    int num2            = matClose(matlabfile2);

    //get pointers to the mxArrays
    x0ptr               = (float*) mxGetPr(x0);
    z0ptr               = (float*) mxGetPr(z0);

    x0_sz               = mxGetN(x0);
    z0_sz               = mxGetN(z0);
}