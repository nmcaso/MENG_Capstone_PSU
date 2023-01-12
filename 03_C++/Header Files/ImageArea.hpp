#include <iostream>

using namespace std;

//include guards
#ifndef _ImageAreaClass_
#define _ImageAreaClass_

class ImageArea{
public:

//pointers
float*             xptr;
float*             yptr;

//array dimensions
int                 x_arr_sz;
int                 y_arr_sz;

//constructor
ImageArea(double x_min, double x_max, double y_min, double y_max, double intval);
~ImageArea();
};

#endif