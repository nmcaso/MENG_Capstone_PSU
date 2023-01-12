#include "Header Files\ImageArea.hpp"

using namespace std;

ImageArea::ImageArea(double x_min, double x_max, double y_min, double y_max, double intval) {

    //constants
    x_arr_sz        = (x_max-x_min)/intval;
    y_arr_sz        = (y_max-y_min)/intval;

    //Make pointers
    xptr            = new float [x_arr_sz];
    yptr            = new float [y_arr_sz];
    
    //Assign Pointers
    for (int i = 0; i < x_arr_sz; ++i) {
        *(xptr+i) = x_min + intval*i;
    } 
    for (int i = 0; i < y_arr_sz; ++i) {
        *(yptr+i) = y_min + intval*i;
    }

}

ImageArea::~ImageArea() {};