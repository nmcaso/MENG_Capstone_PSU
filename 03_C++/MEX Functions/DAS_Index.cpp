#include "mex.hpp"
#include "mexAdapter.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>

#define timenow std::chrono::high_resolution_clock::now()
#define dur std::chrono::duration_cast<std::chrono::duration<double>>
typedef std::chrono::high_resolution_clock::time_point fast_time;
typedef std::chrono::duration<double> tspan;
typedef std::ostringstream str;

using namespace matlab::mex;
using namespace matlab::data;


//0.07764 average time for 400x400x500 index matrix
class MexFunction : public Function{

    private:
        
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;

        double* rfdata;
        uint32_t* delay_matrix;
        double* image;

        ArrayDimensions delay_size;
        ArrayDimensions rf_size;
        ArrayDimensions img_size;
        size_t n_transducers;
        size_t frame_size;
        size_t n_frames;
        size_t n_pixels;
        
    public:

        //class constructor
        MexFunction() {
            matlabPtr = getEngine();
        }

        //function to display error messages
        void displayError(std::string errorMessage) {
            ArrayFactory factory;
            matlabPtr -> feval(u"error", 0, std::vector<Array> ({factory.createScalar(errorMessage) }));
        }

        //function to display things inside of C
        void disp(std::string message) {
            ArrayFactory factory;
            matlabPtr -> feval(u"fprintf",0, std::vector<Array>({factory.createScalar(message)}));
        }

        // //function to check class of a function (replicates the 'isa' function in matlab)
        // TypedArray<bool> isa(const Array obj, std::string type){

        //     ArrayFactory factory;
        //     std::vector<Array> args{obj, factory.createCharArray(type)};
        //     TypedArray<bool> result = matlabPtr->feval(u"isa", args);
        //     return result;

        // }
        
        // MEX Gateway Function
        void operator()(ArgumentList outputs, ArgumentList inputs) {
        
            //Argument parsing - 0.000020 s
            checkArguments (outputs, inputs);

            //initialize the array factory that we will inevitaby need: 0.000003s
            ArrayFactory f;

            TypedArray<uint32_t> delay_temp = std::move(inputs[0]); 
            TypedArray<double> rf_temp      = std::move(inputs[1]);
            
            //set some C++ class properties - 0.000001s
            delay_size  = delay_temp.getDimensions();
            rf_size     = rf_temp.getDimensions();
            frame_size  = rf_size[0];
            n_transducers = delay_size[2];
            n_pixels    = delay_size[0]*delay_size[1];
            img_size = {delay_size[0], delay_size[1]};
            n_frames = 1;
            
            //preallocate the image array: 0.000029s
            if (rf_size.size() == 3) {
                n_frames = rf_size[2];
                img_size.push_back(n_frames);
            }
            TypedArray<double> image_temp = f.createArray<double>(img_size);

            //get pointers
            image = &image_temp[0][0][0].operator double &();
            rfdata = &rf_temp[0][0][0].operator double &();
            delay_matrix = &delay_temp[0][0][0].operator uint32_t &();
            
            //run the loop: 0.059672s (50% of time)
            if (n_frames == 1) {
                indexDAS_singleframe();
            } else {
                indexDAS_multiframe();
            }
            
            // 0.000006s
            outputs[0] = image_temp;

            // time block for profiling
            // fast_time t0 = timenow;
            // fast_time t00 = timenow;
            // tspan d1 = dur(t00 - t0);
            // disp("time: " + std::to_string(d1.count()) + "\n");
            
        }

        //function to parse arguments and issue basic errors
        void checkArguments(ArgumentList outputs, ArgumentList inputs) {
            //check in/out arg numbers
            if (inputs.size() != 2) {
                displayError("Two inputs: DelayMatrix object and Dataset object required.");
            }
            if (outputs.size() > 1) {
                displayError("Invalid number of output arguments; the maximum is 1");
            }
            if (inputs[0].getType() != ArrayType::UINT32) {
                displayError("Bad delay matrix type, expected class uint32");
            }
            if (inputs[0].getDimensions().size()!= 3) {
                displayError("Expected a 3D Array for the delay matrix.");
            }
            if (inputs[1].getType() != ArrayType::DOUBLE) {
                displayError("Bad rfdata type, expected class double");
            }
            if (inputs[1].getDimensions().size() >3 | inputs[1].getDimensions().size() < 2) {
                displayError("Expected a 2 or 3-dimensional rfdata, where the size is [n_samples, n_transducers, n_frames]");
            }
            if (inputs[1].getDimensions().size() == 3 & inputs[1].getDimensions()[2] == 1) {
                displayError("The frame dimension cannot include singleton dimensions! Use squeeze() in MATLAB first to get rid of them.");
            }

        }

        void indexDAS_multiframe() {
            int pixel_frame_size;
            int rf_frame_size;

            for (int frame = 0; frame < n_frames; ++frame) {
               rf_frame_size = frame*frame_size*n_transducers;
               pixel_frame_size = frame*n_pixels;
            
                for (int sens = 0; sens < n_transducers; ++sens) {
                    for (int pixel = 0; pixel < n_pixels; ++ pixel) {
                        *(image + pixel + pixel_frame_size) += *(rfdata + *(delay_matrix + pixel + sens*n_pixels) + rf_frame_size);
                    }
                }
            }

        }

        void indexDAS_singleframe() {

            for (int sens = 0; sens < n_transducers; ++sens) {
                    for (int pixel = 0; pixel < n_pixels; ++ pixel) {
                        *(image + pixel) += *(rfdata + *(delay_matrix + pixel + sens*n_pixels));
                    }
                }
        }

};