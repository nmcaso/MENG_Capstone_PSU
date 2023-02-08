#include "mex.hpp"
#include "mexAdapter.hpp"
#include <string>
#include <vector>
#include <chrono>

using namespace matlab::mex;
using namespace matlab::data;

#define timenow std::chrono::high_resolution_clock::now()
#define dur std::chrono::duration_cast<std::chrono::duration<double>>
typedef std::chrono::high_resolution_clock::time_point fast_time;
typedef std::chrono::duration<double> tspan;
typedef std::ostringstream str;

//a struct for array sizes to pass into the DAS routine itself
typedef struct array_sizes {
    size_t x_dim;                             // number of rows
    size_t y_dim;                             // number of columns
    size_t z_dim;                             // number of sensor elements
    size_t frame_size;                        // rfdata number of rows
} array_sizes;

//a struct for c and dt to pass into the index matrix function
typedef struct physics_parameters {
    double c;
    double dt;
} physics_parameters;

//alternate form of the square root function
float sqrt1(const float &n) 
{
   static union{int i; float f;} u;
   u.i = 0x5F375A86 - (*(int*)&n >> 1);
   return (int(3) - n * u.f * u.f) * n * u.f * 0.5f;
}

class MexFunction : public Function {

    private:

        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;

    public:

        //function to display error messages
        void displayError(std::string errorMessage) {
            ArrayFactory factory;
            matlabPtr -> feval(u"error", 0, std::vector<Array> ({factory.createScalar(errorMessage) }));
        }

        void disp(std::ostringstream& stream) {
            ArrayFactory factory;
            matlabPtr -> feval(u"fprintf",0, std::vector<Array>({factory.createScalar(stream.str())}));
            stream.str("");
        }
        
        //Class constructor ~1.5e-7 seconds
        MexFunction() { 
            
            matlabPtr = getEngine();
        
        }     
    
        //MEX Gateway function
        void operator() (ArgumentList outputs, ArgumentList inputs) { 
        
            //Check arguments (fast)
            checkArguments (outputs,inputs);

            //Generate CPP MEX arrays from the inputs - 1.05e-5 seconds
            Array dataset(inputs[0]);
            Array imagearea(inputs[1]);
            Array sensorarray(inputs[2]);
            double depth = inputs[3][0];
                       
            //Generate typed arrays from the initial array properties - 2.518e-4 seconds (fast)
            TypedArray<double> x_arr    = matlabPtr -> getProperty(imagearea,u"x_arr");
            TypedArray<double> y_arr    = matlabPtr -> getProperty(imagearea,u"y_arr");
            TypedArray<double> x0       = matlabPtr -> getProperty(sensorarray,u"x0");
            TypedArray<double> z0       = matlabPtr -> getProperty(sensorarray,u"z0");
            TypedArray<double> rfdata   = matlabPtr -> getProperty(dataset,u"rfdata");
            TypedArray<double> c_arr    = matlabPtr -> getProperty(dataset,u"c");
            TypedArray<double> dt_arr   = matlabPtr -> getProperty(dataset,u"dt");

            //collect physical parameters of system speed of sound and sample rate - this and other struct - 2.2e-6 seconds
            physics_parameters sys_config = {
                c_arr[0],
                dt_arr[0],
            };

            //collect image size parameters (pixels in x direction, y direction, number of sensors, and frame size)
            array_sizes img_sizes = {
                x_arr.getNumberOfElements(),
                y_arr.getNumberOfElements(),
                x0.getNumberOfElements(),
                rfdata.getNumberOfElements()/x0.getNumberOfElements(),
            };
            
            //Preallocate arrays for the index matrix and image - 5.2e-5 seconds
            ArrayFactory                    f;
            TypedArray<int> index_matrix    = f.createArray<int>({img_sizes.x_dim, img_sizes.y_dim});
            TypedArray<double> image        = f.createArray<double>({img_sizes.x_dim, img_sizes.y_dim});

            //initiate iterators for x0 and z0 (and an integer for the sensor index) - 4.3e-6
            TypedIterator<double> x0_itr = x0.begin();
            TypedIterator<double> z0_itr = z0.begin();
            int sensor_index = 0;

            //Loop over sensors, creating the index matrix for the sensor and subsequently perform DAS by indexing.
            for(; x0_itr != x0.end(); ++x0_itr, ++z0_itr) {
                makeIndexMatrix(index_matrix, x_arr, y_arr, depth, *x0_itr, *z0_itr, sys_config, img_sizes);
                indexDAS(image, index_matrix, rfdata, img_sizes, sensor_index); //moderate - 0.113525
                ++sensor_index;
            }

            //output the image
            outputs[0] = image;
        }
    
        void checkArguments(ArgumentList outputs, ArgumentList inputs) {
            if (inputs.size() != 4) {
            displayError("Expected 4 inputs: dataset, image area, and sensor array objects, and a z-depth");
            }
            if (outputs.size() > 1) {
            displayError("Invalid number of output arguments; the maximum is 1");
            }
        }

        //Generates a delay matrix for virtual image locations in x_arr and y_arr, for current sensor positions x0 and z0
        void makeIndexMatrix(TypedArray<int>& ind_mat, TypedArray<double>& x_arr, TypedArray<double>& y_arr, 
                    double& z, double& x0, double& z0, physics_parameters& rf_params, array_sizes& sizes) {

            //generate pointer for the index matrix
            TypedIterator<int> iptr = ind_mat.begin();

            //Loop across the x array and y arrays (for the constant values of x0 and z0)
            for(auto& x : x_arr) {
                for(auto& y : y_arr) {

                    //calculate the index
                    float q1 = (float)(x - x0);
                    float q2 = (float)(y - z0);
                    // float q3 = (float)(z);
                    int value = (z + sqrt(q1*q1 + q2*q2)) / rf_params.c / rf_params.dt;

                    //correct any bad indices
                    if (value <= 0) {
                        *iptr = 0;
                    } else if (value > sizes.frame_size - 1) {
                        *iptr = sizes.frame_size - 1;
                    } else {
                        *iptr = value;
                    }

                    //increment the pointer
                    ++iptr;
                }
            }
        }

        //Performs Delay and Sum beamforming for a 2D image.
        void indexDAS(TypedArray<double>& img, TypedArray<int>& ind_mat, TypedArray<double>& rf_data, array_sizes& sizes, int& current_sensor) {
            
            //set up 2 pointers for the delay matrix and rfdata
            TypedIterator<double> img_ptr   = img.begin();
            TypedIterator<double> rf_ptr    = rf_data.begin();

            //loop across all values in the delay index matrix
            for(auto& index : ind_mat) {
                
                //add the contributions of the current sensor's beamformed signal arrays to the image
                *img_ptr += *(rf_ptr + index + (int)(sizes.frame_size*current_sensor));
                
                //increment the image pointer
                ++img_ptr;
            }
        }
};