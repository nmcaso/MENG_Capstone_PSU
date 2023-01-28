/* ========================================================================
 * phonebook.cpp
 * example for illustrating how to manipulate structure.
 *
 * takes a (MxN) structure matrix which has first field as 
 * character array(name), and second field as scalar double (phone number). 
 * This function returns a new structure (1x1)containing following fields: 
 * for character array input, it will be (MxN) cell array; 
 * and for numeric double (noncomplex, scalar) input, it will be (MxN)
 * cell array where each field is numeric array of type double.
 *
 * Build : from MATLAB
 *         >> mex phonebook.cpp
 * Usage with example : from MATLAB
 *         >> friends(1).name = 'Jordan Robert';
 *         >> friends(1).phone = 3386;
 *         >> friends(2).name = 'Mary Smith';
 *         >> friends(2).phone = 3912;
 *         >> friends(3).name = 'Stacy Flora';
 *         >> friends(3).phone = 3238;
 *         >> friends(4).name = 'Harry Alpert';
 *         >> friends(4).phone = 3077;
 *         >> phonebook(friends)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2017 The MathWorks, Inc.
 *=======================================================================*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <string>
#include <memory>

using namespace matlab::mex;
using namespace matlab::data;

class MexFunction : public Function {
private:
  std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
public:
  /* Constructor for the class. */
  MexFunction()
  {
    matlabPtr = getEngine();
  }
  
  /* Helper function to print output string on MATLAB command prompt. */
  void displayOnMATLAB(std::ostringstream stream)
  {
    ArrayFactory factory;
    matlabPtr->feval(u"fprintf", 0, std::vector<Array>
            ({ factory.createScalar(stream.str())}));
  }
  
  /* Helper function to generate an error message from given string,
   * and display it over MATLAB command prompt.
   */
  void displayError(std::string errorMessage)
  {
    ArrayFactory factory;
    matlabPtr->feval(u"error", 0, std::vector<Array>({
      factory.createScalar(errorMessage) }));
  }
  
  /* Helper function to information about an empty field in the structure. */
  void emptyFieldInformation(std::string fieldName, size_t index)
  {
    std::ostringstream stream;
    stream<<"Field: "<<std::string(fieldName)<<" of the element at index: "
          <<index+1<<" is empty."<<std::endl;
    displayOnMATLAB(std::move(stream));
  }
  
  /* Helper function to information about an invalid field in the structure. */
  void invalidFieldInformation(std::string fieldName, size_t index)
  {
    std::ostringstream stream;
    stream<<"Field: "<<std::string(fieldName)<<" of the element at index: "
          <<index+1<<" contains wrong value."<<std::endl;
    displayOnMATLAB(std::move(stream));
  }
  
  
  /* This is the gateway routine for the MEX-file. */
  void 
  operator()(ArgumentList outputs, ArgumentList inputs) { 
    
    checkArguments (outputs,inputs);
    
    ArrayFactory factory;   
    StructArray const matlabStructArray = inputs[0];
    checkStructureElements(matlabStructArray);
    auto fields = matlabStructArray.getFieldNames();
    size_t total_num_of_elements = matlabStructArray.getNumberOfElements();   
    CellArray phoneNumberStringArray = 
            factory.createCellArray({ 1,total_num_of_elements });
    CellArray nameStringArray = 
            factory.createCellArray({ 1,total_num_of_elements });
    std::vector<std::string> fieldNames(fields.begin(), fields.end());
    
    StructArray result = factory.createStructArray({ 1,1 },fieldNames );
    
    /* Walk through each structure element. */
    for (size_t entryIndex=0; entryIndex<total_num_of_elements; entryIndex++) {
      /*Copy data from structure array to cell array. */
      Array const structField1 = 
              matlabStructArray[entryIndex][fieldNames[0]];
      Array const structField2 = 
              matlabStructArray[entryIndex][fieldNames[1]];
      phoneNumberStringArray[entryIndex] = structField1;
      nameStringArray[entryIndex] = structField2;
    }
    result[0][fieldNames[0]] = phoneNumberStringArray;
    result[0][fieldNames[1]] = nameStringArray;
    outputs[0] = result;
  }
  
  /* Make sure that the passed structure has valid data. */
  void checkStructureElements(StructArray const & matlabStructArray)
  {
    std::ostringstream stream;
    size_t nfields = matlabStructArray.getNumberOfFields();
    auto fields = matlabStructArray.getFieldNames();
    size_t total_num_of_elements = matlabStructArray.getNumberOfElements();
    std::vector<std::string> fieldNames(fields.begin(), fields.end());
    
    /* Produce error if structure has more than 2 fields. */
    if(nfields != 2) {
      displayError("Struct must consist of 2 entries."
                   "(First: char array, Second: numeric double scalar).");
    }
    
    /* Walk through each structure element. */
    for (size_t entryIndex=0; entryIndex<total_num_of_elements; entryIndex++) {
      const Array structField1 = 
              matlabStructArray[entryIndex][fieldNames[0]];
      const Array structField2 = 
              matlabStructArray[entryIndex][fieldNames[1]];
      
      /* Produce error if name field in structure is empty. */
      if (structField1.isEmpty()) {
        emptyFieldInformation(fieldNames[0],entryIndex);
        displayError("Empty fields are not allowed in this program." 
                     "This field must contain character array.");
      }
      
      /* Produce error if phone number field in structure is empty. */
      if(structField2.isEmpty()) {
        emptyFieldInformation(fieldNames[1],entryIndex);
        displayError("Empty fields are not allowed in this program." 
                     "This field must contain numeric double scalar.");
      }
      
      /* Produce error if name is not a valid character array. */
      if(structField1.getType()!= ArrayType::CHAR) {
        invalidFieldInformation(fieldNames[0],entryIndex);
        displayError("This field must contain character array.");
      }
      /* Produce error if phone number is not a valid double scalar. */
      if (structField2.getType() != ArrayType::DOUBLE
          || structField2.getNumberOfElements() != 1) {
        invalidFieldInformation(fieldNames[1],entryIndex);
        displayError("This field must contain numeric double scalar.");
      }
    }
  }
  
  /* This function makes sure that user has provided structure as input,
   * and is not expecting more than one output in results.
   */
  void checkArguments(ArgumentList outputs, ArgumentList inputs) {
    if (inputs.size() != 1) {
      displayError("One input required.");
    }
    if (outputs.size() > 1) {
      displayError("Too many outputs specified.");
    }
    if (inputs[0].getType() != ArrayType::STRUCT) {
      displayError("Input must be a structure.");
    }
  }
};
