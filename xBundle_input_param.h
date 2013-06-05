//
//  xBundle_input_param.h
//  BundleEx
//
//  Created by arindam bhattacharya on 5/29/13.
//
//

#ifndef BundleEx_xBundle_input_param_h
#define BundleEx_xBundle_input_param_h
#include <string>
using namespace std;
class INPUT_PARAMS {
public:
  string inputFileName;
  bool inputFileUI;
  bool inputTh;
  bool computeHessianBasedComp;
  float sigma;
  int thresh;
  int eigenVecNum;
  bool degugMode;
  INPUT_PARAMS(){
    inputFileUI = false;
    inputTh = false;
    computeHessianBasedComp = false;
    eigenVecNum=0;
    sigma=3.0f;
    degugMode = false;
  }
};

#endif

