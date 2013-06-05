/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
//  Software Guide : BeginLatex
//
//  The following code is an implementation of a small Insight
//  program. It tests including header files and linking with ITK
//  libraries.
//
//  Software Guide : EndLatex

// Software Guide : BeginCodeSnippet


#include <iostream>
using namespace std;

#include "xBundle_input_param.h"
#include "xBundle_hessian.h"




void usage(char* argv[]){
  cout <<"Usage: " << argv[0] <<" Options"<<endl;
  cout <<"Options :" <<endl;
  cout <<"\t -i INPUT_FILENAME // required"<<endl;
  cout <<"\t -thres <int> //required"<<endl;
  cout <<"\t -hessian"<<endl;
  cout <<"\t -sigma <float> default 3.0"<<endl;
  cout <<"\t -eig <int> default 2"<<endl;
  
}

// read the parameters from the users
int readParameters(int argc, char* argv[], INPUT_PARAMS * input)
{
  int ii = 1;
  while (ii < argc) {
    cout << " " << argv[ii] << endl;
    if (strcmp(argv[ii], "-hessian") == 0)
      {
      cout << "Hessian !!!!" << endl;
      input->computeHessianBasedComp=true;
      }
    else if (strcmp(argv[ii], "-i") == 0)
      {
      ii++;
      cout <<"FileName " << argv[ii] << endl;
      input->inputFileUI=true;
      }
    else if (strcmp(argv[ii], "-s") == 0)
      {
      ii++;
      cout <<"sigma " << argv[ii] << endl;
      input->sigma=atof(argv[ii]);
      }
    else if (strcmp(argv[ii], "-thres") == 0)
      {
      ii++;
      input->thresh=atoi(argv[ii]);
      cout <<"thres " <<  input->thresh << endl;
      input->inputTh=true;
      }
    else if (strcmp(argv[ii], "-eig") == 0)
      {
      ii++;
      input->eigenVecNum=atoi(argv[ii]);
      cout <<"eig vec num " <<  input->eigenVecNum << endl;
      }
    else if (strcmp(argv[ii], "-debug") == 0)
      {
      input->degugMode = true;
      }
    else
      {
      cout << "Option[" << argv[ii] << "]does not match anything "<< endl;
      usage(argv);
      return 10;
    }
    ii++;
  }
  if (!input->inputFileUI && !input->inputTh)
    cout <<" no input file file or threshold" << endl;
  return 0;
}


// Compute the required functions on the data
int performComputations (INPUT_PARAMS * input)
{
  if(input->computeHessianBasedComp)
    {
      hessian_based_computation(input);
      cout <<"done!!" << endl;
    }
  return 0;
}

//Main function
int main (int argc, char* argv[])
{
  // read input parameters
  INPUT_PARAMS * input = new INPUT_PARAMS();

  int readParamRes = readParameters(argc, argv, input);
  if (readParamRes !=0)
    cout << "function readParameters error " << readParamRes <<endl;

  int computeRes = performComputations(input);
  cout <<" computeRes " << computeRes << endl;
  return 0;
}

/*
 int main()
 {
 typedef unsigned short PixelType;
 const unsigned int Dimension = 3;
 typedef itk::Image< PixelType, Dimension > ImageType;
 typedef itk::ImageFileReader< ImageType > ReaderType;
 ReaderType::Pointer reader = ReaderType::New();

 //reader->SetFileName( "/Users/arindambhattacharya/Research/internship/data/CFK-Prepreg-klein-570x356x374-1umVS-16bit.mhd" );

 reader->SetFileName( "/Users/arindambhattacharya/Research/internship/data/KUL.mhd" );


 typedef itk::ImageFileWriter< ImageType > WriterType;
 WriterType::Pointer writer = WriterType::New();
 writer->SetFileName( "./kul.nrrd" );
 writer->SetInput( reader->GetOutput() );
 try {
 writer->Update();
 }
 catch( itk::ExceptionObject & err )
 {
 std::cerr << "ExceptionObject caught !" << std::endl; std::cerr << err << std::endl;
 return EXIT_FAILURE;
 }

 std::cout << "Successful run !!" << std::endl;

 return 0;
 }

 */
// Software Guide : EndCodeSnippet

//  Software Guide : BeginLatex
//
//  This code instantiates a $3D$ image\footnote{Also known as a
//  \emph{volume}.} whose pixels are represented with type \code{unsigned
//  short}. The image is then constructed and assigned to a
//  \doxygen{SmartPointer}. Although later in the text we will discuss
//  \code{SmartPointer}'s in detail, for now think of it as a handle on an
//  instance of an object (see section \ref{sec:SmartPointers} for more
//  information).  The \doxygen{Image} class will be described in
//  Section~\ref{sec:ImageSection}.
//
//  Software Guide : EndLatex

