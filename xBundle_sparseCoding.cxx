//
//  xBundle_sparseCoding.cxx
//  BundleEx
//
//  Created by arindam bhattacharya on 6/13/13.
//
//

#include "xBundle_sparseCoding.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkImageRegionIterator.h"

using namespace std;
using namespace itk;


typedef unsigned short PixelType;
const unsigned int Dimension=3;




// Images
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageRegionIterator< ImageType > InIteratorType;






void sparseCoding_computations (const INPUT_PARAMS * input)
{
  cout <<" In sparse coding" << endl;
  /*
   read the image
   iterate thorough it
   */
  
  // read image
  cout <<" fileName " << input->inputFileName<< endl;
  ReaderType::Pointer reader = ReaderType::New(); reader->SetFileName(input->inputFileName.c_str() );
  try
  { reader->Update(); }
  catch ( itk::ExceptionObject &err) {
    std::cout << "ExceptionObject caught !" << std::endl; std::cout << err << std::endl;
    exit(0);
  }
  
  //
  // neighborhood
  typedef itk::ConstantBoundaryCondition<ImageType>  BoundaryConditionType;
  typedef itk::ConstNeighborhoodIterator< ImageType, BoundaryConditionType > NeighborhoodIteratorType;
  
  NeighborhoodIteratorType::RadiusType radius; radius.Fill(input->radius);
  NeighborhoodIteratorType it( radius, reader->GetOutput(),
                              reader->GetOutput()->GetRequestedRegion() );
  //debug
  cout <<"size  " << reader->GetOutput()->GetRequestedRegion().GetSize();
  
  ofstream myfile (input->dataFileName.c_str());
  if (myfile.is_open())
    {
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
      ImageType::IndexType idx = it.GetIndex();
      
      for (int i=0; i< input->linInd;i++)
        {
        if (idx[0]== 0 &&  idx[1]== 0 && idx[2]== 0 )
          {
          cout <<i <<" - "<< it.GetIndex(i)<<" -- " << (int)it.GetPixel(i)<< endl;
          }
        myfile << (int)it.GetPixel(i);
        if (i < input->linInd-1)
          myfile <<",";
        }
      myfile <<"\n";
      }
    myfile.close();
    }
  else cout << "Unable to open file";
  
}

