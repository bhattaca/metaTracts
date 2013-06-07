
#include <iostream>
//#include <vector>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "xBundle_hessian.h"


#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>
using namespace std;
using namespace itk;

typedef unsigned short PixelType;
const unsigned int Dimension =3;
typedef itk::Vector< double, Dimension > DoubleVectorType;


// Images
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageRegionIterator< ImageType > InIteratorType;

// Unsigned short
typedef itk::Matrix< PixelType, Dimension, Dimension > MatrixType;
typedef itk::Image< MatrixType, Dimension > MatrixImageType;
typedef MatrixImageType::RegionType     RegionType;
typedef MatrixImageType::SpacingType    SpacingType;
typedef MatrixImageType::PointType      PointType;


// Double precison
typedef itk::Matrix< double, Dimension, Dimension > DoubleMatrixType;
typedef itk::Image< DoubleMatrixType, Dimension > DoubleMatrixImageType;
typedef DoubleMatrixImageType::RegionType     DoubleRegionType;
typedef DoubleMatrixImageType::Pointer        DoubleMatrixImagePointer;
typedef itk::ImageRegionIterator< DoubleMatrixImageType > DoubleIteratorType;


typedef itk::RGBPixel< unsigned char > RGBPixelType;
typedef itk::Image< RGBPixelType, 3 >  RGBImageType;
typedef RGBImageType::Pointer          RGBImageTypePointer;
typedef itk::ImageRegionIterator< RGBImageType > RGBIteratorType;

// Hessian Filter
typedef HessianRecursiveGaussianImageFilter<ImageType> HessianFilterType;
typedef HessianFilterType::OutputImageType myHessianImageType;
typedef itk::ImageRegionIteratorWithIndex<HessianFilterType::OutputImageType>  HessianIteratorType;


typedef itk::SymmetricEigenAnalysis< DoubleMatrixType, DoubleVectorType, DoubleMatrixType > EigenVectorAnalysisType;



// Helper Functions
void setColor (DoubleMatrixType eigenMatrix, unsigned char * color)
{
  for (int i=0; i<Dimension; i++)
    {
    color[i]=((eigenMatrix[0][i]+1.0)/2.0)*255;
    }
}

// Compute Hessian / eigen analysis
void hessian_based_computation(const INPUT_PARAMS * input)
{
	cout <<"\t Hessian \n";
	//read an image
  ReaderType::Pointer reader = ReaderType::New();
  //reader->SetFileName( input->inputFileName.c_str());
  
  // reader->SetFileName( "/Users/arindambhattacharya/Research/internship/data/CFK-Prepreg-klein-570x356x374-1umVS-16bit.mhd");
  reader->SetFileName( "/Users/arindambhattacharya/Research/internship/data/Datensatz_simuliert.mhd");
  
  reader->Update();
  
  
  //debug
  /*
   ImageType::Pointer image = reader->GetOutput();
   
   ImageType::RegionType region = image->GetLargestPossibleRegion();
   
   ImageType::SizeType size = region.GetSize();
   
   std::cout << "size " << size << std::endl;
   */
  
  //hessian
  HessianFilterType::Pointer hessian = HessianFilterType::New();
  //float sigma = 3.0f;
  
  hessian->SetInput(reader->GetOutput());
  hessian->SetSigma(input->sigma);
  
  hessian->Update();
  
  cout <<"Hessian computation complete."<<endl;
  
  // Eigen Computation
  typedef itk::Matrix< double, Dimension> EigenValuesArrayType;
  typedef itk::Matrix< double, Dimension, Dimension > EigenVectorMatrixType;
  
  // Get info from the hessian image
  HessianFilterType::OutputImageType::Pointer HessianOutputImage = hessian->GetOutput();
  RegionType region = HessianOutputImage->GetLargestPossibleRegion();
  SpacingType spacing = HessianOutputImage->GetSpacing();
  PointType origin = HessianOutputImage->GetOrigin();
  
  // eigen vector /values place holders
  DoubleMatrixType eigenMatrix, myMatrix;
  eigenMatrix[0][0] = 0; eigenMatrix[0][1]=0;eigenMatrix[0][2]=0;
  eigenMatrix[1][0] = 0; eigenMatrix[1][1]=0;eigenMatrix[1][2]=0;
  eigenMatrix[2][0] = 0; eigenMatrix[2][1]=0;eigenMatrix[2][2]=0;
  
  DoubleVectorType eigenVal;
  eigenVal[0]=eigenVal[1]=eigenVal[2]=0;
  
  // Allocate the output image
  DoubleMatrixImagePointer m_Output = DoubleMatrixImageType::New();
  m_Output->SetRegions( region );
  m_Output->SetOrigin( origin );
  m_Output->SetSpacing( spacing );
  m_Output->Allocate();
  m_Output->FillBuffer( eigenMatrix );
  
  
  // test rgb image
  RGBImageTypePointer rgbEigen = RGBImageType::New();
  unsigned char color[3]={0.0, 0.0, 0.0};
  rgbEigen->SetRegions( region );
  rgbEigen->SetOrigin( origin );
  rgbEigen->SetSpacing( spacing );
  rgbEigen->Allocate();
  rgbEigen->FillBuffer( color );
  
  
  // set up iterators
  HessianIteratorType mIt( HessianOutputImage, region );
  DoubleIteratorType oIt( m_Output, region );
  RGBIteratorType rgbIt(rgbEigen, region);
  InIteratorType inIt(reader->GetOutput(),region);
  
  itk::SymmetricSecondRankTensor<double, 3> myTensor;
  mIt.GoToBegin();
  oIt.GoToBegin();
  rgbIt.GoToBegin();
  inIt.GoToBegin();
  
  
  EigenVectorAnalysisType eig;
  eig.SetDimension( Dimension );
  eig.SetOrderEigenMagnitudes( true );
  
  cout <<"All set for eigen computations " << endl;
  int ii=0;
  while( !mIt.IsAtEnd() )
    {
      eig.SetOrderEigenMagnitudes( true );
    // Compute eigen values and eigen matrix
    color[0]= 0;color[1]= 0;color[2]= 0;
    myTensor = mIt.Get();
    myMatrix[0][0] = myTensor[0];
    myMatrix[0][1] = myTensor[1];
    myMatrix[0][2] = myTensor[2];
    myMatrix[1][0] = myTensor[1];
    myMatrix[1][1] = myTensor[3];
    myMatrix[1][2] = myTensor[4];
    myMatrix[2][0] = myTensor[2];
    myMatrix[2][1] = myTensor[4];
    myMatrix[2][2] = myTensor[5];
    eig.ComputeEigenValuesAndVectors( myMatrix, eigenVal, eigenMatrix );
    
    
    //setColor (eigenMatrix, color);
    /*
     if (int (inIt.Get()) > input->thresh){
     for (int i=0; i<Dimension; i++)
     {
     //color[i]=((eigenMatrix[input->eigenVecNum][i]+1.0)/2.0)*255;
     //color[i]= abs(eigenMatrix[2][i])*255;
     //color[i]= int(inIt.Get()/65536.0f * 255);
     }
     }
     else{
     // debug setting everything to original color
     color[0]= int(inIt.Get()/65536.0f * 255);
     color[1]= int(inIt.Get()/65536.0f * 255);
     color[2]=  int(inIt.Get()/65536.0f * 255);
     //color[0]= 0;color[1]= 0;color[2]= 0;
     }
     */
    // compute from tubularness
    bool bug1 = false;
    bool bug2 = false;
    if ( eigenVal[1] < 0.0f  && eigenVal[2] < 0.0f && eigenVal[0] > 0.0f && int (inIt.Get()) > input->thresh )
      {
      bug1=true;
      if (abs(eigenVal[1])/(abs(eigenVal[2])*1.0f) < 0.9)
        {
        
        for (int i=0; i<Dimension; i++)
          {
          //color[i]=((eigenMatrix[input->eigenVecNum][i]+1.0)/2.0)*255;
          color[i]= abs(eigenMatrix[0][i])*255;
          //color[i]= int(inIt.Get()/65536.0f * 255);
          }
        bug2 = true;
        }
      }
    
    
    //    if (ii < 90){
    //      std::cout << " After ComputeEigenValuesAndVectors" << std::endl;
    //      cout <<" coords " << rgbIt.GetIndex() << endl;
    //      cout <<"scalar " << inIt.Get() <<"  th " << input->thresh<< endl;
    //      std::cout << "myTensor " << myTensor << std::endl;
    //      std::cout << myMatrix << std::endl;
    //      std::cout << "eigenVal " << eigenVal << std::endl;
    //      std::cout << "eigenMatrix \n" << eigenMatrix << std::endl;
    //      std::cout << "color " << int(color [0]) << " " << int(color [1]) <<" "<< int(color [2]) << std::endl;
    //      ii++;
    //    }
    
    if (input->degugMode )
      {
      ImageType::IndexType idx = mIt.GetIndex();
      if (idx[0]== 11 &&  idx[1]== 32 && idx[2]== 142 ){
        // set things
        //color[0]= 255;color[1]= 0;color[2]= 0;
        
        cout <<"idx " << idx;
        cout <<"pixel " << HessianOutputImage->GetPixel(idx) << endl;
        cout <<"tubularness  " <<abs(eigenVal[1])/(abs(eigenVal[2])*1.0f) << endl;
        cout <<"bug " << int (bug1) << " "<< int (bug2) << endl;
        std::cout << "myTensor " << myTensor << std::endl;
        cout <<"scalar " << inIt.Get() <<"  th " << input->thresh<< endl;
        std::cout << "color " << int(color [0]) << " " << int(color [1]) <<" "<< int(color [2]) << std::endl;
        std::cout << myMatrix << std::endl;
        std::cout << "eigenVal " << eigenVal << std::endl;
        std::cout << "eigenMatrix \n" << eigenMatrix << std::endl;
        
      }
      
      }
    
    rgbIt.Set(color); // set the color of the eigen matrix
    
    oIt.Set( eigenMatrix );
    ++oIt;
    ++mIt;
    ++rgbIt;
    ++inIt;
    }
  
  typedef itk::ImageFileWriter< RGBImageType > RGBWriterType;
  RGBWriterType::Pointer   writer = RGBWriterType::New();
  writer->SetInput( rgbEigen );
  writer->SetFileName( "EigenMatrix.mhd" );
  try
  {
  writer->Update();
  }
  catch ( itk::ExceptionObject e )
  {
  std::cerr << "Error: " << e << std::endl;
  }
  
}