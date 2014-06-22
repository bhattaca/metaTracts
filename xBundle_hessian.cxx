
#define _USE_MATH_DEFINES
#include <iostream>
#include <ctime>
#include <vector>
#include <queue>
#include <sstream>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <conio.h>

#include "itkPoint.h"
#include "itkindex.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "xBundle_hessian.h"
#include "itkRegionOfInterestImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include "itkBinaryMask3DMeshSource.h"
#include "itkMesh.h"
#include "itkVTKPolyDataWriter.h"


#include "itkBinaryMedianImageFilter.h"


#include "xBundle_fiber_info.h"
#include <numeric>
#include <math.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include "itkConstNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkImageRegionIterator.h"

#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include <vtkSmartPointer.h>
#include <vtkDiscreteMarchingCubes.h>
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"

#include <vtkMetaImageReader.h>
#include <vtkImageAccumulate.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkMaskFields.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtksys/ios/sstream>

// boost inlcudes
#include "boost/tokenizer.hpp"

//#include "itkSampleClassifier.h"
using namespace std;
using namespace itk;
using namespace boost;

//typedef unsigned short PixelType;
typedef float PixelType;
const unsigned int Dimension = 3;
typedef itk::Vector< double, Dimension > DoubleVectorType;


// Images
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::Image< float, Dimension > floatImageType;
typedef itk::Image< bool, Dimension > boolImageType;
typedef itk::Image< int, Dimension > intImageType;
typedef itk::Image< unsigned char, Dimension > ucharImageType;

typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileReader< floatImageType > floatReaderType;
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

//boolean 
typedef boolImageType::Pointer boolImagePointer;
typedef itk::ImageRegionIterator< boolImageType > boolIteratorType;
//int
typedef intImageType::Pointer intImageTypePointer;
typedef itk::ImageRegionIterator<intImageType> intIteratorType;

typedef itk::RGBPixel< unsigned char > RGBPixelType;
//typedef itk::RGBPixel< PixelType > RGBPixelType;
typedef itk::Image< RGBPixelType, 3 >  RGBImageType;
typedef RGBImageType::Pointer          RGBImageTypePointer;
typedef itk::ImageRegionIterator< RGBImageType > RGBIteratorType;

// Hessian Filter
typedef HessianRecursiveGaussianImageFilter<ImageType> HessianFilterType;
typedef HessianFilterType::OutputImageType myHessianImageType;
typedef itk::ImageRegionIteratorWithIndex<HessianFilterType::OutputImageType>  HessianIteratorType;


typedef itk::SymmetricEigenAnalysis< DoubleMatrixType, DoubleVectorType, DoubleMatrixType > EigenVectorAnalysisType;


// Function to correct the orientation of fibers in fiber bundles. 
void correctFiberBundleOrientaion(vector<FIBER> &bundle);
// Compute the average orientattion
void computeAverageOrientation(vector<FIBER> &bundle);

void averageOrientationEndPointsBundle(vector<FIBER> &bundle);
void correctOrientationEndPointsBundle(vector<FIBER> &bundle);


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	return buf;
}

template<typename t>
t dist(vector<t> v1, vector<t> v2)
{
	t tempD = 0.0;
	for (int d = 0; d < 3; d++)
	{
		tempD = tempD + (v1[d] - v2[d])*(v1[d] - v2[d]);
	}
	return tempD;
}

inline float distX(vector <int> v1, ImageType::IndexType v2)
{
	float tempD = 0.0;
	for (int d = 0; d < 3; d++)
	{
		tempD = tempD + (v1[d] - v2[d])*(v1[d] - v2[d]);
	}

	return tempD;
}


int smallFiberNum = 0;
// Helper Functions
void setColor(DoubleMatrixType eigenMatrix, unsigned char * color)
{
	for (int i = 0; i < Dimension; i++)
	{
		color[i] = ((eigenMatrix[0][i] + 1.0) / 2.0) * 255;
	}
}

// find near by points to tempIndex
void findNearPts(
	vector <int> tempIndex,
	const int spacing,
	const vector<int> radius,
	vector <vector<int> > & nearPts)
{
	int i = -radius[2];
	do{
		int j = -radius[1];
		do{
			int k = -radius[0];
			do{
				//
				vector<int> nearPtTemp(3, 0);
				nearPtTemp[0] = tempIndex[0] + spacing*k;
				nearPtTemp[1] = tempIndex[1] + spacing*j;
				nearPtTemp[2] = tempIndex[2] + spacing*i;
				//print_vec(nearPtTemp);
				nearPts.push_back(nearPtTemp);
				//
				k++;
			} while (k <= radius[0]);
			j++;
		} while (j <= radius[1]);
		i++;
	} while (i <= radius[2]);
}




//  print information
void print_info(const DoubleVectorType eigenVal, const DoubleMatrixType eigenMatrix, const  InIteratorType inIt,
	const ImageType::IndexType idx, const unsigned char *color, const INPUT_PARAMS * input, double vness);
//  print information
void print_info_small(const DoubleVectorType eigenVal, const DoubleMatrixType eigenMatrix, const  InIteratorType inIt,
	const ImageType::IndexType idx, const unsigned char *color, const INPUT_PARAMS * input, double vness);




// help functions for computeVesselness
// compute RB
void computeRB(const DoubleVectorType &eigenVal, const bool debug, double & RB)
{
	RB = abs(eigenVal[0]) / sqrt(abs(eigenVal[1] * eigenVal[2]));// deviation from  blob
	if (debug)
	{
		cout << "RB(deviation from blob)=" << RB << endl;
	}
}
// compute RA
void computeRA(const DoubleVectorType &eigenVal, const bool debug, double & RA)
{
	RA = abs(eigenVal[1] / eigenVal[2]); // plate vs line feature (small is good)
	if (debug)
	{
		cout << "RA(plate vs line)=" << RA << endl;
	}
}

// function to compute the second orderness 
void computeSecOrderness(const DoubleVectorType &eigenVal, const bool debug, double &secOness)
{
	double temp = 0.0;
	for (int i = 0; i < 3; i++)
	{
		temp = temp + eigenVal[i] * eigenVal[i];
	}
	secOness = sqrt(temp);

	if (debug)
	{
		cout << "sec orderness=" << secOness << endl;
	}

}
// function to compute the eigen vals and eigen vectors 
void  computeVesselness(const DoubleVectorType &eigenVal, const DoubleMatrixType &eigenMatrix, const bool debug, const float secT, double &vness)
{

	if (eigenVal[1] < 0.0 || eigenVal[2] < 0.0)
	{
		double RA = 0.0, RB = 0.0;
		computeRA(eigenVal, debug, RA);
		computeRB(eigenVal, debug, RB);
		//double secOrderness = 0.0;
		//computeSecOrderness(eigenVal, debug, secOrderness);
		vness = (1.0 - exp(-(RA*RA) / (2 * 0.25)))*(exp(-(RB*RB) / (2 * 0.25)));//*(1.0-exp(-(secOrderness*secOrderness)/(2*250.0)));

		/*		if (false)
		{
		cout <<"RA "<<RA<<" RB "<<RB <<" sec "<< secOrderness <<": ";
		cout <<"vness="<<vness<<",RA="<<(1.0-exp(-(RA*RA)/(0.25)))
		<<",RB="<<(exp(-(RB*RB)/(0.25)))
		<<",s="<<(1.0-exp(-(secOrderness*secOrderness)/(500.0)))<< endl;
		}*/
	}
}

// write to the files
void writeToFile(ofstream &f1, ofstream &f2, const DoubleMatrixType &eigenMatrix,
	const float v, const ImageType::IndexType idx)
{
	if (f1.is_open() && f2.is_open())
	{
		f1 << eigenMatrix[0][0] << "," << eigenMatrix[0][1] << "," << eigenMatrix[0][2] << "\n";
		f2 << idx[0] << "," << idx[1] << "," << idx[2] << endl;
	}
}
//compute the angle
/*
void computeAngle(DoubleMatrixImageType::IndexType idx, DoubleMatrixImageType::IndexType tempIdx,
DoubleMatrixImagePointer &m_Output, const INPUT_PARAMS * input,  int & A )
{
DoubleMatrixType m1,m2;
m1=m_Output->GetPixel(idx);
m2=m_Output->GetPixel(tempIdx);

float dotPdt = 0.0f;
for (int i=0;i<3;i++)
{
dotPdt=dotPdt+(m1[0][i]*m2[0][i]);
}
// if the product is negative then invert it
if (dotPdt < 0.0)
{
dotPdt = 0.0f;
for (int i=0;i<3;i++)
{
dotPdt=dotPdt+(m1[0][i]*m2[0][i]*(-1.0));
}
}

A=acos (dotPdt) * 180.0 / M_PI;

//cout <<"\n["<<m1[0][0]<<" , "<<m1[0][1]<<" , "<<m1[0][2]<<"]\n["
//	<<m2[0][0]<<" , "<<m2[0][1]<<" , "<<m2[0][2]<<"] : "<< A<<endl;

}
*/

//This one takes the m value to take care of the flipping 
//compute the angle
void computeAngle(const DoubleMatrixType &m1, DoubleMatrixImageType::IndexType tempIdx,
	DoubleMatrixImagePointer &m_Output, const INPUT_PARAMS * input, int & A)
{
	DoubleMatrixType m2;

	m2 = m_Output->GetPixel(tempIdx);

	float dotPdt = 0.0f;
	for (int i = 0; i < 3; i++)
	{
		dotPdt = dotPdt + (m1[0][i] * m2[0][i]);
	}
	// if the product is negative then invert it 
	if (dotPdt < 0.0)
	{
		A = 180 - acos(dotPdt) * 180.0 / M_PI;
	}
	else
	{
		A =  acos(dotPdt) * 180.0 / M_PI;
	}
	/*
	if (dotPdt < 0.0)
	{
		dotPdt = 0.0f;
		for (int i = 0; i < 3; i++)
		{
			dotPdt = dotPdt + (m1[0][i] * m2[0][i] * (-1.0));
		}
	}

	A = acos(dotPdt) * 180.0 / M_PI;
	*/
	//cout <<"m1 ["<<m1[0][0]<<" , "<<m1[0][1]<<" , "<<m1[0][2]<<"] m2 ["
	//	<<m2[0][0]<<" , "<<m2[0][1]<<" , "<<m2[0][2]<<"] : "<< " dotpdt "<< dotPdt <<" angle " << A<<endl;

}
// compute the perpendicular distance 
void findPerpendicularDistanceToPlane(const DoubleMatrixType &m, DoubleMatrixImageType::IndexType idx, DoubleMatrixImageType::IndexType tempIdx, float & d)
{
	float vec[3] = { 0.0f, 0.0f, 0.0f };
	// find the first vector
	for (int i = 0; i < 3; i++)
		vec[i] = tempIdx[i] - (idx[i] + m[1][i]);
	// find the dot pdt
	for (int j = 0; j < 3; j++)
		d = d + (m[0][j] * vec[j]);
	//d=abs(d);
}

// compute the horrizontal distance 
void findHorrDistanceToPlane(const DoubleMatrixType &m, DoubleMatrixImageType::IndexType idx, DoubleMatrixImageType::IndexType tempIdx, float & d)
{
	float vec[3] = { 0.0f, 0.0f, 0.0f };
	// find the first vector
	for (int i = 0; i < 3; i++)
		vec[i] = tempIdx[i] - (idx[i] + m[0][i]);
	float d1 = 0.0, d2 = 0.0;
	// find the dot pdt
	for (int j = 0; j < 3; j++)
	{
		d1 = d1 + (m[1][j] * vec[j]);
		d2 = d2 + (m[2][j] * vec[j]);
	}
	d1 = abs(d1); d2 = abs(d2);
	if (d1 < d2)
		d = d1;
	else
		d = d2;
}


void findHorrDistance(const DoubleMatrixType &m, DoubleMatrixImageType::IndexType idx,
	DoubleMatrixImageType::IndexType tempIdx, const float & perpD, float & d)
{
	float newPt[3] = { 0.0 };
	for (int i = 0; i < 3; i++)
	{
		newPt[i] = tempIdx[i] - (perpD*m[0][i]);
	}
	float dtmp = 0.0;
	for (int j = 0; j < 3; j++)
	{
		dtmp = dtmp + (newPt[j] - idx[j])*(newPt[j] - idx[j]);
	}
	d = sqrt(dtmp);
}
//Track  point 
class compare{
public:
	bool operator () (pair<float, vector<int> > &a, pair<float, vector<int> > &b)
	{
		/*	if (a.first > b.first)
		return true;
		else
		return false;*/
		if (a.first < b.first)
			return true;
		else
			return false;

	}
};
// push a point into the pqueue
void addPt2pq(priority_queue< pair<float, vector<int> >, vector<pair<float, vector<int> > >, compare > &pq,
	const ImageType::IndexType &idx)
{
	vector <int> vecIdx(3, 0);
	vecIdx[0] = idx[0];
	vecIdx[1] = idx[1];
	vecIdx[2] = idx[2];
	pair<float, vector<int> > a = make_pair(0.0, vecIdx);
	pq.push(a);
}



// Set the start points 
void setStartPointsFromConfig(
	const vector<int> startLocations,
	vector < vector <int> > & configStartPts
	)
{
	int numOfPoints = startLocations.size() / 3;
	for (int i = 0; i < numOfPoints; i++)
	{
		vector <int> temp (3,0); 
		temp[0] = startLocations[3 * i + 0];
		temp[1] = startLocations[3 * i + 1];
		temp[2] = startLocations[3 * i + 2];
		int spacings = 5;
		vector <int> radius(3, 0);
		radius[0] = 10;
		radius[1] = 10;
		radius[2] = 10;
		findNearPts(temp, spacings, radius, configStartPts);
	}
}

//Function to compute the average orientation at each location.
void computeAverageOrientation( vector<FIBER> & bundle)
{
	int bundleSise = bundle.size();
	for (int n = 0; n < bundleSise; n++)
	{
		const int fiberSize = bundle[n].points.size();
		if (fiberSize < 2)
			continue;
		//temporary vector to store the average orientaions
		vector< vector<float> > tempAvgOrientation(fiberSize, vector<float>(3,0.0));

		//Orientation at first point is the average of the first and second point. 
		vector<float> tempDir(3,0.0);
		for (int i = 0; i < 2; i++)
		{
			for (int d = 0; d < 3; d++)
			{
				tempDir[d] = tempDir[d] + bundle[n].dir[i][d];
			}
		}
		//set the average orientation of the first point in the fiber
		for (int d = 0; d < 3; d++)
		{
			//bundle[n].dir[0][d] = tempDir[d] / 2.0;
			tempAvgOrientation[0][d] = tempDir[d] / 2.0;
		}


		//Orientation of the last point is the average of the last and the second last point
		tempDir.clear();
		tempDir.resize(3, 0.0);
		for (int i = 0; i < 2; i++)
		{
			for (int d = 0; d < 3; d++)
			{
				tempDir[d] = tempDir[d] + bundle[n].dir[fiberSize - i - 1][d];
			}
	
		}
		for (int d = 0; d < 3; d++)
		{
			//bundle[n].dir[fiberSize - 1][d] = tempDir[d] / 2.0;
			tempAvgOrientation[fiberSize - 1][d] = tempDir[d] / 2.0;
		}
		

		//
		for (int p = 1; p < fiberSize - 1; p++)
		{
			vector<float> tempDir(3, 0.0);
			for (int j = -1; j < 2; j++)
			{
				for (int d = 0; d < 3; d++)
				{
					tempDir[d] = tempDir[d] + bundle[n].dir[p + j][d];
				}
			}
			for (int d = 0; d < 3; d++)
			{
				//bundle[n].dir[p][d] = tempDir[d] / 3.0;
				tempAvgOrientation[p][d] = tempDir[d] / 3.0;
			}
		}

		// set the orienataion by copying from tempAvgOrientaion

		for (int i = 0; i < fiberSize; i++)
		{
			for (int d = 0; d < 3; d++)
			{
				bundle[n].dir[i][d] = tempAvgOrientation[i][d];
			}
		}
	}
}

// Function to correct the orientation of fibers in fiber bundles. 
void correctFiberBundleOrientaion(vector<FIBER> &bundle)
{
	int bundleSise = bundle.size();
	for (int n = 0; n < bundleSise; n++)
	{
		int indexOfPointBeforeEndpt1;
		int indexOfPointBeforeEndpt2;
		int fiberSise = bundle[n].points.size();

		if (bundle[n].endPt1Index != 0){
			indexOfPointBeforeEndpt1 = bundle[n].endPt1Index - 1;
			//indexOfPointBeforeEndpt2 = bundle[n].endPt2Index - 1;

			float vecA[3] = { 0.0f, 0.0f, 0.0f };
			// find the first vector
			for (int i = 0; i < 3; i++)
				vecA[i] = bundle[n].points[indexOfPointBeforeEndpt1][i]
				- bundle[n].points[bundle[n].endPt1Index][i];

			float dotPdt = 0.0;
			for (int i = 0; i < 3; i++)
			{
				dotPdt = vecA[i] * bundle[n].dir[bundle[n].endPt1Index][i];
			}

			//DEBUG
			if (bundle[n].id == 219 || bundle[n].id == 235)
			{
				cout << "vecA " << vecA[0] << " " << vecA[1] << " " << vecA[2] << endl;
				cout << "dir at endpt1 " << bundle[n].dir[bundle[n].endPt1Index][0] << " " << bundle[n].dir[bundle[n].endPt1Index][1]
					<< " " << bundle[n].dir[bundle[n].endPt1Index][2] << endl;
				cout << "dotpdt " << dotPdt << endl;
			}
			//cout << "dotPdt " << dotPdt << endl;
			if (dotPdt < 0.0)
			{
				for (int i = 0; i <= bundle[n].endPt1Index; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						bundle[n].dir[i][j] = bundle[n].dir[i][j] * (-1.0);
					}
				}
			}
			else
			{
				// change the second half
				// change the second half 
				int start = bundle[n].endPt1Index + 1;

				for (int i = start; i < fiberSise; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						bundle[n].dir[i][j] = bundle[n].dir[i][j] * (-1.0);
					}
				}
			}
		}
		else
		{

			cout << "bundle[n].endPt1Index " << bundle[n].endPt1Index << endl;
			cout << "bundle[n].endPt2Index " << bundle[n].endPt2Index << endl;
			cout << "tracknum " << bundle[n].id << endl;
			// endPT1 is 0
			indexOfPointBeforeEndpt2 = bundle[n].endPt2Index - 1;
			if (indexOfPointBeforeEndpt2 < 0)
			{
				cerr << "ERROR in function reverse direction.";
				exit(0);
			}
			float vecA[3] = { 0.0f, 0.0f, 0.0f };
			// find the first vector
			for (int i = 0; i < 3; i++)
				vecA[i] = bundle[n].points[indexOfPointBeforeEndpt2][i]
				- bundle[n].points[bundle[n].endPt2Index][i];

			float dotPdt = 0.0;
			for (int i = 0; i < 3; i++)
			{
				dotPdt = vecA[i] * bundle[n].dir[bundle[n].endPt2Index][i];
			}
			cout << "*** DotPdt " << dotPdt << endl;
			if (dotPdt < 0.0)
			{
				for (int i = 1; i < fiberSise; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						bundle[n].dir[i][j] = bundle[n].dir[i][j] * (-1.0);
					}
				}
			}
			else
			{
				//change just the first one 
				for (int j = 0; j < 3; j++)
				{
					bundle[n].dir[0][j] = bundle[n].dir[0][j] * (-1.0);
				}
			}
		}
	}
}

//color a particular track 
void colorTrack(RGBImageTypePointer &rgbEigen, vector<FIBER> &bundle, const int id)
{
	//cout <<"In color track "<< id<<" bundle id "<<bundle[id].id <<endl;
	int m = id;

	for (int i = 0; i < bundle[id].allPoints.size(); i++)
	{
		ImageType::IndexType idx;
		for (int d = 0; d < 3; d++)
		{
			idx[d] = bundle[id].allPoints[i][d];
		}

		unsigned char color[3] = { 0, 0, 0 };
		color[0] = abs(bundle[id].dir[0][0]) * 255;
		color[2] = abs(bundle[id].dir[0][1]) * 255;
		color[1] = abs(bundle[id].dir[0][2]) * 255;
		//if (color[1] = abs(bundle[id].dir[0][0]) < 0.0)
		//{
		//	color[0] = abs(bundle[id].dir[bundle[id].endPt1Index][2]) * 255;
		//}
		//else
		//{
		//	color[2] = abs(bundle[id].dir[bundle[id].endPt1Index][2]) * 255;
		//}
		

		rgbEigen->SetPixel(idx, color);
	}
}


//color a particular track 
// color tracks according to cluster
void colorTrackID(intImageTypePointer &rgbEigen, vector<FIBER> &bundle, const int id)
{
	//cout <<"In color track "<< id<<" bundle id "<<bundle[id].id <<endl;
	int m = id;

	for (int i = 0; i < bundle[id].allPoints.size(); i++)
	{
		ImageType::IndexType idx;
		for (int d = 0; d < 3; d++)
		{
			idx[d] = bundle[id].allPoints[i][d];
		}

		rgbEigen->SetPixel(idx, bundle[id].id);
	}
}

//color a particular track
//color IDS
void colorTrackIDCluster(intImageTypePointer &rgbEigen,
	vector<FIBER> &bundle,
	const int id,
	const int color)
{
	//cout <<"In color track "<< id<<" bundle id "<<bundle[id].id <<endl;
	int m = id;

	for (int i = 0; i < bundle[id].allPoints.size(); i++)
	{
		ImageType::IndexType idx;
		for (int d = 0; d < 3; d++)
		{
			idx[d] = bundle[id].allPoints[i][d];
		}

		rgbEigen->SetPixel(idx, color);
	}
}


void colorTrack(RGBImageTypePointer &rgbEigen,
	vector<FIBER> &bundle,
	const int clusterInd,
	const int id)
{
	
	int m = id;

	unsigned char color[3] = { 0, 0, 0 };
	// the settings of the direction is not accurate 
	if (clusterInd == 1)
	{

		color[0] = 24;
		color[1] = 106;
		color[2] = 82;
	}
	else if (clusterInd == 2)
	{
		//color[0] = 255; color[1] = 128;
		color[0] = 252;
		color[1] = 70;
		color[2] = 12;
	}
	else if (clusterInd == 3)
	{
		//color[0]=255;color[1]=255;
		color[0] = 131;
		color[1] = 67;
		color[2] = 187;

	}
	else if (clusterInd == 4)
	{
		//color[0]=128;
		//color[1]=255;
		//color[0] = 61;
		//color[1] = 2;
		//color[2] = 25;
		color[2] = 255;
	}
	else if (clusterInd == 5)
	{
		//color[1]=255;
		color[0] = 252;
		color[1] = 230;
		color[2] = 146;

	}
	else if (clusterInd == 6)
	{

		//color[1]=255;
		//color[2]=128;
		color[0] = 64;
		color[1] = 191;
		color[2] = 222;
	}
	else if (clusterInd == 7)
	{
		//color[1]=255;
		//color[2]=255;
		color[0] = 128;
		color[1] = 100;
		color[2] = 10;
	}
	else if (clusterInd == 8)
	{
		//color[1]=128;
		//color[2]=255;
		color[0] = 140;
		color[1] = 81;
		color[2] = 10;
	}

	else if (clusterInd == 9)
	{
		color[0] = 216;
		color[1] = 179;
		color[2] = 101;
	}
	else if (clusterInd == 10)
	{
		color[0] = 246;
		color[1] = 232;
		color[2] = 195;
	}
	else if (clusterInd == 12)
	{
		color[0] = 199;
		color[1] = 234;
		color[2] = 229;
	}
	else if (clusterInd == 13)
	{
		color[0] = 90;
		color[1] = 180;
		color[2] = 172;
	}
	else if (clusterInd == 14)
	{
		color[0] = 1;
		color[1] = 102;
		color[2] = 94;
	}
	else if (clusterInd == 11)
	{
		color[1] = 245;
		color[0] = 127;
		color[2] = 127;
	}
	else
	{
		// if non of the clusters match return
		//cerr << "ERROR " << id << endl;
		return;
	}

	cout << "In color track "<<id <<" -- " << bundle[id].id << " cluster " << clusterInd << " color "
		<< int(color[0]) << " " << int(color[1]) << " " << int(color[2])
		<< " total size  " << bundle[id].allPoints.size() << endl;
	for (int i = 0; i < bundle[id].allPoints.size(); i++)
	{
		ImageType::IndexType idx;
		for (int d = 0; d < 3; d++)
		{
			idx[d] = bundle[id].allPoints[i][d];
		}
		rgbEigen->SetPixel(idx, color);
	}
}
typedef itk::RGBPixel< unsigned char > RGBPixelIntType;
//typedef itk::RGBPixel< PixelType > RGBPixelType;
typedef itk::Image< RGBPixelIntType, 3 >  RGBPixelIntImageType;
typedef RGBPixelIntImageType::Pointer          RGBIntImageTypePointer;
typedef itk::ImageRegionIterator< RGBPixelIntImageType > RGBIntIteratorType;

// function to  voxelize space
// rgbEigen  : each set of cluster fibers have the same integer color value
// bundle : fiber cluster information
// possible seed locations
// anisoPnt: result of anisotropic diffusion 
void voxelizeSpace(intImageTypePointer &rgbEigen,
	vector<FIBER> &bundle,
	const boolImagePointer &possibleSeeds,
	const ImageType::Pointer &anisoPnt)
{

	cout << "in voxelize space" << endl;
	// need the actual data
	intIteratorType it(rgbEigen, rgbEigen->GetLargestPossibleRegion());
	it.GoToBegin();
	int spacing = 2;
	vector <int> radius(3, 2);
	vector <vector<int> >nearPts;
	// list of indices and the corresponding colors
	vector <pair <intImageType::IndexType, int > > newColors;
	while (!it.IsAtEnd())
	{
		int color = it.Get();
		intImageType::IndexType idx;
		idx = it.GetIndex();
		//cout << idx[0]<<" "<<idx[1]<<" "<<idx[2]<<endl;

		if (color == 0) // add actual data
		{
			vector <vector<int> >nearPts;
			vector<int> vidx(3, 0);
			intImageType::IndexType idx;
			idx = it.GetIndex();
			if (anisoPnt->GetPixel(idx) > 0.0){
				//if (possibleSeeds->GetPixel(idx)){
				// not in cylinder but must be
				//cout <<"must give a color "<<endl;
				//cout << idx[0]<<" "<<idx[1]<<" "<<idx[2]<<endl;
				vidx[0] = idx[0];
				vidx[1] = idx[1];
				vidx[2] = idx[2];
				nearPts.clear();
				findNearPts(vidx, spacing, radius, nearPts);

				vector <pair<int, int> > hist;
				hist.clear();
				for (int i = 0; i<nearPts.size(); i++){
					idx[0] = nearPts[i][0];
					idx[1] = nearPts[i][1];
					idx[2] = nearPts[i][2];
					// if neighbor is in side the volume and clustered
					if (rgbEigen->GetLargestPossibleRegion().IsInside(idx) &&
						rgbEigen->GetPixel(idx) > 0)
					{
						//check if the id has been seen before
						bool found = false;
						for (int k = 0; k < hist.size(); k++)
						{
							if (hist[k].first == rgbEigen->GetPixel(idx))
							{
								hist[k].second++;
								//cout <<" increasing "<< hist[k].first << endl;
								found = true;
								break;
							}
						}
						if (found == false)
						{
							//cout <<" adding " << rgbEigen->GetPixel(idx) << endl;
							pair <int, int> foo;
							foo = make_pair(rgbEigen->GetPixel(idx), 0);
							hist.push_back(foo);
						}
					}
				}
				//if (!hist.size()==0)// old code
				// hist size must be above a min 
				if (hist.size() > 0)
				{
					int maxVal = -1;
					int maxClusID = 0;
					for (int k = 0; k<hist.size(); k++)
					{
						if (hist[k].second > maxVal)
						{
							maxClusID = hist[k].first;
						}
					}
					idx = it.GetIndex();
					pair<intImageType::IndexType, int> foo;
					foo = make_pair(idx, maxClusID);
					newColors.push_back(foo);
					//rgbEigen->SetPixel(idx, maxClusID);
					//cout <<"color set to "<< maxClusID <<endl;
					//set the surrounding to maxClusID
					/*
					idx = it.GetIndex();
					vidx[0]=idx[0];
					vidx[1]=idx[1];
					vidx[2]=idx[2];
					spacing = 1;
					vector <int> radius1 (3,3);
					nearPts.clear();
					findNearPts(vidx, spacing, radius1, nearPts);
					for (int k=0;k<nearPts.size();k++)
					{
					idx[0]=nearPts[k][0];
					idx[1]=nearPts[k][1];
					idx[2]=nearPts[k][2];
					// if neighbor is in side the volume and clustered
					if (rgbEigen->GetLargestPossibleRegion().IsInside(idx) &&
					rgbEigen->GetPixel(idx) == 0 &&
					anisoPnt->GetPixel(idx) > 20.0)
					{
					rgbEigen->SetPixel(idx, maxClusID);

					}
					}
					*/

				}//hist 
			}
		}
		++it;
	}
	// set all the pixels with new color
	for (int i = 0; i < newColors.size(); i++)
	{
		rgbEigen->SetPixel(newColors[i].first, newColors[i].second);
	}
}
// color the entire volume 
void colorVolume(RGBImageTypePointer &rgbEigen, RGBIntImageTypePointer &rgbVol,
	const boolImagePointer &possibleSeeds, const ImageType::Pointer &anisoPnt)
{
	cerr << "in color volume " << endl;
	RGBIteratorType rgbIt(rgbEigen, rgbEigen->GetLargestPossibleRegion());
	RGBIntIteratorType rgbVolIt(rgbVol, rgbVol->GetLargestPossibleRegion());
	boolIteratorType seedIt(possibleSeeds, rgbEigen->GetLargestPossibleRegion());
	rgbIt.GoToBegin();
	rgbVolIt.GoToBegin();
	seedIt.GoToBegin();
	RGBPixelType color;

	int spacings = 2;
	vector <int> radius(3, 0);
	radius[0] = 3;
	radius[1] = 3;
	radius[2] = 3;
	vector<int> vidx(3, 0);
	vector <vector<int> >nearPts;
	while (!rgbIt.IsAtEnd())
	{
		color = rgbIt.Get();
		rgbVolIt.Set(color);
		//cout << int(color.GetRed()) << int (color.GetGreen()) << int (color.GetBlue()) << endl;
		ImageType::IndexType idx;
		idx = seedIt.GetIndex();
		//if (seedIt.Get() && (int(color.GetRed())==0) && (int(color.GetBlue())==0) && (int(color.GetGreen())==0))
		if ((anisoPnt->GetPixel(idx) > 0.0) && (int(color.GetRed()) == 0) && (int(color.GetBlue()) == 0) && (int(color.GetGreen()) == 0))
		{

			vector <vector<int> >nearPts;
			vector<int> vidx(3, 0);
			//ImageType::IndexType idx;
			idx = seedIt.GetIndex();
			//cout <<"\t "<<idx[0]<<","<<idx[1]<<","<<idx[2]<<endl;
			vidx[0] = idx[0];
			vidx[1] = idx[1];
			vidx[2] = idx[2];
			findNearPts(vidx, spacings, radius, nearPts);
			int colorTemp[3] = { 0, 0, 0 };
			int j = 0;
			for (int i = 0; i<nearPts.size(); i++){
				idx[0] = nearPts[i][0];
				idx[1] = nearPts[i][1];
				idx[2] = nearPts[i][2];
				//cout <<i<<" out of "<<tot <<endl; 
				//cout <<"start " << idx << " value " << aniso_filter->GetOutput()->GetPixel(idx) << endl;
				if (rgbEigen->GetLargestPossibleRegion().IsInside(idx))
				{
					RGBPixelType  px = rgbEigen->GetPixel(idx);
					if ((int(px.GetRed()) > 0) || (int(px.GetGreen()) > 0) || (int(px.GetBlue()) > 0))
					{
						colorTemp[0] = colorTemp[0] + int(px.GetRed());
						colorTemp[1] = colorTemp[1] + int(px.GetGreen());
						colorTemp[2] = colorTemp[2] + int(px.GetBlue());
						j++;
						//cout <<"j "<<j <<endl;
					}
				}
			}
			if (j > 0){
				for (int k = 0; k < 3; k++){

					color[k] = (unsigned char)(int(colorTemp[k] / j*1.0));
				}
				//cout <<endl;
				rgbVolIt.Set(color);
			}

		}
		++rgbIt;
		++seedIt;
		++rgbVolIt;
	}
}

//Function to track individual fibers.
void trackPoint2Debug(DoubleMatrixImagePointer &m_Output,
	const INPUT_PARAMS * input,
	const ImageType::IndexType & trackId, // starting location of the fiber
	int &trackNum, //starts from 1.
	RGBImageTypePointer &rgbEigen,
	const boolImagePointer &possibleSeeds,
	vector<FIBER> &bundle)
{

	RegionType region = m_Output->GetLargestPossibleRegion();

	ImageType::SizeType size = region.GetSize();
	//set up output
	// write an image out 
	//ImageType::Pointer image = ImageType::New();
	//RegionType Outregion;
	///Outregion.SetSize(size);
	//image->SetRegions(Outregion);
	float dark = 0.0; float white = 255.0; float halfwhite = 127.0;
	unsigned char red[3] = { 255, 0, 0 };
	unsigned char blue[3] = { 0, 0, 255 };
	//image->Allocate();
	//image->FillBuffer(dark);
	//typedef  itk::ImageFileWriter< ImageType  > WriterType;
	//WriterType::Pointer writer = WriterType::New();
	//writer->SetFileName("out-color.mhd");
	//
	ImageType::IndexType idx, tempIdx;
	//setup the priority queue
	priority_queue< pair<float, vector<int> >, vector<pair<float, vector<int> > >, compare > pq;

	//addPt2pq(pq,idx);

	//set for finding near point


	vector <int> vecIdx(3, 0);
	vector <int> vecTmpIdx(3, 0);

	int spacings = 1;
	vector <int> radius(3, 0);
	radius[0] = 10;
	radius[1] = 10;
	radius[2] = 10;

	//while the priority queue is not empty
	int nm = 0;
	float OldDir[3] = { 0.0 };
	float CurrDir[3] = { 0.0 };
	//Run the program first in one direction and then the reverse

	//setup a fiber
	FIBER f;
	f.id = trackNum;
	float len = 0.0;
	bool reverseVector = false;

	for (int dir = 0; dir < 2; dir++){

		pq = priority_queue< pair<float, vector<int> >, vector<pair<float, vector<int> > >, compare >();
		bool firstRun = true; // first run in this  direction 
		//set idx to trackId
		for (int d = 0; d < 3; d++)
		{
			idx[d] = trackId[d];
		}

		addPt2pq(pq, idx);

		for (int d = 0; d < 3; d++)
		{
			OldDir[d] = 0.0;
		}
		//DEBUG
		//cout << "idx " << idx[0] << " " << idx[1] << " " << idx[2] << endl;
		//cout <<"starting in a new direction \n size of pq " << pq.size() <<endl;
		//cout  <<"tracking idx " << idx <<endl;
		bool flagExpand = true;
		//while ((!pq.empty()) ){
		//while ( pq.size() < 200  && !pq.empty() && flagExpand){


		while (flagExpand){

			//set Old dir

			//cout <<"---------------------------\nOld Dir ";
			//for (int d=0;d<3;d++)
			//{
			//	cout <<OldDir[d]<<" ";
			//}
			//cout <<"\n";

			//update vecIdx and idx
			vecIdx[0] = pq.top().second[0];
			vecIdx[1] = pq.top().second[1];
			vecIdx[2] = pq.top().second[2];


			//compute the length between vecIdx and idx
			float tempD = 0.0;
			for (int d = 0; d < 3; d++)
			{
				tempD = tempD + (vecIdx[d] - idx[d])*(vecIdx[d] - idx[d]);
			}
			len = len + sqrt(tempD);

			idx[0] = pq.top().second[0];
			idx[1] = pq.top().second[1];
			idx[2] = pq.top().second[2];

			float tempPerpDist = pq.top().first;


			pq.pop();

			vector <vector<int> >nearPts;
			findNearPts(vecIdx, spacings, radius, nearPts);
			// cout <<"size of near points "<< nearPts.size() <<endl;
			DoubleMatrixType m = m_Output->GetPixel(idx);

			// first one way then the other to get the entire thing.
			
			if (dir == 0 && firstRun){
				for (int k = 0; k < 3; k++)
				{
					m[0][k] = m[0][k] * -1.0;
				}

				firstRun = false;
			}
			

			//if true then the vector is dir is reversed in the first run, no need to do it for dir ==1
			//else for dir == 1 must reverse vector
			/*
			if (firstRun && dir == 0){
					if (m[0][2] < 0.0)
					{
						reverseVector = true; 
					}
				
				if (reverseVector){
					for (int d = 0; d < 3; d++)
					{
						m[0][d] = m[0][d] * -1.0;
					}
				}

				firstRun = false;
			}
			if (firstRun  && dir == 1){
				if (!reverseVector)
				{
					for (int d = 0; d < 3; d++)
					{
						m[0][d] = m[0][d] * -1.0;
					}
				}
				firstRun = false;
			}
			*/
			//cout << "dir " << dir << endl;
			//cout << "m " << m[0][0] << " " << m[0][1] << " " << m[0][2] << endl;

			f.pushPoint(vecIdx);
			f.pushPerpDist(tempPerpDist);
			vector <float> tempDir(3, 0.0);

			//Update m 
			float dotPdt = 0.0f;
			for (int i = 0; i < 3; i++)
			{
				dotPdt = dotPdt + (m[0][i] * OldDir[i]);
			}
			// if the product is negative then invert it 
			if (dotPdt < 0.0)
			{
				for (int i = 0; i < 3; i++)
				{
					m[0][i] = -1.0*m[0][i];
				}
			}
			//
			//DEBUG
			//cout << "OldDir ";
			for (int d = 0; d < 3; d++)
			{
				OldDir[d] = m[0][d];
				tempDir[d] = m[0][d];
				//cout << OldDir[d] << " ";
			}
			//cout << endl;

			f.pushDir(tempDir);
			int n = 0;// number of elements added

			//debug donot store
			pq = priority_queue< pair<float, vector<int> >, vector<pair<float, vector<int> > >, compare >();
			//DEBUG
			//cout << "nearPoints size " << nearPts.size() << endl;
			//cout <<" points being added "<< endl;
			//for all near points


			for (int i = 0; i < nearPts.size() && n < 400; i++)
			{
				//set temp Idx as the near point
				for (int j = 0; j < 3; j++)
				{
					tempIdx[j] = nearPts[i][j];
					vecTmpIdx[j] = nearPts[i][j];;
				}
				if (idx == tempIdx)
					continue;
				// if the point is inside the region
				if (region.IsInside(tempIdx))
				{
					float perpDist = 0.0;
					float horrDist = 0.0;

					findPerpendicularDistanceToPlane(m, idx, tempIdx, perpDist);

					findHorrDistance(m, idx, tempIdx, perpDist, horrDist);

					if (perpDist > 0.0 && perpDist < input->cylLength && horrDist < input->cylRadius)
						f.allPoints.push_back(vecTmpIdx);

					// inside the cylinder
					//if (perpDist > 0.5 && perpDist < 4.0 && horrDist < 2.0)
					if (perpDist > 4.5 && perpDist < input->cylLength && horrDist < input->cylRadius)
					{

						//if (dir==0){
						//	image->SetPixel(tempIdx, white);
						//	rgbEigen->SetPixel(tempIdx,red);
						//}
						//else{
						//	image->SetPixel(tempIdx, halfwhite);
						//	rgbEigen->SetPixel(tempIdx,blue);
						//}

						//unsigned char color[3]={trackNum, trackNum, trackNum};
						unsigned char color[3] = { 255, 255, 255 };
						//color[0]=abs(m[0][0])*255;
						//color[1]=abs(m[0][1])*255;
						//color[2]=abs(m[0][2])*255;
						//if (!input->colorTrack)
						//debug




						if (possibleSeeds->GetPixel(tempIdx) == false)
						{
							continue;
						}

						//compare angles
						int angle = 0;
						computeAngle(m, tempIdx, m_Output, input, angle);

						if (angle < 8.0){

							float combinedD = 0.0;
							//combinedD = 0.33*exp(-angle*angle / 5.0) + 0.66*(1 - exp(-perpDist*perpDist / 4.0));
							//DEBUG
							combinedD = 0.33*exp(-angle*angle / 5.0) + 0.66*(1 - exp(-perpDist*perpDist / 10.0));
							pair<float, vector<int> > a = make_pair((combinedD), vecTmpIdx);
							pq.push(a);
							//check if anything is added in this turn
							n++;
							//debug
							/*
							color[0]=combinedD*color[0];
							color[1]=combinedD*color[1];
							color[2]=combinedD*color[2];
							rgbEigen->SetPixel(tempIdx,color);
							*/
							const DoubleMatrixType m2 = m_Output->GetPixel(tempIdx);

							//cout <<"["<<nearPts[i][0]<<","<<nearPts[i][1]<<","<<nearPts[i][2]<<"] :";
							//cout <<"perp dis: "<< perpDist<<" horr "<< horrDist <<","
							//	<<" dist "<<perpDist <<" angle "<< angle;
							////nm++;
							//cout <<"] n "<<n<<endl;
							//image->SetPixel(tempIdx, white);
						}
					}
				}
			}//for nearpoints

			if (n == 0)
			{
				flagExpand = false;
				if (dir == 0){
					f.endPt1Index = f.points.size()-1;
				}
				else
				{
					f.endPt2Index = f.points.size()-1;
				}

				// DEBUG
				//cout << "endpt1 " << idx[0] << " " << idx[1] << " " << idx[2];
				//cout << " Length of fiber " << len << " dir " << dir << endl;
			}

		}//while flagExapand

	}//for 

	f.length = len;
	if (len > input->minFiberLength)
	{

		cout << "trackID " << trackId[0] << " " << trackId[1] << " " << trackId[2];
		cout << " track num " << trackNum << endl;

		bundle.push_back(f);
		trackNum++;
	}
	else
	{
		smallFiberNum++;
	}
}


/*
// track a point
void trackPoint ( DoubleMatrixImagePointer &m_Output,
const INPUT_PARAMS * input, RGBImageTypePointer &rgbEigen)
{
cerr <<"in function track point " << endl;

RegionType region = m_Output->GetLargestPossibleRegion();

ImageType::SizeType size = region.GetSize();
//set up output
// write an image out
ImageType::Pointer image = ImageType::New();
RegionType Outregion;
Outregion.SetSize(size);
image->SetRegions(Outregion);
float dark = 0.0; float white = 255.0;float halfwhite = 127.0;
unsigned char red[3]={255, 0, 0};
unsigned char blue[3]={0, 0, 255};
image->Allocate();
image->FillBuffer(dark);
typedef  itk::ImageFileWriter< ImageType  > WriterType;
WriterType::Pointer writer = WriterType::New();
writer->SetFileName("out-color.mhd");
//
ImageType::IndexType idx, tempIdx;
//idx[0]=11;
//idx[1]=29;
//idx[2]=99;//98;
//48,15,21
//49, 14, 13
//setup the priority queue
priority_queue< pair<float,vector<int> >, vector<pair<float,vector<int> > >, compare > pq;

//addPt2pq(pq,idx);

//set for finding near point
vector <int> vecIdx(3,0);
vector <int> vecTmpIdx(3,0);

int spacings=1;
vector <int> radius(3,0);
radius[0]=10;
radius[1]=10;
radius[2]=10;

//while the priority queue is not empty
int nm=0;
float OldDir[3]={0.0};
float CurrDir[3]={0.0};
//Run the program first in one direction and then the reverse

for (int dir=0;dir<2;dir++){
pq = priority_queue< pair<float,vector<int> >, vector<pair<float,vector<int> > >, compare >();
idx[0]=103;//11;
idx[1]=21;//29
idx[2]=78;//99;//98;
addPt2pq(pq,idx);

for (int d=0;d<3;d++)
{
OldDir[d]=0.0;
}

cerr <<"size of pq " << pq.size() <<endl;
cerr <<"idcx " << idx <<endl;
//while ((!pq.empty()) ){

int t=0;
//maintian a vector of selected points
vector<vector<int> > CandidatePts;

while ( pq.size() < 200  && !pq.empty()){

//set Old dir
cout <<"\n---------------------------\nOld Dir ";
for (int d=0;d<3;d++)
{
cout <<OldDir[d]<<" ";
}
cout <<"\n";
//update vecIdx and idx
vecIdx[0]=pq.top().second[0];
vecIdx[1]=pq.top().second[1];
vecIdx[2]=pq.top().second[2];

idx[0]=pq.top().second[0];
idx[1]=pq.top().second[1];
idx[2]=pq.top().second[2];
float tempPerpDist = pq.top().first;
//put the selected point in the list


//print the current selected point list
cout <<"selected points:\n";
for(int i=0;i<CandidatePts.size();i++)
{
cout <<CandidatePts[i][0]<<","<<CandidatePts[i][1]<<","<<CandidatePts[i][2]<<","<<endl;
}
pq.pop();

vector <vector<int> >nearPts;
findNearPts(vecIdx, spacings, radius, nearPts);
//cout <<"size of near points "<< nearPts.size() <<endl;
DoubleMatrixType m=m_Output->GetPixel(idx);

// first one way then the other to get the entire thing.
if (dir==0)
for (int k=0;k<3;k++)
{
m[0][k]=m[0][k]*-1.0;
}

cout <<"idx "<<idx<<" m "<<m[0][0]<<" , "<<m[0][1]<<" , "<<m[0][2]<<
" perpDist " <<tempPerpDist <<
" size of pq "<<pq.size()<<endl;

//Update m
float dotPdt = 0.0f;
for (int i=0;i<3;i++)
{
dotPdt=dotPdt+(m[0][i]*OldDir[i]);
}
// if the product is negative then invert it
if (dotPdt < 0.0)
{
for (int i=0;i<3;i++)
{
m[0][i]=-1.0*m[0][i];
}
}
//
cout <<"New Dir ";
for (int d=0;d<3;d++)
{
OldDir[d]=m[0][d];
cout <<OldDir[d]<<" ";
}

cout <<"   t " << t ;
cout <<"\n-------------------------\n";

int n=0;// number of elements added
int l=0;
//for all near points
for (int i=0;i<nearPts.size();i++)
{
//set temp Idx as the near point
for (int j=0;j<3;j++)
{
tempIdx[j]=nearPts[i][j];
vecTmpIdx[j]=nearPts[i][j];;
}
if (idx == tempIdx)
continue;
// if the point is inside the region
if(region.IsInside(tempIdx))
{
float perpDist = 0.0;
float horrDist = 0.0;

findPerpendicularDistanceToPlane(m, idx, tempIdx, perpDist);

findHorrDistance(m, idx, tempIdx, perpDist, horrDist);

// inside the cylinder
if (perpDist >0.0 && perpDist < 10.0 && horrDist < 3.0)
{

if (dir==0){
image->SetPixel(tempIdx, white);
rgbEigen->SetPixel(tempIdx,red);
}
else{
image->SetPixel(tempIdx, halfwhite);
rgbEigen->SetPixel(tempIdx,blue);
}

//compare angles
int angle=0;
computeAngle(m, tempIdx, m_Output, input, angle);



if (angle < 3.0  && l < 5){
//check if the point has been selected before
vector<vector<int> >::iterator it;
it = find (CandidatePts.begin(), CandidatePts.end(), vecTmpIdx);

if ((it==CandidatePts.end())){
pair<float,vector<int> > a = make_pair((perpDist), vecTmpIdx);
pq.push(a);
l++;
CandidatePts.push_back(vecTmpIdx);
//debug
const DoubleMatrixType m2=m_Output->GetPixel(tempIdx);
cout <<"["<<nearPts[i][0]<<","<<nearPts[i][1]<<","<<nearPts[i][2]<<"] :";
cout <<"perp dis: "<< perpDist<<" horr "<< horrDist <<","
<<" dist "<<perpDist <<" angle "<< angle;
nm++;
cout <<"] nm "<<nm<<endl;
//image->SetPixel(tempIdx, white);
n++;
}
else
{
cout <<"["<<nearPts[i][0]<<","<<nearPts[i][1]<<","<<nearPts[i][2]<<"] :";
cout <<"["<<vecTmpIdx[0]<<","<<vecTmpIdx[1]<<","<<vecTmpIdx[2]<<"] :";
cout <<" is already selected \n";
}
}
}
}
}//for

if (n==0){
cout <<"no elements added in this iteration \n";
t=t-10;;
}
else{
t=t+10;
}

}//while

}//for
writer->SetInput(image);
writer->Update();
}


*/
// compute reliable hessians
// donot use // outdated code
void computeReliableHessians(DoubleMatrixImagePointer &m_Output, const INPUT_PARAMS * input)
{

	// setup file output
	ofstream datafile("data.csv");
	ofstream indexFile("index.csv");
	datafile << "x,y,z\n";
	DoubleIteratorType oIt(m_Output, m_Output->GetLargestPossibleRegion());
	oIt.GoToBegin();
	DoubleMatrixImageType::RegionType region = m_Output->GetLargestPossibleRegion();
	DoubleMatrixImageType::SizeType size = region.GetSize();
	std::cout << size << std::endl;
	while (!oIt.IsAtEnd())
	{
		DoubleMatrixImageType::IndexType idx = oIt.GetIndex();
		DoubleMatrixImageType::IndexType tempIdx;
		//cout <<idx[0]<<","<<idx[1]<<","<<idx[2]<<endl;
		//settings for finding near points
		vector <int> vecIdx(3, 0);
		vecIdx[0] = idx[0];
		vecIdx[1] = idx[1];
		vecIdx[2] = idx[2];
		int spacings = 1;
		vector <int> radius(3, 0);
		radius[0] = 8;
		radius[1] = 8;
		radius[2] = 8;
		// check vnesss

		vector <vector<int> >nearPts;
		findNearPts(vecIdx, spacings, radius, nearPts);

		int tot = 0;
		int match = 0;
		const DoubleMatrixType m = oIt.Get();
		if (input->degugMode || 1)
		{
			cout << "point " << vecIdx[0] << "," << vecIdx[1] << "," << vecIdx[2] << endl;
			//cout <<"m "<<m <<endl;
		}
		for (int i = 0; i < nearPts.size(); i++)
		{
			//set temp Idx as the near point
			for (int j = 0; j < 3; j++)
			{
				tempIdx[j] = nearPts[i][j];
			}
			if (input->degugMode)
				cout << nearPts[i][0] << "," << nearPts[i][1] << "," << nearPts[i][2] << " is inside " << region.IsInside(tempIdx) << "\n";
			// if the point is inside the region
			if (region.IsInside(tempIdx))
			{
				float perpDist = 0.0;
				float horrDist = 0.0;

				findPerpendicularDistanceToPlane(m, idx, tempIdx, perpDist);

				if (input->degugMode)
					cout << "distance: " << m[0][0] << "," << m[0][1] << "," << m[0][2] << " : " << perpDist << endl;

				findHorrDistanceToPlane(m, idx, tempIdx, horrDist);

				//if(input->degugMode)
				//	cout <<"horr distance: "<<m[0][0]<<","<<m[0][1]<<","<<m[0][2]<<" : "<<m[1][0]<<","<<m[1][1]<<","<<m[1][2]<<" :"<< horrDist <<endl;

				//if(input->degugMode)
				//	cout <<"perp dis: "<< perpDist<<" horr "<< horrDist <<endl;

				// inside the cylinder
				if (perpDist < 8.0 && horrDist < 4.0)
				{
					tot++;
					//compare angles
					int angle = 0;

					//computeAngle(idx, tempIdx, m_Output, input, angle);
					if (angle < 30 || angle  > 150)
					{
						match++;
						if (input->degugMode)
							cout << " angle error " << endl;
					}
				}

			}
		}
		if (input->degugMode)
			cout << "tot " << tot << " match " << match << endl;
		if (input->degugMode)
			cout << "\n";

		if (tot > 0)
		if (match / (tot*1.0) > 0.6){
			writeToFile(datafile, indexFile, m, 0, idx);
			cout << m[0][0] << "," << m[0][1] << "," << m[0][2] << match / (tot*1.0) << endl;
		}
		oIt++;
	}

}


// by part
//ofstream bigdatafile ("bigdata.csv");
//ofstream bigindexFile ("bigindex.csv");
//ofstream numSelectFile ("numSelected.csv");
//typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
//
//void regionHessian(const INPUT_PARAMS * input, FilterType::Pointer &regionFilter, int &n);
//void bypart(const INPUT_PARAMS *input)
//{
//	cout <<"in by parts"<<endl;
//	//
//	//read an image
//	ReaderType::Pointer reader = ReaderType::New(); 
//
//	string part = ("K:\\Internship\\data\\crop-11\\crop-11.mhd");
//	// string prepreg = ("K:\\Internship\\data\\prepreg\\CFK-Prepreg-Cut1.mhd");
//	//string prepreg = ("K:\\Internship\\data\\prepreg\\crop-4\\crop-4.mhd");
//
//	reader->SetFileName(part);
//	reader->Update();
//	cerr <<"status: reading complete."<<endl;
//
//	ImageType::SizeType inSize = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
//	cout <<"size "<< inSize <<endl;
//
//	int boxSize = input->boxSize;
//	int numBox[3]={0};
//	for (int i=0;i<3;i++)
//	{
//		numBox[i] = inSize[i]/boxSize;
//	}
//	cout <<"number of boxes " << numBox[0]<<","<<numBox[1]<<","<<numBox[2]<<endl;
//	if (input->aniso_diff){
//		cerr <<"using the aniso diffusion before the hessian "<<endl;
//	}
//
//	bigdatafile << "x,y,z\n";
//	numSelectFile <<"n \n";
//	ImageType::IndexType start;
//	ImageType::SizeType size;
//	size[0]=boxSize;size[1]=boxSize;size[2]=boxSize;
//	//debug
//	int count =0;
//	for (int z=0;z<numBox[2];z++)
//		for(int y=0;y<numBox[1];y++)
//			for(int x=0;x<numBox[0];x++)
//			{
//				//debug
//				if (count < 3 || 1){
//					start[0]=x*boxSize;
//					start[1]=y*boxSize;
//					start[2]=z*boxSize;		
//					cout <<"start "<<start<<endl;
//
//					FilterType::Pointer regionFilter = FilterType::New();
//					ImageType::RegionType desiredRegion;
//					desiredRegion.SetSize(size);
//					desiredRegion.SetIndex(start);
//
//					regionFilter->SetRegionOfInterest(desiredRegion);
//					regionFilter->SetInput(reader->GetOutput());
//					try
//					{
//						regionFilter->Update();
//					}
//					catch ( itk::ExceptionObject e )
//					{
//						std::cerr << "Error: " << e << std::endl;
//					}
//					int n=0;
//					regionHessian(input, regionFilter,n);
//					numSelectFile << n <<endl;
//					count++;
//				}
//			}
//
//			/*
//			typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
//			FilterType::Pointer regionFilter = FilterType::New();
//
//			ImageType::RegionType desiredRegion;
//			desiredRegion.SetSize(size);
//			desiredRegion.SetIndex(start);
//
//			regionFilter->SetRegionOfInterest(desiredRegion);
//			regionFilter->SetInput(reader->GetOutput());
//			regionFilter->Update();
//			///write output
//			typedef itk::ImageFileWriter< ImageType > WriterType;
//			WriterType::Pointer   writer = WriterType::New();
//			writer->SetInput( regionFilter->GetOutput());
//			writer->SetFileName("test.mhd" );
//			try
//			{
//			writer->Update();
//			}
//			catch ( itk::ExceptionObject e )
//			{
//			std::cerr << "Error: " << e << std::endl;
//			}
//			*/
//
//}

//void regionHessian(const INPUT_PARAMS * input,
//				   FilterType::Pointer &regionFilter,
//				   int &n)
//{
//
//	//aniso diff
//	// put input into
//	typedef GradientAnisotropicDiffusionImageFilter< ImageType,ImageType> aniso_filter_type;
//	aniso_filter_type::Pointer aniso_filter = aniso_filter_type::New();
//	if (input->aniso_diff){
//		aniso_filter->SetInput(regionFilter->GetOutput());
//		int numberOfIterations = 5;
//		aniso_filter->SetNumberOfIterations( numberOfIterations );
//		aniso_filter->SetTimeStep( 0.0625 );
//		aniso_filter->SetConductanceParameter( 3.0 );
//		aniso_filter->Update();
//	}
//
//	//hessian
//	HessianFilterType::Pointer hessian = HessianFilterType::New();
//	if (input->aniso_diff){
//		hessian->SetInput(aniso_filter->GetOutput());
//		//cout <<"using the aniso diffusion before the hessian "<<endl;
//	}
//	else
//	{
//		//cout <<"NO aniso diffusion before the hessian "<<endl;
//		hessian->SetInput(regionFilter->GetOutput());
//	}
//
//
//	hessian->SetSigma(input->sigma);
//	hessian->Update();
//
//	//cerr <<"status: Hessian complete."<<endl;
//
//	// Eigen Computation
//	typedef itk::Matrix< double, Dimension> EigenValuesArrayType;
//	typedef itk::Matrix< double, Dimension, Dimension > EigenVectorMatrixType;
//
//	// Get info from the hessian image
//	HessianFilterType::OutputImageType::Pointer HessianOutputImage = hessian->GetOutput();
//	RegionType region = HessianOutputImage->GetLargestPossibleRegion();
//	SpacingType spacing = HessianOutputImage->GetSpacing();
//	PointType origin = HessianOutputImage->GetOrigin();
//
//	// eigen vector /values place holders
//	DoubleMatrixType eigenMatrix, myMatrix;
//	eigenMatrix[0][0] = 0; eigenMatrix[0][1]=0;eigenMatrix[0][2]=0;
//	eigenMatrix[1][0] = 0; eigenMatrix[1][1]=0;eigenMatrix[1][2]=0;
//	eigenMatrix[2][0] = 0; eigenMatrix[2][1]=0;eigenMatrix[2][2]=0;
//
//	DoubleVectorType eigenVal;
//	eigenVal[0]=eigenVal[1]=eigenVal[2]=0;
//
//	// set up iterators
//
//	InIteratorType inIt(regionFilter->GetOutput(),region);
//
//	itk::SymmetricSecondRankTensor<double, 3> myTensor;
//	HessianIteratorType mIt( HessianOutputImage, region );
//
//	mIt.GoToBegin();
//	inIt.GoToBegin();
//	EigenVectorAnalysisType eig;
//	eig.SetDimension( Dimension );
//	eig.SetOrderEigenMagnitudes( true );
//
//	//cerr <<"status: set for eigen computation." << endl;
//
//	// setup file output
//	//ofstream datafile ("data.csv");
//	//ofstream indexFile ("index.csv");
//	//ofstream vnessFile ("vness.csv");
//	double expThSecOrderness = 40.0;
//	unsigned char color[3]={0};
//	while( !inIt.IsAtEnd() )
//	{
//		eig.SetOrderEigenMagnitudes( true );
//		// Compute eigen values and eigen matrix
//		color[0]= 0;color[1]= 0;color[2]= 0;
//		myTensor = mIt.Get();
//
//		myMatrix[0][0] = myTensor[0];
//		myMatrix[0][1] = myTensor[1];
//		myMatrix[0][2] = myTensor[2];
//		myMatrix[1][0] = myTensor[1];
//		myMatrix[1][1] = myTensor[3];
//		myMatrix[1][2] = myTensor[4];
//		myMatrix[2][0] = myTensor[2];
//		myMatrix[2][1] = myTensor[4];
//		myMatrix[2][2] = myTensor[5];
//
//		// compute the eigen values and eigen vectors
//		eig.ComputeEigenValuesAndVectors( myMatrix, eigenVal, eigenMatrix );
//		// Compute vesseleness 
//		double vness=0.0;
//		computeVesselness(eigenVal, eigenMatrix, false, expThSecOrderness, vness);
//		// if the input is greater than the threhold.
//		ImageType::IndexType idx = inIt.GetIndex();
//
//		if (input->reliableHess == false)
//			if (vness >= input->vness){
//				writeToFile(bigdatafile, bigindexFile, eigenMatrix,vness,idx);
//				n++;
//			}
//			++mIt;
//			++inIt;
//	}
//}

void computeColor(const INPUT_PARAMS * input)
{
	string directions("angle.csv");
	cout << "reading the color " << endl;
	ifstream in(directions.c_str());
	if (!in.is_open()) { cerr << "error in reading,"; exit(0); };

	typedef tokenizer< escaped_list_separator<char> > Tokenizer;
	vector <float> directionsIndices;
	vector< string > vec;
	string line;
	while (getline(in, line))
	{
		Tokenizer tok(line);
		vec.assign(tok.begin(), tok.end());
		//copy(vec.begin(), vec.end(),
		//ostream_iterator<string>(cout, "|"));

		vector<string>::iterator it = vec.begin();
		istringstream os(*it);
		float d;
		os >> d;

		directionsIndices.push_back(d);
	}
	//
	//read an image
	ReaderType::Pointer reader = ReaderType::New();

	string part = ("K:\\Internship\\data\\crop-11\\crop-11.mhd");
	//string prepreg = ("K:\\Internship\\data\\prepreg\\CFK-Prepreg-Cut1.mhd");
	//string prepreg = ("K:\\Internship\\data\\prepreg\\crop-4\\crop-4.mhd");

	reader->SetFileName(part);
	reader->Update();
	cerr << "status: reading complete." << endl;

	ImageType::SizeType inSize = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
	cout << "size " << inSize << endl;

	int boxSize = input->boxSize;
	cout << "boxSize " << boxSize << endl;
	int numBox[3] = { 0 };
	for (int i = 0; i < 3; i++)
	{
		numBox[i] = inSize[i] / boxSize;
	}
	//setup the output image
	ImageType::Pointer image = reader->GetOutput();
	RGBImageTypePointer rgbEigen = RGBImageType::New();
	unsigned char startcolor[3] = { 0, 0, 0 };
	rgbEigen->SetRegions(image->GetLargestPossibleRegion());
	rgbEigen->SetOrigin(image->GetOrigin());
	rgbEigen->SetSpacing(image->GetSpacing());
	rgbEigen->Allocate();
	rgbEigen->FillBuffer(startcolor);

	cout << "number of boxes " << numBox[0] << "," << numBox[1] << "," << numBox[2] << endl;
	ImageType::IndexType start;
	ImageType::SizeType size;
	int count = 0;
	for (int z = 0; z < numBox[2]; z++)
	for (int y = 0; y < numBox[1]; y++)
	for (int x = 0; x < numBox[0]; x++)
	{
		//debug
		if (count < 3 || 1){
			start[0] = x*boxSize;
			start[1] = y*boxSize;
			start[2] = z*boxSize;
			cout << "start " << start << endl;
			//unsigned char color[3]={0.0};
			unsigned char color[3] = { 0.0 };
			for (int d = 0; d < 3; d++)
			{
				color[d] = abs(directionsIndices[3 * count + d]) * 255;
				//color[d]=directionsIndices[ 3*count+d];
			}

			for (int i = 0; i < boxSize; i++)
			for (int j = 0; j < boxSize; j++)
			for (int k = 0; k < boxSize; k++)
			{
				RGBImageType::IndexType idx;
				idx[0] = start[0] + i;
				idx[1] = start[1] + j;
				idx[2] = start[2] + k;
				rgbEigen->SetPixel(idx, color);
			}
			count++;
		}
	}

	//
	typedef itk::ImageFileWriter< RGBImageType > RGBWriterType;
	RGBWriterType::Pointer   writer = RGBWriterType::New();
	writer->SetInput(rgbEigen);

	writer->SetFileName("finalColor.mhd");
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject e)
	{
		std::cerr << "Error: " << e << std::endl;
	}

}

void write_fiber_n(vector<FIBER> &bundle, const int n, const int  clusterIndex, ofstream & fo);

// Compute Hessian / eigen analysis
void hessian_based_computation(INPUT_PARAMS * input,
	vector<FIBER> &bundle)
{
	//read an image
	ReaderType::Pointer reader = ReaderType::New();
	// reader->SetFileName( "/Users/arindambhattacharya/Research/internship/data/CFK-Prepreg-klein-570x356x374-1umVS-16bit.mhd");
	// full data
	//reader->SetFileName( "C:/Users/p41123/Documents/internship/kul/Kul/KUL.mhd");

	//reader->SetFileName( "/Volumes/IAMRNDM/Internship/data/crop-7/crop-7.mhd");
	//reader->SetFileName( "C:\\Users\\p41123\\Documents\\internship\\data\\output\\test-6\\hard-testcase-1.mhd");

	//set reader

	//

	//
	//crop-8
	string crop8 = "K:\\Internship\\data\\crop-8\\crop-8.mhd";
	string crop8_hess_eig0 = "K:\\Internship\\data\\crop-8\\crop-8-hess-eig0.mhd";
	//crop-9
	string crop9 = "K:\\Internship\\data\\crop-9\\crop-9.mhd";
	string crop9_hess_eig0 = "K:\\Internship\\data\\crop-9\\crop-9-hess-eig0.mhd";




	// simulation
	string sim1 = "K:\\Internship\\data\\simulation\\Datensatz_simuliert.mhd";
	string sim1_hess_eig0 = "K:\\Internship\\data\\simulation\\Datensatz_simuliert-hess-eig0.mhd";


	//crop-12
	string crop12 = "K:\\Internship\\data\\crop-12\\crop-12.mhd";
	string crop12_hess_eig0 = "K:\\Internship\\data\\crop-12\\crop-12-hess-eig0.mhd";

	//crop-13
	string crop13 = "D:\\ABhattacharya\\Internship\\data\\crop-13\\crop-13.mhd";

	string out = "D:\\ABhattacharya\\Internship\\data\\crop-13\\crop-out.mhd";

	//crop-10
	string crop10 = "K:\\Internship\\data\\crop-10\\crop-10.mhd";
	string crop10_hess_eig0 = "K:\\Internship\\data\\crop-10\\crop-10-hess-eig0.mhd";
	//string out="K:\\Internship\\data\\crop-10\\output\\tracks-8-2-b.mhd";


	//crop-11
	string crop11 = "K:\\Internship\\data\\crop-11\\crop-11.mhd";
	string crop11_hess_eig0 = "K:\\Internship\\data\\crop-11\\crop-11-hess-eig0.mhd";

	string set2 = "K:\\Internship\\data\\vgtestVolume_686x415x540_u16bit_2um.mhd";
	//string out  = "K:\\Internship\\data\\vgtestVolume_686x415x540_u16bit_2um_out.mhd";
	//string kul = "K:\\Internship\\data\\KUL.mhd";
	//string out="C:\\Users\\p41123\\Documents\\internship\\output\\kul-tracks-8-9.mhd";

	string crop15a = "D:\\ABhattacharya\\Internship\\data\\crop-15-a\\crop-15-a.mhd";
	//string crop13_hess_eig0="K:\\Internship\\data\\crop-13\\crop-13-hess-eig0.mhd";

	string crop16 = "D:\\ABhattacharya\\Internship\\data\\crop-16\\crop-16.mhd";
	//dataset 2
	string crop6 = "D:\\ABhattacharya\\Internship\\data\\crop-6\\crop-6.mhd";
	string crop7 = "D:\\ABhattacharya\\Internship\\data\\crop-7\\crop-7.mhd";

	string kul = "D:\\ABhattacharya\\Internship\\data\\KUL.mhd";
	string glass = "D:\\ABhattacharya\\Internship\\data\\glassfibers\\Layerstructure\\cut2-GF_tool_Platte1-2um-part1.mhd";

	const string A = "D:/ABhattacharya/Internship/data/BrenardData1Crops/crop1/crop1.mhd";
	const char * B = input->inputFileName.c_str();

	
	reader->SetFileName(B);
	reader->Update();
	cerr << "status: reading complete." << endl;

	// aniso diff
	// put input into
	typedef GradientAnisotropicDiffusionImageFilter< ImageType, ImageType> aniso_filter_type;
	aniso_filter_type::Pointer aniso_filter = aniso_filter_type::New();
	if (input->aniso_diff){
		aniso_filter->SetInput(reader->GetOutput());
		int numberOfIterations = 5;
		aniso_filter->SetNumberOfIterations(numberOfIterations);
		aniso_filter->SetTimeStep(0.0625);
		aniso_filter->SetConductanceParameter(3.0);
		aniso_filter->Update();
	}

	//hessian
	HessianFilterType::Pointer hessian = HessianFilterType::New();
	if (input->aniso_diff){
		hessian->SetInput(aniso_filter->GetOutput());
		cout << "Using the aniso diffusion before the hessian " << endl;
	}
	else
	{
		cout << "NO aniso diffusion before the hessian " << endl;
		hessian->SetInput(reader->GetOutput());
	}


	hessian->SetSigma(input->sigma);
	hessian->Update();

	cerr << "status: Hessian complete." << endl;

	// Eigen Computation
	typedef itk::Matrix< double, Dimension> EigenValuesArrayType;
	typedef itk::Matrix< double, Dimension, Dimension > EigenVectorMatrixType;

	//aniso image store
	ImageType::Pointer anisoPnt = aniso_filter->GetOutput();
	// Get info from the hessian image
	HessianFilterType::OutputImageType::Pointer HessianOutputImage = hessian->GetOutput();
	RegionType region = HessianOutputImage->GetLargestPossibleRegion();
	SpacingType spacing = HessianOutputImage->GetSpacing();
	PointType origin = HessianOutputImage->GetOrigin();

	// eigen vector /values place holders
	DoubleMatrixType eigenMatrix, myMatrix;
	eigenMatrix[0][0] = 0; eigenMatrix[0][1] = 0; eigenMatrix[0][2] = 0;
	eigenMatrix[1][0] = 0; eigenMatrix[1][1] = 0; eigenMatrix[1][2] = 0;
	eigenMatrix[2][0] = 0; eigenMatrix[2][1] = 0; eigenMatrix[2][2] = 0;

	DoubleVectorType eigenVal;
	eigenVal[0] = eigenVal[1] = eigenVal[2] = 0;

	// Allocate the output image for the hessian 
	DoubleMatrixImagePointer m_Output = DoubleMatrixImageType::New();
	m_Output->SetRegions(region);
	m_Output->SetOrigin(origin);
	m_Output->SetSpacing(spacing);
	m_Output->Allocate();
	m_Output->FillBuffer(eigenMatrix);

	//possible seeds
	boolImagePointer possibleSeeds = boolImageType::New();
	possibleSeeds->SetRegions(region);
	possibleSeeds->SetOrigin(origin);
	possibleSeeds->SetSpacing(spacing);
	possibleSeeds->Allocate();
	possibleSeeds->FillBuffer(false);



	//start locations 
	intImageTypePointer startLocs = intImageType::New();
	startLocs->SetRegions(region);
	startLocs->SetOrigin(origin);
	startLocs->SetSpacing(spacing);
	startLocs->Allocate();
	startLocs->FillBuffer(1);


	//
	//start locations 
	intImageTypePointer clusterIDS = intImageType::New();
	clusterIDS->SetRegions(region);
	clusterIDS->SetOrigin(origin);
	clusterIDS->SetSpacing(spacing);
	clusterIDS->Allocate();
	clusterIDS->FillBuffer(0);

	//Test output
	RGBImageTypePointer rgbEigen = RGBImageType::New();
	unsigned char color[3] = { 0, 0, 0 };
	rgbEigen->SetRegions(region);
	rgbEigen->SetOrigin(origin);
	rgbEigen->SetSpacing(spacing);
	rgbEigen->Allocate();
	rgbEigen->FillBuffer(color);


	// set up iterators
	HessianIteratorType mIt(HessianOutputImage, region);
	//iterator for storing the hessian
	DoubleIteratorType oIt(m_Output, region);
	RGBIteratorType rgbIt(rgbEigen, region);
	boolIteratorType seedIt(possibleSeeds, region);
	intIteratorType startLocsIt(startLocs, region);

	InIteratorType inIt(aniso_filter->GetOutput(), region);


	itk::SymmetricSecondRankTensor<double, 3> myTensor;
	mIt.GoToBegin();
	oIt.GoToBegin();
	rgbIt.GoToBegin();
	inIt.GoToBegin();
	seedIt.GoToBegin();
	//startLocsIt.GoToBegin();


	EigenVectorAnalysisType eig;
	eig.SetDimension(Dimension);
	eig.SetOrderEigenMagnitudes(true);

	cerr << "status: set for eigen computation." << endl;
	cout << region.GetSize()[0] << "," << region.GetSize()[1] << "," << region.GetSize()[2] << endl;
	int ii = 0;

	// for debug 
	vector <int> tempIndex(3, 0);
	double expThSecOrderness = 40.0;
	//datafile << "x,y,z\n";
	while (!inIt.IsAtEnd())
	{
		eig.SetOrderEigenMagnitudes(true);
		// Compute eigen values and eigen matrix
		color[0] = 0; color[1] = 0; color[2] = 0;

		if (inIt.Get() > input->thresh){

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

			// compute the eigen values and eigen vectors
			eig.ComputeEigenValuesAndVectors(myMatrix, eigenVal, eigenMatrix);

			// Compute vesseleness
			//// vness has been turned off 
			double vness = 0.0;
			//debug 
			//cout <<"eigen val "<< endl;
			//cout <<eigenVal[0]<<","<<eigenVal[1]<<","<<eigenVal[2]<<endl;

			//compute vesselness
			computeVesselness(eigenVal, eigenMatrix, false, expThSecOrderness, vness);
			// if the input is greater than the threhold.
			if (vness > input->vness)
			{
				for (int i = 0; i < Dimension; i++)
				{
					//color[i]=120*abs(eigenMatrix[0][i]);
					//color according to grayscale value
					//color[i] = 0.1*inIt.Get();
				}
				oIt.Set(eigenMatrix);
				seedIt.Set(true);
				unsigned tempCol[3] = { 27, 158, 119 };
				if (!input->colorTrack)
				for (int i = 0; i < Dimension; i++)
				{
					color[i] = 219 * abs(eigenMatrix[0][i]);
				}
				rgbIt.Set(color); // set the color of the eigen matrix
			}
			//else
			//{
			//	for (int i=0; i<Dimension; i++)
			//	{
			//		color[i]=255*abs(eigenMatrix[0][i]);
			//		//color[i]=255;
			//	}
			//}
		}
		ImageType::IndexType id = rgbIt.GetIndex();
		//cerr <<" aniso value " <<aniso_filter->GetOutput()->GetPixel(id);
		//cerr <<" anis val " << anisoPnt->GetPixel(id)<<endl;


		++oIt;
		++mIt;
		++rgbIt;
		++inIt;
		++seedIt;
	}
	cerr << "All hessians completed" << endl;
	/*
	for (int i=0;i<vnessVec.size();i++)
	vnessFile <<vnessVec[i]<<"\n";
	*/
	// track a point 
	if (input->trackPoint){
		ImageType::IndexType idx;
		int ii, jj, kk;
		ii = 162;
		jj = 100;
		kk = 244;


		//crop6 150,175,150, 15 17 15 10 spacing
		//crop 13 200 75 154 20 15 20 spacings 5
		//crop 7 175 175 150 17 15 17 10
		//glass 250 250 250 10 20 20 20 
		idx[0] = ii;
		idx[1] = jj;
		idx[2] = kk;

		vector<int> vidx(3, 0);
		vidx[0] = ii;
		vidx[1] = jj;
		vidx[2] = kk;
		//
		int spacings = 4;
		vector <int> radius(3, 0);
		//radius[0]=20;
		//radius[1]=10;
		//radius[2]=15;
		//test
		radius[0] = 7;//40;
		radius[1] = 7;//22;
		radius[2] = 7;//34;
		//3,3,10
		//20,11,17
		vector <vector<int> >nearPts;
		//DEBUG
		//findNearPts(vidx, spacings, radius, nearPts);

		vector<int> nearPtTemp(3, 0);
		nearPtTemp[0] = 162;
		nearPtTemp[1] = 100;
		nearPtTemp[2] = 68;

		//findNearPts(nearPtTemp, spacings, radius, nearPts);
		/*
		nearPtTemp[0] = 162; nearPtTemp[1] = 70; nearPtTemp[2] = 84;
		findNearPts(nearPtTemp, spacings, radius, nearPts);

		nearPtTemp[0] = 162; nearPtTemp[1] = 70; nearPtTemp[2] = 249;
		findNearPts(nearPtTemp, spacings, radius, nearPts);


		nearPtTemp[0] = 162; nearPtTemp[1] = 70; nearPtTemp[2] = 167;
		findNearPts(nearPtTemp, spacings, radius, nearPts);

		nearPtTemp[0] = 162; nearPtTemp[1] = 70; nearPtTemp[2] = 42;
		findNearPts(nearPtTemp, spacings, radius, nearPts);

		nearPtTemp[0] = 162; nearPtTemp[1] = 115; nearPtTemp[2] = 263;
		findNearPts(nearPtTemp, spacings, radius, nearPts);

		nearPtTemp[0] = 162; nearPtTemp[1] = 116; nearPtTemp[2] = 47;
		findNearPts(nearPtTemp, spacings, radius, nearPts);

		nearPtTemp[0] = 162; nearPtTemp[1] = 27; nearPtTemp[2] = 161;
		findNearPts(nearPtTemp, spacings, radius, nearPts);


		nearPtTemp[0] = 45; nearPtTemp[1] = 129; nearPtTemp[2] = 167;
		findNearPts(nearPtTemp, spacings, radius, nearPts);

		nearPtTemp[0] = 43; nearPtTemp[1] = 125; nearPtTemp[2] = 155;
		findNearPts(nearPtTemp, spacings, radius, nearPts);

		nearPtTemp[0] = 342; nearPtTemp[1] = 71; nearPtTemp[2] = 167;
		findNearPts(nearPtTemp, spacings, radius, nearPts);

		nearPtTemp[0] = 299; nearPtTemp[1] = 47; nearPtTemp[2] = 155;
		findNearPts(nearPtTemp, spacings, radius, nearPts);
		*/
		//
		nearPtTemp[0] = 162;
		nearPtTemp[1] = 70;
		nearPtTemp[2] = 293;

		nearPts.push_back(nearPtTemp);

		nearPtTemp[0] = 162;
		nearPtTemp[1] = 70;
		nearPtTemp[2] = 285;

		nearPts.push_back(nearPtTemp);
		//
		nearPtTemp[0] = 176;
		nearPtTemp[1] = 70;
		nearPtTemp[2] = 64;

		nearPts.push_back(nearPtTemp);

		//read the start points if any from the config file.
		vector <vector<int> > configStartPts;
		if (!input->startLocations.empty())
		{
			setStartPointsFromConfig(input->startLocations, configStartPts);
			cerr << "Total number of start points : " << configStartPts.size() << endl;
		}
		

		cerr << "number of seeds " << nearPts.size() << endl;
		int tot = configStartPts.size();
		//
		int j = 1; // id of the pixel

		cerr << " Compute fibers time: start  " << currentDateTime() << endl;
		for (int i = 0; i< configStartPts.size(); i++){
			idx[0] = configStartPts[i][0];
			idx[1] = configStartPts[i][1];
			idx[2] = configStartPts[i][2];


			cerr << (i*1.0 / configStartPts.size()) * 100.0 << " percent done [" << i << " out of " << configStartPts.size() << "]";
			cerr << " idx " << idx[0] << " " << idx[1] << " " << idx[2] << endl;
				//cout <<i<<" out of "<<tot <<endl; 
			//cout <<"start " << idx << " value " << aniso_filter->GetOutput()->GetPixel(idx) << endl;
			if (region.IsInside(idx))
			if (aniso_filter->GetOutput()->GetPixel(idx) > input->thresh && possibleSeeds->GetPixel(idx))
			{
				trackPoint2Debug( m_Output, input, idx, j, rgbEigen, possibleSeeds, bundle);
				startLocs->SetPixel(idx, 255);
			}
		}
		cerr << " Compute fibers time: end  " << currentDateTime() << endl;
		cerr << "number of small fibers " << smallFiberNum << endl;
		cerr << "number of large fibers " << j - 1 << endl;


		averageOrientationEndPointsBundle(bundle);
		correctOrientationEndPointsBundle(bundle);
		//Correct the orientation of the fibers. 
		//correctFiberBundleOrientaion(bundle);
		//Compute Average Orientation
		//computeAverageOrientation(bundle);


		// all the fibers are tracked 
		cerr << "color track . " << int(input->colorTrack) << endl;
		if (input->colorTrack)
		{
			cout << "color track id " << input->colorTrackID << endl;
			if (input->colorTrackID == -1)
			{
				//color all tracks
				if (input->computeKmeans)
				{
					string cluster("D:\\ABhattacharya\\Internship\\output\\cluster2014.csv");
					cerr << "*********************** read data from cluster.csv***********" << endl;
					ifstream in(cluster.c_str());
					if (!in.is_open()) { cerr << "error in reading cluster \n"; exit(0); };

					typedef tokenizer< escaped_list_separator<char> > Tokenizer;
					vector <int> clusterIndices;
					vector< string > vec;
					string line;
					while (getline(in, line))
					{
						Tokenizer tok(line);
						vec.assign(tok.begin(), tok.end());
						/*
						copy(vec.begin(), vec.end(),
						ostream_iterator<string>(cout, "|"));
						*/
						vector<string>::iterator it = vec.begin();
						istringstream os(*it);
						int d;
						os >> d;
						clusterIndices.push_back(d);
					}
					cerr << "color coding by cluster " << endl;
					cerr << "size " << clusterIndices.size() << endl;

					//color according to cluster
					for (int i = 0; i < bundle.size(); i++)
					{
						if (input->clusterSpecific){
							// Color only specific clusters based on config file
							std::vector<int>::iterator it;
							it = find(input->clusterCols.begin(), input->clusterCols.end(),
								clusterIndices[i]);
							++it;
							if (it != input->clusterCols.end())
								colorTrack(rgbEigen, bundle, clusterIndices[i], i);
						}
						else
						{
							//color all clusters.
							colorTrack(rgbEigen, bundle, clusterIndices[i], i);
						}

					}
					// debug
					// color according to cluster 
					cerr << " in color cluster " << endl;
					ofstream fo("fiber_out.txt");
					for (int i = 0; i < bundle.size(); i++)
					{
						colorTrackIDCluster(clusterIDS, bundle, i, clusterIndices[i]);
						write_fiber_n(bundle, i, clusterIndices[i], fo);
						fo << endl;
						//write_fiber_n(bundle,i,clusterIndices[i], fo);

					}

					cout << " before voxelize  Here" << endl;
					// ***DEBUG***
					/*
					voxelizeSpace(clusterIDS, bundle, possibleSeeds,
					anisoPnt);

					cout << "voxelize complete" << endl;
					*/

					/*
					// ***DEBUG***
					// different approach
					// compute binary image
					// extract surface
					for(int l=0;l<14;l++)
					{
					// create a binary image
					typedef itk::BinaryThresholdImageFilter
					<intImageType,ucharImageType > Binary_FilterType;
					Binary_FilterType::Pointer binary_filter = Binary_FilterType::New();
					binary_filter->SetInput(clusterIDS);
					binary_filter->SetOutsideValue( 0 );
					binary_filter->SetInsideValue( 255 );

					int current_label = l;
					binary_filter->SetLowerThreshold( current_label );
					binary_filter->SetUpperThreshold( current_label );
					binary_filter->Update();
					// binary closing operator
					cout <<"binary closing \n";
					typedef itk::BinaryBallStructuringElement
					< unsigned char, 3> KernelType;
					KernelType ball;
					KernelType::SizeType ballSize;
					ballSize.Fill(10 );
					ball.SetRadius(ballSize);
					ball.CreateStructuringElement();

					typedef itk::BinaryMorphologicalClosingImageFilter
					< ucharImageType,ucharImageType, KernelType > CloseType;
					CloseType::Pointer close = CloseType::New();
					close->SetInput( binary_filter->GetOutput() );
					close->SetKernel( ball );
					//close->SetSafeBorder( false );
					close->Update();
					//

					typedef itk::BinaryMorphologicalOpeningImageFilter
					< ucharImageType,ucharImageType, KernelType > OpenType;
					OpenType::Pointer open = OpenType::New();
					open->SetInput( close->GetOutput() );
					open->SetKernel( ball );

					open->Update();


					//connected components
					typedef itk::ConnectedComponentImageFilter<ucharImageType, ucharImageType >
					ConnectedComponentImageFilterType;

					ConnectedComponentImageFilterType::Pointer connected =
					ConnectedComponentImageFilterType::New ();
					connected->SetInput(open->GetOutput());
					connected->Update();

					std::cout << "Number of objects: " << connected->GetObjectCount() << std::endl;
					//find the largest one
					typedef itk::LabelShapeKeepNObjectsImageFilter< ucharImageType > LabelShapeKeepNObjectsImageFilterType;
					LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter = LabelShapeKeepNObjectsImageFilterType::New();
					labelShapeKeepNObjectsImageFilter->SetInput( connected->GetOutput() );
					labelShapeKeepNObjectsImageFilter->SetBackgroundValue( 0 );
					labelShapeKeepNObjectsImageFilter->SetNumberOfObjects( 1 );
					labelShapeKeepNObjectsImageFilter->SetAttribute( LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);

					//extract mesh
					typedef itk::Mesh<double> MeshType;
					typedef itk::BinaryMask3DMeshSource< ucharImageType, MeshType > MeshSourceType;
					MeshSourceType::Pointer meshSource = MeshSourceType::New();
					meshSource->SetObjectValue(1);
					meshSource->SetInput(labelShapeKeepNObjectsImageFilter->GetOutput());

					try
					{
					meshSource->Update();
					}
					catch( itk::ExceptionObject & exp )
					{
					std::cerr << "Exception thrown during Update() " << std::endl;
					std::cerr << exp << std::endl;
					exit(0);
					}

					// write the mesh out
					std::string filePrefix = "Mesh_";
					typedef itk::VTKPolyDataWriter<MeshType> MeshWriterType;
					MeshWriterType::Pointer MeshWriter = MeshWriterType::New();
					MeshWriter->SetInput(meshSource->GetOutput());
					vtksys_stl::stringstream ss;
					ss << filePrefix << l << ".vtk";
					cout << " writing " << ss.str() << endl;

					MeshWriter->SetFileName(ss.str().c_str());
					MeshWriter->Update();
					}
					*/
					//surface extraction 
					///SURFACE EXTRACTION CURRENT CODE
					/*
					typedef itk::ImageToVTKImageFilter<intImageType> itk2vtkType;
					itk2vtkType::Pointer itk2vtk = itk2vtkType::New();
					// ***DEBUG***
					itk2vtk->SetInput(clusterIDS);
					//itk2vtk->SetInput(labelShapeKeepNObjectsImageFilter->GetOutput());
					itk2vtk->Update();
					double* origin = itk2vtk->GetOutput()->GetOrigin();
					int* extent = itk2vtk->GetOutput()->GetExtent();

					std::cout << extent[0] << " " << extent[1] << " "
					<< extent[2] << " " << extent[3] << " "
					<< extent[4] << " " << extent[5] << std::endl;

					vtkSmartPointer<vtkDiscreteMarchingCubes> discreteCubes = vtkSmartPointer<vtkDiscreteMarchingCubes>::New();

					//DEBUG
					discreteCubes->SetInputData(itk2vtk->GetOutput());
					discreteCubes->GenerateValues(14, 1, 14); //TODO Fixed numbers
					discreteCubes->Update();
					vtkSmartPointer<vtkPolyData> cubes = discreteCubes->GetOutput();
					cout << " here " << cubes->GetNumberOfPoints() << endl;

					//
					std::string filePrefix = "euro_label";
					unsigned int smoothingIterations = 15;
					double passBand = 0.001;
					double featureAngle = 120.0;
					vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
					vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
					vtkSmartPointer<vtkMaskFields> scalarsOff =
					vtkSmartPointer<vtkMaskFields>::New();
					vtkSmartPointer<vtkGeometryFilter> geometry =
					vtkSmartPointer<vtkGeometryFilter>::New();
					vtkSmartPointer<vtkXMLPolyDataWriter> writer =
					vtkSmartPointer<vtkXMLPolyDataWriter>::New();
					smoother->SetInputConnection(discreteCubes->GetOutputPort());
					smoother->SetNumberOfIterations(smoothingIterations);
					smoother->BoundarySmoothingOff();
					smoother->FeatureEdgeSmoothingOff();
					smoother->SetFeatureAngle(featureAngle);
					smoother->SetPassBand(passBand);
					smoother->NonManifoldSmoothingOn();
					smoother->NormalizeCoordinatesOn();
					smoother->Update();
					vtkSmartPointer<vtkThreshold> selector =
					vtkSmartPointer<vtkThreshold>::New();
					selector->SetInputConnection(smoother->GetOutputPort());
					selector->SetInputArrayToProcess(0, 0, 0,
					vtkDataObject::FIELD_ASSOCIATION_CELLS,
					vtkDataSetAttributes::SCALARS);

					// Strip the scalars from the output
					scalarsOff->SetInputConnection(selector->GetOutputPort());
					scalarsOff->CopyAttributeOff(vtkMaskFields::POINT_DATA,
					vtkDataSetAttributes::SCALARS);
					scalarsOff->CopyAttributeOff(vtkMaskFields::CELL_DATA,
					vtkDataSetAttributes::SCALARS);

					geometry->SetInputConnection(scalarsOff->GetOutputPort());

					writer->SetInputConnection(geometry->GetOutputPort());

					for (unsigned int i = 1; i <= 14; i++)
					{
					// see if the label exists, if not skip it


					// select the cells for a given label
					selector->ThresholdBetween(i, i);

					// output the polydata
					vtksys_stl::stringstream ss;
					ss << filePrefix << i << ".vtp";
					cout << " writing " << ss.str() << endl;

					writer->SetFileName(ss.str().c_str());
					writer->Write();

					}

					*/

					//

					//debug
					if (input->colorVol)
					{

						RGBIntImageTypePointer rgbVol = RGBPixelIntImageType::New();
						unsigned char color[3] = { 0, 0, 0 };
						rgbVol->SetRegions(region);
						rgbVol->SetOrigin(origin);
						rgbVol->SetSpacing(spacing);
						rgbVol->Allocate();
						rgbVol->FillBuffer(color);
						cerr << "com[puting color volume" << endl;
						colorVolume(rgbEigen, rgbVol, possibleSeeds, anisoPnt);
						cerr << "writing color volume" << endl;
						typedef itk::ImageFileWriter< RGBPixelIntImageType > RGBWriterType;
						RGBWriterType::Pointer   writer = RGBWriterType::New();
						writer->SetInput(rgbVol);
						writer->SetFileName(input->bundleVolFname);
						try
						{
							writer->Update();
						}
						catch (itk::ExceptionObject e)
						{
							std::cerr << "Error: " << e << std::endl;
						}

					}
				}
				else
					// no kmeans
				{
					if (input->colorSpecificFibers)
					{
						cerr << "Coloring specific fibers." << endl;
						for (int i = 0; i < input->indicesColorSpecificFibers.size(); i++)
						{
							int tempTrackId = input->indicesColorSpecificFibers[i] - 1;
							if (tempTrackId >= bundle.size())
								continue;// the  requested track id does not exist
							colorTrack(rgbEigen, bundle, input->indicesColorSpecificFibers[i] - 1);
							colorTrackID(clusterIDS, bundle, input->indicesColorSpecificFibers[i] - 1);
						}
					}
					else
					{
						for (int i = 0; i < bundle.size(); i++)
						{
							colorTrack(rgbEigen, bundle, i);
							//cout << "*** ClusterIDs Color  " << endl;
							colorTrackID(clusterIDS, bundle, i);
						}
					}
				}
			}
			else
			{
				//bundle[input->colorTrackID].colorFiber();
				cout << "colorTrackid != -1 " << endl;

				colorTrack(rgbEigen, bundle, input->colorTrackID);

			}
		}
		cout << "end hessian comp \n";
	}
	//compute reliable hessians.
	//along a small region the hessians must agree
	//computeReliableHessians(m_Output,input);

	//double sum = std::accumulate(secOnessVec.begin(), secOnessVec.end(), 0.0);
	//double mean = sum / secOnessVec.size();
	//cout <<"sec orderness mean  " << mean <<endl;
	//vector<float>::const_iterator it;
	//it = max_element(secOnessVec.begin(), secOnessVec.end());
	//cout <<"max elemeenyt  " << *it <<endl;
	cerr << "start writting" << endl;
	cout << "********* writting rgb Eigen to D:\\ABhattacharya\\Internship\\output\\outColor.mhd" << endl;
	if (bundle.size() >= 0){
		typedef itk::ImageFileWriter< RGBImageType > RGBWriterType;
		RGBWriterType::Pointer   writer = RGBWriterType::New();
		writer->SetInput(rgbEigen);
		if (!(input->computeKmeans))
			writer->SetFileName("D:\\ABhattacharya\\Internship\\output\\outColor2014.mhd");

		else
			writer->SetFileName("D:\\ABhattacharya\\Internship\\output\\outClusterResult.mhd");
		try
		{
			writer->Update();
		}
		catch (itk::ExceptionObject e)
		{
			std::cerr << "Error: " << e << std::endl;
		}



		//clusterIDS
		if (input->colorTrack){
			typedef itk::ImageFileWriter< intImageType > writerType2;
			writerType2::Pointer   writer_bundle = writerType2::New();
			writer_bundle->SetInput(clusterIDS);

			writer_bundle->SetFileName("D:\\ABhattacharya\\Internship\\output\\bundleID-clusterID2014.mhd");
			try
			{
				writer_bundle->Update();
			}
			catch (itk::ExceptionObject e)
			{
				std::cerr << "Error: " << e << std::endl;
			}
		}



	}
	else
		cerr << "writting nothing " << endl;


	//seed locations
	typedef itk::ImageFileWriter< intImageType > writerType;
	writerType::Pointer   writer_start_loc = writerType::New();
	writer_start_loc->SetInput(startLocs);

	writer_start_loc->SetFileName("D:\\ABhattacharya\\Internship\\output\\startLocs.mhd");
	try
	{
		writer_start_loc->Update();
	}
	catch (itk::ExceptionObject e)
	{
		std::cerr << "Error: " << e << std::endl;
	}


	cerr << "finished writting " << endl;
}

//  print information
void print_info(const DoubleVectorType eigenVal, const DoubleMatrixType eigenMatrix, const  InIteratorType inIt,
	const ImageType::IndexType idx, const unsigned char *color, const INPUT_PARAMS * input, double vness)
{
	// set things
	//color[0]= 255;color[1]= 0;color[2]= 0;
	cout << "\n--------------------\n";
	cout << "idx " << idx;
	//cout <<"pixel " << HessianOutputImage->GetPixel(idx) << endl;
	//cout <<"tubularness  " <<abs(eigenVal[1])/(abs(eigenVal[2])*1.0f) << endl;
	//cout <<"bug " << int (bug1) << " "<< int (bug2) << endl;
	//std::cout << "myTensor " << myTensor << std::endl;
	cout << "scalar " << inIt.Get() << "  th " << input->thresh << endl;
	std::cout << "color " << int(color[0]) << " " << int(color[1]) << " " << int(color[2]) << std::endl;
	//std::cout << myMatrix << std::endl;
	std::cout << "eigenVal " << eigenVal << std::endl;
	std::cout << "eigenMatrix \n" << eigenMatrix << std::endl;
	//computeVesselness(eigenVal, eigenMatrix, input->degugMode, vness);

}
//  print information
void print_info_small(const DoubleVectorType eigenVal, const DoubleMatrixType eigenMatrix, const  InIteratorType inIt,
	const ImageType::IndexType idx, const unsigned char *color, const INPUT_PARAMS * input, double vness)
{
	// set things
	//color[0]= 255;color[1]= 0;color[2]= 0;

	cout << "idx " << idx;
	//cout <<"pixel " << HessianOutputImage->GetPixel(idx) << endl;
	//cout <<"tubularness  " <<abs(eigenVal[1])/(abs(eigenVal[2])*1.0f) << endl;
	//cout <<"bug " << int (bug1) << " "<< int (bug2) << endl;
	//std::cout << "myTensor " << myTensor << std::endl;
	cout << "scalar " << inIt.Get() << "  th " << input->thresh << endl;
	//std::cout << "color " << int(color [0]) << " " << int(color [1]) <<" "<< int(color [2]) << std::endl;
	//std::cout << myMatrix << std::endl;
	std::cout << "eigenVal " << eigenVal << std::endl;
	std::cout << "eigenMatrix " << eigenMatrix[0][0] << "," << eigenMatrix[0][1] << "," << eigenMatrix[0][2] << "   "
		<< eigenMatrix[1][0] << "," << eigenMatrix[1][1] << "," << eigenMatrix[1][2] << "   " <<
		eigenMatrix[2][0] << "," << eigenMatrix[2][1] << "," << eigenMatrix[2][2] << std::endl;
	//computeVesselness(eigenVal, eigenMatrix, input->degugMode, vness);
	cout << "\n--------------------\n";

}




void computeKmeansOutput(const INPUT_PARAMS * input)
{
	string cluster("cluster.csv");
	cout << "read data from cluster.csv" << endl;
	ifstream in(cluster.c_str());
	if (!in.is_open()) { cerr << "error in reading,"; exit(0); };

	typedef tokenizer< escaped_list_separator<char> > Tokenizer;
	vector <int> clusterIndices;
	vector< string > vec;
	string line;
	while (getline(in, line))
	{
		Tokenizer tok(line);
		vec.assign(tok.begin(), tok.end());
		/*
		copy(vec.begin(), vec.end(),
		ostream_iterator<string>(cout, "|"));
		*/
		vector<string>::iterator it = vec.begin();
		istringstream os(*it);
		int d;
		os >> d;
		clusterIndices.push_back(d);
	}

	for (int i = 0; i < clusterIndices.size(); i++)
	{
		cout << " " << clusterIndices[i];
	}
}
// 
// compute the kmeans results by plotting backto the image
// old code for reading back the cluster results

void computeKmeansOutputOLd(const INPUT_PARAMS * input)
{
	string cluster("cluster.csv");
	cout << "read data from cluster.csv" << endl;
	ifstream in(cluster.c_str());
	if (!in.is_open()) { cerr << "error in reading,"; exit(0); };

	typedef tokenizer< escaped_list_separator<char> > Tokenizer;
	vector <int> clusterIndices;
	vector< string > vec;
	string line;
	while (getline(in, line))
	{
		Tokenizer tok(line);
		vec.assign(tok.begin(), tok.end());
		/*
		copy(vec.begin(), vec.end(),
		ostream_iterator<string>(cout, "|"));
		*/
		vector<string>::iterator it = vec.begin();
		istringstream os(*it);
		int d;
		os >> d;
		clusterIndices.push_back(d);
	}
	cerr << "reading of the ClusterIndices done";

	// read the coordinates
	string coords("index.csv");
	cout << "read data from index.csv" << endl;
	ifstream inCoords(coords.c_str());
	if (!inCoords.is_open()) { cerr << "error in reading,"; exit(0); };


	vector <vector<int> > coordsIndices;
	while (getline(inCoords, line))
	{
		Tokenizer tok(line);
		vec.assign(tok.begin(), tok.end());

		//copy(vec.begin(), vec.end(),
		//ostream_iterator<string>(cout, "|"));

		if (vec.size() < 3) continue;
		vector <int> tempCoord;
		for (std::vector<string>::iterator it = vec.begin(); it != vec.end(); ++it) {

			istringstream os(*it);
			int d;
			os >> d;
			tempCoord.push_back(d);
		}
		coordsIndices.push_back(tempCoord);
	}
	cerr << "reading of the coordinates done";
	ReaderType::Pointer reader = ReaderType::New();
	//string crop = "K:\\Internship\\data\\simulation\\Datensatz_simuliert.mhd"; //"K:\\Internship\\data\\crop-8\\crop-8.mhd";
	string crop = "K:\\Internship\\data\\prepreg\\crop-3\\crop-3.mhd";

	reader->SetFileName(crop);
	reader->Update();

	RegionType region = reader->GetOutput()->GetLargestPossibleRegion();
	SpacingType spacing = reader->GetOutput()->GetSpacing();
	PointType origin = reader->GetOutput()->GetOrigin();
	//set output color
	RGBImageTypePointer rgbCluster = RGBImageType::New();
	unsigned char color[3] = { 0.0, 0.0, 0.0 };
	rgbCluster->SetRegions(region);
	rgbCluster->SetOrigin(origin);
	rgbCluster->SetSpacing(spacing);
	rgbCluster->Allocate();
	rgbCluster->FillBuffer(color);
	InIteratorType inIt(reader->GetOutput(), region);
	RGBIteratorType rgbIt(rgbCluster, region);
	inIt.GoToBegin();
	rgbIt.GoToBegin();
	while (!inIt.IsAtEnd())
	{
		ImageType::IndexType idx = inIt.GetIndex();
		color[0] = inIt.Get(); color[1] = inIt.Get(); color[2] = inIt.Get();
		//rgbIt.Set(color);
		inIt++;
		rgbIt++;
	}

	unsigned char red[3] = { 100, 0, 0 };
	unsigned char green[3] = { 0, 100, 0 };
	unsigned char blue[3] = { 0, 0, 100 };
	unsigned char purple[3] = { 100, 0, 100 };
	cerr << "setting color " << endl;

	for (int i = 0; i < coordsIndices.size(); i++)
	{
		ImageType::IndexType idx;
		//cout <<coordsIndices[i][0]<<","<<coordsIndices[i][1]<<","<<coordsIndices[i][2]
		//<<"-- "<<clusterIndices[i]<<endl;

		for (int j = 0; j < 3; j++)
			idx[j] = coordsIndices[i][j];

		if (clusterIndices[i] == 1)
			rgbCluster->SetPixel(idx, red);
		if (clusterIndices[i] == 2)
			rgbCluster->SetPixel(idx, blue);

		if (clusterIndices[i] == 3)
			rgbCluster->SetPixel(idx, green);
		if (clusterIndices[i] == 4)
			rgbCluster->SetPixel(idx, purple);
		/*
		if (clusterIndices[i]==2 || clusterIndices[i]==3)
		rgbCluster->SetPixel(idx,blue);
		*/
		/*
		else if (clusterIndices[i]==2)
		rgbCluster->SetPixel(idx,green);
		else if (clusterIndices[i]==3)
		rgbCluster->SetPixel(idx,green);
		else  if (clusterIndices[i]==4)
		rgbCluster->SetPixel(idx,red);
		else
		cout <<"does not matcjh"<<endl;
		*/
	}


	//write it out
	typedef itk::ImageFileWriter< RGBImageType > RGBWriterType;
	RGBWriterType::Pointer   writer = RGBWriterType::New();
	writer->SetInput(rgbCluster);

	writer->SetFileName("RGBCluster.mhd");
	try
	{
		writer->Update();

	}
	catch (itk::ExceptionObject e)
	{
		std::cerr << "Error: " << e << std::endl;
	}

}

// self similarrity metric computation

void selfSimMetricComputation(const INPUT_PARAMS * input)
{
	cerr << "status: in function xbundle_aniso_diff" << endl;

	floatReaderType::Pointer reader = floatReaderType::New();
	reader->SetFileName("/Volumes/IAMRNDM/Internship/data/kul-crop-4.mhd");

	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &err) {
		std::cout << "ExceptionObject caught !" << std::endl; std::cout << err << std::endl;
		exit(0);
	}
	//
	typedef GradientAnisotropicDiffusionImageFilter< floatImageType, floatImageType> aniso_filter_type;
	aniso_filter_type::Pointer aniso_filter = aniso_filter_type::New();

	// put input into
	aniso_filter->SetInput(reader->GetOutput());
	int numberOfIterations = 5;
	aniso_filter->SetNumberOfIterations(numberOfIterations);
	aniso_filter->SetTimeStep(0.0625);
	aniso_filter->SetConductanceParameter(3.0);
	aniso_filter->Update();

}


// Function to remove very small fibers
void remove_small_fibers(const INPUT_PARAMS  * input, vector<FIBER> &bundle, vector<FIBER> &revisedBundle)
{
	for (int m = 0; m<bundle.size(); m++)
	{

		if (bundle[m].length > input->minFiberLength){
			//
			//cout <<"fiber number "<< bundle[m].id;
			//cout <<"  length of the fiber " << bundle[m].length <<endl;
			//for(int o=0;o<bundle[m].points.size();o++)
			//{
			//	cout <<"["<<o<<"]"<<bundle[m].points[o][0]<<","<<bundle[m].points[o][1]
			//	<<","<<bundle[m].points[o][2]<<" dist " << bundle[m].perDist[o]<<" dir "
			//		<<bundle[m].dir[o][0]<<","
			//		<<bundle[m].dir[o][1]<<","
			//		<<bundle[m].dir[o][2]<<","
			//		<<endl;
			//}

			//
			revisedBundle.push_back(bundle[m]);
		}
	}
}

//help functions : compute the distances between fibers

void compare(vector<FIBER> &bundle, EDGE &e, const int a, const int b, float & d);
//Advanced distance computations
inline void compareAdvanced(vector<FIBER> &bundle, EDGE &e, const int a, const int b, float & d);


//
//Compute the distances between fiberss
//
void compute_distances(const INPUT_PARAMS  * input, vector<FIBER> &bundle, vector<EDGE> &graph)
{
	//DEBUG
	cout << "Compute Distances" << endl;
	int M = bundle.size();
	for (int m = 0; m < M - 1; m++)
	{
		if (int((m*1.0 / M) * 100) % 20 == 0)
			cerr << int((m*1.0 / M) * 100) << " percent done." << endl;

		for (int n = m + 1; n < M; n++)
		{
			/*cout <<"**********************\n";*/

			/*cout << "comparing " << m << "[" << bundle[m].id << "]"
				<< "  with " << n << "[" << bundle[n].id << "] \n";*/
			float d = 0.0;
			EDGE a;
			//compare(bundle,a, m, n, d);
			compareAdvanced(bundle, a, m, n, d);
			//cout << "mean of mean of min  " << d << endl;
			//cout <<"loc 1["<< a.node1_loc <<"] loc 2 [" << a.node2_loc <<"]"<<endl; 

			a.node1 = bundle[m].id;
			a.node2 = bundle[n].id;
			float dotPdt = 0.0f;
			for (int i = 0; i < 3; i++)
			{
				dotPdt = dotPdt + (bundle[m].dir[a.node1_loc][i] * bundle[n].dir[a.node2_loc][i]);
			}

			if (dotPdt < 0.0)
			{
				dotPdt = 0.0;
				for (int i = 0; i < 3; i++)
				{
					dotPdt = dotPdt + (bundle[m].dir[a.node1_loc][i] * (-1.0)*bundle[n].dir[a.node2_loc][i]);
				}
			}
			//debug 
			/*for (int i=0;i<3;i++)
			cout << bundle[m].dir[a.node1_loc][i] <<" , ";
			cout <<endl;
			for (int i=0;i<3;i++)
			cout << bundle[n].dir[a.node2_loc][i] <<" , ";*/

			//cout <<"cos dist "<< dotPdt;
			// for d exp(-x^2/500)  
			// for angle just the cos
			//a.edgeW=exp(-(d)/500)*dotPdt;

			//a.edgeW=dotPdt*exp(-(d)/500);
			a.edgeW = dotPdt;
			///a.edgeW2=exp(-(sqrt(d))/700);
			//a.edgeW2=exp(-sqrt(d)/50.0);
			a.edgeW2 = int(d);
			/*cout <<" w2 "<<a.edgeW2<<endl ;
			cout <<"----------------------------\n";*/
			//a.edgeW=d;
			graph.push_back(a);
		}
	}
}
//
//
//
void write_graph(vector<EDGE> &graph, const string graphFileName)
{
	ofstream f(graphFileName);
	f << "n1,n2,ew,d\n";
	for (int i = 0; i < graph.size(); i++)
	{
		f << graph[i].node1 << "," << graph[i].node2 << "," << graph[i].edgeW << "," << graph[i].edgeW2 << endl;

	}
}


//void read_fibers (vector<FIBER> &bundle, const string bundleInfoFname)
//{
//
//	ifstream in (bundleInfoFname.c_str());
//	if (!in.is_open()) {cerr <<"error in reading,";exit(0);};
//
//	typedef tokenizer< escaped_list_separator<char> > Tokenizer;
//	vector <float> info;  
//	vector< string > vec;
//	string line;
//	while (getline(in,line))
//	{
//		Tokenizer tok(line);
//		vec.assign(tok.begin(),tok.end());
//		/*
//		copy(vec.begin(), vec.end(),
//		ostream_iterator<string>(cout, "|"));
//		*/
//		vector<string>::iterator it = vec.begin();
//		istringstream os(*it);
//		float d;
//		os >> d;
//		info.push_back(d);
//	}
//
//
//	int j=0;
//	int s=int(info[0]);
//		for (int i=0;i<s;i++)
//	{
//		FIBER f;
//		for (int i=0;i<f.points.size();i++)
//		{
//			fo <<f.points[i][0]<<","<<f.points[i][1]<<","<<f.points[i][2]<<",";
//		}
//		for (int i=0;i<f.dir.size();i++)
//		{
//			fo <<f.dir[i][0]<<","<<f.dir[i][1]<<","<<f.dir[i][2]<<",";
//		}
//	}
//}

//write fibers to file
void write_fibers(vector<FIBER> &bundle, const string bundleInfoFname)
{
	ofstream fo(bundleInfoFname);
	fo << bundle.size() << ",";
	for (int i = 0; i < bundle.size(); i++)
	{
		FIBER f = bundle[i];
		for (int i = 0; i < f.points.size(); i++)
		{
			fo << f.points[i][0] << "," << f.points[i][1] << "," << f.points[i][2] << ",";
		}
		/*
		for (int i=0;i<f.dir.size();i++)
		{
		fo <<f.dir[i][0]<<","<<f.dir[i][1]<<","<<f.dir[i][2]<<",";
		}
		*/
	}
}

void write_fiber_n(vector<FIBER> &bundle, const int n, const int  clusterIndex, ofstream & fo)
{
	FIBER f = bundle[n];
	for (int i = 0; i < f.points.size(); i++)
	{
		fo << f.points[i][0] << "," << f.points[i][1] << "," << f.points[i][2] << "," << clusterIndex << ",";
	}
	/*
	for (int i=0;i<f.dir.size();i++)
	{
	fo <<f.dir[i][0]<<","<<f.dir[i][1]<<","<<f.dir[i][2]<<",";
	}
	*/

}






void compare(vector<FIBER> &bundle, EDGE &e, const int a, const int b, float & d)
{
	float minD = 100000000.0;
	int loc1, loc2;
	for (int k = 0; k < bundle[a].points.size(); k++)
	{
		for (int l = 0; l<bundle[b].points.size(); l++)
		{
			float x = dist<int>(bundle[a].points[k], bundle[b].points[l]);
			if (minD>x)
			{
				minD = x;
				loc1 = k;
				loc2 = l;
			}
		}
	}
	d = minD;
	e.node1_loc = loc1;
	e.node2_loc = loc2;
	//cout <<loc1<<","<<loc2<<endl;
}

//compare 2  fibers called by compareAdvanced to compute dst as mentioned in 
//Identifying white matter fiber bundles in DTI data using an automated proximity based fiver clustering
// not used !!
/*
inline void compareFibers (vector<FIBER> &bundle, EDGE &e, const int a, const int b, float & d)
{
cout <<"cmp ";
float tot=0.0;
float combinedShortestDist=10000000.0;
int M=bundle[a].points.size();
int N=bundle[b].points.size();
cout <<"m " << M <<" n " << N << endl;
for (int k=0;k<M;k++)
{
double D=1000000.0;
for (int l=0;l<N;l++)
{
double tempD=sqrt(dist (bundle[a].points[k], bundle[b].points[l]));
if (tempD<D)
D=tempD;
if (tempD<combinedShortestDist)
{
e.node1_loc=k;
e.node2_loc=l;
combinedShortestDist = tempD;
}
}
tot=tot+D;
}
d=tot/M;
}
*/
inline void compareFibers(
	vector<FIBER> &bundle,
	EDGE &e,
	const int a, // fiber index in to the fiber bundle starts from 0.
	const int b, // fiber index
	float & d,
	bool & flagUseMax // if true use the maximum of the distances.
	)
{

	float tot = 0.0;
	float combinedShortestDist = 10000000.0;
	int M = bundle[a].points.size(); //size of fiber a
	int N = bundle[b].points.size(); //suze of fiber b
	int num = 0;
	float total = 0.0;
	int ClosestPointInA = 0;
	int ClosestPointInB = 0;


	for (int k = 0; k < M; k++)
	{	
		float D=0.0;
		int tempNode2_loc; 
		for (int l = 0; l < N; l++)
		{
			float tempD = sqrt(dist(bundle[a].points[k], bundle[b].points[l]));
			if (l == 0)
			{
				D = tempD;
				tempNode2_loc = l;
			}
			else
			{
				if (tempD < D){
					D = tempD;
					tempNode2_loc = l;
				}
			}
		}

		if (k == 0){ 
			combinedShortestDist = D;
			e.node1_loc = 0;
			e.node2_loc = tempNode2_loc;
		}
		else
		{
			if (D<combinedShortestDist)
			{
				combinedShortestDist = D;
				e.node1_loc = k;
				e.node2_loc = tempNode2_loc;
			}
		}
		// threshold
		if (D > 2.0)
		{
			tot = tot + D;
			num++;
		}
	}

	if (num == 0)
		d = 0.0;
	else
	{
		d = tot / num;
	}


	if (abs(d - 0.0) < 0.01)
	{
		d = 0.0;
		//cout << "total " << tot << ", num " << num << " \n";
	}

	//cout << "\nDistance between fiber " << a << " and fiber " << b << " is " << d << endl;

	/*
	int closestEndPt1, closestEndPt2;
	float closestDistance;

	int endpt1A = bundle[a].endPt1Index;
	int endpt2B = bundle[b].endPt2Index;
	int endpt2A = bundle[a].endPt2Index;
	int endpt1B = bundle[b].endPt1Index;
	float dA2B1 = sqrt(dist(bundle[a].points[endpt2A], bundle[b].points[endpt1B]));
	
	float cosAng1 =
		bundle[a].dir[endpt2A][0] * bundle[b].dir[endpt1B][0] +
		bundle[a].dir[endpt2A][1] * bundle[b].dir[endpt1B][1] +
		bundle[a].dir[endpt2A][2] * bundle[b].dir[endpt1B][2];

	float cosAng2 =
		bundle[a].dir[endpt1A][0] * bundle[b].dir[endpt2B][0] +
		bundle[a].dir[endpt1A][1] * bundle[b].dir[endpt2B][1] +
		bundle[a].dir[endpt1A][2] * bundle[b].dir[endpt2B][2];

	if (bundle[a].id == 77 || bundle[a].id == 550 || bundle[a].id == 573)
	{
		cout << "a " << bundle[a].id << " b " << bundle[b].id << endl;
		cout << "cosAngle1 " << cosAng1 << " cosAng2 " << cosAng2 << endl;
		cout << "distance changed from " << d << " to " << dA2B1 << endl;
	}

	if ((cosAng1 < -0.9) && (cosAng1 > -1.0))
	{
		if ((cosAng2 < -0.9) && (cosAng2 > -1.0))
		{
			
			d = dA2B1;
		}
	}

*/

	// Continity  computation
	
	//Compute distances. 
	
	int closestEndPt1, closestEndPt2;
	float closestDistance;

	int endpt1A = bundle[a].endPt1Index;
	int endpt2B = bundle[b].endPt2Index;
	int endpt2A = bundle[a].endPt2Index;
	int endpt1B = bundle[b].endPt1Index;
	float dA1B1 = sqrt(dist(bundle[a].points[endpt1A], bundle[b].points[endpt1B]));
	float dA1B2 = sqrt(dist(bundle[a].points[endpt1A], bundle[b].points[endpt2B]));
	float dA2B1 = sqrt(dist(bundle[a].points[endpt2A], bundle[b].points[endpt1B]));
	float dA2B2 = sqrt(dist(bundle[a].points[endpt2A], bundle[b].points[endpt2B]));

	map<float, pair<int, int>> endPtDistances;
	endPtDistances[dA1B1] = make_pair(endpt1A, endpt1B);
	endPtDistances[dA1B2] = make_pair(endpt1A, endpt2B);
	endPtDistances[dA2B1] = make_pair(endpt2A, endpt1B);
	endPtDistances[dA2B2] = make_pair(endpt2A, endpt2B);


	map<float, pair<int, int> >::iterator it = endPtDistances.begin();
	pair<int, int> p;
	p = it->second;

	//if (bundle[a].id == 77 || bundle[a].id == 550 || bundle[a].id == 573)
	//{
	//	cout << "a " << bundle[a].id << " b " << bundle[b].id << endl;
	//	cout << "closest distance " << it->first << " endpts " << p.first << " and " << p.second << endl;
	//}

	float cosAng1 =
		bundle[a].dir[p.first][0] * bundle[b].dir[p.second][0] +
		bundle[a].dir[p.first][1] * bundle[b].dir[p.second][1] +
		bundle[a].dir[p.first][2] * bundle[b].dir[p.second][2];

	if (cosAng1 < 0.0)
	{
		//compute the other angle. 
		//find the farthest and compute the angle
		map<float, pair<int, int> >::iterator itEnd = endPtDistances.end();
		--itEnd; //find the last element
		p = itEnd->second;
		//cout << "furhtest distance " << itEnd->first << " endpts " << p.first << " and " << p.second << endl;

		float cosAng2 =
			bundle[a].dir[p.first][0] * bundle[b].dir[p.second][0] +
			bundle[a].dir[p.first][1] * bundle[b].dir[p.second][1] +
			bundle[a].dir[p.first][2] * bundle[b].dir[p.second][2];


		if (cosAng2 < 0.0)
		{
			//cout << "Distance updated from " << d << " to " << it->first << endl;
			if ( it->first < d)
				d = it->first;
		}
	}
	
}
// OLD CODE
inline void compareFibersOld(
	vector<FIBER> &bundle,
	EDGE &e,
	const int a, // fiber index in to the fiber bundle starts from 0.
	const int b, //fiber index
	float & d, 
	bool & flagUseMax // if true use the maximum of the distances.
	)
{
	
	float tot = 0.0;
	float combinedShortestDist = 10000000.0;
	int M = bundle[a].points.size(); //size of fiber a
	int N = bundle[b].points.size(); //suze of fiber b
	int num = 0;
	float total = 0.0;
	int ClosestPointInA = 0;
	int ClosestPointInB = 0;

	/*
	for (int k = 0; k < M; k++)
	{
		float D = 999999999.0;
		int  TempClosestPointInB = 0;
		for (int l = 0; l < N; l++)
		{
			float tempD = sqrt(dist(bundle[a].points[k], bundle[b].points[l]));
			if (l == 0)
			{
				D = tempD;
			}
			if (tempD - D < 0.001)
			{
				D = tempD;
				TempClosestPointInB = l;
			}
		} // end for fiber B

		//float cosTheta = bundle[a].dir[k][0] * bundle[b].dir[TempClosestPointInB][0] +
		//	bundle[a].dir[k][1] * bundle[b].dir[TempClosestPointInB][1] +
		//	bundle[a].dir[k][2] * bundle[b].dir[TempClosestPointInB][2];
		//float sinTheta = sqrt(1 - cosTheta*cosTheta);
	

		//cout << " " << D<<"{"<< TempClosestPointInB<<"} ";
		//DEBUG
		total = total + D;

		//total = total + D + D*sinTheta;
		//total = total + D*sinTheta;

		if (D - combinedShortestDist < 0.001)
		{
			combinedShortestDist = D;
			ClosestPointInA = k;
			ClosestPointInB = TempClosestPointInB;
		}
	} // for end 
	d = total / M;
	e.node1_loc = ClosestPointInA;
	e.node2_loc = ClosestPointInB;
	*/

	//DEBUG
	//cout << "\nDistance between fiber " << a << " and fiber " << b << " is " << d << endl;
	//float tempDistanceBetweenEndPoints = sqrt(dist(bundle[a].points[bundle[a].endPt1Index],
	//	bundle[b].points[bundle[b].endPt2Index]));
	//cout << "Distance between endpt " << bundle[a].endPt1Index << " and " << bundle[b].endPt2Index << " is " << tempDistanceBetweenEndPoints << endl;
	


	/*
	// Continity  computation
	//Compute distances. 
	int closestEndPt1, closestEndPt2;
	float closestDistance;
	
	int endpt1A = bundle[a].endPt1Index;
	int endpt2B = bundle[b].endPt2Index;
	int endpt2A = bundle[a].endPt2Index;
	int endpt1B = bundle[b].endPt1Index;
	float dA1B1 = sqrt(dist(bundle[a].points[endpt1A], bundle[b].points[endpt1B]));
	float dA1B2 = sqrt(dist(bundle[a].points[endpt1A], bundle[b].points[endpt2B]));
	float dA2B1 = sqrt(dist(bundle[a].points[endpt2A], bundle[b].points[endpt1B]));
	float dA2B2 = sqrt(dist(bundle[a].points[endpt2A], bundle[b].points[endpt2B]));
	
	map<float, pair<int, int>> endPtDistances;
	endPtDistances[dA1B1] = make_pair(endpt1A, endpt1B);
	endPtDistances[dA1B2] = make_pair(endpt1A, endpt2B);
	endPtDistances[dA2B1] = make_pair(endpt2A, endpt1B);
	endPtDistances[dA2B2] = make_pair(endpt2A, endpt2B);


	map<float, pair<int, int> >::iterator it = endPtDistances.begin();
	pair<int, int> p;
	p = it->second;
	//cout << "closest distance " << it->first << " endpts " << p.first << " and " << p.second << endl;

	float cosAng1 = 
		bundle[a].dir[p.first][0] * bundle[b].dir[p.second][0] +
		bundle[a].dir[p.first][1] * bundle[b].dir[p.second][1] +
		bundle[a].dir[p.first][2] * bundle[b].dir[p.second][2];
	float sinAng1 = sqrt(1 - cosAng1*cosAng1);
	if (sinAng1 < 0.5)
	{
		//compute the other angle. 
		//find the farthest and compute the angle
		map<float, pair<int, int> >::iterator itEnd = endPtDistances.end();
		--itEnd; //find the last element
		p = itEnd->second;
		//cout << "furhtest distance " << it->first << " endpts " << p.first << " and " << p.second << endl;

		float cosAng2 =
			bundle[a].dir[p.first][0] * bundle[b].dir[p.second][0] +
			bundle[a].dir[p.first][1] * bundle[b].dir[p.second][1] +
			bundle[a].dir[p.first][2] * bundle[b].dir[p.second][2];
		float sinAng2 = sqrt(1 - cosAng2*cosAng2);


		if (sinAng2 < 0.5)
		{
			d = it->first;
		}
	}
	*/
	////////////////////////////


	/*
	float cosAng1 = bundle[a].dir[endpt1A][0] * bundle[b].dir[endpt2B][0] +
		bundle[a].dir[endpt1A][1] * bundle[b].dir[endpt2B][1] +
		bundle[a].dir[endpt1A][2] * bundle[b].dir[endpt2B][2];
	float sinAng1 = sqrt(1 - cosAng1*cosAng1);

	cout << "cos angle " << endpt1A << "  and " << endpt2B << " is " << cosAng1
		<< " sinAngle is "<< sinAng1 << endl;
	if (sinAng1 < 0.5)
	{
		float cosAng2 = bundle[a].dir[endpt2A][0] * bundle[b].dir[endpt1B][0] +
			bundle[a].dir[endpt2A][1] * bundle[b].dir[endpt1B][1] +
			bundle[a].dir[endpt2A][2] * bundle[b].dir[endpt1B][2];

		float sinAng2 = sqrt(1 - cosAng2*cosAng2);
		cout << "sinAngle between " << endpt2A << " and " << endpt1B << " is " << sinAng2 << endl;

		if (sinAng2 < 0.5)
		{
			float d1 = sqrt(dist(bundle[a].points[endpt1A],
				bundle[b].points[endpt2B]));
			float d2 = sqrt(dist(bundle[a].points[endpt2A],
				bundle[b].points[endpt1B]));

			if (d1 < d2)
				d = d1;
			else
				d = d2;
			cout << "Distance changed to "<< d << endl;
		}
		else
		{
			//angles dont match 
			flagUseMax = true;
		}
	}
	else
	{
		//angles dont match 
		flagUseMax = true;
	}
	*/





	
	//Original distance computataion.

	for (int k = 0; k < M; k++)
	{
		float D = 1000000.0;
		for (int l = 0; l < N; l++)
		{
			float tempD = sqrt(dist(bundle[a].points[k], bundle[b].points[l]));
			if (tempD < D){
				D = tempD;
				if (tempD < combinedShortestDist)
				{
					e.node1_loc = k;
					e.node2_loc = l;
					combinedShortestDist = tempD;
				}
			}
		}
		// threshold
		if (D > 5.0){
			tot = tot + D;
			num++;
		}
	}

	if (num == 0)
		d = 0.0;
	else
	{
		d = tot / num;
	}


	if (abs(d - 0.0) < 0.01){
		d = 0.0;
		cout << "total " << tot << ", num " << num << " \n";
	}
		
	cout << "\nDistance between fiber " << a << " and fiber " << b << " is " << d << endl;



}

//Advanced distance computations // no parallel
//a,b are the fibers being compared
//d Min mean distance
inline void compareAdvanced(
	vector<FIBER> &bundle,
	EDGE &e,
	const int a,
	const int b,
	float & d)
{
	float mean_min_distAB = 0.0, mean_min_distBA = 0.0;
	bool  flagUseMax = false;
	compareFibers(bundle, e, b, a, mean_min_distBA, flagUseMax);

	compareFibers(bundle, e, a, b, mean_min_distAB, flagUseMax);

	//DEBUG return the minumum
	//d = (mean_min_distAB + mean_min_distBA) / 2.0;
	// INSTEAD OF THE MIN returning the average

	
	if (mean_min_distAB < mean_min_distBA)
	{
		if (!flagUseMax) // if flagUseMax is false then return the minimum
		{
			d = mean_min_distAB;
		}
		else
		{
			d = mean_min_distBA; // return the maximum
		}
	}
	else
	{
		if (!flagUseMax){ d = mean_min_distBA; }
		else { d = mean_min_distAB; } // return the max 
		
	}

	//cout << "mean distance AB " << mean_min_distAB << " mean distance BA " << mean_min_distBA << "  d " << d << endl;
	// cout << "distance returned " << d << " using max ? "<< int(flagUseMax) << endl; 
	//cout << " AB: " << mean_min_distAB <<", BA: " << mean_min_distBA << ", d: " << d <<endl;

}


//
void correctOrientationEndPointsFiber
( vector<FIBER> & bundle,
const int index // fiber index
)
{
	const int fiberSize = bundle[index].points.size();
	int endpt1A, endpt1B, endpt2A, endpt2B, endpt1APrevious;

	if (fiberSize < 2)
		return;
	FIBER f = bundle[index];

	if (f.endPt1Index != 0)
	{
		endpt1A = f.endPt1Index;
		endpt1APrevious = f.endPt1Index - 1; 

		vector <float> vecAB(3, 0.0);
		for (int d = 0; d < 3; d++)
		{
			vecAB[d] = f.points[endpt1APrevious][d] - f.points[endpt1A][d];
		}

		// magnitude compuatation
		float magVecAB = 0.0;
		for (int d = 0; d < 3; d++)
		{
			magVecAB = magVecAB +  vecAB[d] * vecAB[d];
		}

		magVecAB = sqrt(magVecAB);

		if (abs(magVecAB - 0.0) < 0.0001)
			return;

		float dotPdt = 0.0;
		for (int d = 0; d < 3; d++)
		{
			dotPdt = dotPdt + (vecAB[d]/magVecAB) * f.dir[endpt1A][d];
		}

		if (dotPdt < 0.0)
		{
			for (int d = 0; d < 3; d++)
			{
				bundle[index].dir[endpt1A][d] = bundle[index].dir[endpt1A][d] * (-1.0);
			}
			//cout << "id " << f.id << " dotpdt negative" << endl;
		}
	}

	//change endpt2
	endpt2A = f.endPt2Index;
	int endpt2APrevious = endpt2A - 1;
	vector <float> vecAB(3, 0.0);
	for (int d = 0; d < 3; d++)
	{
		vecAB[d] = f.points[endpt2APrevious][d] - f.points[endpt2A][d];
	}
	// magnitude compuatation
	float magVecAB = 0.0;
	for (int d = 0; d < 3; d++)
	{
		magVecAB = magVecAB + vecAB[d] * vecAB[d];
	}

	magVecAB = sqrt(magVecAB);

	if (abs(magVecAB - 0.0) < 0.0001)
		return;

	float dotPdt = 0.0;
	for (int d = 0; d < 3; d++)
	{
		dotPdt = dotPdt + vecAB[d]/magVecAB * f.dir[endpt2A][d];
	}
	if (dotPdt < 0.0)
	{
		for (int d = 0; d < 3; d++)
		{
			bundle[index].dir[endpt2A][d] = bundle[index].dir[endpt2A][d] * (-1.0);
		}
	}
}
void correctOrientationEndPointsBundle(vector<FIBER> &bundle)
{
	const int numberOfFibers = bundle.size();
	for (int f = 0; f < numberOfFibers; f++)
	{
		correctOrientationEndPointsFiber(bundle, f);
	}
}

// Average orientation at the end points fiber
void averageOrientationEndPointsFiber(
	vector<FIBER> &bundle, const int index)
{
	const int fiberSize = bundle[index].points.size();
	const FIBER f = bundle[index];

	if (fiberSize < 2){
		cerr << "Error in averageOrientationEndPointsFiber" << endl;
		exit(0);
	}
	if (f.endPt1Index != 0)
	{
		vector<float> tempDir(3, 0.0);
		// orientaion at 0
		for (int d = 0; d < 3; d++)
		{
			tempDir[d] = tempDir[d] + f.dir[f.endPt1Index][d];
		}
		// Orientation at 1
		for (int d = 0; d < 3; d++)
		{
			tempDir[d] = tempDir[d] + f.dir[f.endPt1Index - 1][d];
		}
		//set endpoint1 
		for (int d = 0; d < 3; d++)
		{
			bundle[index].dir[f.endPt1Index][d] = tempDir[d] / 2.0;
		}
		//DEBUG 
		if (f.id == 77 || f.id == 550 || f.id == 573)
		{
			cout << "id " << f.id << endl;
			cout << f.endPt1Index << " , " << f.endPt1Index - 1 << " || " 
				<< bundle[index].dir[f.endPt1Index][0] << " " << bundle[index].dir[f.endPt1Index][1]
				<< " " << bundle[index].dir[f.endPt1Index][2] << endl;
		}	
	}
	
	if (f.endPt1Index == 0)
	{
		cout << "id " << f.id << endl;
		vector<float> tempDir(3, 0.0);
		for (int d = 0; d < 3; d++)
		{
			tempDir[d] = tempDir[d] + f.dir[f.endPt1Index][d];
			// Note the f.endpt1Index + 1 has the same direction just in the opposite direction
			// so we use the +2
			tempDir[d] = tempDir[d] + f.dir[f.endPt1Index + 2][d]*(-1.0);
		}
		for (int d = 0; d < 3; d++)
		{
			bundle[index].dir[f.endPt1Index][d] = tempDir[d] / 2.0;
		}
	}
	// Set the endpt2
	vector<float> tempDir(3, 0.0);
	for (int d = 0; d < 3; d++)
	{
		tempDir[d] = tempDir[d] + f.dir[f.endPt2Index][d];
		tempDir[d] = tempDir[d] + f.dir[f.endPt2Index - 1][d];
	}
	for (int d = 0; d < 3; d++)
	{
		bundle[index].dir[f.endPt2Index][d] = tempDir[d] / 2.0;
	}
}
//AverageOrientationEndPointsBundle
void averageOrientationEndPointsBundle(vector<FIBER> &bundle)
{
	const int numberOfFibers = bundle.size();
	for (int f = 0; f < numberOfFibers; f++)
	{
		averageOrientationEndPointsFiber(bundle, f);
	}
}
