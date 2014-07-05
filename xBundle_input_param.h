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
#include <vector>
using namespace std;
#include "xBundle_fiber_info.h"
class INPUT_PARAMS {
public:
  string inputFileName;
  bool inputFileUI;
  bool inputTh;
  bool computeHessianBasedComp;
  bool sparseCoding;
  float sigma;
  int thresh;
  int eigenVecNum;
  bool degugMode;
  bool writeDataMode;
  bool aniso_diff;
  bool computeKmeans;
  int radius;
  int linInd;
  int Dimension;
  float vness;
  string dataFileName;
  bool reliableHess;
  bool byParts;
  bool color;
  bool trackPoint;
  int boxSize;
  bool colorTrack;
  int colorTrackID;
  bool colorVol;
  float minFiberLength;
  string graphFileName;
  string graph2FileName;
  string bundleInfoFname;
  string bundleVolFname;
  
  //config
  bool clusterSpecific;
  bool colorSpecificFibers;
  vector <int> indicesColorSpecificFibers;
  vector <int> clusterCols;
  vector <int> indicesFiberInfo;
  vector <int> startLocations;
  vector <int> voxelClusters;


  //cylinder 
  float cylLength;
  float cylRadius;
  



  INPUT_PARAMS()
  {
    inputFileUI = false;
    writeDataMode = false;
    inputTh = false;
    computeHessianBasedComp = false;
    sparseCoding = false;
    eigenVecNum=0;
    sigma=3.0f;
    radius=7;
    Dimension =3;
    dataFileName = "data.csv";
    degugMode = false;
	thresh = 10;
	aniso_diff=true;
	vness = 0.01;
	computeKmeans=false;
	reliableHess = false;
	byParts = false;
	color = false;
	boxSize = 40;
	trackPoint = false;
	colorTrack = false;
	colorTrackID = -1;
	
	minFiberLength = 50.0;
	
	colorVol=false;
	clusterSpecific = false;
	colorSpecificFibers = false; 
	
	//Cylinder
	cylLength = 10.0;
	cylRadius = 3.0;
	
	//outputs D:\\ABhattacharya\\Internship\\output
	/*
	graphFileName = "C:\\Users\\p41123\\Documents\\internship\\output\\graph.csv";
	graph2FileName = "C:\\Users\\p41123\\Documents\\internship\\output\\graphClus.csv";
	bundleInfoFname="C:\\Users\\p41123\\Documents\\internship\\output\\bundle.csv";
	bundleVolFname="C:\\Users\\p41123\\Documents\\internship\\output\\bundleVol.mhd";
	*/
	graphFileName = "D:\\ABhattacharya\\Internship\\output\\graph.csv";
	graph2FileName = "D:\\ABhattacharya\\Internship\\output\\graphClus.csv";
	bundleInfoFname="D:\\ABhattacharya\\Internship\\output\\bundle.csv";
	bundleVolFname="D:\\ABhattacharya\\Internship\\output\\bundleVol.mhd";

 }
};
//K:\Internship\code\itkvtkTest\WindowsBuilds\bundleExBuild\Debug>xBundle.exe -i test -byParts -aniso_diff  -vness 0.8 -s 3
#endif

