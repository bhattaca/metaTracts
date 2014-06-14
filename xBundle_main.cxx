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

#include <thread>
#include <mutex>
#include <iostream>
#include <fstream>
using namespace std;

// boost inlcudes
#include "boost/tokenizer.hpp"
#include <boost/algorithm/string.hpp>
using namespace boost;

#include "xBundle_input_param.h"
#include "xBundle_hessian.h"
#include "xBundle_sparseCoding.h"



vector<EDGE> graph1;

void hello(){

	std::cout << "Hello from thread " << std::this_thread::get_id() << std::endl;
}
void compute_distances_threaded(const INPUT_PARAMS  * input, vector<FIBER> &bundle);


//Function to read the config  file.
//Return True if there are any errors in reading the file.
void readConfigFile(INPUT_PARAMS * input, const string configFileName)
{
	cout << " In function read config file" << endl;
	cout << " Config file name [" << configFileName << "]" << endl;
	ifstream in(configFileName.c_str());
	if (!in.is_open()){ throw "FNAME_ERR"; }
	string line;
	vector <string>  configInfo;

	string s = "exec script1 \"script argument number one\"";
	string separator1("");//dont let quoted arguments escape themselves
	string separator2(",");//split on spaces
	string separator3("\"\'");//let it have quoted arguments

	escaped_list_separator<char> els(separator1, separator2, separator3);

	typedef tokenizer< escaped_list_separator<char> > Tokenizer;


	while (getline(in, line))
	{
		Tokenizer tokens(line, els);
		configInfo.assign(tokens.begin(), tokens.end());
		cerr << "Config Option: " << configInfo[0] << endl;
		if (configInfo[0] == "CLUS_SPECIFIC")
		{
			input->clusterSpecific = true;
			for (int i = 1; i < configInfo.size(); i++)
			{
				input->clusterCols.push_back(atoi(configInfo[i].c_str()));
			}
		}
		else if (configInfo[0] == "FIBER_SPECIFIC")
		{
			input->colorSpecificFibers = true;
			for (int i = 1; i < configInfo.size(); i++)
			{
				input->indicesColorSpecificFibers.push_back(atoi(configInfo[i].c_str()));
			}
		}
		else if (configInfo[0] == "FIBER_INFO")
		{
			for (int i = 1; i < configInfo.size(); i++){
				input->indicesFiberInfo.push_back(atoi(configInfo[i].c_str()));
			}
		}
		else if (configInfo[0] == "START_LOCS")
		{
			if ((configInfo.size() - 1) % 3 != 0)
			{
				throw "Start locations must be multiples of 3.";
			}
			else
			for (int i = 1; i < configInfo.size(); i++)
			{
				input->startLocations.push_back(atoi(configInfo[i].c_str()));
			}
		}
		else if (configInfo[0] == "INPUTFILE_FULLPATH")
		{
			if (configInfo.size() != 2)
				throw "NO FILE NAME FILE NAME FORMAT ERROR";
			else
			{
				input->inputFileName = configInfo[1];
				trim(input->inputFileName);
				input->inputFileUI = true;
				cerr << "InputFileName:" << input->inputFileName << endl;
			}
		}
		else
		{
			cerr << "NONE Of the options match" << endl;
			throw configInfo[0];
		}
			
	}

	//throw "exception";
	//throw (10);
}


void usage(char* argv[]){
	cout << "Usage: " << argv[0] << " Options" << endl;
	cout << "Options :" << endl;
	cout << "\t -i INPUT_FILENAME // required" << endl;
	cout << "\t -thres <int> //required" << endl;
	cout << "\t -hessian. perform hessian based computation" << endl;
	cout << "\t -s <float> default 3.0" << endl;
	cout << "\t -eig <int> default 2" << endl;
	cout << "\t -config : Reads the config file." << endl;
	cout << "\t-byParts do in regions" << endl;
	cout << "\t-colorTrack -1 color all tracks , else color the track id" << endl;

}

// read the parameters from the users
int readParameters(int argc, char* argv[], INPUT_PARAMS * input)
{
	int ii = 1;
	while (ii < argc) {
		cout << "OPTION " << argv[ii] << endl;
		if (strcmp(argv[ii], "-hessian") == 0)
		{
			cout << "Hessian option" << endl;
			input->computeHessianBasedComp = true;
		}
		else if (strcmp(argv[ii], "-color") == 0)
		{
			input->color = true;
		}
		else if (strcmp(argv[ii], "-reliableHess") == 0)
		{
			input->reliableHess = true;
			cout << "Computing reliable hessians " << endl;
		}
		//else if (strcmp(argv[ii], "-byParts") == 0)
		//{
		//	cout <<"compute byparts"<<endl;
		//	input->byParts=true;
		//}
		else if (strcmp(argv[ii], "-track") == 0)
		{
			cout << "Setting option to track points." << endl;
			input->trackPoint = true;
		}
		else if (strcmp(argv[ii], "-sparseCoding") == 0)
		{
			cout << "SparseCoding !!!!" << endl;
			input->sparseCoding = true;
		}
		/*else if (strcmp(argv[ii], "-i") == 0)
		{
			ii++;
			cout << "FileName " << argv[ii] << endl;
			input->inputFileUI = true;
			input->inputFileName = argv[ii];
		}*/
		else if (strcmp(argv[ii], "-colorTrack") == 0)
		{
			ii++;
			input->colorTrack = true;
			input->trackPoint = true;
			input->colorTrackID = atoi(argv[ii]);
		}
		else if (strcmp(argv[ii], "-aniso_diff") == 0)
		{
			cout << "Anisotropic diffusion " << endl;
			input->aniso_diff = true;
		}
		else if (strcmp(argv[ii], "-config") == 0)
		{
			cout << "Reading the config file.\n";

			try
			{
				ii++;
				string configFileName = argv[ii];
				readConfigFile(input, configFileName);
				cerr << "Reading complete" << endl;
			}
			catch (char const * pch)
			{
				cerr << "Threw this exception: " << pch << endl;
				if (strcmp(pch, "FNAME_ERR") == 0)
				{
					cerr << "File Name error." << endl;
				}
			}
			catch (...)
			{
				cerr << "Unknow Exception. Execution halted, cannot read config file. " << endl;
				exit(0);
			}

		}
		else if (strcmp(argv[ii], "-kmeans") == 0)
		{
			cerr << "Kmeans flag" << endl;
			input->computeKmeans = true;
		}
		else if (strcmp(argv[ii], "-colVol") == 0)
		{
			cerr << "Color vol" << endl;
			input->colorVol = true;
		}
		else if (strcmp(argv[ii], "-s") == 0)
		{
			ii++;
			cout << "sigma " << argv[ii] << endl;
			input->sigma = float(atof(argv[ii]));
		}
		else if (strcmp(argv[ii], "-boxSize") == 0)
		{
			ii++;
			input->boxSize = int(atoi(argv[ii]));
			cout << "boxSize " << input->boxSize << endl;
		}
		else if (strcmp(argv[ii], "-thres") == 0)
		{
			ii++;
			input->thresh = atoi(argv[ii]);
			cout << "Thres " << input->thresh << endl;
			input->inputTh = true;
		}
		else if (strcmp(argv[ii], "-vness") == 0)
		{
			ii++;
			input->vness = atof(argv[ii]);
			cout << "Vness " << input->vness << endl;
		}
		else if (strcmp(argv[ii], "-eig") == 0)
		{
			ii++;
			input->eigenVecNum = atoi(argv[ii]);
			cout << "eig vec num " << input->eigenVecNum << endl;
		}
		else if (strcmp(argv[ii], "-radius") == 0)
		{
			ii++;
			input->radius = atoi(argv[ii]);
			input->linInd = (2 * input->radius + 1)*(2 * input->radius + 1)*input->Dimension;
		}
		else if (strcmp(argv[ii], "-writeData") == 0)
		{
			ii++;
			input->dataFileName = argv[ii];
			input->writeDataMode = true;
			cout << "Write data " << endl;
		}
		else if (strcmp(argv[ii], "-debug") == 0)
		{
			input->degugMode = true;
		}
		else
		{
			cout << "Option[" << argv[ii] << "]does not match anything " << endl;
			usage(argv);
			return 10;
		}
		ii++;
	}

	return 0;
}



// Compute the required functions on the data
int performComputations(INPUT_PARAMS * input)
{
	cout << "In Perform Computation." << endl;
	//fiber bundle info
	vector<FIBER> bundle;
	vector<EDGE> graph;

	if (input->computeHessianBasedComp)
	{

		//If there is no input  file then stop computation.
		if (!input->inputFileUI){
			cerr << "Input file could not be read or is corrupted check config file option INPUTFILE_FULLPATH" << endl;
			exit(0);
		}

		//perform hessian based computation 
		hessian_based_computation(input, bundle);

		// printinfo 
		if (!input->indicesFiberInfo.empty())
		{
			int vecLen = input->indicesFiberInfo.size();
			for (int i = 0; i < vecLen; i++)
			{
				bundle[input->indicesFiberInfo[i] - 1].fiberInfo();
			}
		}

		if (input->trackPoint)
		{
			cerr << "original number of fibers " << bundle.size() << endl;
			
			cerr << "Compute distances" << endl;


			if (!input->computeKmeans && bundle.size() > 0){
				compute_distances(input, bundle, graph);
				//compute_distances_threaded(input, bundle);
				//cout <<"write to graph "<<graph1.size()<<endl;
				//write_graph(graph1, input->graph2FileName);
				write_graph(graph, input->graphFileName);
			}
		}
		cerr << "fiber bundles detected!" << endl;

		//for (int m=0;m<bundle.size();m++)
		//{
		//	cout <<"fiber number "<< bundle[m].id <<endl;
		//	cout <<"length of the fiber " << bundle[m].length <<endl;


		//	if (bundle[m].length > 20.0){

		//		revisedBundle.push_back(bundle[m]);

		//		cout <<"points on the fiber "<<endl;
		//		for(int o=0;o<bundle[m].points.size();o++)
		//		{
		//			cout <<"["<<o<<"]"<<bundle[m].points[o][0]<<","<<bundle[m].points[o][1]
		//			<<","<<bundle[m].points[o][2]<<" dist " << bundle[m].perDist[o]<<" dir "
		//				<<bundle[m].dir[o][0]<<","
		//				<<bundle[m].dir[o][1]<<","
		//				<<bundle[m].dir[o][2]<<","
		//				<<endl;
		//		}
		//	}
		//}

		//cout <<" bundle size " << bundle.size() << endl;
		//cout <<"revised bundle " <<revisedBundle.size() <<endl;
		//for (int m=0;m<revisedBundle.size();m++)
		//{
		//	cout <<"fiber number "<< revisedBundle[m].id <<endl;
		//	cout <<"length of the fiber " << revisedBundle[m].length <<endl;
		//	cout <<"points on the fiber "<<endl;
		//	for(int o=0;o<revisedBundle[m].points.size();o++)
		//	{
		//		cout <<"["<<o<<"]"<<revisedBundle[m].points[o][0]<<","<<revisedBundle[m].points[o][1]
		//		<<","<<revisedBundle[m].points[o][2]<<" dist " << revisedBundle[m].perDist[o]<<" dir "
		//			<<revisedBundle[m].dir[o][0]<<","
		//			<<revisedBundle[m].dir[o][1]<<","
		//			<<revisedBundle[m].dir[o][2]<<","
		//			<<endl;
		//	}

		//}
	}

	else if (input->sparseCoding)
	{
		cout << "sparse coding" << endl;
		input->linInd = (2 * input->radius + 1)*(2 * input->radius + 1)*input->Dimension;
		cout << "call sparse coding, radius set to " << input->radius << " single vector is of size " << (2 * input->radius + 1)*(2 * input->radius + 1) * 3 << endl;
		sparseCoding_computations(input);
	}
	/*
	// aniso diffusion now being done inside the hessian computaion
	else if (input->aniso_diff)
	{

	input->linInd=(2*input->radius+1)*(2*input->radius+1)*(2*input->radius+1);
	cout <<" in here call sparse coding, radius set to " << input->radius << " single vector is of size " << (2*input->radius+1)*(2*input->radius+1)*(2*input->radius+1)<< endl;
	// also outputs reults for sparse coding
	xbundle_aniso_diff (input);
	}
	*/
	//else if ( input->byParts)
	//{
	//	cout <<"call by parts"<<endl;
	//	bypart(input);
	//}
	else if (input->color)
	{
		cout << "Color " << endl;
		computeColor(input);
	}
	else if (input->computeKmeans)
	{
		/*		cout <<"compute k means clustering"<< endl;
		input->computeKmeans= true;
		computeKmeansOutput(input);
		*/
	}
	return 0;
}

//Main function
int main(int argc, char* argv[])
{
	// read input parameters
	INPUT_PARAMS * input = new INPUT_PARAMS();

	int readParamRes = readParameters(argc, argv, input);
	if (readParamRes != 0)
		cout << "function readParameters error " << readParamRes << endl;

	int computeRes = performComputations(input);
	cout << " computeRes " << computeRes << endl;
	return 0;
}


inline float dist(vector <int> v1, vector<int> v2)
{
	float tempD = 0.0;
	for (int d = 0; d < 3; d++)
	{
		tempD = tempD + (v1[d] - v2[d])*(v1[d] - v2[d]);
	}

	return tempD;
}

/*
inline void compareFibers (vector<FIBER> &bundle, EDGE &e, const int a, const int b, float & d)
{

float tot=0.0;
float combinedShortestDist=10000000.0;
int M=bundle[a].points.size();
int N=bundle[b].points.size();
int num=0;
for (int k=0;k<M;k++)
{
float D=1000000.0;
for (int l=0;l<N;l++)
{
float tempD=sqrt(dist (bundle[a].points[k], bundle[b].points[l]));

if (tempD<D){
D=tempD;
if (tempD<combinedShortestDist)
{
e.node1_loc=k;
e.node2_loc=l;
combinedShortestDist=tempD;
}
}
}
// threshold
if (D > 5.0){
tot=tot+D;
num++;
}
}
d=tot/num;

if (abs(d-0.0) < 0.01)
cerr <<" distance is 0 "<<endl;
}
//Advanced distance computations
inline void compareAdvanced (vector<FIBER> &bundle, EDGE &e, const int a, const int b, float & d)
{
float mean_min_distAB=0.0,mean_min_distBA=0.0;
compareFibers (bundle, e, b, a, mean_min_distBA);

compareFibers (bundle, e, a, b, mean_min_distAB);

//if (mean_min_distAB < mean_min_distBA)
//debug
if (mean_min_distAB > mean_min_distBA)
d=mean_min_distAB;
else
d=mean_min_distBA;
}
*/

std::mutex g_pages_mutex;

//void func1(vector<FIBER> &bundle, const int m)
//{
//	for (unsigned int n=m+1; n < bundle.size();n++)
//	{
//
//		float d=0.0;
//		EDGE a;
//		compareAdvanced(bundle,a, m, n, d);
//		a.node1=bundle[m].id;
//		a.node2=bundle[n].id;
//		float dotPdt = 0.0f;
//		for (int i=0;i<3;i++)
//		{
//			dotPdt=dotPdt+(bundle[m].dir[a.node1_loc][i]  * bundle[n].dir[a.node2_loc][i]);
//		}
//
//		if (dotPdt < 0.0)
//		{
//			dotPdt=0.0;
//			for (int i=0;i<3;i++)
//			{
//				dotPdt=dotPdt+(bundle[m].dir[a.node1_loc][i]  * (-1.0)*bundle[n].dir[a.node2_loc][i]);
//			}
//		}
//		a.edgeW=dotPdt;
//		a.edgeW2=exp(-d/500.0);
//
//
//		g_pages_mutex.lock();
//		graph1.push_back(a);
//		g_pages_mutex.unlock();
//
//	}
//}
//void compute_distances_threaded(const INPUT_PARAMS  * input, vector<FIBER> &bundle)
//{
//
//	unsigned int M = bundle.size();
//	std::vector<std::thread> threads;
//	std::mutex mx;
//	for (unsigned int m=0; m < M-1;m++)
//	{
//		if (m%100==0)
//			cout <<"m "<<m<<endl;
//		threads.push_back(std::thread(func1, bundle, m));
//	}
//	cout <<"waiting for threads to finish";
//	for(auto& thread : threads){
//		thread.join();
//	}
//}


/*
void func2(vector<FIBER> &bundle, const int Mstart, const int Mend)
{
for (unsigned int m=Mstart; m < Mend;m++){
for (unsigned int n=m+1; n < bundle.size();n++)
{

float d=0.0;
EDGE a;
compareAdvanced(bundle,a, m, n, d);
a.node1=bundle[m].id;
a.node2=bundle[n].id;
double dotPdt = 0.0f;
for (int i=0;i<3;i++)
{
dotPdt=dotPdt+(bundle[m].dir[a.node1_loc][i]  * bundle[n].dir[a.node2_loc][i]);
}

if (dotPdt < 0.0)
{
dotPdt=0.0;
for (int i=0;i<3;i++)
{
dotPdt=dotPdt+(bundle[m].dir[a.node1_loc][i]  * (-1.0)*bundle[n].dir[a.node2_loc][i]);
}
}
a.edgeW=dotPdt;
a.edgeW2=d;
//a.edgeW2=exp(-d/500.0);


g_pages_mutex.lock();

graph1.push_back(a);

g_pages_mutex.unlock();

}

}
}
void compute_distances_threaded(const INPUT_PARAMS  * input, vector<FIBER> &bundle)
{
cout <<"debug "<<endl;
unsigned int M = bundle.size();
std::vector<std::thread> threads;
std::mutex mx;
int x=0;
int numPerT=300;
if (bundle.size() < numPerT)
x=bundle.size();
else
x=numPerT;
int numThreads = int(M/x);
cout <<"number of threads "<<numThreads <<endl;
int t=0,m=0,n=0;
while (n < M)
{
m=t*x;
n=t*x+x-1;
if (n > M-m)
n=t*x+M-m;
t++;
cout <<" start " <<m <<" end "<<n<<endl;
threads.push_back(std::thread(func2, bundle, m,n+1));
}
cout <<"waiting for threads to finish\n";
for(auto& thread : threads){
thread.join();
}
}
*/