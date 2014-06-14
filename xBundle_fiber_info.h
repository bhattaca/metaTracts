#ifndef BundleEx_xBundle_fiber_info_h
#define BundleEx_xBundle_fiber_info_h
#include <iostream>
#include <vector>  
using namespace std;

class FIBER{
public :
	
	unsigned int id;
	float length;
	int endPt1Index;
	int endPt2Index;
	vector <float> seedPtDir;
	vector <vector<int> > points;
	vector <vector<int> > allPoints;
	vector <vector<float> > dir;
	vector <float>  perDist;
	void pushPoint(vector<int> &); // set of points in the fiber
	void pushPerpDist(float &); // set of perpendicular distances
	void pushDir(vector<float> &); // set of directions
	void colorFiber();
	void fiberInfo();
	FIBER()
	{
		id=0;
		length=0;
	}
};

//GRAPH 
class EDGE
{
public:
	int node1;
	int node1_loc;
	int node2_loc;
	int node2;
	double edgeW;// dot pdt
	double edgeW2;//distance
};




#endif // XBUNDLE_HESSIAN_H_INCLUDED