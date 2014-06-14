#include <iostream>
#include <vector>
#include "xBundle_fiber_info.h"
using namespace std;


void FIBER::pushPoint(vector <int> & p)
{
	points.push_back(p);
}
void FIBER::pushPerpDist(float & d)
{
	perDist.push_back(d);
}

void FIBER::pushDir(vector <float> &f)
{
	dir.push_back(f);
}

void FIBER::colorFiber( )
{
	cout <<"in function color single fiber "<< id << endl;
}

void FIBER::fiberInfo()
{
	cout <<"ID: "<<id;
	cout <<"\tlength "<<length;
	cout <<"\tpoints "<<points.size();
	cout <<"\tallpoints "<<allPoints.size();
	cout << "\tendpts " << endPt1Index << " , " << endPt2Index << endl;
	cout <<endl;
	for (int i=0;i<points.size();i++)
	{
		cout <<i<<" : "<< points[i][0]<<","<<points[i][1]<<","<<points[i][2];
		cout <<" dir "<< dir[i][0]<<","<<dir[i][1]<<","<<dir[i][2]<<endl;
	}
}

