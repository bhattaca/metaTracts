#ifndef XBUNDLE_HESSIAN_H_INCLUDED
#define XBUNDLE_HESSIAN_H_INCLUDED


#include "xBundle_input_param.h"
#include "xBundle_fiber_info.h"
void hessian_based_computation(const INPUT_PARAMS * input, vector<FIBER> &bundle);
void hessian_based_computation(INPUT_PARAMS * input, vector<FIBER> &bundle);
void xbundle_aniso_diff (const INPUT_PARAMS * input);
void computeKmeansOutput(const INPUT_PARAMS * input);
//void bypart(const INPUT_PARAMS *input);
void remove_small_fibers(const INPUT_PARAMS  * input, vector<FIBER> &bundle, vector<FIBER> &revisedBundle);
void compute_distances (const INPUT_PARAMS  * input, vector<FIBER> &bundle, vector<EDGE> &graph);
void write_graph( vector<EDGE> &graph,const string graphFileName);
void write_fibers (vector<FIBER> &bundle, const string bundleInfoFname);
void computeColor ( const INPUT_PARAMS * input);

const std::string currentDateTime();
#endif // XBUNDLE_HESSIAN_H_INCLUDED
