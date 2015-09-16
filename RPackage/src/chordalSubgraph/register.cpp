#include <Rcpp.h>
#include "isChordalIterative.h"
#ifdef _MSC_VER
	#undef RcppExport
	#define RcppExport extern "C" __declspec(dllexport)
#endif
extern "C" const char* package_name = "chordalSubgraph";
R_CallMethodDef callMethods[] = 
{
	{"isChordalIterative_igraph", (DL_FUNC)&isChordalIterative_igraph, 1},
	{"isChordalIterative_graphNEL", (DL_FUNC)&isChordalIterative_graphNEL, 1},
	{"isChordalIterative_graphAM", (DL_FUNC)&isChordalIterative_graphAM, 1},
	{NULL, NULL, 0}
};
extern "C" void R_init_Rcpp(DllInfo* info);
RcppExport void R_init_chordalSubgraph(DllInfo *info)
{
	R_init_Rcpp(info);
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
