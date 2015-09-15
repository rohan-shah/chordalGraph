#include <Rcpp.h>
#ifdef _MSC_VER
	#undef RcppExport
	#define RcppExport extern "C" __declspec(dllexport)
#endif
extern "C" const char* package_name = "chordalSubgraph";
R_CallMethodDef callMethods[] = 
{
	//{"createHexagonalLattice", (DL_FUNC)&createHexagonalLattice, 2},
	{NULL, NULL, 0}
};
extern "C" void R_init_Rcpp(DllInfo* info);
RcppExport void R_init_chordalSubgraph(DllInfo *info)
{
	R_init_Rcpp(info);
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
