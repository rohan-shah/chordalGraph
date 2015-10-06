#include <Rcpp.h>
#include "isChordalIterative.h"
#ifdef _MSC_VER
	#undef RcppExport
	#define RcppExport extern "C" __declspec(dllexport)
#endif
extern "C" const char* package_name = "chordalGraph";
SEXP stochasticEnumeration(SEXP nVertices_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp);
SEXP stochasticEnumerationSpecificEdges(SEXP nVertices_sexp, SEXP nEdges_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp);
SEXP stochasticEnumerationNauty(SEXP nVertices_sexp, SEXP budget_sexp, SEXP seed_sexp);
SEXP stochasticEnumerationNautySpecificEdges(SEXP nVertices_sexp, SEXP nEdges_sexp, SEXP budget_sexp, SEXP seed_sexp);
R_CallMethodDef callMethods[] = 
{
	{"isChordalIterative_igraph", (DL_FUNC)&isChordalIterative_igraph, 1},
	{"isChordalIterative_graphNEL", (DL_FUNC)&isChordalIterative_graphNEL, 1},
	{"isChordalIterative_graphAM", (DL_FUNC)&isChordalIterative_graphAM, 1},
	{"stochasticEnumeration", (DL_FUNC)&stochasticEnumeration, 4},
	{"stochasticEnumerationSpecificEdges", (DL_FUNC)&stochasticEnumerationSpecificEdges, 5},
	{"stochasticEnumerationNauty", (DL_FUNC)&stochasticEnumerationNauty, 3},
	{"stochasticEnumerationNautySpecificEdges", (DL_FUNC)&stochasticEnumerationNautySpecificEdges, 4},
	{NULL, NULL, 0}
};
extern "C" void R_init_Rcpp(DllInfo* info);
RcppExport void R_init_chordalGraph(DllInfo *info)
{
	R_init_Rcpp(info);
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
