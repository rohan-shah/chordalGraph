#include <Rcpp.h>
#include <internal.h>
#include "isChordalIterative.h"
#ifdef _MSC_VER
	#undef RcppExport
	#define RcppExport extern "C" __declspec(dllexport)
#endif
extern "C" const char* package_name = "chordalGraph";
SEXP stochasticEnumeration(SEXP nVertices_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp);
SEXP stochasticEnumerationSpecificEdges(SEXP nVertices_sexp, SEXP nEdges_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp);
SEXP stochasticEnumerationNauty(SEXP nVertices_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp);
SEXP stochasticEnumerationNautySpecificEdges(SEXP nVertices_sexp, SEXP nEdges_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp);
SEXP horvitzThompson(SEXP nVertices_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp);
SEXP horvitzThompson2(SEXP nVertices_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp);
SEXP horvitzThompsonSpecificEdges(SEXP nVertices_sexp, SEXP nEdges_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp);
R_CallMethodDef callMethods[] = 
{
	{"isChordalIterative_igraph", (DL_FUNC)&isChordalIterative_igraph, 1},
	{"isChordalIterative_graphNEL", (DL_FUNC)&isChordalIterative_graphNEL, 1},
	{"isChordalIterative_graphAM", (DL_FUNC)&isChordalIterative_graphAM, 1},
	{"stochasticEnumeration", (DL_FUNC)&stochasticEnumeration, 4},
	{"stochasticEnumerationSpecificEdges", (DL_FUNC)&stochasticEnumerationSpecificEdges, 5},
	{"stochasticEnumerationNauty", (DL_FUNC)&stochasticEnumerationNauty, 4},
	{"stochasticEnumerationNautySpecificEdges", (DL_FUNC)&stochasticEnumerationNautySpecificEdges, 5},
	{"horvitzThompson", (DL_FUNC)&horvitzThompson, 4},
	{"horvitzThompson2", (DL_FUNC)&horvitzThompson2, 4},
	{"horvitzThompsonSpecificEdges", (DL_FUNC)&horvitzThompsonSpecificEdges, 5},
	{NULL, NULL, 0}
};
RcppExport void R_init_chordalGraph(DllInfo *info)
{
	std::vector<R_CallMethodDef> callMethodsVector;
	R_CallMethodDef* mpMap2CallMethods = callMethods;
	while(mpMap2CallMethods->name != NULL) mpMap2CallMethods++;
	callMethodsVector.insert(callMethodsVector.begin(), callMethods, mpMap2CallMethods);

#ifdef CUSTOM_STATIC_RCPP
	R_CallMethodDef* RcppStartCallMethods = Rcpp_get_call();
	R_CallMethodDef* RcppEndCallMethods = RcppStartCallMethods;
	while(RcppEndCallMethods->name != NULL) RcppEndCallMethods++;
	callMethodsVector.insert(callMethodsVector.end(), RcppStartCallMethods, RcppEndCallMethods);
#endif
	R_CallMethodDef blank = {NULL, NULL, 0};
	callMethodsVector.push_back(blank);

	R_registerRoutines(info, NULL, &(callMethodsVector[0]), NULL, NULL);
#ifdef CUSTOM_STATIC_RCPP
	init_Rcpp_cache();
#endif
}
