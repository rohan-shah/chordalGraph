#include "stochasticEnumeration.h"
#include <Rcpp.h>
SEXP stochasticEnumeration(SEXP nVertices_sexp, SEXP budget_sexp, SEXP seed_sexp)
{
BEGIN_RCPP
	int seed = Rcpp::as<int>(seed_sexp);
	int nVertices = Rcpp::as<int>(nVertices_sexp);
	int budget = Rcpp::as<int>(budget_sexp);

	boost::mt19937 randomSource;
	randomSource.seed(seed);
	chordalGraph::stochasticEnumerationArgs args(randomSource);
	args.nVertices = nVertices;
	args.budget = budget;

	std::vector<std::string> estimatesAsStrings;
	for (int nEdges = 0; nEdges < ((nVertices-1)*nVertices/2)+1; nEdges++)
	{
		args.nEdges = nEdges;
		chordalGraph::stochasticEnumeration(args);
		estimatesAsStrings.push_back(args.estimate.str());
	}
	SEXP estimatesAsStrings_sexp = Rcpp::wrap(estimatesAsStrings);

	Rcpp::Function mpfrFunction("mpfr");
	Rcpp::RObject result = mpfrFunction(estimatesAsStrings_sexp, Rcpp::Named("prec", 50));
	return result;
END_RCPP
}
SEXP stochasticEnumerationSpecificEdges(SEXP nVertices_sexp, SEXP nEdges_sexp, SEXP budget_sexp, SEXP seed_sexp)
{
BEGIN_RCPP
	int seed = Rcpp::as<int>(seed_sexp);
	int nVertices = Rcpp::as<int>(nVertices_sexp);
	int nEdges = Rcpp::as<int>(nEdges_sexp);
	int budget = Rcpp::as<int>(budget_sexp);

	boost::mt19937 randomSource;
	randomSource.seed(seed);
	chordalGraph::stochasticEnumerationArgs args(randomSource);
	args.nEdges = nEdges;
	args.nVertices = nVertices;
	args.budget = budget;

	chordalGraph::stochasticEnumeration(args);

	std::string estimateAsString = args.estimate.str();
	SEXP estimateAsString_sexp = Rcpp::wrap(estimateAsString);

	Rcpp::Function mpfrFunction("mpfr");
	Rcpp::RObject result = mpfrFunction(estimateAsString_sexp, Rcpp::Named("prec", 50));
	return result;
END_RCPP
}