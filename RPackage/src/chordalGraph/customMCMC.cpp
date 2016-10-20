#include <Rcpp.h>
#include "customMCMC.h"
SEXP customMCMC(SEXP nVertices_sexp, SEXP approximateCounts_sexp, SEXP seed_sexp, SEXP burnIn_sexp, SEXP runSize_sexp)
{
BEGIN_RCPP
	int nVertices;
	try
	{
		nVertices = Rcpp::as<int>(nVertices_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input nVertices must be an integer");
	}

	std::vector<std::string> approximateCounts_string;
	try
	{
		approximateCounts_string = Rcpp::as<std::vector<std::string> >(approximateCounts_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input approximateCounts must be a character vector");
	}
	std::vector<mpfr_class> approximateCounts;
	std::transform(approximateCounts_string.begin(), approximateCounts_string.end(), std::back_inserter(approximateCounts), [](std::string x){return mpfr_class(x);});

	int seed;
	try
	{
		seed = Rcpp::as<int>(seed_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input seed must be an integer");
	}

	int burnIn;
	try
	{
		burnIn = Rcpp::as<int>(burnIn_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input burnIn must be an integer");
	}

	int runSize;
	try
	{
		runSize = Rcpp::as<int>(runSize_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input runSize must be an integer");
	}
	boost::mt19937 randomSource;
	randomSource.seed(seed);

	chordalGraph::mcmcArgs args(randomSource);
	args.nVertices = nVertices;
	args.approximateCounts.swap(approximateCounts);
	args.burnIn = burnIn;
	args.runSize = runSize;
	chordalGraph::customMCMC(args);

	std::vector<std::string> returnValues;
	std::transform(args.estimates.begin(), args.estimates.end(), std::back_inserter(returnValues), [](mpfr_class x){return x.str(10, std::ios_base::dec); });

	return Rcpp::wrap(returnValues);
END_RCPP
}
