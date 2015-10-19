#include "horvitzThompson.h"
#include <Rcpp.h>
SEXP horvitzThompson(SEXP nVertices_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP sampling_sexp)
{
BEGIN_RCPP
	int seed = Rcpp::as<int>(seed_sexp);
	int nVertices = Rcpp::as<int>(nVertices_sexp);
	int budget = Rcpp::as<int>(budget_sexp);
	std::string samplingString = Rcpp::as<std::string>(sampling_sexp);
	chordalGraph::samplingType sampling = chordalGraph::toSamplingType(samplingString);

	boost::mt19937 randomSource;
	randomSource.seed(seed);
	chordalGraph::horvitzThompsonArgs args(randomSource);
	args.nVertices = nVertices;
	args.budget = budget;
	args.sampling = sampling;

	int maxEdges = ((nVertices - 1)*nVertices / 2);

	//Set up text progress bar
	Rcpp::Function txtProgressBar("txtProgressBar");
	Rcpp::Function close("close");
	Rcpp::Function setTxtProgressBar("setTxtProgressBar");
	Rcpp::RObject barHandle = txtProgressBar(Rcpp::Named("style") = 3, Rcpp::Named("min") = 0, Rcpp::Named("max") = maxEdges, Rcpp::Named("initial") = 0);

	//List for samples
	Rcpp::List samples(maxEdges+1);
	std::vector<std::string> estimatesAsStrings;
	std::vector<bool> exactVector;
	for (int nEdges = 0; nEdges < maxEdges+1; nEdges++)
	{
		setTxtProgressBar(barHandle, nEdges);
		args.nEdges = nEdges;
		chordalGraph::horvitzThompson(args);
		estimatesAsStrings.push_back(args.estimate.str());
		exactVector.push_back(args.exact);
	}
	close(barHandle);
	SEXP estimatesAsStrings_sexp = Rcpp::wrap(estimatesAsStrings);

	Rcpp::Function mpfrFunction("mpfr");
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("data") = mpfrFunction(estimatesAsStrings_sexp, Rcpp::Named("prec", 50)), Rcpp::Named("exact") = Rcpp::wrap(exactVector));
	return result;
END_RCPP
}
SEXP horvitzThompsonSpecificEdges(SEXP nVertices_sexp, SEXP nEdges_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP sampling_sexp)
{
BEGIN_RCPP
	int seed = Rcpp::as<int>(seed_sexp);
	int nVertices = Rcpp::as<int>(nVertices_sexp);
	int nEdges = Rcpp::as<int>(nEdges_sexp);
	int budget = Rcpp::as<int>(budget_sexp);
	std::string samplingString = Rcpp::as<std::string>(sampling_sexp);
	chordalGraph::samplingType sampling = chordalGraph::toSamplingType(samplingString);

	boost::mt19937 randomSource;
	randomSource.seed(seed);
	chordalGraph::horvitzThompsonArgs args(randomSource);
	args.nEdges = nEdges;
	args.nVertices = nVertices;
	args.budget = budget;
	args.sampling = sampling;

	chordalGraph::horvitzThompson(args);

	std::string estimateAsString = args.estimate.str();
	SEXP estimateAsString_sexp = Rcpp::wrap(estimateAsString);

	Rcpp::Function mpfrFunction("mpfr");
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("data") = mpfrFunction(estimateAsString_sexp, Rcpp::Named("prec", 50)), Rcpp::Named("exact") = args.exact);

	return result;
END_RCPP
}
