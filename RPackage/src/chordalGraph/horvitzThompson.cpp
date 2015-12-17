#include "horvitzThompson.h"
#include <Rcpp.h>
SEXP horvitzThompson(SEXP nVertices_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP sampling_sexp, SEXP weight_sexp)
{
BEGIN_RCPP
	int seed = Rcpp::as<int>(seed_sexp);
	int nVertices = Rcpp::as<int>(nVertices_sexp);
	int budget = Rcpp::as<int>(budget_sexp);
	std::string samplingString = Rcpp::as<std::string>(sampling_sexp);
	std::string weightString = Rcpp::as<std::string>(weight_sexp);
	chordalGraph::samplingType sampling = chordalGraph::toSamplingType(samplingString);
	chordalGraph::weightType weight = chordalGraph::toWeightType(weightString);

	boost::mt19937 randomSource;
	randomSource.seed(seed);
	chordalGraph::horvitzThompsonArgs args(randomSource);
	args.nVertices = nVertices;
	args.budget = budget;
	args.sampling = sampling;
	args.weights = weight;

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
	std::vector<int> minimumSizeForExactVector;
	for (int nEdges = 0; nEdges < maxEdges+1; nEdges++)
	{
		setTxtProgressBar(barHandle, nEdges);
		args.nEdges = nEdges;
		chordalGraph::horvitzThompson(args);
		estimatesAsStrings.push_back(args.estimate.str());
		exactVector.push_back(args.exact);
		if(args.minimumSizeForExact != -1)
		{
			minimumSizeForExactVector.push_back(args.minimumSizeForExact);
		}
		else minimumSizeForExactVector.push_back(NA_INTEGER);
	}
	close(barHandle);
	SEXP estimatesAsStrings_sexp = Rcpp::wrap(estimatesAsStrings);

	Rcpp::Function mpfrFunction("mpfr");
	SEXP estimatesAsMPFR = mpfrFunction(estimatesAsStrings_sexp, Rcpp::Named("prec", 50));
	SEXP exactVector_sexp = Rcpp::wrap(exactVector);
	SEXP minimumSizeForExactVector_sexp = Rcpp::wrap(minimumSizeForExactVector);
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("data") = estimatesAsMPFR, Rcpp::Named("exact") = exactVector_sexp, Rcpp::Named("minimumSizeForExact") = minimumSizeForExactVector_sexp);
	return result;
END_RCPP
}
SEXP horvitzThompsonSpecificEdges(SEXP nVertices_sexp, SEXP nEdges_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP sampling_sexp, SEXP weight_sexp)
{
BEGIN_RCPP
	int seed = Rcpp::as<int>(seed_sexp);
	int nVertices = Rcpp::as<int>(nVertices_sexp);
	int nEdges = Rcpp::as<int>(nEdges_sexp);
	int budget = Rcpp::as<int>(budget_sexp);
	std::string samplingString = Rcpp::as<std::string>(sampling_sexp);
	std::string weightString = Rcpp::as<std::string>(weight_sexp);
	chordalGraph::samplingType sampling = chordalGraph::toSamplingType(samplingString);
	chordalGraph::weightType weight = chordalGraph::toWeightType(weightString);

	boost::mt19937 randomSource;
	randomSource.seed(seed);
	chordalGraph::horvitzThompsonArgs args(randomSource);
	args.nEdges = nEdges;
	args.nVertices = nVertices;
	args.budget = budget;
	args.sampling = sampling;
	args.weights = weight;

	chordalGraph::horvitzThompson(args);

	std::string estimateAsString = args.estimate.str();
	SEXP estimateAsString_sexp = Rcpp::wrap(estimateAsString);

	Rcpp::Function mpfrFunction("mpfr");
	SEXP estimatesAsMPFR = mpfrFunction(estimateAsString_sexp, Rcpp::Named("prec", 50));
	SEXP exact_sexp = Rcpp::wrap(args.exact);
	SEXP minimumSizeForExact_sexp;
	if(args.minimumSizeForExact != -1) minimumSizeForExact_sexp = Rcpp::wrap(args.minimumSizeForExact);
	else minimumSizeForExact_sexp = Rcpp::wrap(NA_INTEGER);
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("data") = estimatesAsMPFR, Rcpp::Named("exact") = exact_sexp, Rcpp::Named("minimumSizeForExact") = minimumSizeForExact_sexp);
	return result;
END_RCPP
}
