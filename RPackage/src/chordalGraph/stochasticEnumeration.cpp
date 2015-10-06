#include "stochasticEnumeration.h"
#include <Rcpp.h>
void convert(Rcpp::List& samples, const std::vector<chordalGraph::stochasticEnumerationArgs::matrixType>& boostMatrices)
{
	samples = Rcpp::List(boostMatrices.size());
	for (int i = 0; i < (int)boostMatrices.size(); i++)
	{
		const chordalGraph::stochasticEnumerationArgs::matrixType& currentMatrix = boostMatrices[i];
		int size1 = currentMatrix.size1();
		Rcpp::LogicalMatrix rMatrix(size1);
		for (int j = 0; j < size1; j++)
		{
			for (int k = 0; k < size1; k++)
			{
				rMatrix(j, k) = currentMatrix(j, k);
			}
		}
		samples[i] = rMatrix;
	}
}
SEXP stochasticEnumeration(SEXP nVertices_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp)
{
BEGIN_RCPP
	int seed = Rcpp::as<int>(seed_sexp);
	int nVertices = Rcpp::as<int>(nVertices_sexp);
	int budget = Rcpp::as<int>(budget_sexp);
	Rcpp::List options = options_sexp;
	bool reduceChains = false;
	bool outputSamples = false;
	if (options.containsElementNamed("reduceChains"))
	{
		reduceChains = Rcpp::as<bool>(options("reduceChains"));
	}
	if (options.containsElementNamed("outputSamples"))
	{
		outputSamples = Rcpp::as<bool>(options("outputSamples"));
	}
	

	boost::mt19937 randomSource;
	randomSource.seed(seed);
	chordalGraph::stochasticEnumerationArgs args(randomSource);
	args.nVertices = nVertices;
	args.budget = budget;
	args.outputSamples = outputSamples;

	int maxEdges = ((nVertices - 1)*nVertices / 2);

	//Set up text progress bar
	Rcpp::Function txtProgressBar("txtProgressBar");
	Rcpp::Function close("close");
	Rcpp::Function setTxtProgressBar("setTxtProgressBar");
	Rcpp::RObject barHandle = txtProgressBar(Rcpp::Named("style") = 3, Rcpp::Named("min") = 0, Rcpp::Named("max") = maxEdges, Rcpp::Named("initial") = 0);

	//List for samples
	Rcpp::List samples(maxEdges+1);
	std::vector<bool> exactVector;
	std::vector<std::string> estimatesAsStrings;
	for (int nEdges = 0; nEdges < maxEdges+1; nEdges++)
	{
		setTxtProgressBar(barHandle, nEdges);
		args.nEdges = nEdges;
		args.samples.clear();
		if(reduceChains)
		{
			chordalGraph::stochasticEnumeration2(args);
		}
		else
		{
			chordalGraph::stochasticEnumeration(args);
		}
		estimatesAsStrings.push_back(args.estimate.str());
		exactVector.push_back(args.exact);
		if (args.outputSamples)
		{
			Rcpp::List currentEdgesSamples;
			convert(currentEdgesSamples, args.samples);
			samples[nEdges] = currentEdgesSamples;
		}
	}
	close(barHandle);
	SEXP estimatesAsStrings_sexp = Rcpp::wrap(estimatesAsStrings);

	Rcpp::Function mpfrFunction("mpfr");
	Rcpp::List listResult = Rcpp::List::create(Rcpp::Named("data") = mpfrFunction(estimatesAsStrings_sexp, Rcpp::Named("prec", 50)), Rcpp::Named("exact") = Rcpp::wrap(exactVector));
	if (args.outputSamples)
	{
		listResult("samples") = samples;
	}
	return listResult;
END_RCPP
}
SEXP stochasticEnumerationSpecificEdges(SEXP nVertices_sexp, SEXP nEdges_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp)
{
BEGIN_RCPP
	int seed = Rcpp::as<int>(seed_sexp);
	int nVertices = Rcpp::as<int>(nVertices_sexp);
	int nEdges = Rcpp::as<int>(nEdges_sexp);
	int budget = Rcpp::as<int>(budget_sexp);
	Rcpp::List options = options_sexp;
	bool reduceChains = false;
	bool outputSamples = false;
	if (options.containsElementNamed("reduceChains"))
	{
		reduceChains = Rcpp::as<bool>(options("reduceChains"));
	}
	if (options.containsElementNamed("outputSamples"))
	{
		outputSamples = Rcpp::as<bool>(options("outputSamples"));
	}

	boost::mt19937 randomSource;
	randomSource.seed(seed);
	chordalGraph::stochasticEnumerationArgs args(randomSource);
	args.nEdges = nEdges;
	args.nVertices = nVertices;
	args.budget = budget;
	args.outputSamples = outputSamples;

	if(reduceChains)
	{
		chordalGraph::stochasticEnumeration(args);
	}
	else
	{
		chordalGraph::stochasticEnumeration2(args);
	}

	std::string estimateAsString = args.estimate.str();
	SEXP estimateAsString_sexp = Rcpp::wrap(estimateAsString);

	Rcpp::Function mpfrFunction("mpfr");
	Rcpp::List listResult = Rcpp::List::create(Rcpp::Named("data") = mpfrFunction(estimateAsString_sexp, Rcpp::Named("prec", 50)), Rcpp::Named("exact") = Rcpp::wrap(args.exact));

	if (args.outputSamples)
	{
		Rcpp::List samples;
		convert(samples, args.samples);
		listResult("samples") = samples;
	}
	return listResult;
END_RCPP
}
