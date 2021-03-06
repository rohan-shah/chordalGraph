#include "stochasticEnumerationNauty.h"
#include <Rcpp.h>
#include "graphRepresentation.h"
#include "cliqueTree.h"
#include "cliqueTreeAdjacencyMatrix.h"
SEXP stochasticEnumerationNauty(SEXP nVertices_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp)
{
BEGIN_RCPP
	int seed = Rcpp::as<int>(seed_sexp);
	int nVertices = Rcpp::as<int>(nVertices_sexp);
	int budget = Rcpp::as<int>(budget_sexp);
	Rcpp::List options = options_sexp;
	if(!options.containsElementNamed("reduceChains")) throw std::runtime_error("Unable to find option named reduceChains");
	if(!options.containsElementNamed("graphRepresentation")) throw std::runtime_error("Unable to find option named graphRepresentation");
	bool reduceChains = Rcpp::as<bool>(options("reduceChains"));
	std::string graphRepresentationString = Rcpp::as<std::string>(options("graphRepresentation"));
	graphRepresentation representation = toRepresentation(graphRepresentationString);

	boost::mt19937 randomSource;
	randomSource.seed(seed);
	chordalGraph::stochasticEnumerationNautyArgs args(randomSource);
	args.nVertices = nVertices;
	args.budget = budget;

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
		if(reduceChains)
		{
			if(representation == listRepresentation)
			{
				chordalGraph::stochasticEnumerationNautyReduceChains<chordalGraph::cliqueTree>(args);
			}
			else
			{
				chordalGraph::stochasticEnumerationNautyReduceChains<chordalGraph::cliqueTreeAdjacencyMatrix>(args);
			}
		}
		else 
		{
			if(representation == listRepresentation)
			{
				chordalGraph::stochasticEnumerationNauty<chordalGraph::cliqueTree>(args);
			}
			else
			{
				chordalGraph::stochasticEnumerationNauty<chordalGraph::cliqueTreeAdjacencyMatrix>(args);
			}
		}
		estimatesAsStrings.push_back(args.estimate.str());
		exactVector.push_back(args.exact);
		if(args.minimumSizeForExact == -1) minimumSizeForExactVector.push_back(NA_INTEGER);
		else minimumSizeForExactVector.push_back(args.minimumSizeForExact);
	}
	close(barHandle);
	SEXP estimatesAsStrings_sexp = Rcpp::wrap(estimatesAsStrings);
	SEXP minimumSizeForExactVector_sexp = Rcpp::wrap(minimumSizeForExactVector);

	Rcpp::Function mpfrFunction("mpfr");
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("data") = mpfrFunction(estimatesAsStrings_sexp, Rcpp::Named("prec", 50)), Rcpp::Named("exact") = Rcpp::wrap(exactVector), Rcpp::Named("minimumSizeForExact") = minimumSizeForExactVector_sexp);
	return result;
END_RCPP
}
SEXP stochasticEnumerationNautySpecificEdges(SEXP nVertices_sexp, SEXP nEdges_sexp, SEXP budget_sexp, SEXP seed_sexp, SEXP options_sexp)
{
BEGIN_RCPP
	int seed = Rcpp::as<int>(seed_sexp);
	int nVertices = Rcpp::as<int>(nVertices_sexp);
	int nEdges = Rcpp::as<int>(nEdges_sexp);
	int budget = Rcpp::as<int>(budget_sexp);
	Rcpp::List options = options_sexp;
	if(!options.containsElementNamed("reduceChains")) throw std::runtime_error("Unable to find option named reduceChains");
	if(!options.containsElementNamed("graphRepresentation")) throw std::runtime_error("Unable to find option named graphRepresentation");
	bool reduceChains = Rcpp::as<bool>(options("reduceChains"));
	std::string graphRepresentationString = Rcpp::as<std::string>(options("graphRepresentation"));
	graphRepresentation representation = toRepresentation(graphRepresentationString);

	boost::mt19937 randomSource;
	randomSource.seed(seed);
	chordalGraph::stochasticEnumerationNautyArgs args(randomSource);
	args.nEdges = nEdges;
	args.nVertices = nVertices;
	args.budget = budget;

	if(reduceChains)
	{
		if(representation == listRepresentation)
		{
			chordalGraph::stochasticEnumerationNautyReduceChains<chordalGraph::cliqueTree>(args);
		}
		else
		{
			chordalGraph::stochasticEnumerationNautyReduceChains<chordalGraph::cliqueTreeAdjacencyMatrix>(args);
		}
	}
	else
	{
		if(representation == listRepresentation)
		{
			chordalGraph::stochasticEnumerationNauty<chordalGraph::cliqueTree>(args);
		}
		else
		{
			chordalGraph::stochasticEnumerationNauty<chordalGraph::cliqueTreeAdjacencyMatrix>(args);
		}
	}

	std::string estimateAsString = args.estimate.str();
	SEXP estimateAsString_sexp = Rcpp::wrap(estimateAsString);
	SEXP minimumSizeForExact_sexp;
	if(args.minimumSizeForExact == -1) minimumSizeForExact_sexp = Rcpp::wrap(NA_INTEGER);
	else minimumSizeForExact_sexp = Rcpp::wrap(args.minimumSizeForExact);

	Rcpp::Function mpfrFunction("mpfr");
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("data") = mpfrFunction(estimateAsString_sexp, Rcpp::Named("prec", 50)), Rcpp::Named("exact") = args.exact, Rcpp::Named("minimumSizeForExact") = minimumSizeForExact_sexp);

	return result;
END_RCPP
}
