#include "customSymmetricPosteriorInferenceRPackage.h"
#include "customSymmetricPosteriorInference.h"
SEXP customSymmetricPosteriorInference(SEXP outerProductsSum_sexp, SEXP delta_sexp, SEXP dimension_sexp, SEXP dataPoints_sexp, SEXP psi_sexp, SEXP exactCounts_sexp, SEXP burnIn_sexp, SEXP runSize_sexp, SEXP seed_sexp)
{
BEGIN_RCPP
	chordalGraph::customSymmetricPosteriorInferenceArgs args;

	args.dimension = Rcpp::as<int>(dimension_sexp);
	args.dataPoints = Rcpp::as<int>(dataPoints_sexp);

	Rcpp::NumericMatrix outerProductsSum = Rcpp::as<Rcpp::NumericMatrix>(outerProductsSum_sexp);
	args.sampleCovariance.resize(args.dimension, args.dimension, false);
	
	Rcpp::NumericMatrix psi = Rcpp::as<Rcpp::NumericMatrix>(psi_sexp);

	int seed = Rcpp::as<int>(seed_sexp);

	if(psi.nrow() != (int)args.dimension || psi.ncol() != (int)args.dimension)
	{
		throw std::runtime_error("Input psi had the wrong dimensions");
	}
	args.psi.resize(args.dimension, args.dimension, false);
	for(int i = 0; i < (int)args.dimension; i++)
	{
		for(int j = 0; j < (int)args.dimension; j++)
		{
			args.sampleCovariance(i, j) = outerProductsSum(i, j);
			args.psi(i, j) = psi(i, j);
		}
	}
	args.delta = Rcpp::as<int>(delta_sexp);

	Rcpp::CharacterVector exactCounts = Rcpp::as<Rcpp::CharacterVector>(exactCounts_sexp);
	args.exactValues.resize(exactCounts.size());
	for(int i = 0; i <= (int)(args.dimension*(args.dimension-1))/2; i++)
	{
		args.exactValues[i] = mpfr_class(Rcpp::as<std::string>(exactCounts[i]));
	}
	args.burnIn = Rcpp::as<int>(burnIn_sexp);
	args.sampleSize = Rcpp::as<int>(runSize_sexp);
	args.randomSource.seed(seed);
	customSymmetricPosteriorInference(args);

	Rcpp::List graphs(args.results.size());
	Rcpp::NumericVector probabilities(args.results.size());
	for(chordalGraph::customSymmetricPosteriorInferenceArgs::resultsType::iterator i = args.results.begin(); i != args.results.end(); i++)
	{
		Rcpp::IntegerMatrix currentGraphR(args.dimension, args.dimension);
		std::fill(currentGraphR.begin(), currentGraphR.end(), 0);
		const chordalGraph::customSymmetricPosteriorInferenceArgs::graphType& currentGraph = i->first;
		chordalGraph::customSymmetricPosteriorInferenceArgs::graphType::edge_iterator current, end;
		boost::tie(current, end) = boost::edges(currentGraph);
		for(; current != end; current++)
		{
			int source = boost::source(*current, currentGraph), target = boost::target(*current, currentGraph);
			currentGraphR(source, target) = currentGraphR(target, source) = 1;
		}
		graphs(std::distance(args.results.begin(), i)) = currentGraphR;
		probabilities(std::distance(args.results.begin(), i)) = (double)i->second / args.sampleSize;
	}
	Rcpp::List results = Rcpp::List::create(Rcpp::Named("graphs") = graphs, Rcpp::Named("probabilities") = probabilities);
	return results;
END_RCPP
}
