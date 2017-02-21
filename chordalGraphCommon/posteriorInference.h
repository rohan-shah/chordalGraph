#ifndef POSTERIOR_INFERENCE_HEADER_GUARD
#define POSTERIOR_INFERENCE_HEADER_GUARD
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include "cliqueTree.h"
#include "numericType.h"
#include <boost/random/mersenne_twister.hpp>
namespace chordalGraph
{
	struct posteriorInferenceArgs
	{
		typedef moveable_adjacency_matrix<> graphType;
		posteriorInferenceArgs()
			: delta(0), dataPoints(0), sampleSize(0), burnIn(0)
		{}
		std::size_t delta, dataPoints, dimension;
		boost::numeric::ublas::matrix<double> sampleCovariance;
		boost::numeric::ublas::matrix<double> psi;
		std::vector<std::pair<graphType, double> > results;
		std::vector<mpfr_class> exactValues;
		std::size_t sampleSize, burnIn;
		boost::mt19937 randomSource;
	};
	void posteriorInference(posteriorInferenceArgs& args);
}
#endif
