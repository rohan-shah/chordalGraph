#ifndef CONDITIONAL_POISSON_HEADER_GUARD
#define CONDITIONAL_POISSON_HEADER_GUARD
#include <vector>
#include "numericType.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/numeric/ublas/matrix.hpp>
namespace chordalGraph
{
	struct conditionalPoissonArgs
	{
	public:
		conditionalPoissonArgs()
		{}
		std::size_t n;
		std::vector<bool> deterministicInclusion;
		std::vector<numericType> exponentialParameters;
		std::vector<numericType> expExponentialParameters;
		boost::numeric::ublas::matrix<numericType> expNormalisingConstant;
		bool calculateInclusionProbabilities;
	};
	void conditionalPoisson(conditionalPoissonArgs& args, std::vector<int>& indices, std::vector<numericType>& inclusionProbabilities, std::vector<numericType>& weights, boost::mt19937& randomSource);
	void conditionalPoissonInclusionProbabilities(conditionalPoissonArgs& args, std::vector<numericType>& inclusionProbabilities, std::vector<numericType>& weights);
	void calculateExpNormalisingConstants(std::vector<numericType>& expExponentialParameters, std::vector<numericType>& exponentialParameters, boost::numeric::ublas::matrix<numericType>& expNormalisingConstant, int n, int nUnits, std::vector<bool>& ignore);
}
#endif
