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
		conditionalPoissonArgs(boost::mt19937& randomSource)
			: randomSource(randomSource)
		{}
		boost::mt19937& randomSource;
		std::vector<numericType> inclusionProbabilities;
		std::vector<int> indices;
		std::size_t n;
		std::vector<numericType> weights;
		std::vector<bool> deterministicInclusion;
		std::vector<numericType> exponentialParameters;
		std::vector<numericType> expExponentialParameters;
		boost::numeric::ublas::matrix<numericType> expNormalisingConstant;
	};
	void conditionalPoisson(conditionalPoissonArgs& args);
	void calculateExpNormalisingConstants(std::vector<numericType>& expExponentialParameters, std::vector<numericType>& exponentialParameters, boost::numeric::ublas::matrix<numericType>& expNormalisingConstant, int n, int nUnits, std::vector<bool>& ignore);
}
#endif
