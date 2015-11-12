#ifndef SEMIDETERMINISTIC_SAMPLING_HEADER_GUARD
#define SEMIDETERMINISTIC_SAMPLING_HEADER_GUARD
#include <vector>
#include "numericType.h"
#include <boost/random/mersenne_twister.hpp>
namespace chordalGraph
{
	struct semiDeterministicSamplingArgs
	{
	public:
		semiDeterministicSamplingArgs(boost::mt19937& randomSource)
			: randomSource(randomSource)
		{}
		boost::mt19937& randomSource;
		std::vector<numericType> weights;
		std::vector<int> indices;
		std::size_t n;
		std::vector<numericType> inclusionProbabilities;
		std::vector<bool> deterministicInclusion;
		std::vector<int> deterministicIndices;
	};
	void semiDeterministic(semiDeterministicSamplingArgs& args);
}
#endif
