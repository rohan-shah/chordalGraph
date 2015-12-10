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
		semiDeterministicSamplingArgs()
		{}
		std::size_t n;
		std::vector<bool> deterministicInclusion;
		std::vector<int> deterministicIndices;
	};
	void semiDeterministic(semiDeterministicSamplingArgs& args, std::vector<int>& indices, std::vector<numericType>& inclusionProbabilities, std::vector<numericType>& weights, boost::mt19937& randomSource);
}
#endif
