#ifndef SAMPFORD_SAMPLING_HEADER_GUARD
#define SAMPFORD_SAMPLING_HEADER_GUARD
#include <vector>
#include <boost/random/mersenne_twister.hpp>
namespace chordalGraph
{
	struct sampfordBruteForceArgs
	{
	public:
		sampfordBruteForceArgs(boost::mt19937& randomSource)
			: randomSource(randomSource)
		{}
		boost::mt19937& randomSource;
		std::vector<double> weights;
		std::vector<int> indices;
		std::size_t n;
		std::vector<double> accumulated1;
		std::vector<double> accumulated2;
		std::vector<double> inclusionProbabilities;
		std::vector<bool> deterministicInclusion;
		std::vector<int> deterministicIndices;
	};
	void sampfordBruteForce(sampfordBruteForceArgs& args);
}
#endif
