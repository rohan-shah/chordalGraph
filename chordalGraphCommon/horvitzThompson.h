#ifndef CHORDAL_GRAPH_HORVITZ_THOMPSON_HEADER_GUARD
#define CHORDAL_GRAPH_HORVITZ_THOMPSON_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "numericType.h"
namespace chordalGraph
{
	enum samplingType
	{
		sampfordSamplingMultinomial, conditionalPoissonSampling, paretoSampling, semiDeterministicSampling, sampfordSamplingConditionalPoisson, sampfordSamplingFromParetoNaive
	};
	struct horvitzThompsonArgs
	{
	public:
		horvitzThompsonArgs(boost::mt19937& randomSource)
			:randomSource(randomSource)
		{}
		int budget;
		int nVertices;
		int nEdges;
		mpfr_class estimate;
		boost::mt19937& randomSource;
		bool exact;
		int minimumSizeForExact;
		samplingType sampling;
	};
	samplingType toSamplingType(std::string samplingString);
	void horvitzThompson(horvitzThompsonArgs& args);
}
#endif
