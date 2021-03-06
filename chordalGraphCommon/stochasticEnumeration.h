#ifndef CHORDAL_GRAPH_STOCHASTIC_ENUMERATION_HEADER_GUARD
#define CHORDAL_GRAPH_STOCHASTIC_ENUMERATION_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include "numericType.h"
typedef boost::multiprecision::mpz_int mpz_class;
namespace chordalGraph
{
	struct stochasticEnumerationArgs
	{
	public:
		stochasticEnumerationArgs(boost::mt19937& randomSource)
			:randomSource(randomSource)
		{}
		int budget;
		int nVertices;
		int nEdges;
		bool outputSamples;
		typedef boost::numeric::ublas::symmetric_matrix<bool> matrixType;
		std::vector<matrixType> samples;
		mpfr_class estimate;
		boost::mt19937& randomSource;
		bool exact;
		int minimumSizeForExact;
	};
	template<typename cliqueTree> void stochasticEnumeration(stochasticEnumerationArgs& args);
	template<typename cliqueTree> void stochasticEnumerationReduceChains(stochasticEnumerationArgs& args);
}
#endif
