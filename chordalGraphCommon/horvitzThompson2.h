#ifndef CHORDAL_GRAPH_HORVITZ_THOMPSON_2_HEADER_GUARD
#define CHORDAL_GRAPH_HORVITZ_THOMPSON_2_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "horvitzThompson.h"
#include "numericType.h"
namespace chordalGraph
{
	const double auxWeightPower = 1.0/2.0;
	template<typename cliqueTree> void horvitzThompson2ReduceChains(horvitzThompsonArgs& args);
	template<typename cliqueTree> void horvitzThompson2(horvitzThompsonArgs& args);
}
#endif
