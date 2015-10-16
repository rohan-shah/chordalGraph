#ifndef CHORDAL_GRAPH_HORVITZ_THOMPSON_HEADER_GUARD
#define CHORDAL_GRAPH_HORVITZ_THOMPSON_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
typedef boost::multiprecision::static_mpfr_float_50 mpfr_class;
typedef boost::multiprecision::mpz_int mpz_class;
namespace chordalGraph
{
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
	};
	void horvitzThompson(horvitzThompsonArgs& args);
}
#endif
