#ifndef CHORDAL_GRAPH_STOCHASTIC_ENUMERATION_NAUTY_HEADER_GUARD
#define CHORDAL_GRAPH_STOCHASTIC_ENUMERATION_NAUTY_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include "numericType.h"
#include "cliqueTree.h"
namespace chordalGraph
{
	namespace stochasticEnumeration
	{
		struct weightedCliqueTree
		{
		public:
			weightedCliqueTree(weightedCliqueTree&& other)
				: tree(std::move(other.tree)), weight(other.weight)
			{} 
			weightedCliqueTree(const weightedCliqueTree& other)
				: tree(other.tree), weight(other.weight)
			{}
			weightedCliqueTree(int nVertices)
				: tree(nVertices), weight(1)
			{}
			cliqueTree tree;
			int weight;
		};
	}
	struct stochasticEnumerationNautyArgs
	{
	public:
		stochasticEnumerationNautyArgs(boost::mt19937& randomSource)
			:randomSource(randomSource)
		{}
		int budget;
		int nVertices;
		int nEdges;
		mpfr_class estimate;
		boost::mt19937& randomSource;
		bool exact;
		int minimumSizeForExact;
	};
	void stochasticEnumerationNauty(stochasticEnumerationNautyArgs& args);
	void stochasticEnumerationNautyReduceChains(stochasticEnumerationNautyArgs& args);
}
#endif
