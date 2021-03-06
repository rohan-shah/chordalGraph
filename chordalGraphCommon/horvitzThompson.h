#ifndef CHORDAL_GRAPH_HORVITZ_THOMPSON_HEADER_GUARD
#define CHORDAL_GRAPH_HORVITZ_THOMPSON_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "numericType.h"
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
		int minimumSizeForExact;
	};
	namespace horvitzThompsonPrivate
	{
		template<typename cliqueTree> struct weightedCliqueTree
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
			void addVertex()
			{
				tree.addVertex();
			}
			cliqueTree tree;
			numericType weight;
		};
	}

	template<typename cliqueTree> void horvitzThompsonReduceChains(horvitzThompsonArgs& args);
	template<typename cliqueTree> void horvitzThompson(horvitzThompsonArgs& args);
}
#endif
