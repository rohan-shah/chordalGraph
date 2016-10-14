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
	namespace horvitzThompson2Private
	{
		template<typename cliqueTree> struct weightedCliqueTree
		{
		public:
			weightedCliqueTree(weightedCliqueTree&& other)
				: tree(std::move(other.tree)), weight(std::move(other.weight)), auxWeight(std::move(other.auxWeight))
			{} 
			weightedCliqueTree(const weightedCliqueTree& other)
				: tree(other.tree), weight(other.weight), auxWeight(other.auxWeight)
			{}
			weightedCliqueTree(int nVertices)
				: tree(nVertices), weight(1), auxWeight(1)
			{}
			void addVertex()
			{
				tree.addVertex();
			}
			cliqueTree tree;
			numericType weight;
			numericType auxWeight;
		};
		template <typename weightType> struct childNode
		{
		public:
			childNode(int parentIndex, bool includesEdge, weightType weight, weightType auxWeight)
				:weight(weight), auxWeight(auxWeight), value(2*parentIndex + includesEdge)
			{}
			int getParentIndex() const
			{
				return value / 2;
			}
			bool includesEdge() const 
			{
				return (value % 2) == 1;
			}
			weightType weight, auxWeight;
			bool operator<(const childNode& other) const
			{
				return value < other.value;
			}
		private:
			int value;
		};
	}

	template<typename cliqueTree> void horvitzThompson2ReduceChains(horvitzThompsonArgs& args);
	template<typename cliqueTree> void horvitzThompson2(horvitzThompsonArgs& args);
}
#endif
