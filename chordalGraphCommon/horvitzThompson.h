#ifndef CHORDAL_GRAPH_HORVITZ_THOMPSON_HEADER_GUARD
#define CHORDAL_GRAPH_HORVITZ_THOMPSON_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "numericType.h"
#include "cliqueTree.h"
#include "performSampling.h"
namespace chordalGraph
{
	enum weightType
	{
		weightsMultiplicity, weightsAutomorphismGroup
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
	namespace horvitzThompsonPrivate
	{
		struct weightedCliqueTree
		{
		public:
			weightedCliqueTree(weightedCliqueTree&& other)
				: tree(std::move(other.tree)), weight(other.weight), automorphismGroupSize(other.automorphismGroupSize)
			{} 
			weightedCliqueTree(const weightedCliqueTree& other)
				: tree(other.tree), weight(other.weight), automorphismGroupSize(other.automorphismGroupSize)
			{}
			weightedCliqueTree(int nVertices)
				: tree(nVertices), weight(1), automorphismGroupSize(1)
			{}
			cliqueTree tree;
			numericType weight;
			mpz_class automorphismGroupSize;
		};
	}

	samplingType toSamplingType(std::string samplingString);
	weightType toWeightType(std::string weightString);
	void horvitzThompson(horvitzThompsonArgs& args);
}
#endif
