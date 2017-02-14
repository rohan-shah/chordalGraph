#ifndef CUSTOM_MCMC_HEADER_GUARD
#define CUSTOM_MCMC_HEADER_GUARD
#include "numericType.h"
#include <boost/random/mersenne_twister.hpp>
#include "cliqueTree.h"
#include "cliqueTreeAdjacencyMatrix.h"
namespace chordalGraph
{
	typedef cliqueTreeAdjacencyMatrix cliqueTreeType;
	typedef moveable_adjacency_matrix<> graphType;
	typedef std::bitset<MAX_STORAGE_VERTICES> bitsetType;
	struct working
	{
	public:
		working(int nVertices)
			: copied(nVertices), copied2(nVertices), counts1(nVertices), counts2(nVertices), colourVector(nVertices)
		{}
		//Working memory for the addition of edges to the clique tree
		std::list<cliqueTreeType::cliqueTreeGraphType::vertex_descriptor> vertexSequence;
		std::list<cliqueTreeType::externalEdge> edgeSequence;
		std::vector<cliqueTreeType::externalEdge> addEdges;
		std::vector<cliqueTreeType::externalEdge> removeEdges;
		cliqueTreeType::unionMinimalSeparatorsTemporaries unionMinimalSepTemp;
		//Two copies of the current clique tree.
		cliqueTreeAdjacencyMatrix copied, copied2;
		//Count vectors used when removing edges 
		std::vector<int> counts1, counts2;
		//Colour vector for a depth first search
		std::vector<boost::default_color_type> colourVector;
		//Vector of edges to be removed from the graph and clique tree. 
		std::vector<int> edgesToRemove;
		//Probability distribution of the number of edges to remove
		std::vector<double> probabilities;
		//The number of different graphs with that number of extra edges removed. 
		std::vector<int> stateCounts;
		//Stack used to form the removal tree. 
		std::vector<cliqueTreeAdjacencyMatrix::stackEntry> stack;
		//Collection of unique subsets
		std::unordered_set<bitsetType> uniqueSubsets;
		//
		cliqueTreeAdjacencyMatrix::formRemovalTreeTemporaries removalTemporaries;
	};
	struct mcmcArgs
	{
	public:
		mcmcArgs(boost::mt19937& randomSource)
			: randomSource(randomSource), trackEdgeCounts(false)
		{}
		//These are the inputs
		std::vector<mpfr_class> approximateCounts;
		//These are the outputs
		std::vector<mpfr_class> estimates;
		int nVertices;
		boost::mt19937& randomSource;

		std::size_t burnIn;
		std::size_t runSize;
		bool trackEdgeCounts;
		std::vector<int> edgeCounts;
	};
	void customMCMC(mcmcArgs& args);
}
#endif
