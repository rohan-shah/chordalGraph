#ifndef ARMSTRONG_MCMC_HEADER_GUARD
#define ARMSTRONG_MCMC_HEADER_GUARD
#include "numericType.h"
#include <boost/random/mersenne_twister.hpp>
#include "customMCMC.h"
namespace chordalGraph
{
	struct workingArmstrong
	{
	public:
		workingArmstrong(int nVertices)
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
	void armstrongMCMC(mcmcArgs& args);
}
#endif
