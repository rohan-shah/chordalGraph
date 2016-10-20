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
		std::list<cliqueTreeType::cliqueTreeGraphType::vertex_descriptor> vertexSequence;
		std::list<cliqueTreeType::externalEdge> edgeSequence;
		std::vector<cliqueTreeType::externalEdge> addEdges;
		std::vector<cliqueTreeType::externalEdge> removeEdges;
		cliqueTreeType::unionMinimalSeparatorsTemporaries unionMinimalSepTemp;
		cliqueTreeAdjacencyMatrix copied, copied2;
		std::vector<int> counts1, counts2;
		std::vector<boost::default_color_type> colourVector;
		std::vector<int> vertexList;
		std::vector<double> probabilities;
	};
	struct mcmcArgs
	{
	public:
		mcmcArgs(boost::mt19937& randomSource)
			: randomSource(randomSource)
		{}
		//These are the inputs
		std::vector<mpfr_class> approximateCounts;
		//These are the outputs
		std::vector<mpfr_class> estimates;
		int nVertices;
		boost::mt19937& randomSource;

		std::size_t burnIn;
		std::size_t runSize;
	};
	void customMCMC(mcmcArgs& args);
}
#endif
