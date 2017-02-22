#ifndef POSTERIOR_INFERENCE_HEADER_GUARD
#define POSTERIOR_INFERENCE_HEADER_GUARD
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include "cliqueTree.h"
#include "numericType.h"
#include <boost/random/mersenne_twister.hpp>
#include "customMCMCSymmetric.h"
namespace chordalGraph
{
	struct workingPosteriorInference
	{
	public:
		workingPosteriorInference(int nVertices)
			: nVertices(nVertices), copied(nVertices), copied2(nVertices), copiedTree(nVertices), counts1(nVertices), counts2(nVertices), colourVector(nVertices)
		{
			for(int i = 0; i < nVertices; i++) copiedTree.addVertex();
		}
		int nVertices, delta, deltaStar;
		//Working memory for the addition of edges to the clique tree
		std::list<cliqueTreeType::cliqueTreeGraphType::vertex_descriptor> vertexSequence;
		std::list<cliqueTreeType::externalEdge> edgeSequence;
		std::vector<cliqueTreeType::externalEdge> addEdges;
		std::vector<cliqueTreeType::externalEdge> removeEdges;
		cliqueTreeType::unionMinimalSeparatorsTemporaries unionMinimalSepTemp;
		//Two copies of the current clique tree.
		cliqueTreeAdjacencyMatrix copied, copied2;
		cliqueTreeType copiedTree;
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
		std::vector<std::pair<int, int> > presentEdges, absentEdges;
		boost::numeric::ublas::matrix<double> psi, psiStar, psiPart;
		std::vector<mpfr_class> multivariateGammaDelta, multivariateGammaDeltaStar;
	};
	struct posteriorInferenceArgs
	{
		typedef moveable_adjacency_matrix<> graphType;
		posteriorInferenceArgs()
			: delta(0), dataPoints(0), sampleSize(0), burnIn(0)
		{}
		std::size_t delta, dataPoints, dimension;
		boost::numeric::ublas::matrix<double> sampleCovariance;
		boost::numeric::ublas::matrix<double> psi;
		std::vector<std::pair<graphType, double> > results;
		std::vector<mpfr_class> exactValues;
		std::size_t sampleSize, burnIn;
		boost::mt19937 randomSource;
	};
	void posteriorInference(posteriorInferenceArgs& args);
}
#endif
