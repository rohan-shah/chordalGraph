#include "cliqueTree.h"
#include "cliqueTreeAdjacencyMatrix.h"
#include "numericType.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/random/random_number_generator.hpp>
#include <iostream>
namespace chordalGraph
{
	typedef cliqueTreeAdjacencyMatrix cliqueTreeType;
	typedef moveable_adjacency_matrix<> graphType;
	typedef std::bitset<MAX_STORAGE_VERTICES> bitsetType;
	struct working
	{
	public:
		working(int nVertices)
			: copied(nVertices), copiedGraph(nVertices)
		{}
		std::list<cliqueTreeType::cliqueTreeGraphType::vertex_descriptor> vertexSequence;
		std::list<cliqueTreeType::externalEdge> edgeSequence;
		std::vector<cliqueTreeType::externalEdge> addEdges;
		std::vector<cliqueTreeType::externalEdge> removeEdges;
		cliqueTreeType::unionMinimalSeparatorsTemporaries unionMinimalSepTemp;
		cliqueTreeAdjacencyMatrix copied;
		graphType copiedGraph;
	};
	void step(cliqueTreeType& currentTree, graphType& graph, std::vector<mpfr_class>& exactValues, int nVertices, boost::mt19937& randomSource, working& temp, int edgeLimit)
	{
		boost::random::uniform_int_distribution<> randomVertexDist(0, nVertices-1);
		boost::random::bernoulli_distribution<> standardBernoulli;
		boost::random::uniform_real_distribution<> standardUniform;
		int randomVertex1, randomVertex2;
		do
		{
			randomVertex1 = randomVertexDist(randomSource);
			randomVertex2 = randomVertexDist(randomSource);
		} while(randomVertex1 == randomVertex2);
		std::size_t original_edges = boost::num_edges(graph);
		std::pair<graphType::edge_descriptor, bool> existingEdge = boost::edge(randomVertex1, randomVertex2, graph);
		cliqueTreeAdjacencyMatrix& copied = temp.copied;
		copied.makeCopy(currentTree);

		//Here we remove edges
		if(existingEdge.second)
		{
			copied.reset();
			for(int i = 0; i < nVertices; i++) copied.addVertex();

			temp.copiedGraph.clear();
			temp.copiedGraph.num_vertices = nVertices;

			int proposalEdges = 0;

			int minRandomVertex = std::min(randomVertex1, randomVertex2);
			int maxRandomVertex = std::max(randomVertex1, randomVertex2);
			for(int i = 0; i < nVertices; i++)
			{
				bitsetType added, possible;
				//Work out all the possible extra edges for this vertex
				for(int j = i+1; j < nVertices; j++)
				{
					if(boost::edge(i, j, graph).second)
					{
						possible[j] = true;
					}
				}
				if(i == minRandomVertex) possible[maxRandomVertex] = false;

				for(int j = i+1; j < nVertices; j++)
				{
					if(i == minRandomVertex && j == maxRandomVertex) continue;
					if(added[j]) continue;
					if(boost::edge(i, j, graph).second)
					{
						bitsetType newEdges;
						copied.unionMinimalSeparators(i, j, newEdges, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp);
						if((newEdges & ~possible).none())
						{
							copied.addEdge(i, j, newEdges, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp, true);
							boost::add_edge(i, j, temp.copiedGraph);
							for(int k = i + 1; k < nVertices; k++)
							{
								if(newEdges[k]) boost::add_edge(i, k, temp.copiedGraph);
							}
							added |= newEdges;
							proposalEdges += newEdges.count() + 1;
						}
					}
				}
			}
			
			mpfr_class acceptanceValue = exactValues[original_edges] / exactValues[proposalEdges];

			if(acceptanceValue >= 1 || standardUniform(randomSource) <= acceptanceValue.convert_to<double>())
			{
				currentTree.swap(copied);
				graph.swap(temp.copiedGraph);
#ifdef TRACK_GRAPH
				if(boost::num_edges(currentTree.getGraph()) != boost::num_edges(graph)) throw std::runtime_error("Internal error");
#endif
			}
		}
		//Here we add edges
		else
		{
			bitsetType newEdgesVertex1;
			copied.unionMinimalSeparators(randomVertex1, randomVertex2, newEdgesVertex1, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp);
			int increaseInEdges = (int)newEdgesVertex1.count() + 1;
			if((int)original_edges + increaseInEdges <= edgeLimit)
			{
				mpfr_class acceptanceValue = exactValues[original_edges] / exactValues[original_edges + increaseInEdges];
				copied.addEdge(randomVertex1, randomVertex2, newEdgesVertex1, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp, true);
				if(acceptanceValue >= 1 || standardUniform(randomSource) <= acceptanceValue.convert_to<double>())
				{
					newEdgesVertex1[randomVertex2] = true;
					for(int i = 0; i < nVertices; i++)
					{
						if(newEdgesVertex1[i])
						{
							boost::add_edge(randomVertex1, i, graph);
						}
					}
					currentTree.swap(copied);
				}
			}
#ifdef TRACK_GRAPH
			if(boost::num_edges(currentTree.getGraph()) != boost::num_edges(graph)) throw std::runtime_error("Internal error");
#endif
		}
		currentTree.check();
	}
	void main()
	{
		boost::mt19937 randomSource;
		randomSource.seed(1);

		std::vector<mpfr_class> exactValues;
		exactValues.push_back(1);
		exactValues.push_back(36);
		exactValues.push_back(630);
		exactValues.push_back(7140);
		exactValues.push_back(58527);
		exactValues.push_back(364140);
		exactValues.push_back(1741530);
		exactValues.push_back(6317460);
		exactValues.push_back(16933905);

		mpfr_class exactValuesSum = 0;
		for(std::vector<mpfr_class>::iterator i = exactValues.begin(); i != exactValues.end(); i++) exactValuesSum += *i;

		int nVertices = 9;
		int edgeLimit = exactValues.size()-1;
		cliqueTreeType currentTree(nVertices);
		graphType graph(nVertices);
		for(int i = 0; i < nVertices; i++) currentTree.addVertex();

		working temp(nVertices);
		//burn-in
		for(int i = 0; i < 2000; i++)
		{
			step(currentTree, graph, exactValues, nVertices, randomSource, temp, edgeLimit);
		}
		std::vector<std::size_t> counters(exactValues.size(), 0);
		const int sampleSize = 100000;
		for(int i = 0; i < sampleSize; i++)
		{
			step(currentTree, graph, exactValues, nVertices, randomSource, temp, edgeLimit);
			counters[boost::num_edges(graph)]++;
		}
		for(int j = 0; j < (int)exactValues.size(); j++)
		{
			std::cout << counters[j] << std::endl;
		}
		mpfr_class sumFirstSix = 0;
		for(int i = 0; i < 6; i++) sumFirstSix += mpfr_class(counters[i]) / sampleSize;
		mpfr_class estimate = 6 * exactValues.back() * (mpfr_class(counters[edgeLimit]) / sampleSize) / sumFirstSix;
		std::cout << "Estimate: " << estimate.str(10,std::ios_base::dec) << std::endl;
		std::cout << "Exact: " << exactValues.back().str(10, std::ios_base::dec) << std::endl;

	}
}
int main(int argc, char** argv)
{
	chordalGraph::main();
	return 0;
}
