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
	struct working
	{
	public:
		working(int nVertices)
			: counts1(nVertices), counts2(nVertices), colourVector(nVertices)
		{}
		std::list<cliqueTreeType::cliqueTreeGraphType::vertex_descriptor> vertexSequence;
		std::list<cliqueTreeType::externalEdge> edgeSequence;
		std::vector<cliqueTreeType::externalEdge> addEdges;
		std::vector<cliqueTreeType::externalEdge> removeEdges;
		cliqueTreeType::unionMinimalSeparatorsTemporaries unionMinimalSepTemp;
		std::vector<int> counts1, counts2;
		std::vector<boost::default_color_type> colourVector;
		std::vector<int> vertexList;
	};
	void step(cliqueTreeType& currentTree, graphType& graph, std::vector<mpfr_class>& exactValues, int nVertices, boost::mt19937& randomSource, working& temp, int edgeLimit)
	{
		boost::random::uniform_int_distribution<> randomVertexDist(0, nVertices-1);
		boost::random::bernoulli_distribution<> standardBernoulli;
		boost::random::uniform_real_distribution<> standardUniform;
proposeAnother:
		int randomVertex1, randomVertex2;
		do
		{
			randomVertex1 = randomVertexDist(randomSource);
			randomVertex2 = randomVertexDist(randomSource);
		} while(randomVertex1 == randomVertex2);
		std::size_t original_edges = boost::num_edges(graph);
		std::pair<graphType::edge_descriptor, bool> existingEdge = boost::edge(randomVertex1, randomVertex2, graph);

		//Here we remove edges
		if(existingEdge.second)
		{
			int cliqueVertex = -1;
			if(currentTree.canRemoveEdge(randomVertex1, randomVertex2, temp.counts1, cliqueVertex))
			{
				mpfr_class acceptanceValue = (exactValues[original_edges] / exactValues[original_edges - 1]);
				if(acceptanceValue >= 1 || standardUniform(randomSource) <= acceptanceValue.convert_to<double>())
				{
					currentTree.removeEdgeKnownCliqueVertex(randomVertex1, randomVertex2, temp.colourVector, temp.counts2, cliqueVertex);
					boost::remove_edge(randomVertex1, randomVertex2, graph);
				}
			}
		}
		//Here we add edges
		else
		{
			bitsetType newEdgesVertex1;
			currentTree.unionMinimalSeparators(randomVertex1, randomVertex2, newEdgesVertex1, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp);
			int increaseInEdges = 1;
			for(int i = 0; i < nVertices; i++)
			{
				if(newEdgesVertex1[i] && !boost::edge(i, randomVertex1, graph).second) increaseInEdges++;
			}
			if((int)original_edges + increaseInEdges <= edgeLimit && increaseInEdges == 1)
			{
				mpfr_class acceptanceValue = exactValues[original_edges] / exactValues[original_edges + increaseInEdges];
				if(acceptanceValue >= 1 || standardUniform(randomSource) <= acceptanceValue.convert_to<double>())
				{
					currentTree.addEdge(randomVertex1, randomVertex2, newEdgesVertex1, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp, true);
					boost::add_edge(randomVertex1, randomVertex2, graph);
				}
			}
		}
	}
	void main()
	{
		boost::mt19937 randomSource;
		randomSource.seed(2);

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
		exactValues.push_back(33969628);
		exactValues.push_back(59194170);
		exactValues.push_back(94169376);
		exactValues.push_back(137060700);
		exactValues.push_back(181199340);
		exactValues.push_back(216312390);
		exactValues.push_back(234891000);
		exactValues.push_back(237142836);
		exactValues.push_back(227923920);
		exactValues.push_back(204956724);
		exactValues.push_back(170178120);
		exactValues.push_back(130062807);
		exactValues.push_back(92533764);
		exactValues.push_back(62171838);
		exactValues.push_back(39638592);
		exactValues.push_back(23221338);
		exactValues.push_back(12310704);
		exactValues.push_back(5983866);
		exactValues.push_back(2699508);
		exactValues.push_back(1154547);
		exactValues.push_back(430236);
		exactValues.push_back(137718);
		exactValues.push_back(37800);
		exactValues.push_back(10080);
		exactValues.push_back(2100);
		exactValues.push_back(252);
		exactValues.push_back(36);
		exactValues.push_back(1);

		mpfr_class exactValuesSum = 0;
		for(std::vector<mpfr_class>::iterator i = exactValues.begin(); i != exactValues.end(); i++) exactValuesSum += *i;

		int nVertices = 4;
		int edgeLimit = exactValues.size()-1;
		cliqueTreeType currentTree(nVertices);
		graphType graph(nVertices);
		for(int i = 0; i < nVertices; i++)
		{
			currentTree.addVertex();
			boost::add_vertex(graph);
		}

		{
			std::list<cliqueTreeType::cliqueTreeGraphType::vertex_descriptor> vertexSequence;
			std::list<cliqueTreeType::externalEdge> edgeSequence;
			std::vector<cliqueTreeType::externalEdge> addEdges;
			std::vector<cliqueTreeType::externalEdge> removeEdges;
			cliqueTreeType::unionMinimalSeparatorsTemporaries temp;

			bitsetType minimalSeparator;
			/*for(int i = 0; i < nVertices; i++)
			{
				for(int j = i+1; j < nVertices; j++)
				{
					currentTree.addEdge(i, j, minimalSeparator, vertexSequence, edgeSequence, addEdges, removeEdges, temp, false);
					boost::add_edge(i, j, graph);
				}
			}*/
		}

		working temp(nVertices);
		std::vector<std::size_t> counters(exactValues.size(), 0);
		//burn-in
		for(int i = 0; i < 50000; i++)
		{
			step(currentTree, graph, exactValues, nVertices, randomSource, temp, edgeLimit);
			counters[boost::num_edges(graph)]++;
		}
		const int sampleSize = 1000000;
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

		std::cout << "Estimates: " << std::endl;
		for(int i = 0; i < edgeLimit+1; i++)
		{
			mpfr_class estimate = 6 * exactValues[i] * (mpfr_class(counters[i]) / sampleSize) / sumFirstSix;
			std::cout << "[" << i << "] "  << estimate.str(10,std::ios_base::dec) << std::endl;
		}
	}
}
int main(int argc, char** argv)
{
	chordalGraph::main();
	return 0;
}
