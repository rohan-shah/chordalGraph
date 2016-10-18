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
			: copied(nVertices), counts1(nVertices), counts2(nVertices), colourVector(nVertices)
		{}
		std::list<cliqueTreeType::cliqueTreeGraphType::vertex_descriptor> vertexSequence;
		std::list<cliqueTreeType::externalEdge> edgeSequence;
		std::vector<cliqueTreeType::externalEdge> addEdges;
		std::vector<cliqueTreeType::externalEdge> removeEdges;
		cliqueTreeType::unionMinimalSeparatorsTemporaries unionMinimalSepTemp;
		cliqueTreeAdjacencyMatrix copied;
		std::vector<int> counts1, counts2;
		std::vector<boost::default_color_type> colourVector;
		std::vector<int> vertexList;
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

		//Here we remove edges
		if(existingEdge.second)
		{
			cliqueTreeAdjacencyMatrix& copied = temp.copied;
			copied.makeCopy(currentTree);

			if(copied.tryRemoveEdge(randomVertex1, randomVertex2, temp.colourVector, temp.counts1, temp.counts2))
			{
				std::vector<int>& vertexList = temp.vertexList;
				vertexList.clear();
				int totalOtherRemovableEdges = 0;
tryRemoveAnotherEdge:
				for(int i = 0; i < nVertices; i++)
				{
					if(temp.counts2[i] == 1 && temp.counts1[i] != 1)
					{
						if(!copied.tryRemoveEdge(randomVertex1, i, temp.colourVector, temp.counts1, temp.counts2))
						{
							throw std::runtime_error("Internal error");
						}
						totalOtherRemovableEdges--;
						vertexList.push_back(i);
						goto tryRemoveAnotherEdge;
					}
				}
				boost::random::uniform_int_distribution<> extraNumberToRemoveDist(0, totalOtherRemovableEdges);
				int extraToRemove = extraNumberToRemoveDist(randomSource);
				mpfr_class acceptanceValue = (exactValues[original_edges] / exactValues[original_edges - 1 - extraToRemove]) * (totalOtherRemovableEdges + 1);
				while((int)vertexList.size() != totalOtherRemovableEdges)
				{
					bitsetType newEdges;
					copied.addEdge(randomVertex1, vertexList.back(), newEdges, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp, false);
				}
				if(acceptanceValue >= 1 || standardUniform(randomSource) <= acceptanceValue.convert_to<double>())
				{
					currentTree.swap(copied);
					boost::remove_edge(randomVertex1, randomVertex2, graph);
					for(int i = 0; i < totalOtherRemovableEdges; i++)
					{
						boost::remove_edge(randomVertex1, vertexList[i], graph);
					}
				}
			}
		}
		//Here we add edges
		else
		{
			bitsetType newEdgesVertex1;
			currentTree.unionMinimalSeparators(randomVertex1, randomVertex2, newEdgesVertex1, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp);
			int increaseInEdges = (int)newEdgesVertex1.count() + 1;
			if((int)original_edges + increaseInEdges <= edgeLimit)
			{
				mpfr_class acceptanceValue = exactValues[original_edges] / exactValues[original_edges + increaseInEdges];
				if(acceptanceValue >= 1 || standardUniform(randomSource) <= acceptanceValue.convert_to<double>())
				{
					currentTree.addEdge(randomVertex1, randomVertex2, newEdgesVertex1, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp, true);
					newEdgesVertex1[randomVertex2] = true;
					for(int i = 0; i < nVertices; i++)
					{
						if(newEdgesVertex1[i])
						{
							boost::add_edge(randomVertex1, i, graph);
						}
					}
				}
			}
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
		exactValues.push_back(33969628);
		/*exactValues.push_back(59194170);
		exactValues.push_back(94169376);
		exactValues.push_back(137060700);
		exactValues.push_back(181199340);
		exactValues.push_back(216312390);
		exactValues.push_back(234891000);
		exactValues.push_back(237142836);*/

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
