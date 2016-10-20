#include "armstrongMCMC.h"
#include "cliqueTree.h"
#include "cliqueTreeAdjacencyMatrix.h"
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/random/random_number_generator.hpp>
#include <iostream>
namespace chordalGraph
{
	void stepArmstrong(cliqueTreeType& currentTree, graphType& graph, std::vector<mpfr_class>& exactValues, int nVertices, boost::mt19937& randomSource, working& temp, int edgeLimit)
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
	void armstrongMCMC(mcmcArgs& args)
	{
		//Check that we have at least five vertices. Can still run the MCMC with 4, but there's nothing really to estimate, because using the MCMC assumes that the first six count values are known. Which is all of them in the case of four vertices. 
		int nVertices = args.nVertices;
		if(nVertices <= 4)
		{
			throw std::runtime_error("Input nVertices must be at least 5");
		}
		//The case where the edgeLimit is 5 or smaller is not interesting. 
		int edgeLimit = args.approximateCounts.size()-1;
		if(edgeLimit > ((nVertices + 1) * nVertices) /2 || edgeLimit <= 5)
		{
			throw std::runtime_error("Input edgeLimit out of the valid range or smaller than 6");
		}
		//Set up clique tree and graph
		cliqueTreeType currentTree(nVertices);
		graphType graph(nVertices);
		for(int i = 0; i < nVertices; i++) currentTree.addVertex();

		working temp(nVertices);

		std::size_t burnIn = args.burnIn;
		//burn-in
		for(std::size_t i = 0; i < burnIn; i++)
		{
			stepArmstrong(currentTree, graph, args.approximateCounts, nVertices, args.randomSource, temp, edgeLimit);
		}
		std::vector<std::size_t> counters(args.approximateCounts.size(), 0);
		for(std::size_t i = 0; i < args.runSize; i++)
		{
			stepArmstrong(currentTree, graph, args.approximateCounts, nVertices, args.randomSource, temp, edgeLimit);
			counters[boost::num_edges(graph)]++;
		}
		mpfr_class sumFirstSix = 0;
		for(int i = 0; i < 6; i++) sumFirstSix += mpfr_class(counters[i]) / args.runSize;

		args.estimates.resize(args.approximateCounts.size());
		for(std::size_t i = 0; i < args.approximateCounts.size(); i++)
		{
			args.estimates[i] = 6 * args.approximateCounts[i] * (mpfr_class(counters[i]) / args.runSize) / sumFirstSix;
		}
	}
}
