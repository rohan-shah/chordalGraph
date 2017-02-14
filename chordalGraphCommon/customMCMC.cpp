#include "customMCMC.h"
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <iostream>
namespace chordalGraph
{
	void stepCustom(cliqueTreeType& currentTree, graphType& graph, std::vector<mpfr_class>& exactValues, int nVertices, boost::mt19937& randomSource, working& temp, int edgeLimit)
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
		//Here we remove edges
		if(existingEdge.second)
		{
			int cliqueVertex = -1;
			//Can we remove this edge?
			if(currentTree.canRemoveEdge(randomVertex1, randomVertex2, temp.counts1, cliqueVertex))
			{
				//Form tree of removable edge subsets and work out the counts. 
				currentTree.formRemovalTree(temp.stateCounts, copied, randomVertex1, randomVertex2, temp.uniqueSubsets, temp.removalTemporaries);
				//The maximum number of other removable edges, in addition to edge (u, v).  
				int maximumOtherRemovableEdges = 0;
				for(int i = 0; i < nVertices; i++) 
				{
					if(temp.stateCounts[i] == 0) break;
					maximumOtherRemovableEdges = i;
				}
				//The number of edges to actually remove for the proposal. 
				int extraToRemove = 0;
				//Compute probabilities and normalizing constant. 
				temp.probabilities.clear();
				double sum1 = 0;
				for(int i = 0; i < maximumOtherRemovableEdges + 1; i++)
				{
					temp.probabilities.push_back(mpfr_class(exactValues[original_edges] / exactValues[original_edges - i - 1]).convert_to<double>());
					sum1 += temp.probabilities.back();
				}
				if(maximumOtherRemovableEdges > 0)
				{
					boost::random::discrete_distribution<> extraNumberToRemoveDist(temp.probabilities.begin(), temp.probabilities.end());
					extraToRemove = extraNumberToRemoveDist(randomSource);
				}
				//The acceptance probability for the Metropolis-Hasting proposal. 
				mpfr_class acceptanceProbability = 0;
				//If we actually only remove a single edge then things are more complicated - There are two ways this can happen, the second way corresponds to swapping randomVertex1 and randomVertex2
				if(extraToRemove == 0)
				{
					//Form tree of removable edge subsets and work out the counts. 
					currentTree.formRemovalTree(temp.stateCounts, copied, randomVertex2, randomVertex1, temp.uniqueSubsets, temp.removalTemporaries);
					
					int maximumOtherRemovableEdges2 = 0;
					for(int i = 0; i < nVertices; i++)
					{
						if(temp.stateCounts[i] == 0) break;
						maximumOtherRemovableEdges2 = i;
					}
					double sum2 = 0;
					for(int i = 0; i < maximumOtherRemovableEdges2 + 1; i++)
					{
						sum2 += mpfr_class(exactValues[original_edges] / exactValues[original_edges - i - 1]).convert_to<double>();
					}
					acceptanceProbability = 1/(0.5 * (1/sum1 + 1/sum2));
				}
				else acceptanceProbability = sum1 * temp.stateCounts[extraToRemove];
				if(acceptanceProbability >= 1 || standardUniform(randomSource) <= acceptanceProbability.convert_to<double>())
				{
					currentTree.removeEdgeKnownCliqueVertex(randomVertex1, randomVertex2, temp.colourVector, temp.counts2, cliqueVertex);
					boost::remove_edge(randomVertex1, randomVertex2, graph);
					if(extraToRemove != 0)
					{
						boost::random::uniform_int_distribution<> randomSubset(0, temp.stateCounts[extraToRemove] - 1);
						int index = randomSubset(randomSource);
						bitsetType chosenSubset;
						for(std::unordered_set<bitsetType>::iterator i = temp.uniqueSubsets.begin(); i != temp.uniqueSubsets.end(); i++)
						{
							if((int)i->count() == extraToRemove+1)
							{
								if(index == 0)
								{
									chosenSubset = *i;
									break;
								}
								index--;
							}
						}
						//This is 1 rather than 0, because the edge randomVertex1, randomVertex2 is already deleted. 
						while(chosenSubset.count() > 1)
						{
							currentTree.canRemoveEdge(randomVertex1, randomVertex2, temp.counts1, cliqueVertex);
							for(int i = 0; i < nVertices; i++)
							{
								if(chosenSubset[i] && temp.counts1[i] == 1)
								{
									chosenSubset[i] = false;
									currentTree.tryRemoveEdge(randomVertex1, i, temp.colourVector, temp.counts2);
									boost::remove_edge(randomVertex1, i, graph);
								}
							}
						}
					}
				}
			}
		}
		//Here we add edges
		else
		{
			cliqueTreeAdjacencyMatrix& copied2 = temp.copied2;
			bitsetType newEdgesVertex1;
			currentTree.unionMinimalSeparators(randomVertex1, randomVertex2, newEdgesVertex1, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp);
			int increaseInEdges = 1;
			for(int i = 0; i < nVertices; i++)
			{
				if(newEdgesVertex1[i] && !boost::edge(i, randomVertex1, graph).second) increaseInEdges++;
			}
			if((int)original_edges + increaseInEdges <= edgeLimit)
			{
				mpfr_class acceptanceProbability;
				if(increaseInEdges != 1)
				{
					copied.makeCopy(currentTree);
					//Actually add the edges to the copy
					copied.addEdge(randomVertex1, randomVertex2, newEdgesVertex1, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp, true);
					copied.formRemovalTree(temp.stateCounts, copied2, randomVertex1, randomVertex2, temp.uniqueSubsets, temp.removalTemporaries);
					//Work out how many more edges can be removed from the original
					int otherBackwardsCanRemove = 0;
					for(int i = 0; i < nVertices; i++) 
					{
						if(temp.stateCounts[i] == 0) break;
						otherBackwardsCanRemove = i;
					}

					double sum = 0;
					for(int i = 0; i < otherBackwardsCanRemove + 1; i++)
					{
						sum += mpfr_class(exactValues[original_edges + increaseInEdges] / exactValues[original_edges + increaseInEdges - i - 1]).convert_to<double>();
					}
					acceptanceProbability = 1/(sum * temp.stateCounts[increaseInEdges - 1]);
				}
				else
				{
					copied.makeCopy(currentTree);
					//Actually add the edges to the copy
					copied.addEdge(randomVertex1, randomVertex2, newEdgesVertex1, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp, true);
					//Form removal tree with the vertices the same way round
					copied.formRemovalTree(temp.stateCounts, copied2, randomVertex1, randomVertex2, temp.uniqueSubsets, temp.removalTemporaries);
					//Work out how many more edges can be removed from the original
					int otherBackwardsCanRemove = 0;
					for(int i = 0; i < nVertices; i++) 
					{
						if(temp.stateCounts[i] == 0) break;
						otherBackwardsCanRemove = i;
					}

					double sum1 = 0;
					for(int i = 0; i < otherBackwardsCanRemove + 1; i++)
					{
						sum1 += mpfr_class(exactValues[original_edges + increaseInEdges] / exactValues[original_edges + increaseInEdges - i - 1]).convert_to<double>();
					}
					//Form removal tree with the vertices reversed
					copied.formRemovalTree(temp.stateCounts, copied2, randomVertex2, randomVertex1, temp.uniqueSubsets, temp.removalTemporaries);
					//Work out how many more edges can be removed from the original
					int otherBackwardsCanRemove2 = 0;
					for(int i = 0; i < nVertices; i++) 
					{
						if(temp.stateCounts[i] == 0) break;
						otherBackwardsCanRemove2 = i;
					}

					double sum2 = 0;
					for(int i = 0; i < otherBackwardsCanRemove2 + 1; i++)
					{
						sum2 += mpfr_class(exactValues[original_edges + increaseInEdges] / exactValues[original_edges + increaseInEdges - i - 1]).convert_to<double>();
					}
					acceptanceProbability = 0.5 * (1/sum1 + 1/sum2);
				}
				if (acceptanceProbability >= 1 || standardUniform(randomSource) <= acceptanceProbability.convert_to<double>())
				{
					currentTree.swap(copied);
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
#ifndef NDEBUG
		currentTree.check();
#endif
	}
	void customMCMC(mcmcArgs& args)
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
			stepCustom(currentTree, graph, args.approximateCounts, nVertices, args.randomSource, temp, edgeLimit);
		}
		std::vector<std::size_t> counters(args.approximateCounts.size(), 0);
		if(args.trackEdgeCounts) args.edgeCounts.clear();
		for(std::size_t i = 0; i < args.runSize; i++)
		{
			stepCustom(currentTree, graph, args.approximateCounts, nVertices, args.randomSource, temp, edgeLimit);
			std::size_t nEdges = boost::num_edges(graph);
			counters[nEdges]++;
			if(args.trackEdgeCounts) args.edgeCounts.push_back((int)nEdges);
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
