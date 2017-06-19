#include "customSymmetricPosteriorInference.h"
#include <boost/graph/depth_first_search.hpp>
#include "cliqueTreeAdjacencyMatrix.h"
#include "customMCMCSymmetric.h"
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/random/discrete_distribution.hpp>
#include "posteriorInference.h"
namespace chordalGraph
{
	void posteriorInferenceStepRemoveEdges(cliqueTreeAdjacencyMatrix& currentTree, customSymmetricPosteriorInferenceArgs::graphType& graph, std::vector<mpfr_class>& exactValues, boost::mt19937& randomSource, workingCustomSymmetricPosteriorInference& working)
	{
		boost::random::bernoulli_distribution<> standardBernoulli;
		boost::random::uniform_real_distribution<> standardUniform;
		int original_edges = (int)boost::num_edges(graph);
		int nVertices = working.nVertices;
		int maxEdges = (nVertices * (nVertices - 1)) / 2;

		cliqueTreeAdjacencyMatrix& copiedTree = working.copiedTree;
		int delta = working.delta;
		int deltaStar = working.deltaStar;

		boost::random::uniform_int_distribution<> randomEdgeDist(0, (int)working.presentEdges.size()-1);
		int edgeIndex = randomEdgeDist(randomSource);
		int randomVertex1 = working.presentEdges[edgeIndex].first, randomVertex2 = working.presentEdges[edgeIndex].second;
		//Standardise order
		if(standardBernoulli(randomSource)) std::swap(randomVertex1, randomVertex2);

		int originalCliqueVertex;
		//Can we remove this edge?
		if(currentTree.canRemoveEdge(randomVertex1, randomVertex2, working.counts1, originalCliqueVertex))
		{
			//Form tree of removable edge subsets and work out the counts. 
			currentTree.formRemovalTree(working.stateCounts, copiedTree, randomVertex1, randomVertex2, working.uniqueSubsets, working.removalTemporaries);
			//The maximum number of other removable edges, in addition to edge (u, v).  
			int maximumOtherRemovableEdges = 0;
			for(int i = 0; i < nVertices; i++) 
			{
				if(working.stateCounts[i] == 0) break;
				maximumOtherRemovableEdges = i;
			}
			//The number of edges to actually remove for the proposal. 
			int extraToRemove = 0;
			//Compute probabilities and normalizing constant. 
			working.probabilities.clear();
			double sum1 = 0;
			for(int i = 0; i < maximumOtherRemovableEdges + 1; i++)
			{
				working.probabilities.push_back(mpfr_class(exactValues[original_edges] / exactValues[original_edges - i - 1]).convert_to<double>());
				sum1 += working.probabilities.back();
			}
			if(maximumOtherRemovableEdges > 0)
			{
				boost::random::discrete_distribution<> extraNumberToRemoveDist(working.probabilities.begin(), working.probabilities.end());
				extraToRemove = extraNumberToRemoveDist(randomSource);
			}
			//The acceptance probability for the Metropolis-Hasting proposal. 
			mpfr_class acceptanceProbability = 0;
			//If we actually only remove a single edge then things are more complicated - There are two ways this can happen, the second way corresponds to swapping randomVertex1 and randomVertex2
			if(extraToRemove == 0)
			{
				//Form tree of removable edge subsets and work out the counts. 
				currentTree.formRemovalTree(working.stateCounts, copiedTree, randomVertex2, randomVertex1, working.uniqueSubsets, working.removalTemporaries);
				
				int maximumOtherRemovableEdges2 = 0;
				for(int i = 0; i < nVertices; i++)
				{
					if(working.stateCounts[i] == 0) break;
					maximumOtherRemovableEdges2 = i;
				}
				double sum2 = 0;
				for(int i = 0; i < maximumOtherRemovableEdges2 + 1; i++)
				{
					sum2 += mpfr_class(exactValues[original_edges] / exactValues[original_edges - i - 1]).convert_to<double>();
				}
				if(original_edges == 1)
				{
					acceptanceProbability = (2.0 / maxEdges) * (exactValues[1] / exactValues[0]);
				}
				else if(original_edges == maxEdges)
				{
					acceptanceProbability = (double)original_edges/(double)((maxEdges - (int)(original_edges - 1)) * (1/sum1 + 1/sum2));
				}
				else
				{
					acceptanceProbability = 2 * (double)original_edges/(double)((maxEdges - (int)(original_edges - 1)) * (1/sum1 + 1/sum2));
				}
			}
			else 
			{
				if(original_edges == maxEdges)
				{
					acceptanceProbability = 0.5*((double)original_edges/(double)(maxEdges - (int)(original_edges - 1 - extraToRemove))) * sum1 * working.stateCounts[extraToRemove];
				}
				else
				{
					acceptanceProbability = ((double)original_edges/(double)(maxEdges - (int)(original_edges - 1 - extraToRemove))) * sum1 * working.stateCounts[extraToRemove];
				}
			}
			//If we're removing more than one, edge, work out which edges we're actually erasing. 
			bitsetType chosenSubset;
			if(extraToRemove != 0)
			{
				boost::random::uniform_int_distribution<> randomSubset(0, working.stateCounts[extraToRemove] - 1);
				int indexWithinThatEdgeCount = randomSubset(randomSource);
				for(std::unordered_set<bitsetType>::iterator i = working.uniqueSubsets.begin(); i != working.uniqueSubsets.end(); i++)
				{
					if((int)i->count() == extraToRemove+1)
					{
						if(indexWithinThatEdgeCount == 0)
						{
							chosenSubset = *i;
							break;
						}
						indexWithinThatEdgeCount--;
					}
				}
			}
			//Fast computation of the ratio of the h functions
			mpfr_class extraFactor2 = 1;
			{
				copiedTree.makeCopy(currentTree);
				extraFactor2 = getHRatio(copiedTree, originalCliqueVertex, randomVertex1, randomVertex2, working.psi, working.psiPart, nVertices, delta) / getHRatio(copiedTree, originalCliqueVertex, randomVertex1, randomVertex2, working.psiStar, working.psiPart, nVertices, deltaStar);
				cliqueTreeAdjacencyMatrix::removeReversal reverse;
				copiedTree.removeEdgeKnownCliqueVertex(randomVertex1, randomVertex2, working.colourVector, working.counts2, originalCliqueVertex, reverse);
				if(extraToRemove != 0)
				{
					//This is 1 rather than 0, because the edge randomVertex1, randomVertex2 is already deleted. 
					bitsetType copiedChosenSubset = chosenSubset;
					while(copiedChosenSubset.count() > 1)
					{
						for(int i = 0; i < nVertices; i++)
						{
							if(copiedChosenSubset[i] && working.counts2[i] == 1)
							{
								int newCliqueVertex;
								copiedTree.canRemoveEdge(randomVertex1, i, working.counts2, newCliqueVertex);
								copiedChosenSubset[i] = false;
								extraFactor2 *= getHRatio(copiedTree, newCliqueVertex, randomVertex1, i, working.psi, working.psiPart, nVertices, delta) / getHRatio(copiedTree, newCliqueVertex, randomVertex1, i, working.psiStar, working.psiPart, nVertices, deltaStar);
								cliqueTreeAdjacencyMatrix::removeReversal reverse;
								copiedTree.tryRemoveEdge(randomVertex1, i, working.colourVector, working.counts2, reverse);
							}
						}
					}
				}
			}
#ifndef NDEBUG
			//Slow computation of the ratio of the h functions
			mpfr_class extraFactor = (h(copiedTree, working.delta, working.psi, working.multivariateGammaDelta, working.psiPart, nVertices, working.colourVector) * h(currentTree, deltaStar, working.psiStar, working.multivariateGammaDeltaStar, working.psiPart, nVertices, working.colourVector)) / (h(currentTree, delta, working.psi, working.multivariateGammaDelta, working.psiPart, nVertices, working.colourVector) * h(copiedTree, deltaStar, working.psiStar, working.multivariateGammaDeltaStar, working.psiPart, nVertices, working.colourVector));
			if(std::fabs(1 - (extraFactor * extraFactor2).convert_to<double>()) > 1e-5)
			{
				throw std::runtime_error("Internal error");
			}
#endif
			acceptanceProbability /= extraFactor2;
			if(acceptanceProbability >= 1 || standardUniform(randomSource) <= acceptanceProbability.convert_to<double>())
			{
				currentTree.swap(copiedTree);
				boost::remove_edge(randomVertex1, randomVertex2, graph);
				{
					std::swap(working.presentEdges[edgeIndex], working.presentEdges.back());
					working.presentEdges.pop_back();
					int minVertex = std::min(randomVertex1, randomVertex2), maxVertex = std::max(randomVertex1, randomVertex2);
					working.absentEdges.push_back(std::make_pair(minVertex, maxVertex));
				}
				if(extraToRemove != 0)
				{
					for(int i = 0; i < nVertices; i++)
					{
						if(chosenSubset[i] && i != randomVertex1 && i != randomVertex2)
						{
							boost::remove_edge(randomVertex1, i, graph);
							{
								int minVertex = std::min(randomVertex1, i), maxVertex = std::max(randomVertex1, i);
								int edgeCounter = (int)std::distance(working.presentEdges.begin(), std::find(working.presentEdges.begin(), working.presentEdges.end(), std::make_pair(minVertex, maxVertex)));
								std::swap(working.presentEdges[edgeCounter], working.presentEdges.back());
								working.presentEdges.pop_back();
								working.absentEdges.push_back(std::make_pair(minVertex, maxVertex));
							}
						}
					}
				}
			}
		}
	}
	void posteriorInferenceStepAddEdges(cliqueTreeAdjacencyMatrix& currentTree, customSymmetricPosteriorInferenceArgs::graphType& graph, std::vector<mpfr_class>& exactValues, boost::mt19937& randomSource, workingCustomSymmetricPosteriorInference& working)
	{
		boost::random::bernoulli_distribution<> standardBernoulli;
		boost::random::uniform_real_distribution<> standardUniform;
		int original_edges = (int)boost::num_edges(graph);
		int nVertices = working.nVertices;
		int maxEdges = (nVertices * (nVertices - 1)) / 2;

		cliqueTreeAdjacencyMatrix& copiedTree = working.copiedTree;
		int delta = working.delta;
		int deltaStar = working.deltaStar;

		boost::random::uniform_int_distribution<> randomEdgeDist(0, (int)working.absentEdges.size()-1);
		int edgeIndex = randomEdgeDist(randomSource);
		int randomVertex1 = working.absentEdges[edgeIndex].first, randomVertex2 = working.absentEdges[edgeIndex].second;
		if(standardBernoulli(randomSource)) std::swap(randomVertex1, randomVertex2);

		cliqueTreeAdjacencyMatrix& copied2 = working.copiedTree2;
		bitsetType newEdgesVertex1;
		currentTree.unionMinimalSeparators(randomVertex1, randomVertex2, newEdgesVertex1, working.vertexSequence, working.edgeSequence, working.addEdges, working.removeEdges, working.unionMinimalSepTemp);
		std::list<cliqueTreeType::cliqueTreeGraphType::vertex_descriptor> vertexSequenceCopied = working.vertexSequence;
		int increaseInEdges = 1;
		for(int i = 0; i < nVertices; i++)
		{
			if(newEdgesVertex1[i] && !boost::edge(i, randomVertex1, graph).second) increaseInEdges++;
		}
		if(original_edges + increaseInEdges <= maxEdges)
		{
			mpfr_class acceptanceProbability;
			if(increaseInEdges != 1)
			{
				copiedTree.makeCopy(currentTree);
				//Actually add the edges to the copy
				copiedTree.addEdge(randomVertex1, randomVertex2, newEdgesVertex1, working.vertexSequence, working.edgeSequence, working.addEdges, working.removeEdges, working.unionMinimalSepTemp, true);
				copiedTree.formRemovalTree(working.stateCounts, copied2, randomVertex1, randomVertex2, working.uniqueSubsets, working.removalTemporaries);
				//Work out how many more edges can be removed from the original
				int otherBackwardsCanRemove = 0;
				for(int i = 0; i < nVertices; i++) 
				{
					if(working.stateCounts[i] == 0) break;
					otherBackwardsCanRemove = i;
				}

				double sum = 0;
				for(int i = 0; i < otherBackwardsCanRemove + 1; i++)
				{
					sum += mpfr_class(exactValues[original_edges + increaseInEdges] / exactValues[original_edges + increaseInEdges - i - 1]).convert_to<double>();
				}
				if(increaseInEdges + original_edges == maxEdges)
				{
					acceptanceProbability = 2.0/(((double)(original_edges+increaseInEdges)/(double)(maxEdges - (int)original_edges)) * sum * working.stateCounts[increaseInEdges - 1]);
				}
				else
				{
					acceptanceProbability = 1.0/(((double)(original_edges+increaseInEdges)/(double)(maxEdges - (int)original_edges)) * sum * working.stateCounts[increaseInEdges - 1]);
				}
			}
			else
			{
				copiedTree.makeCopy(currentTree);
				//Actually add the edges to the copy
				copiedTree.addEdge(randomVertex1, randomVertex2, newEdgesVertex1, working.vertexSequence, working.edgeSequence, working.addEdges, working.removeEdges, working.unionMinimalSepTemp, true);
				//Form removal tree with the vertices the same way round
				copiedTree.formRemovalTree(working.stateCounts, copied2, randomVertex1, randomVertex2, working.uniqueSubsets, working.removalTemporaries);
				//Work out how many more edges can be removed from the original
				int otherBackwardsCanRemove = 0;
				for(int i = 0; i < nVertices; i++) 
				{
					if(working.stateCounts[i] == 0) break;
					otherBackwardsCanRemove = i;
				}

				double sum1 = 0;
				for(int i = 0; i < otherBackwardsCanRemove + 1; i++)
				{
					sum1 += mpfr_class(exactValues[original_edges + increaseInEdges] / exactValues[original_edges + increaseInEdges - i - 1]).convert_to<double>();
				}
				//Form removal tree with the vertices reversed
				copiedTree.formRemovalTree(working.stateCounts, copied2, randomVertex2, randomVertex1, working.uniqueSubsets, working.removalTemporaries);
				//Work out how many more edges can be removed from the original
				int otherBackwardsCanRemove2 = 0;
				for(int i = 0; i < nVertices; i++) 
				{
					if(working.stateCounts[i] == 0) break;
					otherBackwardsCanRemove2 = i;
				}

				double sum2 = 0;
				for(int i = 0; i < otherBackwardsCanRemove2 + 1; i++)
				{
					sum2 += mpfr_class(exactValues[original_edges + increaseInEdges] / exactValues[original_edges + increaseInEdges - i - 1]).convert_to<double>();
				}
				if(original_edges == 0)
				{
					acceptanceProbability = (maxEdges / 2.0) * (exactValues[0] / exactValues[1]);
				}
				else if(original_edges == maxEdges - 1)
				{
					acceptanceProbability = ((double)(maxEdges - original_edges)/(original_edges+1.0)) * (1/sum1 + 1/sum2);
				}
				else
				{
					acceptanceProbability = 0.5*((double)(maxEdges - original_edges)/(original_edges+1.0)) * (1/sum1 + 1/sum2);
				}
			}
			mpfr_class extraFactor2 = 1;
			cliqueTreeAdjacencyMatrix& copiedTree2 = working.copiedTree2;
			copiedTree2.makeCopy(currentTree);
			std::list<cliqueTreeType::cliqueTreeGraphType::vertex_descriptor>::reverse_iterator i = vertexSequenceCopied.rbegin();
			if(i == vertexSequenceCopied.rend())
			{
				//In this case we're joining two separate connected components. 
				int cliqueVertex = -1;
				//Find clique tree node that contains these two vertices. There can only be one.
				for(int otherVertex2 = 0; otherVertex2 < nVertices; otherVertex2++)
				{
					bitsetType contents = boost::get(boost::vertex_name, copiedTree.getCliqueGraph(), otherVertex2).contents;
					if(contents[randomVertex1] && contents[randomVertex2])
					{
						cliqueVertex = otherVertex2; break;
					}
				}
				extraFactor2 *= getHRatio(copiedTree, cliqueVertex, randomVertex1, randomVertex2, working.psi, working.psiPart, nVertices, delta) / getHRatio(copiedTree, cliqueVertex, randomVertex1, randomVertex2, working.psiStar, working.psiPart, nVertices, deltaStar);
			}
			else
			{
				bitsetType tmpNoExtraEdges;
				bitsetType previousVertexContents = boost::get(boost::vertex_name, currentTree.getCliqueGraph(), *i).contents;
				bitsetType stillPossible = previousVertexContents;
				stillPossible.flip();
				i++;
				for(; i != vertexSequenceCopied.rend(); i++)
				{
					bitsetType newVertices = boost::get(boost::vertex_name, currentTree.getCliqueGraph(), *i).contents;
					if(newVertices[randomVertex2])
					{
						newVertices &= previousVertexContents & stillPossible;
						newVertices[randomVertex2] = true;
					}
					else newVertices &= previousVertexContents & stillPossible;
					//This is the vertex to which we are going to add an edge
					for(int otherVertex = 0; otherVertex < nVertices; otherVertex++)
					{
						//If we're actually adding the final edge {randomVertex1, randomVertex2}, then this has to be done *last*. Obviously. Everything else is a pre-requisite for being to add that edge, after all.
						if(newVertices[otherVertex] && otherVertex != randomVertex2)
						{
							copiedTree2.addEdge(randomVertex1, otherVertex, tmpNoExtraEdges, working.vertexSequence2, working.edgeSequence, working.addEdges, working.removeEdges, working.unionMinimalSepTemp, false);
							int cliqueVertex = -1;
							//Find clique tree node that contains these two vertices. There can only be one.
							for(int otherVertex2 = 0; otherVertex2 < nVertices; otherVertex2++)
							{
								bitsetType contents = boost::get(boost::vertex_name, copiedTree2.getCliqueGraph(), otherVertex2).contents;
								if(contents[randomVertex1] && contents[otherVertex])
								{
									cliqueVertex = otherVertex2; break;
								}
							}
							//Now that we know the clique vertex and the edge we're adding, compute the contribution from this change. 
							extraFactor2 *= getHRatio(copiedTree2, cliqueVertex, randomVertex1, otherVertex, working.psi, working.psiPart, nVertices, delta) / getHRatio(copiedTree2, cliqueVertex, randomVertex1, otherVertex, working.psiStar, working.psiPart, nVertices, deltaStar);
						}
					}
					//Add edge {randomVertex1, randomVertex2} last, if it is in fact being added this iteration.
					if(newVertices[randomVertex2])
					{
						copiedTree2.addEdge(randomVertex1, randomVertex2, tmpNoExtraEdges, working.vertexSequence2, working.edgeSequence, working.addEdges, working.removeEdges, working.unionMinimalSepTemp, false);
						int cliqueVertex = -1;
						//Find clique tree node that contains these two vertices. There can only be one.
						for(int otherVertex2 = 0; otherVertex2 < nVertices; otherVertex2++)
						{
							bitsetType contents = boost::get(boost::vertex_name, copiedTree2.getCliqueGraph(), otherVertex2).contents;
							if(contents[randomVertex1] && contents[randomVertex2])
							{
								cliqueVertex = otherVertex2; break;
							}
						}
						//Now that we know the clique vertex and the edge we're adding, compute the contribution from this change. 
						extraFactor2 *= getHRatio(copiedTree2, cliqueVertex, randomVertex1, randomVertex2, working.psi, working.psiPart, nVertices, delta) / getHRatio(copiedTree2, cliqueVertex, randomVertex1, randomVertex2, working.psiStar, working.psiPart, nVertices, deltaStar);
					}
					previousVertexContents = boost::get(boost::vertex_name, currentTree.getCliqueGraph(), *i).contents;
					newVertices.flip();
					stillPossible &= newVertices;
				}
			}
			//Slow computation of the ratio of the h functions
#ifndef NDEBUG
			mpfr_class extraFactor = (h(copiedTree, working.delta, working.psi, working.multivariateGammaDelta, working.psiPart, nVertices, working.colourVector) * h(currentTree, deltaStar, working.psiStar, working.multivariateGammaDeltaStar, working.psiPart, nVertices, working.colourVector)) / (h(currentTree, delta, working.psi, working.multivariateGammaDelta, working.psiPart, nVertices, working.colourVector) * h(copiedTree, deltaStar, working.psiStar, working.multivariateGammaDeltaStar, working.psiPart, nVertices, working.colourVector));
			if(std::fabs(1 - (extraFactor / extraFactor2).convert_to<double>()) > 1e-5)
			{
				throw std::runtime_error("Internal error");
			}
#endif
			acceptanceProbability *= extraFactor2;
			if (acceptanceProbability >= 1 || standardUniform(randomSource) <= acceptanceProbability.convert_to<double>())
			{
				currentTree.swap(copiedTree);
				newEdgesVertex1[randomVertex2] = true;
				for(int i = 0; i < nVertices; i++)
				{
					if(newEdgesVertex1[i] && i != randomVertex1)
					{
						bool addedNewEdge = boost::add_edge(randomVertex1, i, graph).second;
						if(addedNewEdge)
						{
							int minVertex = std::min(randomVertex1, i), maxVertex = std::max(randomVertex1, i);
							int edgeCounter = (int)std::distance(working.absentEdges.begin(), std::find(working.absentEdges.begin(), working.absentEdges.end(), std::make_pair(minVertex, maxVertex)));
							std::swap(working.absentEdges[edgeCounter], working.absentEdges.back());
							working.absentEdges.pop_back();
							working.presentEdges.push_back(std::make_pair(minVertex, maxVertex));
						}
					}
				}
			}
		}
	}
	void posteriorInferenceStep(cliqueTreeAdjacencyMatrix& currentTree, customSymmetricPosteriorInferenceArgs::graphType& graph, std::vector<mpfr_class>& exactValues, boost::mt19937& randomSource, workingCustomSymmetricPosteriorInference& working)
	{
		boost::random::bernoulli_distribution<> standardBernoulli;
		int original_edges = (int)boost::num_edges(graph);
		int nVertices = working.nVertices;
		int maxEdges = (nVertices * (nVertices - 1)) / 2;

		//Here we remove edges
		if(original_edges == maxEdges || (original_edges > 0 && standardBernoulli(randomSource)))
		{
			posteriorInferenceStepRemoveEdges(currentTree, graph, exactValues, randomSource, working);
		}
		//Here we add edges
		else
		{
			posteriorInferenceStepAddEdges(currentTree, graph, exactValues, randomSource, working);
		}
#ifndef NDEBUG
		currentTree.check();
		assert(working.presentEdges.size() == boost::num_edges(graph));
		assert((int)working.presentEdges.size() + (int)working.absentEdges.size() == maxEdges);
#endif
	}
	void customSymmetricPosteriorInference(customSymmetricPosteriorInferenceArgs& args)
	{
		if(!((args.sampleSize > 0) ^ (args.uniqueGraphsLimit > 0)))
		{
			throw std::runtime_error("Only one of sampliSize and uniqueGraphsLimit can be specified");
		}
		cliqueTreeAdjacencyMatrix currentTree(args.dimension);
		customSymmetricPosteriorInferenceArgs::graphType graph(args.dimension);
		for(std::size_t i = 0; i < args.dimension; i++) currentTree.addVertex();

		workingCustomSymmetricPosteriorInference working(args.dimension);
		for(int i = 0; i < (int)args.dimension; i++)
		{
			for(int j = i + 1; j < (int)args.dimension; j++) working.absentEdges.push_back(std::make_pair(i, j));
		}
		working.nVertices = args.dimension;
		working.delta = args.delta;
		working.deltaStar = args.delta + args.dataPoints - 1;
		working.psi = args.psi;
		working.psiStar = args.psi + args.sampleCovariance;
		working.psiPart.resize(args.dimension, args.dimension, false);
		working.multivariateGammaDelta.resize(args.dimension + 1);
		working.multivariateGammaDeltaStar.resize(args.dimension + 1);
		
		for(std::size_t i = 0; i < args.burnIn; i++)
		{
			posteriorInferenceStep(currentTree, graph, args.exactValues, args.randomSource, working);
		}
		args.results.clear();
		if(args.sampleSize > 0)
		{
			for(std::size_t i = 0; i < args.sampleSize; i++)
			{
				posteriorInferenceStep(currentTree, graph, args.exactValues, args.randomSource, working);
				customSymmetricPosteriorInferenceArgs::resultsType::iterator existingGraph = args.results.find(graph);
				if(existingGraph == args.results.end())
				{
					args.results[graph] = 1;
				}
				else existingGraph->second++;
			}
		}
		else
		{
			args.sampleSize = 0;
			while(args.results.size() != args.uniqueGraphsLimit)
			{
				posteriorInferenceStep(currentTree, graph, args.exactValues, args.randomSource, working);
				customSymmetricPosteriorInferenceArgs::resultsType::iterator existingGraph = args.results.find(graph);
				if(existingGraph == args.results.end())
				{
					args.results[graph] = 1;
				}
				else existingGraph->second++;
				args.sampleSize++;
			}
		}
	}

}
