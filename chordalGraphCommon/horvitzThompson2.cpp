#include "horvitzThompson2.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions.hpp>
#include "nauty.h"
#include "cliqueTree.h"
#include "cliqueTreeAdjacencyMatrix.h"
#include "sampford.h"
#include "childNode.h"
namespace chordalGraph
{
	template<typename cliqueTree> void horvitzThompson2(horvitzThompsonArgs& args)
	{
		args.exact = true;
		args.minimumSizeForExact = -1;
		args.estimate = 0;

		//Temporary data that's used in cliqueTree calls
		typename cliqueTree::unionMinimalSeparatorsTemporaries temp;

		boost::random_number_generator<boost::mt19937> generator(args.randomSource);
		boost::random::bernoulli_distribution<> standardBernoulli(0.5);
		//Number of edges either present (or to be added later)
		std::vector<int> nEdges(args.budget);
		std::vector<horvitzThompson2Private::weightedCliqueTree<cliqueTree> > cliqueTrees;
		cliqueTrees.reserve(args.budget);

		//At any point we will have a bunch of conditions on what other edges are required to be present (in order to maintain chordality)
		//and what edges have already been marked as present
		//This vector tracks those conditions. Entry 0 contains the conditions for the first sample, then the following entry contains
		//the conditions for the next sample, etc
		std::vector<bitsetType> conditions(args.budget);

		//The vertex we're currently considering
		int currentVertex = 0;
		//The edge we're currently considering
		int currentEdge = 0;
		//We start off with one sample
		{
			horvitzThompson2Private::weightedCliqueTree<cliqueTree> initialTree(args.nVertices);
			initialTree.tree.addVertex();
			cliqueTrees.push_back(initialTree);
			nEdges[0] = 0;
		}
		//It has no conditions
		conditions[0] = 0;

		//The number of total edges, if the new edge is added. We only need one entry
		//per sample, not two, because if the new edge is not added, then the number of edges
		//stays the same (and this is recorded in nEdges)
		std::vector<int> possibilityEdges(args.budget);

		//The conditions, for the new set of samples
		std::vector<bitsetType> newConditions(args.budget);

		std::vector<horvitzThompson2Private::weightedCliqueTree<cliqueTree> > newCliqueTrees;
		newCliqueTrees.reserve(args.budget);
		std::vector<int> newNEdges(args.budget);

		std::vector<std::list<typename cliqueTree::cliqueTreeGraphType::vertex_descriptor> > vertexSequence(args.budget);
		std::vector<std::list<typename cliqueTree::externalEdge> > edgeSequence(args.budget);
		std::vector<std::vector<typename cliqueTree::externalEdge> > removeEdges(args.budget);
		std::vector<std::vector<typename cliqueTree::externalEdge> > addEdges(args.budget);

		std::vector<bitsetType> unionMinimalSeparators(args.budget);
		//Vector used to shuffle indices
		typedef horvitzThompson2Private::childNode<mpfr_class> childNodeType;
		std::vector<childNodeType> childNodes;
		childNodes.reserve(args.budget*2);

		sampling::sampfordFromParetoNaiveArgs samplingArgs;
		//Inclusion probabilities and indices that result from calling one of the sampling functions. These are taken out of the relevant argument struct
		std::vector<int>& indices = samplingArgs.indices;
		std::vector<numericType>& weights = samplingArgs.weights;
		std::vector<numericType>& rescaledWeights = samplingArgs.rescaledWeights;

		//Nauty variables
		std::vector<int> lab, ptn, orbits;
		lab.reserve(args.nVertices);
		ptn.reserve(args.nVertices);
		std::vector<graph> nautyGraph;
		boost::numeric::ublas::matrix<graph, boost::numeric::ublas::row_major> cannonicalNautyGraphs(2 * args.budget, args.nVertices * SETWORDSNEEDED(args.nVertices));

		//Used to count the number of distinct graphs (up to isomorphism)
		std::vector<bool> alreadyConsidered(2*args.budget);
		//Continue while there are samples left.
		while (currentVertex < args.nVertices - 1)
		{
			currentVertex++;
			currentEdge = 0;
			//There are no conditions when we move to the next vertex
			std::fill(conditions.begin(), conditions.end(), 0);
			//Add the extra vertex
			std::for_each(cliqueTrees.begin(), cliqueTrees.end(), std::mem_fun_ref(&horvitzThompson2Private::weightedCliqueTree<cliqueTree>::addVertex));
			while(currentEdge < currentVertex)
			{
				//Clear vector of indices of possibilities
				childNodes.clear();
				//Remaining edges, including the one we're just about to consider. 
				int nRemainingEdges = currentVertex - currentEdge + ((args.nVertices - currentVertex - 1)* (args.nVertices - 2 - currentVertex) / 2) + (args.nVertices - currentVertex - 1) * (currentVertex + 1);
				//Clear data structures;
				std::for_each(vertexSequence.begin(), vertexSequence.end(), std::mem_fun_ref(&std::list<typename cliqueTree::cliqueTreeGraphType::vertex_descriptor>::clear));
				std::for_each(edgeSequence.begin(), edgeSequence.end(), std::mem_fun_ref(&std::list<typename cliqueTree::externalEdge>::clear));
				std::for_each(removeEdges.begin(), removeEdges.end(), std::mem_fun_ref(&std::vector<typename cliqueTree::externalEdge>::clear));
				std::for_each(addEdges.begin(), addEdges.end(), std::mem_fun_ref(&std::vector<typename cliqueTree::externalEdge>::clear));
				std::fill(unionMinimalSeparators.begin(), unionMinimalSeparators.end(), 0);
				for (int sampleCounter = 0; sampleCounter < (int)cliqueTrees.size(); sampleCounter++)
				{
					horvitzThompson2Private::weightedCliqueTree<cliqueTree>& currentCliqueTree = cliqueTrees[sampleCounter];
					bitsetType copiedConditions = conditions[sampleCounter];
					//maximum possible number of edges
					int maxEdges = nEdges[sampleCounter] + nRemainingEdges - (int)(copiedConditions & ~bitsetType((1ULL << (currentEdge+1)) - 1)).count();
					//Do we need this edge to make up the numbers?
					bool requiresEdge = maxEdges == args.nEdges;
					if (maxEdges < args.nEdges)
					{
						//Too few edges to reach the target number. No children. 
					}
					//The first condition corresponds to the case where adding ALL edges past this point gives us the target number
					//The second corresponds to the case where we've already hit the target
					else if ((currentEdge == 0 && requiresEdge) || nEdges[sampleCounter] == args.nEdges)
					{
						args.estimate += currentCliqueTree.weight;
					}
					else if (currentVertex == args.nVertices)
					{
						//Reached the end without the right number of edges (because that case is handled by the previous condition)
					}
					else if (copiedConditions[currentEdge])
					{
						//The case where this edge being off is not allowed.
						//Note that in this case we add an edge without calling unionMinimalSeparator, so we need to be careful later on
						possibilityEdges[sampleCounter] = nEdges[sampleCounter];
						
						//Add correct index to possibilities vector
						childNodes.push_back(childNodeType(sampleCounter, true, currentCliqueTree.weight, currentCliqueTree.auxWeight));
					}
					//This edge could be either present or absent, without further information
					else
					{
						currentCliqueTree.tree.unionMinimalSeparators(currentVertex, currentEdge, unionMinimalSeparators[sampleCounter], vertexSequence[sampleCounter], edgeSequence[sampleCounter], addEdges[sampleCounter], removeEdges[sampleCounter], temp);
						//If we need to add edges that were already considered (and therefore, must have already been rejected), then this edge CANNOT be present
						if ((unionMinimalSeparators[sampleCounter] & (~copiedConditions) & bitsetType((1ULL << currentEdge) - 1)).any())
						{
							//If we don't need this edge to make up the numbers, then mark it as off and 
							//continue to the next edge
							if (!requiresEdge)
							{
								childNodes.push_back(childNodeType(sampleCounter, false, currentCliqueTree.weight, currentCliqueTree.auxWeight));
							}
							//If we do, then it's impossible to reach the target and there are no children. 
						}
						else
						{
							//Here we *may* have an actual branch, both taking and not taking the edge are possible. 
							//We need to remove the ones that we've already conditioned on
							//And also remove the one edge that's already present
							bitsetType additionalEdges = unionMinimalSeparators[sampleCounter] & (~copiedConditions) & ~bitsetType((1ULL << currentEdge) - 1);
							int nAdditionalEdges = (int)additionalEdges.count();
							//If this would push us over the edge limit, then really this
							//edge can only be missing
							if (nEdges[sampleCounter] + nAdditionalEdges + 1 > args.nEdges)
							{
								childNodes.push_back(childNodeType(sampleCounter, false, currentCliqueTree.weight, currentCliqueTree.auxWeight));
							}
							else if (requiresEdge)
							{
								//If we need this edge to make up the numbers, don't consider the case where it's missing. 
								possibilityEdges[sampleCounter] = nEdges[sampleCounter] + 1 + nAdditionalEdges;
								childNodes.push_back(childNodeType(sampleCounter, true, currentCliqueTree.weight, currentCliqueTree.auxWeight * pow(auxWeightPower, nAdditionalEdges)));
							}
							else
							{
								//Here both present and absent are allowed.
								//Increase the number of total edges in the case that the edge is present
								possibilityEdges[sampleCounter] = nEdges[sampleCounter] + 1 + nAdditionalEdges;

								//Add correct indices to possibilities vector
								childNodes.push_back(childNodeType(sampleCounter, false, currentCliqueTree.weight, currentCliqueTree.auxWeight));
								childNodes.push_back(childNodeType(sampleCounter, true, currentCliqueTree.weight, currentCliqueTree.auxWeight * pow(auxWeightPower, nAdditionalEdges)));
							}
						}
					}
				}
				if(currentEdge == currentVertex - 1)
				{
					//Work out how many different graphs we have, up to isomorphism
					//To start with, get out cannonical representations
					for (int childCounter = 0; childCounter < (int)childNodes.size(); childCounter++)
					{
						childNodeType& currentChild = childNodes[childCounter];
						if(!currentChild.includesEdge())
						{
							cliqueTrees[currentChild.getParentIndex()].tree.convertToNauty(lab, ptn, orbits, nautyGraph, &(cannonicalNautyGraphs(childCounter, 0)));
						}
						else
						{
							cliqueTrees[currentChild.getParentIndex()].tree.convertToNautyWithEdge(lab, ptn, orbits, nautyGraph, &(cannonicalNautyGraphs(childCounter, 0)), currentEdge, currentVertex);
						}
					}
					//Work out which graphs are isomorphic to a graph earlier on in the set of samples. The weight for those graphs are added to the earlier one.
					std::fill(alreadyConsidered.begin(), alreadyConsidered.end(), false);
					for (int childCounter = 0; childCounter < (int)childNodes.size(); childCounter++)
					{
						if(!alreadyConsidered[childCounter])
						{
							mpfr_class weightOther = 0, weightOtherAux = 0;
							for (int childCounter2 = childCounter+1; childCounter2 < (int)childNodes.size(); childCounter2++)
							{
								if(!alreadyConsidered[childCounter2])
								{
									int m = SETWORDSNEEDED(currentVertex);
									int memcmpResult = memcmp(&(cannonicalNautyGraphs(childCounter, 0)), &(cannonicalNautyGraphs(childCounter2, 0)), m*currentVertex*sizeof(graph));
									if(memcmpResult == 0)
									{
										weightOther += childNodes[childCounter2].weight;
										weightOtherAux += childNodes[childCounter2].auxWeight;
										alreadyConsidered[childCounter2] = true;
									}
								}
							}
							childNodes[childCounter].weight += weightOther;
							childNodes[childCounter].auxWeight += weightOtherAux;
						}
					}
					std::vector<bool>::reverse_iterator alreadyConsideredIterator = std::vector<bool>::reverse_iterator(alreadyConsidered.begin() + childNodes.size());
					std::vector<childNodeType>::reverse_iterator toErase = childNodes.rbegin();
					for(std::vector<childNodeType>::reverse_iterator i = childNodes.rbegin(); i != childNodes.rend(); i++, alreadyConsideredIterator++)
					{
						if(*alreadyConsideredIterator)
						{
							std::swap(*i, *toErase);
							toErase++;
						}
					}
					childNodes.erase(toErase.base(), childNodes.end()); 
					std::sort(childNodes.begin(), childNodes.end());
				}

				if (nRemainingEdges == 1)
				{
					for(std::vector<childNodeType>::iterator i = childNodes.begin(); i != childNodes.end(); i++)
					{
						args.estimate += i->weight;
					}
					return;
				}

				int toTake = std::min((int)childNodes.size(), args.budget);
				if(toTake != (int)childNodes.size() || !args.exact)
				{
					args.exact = false;
					args.minimumSizeForExact = -1;
				}
				else args.minimumSizeForExact = std::max(args.minimumSizeForExact, toTake);
				newCliqueTrees.clear();
				//If we're taking an exhaustive sample, then skip the resampling-without-replacement section. 
				if(toTake == (int)childNodes.size())
				{
					rescaledWeights.resize(toTake);
					std::fill(rescaledWeights.begin(), rescaledWeights.end(), 1);
					indices.resize(toTake);
					for(int i = 0; i < toTake; i++) indices[i] = i;
				}
				else
				{
					weights.clear();
					for(int i = 0; i < (int)childNodes.size(); i++) weights.push_back(childNodes[i].auxWeight);

					samplingArgs.n = toTake;
					sampling::sampfordFromParetoNaive(samplingArgs, args.randomSource);
				}
				std::sort(indices.begin(), indices.end());
				//Now actually start making copies
				int i = 0;
				while (i < toTake)
				{
					childNodeType& currentChildNode = childNodes[indices[i]];
					int parentIndex = currentChildNode.getParentIndex();
					int copyCount = 1;
					if(i + 1 < toTake && childNodes[indices[i + 1]].getParentIndex() == parentIndex)
					{
						copyCount = 2;
					}
					int newIndex = (int)newCliqueTrees.size();
					newConditions[newIndex] = conditions[parentIndex];
					if (copyCount == 1)
					{
						newCliqueTrees.push_back(std::move(cliqueTrees[parentIndex]));
						newCliqueTrees.back().weight = currentChildNode.weight;
						newCliqueTrees.back().auxWeight = currentChildNode.auxWeight;
						if(childNodes[indices[i]].includesEdge())
						{
							newConditions[newIndex][currentEdge] = 1;
							newNEdges[newIndex] = possibilityEdges[parentIndex];
							if(!conditions[parentIndex][currentEdge])
							{
								newConditions[newIndex] |= unionMinimalSeparators[parentIndex];
								newCliqueTrees[newIndex].tree.addEdge(currentVertex, currentEdge, unionMinimalSeparators[parentIndex], vertexSequence[parentIndex], edgeSequence[parentIndex], addEdges[parentIndex], removeEdges[parentIndex], temp, true);
							}
						}
						else
						{
							newNEdges[newIndex] = nEdges[parentIndex];
							newConditions[newIndex] = conditions[parentIndex];
						}
						newCliqueTrees[newIndex].weight /= rescaledWeights[indices[i]];
						newCliqueTrees[newIndex].auxWeight /= rescaledWeights[indices[i]];
						i++;
					}
					else if (copyCount == 2)
					{
						newCliqueTrees.push_back(cliqueTrees[parentIndex]);
						newCliqueTrees.back().weight = currentChildNode.weight;
						newCliqueTrees.back().auxWeight = currentChildNode.auxWeight;
						newCliqueTrees.push_back(std::move(cliqueTrees[parentIndex]));
						newCliqueTrees.back().weight = childNodes[indices[i + 1]].weight;
						newCliqueTrees.back().auxWeight = childNodes[indices[i + 1]].auxWeight;
						
						newConditions[newIndex+1] = conditions[parentIndex];
						newConditions[newIndex+1][currentEdge] = 1;
						newConditions[newIndex+1] |= unionMinimalSeparators[parentIndex];
						newNEdges[newIndex+1] = possibilityEdges[parentIndex];
						newCliqueTrees[newIndex+1].tree.addEdge(currentVertex, currentEdge, unionMinimalSeparators[parentIndex], vertexSequence[parentIndex], edgeSequence[parentIndex], addEdges[parentIndex], removeEdges[parentIndex], temp, true);
						
						newNEdges[newIndex] = nEdges[parentIndex];
						newConditions[newIndex] = conditions[parentIndex];
						
						newCliqueTrees[newIndex].weight /= rescaledWeights[indices[i]];
						newCliqueTrees[newIndex+1].weight /= rescaledWeights[indices[i+1]];
						newCliqueTrees[newIndex].auxWeight /= rescaledWeights[indices[i]];
						newCliqueTrees[newIndex+1].auxWeight /= rescaledWeights[indices[i+1]];
						i += 2;
					}
				}
				newCliqueTrees.swap(cliqueTrees);
				nEdges.swap(newNEdges);
				conditions.swap(newConditions);
				currentEdge++;
				//If we've reached the end prematurely, break. 
				if (cliqueTrees.size() == 0) return;
			}
		}
	}
	template void horvitzThompson2<cliqueTree>(horvitzThompsonArgs& args);
	template void horvitzThompson2<cliqueTreeAdjacencyMatrix>(horvitzThompsonArgs& args);
}
