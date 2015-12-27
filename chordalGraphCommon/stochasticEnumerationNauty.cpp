#include "stochasticEnumerationNauty.h"
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/math/special_functions.hpp>
#include "nauty.h"
#include "cliqueTree.h"
#include "cliqueTreeAdjacencyMatrix.h"
namespace chordalGraph
{
	template<typename cliqueTree> void stochasticEnumerationNauty(stochasticEnumerationNautyArgs& args)
	{
		args.exact = true;
		args.minimumSizeForExact = -1;
		args.estimate = 0;
		mpfr_class multiple = 1;

		//Temporary data that's used in cliqueTree calls
		typename cliqueTree::unionMinimalSeparatorsTemporaries temp;

		boost::random_number_generator<boost::mt19937> generator(args.randomSource);
		//Number of edges either present (or to be added later)
		std::vector<int> nEdges(args.budget);
		std::vector<stochasticEnumerationNautyPrivate::weightedCliqueTree<cliqueTree> > cliqueTrees;
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
			stochasticEnumerationNautyPrivate::weightedCliqueTree<cliqueTree> initialTree(args.nVertices);
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

		//This tells us how many children of a sample were taken.
		//This is important, because if it's only one child then we don't need to make a copy of the clique tree, 
		//we can use the move constructor instead. 
		std::vector<int> copyCounts(args.budget);


		std::vector<stochasticEnumerationNautyPrivate::weightedCliqueTree<cliqueTree> > newCliqueTrees;
		newCliqueTrees.reserve(args.budget);
		std::vector<int> newNEdges(args.budget);

		std::vector<std::list<typename cliqueTree::cliqueTreeGraphType::vertex_descriptor> > vertexSequence(args.budget);
		std::vector<std::list<typename cliqueTree::externalEdge> > edgeSequence(args.budget);
		std::vector<std::vector<typename cliqueTree::externalEdge> > removeEdges(args.budget);
		std::vector<std::vector<typename cliqueTree::externalEdge> > addEdges(args.budget);

		std::vector<bitsetType> unionMinimalSeparators(args.budget);
		//Vector used to shuffle indices
		std::vector<int> shuffleVector;
		shuffleVector.reserve(args.budget);

		//Nauty variables
		std::vector<int> lab, ptn, orbits;
		lab.reserve(args.nVertices);
		ptn.reserve(args.nVertices);
		std::vector<graph> nautyGraph;
		std::vector<std::vector<graph> > cannonicalNautyGraphs(args.budget);

		//Used to count the number of distinct graphs (up to isomorphism)
		std::vector<bool> alreadyConsidered(args.budget);
		//Continue while there are samples left.
		while (currentVertex < args.nVertices - 1)
		{
			currentVertex++;
			currentEdge = 0;
			//Work out how many different graphs we have, up to isomorphism
			//To start with, get out cannonical representations
			for (int sampleCounter = 0; sampleCounter < (int)cliqueTrees.size(); sampleCounter++)
			{
				cliqueTrees[sampleCounter].tree.convertToNauty(lab, ptn, orbits, nautyGraph, cannonicalNautyGraphs[sampleCounter]);
			}
			//Work out which graphs are isomorphic to a graph earlier on in the set of samples. The weight for those graphs are added to the earlier one. 
			std::fill(alreadyConsidered.begin(), alreadyConsidered.end(), false);
			for(int sampleCounter = 0; sampleCounter < (int)cliqueTrees.size(); sampleCounter++)
			{
				if(!alreadyConsidered[sampleCounter])
				{
					int weightOther = 0;
					for(int sampleCounter2 = sampleCounter+1; sampleCounter2 < (int)cliqueTrees.size(); sampleCounter2++)
					{
						if(!alreadyConsidered[sampleCounter2])
						{
							int m = SETWORDSNEEDED(currentVertex);
							int memcmpResult = memcmp(&(cannonicalNautyGraphs[sampleCounter][0]), &(cannonicalNautyGraphs[sampleCounter2][0]), m*currentVertex*sizeof(graph));
							if(memcmpResult == 0)
							{
								weightOther+= cliqueTrees[sampleCounter2].weight;
								alreadyConsidered[sampleCounter2] = true;
							}
						}
					}
					cliqueTrees[sampleCounter].weight += weightOther;
				}
			}
			//There are no conditions when we move to the next vertex
			std::fill(conditions.begin(), conditions.end(), 0);
			//Add the extra vertex
			std::for_each(cliqueTrees.begin(), cliqueTrees.end(), std::mem_fun_ref(&stochasticEnumerationNautyPrivate::weightedCliqueTree<cliqueTree>::addVertex));

			while(currentEdge < currentVertex)
			{
				//Clear vector of indices of possibilities
				shuffleVector.clear();
				//count the number of nodes which have cost 1
				int knownToBeChordalWeight = 0;
				int nRemainingEdges = currentVertex - currentEdge + ((args.nVertices - currentVertex - 1)* (args.nVertices - 2 - currentVertex) / 2) + (args.nVertices - currentVertex - 1) * (currentVertex + 1);
				//Clear data structures;
				std::for_each(vertexSequence.begin(), vertexSequence.end(), std::mem_fun_ref(&std::list<typename cliqueTree::cliqueTreeGraphType::vertex_descriptor>::clear));
				std::for_each(edgeSequence.begin(), edgeSequence.end(), std::mem_fun_ref(&std::list<typename cliqueTree::externalEdge>::clear));
				std::for_each(removeEdges.begin(), removeEdges.end(), std::mem_fun_ref(&std::vector<typename cliqueTree::externalEdge>::clear));
				std::for_each(addEdges.begin(), addEdges.end(), std::mem_fun_ref(&std::vector<typename cliqueTree::externalEdge>::clear));
				std::fill(unionMinimalSeparators.begin(), unionMinimalSeparators.end(), 0);
				for (int sampleCounter = 0; sampleCounter < (int)cliqueTrees.size(); sampleCounter++)
				{
					if(alreadyConsidered[sampleCounter] && currentEdge == 0) continue;
					//Remaining edges, including the one we're just about to consider. 
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
						knownToBeChordalWeight += cliqueTrees[sampleCounter].weight;
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
						shuffleVector.push_back(2 * sampleCounter + 1);
					}
					//This edge could be either present or absent, without further information
					else
					{
						cliqueTrees[sampleCounter].tree.unionMinimalSeparators(currentVertex, currentEdge, unionMinimalSeparators[sampleCounter], vertexSequence[sampleCounter], edgeSequence[sampleCounter], addEdges[sampleCounter], removeEdges[sampleCounter], temp);
						//If we need to add edges that were already considered (and therefore, must have already been rejected), then this edge CANNOT be present
						if ((unionMinimalSeparators[sampleCounter] & (~copiedConditions) & bitsetType((1ULL << currentEdge) - 1)).any())
						{
							//If we don't need this edge to make up the numbers, then mark it as off and 
							//continue to the next edge
							if (!requiresEdge)
							{
								shuffleVector.push_back(2 * sampleCounter);
							}
							//If we do, then it's impossible to reach the target and there are no children. 
						}
						else
						{
							//We need to remove the ones that we've already conditioned on
							//And also remove the one edge that's already present
							bitsetType additionalEdges = unionMinimalSeparators[sampleCounter] & (~copiedConditions) & ~bitsetType((1ULL << currentEdge) - 1);
							int nAdditionalEdges = (int)additionalEdges.count();
							//If this would push us over the edge limit, then really this
							//edge can only be missing
							if (nEdges[sampleCounter] + nAdditionalEdges + 1 > args.nEdges)
							{
								//Add correct sampleCounter to possibilities vector
								shuffleVector.push_back(2 * sampleCounter);
							}
							else if (requiresEdge)
							{
								//If we need this edge to make up the numbers, don't consider the case where it's missing.
								possibilityEdges[sampleCounter] = nEdges[sampleCounter] + 1 + nAdditionalEdges;
								shuffleVector.push_back(2 * sampleCounter + 1);
							}
							else
							{
								//Here both present and absent are allowed.
								//Increase the number of total edges in the case that the edge is present
								possibilityEdges[sampleCounter] = nEdges[sampleCounter] + 1 + nAdditionalEdges;

								//Add correct sampleCounter to possibilities vector
								shuffleVector.push_back(2 * sampleCounter);
								shuffleVector.push_back(2 * sampleCounter + 1);
							}
						}
					}
				}
				args.estimate += multiple * knownToBeChordalWeight;

				if (nRemainingEdges == 1)
				{
					knownToBeChordalWeight = 0;
					for(std::vector<int>::iterator i = shuffleVector.begin(); i != shuffleVector.end(); i++)
					{
						knownToBeChordalWeight += cliqueTrees[(*i) / 2].weight;
					}
					args.estimate += multiple * knownToBeChordalWeight;
					return;
				}

				boost::range::random_shuffle(shuffleVector, generator);
				int toTake = std::min((int)shuffleVector.size(), args.budget);

				if(toTake != (int)shuffleVector.size() || !args.exact)
				{
					args.exact = false;
					args.minimumSizeForExact = -1;
				}
				else args.minimumSizeForExact = std::max(args.minimumSizeForExact, toTake);
				//Ratio of vertices examined to not examined
				multiple *= (double)shuffleVector.size() / (double)toTake;

				std::fill(copyCounts.begin(), copyCounts.end(), 0);
				//Now work out how many copies (0, 1, 2) are taken of each sample
				for (int i = 0; i < toTake; i++)
				{
					copyCounts[shuffleVector[i]/2]++;
				}

				//Now actually start making copies
				newCliqueTrees.clear();
				for (int i = 0; i < toTake; i++)
				{
					int originalIndex = shuffleVector[i] / 2;
					newConditions[i] = conditions[originalIndex];

					if (copyCounts[originalIndex] == 1)
					{
						newCliqueTrees.push_back(std::move(cliqueTrees[originalIndex]));
					}
					else if (copyCounts[originalIndex] == 2)
					{
						copyCounts[originalIndex]--;
						newCliqueTrees.push_back(cliqueTrees[originalIndex]);
					}
					//This should never happen. 
					else
					{
						throw std::runtime_error("Internal error");
					}
					//Did we take the case where we added an edge?
					if (shuffleVector[i] % 2 == 1)
					{
						newConditions[i][currentEdge] = 1;
						newNEdges[i] = possibilityEdges[originalIndex];
						if (!conditions[originalIndex][currentEdge])
						{
							newConditions[i] |= unionMinimalSeparators[originalIndex];
							newCliqueTrees[i].tree.addEdge(currentVertex, currentEdge, unionMinimalSeparators[originalIndex], vertexSequence[originalIndex], edgeSequence[originalIndex], addEdges[originalIndex], removeEdges[originalIndex], temp, true);
						}
					}
					else
					{
						newNEdges[i] = nEdges[originalIndex];
						newConditions[i] = conditions[originalIndex];
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
	template void stochasticEnumerationNauty<cliqueTree>(stochasticEnumerationNautyArgs& args);
	template void stochasticEnumerationNauty<cliqueTreeAdjacencyMatrix>(stochasticEnumerationNautyArgs& args);
}
