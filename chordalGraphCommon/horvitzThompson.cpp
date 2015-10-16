#include "horvitzThompson.h"
#include "cliqueTree.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions.hpp>
#include "nauty.h"
#include "sampford.h"
namespace chordalGraph
{
	namespace horvitzThompsonPrivate
	{
		struct weightedCliqueTree
		{
		public:
			weightedCliqueTree(weightedCliqueTree&& other)
				: tree(std::move(other.tree)), weight(other.weight)
			{} 
			weightedCliqueTree(const weightedCliqueTree& other)
				: tree(other.tree), weight(other.weight)
			{}
			weightedCliqueTree(int nVertices)
				: tree(nVertices), weight(1)
			{}
			cliqueTree tree;
			double weight;
		};
	}
	using horvitzThompsonPrivate::weightedCliqueTree;
	void horvitzThompson(horvitzThompsonArgs& args)
	{
		args.exact = true;
		args.estimate = 0;

		//Temporary data that's used in cliqueTree calls
		cliqueTree::unionMinimalSeparatorsTemporaries temp;

		boost::random_number_generator<boost::mt19937> generator(args.randomSource);
		boost::random::bernoulli_distribution<> standardBernoulli(0.5);
		//Number of edges either present (or to be added later)
		std::vector<int> nEdges(args.budget);
		std::vector<weightedCliqueTree> cliqueTrees;
		cliqueTrees.reserve(args.budget);

		//At any point we will have a bunch of conditions on what other edges are required to be present (in order to maintain chordality)
		//and what edges have already been marked as present
		//This vector tracks those conditions. Entry 0 contains the conditions for the first sample, then the following entry contains
		//the conditions for the next sample, etc
		std::vector<bitsetType> conditions(args.budget);

		//The vertex we're currently considering for a certain sample
		std::vector<int> currentVertex;
		//The edge we're currently considering
		std::vector<int> currentEdge;
		//We start off with one sample
		{
			weightedCliqueTree initialTree(args.nVertices);
			initialTree.tree.addVertex();
			initialTree.tree.addVertex();
			cliqueTrees.push_back(initialTree);
			nEdges[0] = 0;

			//It has zero vertices and zero edges
			currentEdge.push_back(0);
			currentVertex.push_back(1);
		}
		//It has no conditions
		conditions[0] = 0;

		//The number of total edges, if the new edge is added. We only need one entry
		//per sample, not two, because if the new edge is not added, then the number of edges
		//stays the same (and this is recorded in nEdges)
		std::vector<int> possibilityEdges(args.budget);

		//The conditions, for the new set of samples
		std::vector<bitsetType> newConditions(args.budget);

		std::vector<int> newCurrentVertex, newCurrentEdge;

		//This tells us how many children of a sample were taken.
		//This is important, because if it's only one child then we don't need to make a copy of the clique tree, 
		//we can use the move constructor instead. 
		std::vector<int> copyCounts(args.budget);

		std::vector<weightedCliqueTree> newCliqueTrees;
		newCliqueTrees.reserve(args.budget);
		std::vector<int> newNEdges(args.budget);

		std::vector<std::list<cliqueTree::cliqueTreeGraphType::vertex_descriptor> > vertexSequence(args.budget);
		std::vector<std::list<cliqueTree::externalEdge> > edgeSequence(args.budget);
		std::vector<std::vector<cliqueTree::externalEdge> > removeEdges(args.budget);
		std::vector<std::vector<cliqueTree::externalEdge> > addEdges(args.budget);

		std::vector<bitsetType> unionMinimalSeparators(args.budget);
		//Vector used to shuffle indices
		std::vector<int> hasChildren;
		hasChildren.reserve(args.budget);
		sampfordBruteForceArgs sampfordArgs(args.randomSource);

		//Nauty variables
		std::vector<int> lab, ptn, orbits;
		lab.reserve(args.nVertices);
		ptn.reserve(args.nVertices);
		std::vector<graph> nautyGraph;
		std::vector<std::vector<graph> > cannonicalNautyGraphs(args.budget);

		//Used to count the number of distinct graphs (up to isomorphism)
		std::vector<bool> alreadyConsidered(args.budget);
		//Continue while there are samples left.
		while (currentVertex.size() > 0)
		{
			//Work out how many different graphs we have, up to isomorphism
			//To start with, get out cannonical representations
			for (int sampleCounter = 0; sampleCounter < (int)currentVertex.size(); sampleCounter++)
			{
				if(currentEdge[sampleCounter] == 0)
				{
					cliqueTrees[sampleCounter].tree.convertToNauty(lab, ptn, orbits, nautyGraph, cannonicalNautyGraphs[sampleCounter]);
				}
			}
			//Work out which graphs are isomorphic to a graph earlier on in the set of samples. The weight for those graphs are added to the earlier one. 
			std::fill(alreadyConsidered.begin(), alreadyConsidered.end(), false);
			for(int sampleCounter = 0; sampleCounter < (int)currentVertex.size(); sampleCounter++)
			{
				if(!alreadyConsidered[sampleCounter] && currentEdge[sampleCounter] == 0)
				{
					double weightOther = 0;
					for(int sampleCounter2 = sampleCounter+1; sampleCounter2 < (int)currentVertex.size(); sampleCounter2++)
					{
						if(!alreadyConsidered[sampleCounter2] && currentEdge[sampleCounter2] == 0 && currentVertex[sampleCounter] == currentVertex[sampleCounter2])
						{
							int m = SETWORDSNEEDED(currentVertex[sampleCounter]);
							int memcmpResult = memcmp(&(cannonicalNautyGraphs[sampleCounter][0]), &(cannonicalNautyGraphs[sampleCounter2][0]), m*currentVertex[sampleCounter]*sizeof(graph));
							if(memcmpResult == 0)
							{
								weightOther += cliqueTrees[sampleCounter2].weight;
								cliqueTrees[sampleCounter2].weight = 0;
								alreadyConsidered[sampleCounter2] = true;
							}
						}
					}
					cliqueTrees[sampleCounter].weight += weightOther;
				}
			}

			//Clear vector of indices of possibilities
			hasChildren.clear();
			sampfordArgs.weights.clear();
			double sumHasChildrenWeights = 0;
			for (int sampleCounter = 0; sampleCounter < (int)currentVertex.size(); sampleCounter++)
			{
				if(alreadyConsidered[sampleCounter]) continue;
				while (true)
				{
					//Clear data structures;
					vertexSequence[sampleCounter].clear();
					edgeSequence[sampleCounter].clear();
					removeEdges[sampleCounter].clear();
					addEdges[sampleCounter].clear();
					unionMinimalSeparators[sampleCounter] = 0;
					//Remaining edges, including the one we're just about to consider. 
					int& sampleCurrentVertex = currentVertex[sampleCounter];
					int& sampleCurrentEdge = currentEdge[sampleCounter];
					int nRemainingEdges = sampleCurrentVertex - sampleCurrentEdge + ((args.nVertices - sampleCurrentVertex - 1)* (args.nVertices - 2 - sampleCurrentVertex) / 2) + (args.nVertices - sampleCurrentVertex - 1) * (sampleCurrentVertex + 1);
					bitsetType copiedConditions = conditions[sampleCounter];
					//maximum possible number of edges
					int maxEdges = nEdges[sampleCounter] + nRemainingEdges - (int)(copiedConditions & ~bitsetType((1ULL << (sampleCurrentEdge+1)) - 1)).count();
					//Do we need this edge to make up the numbers?
					bool requiresEdge = maxEdges == args.nEdges;
					if (maxEdges < args.nEdges)
					{
						//Too few edges to reach the target number. No children. 
						break;
					}
					//The first condition corresponds to the case where adding ALL edges past this point gives us the target number
					//The second corresponds to the case where we've already hit the target
					else if ((sampleCurrentEdge == 0 && requiresEdge) || nEdges[sampleCounter] == args.nEdges)
					{
						args.estimate +=cliqueTrees[sampleCounter].weight;
						break;
					}
					else if (sampleCurrentVertex == args.nVertices)
					{
						//Reached the end without the right number of edges (because that case is handled by the previous condition)
						break;
					}
					else if (copiedConditions[sampleCurrentEdge])
					{
						//We already conditioned on this edge being in, so continue
					}
					//This edge could be either present or absent, without further information
					else
					{
						cliqueTrees[sampleCounter].tree.unionMinimalSeparators(sampleCurrentVertex, sampleCurrentEdge, unionMinimalSeparators[sampleCounter], vertexSequence[sampleCounter], edgeSequence[sampleCounter], addEdges[sampleCounter], removeEdges[sampleCounter], temp);
						//If we need to add edges that were already considered (and therefore, must have already been rejected), then this edge CANNOT be present
						if ((unionMinimalSeparators[sampleCounter] & (~copiedConditions) & bitsetType((1ULL << sampleCurrentEdge) - 1)).any())
						{
							//If we don't need this edge to make up the numbers, then mark it as off and 
							//continue to the next edge
							if (!requiresEdge)
							{
								//we don't have to do anything here
							}
							//If we do, then it's impossible to reach the target and there are no children. 
							else break;
						}
						else
						{
							//Here we *may* have an actual branch, both taking and not taking the edge are possible. 
							//We need to remove the ones that we've already conditioned on
							//And also remove the one edge that's already present
							bitsetType additionalEdges = unionMinimalSeparators[sampleCounter] & (~copiedConditions) & ~bitsetType((1ULL << sampleCurrentEdge) - 1);
							int nAdditionalEdges = (int)additionalEdges.count();
							//If this would push us over the edge limit, then really this
							//edge can only be missing
							if (nEdges[sampleCounter] + nAdditionalEdges + 1 > args.nEdges)
							{
							}
							else if (requiresEdge)
							{
								//If we need this edge to make up the numbers, don't consider the case where it's missing. 
								nEdges[sampleCounter] = nEdges[sampleCounter] + 1 + nAdditionalEdges;
								cliqueTrees[sampleCounter].tree.addEdge(sampleCurrentVertex, sampleCurrentEdge, unionMinimalSeparators[sampleCounter], vertexSequence[sampleCounter], edgeSequence[sampleCounter], addEdges[sampleCounter], removeEdges[sampleCounter], temp, true);
								conditions[sampleCounter][sampleCurrentEdge] = true;
								conditions[sampleCounter] |= unionMinimalSeparators[sampleCounter];
							}
							else
							{
								//Here both present and absent are allowed.
								//Increase the number of total edges in the case that the edge is present
								possibilityEdges[sampleCounter] = nEdges[sampleCounter] + 1 + nAdditionalEdges;

								//Add correct indices to possibilities vector
								hasChildren.push_back(sampleCounter);
								sampfordArgs.weights.push_back(cliqueTrees[sampleCounter].weight);
								sampfordArgs.weights.push_back(cliqueTrees[sampleCounter].weight);
								sumHasChildrenWeights += cliqueTrees[sampleCounter].weight;
								break;
							}
						}
					}
					sampleCurrentEdge++;
					if (sampleCurrentEdge == sampleCurrentVertex)
					{
						sampleCurrentVertex++;
						sampleCurrentEdge = 0;
						conditions[sampleCounter] = 0;
						cliqueTrees[sampleCounter].tree.addVertex();
					}
				}
			}
			int toTake = std::min(2*(int)hasChildren.size(), args.budget);
			if(toTake != 2*(int)hasChildren.size()) args.exact = false;
			//If we're taking an exhaustive sample, then skip the resampling-without-replacement section. 
			if(toTake == 2*(int)hasChildren.size())
			{
				std::fill(copyCounts.begin(), copyCounts.end(), 2);
				sampfordArgs.inclusionProbabilities.resize(toTake);
				std::fill(sampfordArgs.inclusionProbabilities.begin(), sampfordArgs.inclusionProbabilities.end(), 1);
			}
			else
			{
				std::fill(copyCounts.begin(), copyCounts.end(), 0);
				sampfordArgs.n = toTake;
				sampfordBruteForce(sampfordArgs);
				for(int i = 0; i < (int)toTake; i++)
				{
					copyCounts[sampfordArgs.indices[i]/2]++;
				}
			}

			//Now actually start making copies
			newCliqueTrees.clear();
			newCurrentVertex.clear();
			newCurrentEdge.clear();
			for (int i = 0; i < (int)hasChildren.size(); i++)
			{
				int originalIndex = hasChildren[i];
				int newIndex = newCliqueTrees.size();

				if (copyCounts[i] == 1)
				{
					newCliqueTrees.push_back(std::move(cliqueTrees[originalIndex]));
					if(standardBernoulli(args.randomSource))
					{
						newConditions[newIndex] = conditions[originalIndex];
						newConditions[newIndex][currentEdge[originalIndex]] = 1;
						newConditions[newIndex] |= unionMinimalSeparators[originalIndex];
						newNEdges[newIndex] = possibilityEdges[originalIndex];
						newCliqueTrees[newIndex].tree.addEdge(currentVertex[originalIndex], currentEdge[originalIndex], unionMinimalSeparators[originalIndex], vertexSequence[originalIndex], edgeSequence[originalIndex], addEdges[originalIndex], removeEdges[originalIndex], temp, true);
					}
					else
					{
						newNEdges[newIndex] = nEdges[originalIndex];
						newConditions[newIndex] = conditions[originalIndex];
					}
					newCurrentVertex.push_back(currentVertex[originalIndex]);
					newCurrentEdge.push_back(currentEdge[originalIndex]);
					newCurrentEdge[newIndex]++;
					if (newCurrentEdge[newIndex] == newCurrentVertex[newIndex])
					{
						newCurrentEdge[newIndex] = 0;
						newCurrentVertex[newIndex]++;
						newCliqueTrees[newIndex].tree.addVertex();
						newConditions[newIndex] = 0;
					}
					newCliqueTrees[newIndex].weight /= sampfordArgs.inclusionProbabilities[2*i];
				}
				else if (copyCounts[i] == 2)
				{
					newCliqueTrees.push_back(cliqueTrees[originalIndex]);
					newCliqueTrees.push_back(std::move(cliqueTrees[originalIndex]));
					
					newConditions[newIndex] = conditions[originalIndex];
					newConditions[newIndex][currentEdge[originalIndex]] = 1;
					newConditions[newIndex] |= unionMinimalSeparators[originalIndex];
					newNEdges[newIndex] = possibilityEdges[originalIndex];
					newCliqueTrees[newIndex].tree.addEdge(currentVertex[originalIndex], currentEdge[originalIndex], unionMinimalSeparators[originalIndex], vertexSequence[originalIndex], edgeSequence[originalIndex], addEdges[originalIndex], removeEdges[originalIndex], temp, true);
					
					newNEdges[newIndex+1] = nEdges[originalIndex];
					newConditions[newIndex+1] = conditions[originalIndex];
					
					newCurrentVertex.push_back(currentVertex[originalIndex]);
					newCurrentVertex.push_back(currentVertex[originalIndex]);
					
					newCurrentEdge.push_back(currentEdge[originalIndex]);
					newCurrentEdge.push_back(currentEdge[originalIndex]);

					newCurrentEdge[newIndex]++;
					newCurrentEdge[newIndex+1]++;

					if (newCurrentEdge[newIndex] == newCurrentVertex[newIndex])
					{
						newCurrentEdge[newIndex] = 0; newCurrentEdge[newIndex+1] = 0;
						newCurrentVertex[newIndex]++; newCurrentVertex[newIndex+1]++;
						newCliqueTrees[newIndex].tree.addVertex(); newCliqueTrees[newIndex+1].tree.addVertex();
						newConditions[newIndex] = 0; newConditions[newIndex+1] = 0;
					}
					newCliqueTrees[newIndex].weight /= sampfordArgs.inclusionProbabilities[2*i];
					newCliqueTrees[newIndex+1].weight = newCliqueTrees[newIndex].weight;
				}
				else if(copyCounts[i] == 0)
				{
				}
				//This should never happen. 
				else
				{
					throw std::runtime_error("Internal error");
				}
			}
			newCliqueTrees.swap(cliqueTrees);
			nEdges.swap(newNEdges);
			conditions.swap(newConditions);
			currentVertex.swap(newCurrentVertex);
			currentEdge.swap(newCurrentEdge);
			//If we've reached the end prematurely, break. 
			if (cliqueTrees.size() == 0) return;
		}
	}
}
