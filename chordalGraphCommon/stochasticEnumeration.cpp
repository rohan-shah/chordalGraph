#include "stochasticEnumeration.h"
#include "cliqueTree.h"
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/math/special_functions.hpp>
namespace chordalGraph
{
	void stochasticEnumeration(stochasticEnumerationArgs& args)
	{
		args.estimate = 0;
		mpfr_class multiple = 1;

		//Temporary data that's used in cliqueTree calls
		cliqueTree::unionMinimalSeparatorsTemporaries temp;

		boost::random_number_generator<boost::mt19937> generator(args.randomSource);
		//Number of edges either present (or to be added later)
		std::vector<int> nEdges(args.budget);
		std::vector<cliqueTree> cliqueTrees;
		cliqueTrees.reserve(args.budget);

		//At any point we will have a bunch of conditions on what other edges are required to be present (in order to maintain chordality)
		//This vector tracks those conditions. Entry 0 contains the conditions for the first sample, then the following entry contains
		//the conditions for the next sample, etc
		std::vector<bitsetType> conditions(args.budget);

		//At the start of the algorithm vertex 0 has already been added, so we are
		//looking at edges from vertex 1
		int currentVertex = 0;
		//The first edge considered is the one with vertex 0
		int currentEdge = 0;
		//We start off with one sample
		{
			cliqueTree initialTree(args.nVertices);
			initialTree.addVertex();
			cliqueTrees.push_back(initialTree);
			nEdges[0] = 0;
		}
		//It has no conditions
		conditions[0] = 0;

		//The number of total edges, if the new edge is added. We only need one entry
		//per sample, not two, because if the new edge is not added, then the number of edges
		//stays the same (and this is recorded in nEdges)
		std::vector<int> possibilityEdges(args.budget);

		//The conditions, if a new edge is added. We only need one entry per sample because if a new edge is not added,
		//the conditions will stay the same. 
		std::vector<bitsetType> newConditions(args.budget);
		//Existing edges for the current vertex
		std::vector<bitsetType> existingEdges(args.budget);
		std::vector<bitsetType> newExistingEdges(args.budget);

		//This tells us how many children of a sample were taken.
		//This is important, because if it's only one child then we don't need to make a copy of the clique tree, 
		//we can use the move constructor instead. 
		std::vector<int> copyCounts(args.budget);


		std::vector<cliqueTree> newCliqueTrees;
		newCliqueTrees.reserve(args.budget);
		std::vector<int> newNEdges(args.budget);

		std::vector<std::list<cliqueTree::cliqueTreeGraphType::vertex_descriptor> > vertexSequence(args.budget);
		std::vector<std::list<cliqueTree::externalEdge> > edgeSequence(args.budget);
		std::vector<std::vector<cliqueTree::externalEdge> > removeEdges(args.budget);
		std::vector<std::vector<cliqueTree::externalEdge> > addEdges(args.budget);

		std::vector<bitsetType> unionMinimalSeparators(args.budget);
		//Vector used to shuffle indices
		std::vector<int> shuffleVector;
		shuffleVector.reserve(args.budget);
		//Continue until we reach the vertex limit.
		while (currentVertex < args.nVertices-1)
		{
			currentVertex++;
			currentEdge = 0;
			//There are no conditions when we move to the next vertex
			std::fill(conditions.begin(), conditions.end(), 0);
			std::fill(existingEdges.begin(), existingEdges.end(), 0);
			//Add the extra vertex
			std::for_each(cliqueTrees.begin(), cliqueTrees.end(), std::mem_fun_ref(&cliqueTree::addVertex));
			while (currentEdge < currentVertex)
			{
				//Remaining edges, including the one we're just about to consider. 
				int nRemainingEdges = currentVertex - currentEdge + ((args.nVertices - currentVertex- 1)* (args.nVertices-2 - currentVertex) / 2) + (args.nVertices - currentVertex-1) * (currentVertex+1);

				//Clear data structures;
				std::for_each(vertexSequence.begin(), vertexSequence.end(), std::mem_fun_ref(&std::list<cliqueTree::cliqueTreeGraphType::vertex_descriptor>::clear));
				std::for_each(edgeSequence.begin(), edgeSequence.end(), std::mem_fun_ref(&std::list<cliqueTree::externalEdge>::clear));
				std::for_each(removeEdges.begin(), removeEdges.end(), std::mem_fun_ref(&std::vector<cliqueTree::externalEdge>::clear));
				std::for_each(addEdges.begin(), addEdges.end(), std::mem_fun_ref(&std::vector<cliqueTree::externalEdge>::clear));

				std::fill(unionMinimalSeparators.begin(), unionMinimalSeparators.end(), 0);
				//Clear vector of indices of possibilities
				shuffleVector.clear();
				//Make a list of all the possibiltiies
				int knownToBeChordal = 0;
				for (std::vector<cliqueTree>::iterator i = cliqueTrees.begin(); i != cliqueTrees.end(); i++)
				{
					int index = (int)std::distance(cliqueTrees.begin(), i);
					bitsetType copiedConditions = conditions[index];
					//maximum possible number of edges
					int maxEdges = nEdges[index] + nRemainingEdges - (int)(copiedConditions & ~bitsetType((1ULL << (currentEdge)) - 1)).count();
					//Do we need this edge to make up the numbers?
					bool requiresEdge = maxEdges == args.nEdges;
					if (maxEdges < args.nEdges)
					{
						//Too few edges to reach the target number
					}
					else if ((currentEdge == 0 && requiresEdge) || nEdges[index] == args.nEdges)
					{
						//Already known to be chordal
						knownToBeChordal++;
					}
					else if (copiedConditions[currentEdge])
					{
						//The case where this edge being off is not allowed.
						//Note that in this case we add an edge without calling unionMinimalSeparator, so we need to be careful later on
						possibilityEdges[index] = nEdges[index];

						//Add correct index to possibilities vector
						shuffleVector.push_back(2 * index + 1);
					}
					//This edge could be either present or absent, without further information
					else
					{
						cliqueTrees[index].unionMinimalSeparators(currentVertex, currentEdge, unionMinimalSeparators[index], vertexSequence[index], edgeSequence[index], addEdges[index], removeEdges[index], temp);
						//If we need to add edges that were already considered (and therefore, must have already been rejected), then this edge CANNOT be present
						//Of course *one* vertex of the minimal separator must correspond to an already added edge, and that's ok. 
						if ((unionMinimalSeparators[index] & ~existingEdges[index] & bitsetType((1ULL << currentEdge) - 1)).any())
						{
							if (!requiresEdge)
							{
								//Add correct index to possibilities vector
								shuffleVector.push_back(2 * index);
							}
						}
						else
						{
							//We need to remove the ones that we've already conditioned on
							//And also remove the one edge that's already present
							bitsetType additionalEdges = unionMinimalSeparators[index] & (~copiedConditions) & ~bitsetType((1ULL << currentEdge) - 1);
							int nAdditionalEdges = (int)additionalEdges.count();
							//If this would push us over the edge limit, then really this
							//edge can only be missing
							if (nEdges[index] + nAdditionalEdges + 1 > args.nEdges)
							{
								//Add correct index to possibilities vector
								shuffleVector.push_back(2 * index);
							}
							else if (requiresEdge)
							{
								//If we need this edge to make up the numbers, don't consider the case where it's missing. 
								possibilityEdges[index] = nEdges[index] + 1 + nAdditionalEdges;
								shuffleVector.push_back(2 * index + 1);
							}
							else
							{
								//Here both present and absent are allowed.
								//Increase the number of total edges in the case that the edge is present
								possibilityEdges[index] = nEdges[index] + 1 + nAdditionalEdges;

								//Add correct index to possibilities vector
								shuffleVector.push_back(2 * index);
								shuffleVector.push_back(2 * index + 1);
							}
						}
					}
				}
				args.estimate += multiple * knownToBeChordal;
				boost::range::random_shuffle(shuffleVector, generator);
				int toTake = std::min((int)shuffleVector.size(), args.budget);

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
					newExistingEdges[i] = existingEdges[originalIndex];
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
						newExistingEdges[i][currentEdge] = 1;
						newNEdges[i] = possibilityEdges[originalIndex];
						//In these cases we are adding edges without calling unionMinimalSeparator / addEdge.
						//In the first case this is because we already conditioned on this edge being present
						//(therefore it's covered by a previous call to addEdge)
						if (conditions[originalIndex][currentEdge])
						{
							newConditions[i] = conditions[originalIndex];
						}
						else
						{
							newCliqueTrees[i].addEdge(currentVertex, currentEdge, unionMinimalSeparators[originalIndex], vertexSequence[originalIndex], edgeSequence[originalIndex], addEdges[originalIndex], removeEdges[originalIndex], temp, true);
							newConditions[i] = conditions[originalIndex] | (unionMinimalSeparators[originalIndex] & (~bitsetType((1ULL << (currentEdge+1)) - 1)) & (~newExistingEdges[i]));
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
				existingEdges.swap(newExistingEdges);
				currentEdge++;
				//If we've reached the end prematurely, break. 
				if (cliqueTrees.size() == 0) return;
				//If nRemainingEdges == 1 then everything in shuffleVector represents a valid chordal graph. So
				//Count these final ones and then return. 
				if (nRemainingEdges == 1)
				{
					args.estimate += multiple * shuffleVector.size();
					return;
				}
			}
		}
	}
}