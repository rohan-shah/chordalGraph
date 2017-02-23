#include "posteriorInference.h"
#include <boost/graph/depth_first_search.hpp>
#include "cliqueTreeAdjacencyMatrix.h"
#include "customMCMCSymmetric.h"
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/numeric/ublas/lu.hpp>
namespace chordalGraph
{
	typedef cliqueTreeAdjacencyMatrix cliqueTreeType;
	mpfr_class multivariateGammaFunction(int m, double alpha)
	{
		mpfr_class result = boost::multiprecision::pow(boost::math::constants::pi<mpfr_class>(), m* (m-1.0)/4.0);
		for(int i = 1; i < m+1; i++)
		{
			result *= boost::math::tgamma<mpfr_class>(alpha - (i - 1.0)/2.0);
		}
		return result;
	}
	void extractSubmatrix(boost::numeric::ublas::matrix<double>& psi, boost::numeric::ublas::matrix<double>& submatrix, bitsetType contents, int dimension)
	{
		std::vector<int> indices;
		for(int i = 0; i < dimension; i++)
		{
			if(contents[i]) indices.push_back(i);
		}
		for(std::vector<int>::iterator i = indices.begin(); i != indices.end(); i++)
		{
			for(std::vector<int>::iterator j = indices.begin(); j != indices.end(); j++)
			{
				submatrix(std::distance(indices.begin(), i), std::distance(indices.begin(), j)) = psi(*i, *j);
			}
		}
	}
	double getDeterminantSign(boost::numeric::ublas::permutation_matrix<std::size_t>& pm)
	{
		int sign = 1;
		for(int i = 0; i < (int)pm.size(); i++)
		{
			if(i != (int)pm(i)) sign *= -1;
		}
		return sign;
	}
	double getDeterminant(boost::numeric::ublas::matrix<double>& square)
	{
		boost::numeric::ublas::permutation_matrix<std::size_t> pm(square.size1());
		double det = 1.0;
		if(boost::numeric::ublas::lu_factorize(square, pm)) return 0;

		for(int i = 0; i < (int)square.size1(); i++)
		{
			det *= square(i, i);
		}
		return det * getDeterminantSign(pm);
	}
	struct computeHHelper
	{
		computeHHelper(std::vector<mpfr_class>& multivariateGamma, boost::numeric::ublas::matrix<double>& psi, int delta, boost::numeric::ublas::matrix<double>& psiPart, int dimension, mpfr_class& numerator, mpfr_class& denominator)
			: multivariateGamma(multivariateGamma), psi(psi), delta(delta), psiPart(psiPart), dimension(dimension), numerator(numerator), denominator(denominator)
		{}
		void initialize_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void finish_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void start_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{
			//Vertex zero hasn't been handled yet.
			bitsetType vertexContents = boost::get(boost::vertex_name, g, v).contents;
			int vertexCount = vertexContents.count();
			mpfr_class& vertexGamma = multivariateGamma[vertexCount];
			if(vertexGamma == 0)
			{
				vertexGamma = multivariateGammaFunction(vertexCount, (delta + vertexCount - 1.0) / 2.0);
			}
			psiPart.resize(vertexCount, vertexCount, false);
			extractSubmatrix(psi, psiPart, vertexContents, dimension);
			double determinant = getDeterminant(psiPart);
			numerator *= (boost::multiprecision::pow(mpfr_class(determinant / (1ULL << vertexCount)), (delta + vertexCount - 1.0) / 2.0) / vertexGamma);
		}
		void examine_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void tree_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{
			int source = boost::source(v, g), target = boost::target(v, g);
			bitsetType sourceContents = boost::get(boost::vertex_name, g, source).contents;
			bitsetType targetContents = boost::get(boost::vertex_name, g, target).contents;
			bitsetType intersectionContents = sourceContents & targetContents;
			//Target goes into numerator
			int targetCount = targetContents.count();
			mpfr_class& targetGamma = multivariateGamma[targetCount];
			if(targetGamma == 0)
			{
				targetGamma = multivariateGammaFunction(targetCount, (delta + targetCount - 1.0) / 2.0);
			}
			psiPart.resize(targetCount, targetCount, false);
			extractSubmatrix(psi, psiPart, targetContents, dimension);
			double determinant = getDeterminant(psiPart);
			numerator *= (boost::multiprecision::pow(mpfr_class(determinant / (1ULL << targetCount)), (delta + targetCount - 1.0) / 2.0) / targetGamma);

			//Intersection goes into denominator
			int intersectionCount = intersectionContents.count();
			mpfr_class& intersectionGamma = multivariateGamma[intersectionCount];
			if(intersectionGamma == 0)
			{
				intersectionGamma = multivariateGammaFunction(intersectionCount, (delta + intersectionCount - 1.0) / 2.0);
			}
			psiPart.resize(intersectionCount, intersectionCount, false);
			extractSubmatrix(psi, psiPart, intersectionContents, dimension);
			determinant = getDeterminant(psiPart);
			denominator *= (boost::multiprecision::pow(mpfr_class(determinant / (1ULL << intersectionCount)), (delta + intersectionCount - 1.0) / 2.0) / intersectionGamma);
		}
		void back_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void forward_or_cross_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void discover_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		std::vector<mpfr_class>& multivariateGamma;
		boost::numeric::ublas::matrix<double>& psi;
		int delta;
		boost::numeric::ublas::matrix<double>& psiPart;
		int dimension;
		mpfr_class& numerator;
		mpfr_class& denominator;
	};
	mpfr_class h(cliqueTreeType& tree, int delta, boost::numeric::ublas::matrix<double>& psi, std::vector<mpfr_class>& multivariateGamma, boost::numeric::ublas::matrix<double>& psiPart, int dimension, std::vector<boost::default_color_type>& colourVector)
	{
		mpfr_class numerator = 1, denominator = 1;
		computeHHelper helper(multivariateGamma, psi, delta, psiPart, dimension, numerator, denominator);
		boost::iterator_property_map<std::vector<boost::default_color_type>::iterator, boost::identity_property_map> colourMap(colourVector.begin());
		std::fill(colourVector.begin(), colourVector.end(), boost::default_color_type::white_color);
		boost::depth_first_search(tree.getCliqueGraph(), helper, colourMap);
		
		return numerator / denominator;
	}
	void posteriorInferenceStep(cliqueTreeType& currentTree, posteriorInferenceArgs::graphType& graph, std::vector<mpfr_class>& exactValues, boost::mt19937& randomSource, workingPosteriorInference& working)
	{
		boost::random::bernoulli_distribution<> standardBernoulli;
		boost::random::uniform_real_distribution<> standardUniform;
		int original_edges = (int)boost::num_edges(graph);

		cliqueTreeAdjacencyMatrix& copiedTree = working.copiedTree;
		int nVertices = working.nVertices;
		int maxEdges = (nVertices * (nVertices - 1)) / 2;
		int delta = working.delta;
		int deltaStar = working.deltaStar;
		//Here we remove edges
		if(original_edges == maxEdges || (original_edges > 0 && standardBernoulli(randomSource)))
		{
			boost::random::uniform_int_distribution<> randomEdgeDist(0, (int)working.presentEdges.size()-1);
			int edgeIndex = randomEdgeDist(randomSource);
			int randomVertex1 = working.presentEdges[edgeIndex].first, randomVertex2 = working.presentEdges[edgeIndex].second;
			if(standardBernoulli(randomSource)) std::swap(randomVertex1, randomVertex2);

			int cliqueVertex = -1;
			//Can we remove this edge?
			if(currentTree.canRemoveEdge(randomVertex1, randomVertex2, working.counts1, cliqueVertex))
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
				//To add the bit relating to h(...), we need to create an update clique tree. 
				{
					copiedTree.makeCopy(currentTree);
					copiedTree.removeEdgeKnownCliqueVertex(randomVertex1, randomVertex2, working.colourVector, working.counts2, cliqueVertex);
					if(extraToRemove != 0)
					{
						//This is 1 rather than 0, because the edge randomVertex1, randomVertex2 is already deleted. 
						bitsetType copiedChosenSubset = chosenSubset;
						while(copiedChosenSubset.count() > 1)
						{
							copiedTree.canRemoveEdge(randomVertex1, randomVertex2, working.counts1, cliqueVertex);
							for(int i = 0; i < nVertices; i++)
							{
								if(copiedChosenSubset[i] && working.counts1[i] == 1)
								{
									copiedChosenSubset[i] = false;
									copiedTree.tryRemoveEdge(randomVertex1, i, working.colourVector, working.counts2);
								}
							}
						}
					}
				}
				mpfr_class extraFactor = (h(copiedTree, working.delta, working.psi, working.multivariateGammaDelta, working.psiPart, nVertices, working.colourVector) * h(currentTree, deltaStar, working.psiStar, working.multivariateGammaDeltaStar, working.psiPart, nVertices, working.colourVector)) / (h(currentTree, delta, working.psi, working.multivariateGammaDelta, working.psiPart, nVertices, working.colourVector) * h(copiedTree, deltaStar, working.psiStar, working.multivariateGammaDeltaStar, working.psiPart, nVertices, working.colourVector));
				acceptanceProbability *= extraFactor;
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
		//Here we add edges
		else
		{
			boost::random::uniform_int_distribution<> randomEdgeDist(0, (int)working.absentEdges.size()-1);
			int edgeIndex = randomEdgeDist(randomSource);
			int randomVertex1 = working.absentEdges[edgeIndex].first, randomVertex2 = working.absentEdges[edgeIndex].second;
			if(standardBernoulli(randomSource)) std::swap(randomVertex1, randomVertex2);

			cliqueTreeAdjacencyMatrix& copied2 = working.copied2;
			bitsetType newEdgesVertex1;
			currentTree.unionMinimalSeparators(randomVertex1, randomVertex2, newEdgesVertex1, working.vertexSequence, working.edgeSequence, working.addEdges, working.removeEdges, working.unionMinimalSepTemp);
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
				mpfr_class extraFactor = (h(copiedTree, working.delta, working.psi, working.multivariateGammaDelta, working.psiPart, nVertices, working.colourVector) * h(currentTree, deltaStar, working.psiStar, working.multivariateGammaDeltaStar, working.psiPart, nVertices, working.colourVector)) / (h(currentTree, delta, working.psi, working.multivariateGammaDelta, working.psiPart, nVertices, working.colourVector) * h(copiedTree, deltaStar, working.psiStar, working.multivariateGammaDeltaStar, working.psiPart, nVertices, working.colourVector));
				acceptanceProbability *= extraFactor;
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
#ifndef NDEBUG
		currentTree.check();
		assert(working.presentEdges.size() == boost::num_edges(graph));
		assert((int)working.presentEdges.size() + (int)working.absentEdges.size() == maxEdges);
#endif
	}
	bool operator==(const graphType& first, const graphType& second)
	{
		return std::equal(first.m_matrix.begin(), first.m_matrix.end(), second.m_matrix.begin());
	}
	void posteriorInference(posteriorInferenceArgs& args)
	{

		cliqueTreeType currentTree(args.dimension);
		posteriorInferenceArgs::graphType graph(args.dimension);
		for(std::size_t i = 0; i < args.dimension; i++) currentTree.addVertex();

		workingPosteriorInference working(args.dimension);
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
		for(std::size_t i = 0; i < args.sampleSize; i++)
		{
			posteriorInferenceStep(currentTree, graph, args.exactValues, args.randomSource, working);
			posteriorInferenceArgs::resultsType::iterator existingGraph = args.results.find(graph);
			if(existingGraph == args.results.end())
			{
				args.results[graph] = 1;
			}
			else existingGraph->second++;
		}
	}

}
