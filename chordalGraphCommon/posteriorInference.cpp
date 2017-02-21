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
			result *= boost::math::tgamma(alpha - (i - 1.0)/2.0);
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
		computeHHelper(boost::numeric::ublas::matrix<mpfr_class>& multivariateGamma, boost::numeric::ublas::matrix<double>& psi, int delta, boost::numeric::ublas::matrix<double>& temp, int dimension, mpfr_class& numerator, mpfr_class& denominator)
			: multivariateGamma(multivariateGamma), psi(psi), delta(delta), temp(temp), dimension(dimension), numerator(numerator), denominator(denominator)
		{}
		void initialize_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void finish_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void start_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
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
			mpfr_class& targetGamma = multivariateGamma(targetCount, (delta + targetCount - 1.0) / 2.0);
			if(targetGamma == 0)
			{
				targetGamma = multivariateGammaFunction(targetCount, (delta + targetCount - 1.0) / 2.0);
			}
			temp.resize(targetCount, targetCount, false);
			extractSubmatrix(psi, temp, targetContents, dimension);
			double determinant = getDeterminant(temp);
			numerator *= (boost::multiprecision::pow(mpfr_class(determinant), (delta + targetCount - 1.0) / 2.0) / targetGamma);

			//Intersection goes into denominator
			int intersectionCount = intersectionContents.count();
			mpfr_class& intersectionGamma = multivariateGamma(intersectionCount, (delta + intersectionCount - 1.0) / 2.0);
			if(intersectionGamma == 0)
			{
				intersectionGamma = multivariateGammaFunction(intersectionCount, (delta + intersectionCount - 1.0) / 2.0);
			}
			temp.resize(intersectionCount, intersectionCount, false);
			extractSubmatrix(psi, temp, intersectionContents, dimension);
			determinant = getDeterminant(temp);
			denominator *= (boost::multiprecision::pow(mpfr_class(determinant), (delta + intersectionCount - 1.0) / 2.0) / intersectionGamma);
		}
		void back_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void forward_or_cross_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void discover_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		boost::numeric::ublas::matrix<mpfr_class>& multivariateGamma;
		boost::numeric::ublas::matrix<double>& psi;
		int delta;
		boost::numeric::ublas::matrix<double>& temp;
		int dimension;
		mpfr_class& numerator;
		mpfr_class& denominator;
	};
	mpfr_class h(cliqueTreeType& tree, int delta, boost::numeric::ublas::matrix<double>& psi, boost::numeric::ublas::matrix<mpfr_class>& multivariateGamma, boost::numeric::ublas::matrix<double>& temp, int dimension, std::vector<boost::default_color_type>& colourVector)
	{
		mpfr_class numerator = 1, denominator = 1;
		computeHHelper helper(multivariateGamma, psi, delta, temp, dimension, numerator, denominator);
		boost::iterator_property_map<std::vector<boost::default_color_type>::iterator, boost::identity_property_map> colourMap(colourVector.begin());
		std::fill(colourVector.begin(), colourVector.end(), boost::default_color_type::white_color);
		boost::depth_first_visit(tree.getCliqueGraph(), 0, helper, colourMap);
		
		//Vertex zero hasn't been handled yet.
		bitsetType vertexZeroContents = boost::get(boost::vertex_name, tree.getCliqueGraph(), 0).contents;
		int vertexZeroCount = vertexZeroContents.count();
		mpfr_class& vertexZeroGamma = multivariateGamma(vertexZeroCount, (delta + vertexZeroCount - 1.0) / 2.0);
		if(vertexZeroGamma == 0)
		{
			vertexZeroGamma = multivariateGammaFunction(vertexZeroCount, (delta + vertexZeroCount - 1.0) / 2.0);
		}
		temp.resize(vertexZeroCount, vertexZeroCount, false);
		extractSubmatrix(psi, temp, vertexZeroContents, dimension);
		double determinant = getDeterminant(temp);
		numerator *= (boost::multiprecision::pow(mpfr_class(determinant), (delta + vertexZeroCount - 1.0) / 2.0) / vertexZeroGamma);
		return numerator / denominator;
	}
	void posteriorInferenceStep(boost::numeric::ublas::matrix<mpfr_class>& multivariateGamma, boost::numeric::ublas::matrix<double>& psiStar, int deltaStar, cliqueTreeType& currentTree, posteriorInferenceArgs::graphType& graph, std::vector<mpfr_class>& exactValues, int nVertices, boost::mt19937& randomSource, workingCustomSymmetric& temp)
	{
		boost::random::bernoulli_distribution<> standardBernoulli;
		boost::random::uniform_real_distribution<> standardUniform;
		int original_edges = (int)boost::num_edges(graph);

		cliqueTreeAdjacencyMatrix& copiedTree = temp.copiedTree;
		int maxEdges = (nVertices * (nVertices - 1)) / 2;
		//Here we remove edges
		if(original_edges == maxEdges || (original_edges > 0 && standardBernoulli(randomSource)))
		{
			boost::random::uniform_int_distribution<> randomEdgeDist(0, (int)temp.presentEdges.size()-1);
			int edgeIndex = randomEdgeDist(randomSource);
			int randomVertex1 = temp.presentEdges[edgeIndex].first, randomVertex2 = temp.presentEdges[edgeIndex].second;
			if(standardBernoulli(randomSource)) std::swap(randomVertex1, randomVertex2);

			int cliqueVertex = -1;
			//Can we remove this edge?
			if(currentTree.canRemoveEdge(randomVertex1, randomVertex2, temp.counts1, cliqueVertex))
			{
				//Form tree of removable edge subsets and work out the counts. 
				currentTree.formRemovalTree(temp.stateCounts, copiedTree, randomVertex1, randomVertex2, temp.uniqueSubsets, temp.removalTemporaries);
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
					currentTree.formRemovalTree(temp.stateCounts, copiedTree, randomVertex2, randomVertex1, temp.uniqueSubsets, temp.removalTemporaries);
					
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
						acceptanceProbability = 0.5*((double)original_edges/(double)(maxEdges - (int)(original_edges - 1 - extraToRemove))) * sum1 * temp.stateCounts[extraToRemove];
					}
					else
					{
						acceptanceProbability = ((double)original_edges/(double)(maxEdges - (int)(original_edges - 1 - extraToRemove))) * sum1 * temp.stateCounts[extraToRemove];
					}
				}
				//To add the bit relating to h(...), we need to create an update clique tree. 
				posteriorInferenceArgs::graphType& copiedGraph = temp.copiedGraph;
				copiedTree.makeCopy(currentTree);
				copiedTree.removeEdgeKnownCliqueVertex(randomVertex1, randomVertex2, temp.colourVector, temp.counts2, cliqueVertex);
				boost::remove_edge(randomVertex1, randomVertex2, copiedGraph);
				{
					std::swap(temp.presentEdges[edgeIndex], temp.presentEdges.back());
					temp.presentEdges.pop_back();
					int minVertex = std::min(randomVertex1, randomVertex2), maxVertex = std::max(randomVertex1, randomVertex2);
					temp.absentEdges.push_back(std::make_pair(minVertex, maxVertex));
				}
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
						copiedTree.canRemoveEdge(randomVertex1, randomVertex2, temp.counts1, cliqueVertex);
						for(int i = 0; i < nVertices; i++)
						{
							if(chosenSubset[i] && temp.counts1[i] == 1)
							{
								chosenSubset[i] = false;
								currentTree.tryRemoveEdge(randomVertex1, i, temp.colourVector, temp.counts2);
								boost::remove_edge(randomVertex1, i, copiedGraph);
								{
									int minVertex = std::min(randomVertex1, i), maxVertex = std::max(randomVertex1, i);
									int edgeCounter = (int)std::distance(temp.presentEdges.begin(), std::find(temp.presentEdges.begin(), temp.presentEdges.end(), std::make_pair(minVertex, maxVertex)));
									std::swap(temp.presentEdges[edgeCounter], temp.presentEdges.back());
									temp.presentEdges.pop_back();
									temp.absentEdges.push_back(std::make_pair(minVertex, maxVertex));
								}
							}
						}
					}
				}
				if(acceptanceProbability >= 1 || standardUniform(randomSource) <= acceptanceProbability.convert_to<double>())
				{
					graph.swap(copiedGraph);
					currentTree.swap(copiedTree);
				}
			}
		}
		//Here we add edges
		else
		{
			boost::random::uniform_int_distribution<> randomEdgeDist(0, (int)temp.absentEdges.size()-1);
			int edgeIndex = randomEdgeDist(randomSource);
			int randomVertex1 = temp.absentEdges[edgeIndex].first, randomVertex2 = temp.absentEdges[edgeIndex].second;
			if(standardBernoulli(randomSource)) std::swap(randomVertex1, randomVertex2);

			cliqueTreeAdjacencyMatrix& copied2 = temp.copied2;
			bitsetType newEdgesVertex1;
			currentTree.unionMinimalSeparators(randomVertex1, randomVertex2, newEdgesVertex1, temp.vertexSequence, temp.edgeSequence, temp.addEdges, temp.removeEdges, temp.unionMinimalSepTemp);
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
					if(increaseInEdges + original_edges == maxEdges)
					{
						acceptanceProbability = 2.0/(((double)(original_edges+increaseInEdges)/(double)(maxEdges - (int)original_edges)) * sum * temp.stateCounts[increaseInEdges - 1]);
					}
					else
					{
						acceptanceProbability = 1.0/(((double)(original_edges+increaseInEdges)/(double)(maxEdges - (int)original_edges)) * sum * temp.stateCounts[increaseInEdges - 1]);
					}
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
				if (acceptanceProbability >= 1 || standardUniform(randomSource) <= acceptanceProbability.convert_to<double>())
				{
					currentTree.swap(copied);
					newEdgesVertex1[randomVertex2] = true;
					for(int i = 0; i < nVertices; i++)
					{
						if(newEdgesVertex1[i] && i != randomVertex1)
						{
							bool addedNewEdge = boost::add_edge(randomVertex1, i, graph).second;
							if(addedNewEdge)
							{
								int minVertex = std::min(randomVertex1, i), maxVertex = std::max(randomVertex1, i);
								int edgeCounter = (int)std::distance(temp.absentEdges.begin(), std::find(temp.absentEdges.begin(), temp.absentEdges.end(), std::make_pair(minVertex, maxVertex)));
								std::swap(temp.absentEdges[edgeCounter], temp.absentEdges.back());
								temp.absentEdges.pop_back();
								temp.presentEdges.push_back(std::make_pair(minVertex, maxVertex));
							}
						}
					}
				}
			}
		}
#ifndef NDEBUG
		currentTree.check();
		assert(temp.presentEdges.size() == boost::num_edges(graph));
		assert((int)temp.presentEdges.size() + (int)temp.absentEdges.size() == maxEdges);
#endif
	}
	void posteriorInference(posteriorInferenceArgs& args)
	{
		boost::numeric::ublas::matrix<mpfr_class> multivariateGamma(args.dimension, (args.delta + args.dimension)/2, 0);
		boost::numeric::ublas::matrix<double> psiStar = args.psi + args.sampleCovariance;
		int deltaStar = args.delta + args.dataPoints - 1;

		cliqueTreeType currentTree(args.dimension);
		posteriorInferenceArgs::graphType graph(args.dimension);
		for(std::size_t i = 0; i < args.dimension; i++) currentTree.addVertex();

		workingCustomSymmetric temp(args.dimension);
		for(int i = 0; i < (int)args.dimension; i++)
		{
			for(int j = i + 1; j < (int)args.dimension; j++) temp.absentEdges.push_back(std::make_pair(i, j));
		}
		
		for(std::size_t i = 0; i < args.burnIn; i++)
		{
			posteriorInferenceStep(multivariateGamma, psiStar, deltaStar, currentTree, graph, args.exactValues, args.dimension, args.randomSource, temp);
		}
		for(std::size_t i = 0; i < args.sampleSize; i++)
		{
//			posteriorInferenceStep();
		}
	}
}
