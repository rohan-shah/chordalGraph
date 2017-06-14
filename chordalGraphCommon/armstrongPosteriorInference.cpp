#include "armstrongPosteriorInference.h"
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
namespace chordalGraph
{
	void armstrongPosteriorInferenceStep(cliqueTreeType& currentTree, armstrongPosteriorInferenceArgs::graphType& graph, std::vector<mpfr_class>& exactValues, boost::mt19937& randomSource, workingArmstrongPosteriorInference& working)
	{
		cliqueTreeAdjacencyMatrix& copiedTree = working.copiedTree;
		int delta = working.delta;
		int deltaStar = working.deltaStar;
		int nVertices = working.nVertices;
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
			int cliqueVertex = -1;
			if(currentTree.canRemoveEdge(randomVertex1, randomVertex2, working.counts1, cliqueVertex))
			{
				copiedTree.makeCopy(currentTree);
				cliqueTreeAdjacencyMatrix::removeReversal reverse;
				copiedTree.removeEdgeKnownCliqueVertex(randomVertex1, randomVertex2, working.colourVector, working.counts2, cliqueVertex, reverse);
				mpfr_class acceptanceValue = (exactValues[original_edges] / exactValues[original_edges - 1]);
				//mpfr_class extraFactor = (h(copiedTree, delta, working.psi, working.multivariateGammaDelta, working.psiPart, nVertices, working.colourVector) * h(currentTree, deltaStar, working.psiStar, working.multivariateGammaDeltaStar, working.psiPart, nVertices, working.colourVector)) / (h(currentTree, delta, working.psi, working.multivariateGammaDelta, working.psiPart, nVertices, working.colourVector) * h(copiedTree, deltaStar, working.psiStar, working.multivariateGammaDeltaStar, working.psiPart, nVertices, working.colourVector));
				mpfr_class extraFactor2 = getHRatio(currentTree, cliqueVertex, randomVertex1, randomVertex2, working.psi, working.psiPart, nVertices, delta) / getHRatio(currentTree, cliqueVertex, randomVertex1, randomVertex2, working.psiStar, working.psiPart, nVertices, deltaStar);
				/*if(std::fabs(1 - (extraFactor * extraFactor2).convert_to<double>()) > 1e-5)
				{
					throw std::runtime_error("Internal error");
				}
				acceptanceValue *= extraFactor;*/
				acceptanceValue /= extraFactor2;
				if(acceptanceValue >= 1 || standardUniform(randomSource) <= acceptanceValue.convert_to<double>())
				{
					currentTree.swap(copiedTree);
					boost::remove_edge(randomVertex1, randomVertex2, graph);
				}
			}
		}
		//Here we add edges
		else
		{
			bitsetType newEdgesVertex1;
			currentTree.unionMinimalSeparators(randomVertex1, randomVertex2, newEdgesVertex1, working.vertexSequence, working.edgeSequence, working.addEdges, working.removeEdges, working.unionMinimalSepTemp);
			int increaseInEdges = 1;
			for(int i = 0; i < nVertices; i++)
			{
				if(newEdgesVertex1[i] && !boost::edge(i, randomVertex1, graph).second) increaseInEdges++;
			}
			if(increaseInEdges == 1)
			{
				copiedTree.makeCopy(currentTree);
				copiedTree.addEdge(randomVertex1, randomVertex2, newEdgesVertex1, working.vertexSequence, working.edgeSequence, working.addEdges, working.removeEdges, working.unionMinimalSepTemp, true);
				//Lazy way of working out clique vertex. 
				int cliqueVertex = -1;
				copiedTree.canRemoveEdge(randomVertex1, randomVertex2, working.counts1, cliqueVertex);
				mpfr_class acceptanceValue = exactValues[original_edges] / exactValues[original_edges + increaseInEdges];
				//mpfr_class extraFactor = (h(copiedTree, delta, working.psi, working.multivariateGammaDelta, working.psiPart, nVertices, working.colourVector) * h(currentTree, deltaStar, working.psiStar, working.multivariateGammaDeltaStar, working.psiPart, nVertices, working.colourVector)) / (h(currentTree, delta, working.psi, working.multivariateGammaDelta, working.psiPart, nVertices, working.colourVector) * h(copiedTree, deltaStar, working.psiStar, working.multivariateGammaDeltaStar, working.psiPart, nVertices, working.colourVector));
				mpfr_class extraFactor2 = getHRatio(copiedTree, cliqueVertex, randomVertex1, randomVertex2, working.psi, working.psiPart, nVertices, delta) / getHRatio(copiedTree, cliqueVertex, randomVertex1, randomVertex2, working.psiStar, working.psiPart, nVertices, deltaStar);
				/*if(std::fabs(1 - (extraFactor / extraFactor2).convert_to<double>()) > 1e-5)
				{
					throw std::runtime_error("Internal error");
				}
				acceptanceValue *= extraFactor;*/
				acceptanceValue *= extraFactor2;
				if(acceptanceValue >= 1 || standardUniform(randomSource) <= acceptanceValue.convert_to<double>())
				{
					currentTree.swap(copiedTree);
					boost::add_edge(randomVertex1, randomVertex2, graph);
				}
			}
		}
	}
	void armstrongPosteriorInference(armstrongPosteriorInferenceArgs& args)
	{

		cliqueTreeType currentTree(args.dimension);
		armstrongPosteriorInferenceArgs::graphType graph(args.dimension);
		for(std::size_t i = 0; i < args.dimension; i++) currentTree.addVertex();

		workingArmstrongPosteriorInference working(args.dimension);
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
			armstrongPosteriorInferenceStep(currentTree, graph, args.exactValues, args.randomSource, working);
		}
		args.results.clear();
		for(std::size_t i = 0; i < args.sampleSize; i++)
		{
			armstrongPosteriorInferenceStep(currentTree, graph, args.exactValues, args.randomSource, working);
			armstrongPosteriorInferenceArgs::resultsType::iterator existingGraph = args.results.find(graph);
			if(existingGraph == args.results.end())
			{
				args.results[graph] = 1;
			}
			else existingGraph->second++;
		}
	}

}
