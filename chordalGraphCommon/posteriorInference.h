#ifndef POSTERIOR_INFERENCE_HEADER_GUARD
#define POSTERIOR_INFERENCE_HEADER_GUARD
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include "cliqueTree.h"
#include "numericType.h"
#include <boost/random/mersenne_twister.hpp>
#include "customMCMCSymmetric.h"
#include <unordered_map>
#include <boost/numeric/ublas/lu.hpp>
namespace chordalGraph
{
	mpfr_class multivariateGammaFunction(int m, double alpha);
	void extractSubmatrix(boost::numeric::ublas::matrix<double>& psi, boost::numeric::ublas::matrix<double>& submatrix, bitsetType contents, int dimension);
	double getDeterminantSign(boost::numeric::ublas::permutation_matrix<std::size_t>& pm);
	double getDeterminant(boost::numeric::ublas::matrix<double>& square);
	struct computeHHelper
	{
		computeHHelper(std::vector<mpfr_class>& multivariateGamma, boost::numeric::ublas::matrix<double>& psi, int delta, boost::numeric::ublas::matrix<double>& psiPart, int dimension, mpfr_class& numerator, mpfr_class& denominator)
			: multivariateGamma(multivariateGamma), psi(psi), delta(delta), psiPart(psiPart), dimension(dimension), numerator(numerator), denominator(denominator)
		{}
		void initialize_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void finish_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void start_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g);
		void examine_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void tree_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g);
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
	mpfr_class h(cliqueTreeType& tree, int delta, boost::numeric::ublas::matrix<double>& psi, std::vector<mpfr_class>& multivariateGamma, boost::numeric::ublas::matrix<double>& psiPart, int dimension, std::vector<boost::default_color_type>& colourVector);
	bool operator==(const graphType& first, const graphType& second);
}
#endif
