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
#include "cholesky.hpp"
namespace chordalGraph
{
	mpfr_class getHRatio(cliqueTreeAdjacencyMatrix& tree, int cliqueVertex, int vertex1, int vertex2, boost::numeric::ublas::matrix<double>& psi, boost::numeric::ublas::matrix<double>& psiPart, int nVertices, int delta)
	{
		const double twoSqrtPi = 3.544907701811032054596334966682290365595098912244774256427;
		bitsetType cliqueVertexContents = boost::get(boost::vertex_name, tree.getCliqueGraph(), cliqueVertex).contents;
		int contentsCount = cliqueVertexContents.count();
		psiPart.resize(contentsCount, contentsCount, false);
		extractSubmatrixPair(psi, psiPart, cliqueVertexContents, nVertices, vertex1, vertex2);
		if(contentsCount != 2)
		{
			boost::numeric::ublasx::cholesky_decompose(psiPart);
			mpfr_class squareMatrixDet = psiPart(contentsCount - 1, contentsCount - 1) * psiPart(contentsCount - 1, contentsCount - 1) * psiPart(contentsCount - 2, contentsCount - 2)* psiPart(contentsCount - 2, contentsCount - 2);
			mpfr_class squareMatrixDetPower = boost::multiprecision::pow(squareMatrixDet, (delta + cliqueVertexContents.count() - 1)/2.0);
			mpfr_class otherDet1 = psiPart(contentsCount - 2, contentsCount - 2) * psiPart(contentsCount - 2, contentsCount - 2);
			mpfr_class otherDet1Power = boost::multiprecision::pow(otherDet1, (delta + cliqueVertexContents.count() - 2)/2.0);
			mpfr_class otherDet2 = psiPart(contentsCount-1, contentsCount - 2) * psiPart(contentsCount-1, contentsCount - 2) + psiPart(contentsCount-1, contentsCount-1) * psiPart(contentsCount-1, contentsCount-1);
			mpfr_class otherDet2Power = boost::multiprecision::pow(otherDet2, (delta + cliqueVertexContents.count() - 2)/2.0);
			mpfr_class ratio = (squareMatrixDetPower * boost::math::tgamma<mpfr_class>((delta + cliqueVertexContents.count() - 2)/2.0)) / (otherDet1Power * otherDet2Power * boost::math::tgamma<mpfr_class>((delta + cliqueVertexContents.count() - 1)/2.0) * twoSqrtPi);
			return ratio;
		}
		else
		{
			const double sqrtPi = 1.772453850905516027298167483341145182797549456122387128213;
			mpfr_class detSquare = psiPart(0, 0) * psiPart(1,1) - psiPart(1, 0) * psiPart(0, 1);
			mpfr_class detSquarePower = boost::multiprecision::pow(detSquare/4, (delta + 1)/2.0);
			mpfr_class detPowers = boost::multiprecision::pow(mpfr_class(psiPart(0,0) * psiPart(1, 1) / 4), delta/2.0);
			mpfr_class gamma1 = boost::math::tgamma<mpfr_class>(delta/2.0);
			mpfr_class gamma2 = sqrtPi * boost::math::tgamma<mpfr_class>((delta+1)/2.0) * boost::math::tgamma<mpfr_class>(delta/2.0);
			mpfr_class ratio = detSquarePower * gamma1 * gamma1 / (detPowers * gamma2);
			return ratio;
		}
	}
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
	void extractSubmatrixPair(boost::numeric::ublas::matrix<double>& psi, boost::numeric::ublas::matrix<double>& submatrix, bitsetType contents, int dimension, int first, int second)
	{
		std::vector<int> indices;
		for(int i = 0; i < dimension; i++)
		{
			if(contents[i] && i != first && i != second) indices.push_back(i);
		}
		indices.push_back(first);
		indices.push_back(second);
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
	void computeHHelper::start_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
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
	void computeHHelper::tree_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
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
	mpfr_class h(cliqueTreeType& tree, int delta, boost::numeric::ublas::matrix<double>& psi, std::vector<mpfr_class>& multivariateGamma, boost::numeric::ublas::matrix<double>& psiPart, int dimension, std::vector<boost::default_color_type>& colourVector)
	{
		mpfr_class numerator = 1, denominator = 1;
		computeHHelper helper(multivariateGamma, psi, delta, psiPart, dimension, numerator, denominator);
		boost::iterator_property_map<std::vector<boost::default_color_type>::iterator, boost::identity_property_map> colourMap(colourVector.begin());
		std::fill(colourVector.begin(), colourVector.end(), boost::default_color_type::white_color);
		boost::depth_first_search(tree.getCliqueGraph(), helper, colourMap);
		
		return numerator / denominator;
	}
	bool operator==(const graphType& first, const graphType& second)
	{
		return std::equal(first.m_matrix.begin(), first.m_matrix.end(), second.m_matrix.begin());
	}
}
