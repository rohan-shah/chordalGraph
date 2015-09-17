#include "isChordalIterative.h"
#include "convertGraph.h"
SEXP isChordalIterative(SEXP graph_sexp, R_GRAPH_TYPE type)
{
BEGIN_RCPP
	::chordalGraph::cliqueTree::graphType graph;
	convertGraph(graph_sexp, graph, type);

	int nVertices = (int)boost::num_vertices(graph);
	::chordalGraph::cliqueTree iterativeTree(nVertices);
	for (int i = 0; i < nVertices; i++)
	{
		::chordalGraph::bitsetType edges;
		for (int j = 0; j < i; j++)
		{
			if (boost::edge(i, j, graph).second)
			{
				edges[j] = true;
			}
		}
		bool isChordal = iterativeTree.tryAddVertexWithEdges(edges);
		if (!isChordal) return Rcpp::wrap(false);
	}
	return Rcpp::wrap(true);
END_RCPP
}
SEXP isChordalIterative_igraph(SEXP graph_sexp)
{
	return isChordalIterative(graph_sexp, IGRAPH);
}
SEXP isChordalIterative_graphNEL(SEXP graph_sexp)
{
	return isChordalIterative(graph_sexp, GRAPHNEL);
}
SEXP isChordalIterative_graphAM(SEXP graph_sexp)
{
	return isChordalIterative(graph_sexp, GRAPHAM);
}