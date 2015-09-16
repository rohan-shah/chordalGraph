#include "convertGraph.h"
void convertGraph(SEXP graph_sexp, ::chordalGraph::cliqueTree::graphType& graphRef, R_GRAPH_TYPE graphType)
{
	if(graphType == IGRAPH)
	{
		convertGraphIGraph(graph_sexp, graphRef);
	}
	else if(graphType == GRAPHAM)
	{
		convertGraphAM(graph_sexp, graphRef);	
	}
	else if(graphType == GRAPHNEL)
	{
		convertGraphNEL(graph_sexp, graphRef);
	}
	else
	{
		throw std::runtime_error("Internal error");
	}
}