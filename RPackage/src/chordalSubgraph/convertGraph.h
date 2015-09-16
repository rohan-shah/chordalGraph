#ifndef CONVERT_GRAPH_HEADER_GUARD
#define CONVERT_GRAPH_HEADER_GUARD
#include "cliqueTree.h"
#include <Rcpp.h>
enum R_GRAPH_TYPE
{
	IGRAPH, GRAPHNEL, GRAPHAM
};
void convertGraph(SEXP graph_sexp, ::chordalSubgraph::cliqueTree::graphType& graphRef, R_GRAPH_TYPE graphType);
void convertGraphIGraph(SEXP graph_sexp, ::chordalSubgraph::cliqueTree::graphType& graphRef);
void convertGraphAM(SEXP graph_sexp, ::chordalSubgraph::cliqueTree::graphType& graphRef);
void convertGraphNEL(SEXP graph_sexp, ::chordalSubgraph::cliqueTree::graphType& graphRef);
#endif