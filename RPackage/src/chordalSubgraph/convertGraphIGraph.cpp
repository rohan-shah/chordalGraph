#include "convertGraph.h"
void convertGraphIGraph(SEXP graph_sexp, ::chordalSubgraph::cliqueTree::graphType& graphRef)
{
	//Convert graph object
	Rcpp::List graph;
	try
	{
		graph = Rcpp::as<Rcpp::List>(graph_sexp);
	}
	catch(Rcpp::not_compatible&)
	{
		throw std::runtime_error("Unable to convert input graph to a list");
	}
	int nVertices = Rcpp::as<int>(graph(0));
	bool directed = Rcpp::as<bool>(graph(1));
	Rcpp::IntegerVector edgesVertex1 = Rcpp::as<Rcpp::IntegerVector>(graph(2));
	Rcpp::IntegerVector edgesVertex2 = Rcpp::as<Rcpp::IntegerVector>(graph(3));
	if(edgesVertex1.size() != edgesVertex2.size())
	{
		throw std::runtime_error("Internal error");
	}
	if(directed)
	{
		throw std::runtime_error("Input graph must be undirected");
	}

	//Construct graph
	graphRef = ::chordalSubgraph::cliqueTree::graphType(nVertices);
	for(int i = 0; i < edgesVertex1.size(); i++)
	{
		boost::add_edge(edgesVertex1(i), edgesVertex2(i), graphRef);
	}
}
