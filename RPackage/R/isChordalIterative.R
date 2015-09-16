isChordalIterative <- function(graph)
{
	if(class(graph) == "igraph")
	{
		if(igraph::is.directed(graph))
		{
			stop("Input `graph' must be undirected")
		}
		estimate <- .Call("isChordalIterative_igraph", graph, PACKAGE="chordalSubgraph")
	}
	else if(class(graph) == "graphNEL")
	{
		estimate <- .Call("isChordalIterative_graphNEL", graph, PACKAGE="chordalSubgraph")
	}
	else if(class(graph) == "graphAM")
	{
		estimate <- .Call("isChordalIterative_graphAM", graph, PACKAGE="chordalSubgraph")
	}
	else 
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
}