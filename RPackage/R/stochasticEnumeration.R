stochasticEnumeration <- function(nVertices, budget, seed, nEdges, options = list(reduceChains = FALSE, graphRepresentation = "list"))
{
	if(missing(nVertices) || missing(budget) || missing(seed))
	{
		stop("Inputs nVertices, budget and seed are required")
	}
	if(length(nVertices) != 1 || length(budget) != 1 || length(seed) != 1)
	{
		stop("Inputs nVertices, budget and seed must be single integers")
	}
	if(!is.numeric(nVertices) || !is.numeric(budget) || !is.numeric(seed))
	{
		stop("Inputs nVertices, budget and seed must be single integers")
	}
	if(abs(nVertices - as.integer(nVertices)) > 1e-3 || abs(budget - as.integer(budget)) > 1e-3 || abs(seed - as.integer(seed)) > 1e-3)
	{
		stop("Inputs nVertices, budget and seed must be single integers")	
	}
	if(nVertices < 1 || budget < 1)
	{
		stop("Inputs nVertices and budget must be positive")
	}
	outputSamples <- FALSE
	if(!is.list(options))
	{
		stop("Input options must be a list")
	}
	if(any(!(names(options) %in% c("reduceChains", "outputSamples", "graphRepresentation"))))
	{
		stop("The only valid options for stochasticEnumeration are \"reduceChains\", \"outputSamples\" and \"graphRepresentation\"")
	}
	if(any(!(c("reduceChains", "graphRepresentation") %in% names(options))))
	{
		stop("Options \"reduceChains\" and \"graphRepresentation\" are required")
	}
	if(!(options$graphRepresentation %in% c("list", "matrix")))
	{
		stop("options$graphRepresentation must be either \"list\" or \"matrix\"")
	}
	if(missing(nEdges))
	{
		start <- Sys.time()
		result <- .Call("stochasticEnumeration", nVertices, budget, seed, options, PACKAGE="chordalGraph")
		end <- Sys.time()
		s4Result <- new("estimatedChordalCounts", data = result$data, call = match.call(), start = start, end = end, options = options, samples = NULL, exact = result$exact, minimumSizeForExact = result$minimumSizeForExact)
		if(outputSamples)
		{
			s4Result@samples <- result$samples
		}
		return(s4Result)
	}
	else
	{
		if(length(nEdges) != 1 || !is.numeric(nEdges) || abs(nEdges - as.integer(nEdges)) > 1e-3)
		{
			stop("Input nEdges must be a single integer")
		}
		if(nEdges > ((nVertices-1)*nVertices/2)+1 || nEdges < 0)
		{
			stop("Input nEdges must be in range [0, ((nVertices-1)*nVertices/2)+1]")
		}
		start <- Sys.time()
		result <- .Call("stochasticEnumerationSpecificEdges", nVertices, nEdges, budget, seed, options, PACKAGE="chordalGraph")
		end <- Sys.time()
		s4Result <- new("estimatedChordalCount", data = result$data, call = match.call(), start = start, end = end, options = options, samples = NULL, exact = result$exact, minimumSizeForExact = result$minimumSizeForExact)
		if(outputSamples)
		{
			s4Result@samples <- result$samples
		}
		return(s4Result)
	}
}
setClassUnion("listOrNULL", c("list", "NULL"))
setClass("estimatedChordalCount", slots = list(data = "mpfr", call = "call", start = "POSIXct", end = "POSIXct", options = "list", exact = "logical", samples = "listOrNULL", minimumSizeForExact = "integer"))
setClass("estimatedChordalCounts", slots = list(data = "mpfr", call = "call", start = "POSIXct", end = "POSIXct", options = "list", exact = "logical", samples = "listOrNULL", minimumSizeForExact = "integer"))
