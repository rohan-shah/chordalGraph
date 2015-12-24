stochasticEnumerationNauty <- function(nVertices, budget, seed, nEdges, options = list())
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
	if(!is.list(options))
	{
		stop("Input options must be a list")
	}
	if(any(!(names(options) %in% c("reduceChains"))))
	{
		stop("Input options contained invalid entries")
	}
	if(missing(nEdges))
	{
		start <- Sys.time()
		result <- .Call("stochasticEnumerationNauty", nVertices, budget, seed, options, PACKAGE="chordalGraph")
		end <- Sys.time()
		s4Result <- new("estimatedChordalCounts", data = result$data, call = match.call(), start = start, end = end, samples = NULL, options = options, exact = result$exact, minimumSizeForExact = result$minimumSizeForExact)
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
		result <- .Call("stochasticEnumerationNautySpecificEdges", nVertices, nEdges, budget, seed, options, PACKAGE="chordalGraph")
		end <- Sys.time()
		s4Result <- new("estimatedChordalCount", data = result$data, call = match.call(), start = start, end = end, samples = NULL, options = options, exact = result$exact, minimumSizeForExact = result$minimumSizeForExact)
		return(s4Result)
	}
}
