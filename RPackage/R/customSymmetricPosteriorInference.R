setClass("customPosterior", slots = list(data = "ANY", delta = "numeric", psi = "matrix", exactCounts = "character", call = "call", start = "POSIXct", end = "POSIXct", seed = "integer", burnin = "integer", runSize = "integer", probabilities = "ANY", graphs = "ANY", uniqueGraphsLimit = "integer"))
customSymmetricPosteriorInference <- function(data, delta, psi, exactCounts, seed, burnin, runSize = 0, uniqueGraphsLimit = 0)
{
	if((uniqueGraphsLimit <= 0 && runSize <= 0) || (uniqueGraphsLimit > 0 && runSize > 0))
	{
		stop("Only one of inputs uniqueGraphsLimit and runSize can be specified")
	}
	mean <- apply(data, 2, mean)
	differencesFromMean <- t(apply(data, 1, function(x) x - mean))
	outerProducts <- lapply(1:nrow(data), function(x) outer(differencesFromMean[x,], differencesFromMean[x,]))
	outerProductsSum <- Reduce('+', outerProducts)
	if(class(exactCounts) == "mpfr")
	{
		exactCounts <- format(exactCounts)
	}
	start <- Sys.time()
	result <- .Call("customSymmetricPosteriorInference", outerProductsSum, delta, ncol(data), nrow(data), psi, as.character(exactCounts), burnin, runSize, seed, uniqueGraphsLimit, PACKAGE="chordalGraph")
	end <- Sys.time()
	return(new("customPosterior", data = data, delta = delta, psi = psi, exactCounts = as.character(exactCounts), call = match.call(), start = start, end = end, seed = as.integer(seed), burnin = as.integer(burnin), runSize = as.integer(result$sampleSize), probabilities = result$probabilities, graphs = result$graphs, uniqueGraphsLimit = as.integer(uniqueGraphsLimit)))
}
