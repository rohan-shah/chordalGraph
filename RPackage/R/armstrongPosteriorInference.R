setClass("armstrongPosterior", slots = list(data = "ANY", delta = "numeric", psi = "matrix", exactCounts = "character", call = "call", start = "POSIXct", end = "POSIXct", seed = "integer", burnin = "integer", runSize = "integer", probabilities = "ANY", graphs = "ANY"))
armstrongPosteriorInference <- function(data, delta, psi, exactCounts, seed, burnin, runSize)
{
	mean <- apply(data, 2, mean)
	differencesFromMean <- t(apply(data, 1, function(x) x - mean))
	outerProducts <- lapply(1:nrow(data), function(x) outer(differencesFromMean[x,], differencesFromMean[x,]))
	outerProductsSum <- Reduce('+', outerProducts)
	if(class(exactCounts) == "mpfr")
	{
		exactCounts <- format(exactCounts)
	}
	start <- Sys.time()
	result <- .Call("armstrongPosteriorInference", outerProductsSum, delta, ncol(data), nrow(data), psi, as.character(exactCounts), burnin, runSize, seed, PACKAGE="chordalGraph")
	end <- Sys.time()
	return(new("armstrongPosterior", data = data, delta = delta, psi = psi, exactCounts = as.character(exactCounts), call = match.call(), start = start, end = end, seed = as.integer(seed), burnin = as.integer(burnin), runSize = as.integer(runSize), probabilities = result$probabilities, graphs = result$graphs))
}
