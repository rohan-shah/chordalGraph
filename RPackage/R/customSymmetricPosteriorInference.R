customSymmetricPosteriorInference <- function(data, delta, psi, exactCounts, seed, burnin, runSize)
{
	mean <- apply(data, 2, mean)
	differencesFromMean <- t(apply(data, 1, function(x) x - mean))
	outerProducts <- lapply(1:nrow(data), function(x) outer(differencesFromMean[x,], differencesFromMean[x,]))
	outerProductsSum <- Reduce('+', outerProducts)
	result <- .Call("customSymmetricPosteriorInference", outerProductsSum, delta, ncol(data), nrow(data), psi, as.character(exactCounts), burnin, runSize, PACKAGE="chordalGraph")
	return(result)
}
