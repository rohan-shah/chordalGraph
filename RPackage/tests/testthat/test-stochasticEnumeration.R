context("Test stochasticEnumeration function")

test_that("Stochastic enumeration function gives unbiased results for 6 x 6 graph",
{
	data(exact6, envir = environment())
	for(reduceChains in c(TRUE, FALSE))
	{
		nReps <- 1000
		results <- matrix(data=NA, nrow = nReps, ncol = 16)
		for(i in 1:nReps)
		{
			capture.output(results[i,] <- as.numeric(stochasticEnumeration(nVertices = 6, seed = i, budget = 45, options = list(reduceChains = reduceChains, graphRepresentation = "list"))@data))
		}
		means <- apply(results, 2, mean)
		for(edgeCount in 1:16)
		{
			expect_equal(means[edgeCount], as.numeric(exact6@data[edgeCount]), tolerance = 0.015)
		}
	}
})
test_that("Stochastic enumeration function gives the same results for different values of graphRepresentation",
{
	for(reduceChains in c(TRUE, FALSE))
	{
		nReps <- 10
		for(i in 1:nReps)
		{
			capture.output(resultList <- stochasticEnumeration(nVertices = 6, seed = i, budget = 45, options = list(reduceChains = reduceChains, graphRepresentation = "list")))
			capture.output(resultMatrix <- stochasticEnumeration(nVertices = 6, seed = i, budget = 45, options = list(reduceChains = reduceChains, graphRepresentation = "list")))
			resultList@call <- resultMatrix@call <- call("list")
			resultList@options <- resultMatrix@options <- list()
			resultList@start <- resultList@end <- resultMatrix@start <- resultMatrix@end
			expect_identical(resultList, resultMatrix)
		}
	}
})
