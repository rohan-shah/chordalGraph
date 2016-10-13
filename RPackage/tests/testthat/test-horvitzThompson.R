context("Test horvitzThompson function")
test_that("Horvitz Thompson algorithm gives unbiased results for 5 x 5 graph",
{
	data(exact5, envir = environment())
	nReps <- 600
	results <- matrix(data=NA, nrow = nReps, ncol = 11)
	for(reduceChains in c(TRUE, FALSE))
	{
		for(i in 1:nReps)
		{
			capture.output(results[i,] <- as.numeric(horvitzThompson(nVertices = 5, seed = i, budget = 10, options = list(reduceChains = reduceChains, graphRepresentation = "matrix"))@data))
		}
		means <- apply(results, 2, mean)
		for(edgeCount in 1:11)
		{
			expect_equal(means[edgeCount], as.numeric(exact5@data[edgeCount]), tolerance = 0.01)
		}
	}
})
test_that("Horvitz Thompson algorithm gives identical results for different values of graphRepresentation",
{
	nReps <- 10
	for(reduceChains in c(TRUE, FALSE))
	{
		for(i in 1:nReps)
		{
			capture.output(resultList <- horvitzThompson(nVertices = 5, seed = i, budget = 10, options = list(reduceChains = reduceChains, graphRepresentation = "list")))
			capture.output(resultMatrix <- horvitzThompson(nVertices = 5, seed = i, budget = 10,  options = list(reduceChains = reduceChains, graphRepresentation = "matrix")))
			resultList@call <- resultMatrix@call <- call("list")
			resultList@options <- resultMatrix@options <- list()
			resultList@start <- resultList@end <- resultMatrix@start <- resultMatrix@end
			expect_identical(resultList, resultMatrix)
		}
	}
})

test_that("Horvitz Thompson algorithm gives unbiased results for 6 x 6 graph",
{
	data(exact6, envir = environment())
	nReps <- 600
	for(reduceChains in c(TRUE, FALSE))
	{
		results <- matrix(data=NA, nrow = nReps, ncol = 16)
		for(i in 1:nReps)
		{
			capture.output(results[i,] <- as.numeric(horvitzThompson(nVertices = 6, seed = i, budget = 30, options = list(reduceChains = reduceChains, graphRepresentation = "matrix"))@data))
		}
		means <- apply(results, 2, mean)
		for(edgeCount in 1:16)
		{
			expect_equal(means[edgeCount], as.numeric(exact6@data[edgeCount]), tolerance = 0.01)
		}
	}
})