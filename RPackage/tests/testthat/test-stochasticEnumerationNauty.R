context("Test stochasticEnumerationNauty function")

test_that("Gives unbiased results for 6 x 6 graph",
{
	data(exact6, envir = environment())
	nReps <- 800
	for(reduceChains in c(TRUE, FALSE))
	{
		results <- matrix(data=NA, nrow = nReps, ncol = 16)
		for(i in 1:nReps)
		{
			capture.output(results[i,] <- as.numeric(stochasticEnumerationNauty(nVertices = 6, seed = i, budget = 45, options = list(reduceChains = reduceChains))@data))
		}
		means <- apply(results, 2, mean)
		for(edgeCount in 1:16)
		{
			expect_equal(means[edgeCount], as.numeric(exact6@data[edgeCount]), tolerance = 0.015)
		}
	}
})
