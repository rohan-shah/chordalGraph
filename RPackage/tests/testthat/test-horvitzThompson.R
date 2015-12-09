context("Test horvitzThompson function")
test_that("ConditionalPoisson method gives unbiased results for 5 x 5 graph",
{
	data(exact5, envir = environment())
	nReps <- 400
	results <- matrix(data=NA, nrow = nReps, ncol = 11)
	for(i in 1:nReps)
	{
		capture.output(results[i,] <- as.numeric(horvitzThompson(nVertices = 5, seed = i, budget = 10, sampling = "conditionalPoisson")@data))
	}
	means <- apply(results, 2, mean)
	for(edgeCount in 1:11)
	{
		expect_equal(means[edgeCount], as.numeric(exact5@data[edgeCount]), tolerance = 0.025)
	}
})
test_that("All methods except sampfordMultinomial, conditionalPoisson and pareto give unbiased results for 6 x 6 graph",
{
	data(exact6, envir = environment())
	nReps <- 400
	methods <- setdiff(chordalGraph::samplingMethods, c("sampfordMultinomial", "conditionalPoisson", "pareto"))
	for(method in methods)
	{
		results <- matrix(data=NA, nrow = nReps, ncol = 16)
		for(i in 1:nReps)
		{
			capture.output(results[i,] <- as.numeric(horvitzThompson(nVertices = 6, seed = i, budget = 40, sampling = method)@data))
		}
		means <- apply(results, 2, mean)
		for(edgeCount in 1:16)
		{
			expect_equal(means[edgeCount], as.numeric(exact6@data[edgeCount]), tolerance = 0.025)
		}
	}
})