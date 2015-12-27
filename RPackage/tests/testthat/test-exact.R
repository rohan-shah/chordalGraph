context("Test exact data calculations are correct")
exactResults <- list()
exactResults[["3"]] <- c(1,3,3,1)
exactResults[["4"]] <- c(1,6,15,20,12,6,1)
exactResults[["5"]] <- c(1,10,45,120,195,180,140,90,30,10,1)
exactResults[["6"]] <- c(1,15,105,455,1320,2526,3085,3255,3000,2235,1206,615,260,60,15,1)

test_that("Exact data calculations are correct using stochastic enumeration",
{
	for(nVertices in as.character(3:6))
	{
		for(reduceChains in c(TRUE, FALSE))
		{
			capture.output(exactList <- stochasticEnumeration(nVertices = as.integer(nVertices), budget = 1000000, seed = 1, options = list(reduceChains = reduceChains, graphRepresentation = "list")))
			capture.output(exactMatrix <- stochasticEnumeration(nVertices = as.integer(nVertices), budget = 1000000, seed = 1, options = list(reduceChains = reduceChains, graphRepresentation = "matrix")))
			expect_that(Rmpfr::all.equal(exactList@data, exactResults[[nVertices]]), is_true())
			#Remove the bits that are different
			exactMatrix@call <- exactList@call <- call("list")
			exactMatrix@options <- exactList@options <- list()
			exactMatrix@start <- exactMatrix@end <- exactList@start <- exactList@end
			expect_identical(exactMatrix, exactList)
		}
	}
})
test_that("Exact data calculations are correct using stochastic enumeration nauty",
{
	startTime <- Sys.time()
	for(nVertices in as.character(3:6))
	{
		for(reduceChains in c(TRUE, FALSE))
		{
			capture.output(exactList <- stochasticEnumerationNauty(nVertices = as.integer(nVertices), budget = 1000000, seed = 1, options = list(reduceChains = reduceChains, graphRepresentation = "list")))
			capture.output(exactMatrix <- stochasticEnumerationNauty(nVertices = as.integer(nVertices), budget = 1000000, seed = 1, options = list(reduceChains = reduceChains, graphRepresentation = "matrix")))
			expect_that(Rmpfr::all.equal(exactList@data, exactResults[[nVertices]]), is_true())
			#Remove the bits that are different
			exactMatrix@call <- exactList@call <- call("list")
			exactMatrix@options <- exactList@options <- list()
			exactMatrix@start <- exactMatrix@end <- exactList@start <- exactList@end
			expect_identical(exactMatrix, exactList)
		}
	}
})
test_that("Exact data calculations are correct using Horvitz Thompson estimation",
{
	for(nVertices in as.character(3:6))
	{
		for(reduceChains in c(TRUE, FALSE))
		{
			for(method in chordalGraph:::samplingMethods)
			{
				capture.output(exactList <- horvitzThompson(nVertices = as.integer(nVertices), budget = 1000000, seed = 1, sampling = method, options = list(reduceChains = reduceChains, graphRepresentation = "list")))
				capture.output(exactMatrix <- horvitzThompson(nVertices = as.integer(nVertices), budget = 1000000, seed = 1, sampling = method, options = list(reduceChains = reduceChains, graphRepresentation = "matrix")))
				expect_that(Rmpfr::all.equal(exactList@data, exactResults[[nVertices]]), is_true())
				#Remove the bits that are different
				exactMatrix@call <- exactList@call <- call("list")
				exactMatrix@options <- exactList@options <- list()
				exactMatrix@start <- exactMatrix@end <- exactList@start <- exactList@end
				expect_identical(exactMatrix, exactList)
			}
		}
	}
})
rm(exactResults)
