context("Test exact data calculations are correct")

test_that("Exact data calculations are correct using stochastic enumeration",
{
	for(reduceChains in c(TRUE, FALSE))
	{
		capture.output(exact3 <- stochasticEnumeration(nVertices = 3, budget = 1000000, seed = 1, options = list(reduceChains = reduceChains)))
		expect_that(Rmpfr::all.equal(exact3@data, c(1,3,3,1)), is_true())

		capture.output(exact4 <- stochasticEnumeration(nVertices = 4, budget = 1000000, seed = 1, options = list(reduceChains = reduceChains)))
		expect_that(Rmpfr::all.equal(exact4@data, c(1,6,15,20,12,6,1)), is_true())

		capture.output(exact5 <- stochasticEnumeration(nVertices = 5, budget = 1000000, seed = 1, options = list(reduceChains = reduceChains)))
		expect_that(Rmpfr::all.equal(exact5@data, c(1,10,45,120,195,180,140,90,30,10,1)), is_true())

		capture.output(exact6 <- stochasticEnumeration(nVertices = 6, budget = 1000000, seed = 1, options = list(reduceChains = reduceChains)))
		expect_that(Rmpfr::all.equal(exact6@data, c(1,15,105,455,1320,2526,3085,3255,3000,2235,1206,615,260,60,15,1)), is_true())
	}
})
test_that("Exact data calculations are correct using stochastic enumeration nauty",
{
	for(reduceChains in c(TRUE, FALSE))
	{
		capture.output(exact3 <- stochasticEnumerationNauty(nVertices = 3, budget = 1000000, seed = 1, options = list(reduceChains = reduceChains)))
		expect_that(Rmpfr::all.equal(exact3@data, c(1,3,3,1)), is_true())

		capture.output(exact4 <- stochasticEnumerationNauty(nVertices = 4, budget = 1000000, seed = 1, options = list(reduceChains = reduceChains)))
		expect_that(Rmpfr::all.equal(exact4@data, c(1,6,15,20,12,6,1)), is_true())

		capture.output(exact5 <- stochasticEnumerationNauty(nVertices = 5, budget = 1000000, seed = 1, options = list(reduceChains = reduceChains)))
		expect_that(Rmpfr::all.equal(exact5@data, c(1,10,45,120,195,180,140,90,30,10,1)), is_true())

		capture.output(exact6 <- stochasticEnumerationNauty(nVertices = 6, budget = 1000000, seed = 1, options = list(reduceChains = reduceChains)))
		expect_that(Rmpfr::all.equal(exact6@data, c(1,15,105,455,1320,2526,3085,3255,3000,2235,1206,615,260,60,15,1)), is_true())
	}
})
test_that("Exact data calculations are correct using Horvitz Thompson estimation",
{
	for(method in chordalGraph:::samplingMethods)
	{
		capture.output(exact3 <- horvitzThompson(nVertices = 3, budget = 1000000, seed = 1, sampling = method))
		expect_that(Rmpfr::all.equal(exact3@data, c(1,3,3,1)), is_true())
	
		capture.output(exact4 <- horvitzThompson(nVertices = 4, budget = 1000000, seed = 1, sampling = method))
		expect_that(Rmpfr::all.equal(exact4@data, c(1,6,15,20,12,6,1)), is_true())
	
		capture.output(exact5 <- horvitzThompson(nVertices = 5, budget = 1000000, seed = 1, sampling = method))
		expect_that(Rmpfr::all.equal(exact5@data, c(1,10,45,120,195,180,140,90,30,10,1)), is_true())
	
		capture.output(exact6 <- horvitzThompson(nVertices = 6, budget = 1000000, seed = 1, sampling = method))
		expect_that(Rmpfr::all.equal(exact6@data, c(1,15,105,455,1320,2526,3085,3255,3000,2235,1206,615,260,60,15,1)), is_true())
	}
})
