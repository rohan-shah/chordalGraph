library(chordalGraph)
vertices <- 34
maxEdges <- vertices * (vertices - 1)/2
exactCounts <- mpfr(c(1, choose(maxEdges, 1), choose(maxEdges, 2), choose(maxEdges, 3), choose(maxEdges, 4) - choose(vertices, 4)*3, choose(maxEdges, 5) - choose(vertices, 5) * 12 - choose(vertices, 4)*3*(maxEdges - 6)), precBits = 64)
currentCounts <- exactCounts
sampleSize <- 1000000
burnin <- 10000
maxEdges <- 561

#Run iterative MCMC
currentCounts <- exactCounts
customSymmetricCounts <- list()
for(edgeCounter in 6:maxEdges)
{
	last <- tail(currentCounts, 1)
	secondLast <- tail(currentCounts, 2)[1]
	estimatedCounts <- c(currentCounts, 0.5 * (last^2 / secondLast))
	estimatedCounts[1:6] <- exactCounts
	customSymmetricResults <- customMCMCSymmetric(vertices, estimatedCounts, seed = edgeCounter, burnin, sampleSize)
	customSymmetricCounts[[edgeCounter]] <- customSymmetricResults$estimates
	currentCounts <- customSymmetricResults$estimates#c(currentCounts, tail(customResults$estimates, 1))
	if(any(currentCounts == 0)) stop("Estimated a value of 0")
	cat(edgeCounter, " / ", maxEdges, "\n", sep="")
}

save(customSymmetricCounts, file = "results.RData")
