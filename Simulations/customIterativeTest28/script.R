library(chordalGraph)
vertices <- as.integer(Sys.getenv("VERTICES"))
cat("vertices = ", vertices, "\n", sep="")

maxEdges <- vertices * (vertices - 1)/2
exactCounts <- mpfr(c(1, choose(maxEdges, 1), choose(maxEdges, 2), choose(maxEdges, 3), choose(maxEdges, 4) - choose(vertices, 4)*3, choose(maxEdges, 5) - choose(vertices, 5) * 12 - choose(vertices, 4)*3*(maxEdges - 6)), precBits = 64)
currentCounts <- exactCounts
sampleSize <- 2000000
burnin <- 10000

filename <- paste0("results", vertices, ".RData")

#Run iterative MCMC
if(!file.exists(filename))
{
	currentCounts <- exactCounts
	customSymmetricCounts <- list()
	start <- 6
} else
{
	load(filename)
	start <- length(customSymmetricCounts)+1
	currentCounts <- customSymmetricCounts[[length(customSymmetricCounts)]]
	if(start == maxEdges+1) stop("Already finished")
}
for(edgeCounter in start:maxEdges)
{
	last <- tail(currentCounts, 1)
	secondLast <- tail(currentCounts, 2)[1]
	estimatedCounts <- c(currentCounts, 0.5 * (last^2 / secondLast))
	estimatedCounts[1:6] <- exactCounts
	customSymmetricResults <- customMCMCSymmetric(vertices, estimatedCounts, seed = edgeCounter, burnin, sampleSize)
	customSymmetricCounts[[edgeCounter]] <- customSymmetricResults$estimates
	currentCounts <- customSymmetricResults$estimates#c(currentCounts, tail(customResults$estimates, 1))
	if(any(is.na(currentCounts)) || any(currentCounts == 0))
	{
		save(customSymmetricCounts, file = filename)
		stop("Estimated a value of 0")
	}
	cat(edgeCounter, " / ", maxEdges, "\n", sep="")
	save(customSymmetricCounts, file = filename)
}

save(customSymmetricCounts, file = filename)
