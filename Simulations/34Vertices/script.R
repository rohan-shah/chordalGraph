library(chordalGraph)
vertices <- 34
maxEdges <- vertices * (vertices - 1)/2
exactCounts <- mpfr(c(1, choose(maxEdges, 1), choose(maxEdges, 2), choose(maxEdges, 3), choose(maxEdges, 4) - choose(vertices, 4)*3, choose(maxEdges, 5) - choose(vertices, 5) * 12 - choose(vertices, 4)*3*(maxEdges - 6)), precBits = 64)
currentCounts <- exactCounts
sampleSize <- 200000
burnin <- 10000
maxEdges <- 150

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
	cat(edgeCounter, " / ", maxEdges, "\n", sep="")
}

#Run more iterations to stabilise bad count data.
approximateCounts <- customSymmetricCounts[[maxEdges]]
approximateCounts[1:6] <- exactCounts
for(i in 1:6)
{
	customSymmetricResults <- customMCMCSymmetric(vertices, approximateCounts = approximateCounts, seed = 1000+i, burnin, sampleSize)
	approximateCounts <- customSymmetricResults$estimates
	approximateCounts[1:6] <- exactCounts
}
approximateCounts <- customSymmetricResults$estimates
approximateCounts[1:6] <- exactCounts
customSymmetricResults <- customMCMCSymmetric(vertices, approximateCounts = approximateCounts, seed = 1000+i, burnin, sampleSize)

pdf("./customSymmetric34.pdf")
plot(1:sampleSize, customSymmetricResults$edgeCounts, type="l", main = "Sample path", xlab = "Step", ylab = "Edges")
dev.off()
pdf("./customSymmetricHist34.pdf")
hist(customSymmetricResults$edgeCounts, breaks = -1:maxEdges, main = "MCMC results", xlab = "Edges")
dev.off()

armstrongResults <- armstrongMCMC(vertices, approximateCounts = approximateCounts, seed = edgeCounter, burnin, sampleSize)
pdf("./armstrong34.pdf")
plot(1:sampleSize, armstrongResults$edgeCounts, type="l", main = "Sample path", xlab = "Step", ylab = "Edges")
dev.off()
pdf("./armstrongHist34.pdf")
hist(armstrongResults$edgeCounts, breaks = -1:maxEdges, main = "MCMC results", xlab = "Edges")
dev.off()

#There are numbers of edges for which only one graph is generated!!
table(armstrongResults$edgeCounts)






#Validate these results by looking at 20 edges - The results should be *highly* accurate, for both methods. 
approximateCounts <- customSymmetricCounts[[20]]
approximateCounts[1:6] <- exactCounts
for(i in 1:6)
{
	customSymmetricResults <- customMCMCSymmetric(vertices, approximateCounts = approximateCounts, seed = 2000+i, burnin, sampleSize)
	approximateCounts <- customSymmetricResults$estimates
	approximateCounts[1:6] <- exactCounts
}
approximateCounts <- customSymmetricResults$estimates
approximateCounts[1:6] <- exactCounts
customSymmetricResults <- customMCMCSymmetric(vertices, approximateCounts = approximateCounts, seed = 2000+i, burnin, sampleSize)

pdf("./customSymmetric34_20.pdf")
plot(1:sampleSize, customSymmetricResults$edgeCounts, type="l", main = "Sample path", xlab = "Step", ylab = "Edges")
dev.off()
pdf("./customSymmetricHist34_20.pdf")
hist(customSymmetricResults$edgeCounts, breaks = -1:20, main = "MCMC results", xlab = "Edges")
dev.off()

armstrongResults <- armstrongMCMC(vertices, approximateCounts = approximateCounts, seed = edgeCounter, burnin, sampleSize*10)
pdf("./armstrong34_20.pdf")
plot(1:(sampleSize*10), armstrongResults$edgeCounts, type="l", main = "Sample path", xlab = "Step", ylab = "Edges")
dev.off()
pdf("./armstrongHist34_20.pdf")
hist(armstrongResults$edgeCounts, breaks = -1:20, main = "MCMC results", xlab = "Edges")
dev.off()







save(customSymmetricCounts, file = "results.RData")
