library(Rmpfr)
load("./results27.RData")

library(tikzDevice)
tikz("./estimateCounts27.tex")
plot(log10(as.numeric(customSymmetricCounts[[length(customSymmetricCounts)]])), xlab = "k", ylab = "$A_{27, k}$")
dev.off()

pdf("./edgeDistribution27.pdf")
hist(customSymmetricResults$edgeCounts, xlab = "k")
dev.off()
