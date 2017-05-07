library(Rmpfr)
load("./results28.RData")

pdf("./estimateCounts28.pdf")
plot(log10(as.numeric(customSymmetricCounts[[length(customSymmetricCounts)]])), xlab = "k", ylab = "A_{28, k}")
dev.off()

library(tikzDevice)
tikz("./estimateCounts28.tex")
plot(log10(as.numeric(customSymmetricCounts[[length(customSymmetricCounts)]])), xlab = "k", ylab = "A_{28, k}")
dev.off()

