library(Rmpfr)
library(ggplot2)
load("./results27.RData")

library(tikzDevice)
tikz("./estimateCounts27.tex", standAlone = TRUE)
plot(log10(as.numeric(customSymmetricCounts[[length(customSymmetricCounts)]])), xlab = "k", ylab = "$\\log_{10}\\left(A_{27, k}\\right)$", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

pdf("./edgeDistribution27.pdf")
data <- data.frame(value = customSymmetricResults$edgeCounts)
ggplot(data = data, mapping = aes(x = value)) + geom_histogram(bins = max(customSymmetricResults$edgeCounts)+1) + theme_bw() + xlab("k") + ylab("Count") + theme(axis.text.x = element_text(size = 18), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 20))
dev.off()


empty <- customSymmetricResults$edgeCounts == 0
complete <- customSymmetricResults$edgeCounts == max(customSymmetricResults$edgeCounts)
traverses <- 0
start <- head(which(empty | complete), 1)
current <- start
while(TRUE)
{
	browser()
	if(empty[current])
	{
		newCurrent <- head(which(complete[current:length(complete)]), 1) + (current - 1)
		if(length(newCurrent) == 0) break
		traverses <- traverses + 1
		current <- newCurrent
	}
	else if(complete[current])
	{
		newCurrent <- head(which(empty[current:length(complete)]), 1) + (current - 1)
		if(length(newCurrent) == 0) break
		traverses <- traverses + 1
		current <- newCurrent
	}
	else stop("Internal error")
}