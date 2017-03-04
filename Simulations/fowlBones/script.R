library(SIN)
library(chordalGraph)
data(fowlbones)
outerProductsSum <- fowlbones$corr * 276 * outer(fowlbones$stddev, fowlbones$stddev)
delta <- 5
exactCounts <- c(1, 15, 105, 455, 1320, 2526, 3085, 3255, 3000, 2235, 1206, 615, 260, 60, 15, 1)
burnin <- 1000
runSize <- 1000000

result <- .Call("customSymmetricPosteriorInference", outerProductsSum, delta, 6, 276, diag(6), as.character(exactCounts), burnin, runSize, PACKAGE = "chordalGraph")
sum <- matrix(0, 6, 6)
for(i in 1:length(result$graphs))
{
	sum <- sum + result$graphs[[i]] * result$probabilities[i]
}
