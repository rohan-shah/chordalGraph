customMCMCSymmetric <- function(nVertices, approximateCounts, seed, burnIn, runSize)
{
	if(nVertices < 5)
	{
		stop("Input nVertices must be at least 5")
	}
	edgeLimit <- length(approximateCounts)-1
	if(edgeLimit > (nVertices * (nVertices + 1))/2 || edgeLimit <= 5)
	{
		stop("Input edgeLimit out of the valid range or smaller than 6")
	}
	if(isS4(approximateCounts))
	{
		approximateCountsCharacter <- as(approximateCounts, "character")
	}
	else
	{
		approximateCountsCharacter <- as.character(approximateCounts)
	}
	result <- .Call("customMCMCSymmetric", nVertices, approximateCountsCharacter, seed, burnIn, runSize, PACKAGE="chordalGraph")
	result$estimates <- mpfr(result$estimates)
	return(result)
}
