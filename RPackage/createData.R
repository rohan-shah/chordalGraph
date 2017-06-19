library(chordalGraph)
for(i in 3:10)
{
	if(!file.exists(file.path("data", paste0("exact", i, ".csv"))))
	{
		val <- stochasticEnumerationNauty(nVertices = i, budget = 1000000, seed = 1)
		if(any(!(val@exact))) stop("Non exact computation")
		assign(x=paste0("exact", i), val="val")
		if(!file.exists("data")) dir.create("data")
		write.table(formatMpfr(val@data), file = file.path("data", paste0("exact", i, ".csv")), col.names=FALSE, row.names=FALSE)
	}
}
