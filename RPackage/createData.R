library(chordalGraph)
for(i in 7:9)
{
	if(!file.exists(file.path("data", paste0("exact", i, ".RData"))))
	{
		assign(x=paste0("exact", i), val=stochasticEnumerationNauty(nVertices = i, budget = 1000000, seed = 1))
		if(!file.exists("data")) dir.create("data")
		save(list = paste0("exact", i), file = file.path("data", paste0("exact", i, ".RData")))
	}
}
#if(!file.exists("data", "exact9.RData"))
#{
	#exact9 <- stochasticEnumerationNauty(nVertices = 9, budget = 1000000, seed = 1)
	#if(!file.exists("data")) dir.create("data")
	#save(exact9, file = file.path("data", "exact9.RData"))
#}