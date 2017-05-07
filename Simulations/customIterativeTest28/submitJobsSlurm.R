for(vertices in 9:28)
{
	maxEdges <- vertices * (vertices - 1) / 2
	filename <- paste0("results", vertices, ".RData")
	shouldSubmit <- TRUE
	if(file.exists(filename))
	{
		load(filename)
		if(length(customSymmetricCounts) == maxEdges) shouldSubmit <- FALSE
	}
	if(shouldSubmit)
	{
		system2(command = "sbatch", args = c(paste0("--export=VERTICES=", vertices), "submitScriptSlurm.sh"), wait=TRUE)
	}
}
