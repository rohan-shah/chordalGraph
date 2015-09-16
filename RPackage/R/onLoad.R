.onLoad <- function(libname, pkgname)
{
	library.dynam(package="chordalSubgraph", chname="chordalSubgraph", lib.loc = .libPaths())
}