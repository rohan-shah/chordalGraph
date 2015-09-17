.onLoad <- function(libname, pkgname)
{
	library.dynam(package="chordalGraph", chname="chordalGraph", lib.loc = .libPaths())
}