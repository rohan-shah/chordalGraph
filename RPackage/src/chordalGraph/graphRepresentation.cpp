#include "graphRepresentation.h"
#include <stdexcept>
graphRepresentation toRepresentation(std::string str)
{
	if(str == "matrix")
	{
		return matrixRepresentation;
	}
	else if(str == "list")
	{
		return listRepresentation;
	}
	else
	{
		throw std::runtime_error("Input graphRepresentation must be either \"list\" or \"matrix\"");
	}
}
