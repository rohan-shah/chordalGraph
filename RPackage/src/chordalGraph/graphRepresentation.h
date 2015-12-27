#ifndef GRAPH_REPRESENTATION_HEADER_GUARD
#define GRAPH_REPRESENTATION_HEADER_GUARD
#include <string>
enum graphRepresentation
{
	listRepresentation, matrixRepresentation
};
graphRepresentation toRepresentation(std::string str);
#endif
