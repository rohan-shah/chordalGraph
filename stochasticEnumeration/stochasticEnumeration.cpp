#include "cliqueTree.h"
#include "stochasticEnumeration.h"
namespace chordalGraph
{
	void main()
	{
		boost::mt19937 randomSource;
		randomSource.seed(1);
		stochasticEnumerationArgs args(randomSource);
		args.budget = 10000;
		args.nEdges = 6;
		args.nVertices = 6;
		stochasticEnumeration(args);
		std::cout << args.estimate.str() << std::endl;
	}
}
int main(int argc, char** argv)
{
	chordalGraph::main();
	return 0;
}