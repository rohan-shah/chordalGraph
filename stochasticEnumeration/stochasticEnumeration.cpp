#include "cliqueTree.h"
#include "stochasticEnumeration.h"
#include <boost/program_options.hpp>
namespace chordalGraph
{
	int main(int argc, char** argv)
	{
		boost::program_options::options_description options("Usage");
		options.add_options()
			("seed", boost::program_options::value<int>(), "(int) The random seed to use. ")
			("budget", boost::program_options::value<int>(), "(int) The budget to use. ")
			("nVertices", boost::program_options::value<int>(), "(int) The number of vertices. ")
			("help", "Display this message");
		boost::program_options::variables_map variableMap;
		try
		{
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), variableMap);
		}
		catch (boost::program_options::error& ee)
		{

			std::cerr << "Error parsing command line arguments: " << ee.what() << std::endl << std::endl;
			std::cerr << options << std::endl;
			return -1;
		}
		if (variableMap.count("help") > 0)
		{
			std::cout << "Compute the number of chordal graphs with every number of edges, for the given number of vertices" << std::endl;
			std::cout << options << std::endl;
			return 0;
		}
		if (variableMap.count("seed") != 1 || variableMap.count("nVertices") != 1 || variableMap.count("budget") != 1)
		{
			std::cout << "Inputs seed, nVertices and budget are required" << std::endl;
			return 0;
		}
		int seed = variableMap["seed"].as<int>();
		int budget = variableMap["budget"].as<int>();
		int nVertices = variableMap["nVertices"].as<int>();

		if (budget < 1)
		{
			std::cout << "Input budget must be at least 1" << std::endl;
			return 0;
		}
		if (nVertices < 1)
		{
			std::cout << "Input nVertices must be at least 1" << std::endl;
		}
		boost::mt19937 randomSource;
		randomSource.seed(seed);
		stochasticEnumerationArgs args(randomSource);
		args.budget = budget;
		args.nVertices = nVertices;
		for (int nEdges = 0; nEdges < ((args.nVertices-1)*args.nVertices/2)+1; nEdges++)
		{
			args.nEdges = nEdges;
			stochasticEnumeration(args);
			std::cout << args.estimate.str() << std::endl;
		}
		return 0;
	}
}
int main(int argc, char** argv)
{
	return chordalGraph::main(argc, argv);
}