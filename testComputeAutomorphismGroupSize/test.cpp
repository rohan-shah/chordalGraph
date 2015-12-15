#include "nauty.h"
#include <boost/program_options.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <iostream>
namespace chordalGraph
{
	std::size_t product = 1;
	void userlevelproc(int* lab, int* ptn, int level, int* orbits, statsblk* stats, int tv, int index, int tcellsize, int numcells, int childcount, int n)
	{
		product *= index;
	}
	void main(int argc, char** argv)
	{
		boost::program_options::options_description options("Usage");
		options.add_options()
			("graphFile", boost::program_options::value<std::string>(), "(path) The path to a graph file, in graphml format. ")
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
			return;
		}
		if (variableMap.count("help") > 0)
		{
			std::cout << "Compute the size of the automorphism group" << std::endl;
			std::cout << options << std::endl;
			return;
		}

		std::ifstream graphStream(variableMap["graphFile"].as<std::string>());
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graphType;
		graphType boostGraph;
		boost::dynamic_properties properties;
		try
		{
			boost::read_graphml(graphStream, boostGraph, properties, 0);
		}
		catch(std::exception& exp)
		{
			std::cout << "Error calling boost::read_graphml: " << exp.what() << std::endl;
			return;
		}
		std::size_t nVertices = boost::num_vertices(boostGraph);
		std::vector<int> lab;
		std::vector<int> ptn;
		std::vector<int> orbits;
		std::vector<graph> nautyGraph;
		std::vector<graph> cannonicalNautyGraph;

		static DEFAULTOPTIONS_GRAPH(nautyOptions);
		statsblk stats;
		//nautyOptions.getcanon = true;
		nautyOptions.userlevelproc = userlevelproc;

		int n = nVertices;
		int m = SETWORDSNEEDED(n);
		nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
		lab.resize(n);
		ptn.resize(n);
		orbits.resize(n);
		nautyGraph.resize(n * m);
		EMPTYGRAPH(&(nautyGraph[0]), m, n);
		graphType::edge_iterator current, end;
		boost::tie(current, end) = boost::edges(boostGraph);
		for(; current != end; current++)
		{
			ADDONEEDGE(&(nautyGraph[0]), (int)boost::source(*current, boostGraph), (int)boost::target(*current, boostGraph), m);
		}
		//cannonicalNautyGraph.resize(n * m);
		densenauty(&(nautyGraph[0]), &(lab[0]), &(ptn[0]), &(orbits[0]), &nautyOptions, &stats, m, n, &(cannonicalNautyGraph[0]));
		std::cout << product << std::endl;
	}
}
int main(int argc, char** argv)
{
	chordalGraph::main(argc, argv);
	return 0;
}
