#include "cliqueTreeAdjacencyMatrix.h"
#include <iostream>
namespace chordalGraph
{
	int main(int argc, char** argv)
	{
		std::vector<int> edges{0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 1, 4, 2, 4, 3, 4};
		cliqueTreeAdjacencyMatrix::graphType graph(5);
		for(int i = 0; i < (int)edges.size(); i+=2)
		{
			boost::add_edge(edges[i], edges[i + 1], graph);
		}
		cliqueTreeAdjacencyMatrix cliqueTree(5);
		bitsetType edgesBitset;
		cliqueTreeAdjacencyMatrix::unionMinimalSeparatorsTemporaries temp;
	
		//Vertex 0
		cliqueTree.addVertex();
	
		//Vertex 1
		edgesBitset.reset();
		edgesBitset[0] = true;
		cliqueTree.tryAddVertexWithEdges(edgesBitset, temp);

		//Vertex 2
		edgesBitset.reset();
		edgesBitset[0] = edgesBitset[1] = true;
		cliqueTree.tryAddVertexWithEdges(edgesBitset, temp);

		//Vertex 3
		edgesBitset.reset();
		edgesBitset[0] = edgesBitset[1] = true;
		cliqueTree.tryAddVertexWithEdges(edgesBitset, temp);

		//Vertex 4
		edgesBitset.reset();
		edgesBitset[0] = edgesBitset[1] = edgesBitset[2] = edgesBitset[3] = true;
		cliqueTree.tryAddVertexWithEdges(edgesBitset, temp);

		std::vector<int> outputCounts, counts1, counts2;
		std::vector<cliqueTreeAdjacencyMatrix::stackEntry> stack;
		std::vector<boost::default_color_type> colourVector;

		std::unordered_set<bitsetType> subsets;
		cliqueTreeAdjacencyMatrix copied(5);
		cliqueTreeAdjacencyMatrix::formRemovalTreeTemporaries removalTemporaries;
		cliqueTree.formRemovalTree(outputCounts, copied, 4, 3, subsets, removalTemporaries);
		return 0;
	}
}
int main(int argc, char** argv)
{
	return chordalGraph::main(argc, argv);
}
