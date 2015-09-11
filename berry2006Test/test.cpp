#include "cliqueTree.h"
namespace chordalSubgraph
{
	void main()
	{
		cliqueTree::graphType baseGraph(7);
		boost::add_edge(1, 0, baseGraph);
		boost::add_edge(2, 1, baseGraph);
		boost::add_edge(3, 0, baseGraph);
		boost::add_edge(4, 2, baseGraph);
		boost::add_edge(4, 3, baseGraph);
		boost::add_edge(5, 2, baseGraph);
		boost::add_edge(5, 3, baseGraph);
		boost::add_edge(5, 4, baseGraph);
		boost::add_edge(6, 0, baseGraph);
		boost::add_edge(6, 1, baseGraph);

		cliqueTree maximalChordalSubgraph(7);
		cliqueTree minimalTriangulation(7);
		bitsetType vertices;
		std::list<cliqueTree::cliqueTreeGraphType::edge_descriptor> edgeSequence;
		std::list<cliqueTree::cliqueTreeGraphType::vertex_descriptor> vertexSequence;
		for (int i = 0; i < 7; i++)
		{
			maximalChordalSubgraph.addVertex();
			maximalChordalSubgraph.check();
			//minimalTriangulation.addVertex();
			//minimalTriangulation.check();
			cliqueTree::graphType::out_edge_iterator current, end;
			boost::tie(current, end) = boost::out_edges(i, baseGraph);
			while (current != end)
			{
				int targetVertex = (int)boost::target(*current, baseGraph);
				if (targetVertex < i)
				{
					maximalChordalSubgraph.unionMinimalSeparators(i, targetVertex, vertices, vertexSequence, edgeSequence);
					//vertices contains vertices on chordless paths between the two. 
					//We can only proceed if there are no vertices on the path (different components), 
					//or the set forms a complete graph / edges to those vertices already exist
					for (int j = 0; j < i; j++)
					{
						if (boost::edge(i, j, maximalChordalSubgraph.getGraph()).second)
						{
							vertices[j] = false;
						}
					}
					if (vertices.none())
					{
						maximalChordalSubgraph.addEdge(targetVertex, i, vertices, vertexSequence, edgeSequence, true);
						maximalChordalSubgraph.check();
					}
					//minimalTriangulation.unionMinimalSeparators(i, targetVertex, vertices, vertexSequence, edgeSequence);
					//minimalTriangulation.addEdge(i, targetVertex, vertices, vertexSequence, edgeSequence, true);
					//minimalTriangulation.check()
				}
				current++;
			}
		}
	}
}
int main(int argc, char** argv)
{
	chordalSubgraph::main();
	return 0;
}