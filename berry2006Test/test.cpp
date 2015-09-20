#include "cliqueTree.h"
namespace chordalGraph
{
	void main()
	{
		cliqueTree::unionMinimalSeparatorsTemporaries temp;
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
#ifdef TRACK_GRAPH
		cliqueTree minimalTriangulation(7);
#endif
		bitsetType vertices;
		std::list<cliqueTree::externalEdge> edgeSequence;
		std::list<cliqueTree::cliqueTreeGraphType::vertex_descriptor> vertexSequence;
		std::vector<cliqueTree::externalEdge> addEdges;
		std::vector<cliqueTree::externalEdge> removeEdges;
		for (int i = 0; i < 7; i++)
		{
			maximalChordalSubgraph.addVertex();
			maximalChordalSubgraph.check();
#ifdef TRACK_GRAPH
			minimalTriangulation.addVertex();
			minimalTriangulation.check();
#endif
			cliqueTree::graphType::out_edge_iterator current, end;
			boost::tie(current, end) = boost::out_edges(i, baseGraph);
			while (current != end)
			{
				int targetVertex = (int)boost::target(*current, baseGraph);
				if (targetVertex < i)
				{
					maximalChordalSubgraph.unionMinimalSeparators(i, targetVertex, vertices, vertexSequence, edgeSequence, addEdges, removeEdges, temp);
					//vertices contains vertices on chordless paths between the two. 
					//We can only proceed if there are no vertices on the path (different components), 
					//or the set forms a complete graph / edges to those vertices already exist
					if (vertices.none() || vertices.count() == 1)
					{
						maximalChordalSubgraph.addEdge(targetVertex, i, vertices, vertexSequence, edgeSequence, addEdges, removeEdges, temp, true);
						maximalChordalSubgraph.check();
					}
					//Don't try and add an edge that's already been added into the triangulation
#ifdef TRACK_GRAPH
					if (!boost::edge(i, targetVertex, minimalTriangulation.getGraph()).second)
					{
						minimalTriangulation.unionMinimalSeparators(i, targetVertex, vertices, vertexSequence, edgeSequence, addEdges, removeEdges, temp);
						minimalTriangulation.addEdge(i, targetVertex, vertices, vertexSequence, edgeSequence, addEdges, removeEdges, temp, true);
					}
					minimalTriangulation.check();
#endif
				}
				current++;
			}
		}
	}
}
int main(int argc, char** argv)
{
	chordalGraph::main();
	return 0;
}