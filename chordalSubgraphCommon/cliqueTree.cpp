#include "cliqueTree.h"
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
namespace chordalSubgraph
{
	const cliqueTree::graphType& cliqueTree::getGraph() const
	{
		return graph;
	}
	const cliqueTree::cliqueTreeGraphType& cliqueTree::getCliqueGraph() const
	{
		return cliqueGraph;
	}
	void cliqueTree::addVertex()
	{
		int previousVertexCount = (int)boost::num_vertices(graph);
		boost::add_vertex(graph);

		bitsetType newBitset(0);
		newBitset[previousVertexCount] = true;
		cliqueVertex newVertex(newBitset);
		int newCliqueVertexCount = (int)boost::add_vertex(newVertex, cliqueGraph);

		componentIDs.push_back(previousVertexCount);
		verticesToCliqueVertices.push_back(newCliqueVertexCount);
	}
	void cliqueTree::addVertexWithEdges(const bitsetType& involvedEdges)
	{
		bitsetType copiedInvolvedEdges = involvedEdges;

		int nVertices = (int)boost::num_vertices(graph);
		int nCliqueVertices = (int)boost::num_vertices(cliqueGraph);
		if (nVertices == 0)
		{
			cliqueVertex newVertex;
			newVertex.contents[0] = true;
			boost::add_vertex(newVertex, cliqueGraph);

			boost::add_vertex(graph);
			verticesToCliqueVertices.push_back(0);

			componentIDs.push_back(0);
			return;
		}
		//We begin by adding the vertex to the clique graph
		cliqueVertex newCliqueVertex;
		newCliqueVertex.contents[nVertices] = true;
		boost::add_vertex(newCliqueVertex, cliqueGraph);

		//and to the graph
		boost::add_vertex(graph);

		//Initially the added vertex is in its own connected component
		componentIDs.push_back(nVertices);

		//It's represented by the new clique vertex in the clique graph
		verticesToCliqueVertices.push_back(nCliqueVertices);

		//If the new vertex isn't connected to anything else, take a short-cut
		if (copiedInvolvedEdges.none())	return;

		//Now the general case

		//So now we have one more clique vertex and one more graph vertex
		nCliqueVertices++;
		nVertices++;

		bitsetType unionMinimalSeparatorBitset;
		std::list<cliqueTreeGraphType::vertex_descriptor> vertexSequence;
		std::list<cliqueTreeGraphType::edge_descriptor> edgeSequence;
		//This is nVertices-1 because we the new vertex can only have edges with the other
		//pre-existing vertices
		for (int i = 0; i < nVertices-1; i++)
		{
			if (copiedInvolvedEdges[i])
			{
				unionMinimalSeparatorBitset.reset();
				unionMinimalSeparators(nVertices, i, unionMinimalSeparatorBitset, vertexSequence, edgeSequence);
				//If we need to add something that we weren't going to, then we don't
				//have a chordal graph. So throw an error. 
				if ((unionMinimalSeparatorBitset & (~copiedInvolvedEdges)).any())
				{
					throw std::runtime_error("Trying to construct a clique tree for a non-chordal graph");
				}
				addEdge(nVertices, i, unionMinimalSeparatorBitset, vertexSequence, edgeSequence, true);
				//We can now safely ignore the edges in unionMinimalSeparatorBitset
				//That is, they've been added so don't try and add them again
				copiedInvolvedEdges = copiedInvolvedEdges & (~unionMinimalSeparatorBitset);
			}
		}
	}
	void cliqueTree::addEdge(int vertexForExtraEdges, int v, bitsetType& unionMinimalSeparatorBitset, std::list<cliqueTreeGraphType::vertex_descriptor>& vertexSequence, std::list<cliqueTreeGraphType::edge_descriptor>& edgeSequence, bool hasPrecomputedUnionMinimalSeparator)
	{
		if (!hasPrecomputedUnionMinimalSeparator)
		{
			unionMinimalSeparatorBitset.reset();
			vertexSequence.clear();
			edgeSequence.clear();
			unionMinimalSeparators(vertexForExtraEdges, v, unionMinimalSeparatorBitset, vertexSequence, edgeSequence);
		}
		cliqueVertex& extraEdgesCliqueVertex = boost::get(boost::vertex_name, cliqueGraph, verticesToCliqueVertices[vertexForExtraEdges]);
		cliqueVertex& vCliqueVertex = boost::get(boost::vertex_name, cliqueGraph, verticesToCliqueVertices[v]);

		int extraEdgesComponent = componentIDs[vertexForExtraEdges];
		int vComponent = componentIDs[v];
		//Note that the actual edge in the graph is added at the end of the function
		//First the case where they're in different connected components
		if (vComponent != extraEdgesComponent)
		{
			//Both the new vertex which was added, and the one its being joined with,
			//are contained in cliques of size 1 in the clique tree. 
			if (extraEdgesCliqueVertex.contents.count() == 1 && vCliqueVertex.contents.count() == 1)
			{
				//If so we delete one clique vertex
				boost::remove_vertex(verticesToCliqueVertices[v], cliqueGraph);
				//Get reference again, seeing as we removed a vertex. 
				cliqueVertex& extraEdgesCliqueVertex = boost::get(boost::vertex_name, cliqueGraph, verticesToCliqueVertices[vertexForExtraEdges]);
				//The other clique vertex now represents both vertices
				extraEdgesCliqueVertex.contents[v] = true;
				//update verticesToCliqueVertices
				verticesToCliqueVertices[v] = verticesToCliqueVertices[vertexForExtraEdges];
				//update componentIDs
				componentIDs[v] = componentIDs[vertexForExtraEdges];
				//add corresponding edge to the graph is deferred to end of function
			}
			else if (extraEdgesCliqueVertex.contents.count() == 1)
			{
				//Here the number of clique vertices is unchanged
				//Add v to the clique vertex of size 1
				extraEdgesCliqueVertex.contents[v] = true;

				//Add an extra edge having vertex v
				cliqueEdge newEdge;
				newEdge.contents[v] = true;
				boost::add_edge(verticesToCliqueVertices[vertexForExtraEdges], verticesToCliqueVertices[v], newEdge, cliqueGraph);

				//Update componentIDs. If we do it this way around it's one change
				componentIDs[vertexForExtraEdges] = componentIDs[v];
				//No change to verticesToCliqueVertices
				//add corresponding edge to the graph is deferred to end of function
			}
			//Essentially the same as the previous case
			else if (vCliqueVertex.contents.count() == 1)
			{
				//Add vertexForExtraEdges to the clique vertex of size 1
				vCliqueVertex.contents[vertexForExtraEdges] = true;

				//Add an extra edge having vertex v
				cliqueEdge newEdge;
				newEdge.contents[vertexForExtraEdges] = true;
				boost::add_edge(verticesToCliqueVertices[vertexForExtraEdges], verticesToCliqueVertices[v], newEdge, cliqueGraph);

				//Update componentIDs. If we do it this way around it's one change
				componentIDs[v] = componentIDs[vertexForExtraEdges];
				//No change to verticesToCliqueVertices
				//add corresponding edge to the graph is deferred to end of function
			}
			//General case, clique vertices representing both vertices have size > 1
			else
			{
				//Here we need a new clique vertex
				cliqueVertex newVertex;
				newVertex.contents[v] = newVertex.contents[vertexForExtraEdges] = true;
				int newVertexId = (int)boost::add_vertex(newVertex, cliqueGraph);

				//Add extra edges
				cliqueEdge newEdge1, newEdge2;
				newEdge1.contents[v] = true;
				newEdge2.contents[vertexForExtraEdges] = true;
				boost::add_edge(verticesToCliqueVertices[v], newVertexId, newEdge1, cliqueGraph);
				boost::add_edge(verticesToCliqueVertices[vertexForExtraEdges], newVertexId, newEdge2, cliqueGraph);

				//Update componentIDs. Slightly more complicated this time
				std::replace(componentIDs.begin(), componentIDs.end(), componentIDs[v], componentIDs[vertexForExtraEdges]);
				//No change to verticesToCliqueVertices
				//add corresponding edge to the graph is deferred to end of function
			}
		}
		//Now the case where the pair of vertices were already in the same connected component
		else
		{
			//Add vertex to every clique tree node and vertex on the path
			for (std::list<cliqueTreeGraphType::edge_descriptor>::iterator i = edgeSequence.begin(); i != edgeSequence.end(); i++)
			{
				boost::get(boost::edge_name, cliqueGraph, *i).contents[vertexForExtraEdges] = true;
			}
			for (std::list<cliqueTreeGraphType::vertex_descriptor>::iterator i = vertexSequence.begin(); i != vertexSequence.end(); i++)
			{
				boost::get(boost::vertex_name, cliqueGraph, *i).contents[vertexForExtraEdges] = true;
			}
			//Now we go through and check for vertices that lie in only one tree node on the path
			//This is simple, as every vertex induces a connected sub-tree. So we only need to look
			//at the vertices to the left and to the right. 
			std::list<cliqueTreeGraphType::vertex_descriptor>::iterator current = std::next(vertexSequence.begin()), end = std::prev(vertexSequence.end());
			int nVertices = (int)boost::num_vertices(graph);
			for (; current != end; current++)
			{
				std::list<cliqueTreeGraphType::vertex_descriptor>::iterator previous = std::prev(current), next = std::next(current);
				cliqueVertex& currentCliqueVertex = boost::get(boost::vertex_name, cliqueGraph, *current);
				cliqueVertex& nextCliqueVertex = boost::get(boost::vertex_name, cliqueGraph, *next);
				cliqueVertex& previousCliqueVertex = boost::get(boost::vertex_name, cliqueGraph, *previous);
				bitsetType toSplitOut = currentCliqueVertex.contents & (~nextCliqueVertex.contents) & (~previousCliqueVertex.contents);
				if (toSplitOut.any())
				{
					cliqueVertex newCliqueVertex;
					//new clique vertex has the contents of the old clique vertex (before we added the extra vertex in the loops above)
					newCliqueVertex.contents = currentCliqueVertex.contents;
					newCliqueVertex.contents[vertexForExtraEdges] = false;
					int newCliqueVertexId = (int)boost::add_vertex(newCliqueVertex, cliqueGraph);

					//Remove the ones that are going to be split out from the current clique vertex
					currentCliqueVertex.contents &= ~toSplitOut;
					cliqueEdge newCliqueEdge;
					newCliqueEdge.contents = newCliqueVertex.contents & currentCliqueVertex.contents;
					//update the mapping of vertices to cliqueVertices
					for (int i = 0; i < nVertices; i++)
					{
						if (toSplitOut[i])
						{
							verticesToCliqueVertices[i] = newCliqueVertexId;
						}
					}
				}
			}
			//Now we check for non-maximal cliques
			//When we find one we contract that edge (and delete it from the collection of currentPath edges
			std::list<cliqueTreeGraphType::edge_descriptor>::iterator currentPathEdge = edgeSequence.begin();
			std::list<cliqueTreeGraphType::vertex_descriptor>::iterator currentPathVertex = vertexSequence.begin();
			while(currentPathEdge != edgeSequence.end())
			{
				int firstVertexIndex = (int)boost::source(*currentPathEdge, cliqueGraph);
				int secondVertexIndex = (int)boost::target(*currentPathEdge, cliqueGraph);
				cliqueVertex& firstVertex = boost::get(boost::vertex_name, cliqueGraph, firstVertexIndex);
				cliqueVertex& secondVertex = boost::get(boost::vertex_name, cliqueGraph, secondVertexIndex);
				bool contract = false;
				int smallerVertexIndex, biggerVertexIndex;
				if ((firstVertex.contents & (~secondVertex.contents)).none())
				{
					contract = true;
					smallerVertexIndex = firstVertexIndex;
					biggerVertexIndex = secondVertexIndex;
				}
				else if ((secondVertex.contents & (~firstVertex.contents)).none())
				{
					contract = true;
					smallerVertexIndex = secondVertexIndex;
					biggerVertexIndex = firstVertexIndex;
				}
				if (contract)
				{
					cliqueTreeGraphType::out_edge_iterator currentEdgeIterator, endEdgeIterator;
					boost::tie(currentEdgeIterator, endEdgeIterator) = boost::out_edges(smallerVertexIndex, cliqueGraph);
					while (currentEdgeIterator != endEdgeIterator)
					{
						cliqueEdge newEdge;
						newEdge.contents = boost::get(boost::vertex_name, cliqueGraph, currentEdgeIterator->m_target).contents & boost::get(boost::vertex_name, cliqueGraph, biggerVertexIndex).contents;
						boost::add_edge(currentEdgeIterator->m_target, biggerVertexIndex, newEdge, cliqueGraph);
					}
					boost::remove_vertex(smallerVertexIndex, cliqueGraph);
				}
				else
				{
					currentPathVertex++;
					currentPathEdge++;
				}
			}
		}
		boost::add_edge(vertexForExtraEdges, v, graph);
	}
	void cliqueTree::unionMinimalSeparators(int u, int v, bitsetType& vertices, std::list<cliqueTreeGraphType::vertex_descriptor>& vertexSequence, std::list<cliqueTreeGraphType::edge_descriptor>& edgeSequence)
	{
		if (u == v)
		{
			throw std::runtime_error("Cannot have u == v in call to unionMinimalSeparators");
		}
		vertices.reset();
		int nCliqueVertices = (int)boost::num_vertices(cliqueGraph), nCliqueEdges = (int)boost::num_edges(cliqueGraph);
		cliqueTreeGraphType::vertex_descriptor uVertex, vVertex;
		cliqueTreeGraphType::vertex_iterator current, end;
		boost::tie(current, end) = boost::vertices(cliqueGraph);
		while (current != end)
		{
			if (boost::get(boost::vertex_name, cliqueGraph, *current).contents[u])
			{
				uVertex = *current;
			}
			if (boost::get(boost::vertex_name, cliqueGraph, *current).contents[v])
			{
				vVertex = *current;
			}
			current++;
		}
		//Work out edges that are on the path from uVertex to vVertex, by doing a breadth first search. 
		//Added +1 to ensure we don't try and access an entry of a zero-length vector later on
		std::vector<cliqueTreeGraphType::edge_descriptor> predecessorEdges(nCliqueEdges+1);
		//This represents the vertices of the path
		//Added +1 to ensure we don't try and access an entry of a zero-length vector later on
		std::vector<cliqueTreeGraphType::vertex_descriptor> predecessorVertices(nCliqueVertices+1);
		
		typedef boost::color_traits<boost::default_color_type> Color;
		std::vector<boost::default_color_type> colorMap(nCliqueVertices, Color::white());
		//start from vertex u, go to vertex v
		boost::breadth_first_search(cliqueGraph, uVertex,
			boost::visitor(
				boost::make_bfs_visitor(
					std::make_pair(
						boost::record_edge_predecessors(&predecessorEdges[0], boost::on_tree_edge()),
						boost::record_predecessors(&predecessorVertices[0], boost::on_tree_edge())
					)
				)
				).color_map(&colorMap[0])
		);
		edgeSequence.clear();
		vertexSequence.clear();
		if (colorMap[vVertex] == Color::black())
		{
			//Same connected component, so there is a minimal separator
			//Now we can get out the sequence of edges / vertices, starting at v and going through to u
			vertexSequence.push_back(vVertex);
			//Note there is one more element in vertexSequence than there is in edgeSequence
			cliqueTreeGraphType::vertex_descriptor currentVertex = vVertex;
			while (true)
			{
				edgeSequence.push_back(predecessorEdges[currentVertex]);
				vertexSequence.push_back(predecessorVertices[currentVertex]);
				if (predecessorVertices[currentVertex] == uVertex)
				{
					break;
				}
				currentVertex = predecessorVertices[currentVertex];
			}
			//Now start rearranging the clique tree if one of the vertex subsets used to represent one of the edges is a subset of the vertex subset used
			//to represent an adjacent edge
			while (true)
			{
				bool rearranged = false;
				std::list<cliqueTreeGraphType::edge_descriptor>::iterator edgeSequenceIterator = edgeSequence.begin();
				std::list<cliqueTreeGraphType::vertex_descriptor>::iterator vertexSequenceIterator = std::next(vertexSequence.begin());
				for (; std::next(edgeSequenceIterator) != edgeSequence.end(); edgeSequenceIterator++,vertexSequenceIterator++)
				{
					cliqueTreeGraphType::edge_descriptor currentEdge = *edgeSequenceIterator, nextEdge = *std::next(edgeSequenceIterator);
					bitsetType currentEdgeSet = boost::get(boost::edge_name, cliqueGraph, currentEdge).contents;
					bitsetType nextEdgeSet = boost::get(boost::edge_name, cliqueGraph, nextEdge).contents;
					//this is the clique tree vertex shared between these edges
					cliqueTreeGraphType::vertex_descriptor sharedVertex = *vertexSequenceIterator;
					//currentEdgeSet is contained within nextEdgeSet
					if ((currentEdgeSet&(~nextEdgeSet)).none())
					{
						//Remove the edge associated with the smaller subset
						boost::remove_edge(currentEdge, cliqueGraph);
						//Add an edge which skips over
						cliqueTreeGraphType::edge_descriptor newlyAddedEdge = boost::add_edge(*std::prev(vertexSequenceIterator), *std::next(vertexSequenceIterator), cliqueGraph).first;
						boost::get(boost::edge_name, cliqueGraph, newlyAddedEdge).contents = currentEdgeSet;

						vertexSequence.erase(vertexSequenceIterator);
						edgeSequence.erase(edgeSequenceIterator);
						rearranged = true;
						break;
					}
					//nextEdgeSet is contained within currentEdgeSet
					else if ((nextEdgeSet&(~currentEdgeSet)).none())
					{
						//Remove the edge associated with the smaller subset
						boost::remove_edge(nextEdge, cliqueGraph);
						//Add an edge which skips over
						cliqueTreeGraphType::edge_descriptor newlyAddedEdge = boost::add_edge(*std::prev(vertexSequenceIterator), *std::next(vertexSequenceIterator), cliqueGraph).first;
						boost::get(boost::edge_name, cliqueGraph, newlyAddedEdge).contents = currentEdgeSet;

						vertexSequence.erase(vertexSequenceIterator);
						edgeSequence.erase(edgeSequenceIterator);
						rearranged = true;
						break;
					}
				}
				if (!rearranged) break;
			}
			//Now take the union of everything that's left in the edge sequence
			for (std::list<cliqueTreeGraphType::edge_descriptor>::iterator edgeSequenceIterator = edgeSequence.begin(); edgeSequenceIterator != edgeSequence.end(); edgeSequenceIterator++)
			{
				bitsetType relevantSubset = boost::get(boost::edge_name, cliqueGraph, *edgeSequenceIterator).contents;
				vertices |= relevantSubset;
			}
			return;
		}
		else
		{
			//Different connected components, so vertices remains empty
			return;
		}
	}
	void cliqueTree::check() const
	{
		cliqueTreeGraphType::edge_iterator current, end;
		boost::tie(current, end) = boost::edges(cliqueGraph);
		for (; current != end; current++)
		{
			cliqueTreeGraphType::vertex_descriptor source = boost::source(*current, cliqueGraph);
			cliqueTreeGraphType::vertex_descriptor target = boost::target(*current, cliqueGraph);
			bitsetType sourceSet = boost::get(boost::vertex_name, cliqueGraph, source).contents;
			bitsetType targetSet = boost::get(boost::vertex_name, cliqueGraph, target).contents;
			bitsetType edgeSet = boost::get(boost::edge_name, cliqueGraph, *current).contents;
			if (edgeSet != (sourceSet & targetSet))
			{
				throw std::runtime_error("Found an edge with a vertex set not equal to the vertex sets of its clique tree nodes");
			}
		}
	}
}