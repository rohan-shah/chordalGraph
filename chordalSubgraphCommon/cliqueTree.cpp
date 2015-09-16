#include "cliqueTree.h"
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
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
		int newCliqueVertexId = (int)boost::add_vertex(newVertex, cliqueGraph);

		componentIDs.push_back(previousVertexCount);
		verticesToCliqueVertices.push_back(newCliqueVertexId);
	}
	bool cliqueTree::tryAddVertexWithEdges(const bitsetType& involvedEdges)
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
			return true;
		}
		addVertex();

		//If the new vertex isn't connected to anything else, take a short-cut
		if (copiedInvolvedEdges.none())	return true;

		//Now the general case
		bitsetType unionMinimalSeparatorBitset;
		std::list<cliqueTreeGraphType::vertex_descriptor> vertexSequence;
		std::list<cliqueTreeGraphType::edge_descriptor> edgeSequence;
		for (int i = 0; i < nVertices; i++)
		{
			if (copiedInvolvedEdges[i])
			{
				unionMinimalSeparatorBitset.reset();
				unionMinimalSeparators(nVertices, i, unionMinimalSeparatorBitset, vertexSequence, edgeSequence);
				//If we need to add something that we weren't going to, then we don't
				//have a chordal graph. So throw an error. 
				if ((unionMinimalSeparatorBitset & (~involvedEdges)).any())
				{
					return false;
				}
				addEdge(nVertices, i, unionMinimalSeparatorBitset, vertexSequence, edgeSequence, true);
				//We can now safely ignore the edges in unionMinimalSeparatorBitset
				//That is, they've been added so don't try and add them again
				copiedInvolvedEdges = copiedInvolvedEdges & (~unionMinimalSeparatorBitset);
				copiedInvolvedEdges[i] = false;
			}
		}
		return true;
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
				std::size_t nVertices = boost::num_vertices(graph);
				//If so we delete one clique vertex. This may invalidate references.
				boost::remove_vertex(verticesToCliqueVertices[v], cliqueGraph);
				//Compensate for removing the vertex by decreasing by 1 some of the entries
				//in verticesToCliqueVertices
				int deleted = verticesToCliqueVertices[v];
				for (int i = 0; i < (int)nVertices; i++)
				{
					if (verticesToCliqueVertices[i] >= deleted) verticesToCliqueVertices[i]--;
				}
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
				//This may invalidate references
				int newVertexId = (int)boost::add_vertex(newVertex, cliqueGraph);

				//Add extra edges
				cliqueEdge newEdge1, newEdge2;
				newEdge1.contents[v] = true;
				newEdge2.contents[vertexForExtraEdges] = true;
				boost::add_edge(verticesToCliqueVertices[v], newVertexId, newEdge1, cliqueGraph);
				boost::add_edge(verticesToCliqueVertices[vertexForExtraEdges], newVertexId, newEdge2, cliqueGraph);

				//Update componentIDs. Slightly more complicated this time. Note that std::replace takes inputs by reference.
				int oldValue = componentIDs[v], newValue = componentIDs[vertexForExtraEdges];
				std::replace(componentIDs.begin(), componentIDs.end(), oldValue, newValue);
				//No change to verticesToCliqueVertices
				//add corresponding edge to the graph is deferred to end of function
			}
		}
		//Now the case where the pair of vertices were already in the same connected component
		else
		{
			int nVertices = (int)boost::num_vertices(graph);
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
			{
				std::list<cliqueTreeGraphType::vertex_descriptor>::iterator current = std::next(vertexSequence.begin()), end = std::prev(vertexSequence.end());
				std::list<cliqueTreeGraphType::edge_descriptor>::iterator currentPathEdge = edgeSequence.begin();
				for (; current != end; current++, currentPathEdge++)
				{
					std::list<cliqueTreeGraphType::vertex_descriptor>::iterator previous = std::prev(current), next = std::next(current);
					cliqueVertex& currentCliqueVertex = boost::get(boost::vertex_name, cliqueGraph, *current);
					cliqueVertex& nextCliqueVertex = boost::get(boost::vertex_name, cliqueGraph, *next);
					cliqueVertex& previousCliqueVertex = boost::get(boost::vertex_name, cliqueGraph, *previous);
					bitsetType toSplitOut = currentCliqueVertex.contents & (~nextCliqueVertex.contents) & (~previousCliqueVertex.contents);
					//If the current clique vertex contains anything that's not contained in the previous or next vertices on the path, then we need to split something out.
					if (toSplitOut.any())
					{
						currentCliqueVertex.contents[vertexForExtraEdges] = false;

						cliqueVertex newCliqueVertex;
						//new clique vertex has the contents of the old clique vertex, restricted to those that are not contained in the neighbours
						//And it also contains vertexForExtraEdges. 
						newCliqueVertex.contents = currentCliqueVertex.contents;
						newCliqueVertex.contents &= ~toSplitOut;
						newCliqueVertex.contents[vertexForExtraEdges] = true;

						//Copy the contents of the vertices we have references to, as adding a vertex may invalidate the references.
						bitsetType nextCliqueVertexContents = nextCliqueVertex.contents;
						bitsetType previousCliqueVertexContents = previousCliqueVertex.contents;
						bitsetType currentCliqueVertexContents = currentCliqueVertex.contents;
						//This may invalidate references
						int newCliqueVertexId = (int)boost::add_vertex(newCliqueVertex, cliqueGraph);

						//Set up new edges between previous / next and new vertex
						cliqueEdge newCliqueEdgePrevious, newCliqueEdgeNext;
						newCliqueEdgePrevious.contents = newCliqueVertex.contents & previousCliqueVertexContents;
						newCliqueEdgeNext.contents = newCliqueVertex.contents & nextCliqueVertexContents;
						cliqueTreeGraphType::edge_descriptor newCliqueEdgePreviousDescriptor = boost::add_edge(newCliqueVertexId, *previous, newCliqueEdgePrevious, cliqueGraph).first;
						cliqueTreeGraphType::edge_descriptor newCliqueEdgeNextDescriptor = boost::add_edge(newCliqueVertexId, *next, newCliqueEdgeNext, cliqueGraph).first;
						//delete old edges
						boost::remove_edge(*current, *previous, cliqueGraph);
						boost::remove_edge(*current, *next, cliqueGraph);

						//Add edge between new path vertex and old path vertex
						cliqueEdge oldAndNewEdge;
						oldAndNewEdge.contents = newCliqueVertex.contents & currentCliqueVertexContents;
						boost::add_edge(newCliqueVertexId, *current, oldAndNewEdge, cliqueGraph);

						//Update path
						*current = newCliqueVertexId;
						*currentPathEdge = newCliqueEdgePreviousDescriptor;
						*std::next(currentPathEdge) = newCliqueEdgeNextDescriptor;

						//update the mapping of vertices to cliqueVertices
						verticesToCliqueVertices[vertexForExtraEdges] = newCliqueVertexId;
					}
				}
				//Special cases for the first and last vertices on the path. 

				//First
				int firstVertexId = (int)*vertexSequence.begin(), secondVertexId = (int)*std::next(vertexSequence.begin());
				cliqueVertex& firstVertex = boost::get(boost::vertex_name, cliqueGraph, firstVertexId);
				cliqueVertex& secondVertex = boost::get(boost::vertex_name, cliqueGraph, secondVertexId);
				
				//make a copy of the minimal separators that contains both end-points
				bitsetType unionMinimalSeparatorBitsetCopy = unionMinimalSeparatorBitset;
				unionMinimalSeparatorBitsetCopy[vertexForExtraEdges] = unionMinimalSeparatorBitsetCopy[v] = true;
				if ((firstVertex.contents & ~unionMinimalSeparatorBitsetCopy).any())
				{
					cliqueVertex newVertex;
					newVertex.contents = firstVertex.contents & unionMinimalSeparatorBitsetCopy;
					newVertex.contents[vertexForExtraEdges] = true;
					firstVertex.contents[vertexForExtraEdges] = false;

					//This addition of a vertex may invalidate the references firstVertex and secondVertex. So 
					//make a copy of their contents.
					bitsetType firstVertexContents = firstVertex.contents;
					bitsetType secondVertexContents = secondVertex.contents;
					int newVertexId = (int)boost::add_vertex(newVertex, cliqueGraph);

					//New edges
					cliqueEdge newEdge;
					newEdge.contents = firstVertexContents & newVertex.contents;
					boost::add_edge(firstVertexId, newVertexId, newEdge, cliqueGraph);
					newEdge.contents = secondVertexContents & newVertex.contents;
					cliqueTreeGraphType::edge_descriptor newEdgeDescriptor = boost::add_edge(secondVertexId, newVertexId, newEdge, cliqueGraph).first;
					
					//Update path
					*edgeSequence.begin() = newEdgeDescriptor;
					*vertexSequence.begin() = newVertexId;

					//Remove existing edge
					boost::remove_edge(firstVertexId, secondVertexId, cliqueGraph);
					//Update verticesToCliqueVertices
					verticesToCliqueVertices[vertexForExtraEdges] = newVertexId;
				}
				//last
				int lastVertexId = (int)*vertexSequence.rbegin(), secondLastVertexId = (int)*std::next(vertexSequence.rbegin());
				cliqueVertex& lastVertex = boost::get(boost::vertex_name, cliqueGraph, lastVertexId);
				cliqueVertex& secondLastVertex = boost::get(boost::vertex_name, cliqueGraph, secondLastVertexId);
				if ((lastVertex.contents & ~unionMinimalSeparatorBitsetCopy).any())
				{
					cliqueVertex newVertex;
					newVertex.contents = lastVertex.contents & unionMinimalSeparatorBitsetCopy;
					newVertex.contents[vertexForExtraEdges] = true;
					//lastVertex.contents[vertexForExtraEdges] = false;
					bitsetType lastVertexContents = lastVertex.contents;
					bitsetType secondLastVertexContents = secondLastVertex.contents;
					int newVertexId = (int)boost::add_vertex(newVertex, cliqueGraph);

					//New edges
					cliqueEdge newEdge;
					newEdge.contents = lastVertexContents & newVertex.contents;
					boost::add_edge(lastVertexId, newVertexId, newEdge, cliqueGraph);
					newEdge.contents = secondLastVertexContents & newVertex.contents;
					cliqueTreeGraphType::edge_descriptor newEdgeDescriptor = boost::add_edge(secondLastVertexId, newVertexId, newEdge, cliqueGraph).first;

					//Update path
					*edgeSequence.begin() = newEdgeDescriptor;
					*vertexSequence.begin() = newVertexId;

					//Remove existing edge
					boost::remove_edge(lastVertexId, secondLastVertexId, cliqueGraph);
					//Update verticesToCliqueVertices
					verticesToCliqueVertices[vertexForExtraEdges] = newVertexId;
				}
			}
			//Now we check for non-maximal cliques
			//When we find one we contract that edge (and delete it from the collection of currentPath edges
			std::list<cliqueTreeGraphType::edge_descriptor>::iterator currentPathEdge = edgeSequence.begin();
			std::list<cliqueTreeGraphType::vertex_descriptor>::iterator currentPathVertex = vertexSequence.begin();
			while (currentPathEdge != edgeSequence.end())
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
						if (currentEdgeIterator->m_target != biggerVertexIndex)
						{
							cliqueEdge newEdge;
							newEdge.contents = boost::get(boost::vertex_name, cliqueGraph, currentEdgeIterator->m_target).contents & boost::get(boost::vertex_name, cliqueGraph, biggerVertexIndex).contents;
							boost::add_edge(currentEdgeIterator->m_target, biggerVertexIndex, newEdge, cliqueGraph);
						}
						currentEdgeIterator++;
					}
					//Compensate for removing vertex and update verticesToCliqueVertices
					bitsetType deletedVertexContents = boost::get(boost::vertex_name, cliqueGraph, smallerVertexIndex).contents;
					int deleted = smallerVertexIndex;
					for (int i = 0; i < (int)nVertices; i++)
					{
						if (deletedVertexContents[i]) verticesToCliqueVertices[i] = biggerVertexIndex;
						if (verticesToCliqueVertices[i] >= deleted) verticesToCliqueVertices[i]--;
					}
					//Clear edges
					boost::clear_vertex(smallerVertexIndex, cliqueGraph);
					//remove vertex
					boost::remove_vertex(smallerVertexIndex, cliqueGraph);
					currentPathEdge = edgeSequence.erase(currentPathEdge);
					currentPathVertex = vertexSequence.erase(currentPathVertex);
				}
				else
				{
					currentPathVertex++;
					currentPathEdge++;
				}
			}
			//Add all the other required edges
			for (int i = 0; i < nVertices; i++)
			{
				if(unionMinimalSeparatorBitset[i] && !boost::edge(vertexForExtraEdges, i, graph).second) boost::add_edge(vertexForExtraEdges, i, graph);
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
		std::vector<cliqueTreeGraphType::edge_descriptor> predecessorEdges(nCliqueVertices+1);
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
			//We need to go from the start to the *last* element in the sequence that contains v. 
			while (true)
			{
				if (vertexSequence.size() > 1)
				{
					int firstCliqueVertex = (int)*vertexSequence.begin();
					int secondCliqueVertex = (int)*std::next(vertexSequence.begin());
					if (boost::get(boost::vertex_name, cliqueGraph, secondCliqueVertex).contents[v])
					{
						vertexSequence.pop_front();
						edgeSequence.pop_front();
						continue;
					}
				}
				break;
			}
			
			//Similarly from the end to the earliest element in the sequence that contains u. 
			while (true)
			{
				if (vertexSequence.size() > 1)
				{
					int firstCliqueVertex = (int)*vertexSequence.rbegin();
					int secondCliqueVertex = (int)*std::next(vertexSequence.rbegin());
					if (boost::get(boost::vertex_name, cliqueGraph, secondCliqueVertex).contents[u])
					{
						vertexSequence.pop_back();
						edgeSequence.pop_back();
						continue;
					}
				}
				break;
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
						boost::get(boost::edge_name, cliqueGraph, newlyAddedEdge).contents = currentEdgeSet & nextEdgeSet;

						vertexSequence.erase(vertexSequenceIterator);
						edgeSequence.insert(edgeSequenceIterator, newlyAddedEdge);
						edgeSequence.erase(std::next(edgeSequenceIterator));
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
						boost::get(boost::edge_name, cliqueGraph, newlyAddedEdge).contents = currentEdgeSet & nextEdgeSet;

						vertexSequence.erase(vertexSequenceIterator);
						edgeSequence.insert(edgeSequenceIterator, newlyAddedEdge);
						edgeSequence.erase(std::next(edgeSequenceIterator));
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
	//Used to filter clique tree by a graph vertex (the resulting graph should be connected)
	struct filterByVertex
	{
		filterByVertex(const cliqueTree::cliqueTreeGraphType* cliqueGraph)
			:cliqueGraph(cliqueGraph)
		{}
		filterByVertex()
			:cliqueGraph(NULL), vertex(-1)
		{}
		filterByVertex(const filterByVertex& other)
			:cliqueGraph(other.cliqueGraph), vertex(other.vertex)
		{}
		filterByVertex& operator=(const filterByVertex& other)
		{
			cliqueGraph = other.cliqueGraph;
			vertex = other.vertex;
			return *this;
		}
		const cliqueTree::cliqueTreeGraphType* cliqueGraph;
		int vertex;
		template <typename Vertex> bool operator()(const Vertex& v) const 
		{
			return boost::get(boost::vertex_name, *cliqueGraph, v).contents[vertex];
		}
	};
	void cliqueTree::check() const
	{
		std::size_t nVertices = boost::num_vertices(graph);
		std::size_t nCliqueVertices = boost::num_vertices(cliqueGraph);
		std::size_t nCliqueEdges = boost::num_edges(cliqueGraph);
		{
			cliqueTreeGraphType::vertex_iterator current, end;
			boost::tie(current, end) = boost::vertices(cliqueGraph);
			bitsetType unionCliqueVertices;
			for (; current != end; current++)
			{
				bitsetType contents = boost::get(boost::vertex_name, cliqueGraph, *current).contents;
				unionCliqueVertices |= contents;
				for (int i = 0; i < nVertices; i++)
				{
					//Check that every clique is in fact a clique
					if (contents[i])
					{
						for (int j = 0; j < nVertices; j++)
						{
							if (contents[j] & (i != j))
							{
								if (!boost::edge(i, j, graph).second)
								{
									throw std::runtime_error("Clique graph node was not actually a clique");
								}
							}
						}
					}
					//Check that every clique is maximal
					if (!contents[i])
					{
						bool foundMissingEdge = false;
						for (int j = 0; j < nVertices; j++)
						{
							if (contents[j])
							{
								if (!boost::edge(i, j, graph).second)
								{
									foundMissingEdge = true;
									break;
								}
							}
						}
						if (!foundMissingEdge)
						{
							throw std::runtime_error("Clique graph node was not maximal");
						}
					}
				}
			}
			//Every vertex appears in at least one clique vertex
			for (int i = 0; i < nVertices; i++)
			{
				if (!unionCliqueVertices[i])
				{
					throw std::runtime_error("A vertex was not present in any clique vertex");
				}
			}
		}
		
		//Check that vertices induce connected subtrees.
		std::vector<int> connectedComponents(nCliqueVertices);
		for (int i = 0; i < nVertices; i++)
		{
			filterByVertex filterObject(&cliqueGraph);
			filterObject.vertex = i;
			boost::filtered_graph<cliqueTree::cliqueTreeGraphType, boost::keep_all, filterByVertex> filteredGraph(cliqueGraph, boost::keep_all(), filterObject);
			int componentsCount = boost::connected_components(filteredGraph, &(connectedComponents[0]));
			if (componentsCount > 1)
			{
				throw std::runtime_error("Vertex induced a disconnected tree of the clique graph");
			}
		}

		//Sets for the edges are the intersections of the sets for the vertices
		{
			cliqueTreeGraphType::vertex_iterator firstVertex, secondVertex, end;
			boost::tie(firstVertex, end) = boost::vertices(cliqueGraph);
			for (; firstVertex != end; firstVertex++)
			{
				for (secondVertex = firstVertex + 1; secondVertex != end; secondVertex++)
				{
					bitsetType firstVertexSet = boost::get(boost::vertex_name, cliqueGraph, *firstVertex).contents;
					bitsetType secondVertexSet = boost::get(boost::vertex_name, cliqueGraph, *secondVertex).contents;
					bool hasEdge;
					cliqueTreeGraphType::edge_descriptor potentialEdge;
					boost::tie(potentialEdge, hasEdge) = boost::edge(*firstVertex, *secondVertex, cliqueGraph);
					if (hasEdge)
					{
						bitsetType edgeSet = boost::get(boost::edge_name, cliqueGraph, potentialEdge).contents;
						if (edgeSet != (firstVertexSet & secondVertexSet))
						{
							throw std::runtime_error("Found an edge with a vertex set not equal to the intersection of the vertex sets of its clique tree nodes");
						}
					}
				}
			}
		}
		//Intersection of two clique nodes is contained in the path between them.
		{
			cliqueTreeGraphType::vertex_iterator firstVertex, end;
			boost::tie(firstVertex, end) = boost::vertices(cliqueGraph);

			std::vector<cliqueTreeGraphType::edge_descriptor> predecessorEdges(nCliqueVertices + 1);
			std::vector<cliqueTreeGraphType::vertex_descriptor> predecessorVertices(nCliqueVertices + 1);

			typedef boost::color_traits<boost::default_color_type> Color;
			std::vector<boost::default_color_type> colorMap(nCliqueVertices, Color::white());
			for (; firstVertex != end; firstVertex++)
			{
				std::fill(colorMap.begin(), colorMap.end(), Color::white());
				//start from vertex u, go to vertex v
				boost::breadth_first_search(cliqueGraph, *firstVertex,
					boost::visitor(
						boost::make_bfs_visitor(
							std::make_pair(
								boost::record_edge_predecessors(&predecessorEdges[0], boost::on_tree_edge()),
								boost::record_predecessors(&predecessorVertices[0], boost::on_tree_edge())
							)
						)
					).color_map(&colorMap[0])
				);
				for (cliqueTreeGraphType::vertex_iterator secondVertex = firstVertex + 1; secondVertex != end; secondVertex++)
				{
					//There is only a path if they're in the same connected component
					if (colorMap[*secondVertex] == Color::black())
					{
						bitsetType firstVertexSet = boost::get(boost::vertex_name, cliqueGraph, *firstVertex).contents;
						bitsetType secondVertexSet = boost::get(boost::vertex_name, cliqueGraph, *secondVertex).contents;
						bitsetType intersectionSet = firstVertexSet & secondVertexSet;

						cliqueTreeGraphType::vertex_descriptor currentVertex = *secondVertex;
						while (true)
						{
							cliqueTreeGraphType::edge_descriptor existingEdge;
							bool hasExistingEdge;

							boost::tie(existingEdge, hasExistingEdge) = boost::edge(currentVertex, predecessorVertices[currentVertex], cliqueGraph);
							//If there is no path between them, break
							if (!hasExistingEdge) break;

							bitsetType existingEdgeSet = boost::get(boost::edge_name, cliqueGraph, existingEdge).contents;
							if ((existingEdgeSet & intersectionSet) != intersectionSet)
							{
								throw std::runtime_error("There was a pair of vertices so that the path between them did not contain the intersection");
							}
							currentVertex = predecessorVertices[currentVertex];
							if (currentVertex == *firstVertex)
							{
								break;
							}
						}
					}
				}
			}
		}
	}
}