#include "cliqueTree.h"
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
namespace chordalGraph
{
#ifdef TRACK_GRAPH
	const cliqueTree::graphType& cliqueTree::getGraph() const
	{
		return graph;
	}
#endif
	const cliqueTree::cliqueTreeGraphType& cliqueTree::getCliqueGraph() const
	{
		return cliqueGraph;
	}
	cliqueTree::cliqueTree(int maximumVertices)
		:nVertices(0)
	{
		if (maximumVertices > MAX_STORAGE_VERTICES)
		{
			throw std::runtime_error("Requested number of vertices exceeded the amount of storage available");
		}
#ifdef TRACK_GRAPH
#ifdef USE_ADJACENCY_MATRIX_FOR_GRAPH
		graph = graphType(maximumVertices);
#endif
#endif
	}
	void cliqueTree::addVertex()
	{
		int previousVertexCount = nVertices;
		nVertices++;
#ifdef TRACK_GRAPH
		boost::add_vertex(graph);
#endif

		bitsetType newBitset(0);
		newBitset[previousVertexCount] = true;
		cliqueVertex newVertex(newBitset);
		int newCliqueVertexId = (int)boost::add_vertex(newVertex, cliqueGraph);

		componentIDs.push_back(previousVertexCount);
		verticesToCliqueVertices.push_back(newCliqueVertexId);
	}
	bool cliqueTree::tryAddVertexWithEdges(const bitsetType& involvedEdges, unionMinimalSeparatorsTemporaries& temp)
	{
		bitsetType copiedInvolvedEdges = involvedEdges;

		if (nVertices == 0)
		{
			cliqueVertex newVertex;
			newVertex.contents[0] = true;
			boost::add_vertex(newVertex, cliqueGraph);

#ifdef TRACK_GRAPH
			boost::add_vertex(graph);
#endif
			nVertices = 1;
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
		std::list<externalEdge> edgeSequence;
		std::vector<externalEdge> addEdges;
		std::vector<externalEdge> removeEdges;
		for (int i = 0; i < nVertices; i++)
		{
			if (copiedInvolvedEdges[i])
			{
				unionMinimalSeparatorBitset.reset();
				unionMinimalSeparators(nVertices, i, unionMinimalSeparatorBitset, vertexSequence, edgeSequence, addEdges, removeEdges, temp);
				//If we need to add something that we weren't going to, then we don't
				//have a chordal graph. So throw an error. 
				if ((unionMinimalSeparatorBitset & (~involvedEdges)).any())
				{
					return false;
				}
				addEdge(nVertices, i, unionMinimalSeparatorBitset, vertexSequence, edgeSequence, addEdges, removeEdges, temp, true);
				//We can now safely ignore the edges in unionMinimalSeparatorBitset
				//That is, they've been added so don't try and add them again
				copiedInvolvedEdges = copiedInvolvedEdges & (~unionMinimalSeparatorBitset);
				copiedInvolvedEdges[i] = false;
			}
		}
		return true;
	}
	void cliqueTree::addEdge(int vertexForExtraEdges, int v, bitsetType& unionMinimalSeparatorBitset, std::list<cliqueTreeGraphType::vertex_descriptor>& vertexSequence, std::list<externalEdge>& edgeSequence, std::vector<externalEdge>& addEdges, std::vector<externalEdge>& removeEdges, unionMinimalSeparatorsTemporaries& temp, bool hasPrecomputedUnionMinimalSeparator)
	{
		if (!hasPrecomputedUnionMinimalSeparator)
		{
			unionMinimalSeparatorBitset.reset();
			vertexSequence.clear();
			edgeSequence.clear();
			unionMinimalSeparators(vertexForExtraEdges, v, unionMinimalSeparatorBitset, vertexSequence, edgeSequence, addEdges, removeEdges, temp);
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
				boost::add_edge(verticesToCliqueVertices[vertexForExtraEdges], verticesToCliqueVertices[v], cliqueGraph);

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
				boost::add_edge(verticesToCliqueVertices[vertexForExtraEdges], verticesToCliqueVertices[v], cliqueGraph);

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
				boost::add_edge(verticesToCliqueVertices[v], newVertexId, cliqueGraph);
				boost::add_edge(verticesToCliqueVertices[vertexForExtraEdges], newVertexId, cliqueGraph);

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
			//Add every edge that is marked to be added.
			for (std::vector<externalEdge>::iterator i = addEdges.begin(); i != addEdges.end(); i++)
			{
				boost::add_edge(i->source, i->target, cliqueGraph);
			}
			//Remove every edge that has been marked for removal.
			for (std::vector<externalEdge>::iterator i = removeEdges.begin(); i != removeEdges.end(); i++)
			{
				boost::remove_edge(i->source, i->target, cliqueGraph);
			}
			std::vector<cliqueTreeGraphType::edge_descriptor> edgeSequenceInternalEdges;
			//Add vertex to every clique tree node and vertex on the path. 
			//Also make a copy of the edge_descriptor, now that we know that everything has been added.
			for (std::list<externalEdge>::iterator i = edgeSequence.begin(); i != edgeSequence.end(); i++)
			{
				cliqueTreeGraphType::edge_descriptor currentEdge = boost::edge(i->source, i->target, cliqueGraph).first;
				edgeSequenceInternalEdges.push_back(currentEdge);
				//After we've called addTo, we know we can treat it as a regular graph edge
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
				std::vector<cliqueTreeGraphType::edge_descriptor>::iterator currentPathEdge = edgeSequenceInternalEdges.begin();
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

						//This may invalidate references
						int newCliqueVertexId = (int)boost::add_vertex(newCliqueVertex, cliqueGraph);

						//Set up new edges between previous / next and new vertex
						cliqueTreeGraphType::edge_descriptor newCliqueEdgePreviousDescriptor = boost::add_edge(newCliqueVertexId, *previous, cliqueGraph).first;
						cliqueTreeGraphType::edge_descriptor newCliqueEdgeNextDescriptor = boost::add_edge(newCliqueVertexId, *next, cliqueGraph).first;
						//delete old edges
						boost::remove_edge(*current, *previous, cliqueGraph);
						boost::remove_edge(*current, *next, cliqueGraph);

						//Add edge between new path vertex and old path vertex
						boost::add_edge(newCliqueVertexId, *current, cliqueGraph);

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
				
				//make a copy of the minimal separators that contains both end-points
				bitsetType unionMinimalSeparatorBitsetCopy = unionMinimalSeparatorBitset;
				unionMinimalSeparatorBitsetCopy[vertexForExtraEdges] = unionMinimalSeparatorBitsetCopy[v] = true;
				if ((firstVertex.contents & ~unionMinimalSeparatorBitsetCopy).any())
				{
					cliqueVertex newVertex;
					newVertex.contents = firstVertex.contents & unionMinimalSeparatorBitsetCopy;
					newVertex.contents[vertexForExtraEdges] = true;
					firstVertex.contents[vertexForExtraEdges] = false;

					//This addition of a vertex may invalidate the reference firstVertex.
					int newVertexId = (int)boost::add_vertex(newVertex, cliqueGraph);

					//New edges
					boost::add_edge(firstVertexId, newVertexId, cliqueGraph);
					cliqueTreeGraphType::edge_descriptor newEdgeDescriptor = boost::add_edge(secondVertexId, newVertexId, cliqueGraph).first;
					
					//Update path
					*edgeSequenceInternalEdges.begin() = newEdgeDescriptor;
					*vertexSequence.begin() = newVertexId;

					//Remove existing edge
					boost::remove_edge(firstVertexId, secondVertexId, cliqueGraph);
					//Update verticesToCliqueVertices
					verticesToCliqueVertices[vertexForExtraEdges] = newVertexId;
				}
				//last
				int lastVertexId = (int)*vertexSequence.rbegin(), secondLastVertexId = (int)*std::next(vertexSequence.rbegin());
				cliqueVertex& lastVertex = boost::get(boost::vertex_name, cliqueGraph, lastVertexId);
				if ((lastVertex.contents & ~unionMinimalSeparatorBitsetCopy).any())
				{
					cliqueVertex newVertex;
					newVertex.contents = lastVertex.contents & unionMinimalSeparatorBitsetCopy;
					newVertex.contents[vertexForExtraEdges] = true;
					int newVertexId = (int)boost::add_vertex(newVertex, cliqueGraph);

					//New edges
					boost::add_edge(lastVertexId, newVertexId, cliqueGraph);
					cliqueTreeGraphType::edge_descriptor newEdgeDescriptor = boost::add_edge(secondLastVertexId, newVertexId, cliqueGraph).first;

					//Update path
					*edgeSequenceInternalEdges.begin() = newEdgeDescriptor;
					*vertexSequence.begin() = newVertexId;

					//Remove existing edge
					boost::remove_edge(lastVertexId, secondLastVertexId, cliqueGraph);
					//Update verticesToCliqueVertices
					verticesToCliqueVertices[vertexForExtraEdges] = newVertexId;
				}
			}
			//Now we check for non-maximal cliques
			//When we find one we contract that edge (and delete it from the collection of currentPath edges
			std::vector<cliqueTreeGraphType::edge_descriptor>::iterator currentPathEdge = edgeSequenceInternalEdges.begin();
			std::list<cliqueTreeGraphType::vertex_descriptor>::iterator currentPathVertex = vertexSequence.begin();
			while (currentPathEdge != edgeSequenceInternalEdges.end())
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
						if ((int)currentEdgeIterator->m_target != biggerVertexIndex)
						{
							boost::add_edge(currentEdgeIterator->m_target, biggerVertexIndex, cliqueGraph);
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
					currentPathEdge = edgeSequenceInternalEdges.erase(currentPathEdge);
					currentPathVertex = vertexSequence.erase(currentPathVertex);
				}
				else
				{
					currentPathVertex++;
					currentPathEdge++;
				}
			}
#ifdef TRACK_GRAPH
			//Add all the other required edges
			for (int i = 0; i < nVertices; i++)
			{
				if(unionMinimalSeparatorBitset[i] && !boost::edge(vertexForExtraEdges, i, graph).second) boost::add_edge(vertexForExtraEdges, i, graph);
			}
#endif
		}
#ifdef TRACK_GRAPH
		boost::add_edge(vertexForExtraEdges, v, graph);
#endif
	}
	void cliqueTree::unionMinimalSeparators(int u, int v, bitsetType& vertices, std::list<cliqueTreeGraphType::vertex_descriptor>& vertexSequence, std::list<externalEdge>& edgeSequence, std::vector<externalEdge>& addEdges, std::vector<externalEdge>& removeEdges, unionMinimalSeparatorsTemporaries& temp)
	{
		if (u == v)
		{
			throw std::runtime_error("Cannot have u == v in call to unionMinimalSeparators");
		}
		vertices.reset();
		addEdges.clear();
		removeEdges.clear();
		edgeSequence.clear();
		vertexSequence.clear();
		int nCliqueVertices = (int)boost::num_vertices(cliqueGraph);
		cliqueTreeGraphType::vertex_descriptor uVertex = verticesToCliqueVertices[u], vVertex = verticesToCliqueVertices[v];
		if(uVertex == vVertex) return;
		//Work out edges that are on the path from uVertex to vVertex, by doing a breadth first search. 
		//Added +1 to ensure we don't try and access an entry of a zero-length vector later on
		//replaced with temp.predecessorEdges
		//std::vector<cliqueTreeGraphType::edge_descriptor> predecessorEdges(nCliqueVertices+1);
		temp.predecessorEdges.resize(nCliqueVertices + 1);
		//This represents the vertices of the path
		//Added +1 to ensure we don't try and access an entry of a zero-length vector later on
		//replaced with temp.predecessorVertices
		//std::vector<cliqueTreeGraphType::vertex_descriptor> predecessorVertices(nCliqueVertices+1);
		temp.predecessorVertices.resize(nCliqueVertices + 1);
		
		temp.colorMap.resize(nCliqueVertices);
		typedef boost::color_traits<boost::default_color_type> Color;
		std::fill(temp.colorMap.begin(), temp.colorMap.end(), Color::white());
		boost::iterator_property_map<std::vector<boost::default_color_type>::iterator, boost::identity_property_map> colorMap(temp.colorMap.begin());
		//start from vertex u, go to vertex v
		boost::depth_first_visit(cliqueGraph, uVertex,
			boost::make_dfs_visitor(
					std::make_pair(
						boost::record_edge_predecessors(&temp.predecessorEdges[0], boost::on_tree_edge()),
						boost::record_predecessors(&temp.predecessorVertices[0], boost::on_tree_edge())
					)
				), colorMap
		);
		if (temp.colorMap[vVertex] == Color::black())
		{
			//Same connected component, so there is a minimal separator
			//Now we can get out the sequence of edges / vertices, starting at v and going through to u
			vertexSequence.push_back(vVertex);
			//Note there is one more element in vertexSequence than there is in edgeSequence
			cliqueTreeGraphType::vertex_descriptor currentVertex = vVertex;
			while (true)
			{
				cliqueTreeGraphType::edge_descriptor internalEdge = temp.predecessorEdges[currentVertex];
				externalEdge internalToExternal((int)boost::source(internalEdge, cliqueGraph), (int)boost::target(internalEdge, cliqueGraph));
				edgeSequence.push_back(internalToExternal);
				vertexSequence.push_back(temp.predecessorVertices[currentVertex]);
				if (temp.predecessorVertices[currentVertex] == uVertex)
				{
					break;
				}
				currentVertex = temp.predecessorVertices[currentVertex];
			}
			//We need to go from the start to the *last* element in the sequence that contains v. 
			while (true)
			{
				if (vertexSequence.size() > 1)
				{
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
				std::list<externalEdge>::iterator edgeSequenceIterator = edgeSequence.begin();
				std::list<cliqueTreeGraphType::vertex_descriptor>::iterator vertexSequenceIterator = std::next(vertexSequence.begin());
				for (; std::next(edgeSequenceIterator) != edgeSequence.end(); edgeSequenceIterator++,vertexSequenceIterator++)
				{
					externalEdge currentEdge = *edgeSequenceIterator, nextEdge = *std::next(edgeSequenceIterator);
					bitsetType currentEdgeSet = boost::get(boost::vertex_name, cliqueGraph, currentEdge.target).contents & boost::get(boost::vertex_name, cliqueGraph, currentEdge.source).contents;
					bitsetType nextEdgeSet = boost::get(boost::vertex_name, cliqueGraph, nextEdge.target).contents & boost::get(boost::vertex_name, cliqueGraph, nextEdge.source).contents;
					//this is the clique tree vertex shared between these edges
					cliqueTreeGraphType::vertex_descriptor sharedVertex = *vertexSequenceIterator;
					//currentEdgeSet is contained within nextEdgeSet
					if ((currentEdgeSet&(~nextEdgeSet)).none())
					{
						//Remove the edge associated with the smaller subset
						removeEdges.push_back(currentEdge);
						//Add an edge which skips over
						externalEdge edgeToAdd((int)*std::prev(vertexSequenceIterator), (int)*std::next(vertexSequenceIterator));
						addEdges.push_back(edgeToAdd);

						vertexSequence.erase(vertexSequenceIterator);
						edgeSequence.insert(edgeSequenceIterator, edgeToAdd);
						edgeSequence.erase(std::next(edgeSequenceIterator));
						edgeSequence.erase(edgeSequenceIterator);
						rearranged = true;
						break;
					}
					//nextEdgeSet is contained within currentEdgeSet
					else if ((nextEdgeSet&(~currentEdgeSet)).none())
					{
						//Remove the edge associated with the smaller subset
						removeEdges.push_back(nextEdge);
						//Add an edge which skips over
						externalEdge edgeToAdd((int)*std::prev(vertexSequenceIterator), (int)*std::next(vertexSequenceIterator));
						addEdges.push_back(edgeToAdd);

						vertexSequence.erase(vertexSequenceIterator);
						edgeSequence.insert(edgeSequenceIterator, edgeToAdd);
						edgeSequence.erase(std::next(edgeSequenceIterator));
						edgeSequence.erase(edgeSequenceIterator);
						rearranged = true;
						break;
					}
				}
				if (!rearranged) break;
			}
			//Now take the union of everything that's left in the edge sequence
			for (std::list<externalEdge>::iterator edgeSequenceIterator = edgeSequence.begin(); edgeSequenceIterator != edgeSequence.end(); edgeSequenceIterator++)
			{
				bitsetType relevantSubset = boost::get(boost::vertex_name, cliqueGraph, edgeSequenceIterator->source).contents & boost::get(boost::vertex_name, cliqueGraph, edgeSequenceIterator->target).contents;
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
	int cliqueTree::getNVertices() const
	{
		return nVertices;
	}
	void cliqueTree::check() const
	{
		std::size_t nCliqueVertices = boost::num_vertices(cliqueGraph);
		{
			cliqueTreeGraphType::vertex_iterator current, end;
			boost::tie(current, end) = boost::vertices(cliqueGraph);
			bitsetType unionCliqueVertices;
			for (; current != end; current++)
			{
				bitsetType contents = boost::get(boost::vertex_name, cliqueGraph, *current).contents;
				unionCliqueVertices |= contents;
#ifdef TRACK_GRAPH
				for (int i = 0; i < (int)nVertices; i++)
				{
					//Check that every clique is in fact a clique
					if (contents[i])
					{
						for (int j = 0; j < (int)nVertices; j++)
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
						for (int j = 0; j < (int)nVertices; j++)
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
#endif
			}
			//Every vertex appears in at least one clique vertex
			for (int i = 0; i < (int)nVertices; i++)
			{
				if (!unionCliqueVertices[i])
				{
					throw std::runtime_error("A vertex was not present in any clique vertex");
				}
			}
		}
		
		//Check that vertices induce connected subtrees.
		std::vector<int> connectedComponents(nCliqueVertices);
		for (int i = 0; i < (int)nVertices; i++)
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
	}
#ifdef HAS_NAUTY
	void cliqueTree::convertToNauty(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<::graph>& nautyGraph, std::vector<::graph>& cannonicalNautyGraph)
	{
		static DEFAULTOPTIONS_GRAPH(options);
		statsblk stats;
		options.getcanon = true;

		int n = nVertices;
		int m = SETWORDSNEEDED(n);
		nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
		lab.resize(n);
		ptn.resize(n);
		orbits.resize(n);
		nautyGraph.resize(n * m);
		EMPTYGRAPH(&(nautyGraph[0]), m, n);
		cliqueTree::cliqueTreeGraphType::vertex_iterator current, end;
		boost::tie(current, end) = boost::vertices(cliqueGraph);
		for(; current != end; current++)
		{
			bitsetType bitset = boost::get(boost::vertex_name, cliqueGraph, *current).contents;
			for(int i = 0; i < n; i++)
			{
				for(int j = 0; j < n; j++)
				{
					if(bitset[i] && bitset[j])
					{
						ADDONEEDGE(&(nautyGraph[0]), i, j, m);
					}
				}
			}
		}
		cannonicalNautyGraph.resize(n * m);
		densenauty(&(nautyGraph[0]), &(lab[0]), &(ptn[0]), &(orbits[0]), &options, &stats, m, n, &(cannonicalNautyGraph[0]));
	}
	void cliqueTree::convertToNautyWithEdge(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<::graph>& nautyGraph, std::vector<::graph>& cannonicalNautyGraph, int v1, int v2)
	{
		static DEFAULTOPTIONS_GRAPH(options);
		statsblk stats;
		options.getcanon = true;

		int n = nVertices;
		int m = SETWORDSNEEDED(n);
		nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
		lab.resize(n);
		ptn.resize(n);
		orbits.resize(n);
		nautyGraph.resize(n * m);
		EMPTYGRAPH(&(nautyGraph[0]), m, n);
		cliqueTree::cliqueTreeGraphType::vertex_iterator current, end;
		boost::tie(current, end) = boost::vertices(cliqueGraph);
		for(; current != end; current++)
		{
			bitsetType bitset = boost::get(boost::vertex_name, cliqueGraph, *current).contents;
			for(int i = 0; i < n; i++)
			{
				for(int j = 0; j < n; j++)
				{
					if(bitset[i] && bitset[j])
					{
						ADDONEEDGE(&(nautyGraph[0]), i, j, m);
					}
				}
			}
		}
		ADDONEEDGE(&(nautyGraph[0]), v1, v2, m);
		cannonicalNautyGraph.resize(n * m);
		densenauty(&(nautyGraph[0]), &(lab[0]), &(ptn[0]), &(orbits[0]), &options, &stats, m, n, &(cannonicalNautyGraph[0]));
	}/*	void cliqueTree::userlevelproc(int* lab, int* ptn, int level, int* orbits, statsblk* stats, int tv, int index, int tcellsize, int numcells, int childcount, int n)
	{
		*(mpz_class*)&(lab[-4]) *= index;
	}
	void cliqueTree::convertToNautyAndCountAutomorphisms(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<graph>& nautyGraph, std::vector<graph>& cannonicalNautyGraph, mpz_class& automorphismCount)
	{
#ifdef TRACK_GRAPH
		std::size_t nVertices = boost::num_vertices(graph);
#endif
		static DEFAULTOPTIONS_GRAPH(options);
		statsblk stats;
		options.getcanon = true;
		options.userlevelproc = &cliqueTree::userlevelproc;

		int n = nVertices;
		int m = SETWORDSNEEDED(n);
		nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
		lab.resize(n+4);
		ptn.resize(n);
		orbits.resize(n);
		nautyGraph.resize(n * m);
		EMPTYGRAPH(&(nautyGraph[0]), m, n);
		cliqueTree::cliqueTreeGraphType::vertex_iterator current, end;
		boost::tie(current, end) = boost::vertices(cliqueGraph);
		for(; current != end; current++)
		{
			bitsetType bitset = boost::get(boost::vertex_name, cliqueGraph, *current).contents;
			for(int i = 0; i < n; i++)
			{
				for(int j = 0; j < n; j++)
				{
					if(bitset[i] && bitset[j])
					{
						ADDONEEDGE(&(nautyGraph[0]), i, j, m);
					}
				}
			}
		}
		cannonicalNautyGraph.resize(n * m);
		automorphismCount = 1;
		*(mpz_class**)&(lab[0]) = &automorphismCount;
		densenauty(&(nautyGraph[0]), &(lab[4]), &(ptn[0]), &(orbits[0]), &options, &stats, m, n, &(cannonicalNautyGraph[0]));
	}*/
#endif
}
