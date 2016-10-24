#include "cliqueTreeAdjacencyMatrix.h"
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
namespace chordalGraph
{
#ifdef TRACK_GRAPH
	const cliqueTreeAdjacencyMatrix::graphType& cliqueTreeAdjacencyMatrix::getGraph() const
	{
		return graph;
	}
#endif
	struct removeEdgeHelper
	{
		removeEdgeHelper(int nVertices, int* componentIDs, int valueToInsert)
			:nVertices(nVertices), componentIDs(componentIDs), valueToInsert(valueToInsert)
		{}
		void initialize_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void finish_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void start_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void examine_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void tree_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void back_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void forward_or_cross_edge(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::edge_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{}
		void discover_vertex(cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_descriptor v, const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& g)
		{
			bitsetType contents = boost::get(boost::vertex_name, g, v).contents;
			for(int i = 0; i < nVertices; i++)
			{
				if(contents[i])
				{
					componentIDs[i] = valueToInsert;
				}
			}
		}
		int nVertices;
		int* componentIDs;
		int valueToInsert;
	};
	bool cliqueTreeAdjacencyMatrix::canRemoveEdge(int u, int v, std::vector<int>& counts, int& cliqueVertex)
	{
		cliqueVertex = -1;
		std::fill(counts.begin(), counts.end(), 0);
		for(int current = 0; current < nVertices; current++)
		{
			bitsetType contents = boost::get(boost::vertex_name, cliqueGraph, current).contents;
			if(contents[u])
			{
				if(contents[v])
				{
					cliqueVertex = current;
				}
				for(int other = 0; other < nVertices; other++)
				{
					if(contents[other]) counts[other]++;
				}
			}
		}
		if(u == v) return false;
		if(counts[v] != 1) return false;
		return true;
	}
	bool cliqueTreeAdjacencyMatrix::tryRemoveEdge(int u, int v, std::vector<boost::default_color_type>& colourVector, std::vector<int>& countsAfter)
	{
		if(u == v) return false;
		int cliqueVertex = -1;
		if(!canRemoveEdge(u, v, countsAfter, cliqueVertex)) return false;
		removeEdgeKnownCliqueVertex(u, v, colourVector, countsAfter, cliqueVertex);
		return true;
	}
	bool cliqueTreeAdjacencyMatrix::removeEdgeKnownCliqueVertex(int u, int v, std::vector<boost::default_color_type>& colourVector, std::vector<int>& counts2, int cliqueVertex)
	{
		bitsetType foundVertexContents = boost::get(boost::vertex_name, cliqueGraph, cliqueVertex).contents;
		bitsetType foundVertexContentsMinusU = foundVertexContents, foundVertexContentsMinusV = foundVertexContents;
		foundVertexContentsMinusU[u] = false;
		foundVertexContentsMinusV[v] = false;
	
		bool removeVertexMinusU = false, removeVertexMinusV = false;
		int largerVertexMinusU = -1, largerVertexMinusV = -1;
		for(int other = 0; other < nVertices; other++)
		{
			if(other == cliqueVertex) continue;
			std::pair<cliqueTreeGraphType::edge_descriptor, bool> p = boost::edge(cliqueVertex, other, cliqueGraph);
			if(p.second)
			{
				bitsetType targetContents = boost::get(boost::vertex_name, cliqueGraph, other).contents;
				if((~targetContents & foundVertexContentsMinusU).none())
				{
					largerVertexMinusU = other;
					removeVertexMinusU = true;
				}
				if((~targetContents & foundVertexContentsMinusV).none())
				{
					largerVertexMinusV = other;
					removeVertexMinusV = true;
				}
			}
		}
		//If there are only two vertices in this clique, then we're actually breaking a connected component
		if(foundVertexContents.count() == 2)
		{
			//Both U and V are connected to other stuff
			if(removeVertexMinusU && removeVertexMinusV)
			{
				remainingCliqueTreeVertices.push_back(cliqueVertex);
				boost::get(boost::vertex_name, cliqueGraph, cliqueVertex).contents.reset();
				verticesToCliqueVertices[v] = largerVertexMinusU;
				verticesToCliqueVertices[u] = largerVertexMinusV;

				for(int other = 0; other < nVertices; other++)
				{
					if(other == cliqueVertex) continue;
					std::pair<cliqueTreeGraphType::edge_descriptor, bool> p = boost::edge(cliqueVertex, other, cliqueGraph);
					if(p.second)
					{
						bitsetType targetContents = boost::get(boost::vertex_name, cliqueGraph, other).contents;
						if(targetContents[v] && other != largerVertexMinusU)
						{
							boost::add_edge(other, largerVertexMinusU, cliqueGraph);
							boost::remove_edge(other, cliqueVertex, cliqueGraph);
						}
						if(targetContents[u] && other != largerVertexMinusV)
						{
							boost::add_edge(other, largerVertexMinusV, cliqueGraph);
							boost::remove_edge(other, cliqueVertex, cliqueGraph);
						}
					}
				}
				boost::clear_vertex(cliqueVertex, cliqueGraph);

				int newComponentSize = *std::max_element(componentIDs.begin(), componentIDs.end())+1;
				removeEdgeHelper visitor(nVertices, &(componentIDs[0]), newComponentSize);
				colourVector.resize(nMaxVertices);
				boost::iterator_property_map<std::vector<boost::default_color_type>::iterator, boost::identity_property_map> colorMap(colourVector.begin());
				std::fill(colourVector.begin(), colourVector.end(), boost::default_color_type::white_color);
				boost::depth_first_visit(cliqueGraph, largerVertexMinusU, visitor, colorMap);

			}
			//Here only v is connected to other stuff, so u ends up by itself. 
			else if(removeVertexMinusU)
			{
				int newComponentSize = *std::max_element(componentIDs.begin(), componentIDs.end())+1;
				bitsetType& firstBitset = boost::get(boost::vertex_name, cliqueGraph, cliqueVertex).contents;
				firstBitset.reset();
				firstBitset[u] = true;
				componentIDs[u] = newComponentSize;

				verticesToCliqueVertices[u] = cliqueVertex;
				verticesToCliqueVertices[v] = largerVertexMinusU;
		
				for(int other = 0; other < nVertices; other++)
				{
					if(other == cliqueVertex) continue;
					std::pair<cliqueTreeGraphType::edge_descriptor, bool> p = boost::edge(cliqueVertex, other, cliqueGraph);
					if(p.second)
					{
						bitsetType targetContents = boost::get(boost::vertex_name, cliqueGraph, other).contents;
						if(targetContents[v] && other != largerVertexMinusU)
						{
							boost::add_edge(other, largerVertexMinusU, cliqueGraph);
							boost::remove_edge(other, cliqueVertex, cliqueGraph);
						}
					}
				}
				boost::clear_vertex(cliqueVertex, cliqueGraph);
			}
			//Here only u is connected to other stuff, so v ends up by itself. 
			else if(removeVertexMinusV)
			{
				int newComponentSize = *std::max_element(componentIDs.begin(), componentIDs.end())+1;
				bitsetType& firstBitset = boost::get(boost::vertex_name, cliqueGraph, cliqueVertex).contents;
				firstBitset.reset();
				firstBitset[v] = true;
				componentIDs[v] = newComponentSize;

				verticesToCliqueVertices[v] = cliqueVertex;
				verticesToCliqueVertices[u] = largerVertexMinusV;

				for(int other = 0; other < nVertices; other++)
				{
					if(other == cliqueVertex) continue;
					std::pair<cliqueTreeGraphType::edge_descriptor, bool> p = boost::edge(cliqueVertex, other, cliqueGraph);
					if(p.second)
					{
						bitsetType targetContents = boost::get(boost::vertex_name, cliqueGraph, other).contents;
						if(targetContents[u] && other != largerVertexMinusV)
						{
							boost::add_edge(other, largerVertexMinusV, cliqueGraph);
							boost::remove_edge(other, cliqueVertex, cliqueGraph);
						}
					}
				}
				boost::clear_vertex(cliqueVertex, cliqueGraph);
			}
			//In this case we're breaking a single isolated link.
			else
			{
				int newComponentSize = *std::max_element(componentIDs.begin(), componentIDs.end())+1;

				boost::clear_vertex(cliqueVertex, cliqueGraph);
				bitsetType& firstBitset = boost::get(boost::vertex_name, cliqueGraph, cliqueVertex).contents;
				firstBitset.reset();
				firstBitset[u] = true;
				componentIDs[u] = newComponentSize;
				verticesToCliqueVertices[u] = cliqueVertex;

				int newCliqueVertex = remainingCliqueTreeVertices.back();
				remainingCliqueTreeVertices.pop_back();
				bitsetType& secondBitset = boost::get(boost::vertex_name, cliqueGraph, newCliqueVertex).contents;
				secondBitset.reset();
				secondBitset[v] = true;
				componentIDs[v] = newComponentSize+1;
				verticesToCliqueVertices[v] = newCliqueVertex;
			}
		}
		else if(removeVertexMinusU && removeVertexMinusV)
		{
			verticesToCliqueVertices[v] = largerVertexMinusU;
			verticesToCliqueVertices[u] = largerVertexMinusV;
			for(int other = 0; other < nVertices; other++)
			{
				if(other == cliqueVertex) continue;
				std::pair<cliqueTreeGraphType::edge_descriptor, bool> p = boost::edge(cliqueVertex, other, cliqueGraph);
				if(p.second)
				{
					bitsetType targetContents = boost::get(boost::vertex_name, cliqueGraph, other).contents;
					if(targetContents[v] && other != largerVertexMinusU)
					{
						boost::add_edge(other, largerVertexMinusU, cliqueGraph);
						boost::remove_edge(other, cliqueVertex, cliqueGraph);
					}
					else if(targetContents[u] && other != largerVertexMinusV)
					{
						boost::add_edge(other, largerVertexMinusV, cliqueGraph);
						boost::remove_edge(other, cliqueVertex, cliqueGraph);
					}
					//Otherwise it doesn't actualy matter which clique it's connected to
					else if(other != largerVertexMinusV && other != largerVertexMinusU)
					{
						boost::add_edge(other, largerVertexMinusV, cliqueGraph);
						boost::remove_edge(other, cliqueVertex, cliqueGraph);
					}
				}
			}
			for(int other = 0; other < nVertices; other++)
			{
				if(verticesToCliqueVertices[other] == cliqueVertex) verticesToCliqueVertices[other] = largerVertexMinusU;
			}
			boost::add_edge(largerVertexMinusU, largerVertexMinusV, cliqueGraph);
			remainingCliqueTreeVertices.push_back(cliqueVertex);
			boost::clear_vertex(cliqueVertex, cliqueGraph);
			cliqueGraph.num_vertices--;
			boost::get(boost::vertex_name, cliqueGraph, cliqueVertex).contents.reset();
		}
		else if(removeVertexMinusU)
		{
			verticesToCliqueVertices[v] = largerVertexMinusU;
			verticesToCliqueVertices[u] = cliqueVertex;
			for(int other = 0; other < nVertices; other++)
			{
				if(other == cliqueVertex) continue;
				std::pair<cliqueTreeGraphType::edge_descriptor, bool> p = boost::edge(cliqueVertex, other, cliqueGraph);
				if(p.second)
				{
					bitsetType targetContents = boost::get(boost::vertex_name, cliqueGraph, other).contents;
					if(targetContents[v] && other != largerVertexMinusU)
					{
						boost::add_edge(other, largerVertexMinusU, cliqueGraph);
						boost::remove_edge(other, cliqueVertex, cliqueGraph);
					}
				}
			}
			boost::get(boost::vertex_name, cliqueGraph, cliqueVertex).contents[v] = false;
		}
		else if(removeVertexMinusV)
		{
			verticesToCliqueVertices[u] = largerVertexMinusV;
			verticesToCliqueVertices[v] = cliqueVertex;
			for(int other = 0; other < nVertices; other++)
			{
				if(other == cliqueVertex) continue;
				std::pair<cliqueTreeGraphType::edge_descriptor, bool> p = boost::edge(cliqueVertex, other, cliqueGraph);
				if(p.second)
				{
					bitsetType targetContents = boost::get(boost::vertex_name, cliqueGraph, other).contents;
					if(targetContents[u] && other != largerVertexMinusV)
					{
						boost::add_edge(other, largerVertexMinusV, cliqueGraph);
						boost::remove_edge(other, cliqueVertex, cliqueGraph);
					}
				}
			}
			boost::get(boost::vertex_name, cliqueGraph, cliqueVertex).contents[u] = false;
		}
		else
		{
			int newCliqueVertex = remainingCliqueTreeVertices.back();
			remainingCliqueTreeVertices.pop_back();
			cliqueGraph.num_vertices++;
			for(int other = 0; other < nVertices; other++)
			{
				if(other == cliqueVertex) continue;
				std::pair<cliqueTreeGraphType::edge_descriptor, bool> p = boost::edge(cliqueVertex, other, cliqueGraph);
				if(p.second)
				{
					bitsetType targetContents = boost::get(boost::vertex_name, cliqueGraph, other).contents;
					if(targetContents[v])
					{
						boost::add_edge(other, newCliqueVertex, cliqueGraph);
						boost::remove_edge(other, cliqueVertex, cliqueGraph);
					}
				}
			}
			boost::add_edge(newCliqueVertex, cliqueVertex, cliqueGraph);
			boost::get(boost::vertex_name, cliqueGraph, newCliqueVertex).contents = foundVertexContentsMinusU;
			boost::get(boost::vertex_name, cliqueGraph, cliqueVertex).contents = foundVertexContentsMinusV;
			verticesToCliqueVertices[u] = cliqueVertex;
			verticesToCliqueVertices[v] = newCliqueVertex;
		}
#ifdef TRACK_GRAPH
		boost::remove_edge(u, v, graph);
#endif
		std::fill(counts2.begin(), counts2.end(), 0);
		for(int other = 0; other < nVertices; other++)
		{
			bitsetType contents = boost::get(boost::vertex_name, cliqueGraph, other).contents;
			if(contents[u])
			{
				for(int other2 = 0; other2 < nVertices; other2++)
				{
					if(contents[other2])
					{
						counts2[other2]++;
					}
				}
			}
		}
#ifndef NDEBUG
		check();
#endif
		return true;
	}
	void cliqueTreeAdjacencyMatrix::makeCopy(const cliqueTreeAdjacencyMatrix& other)
	{
		cliqueGraph.makeCopy(other.cliqueGraph); 
#ifdef TRACK_GRAPH
		graph.makeCopy(other.graph);
#endif
		nVertices = other.nVertices;
		nMaxVertices = other.nMaxVertices;
		verticesToCliqueVertices = other.verticesToCliqueVertices;
		componentIDs = other.componentIDs;
		remainingCliqueTreeVertices = other.remainingCliqueTreeVertices;
	}
	const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType& cliqueTreeAdjacencyMatrix::getCliqueGraph() const
	{
		return cliqueGraph;
	}
	cliqueTreeAdjacencyMatrix::cliqueTreeAdjacencyMatrix(int maximumVertices)
		:nVertices(0), nMaxVertices(maximumVertices)
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
		cliqueGraph = cliqueTreeGraphType(maximumVertices);
		remainingCliqueTreeVertices.resize(maximumVertices);
		for(int i = 0; i < maximumVertices; i++) remainingCliqueTreeVertices[maximumVertices - 1 - i] = i;
	}
	
	void cliqueTreeAdjacencyMatrix::addVertex()
	{
		int previousVertexCount = nVertices;
		nVertices++;
#ifdef TRACK_GRAPH
#ifndef USE_ADJACENCY_MATRIX_FOR_GRAPH
		boost::add_vertex(graph);
#else
		graph.num_vertices++;
#endif
#endif

		bitsetType newBitset(0);
		newBitset[previousVertexCount] = true;
		int newCliqueVertexId = *remainingCliqueTreeVertices.rbegin();
		remainingCliqueTreeVertices.pop_back();
		cliqueGraph.num_vertices++;
		cliqueVertex& newVertex = boost::get(boost::vertex_name, cliqueGraph, newCliqueVertexId);
		newVertex.contents = newBitset;

		componentIDs.push_back(previousVertexCount);
		verticesToCliqueVertices.push_back(newCliqueVertexId);
	}
	bool cliqueTreeAdjacencyMatrix::tryAddVertexWithEdges(const bitsetType& involvedEdges, unionMinimalSeparatorsTemporaries& temp)
	{
		bitsetType copiedInvolvedEdges = involvedEdges;

		if (nVertices == 0)
		{
			int newCliqueVertexId = *remainingCliqueTreeVertices.rbegin();
			remainingCliqueTreeVertices.pop_back();
			cliqueGraph.num_vertices++;
			cliqueVertex& newVertex = boost::get(boost::vertex_name, cliqueGraph, newCliqueVertexId);
			newVertex.contents = 0;
			newVertex.contents[0] = true;
#ifdef TRACK_GRAPH
#ifndef USE_ADJACENCY_MATRIX_FOR_GRAPH
			boost::add_vertex(graph);
#endif
#endif
			nVertices = 1;
			verticesToCliqueVertices.push_back(newCliqueVertexId);

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
		for (int i = 0; i < nVertices-1; i++)
		{
			if (copiedInvolvedEdges[i])
			{
				unionMinimalSeparatorBitset.reset();
				unionMinimalSeparators(nVertices-1, i, unionMinimalSeparatorBitset, vertexSequence, edgeSequence, addEdges, removeEdges, temp);
				//If we need to add something that we weren't going to, then we don't
				//have a chordal graph. So throw an error. 
				if ((unionMinimalSeparatorBitset & (~involvedEdges)).any())
				{
					return false;
				}
				addEdge(nVertices-1, i, unionMinimalSeparatorBitset, vertexSequence, edgeSequence, addEdges, removeEdges, temp, true);
				//We can now safely ignore the edges in unionMinimalSeparatorBitset
				//That is, they've been added so don't try and add them again
				copiedInvolvedEdges = copiedInvolvedEdges & (~unionMinimalSeparatorBitset);
				copiedInvolvedEdges[i] = false;
			}
		}
		return true;
	}
	void cliqueTreeAdjacencyMatrix::addEdge(int vertexForExtraEdges, int v, bitsetType& unionMinimalSeparatorBitset, std::list<cliqueTreeGraphType::vertex_descriptor>& vertexSequence, std::list<externalEdge>& edgeSequence, std::vector<externalEdge>& addEdges, std::vector<externalEdge>& removeEdges, unionMinimalSeparatorsTemporaries& temp, bool hasPrecomputedUnionMinimalSeparator)
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
				int deleted = verticesToCliqueVertices[v];
				boost::clear_vertex(deleted, cliqueGraph);
				boost::get(boost::vertex_name, cliqueGraph, deleted).contents.reset();
				cliqueGraph.num_vertices--;
				remainingCliqueTreeVertices.push_back(deleted);
				
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
				int newVertexId = *remainingCliqueTreeVertices.rbegin();
				remainingCliqueTreeVertices.pop_back();
				cliqueGraph.num_vertices++;
				cliqueVertex& newVertex = boost::get(boost::vertex_name, cliqueGraph, newVertexId);
				newVertex.contents = 0;
				newVertex.contents[v] = newVertex.contents[vertexForExtraEdges] = true;

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

						int newCliqueVertexId = *remainingCliqueTreeVertices.rbegin();
						remainingCliqueTreeVertices.pop_back();
						cliqueGraph.num_vertices++;
						cliqueVertex& newCliqueVertex = boost::get(boost::vertex_name, cliqueGraph, newCliqueVertexId);
						//new clique vertex has the contents of the old clique vertex, restricted to those that are not contained in the neighbours
						//And it also contains vertexForExtraEdges. 
						newCliqueVertex.contents = currentCliqueVertex.contents;
						newCliqueVertex.contents &= ~toSplitOut;
						newCliqueVertex.contents[vertexForExtraEdges] = true;

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
					int newCliqueVertexId = *remainingCliqueTreeVertices.rbegin();
					remainingCliqueTreeVertices.pop_back();
					cliqueGraph.num_vertices++;
					cliqueVertex& newVertex = boost::get(boost::vertex_name, cliqueGraph, newCliqueVertexId);
					newVertex.contents = firstVertex.contents & unionMinimalSeparatorBitsetCopy;
					newVertex.contents[vertexForExtraEdges] = true;
					firstVertex.contents[vertexForExtraEdges] = false;

					//New edges
					boost::add_edge(firstVertexId, newCliqueVertexId, cliqueGraph);
					cliqueTreeGraphType::edge_descriptor newEdgeDescriptor = boost::add_edge(secondVertexId, newCliqueVertexId, cliqueGraph).first;
					
					//Update path
					*edgeSequenceInternalEdges.begin() = newEdgeDescriptor;
					*vertexSequence.begin() = newCliqueVertexId;

					//Remove existing edge
					boost::remove_edge(firstVertexId, secondVertexId, cliqueGraph);
					//Update verticesToCliqueVertices
					verticesToCliqueVertices[vertexForExtraEdges] = newCliqueVertexId;
				}
				//last
				int lastVertexId = (int)*vertexSequence.rbegin(), secondLastVertexId = (int)*std::next(vertexSequence.rbegin());
				cliqueVertex& lastVertex = boost::get(boost::vertex_name, cliqueGraph, lastVertexId);
				if ((lastVertex.contents & ~unionMinimalSeparatorBitsetCopy).any())
				{
					int newCliqueVertexId = *remainingCliqueTreeVertices.rbegin();
					remainingCliqueTreeVertices.pop_back();
					cliqueGraph.num_vertices++;
					cliqueVertex& newVertex = boost::get(boost::vertex_name, cliqueGraph, newCliqueVertexId);
					newVertex.contents = lastVertex.contents & unionMinimalSeparatorBitsetCopy;
					newVertex.contents[vertexForExtraEdges] = true;

					//New edges
					boost::add_edge(lastVertexId, newCliqueVertexId, cliqueGraph);
					cliqueTreeGraphType::edge_descriptor newEdgeDescriptor = boost::add_edge(secondLastVertexId, newCliqueVertexId, cliqueGraph).first;

					//Update path
					*edgeSequenceInternalEdges.begin() = newEdgeDescriptor;
					*vertexSequence.begin() = newCliqueVertexId;

					//Remove existing edge
					boost::remove_edge(lastVertexId, secondLastVertexId, cliqueGraph);
					//Update verticesToCliqueVertices
					verticesToCliqueVertices[vertexForExtraEdges] = newCliqueVertexId;
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
					for (int i = 0; i < (int)nVertices; i++)
					{
						if (deletedVertexContents[i]) verticesToCliqueVertices[i] = biggerVertexIndex;
					}
					//Clear edges
					boost::clear_vertex(smallerVertexIndex, cliqueGraph);
					boost::get(boost::vertex_name, cliqueGraph, smallerVertexIndex).contents.reset();
					cliqueGraph.num_vertices--;
					remainingCliqueTreeVertices.push_back(smallerVertexIndex);
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
	int cliqueTreeAdjacencyMatrix::maxVertices() const
	{
		return nMaxVertices;
	}
	void cliqueTreeAdjacencyMatrix::unionMinimalSeparators(int u, int v, bitsetType& vertices, std::list<cliqueTreeGraphType::vertex_descriptor>& vertexSequence, std::list<externalEdge>& edgeSequence, std::vector<externalEdge>& addEdges, std::vector<externalEdge>& removeEdges, unionMinimalSeparatorsTemporaries& temp)
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
		cliqueTreeGraphType::vertex_descriptor uVertex = verticesToCliqueVertices[u], vVertex = verticesToCliqueVertices[v];
		if(uVertex == vVertex) return;
		//Work out edges that are on the path from uVertex to vVertex, by doing a breadth first search. 
		//Added +1 to ensure we don't try and access an entry of a zero-length vector later on
		//replaced with temp.predecessorEdges
		//std::vector<cliqueTreeGraphType::edge_descriptor> predecessorEdges(nCliqueVertices+1);
		temp.predecessorEdges.resize(maxVertices());
		//This represents the vertices of the path
		//Added +1 to ensure we don't try and access an entry of a zero-length vector later on
		//replaced with temp.predecessorVertices
		//std::vector<cliqueTreeGraphType::vertex_descriptor> predecessorVertices(nCliqueVertices+1);
		temp.predecessorVertices.resize(maxVertices());
		
		temp.colorMap.resize(maxVertices());
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
					bitsetType currentEdgeSet = boost::get(boost::vertex_name, cliqueGraph, currentEdge.source).contents & boost::get(boost::vertex_name, cliqueGraph, currentEdge.target).contents;
					bitsetType nextEdgeSet = boost::get(boost::vertex_name, cliqueGraph, nextEdge.source).contents & boost::get(boost::vertex_name, cliqueGraph, nextEdge.target).contents;
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
	struct filterByVertexAdjacency
	{
		filterByVertexAdjacency(const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType* cliqueGraph)
			:cliqueGraph(cliqueGraph)
		{}
		filterByVertexAdjacency()
			:cliqueGraph(NULL), vertex(-1)
		{}
		filterByVertexAdjacency(const filterByVertexAdjacency& other)
			:cliqueGraph(other.cliqueGraph), vertex(other.vertex)
		{}
		filterByVertexAdjacency& operator=(const filterByVertexAdjacency& other)
		{
			cliqueGraph = other.cliqueGraph;
			vertex = other.vertex;
			return *this;
		}
		const cliqueTreeAdjacencyMatrix::cliqueTreeGraphType* cliqueGraph;
		int vertex;
		template <typename Vertex> bool operator()(const Vertex& v) const 
		{
			return boost::get(boost::vertex_name, *cliqueGraph, v).contents[vertex];
		}
	};
	int cliqueTreeAdjacencyMatrix::getNVertices() const
	{
		return nVertices;
	}
	void cliqueTreeAdjacencyMatrix::check() const
	{
		{
			cliqueTreeGraphType::vertex_iterator current, end;
			boost::tie(current, end) = boost::vertices(cliqueGraph);
			bitsetType unionCliqueVertices;
			for (; current != end; current++)
			{
				bitsetType contents = boost::get(boost::vertex_name, cliqueGraph, *current).contents;
				if(contents.none()) continue;
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
		{
			cliqueTreeGraphType::vertex_iterator current, end;
			boost::tie(current, end) = boost::vertices(cliqueGraph);
			for (; current != end; current++)
			{
				bitsetType contents = boost::get(boost::vertex_name, cliqueGraph, *current).contents;
				cliqueTreeGraphType::vertex_iterator current2, end2;
				boost::tie(current2, end2) = boost::vertices(cliqueGraph);
				for (; current2 != end2; current2++)
				{
					std::pair<cliqueTreeGraphType::edge_descriptor, bool> p = boost::edge(*current, *current2, cliqueGraph);
					bitsetType contents2 = boost::get(boost::vertex_name, cliqueGraph, *current2).contents;
					if(p.second)
					{
						if(current == current2) throw std::runtime_error("Self edges not allowed");
						if((contents & contents2).none())
						{
							throw std::runtime_error("Edge between clique tree vertices without common contents");
						}

					}
				}
			}
		}
		//Check that vertices induce connected subtrees.
		std::vector<int> connectedComponents(nMaxVertices);
		std::vector<boost::default_color_type> colours(nMaxVertices);
		for (int i = 0; i < (int)nVertices; i++)
		{
			filterByVertexAdjacency filterObject(&cliqueGraph);
			filterObject.vertex = i;
			boost::filtered_graph<cliqueTreeAdjacencyMatrix::cliqueTreeGraphType, boost::keep_all, filterByVertexAdjacency> filteredGraph(cliqueGraph, boost::keep_all(), filterObject);
			int componentsCount = boost::connected_components(filteredGraph, &(connectedComponents[0]), boost::color_map(&(colours[0])));
			if (componentsCount > 1)
			{
				throw std::runtime_error("Vertex induced a disconnected tree of the clique graph");
			}
		}
		//Check that verticesToCliqueVertices structure is correct
		for (int i = 0; i < (int)nVertices; i++)
		{
			if (!boost::get(boost::vertex_name, cliqueGraph, verticesToCliqueVertices[i]).contents[i])
			{
				throw std::runtime_error("Member verticesToCliqueVertices contained incorrect data");
			}
		}
#ifdef TRACK_GRAPH
		//Check that the componentIDs entry is correct. 
		boost::connected_components(graph, &(connectedComponents[0]));
		for(int i = 0; i < (int)nVertices; i++)
		{
			for(int j = i + 1; j < (int)nVertices; j++)
			{
				if(connectedComponents[i] == connectedComponents[j])
				{
					if(componentIDs[i] != componentIDs[j]) throw std::runtime_error("Entry componentIDs was incorrect");
				}
				else
				{
					if(componentIDs[i] == componentIDs[j]) throw std::runtime_error("Entry componentIDs was incorrect");
				}
			}
		}
#endif
	}
	void cliqueTreeAdjacencyMatrix::reset()
	{
		nVertices = 0;
		cliqueGraph.clear();
		for(int i = 0; i < nMaxVertices; i++) 
		{
			bitsetType& currentVertexContents = boost::get(boost::vertex_name, cliqueGraph, i).contents;
			currentVertexContents.reset();
		}
#ifdef TRACK_GRAPH
#ifndef USE_ADJACENCY_MATRIX_FOR_GRAPH
		graph = graphType();
#else 
		graph.clear();
#endif
#endif
		remainingCliqueTreeVertices.resize(nMaxVertices);
		componentIDs.clear();
		for(int i = 0; i < nMaxVertices; i++) remainingCliqueTreeVertices[nMaxVertices - 1 - i] = i;
		verticesToCliqueVertices.clear();
	}
	void cliqueTreeAdjacencyMatrix::swap(cliqueTreeAdjacencyMatrix& other)
	{
		cliqueGraph.swap(other.cliqueGraph);
#ifdef TRACK_GRAPH
		graph.swap(other.graph);
#endif
		std::swap(nVertices, other.nVertices);
		std::swap(nMaxVertices, other.nMaxVertices);
		verticesToCliqueVertices.swap(other.verticesToCliqueVertices);
		componentIDs.swap(other.componentIDs);
		remainingCliqueTreeVertices.swap(other.remainingCliqueTreeVertices);
	}
	cliqueTreeAdjacencyMatrix& cliqueTreeAdjacencyMatrix::operator=(cliqueTreeAdjacencyMatrix&& other)
	{
		cliqueGraph = std::move(other.cliqueGraph); 
#ifdef TRACK_GRAPH
		graph = std::move(other.graph); 
#endif
		nVertices = other.nVertices;
		nMaxVertices = other.nMaxVertices;
		verticesToCliqueVertices = std::move(other.verticesToCliqueVertices);
		componentIDs = std::move(other.componentIDs);
		remainingCliqueTreeVertices = std::move(other.remainingCliqueTreeVertices);
		return *this;
	}
#ifdef HAS_NAUTY
	void cliqueTreeAdjacencyMatrix::convertToNauty(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<::graph>& nautyGraph, std::vector<::graph>& cannonicalNautyGraph)
	{
		int n = nVertices;
		int m = SETWORDSNEEDED(n);
		cannonicalNautyGraph.resize(n * m);
		convertToNauty(lab, ptn, orbits, nautyGraph, &(cannonicalNautyGraph[0]));
	}
	void cliqueTreeAdjacencyMatrix::convertToNauty(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<::graph>& nautyGraph, ::graph* cannonicalNautyGraph)
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
		cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_iterator current, end;
		boost::tie(current, end) = boost::vertices(cliqueGraph);
		for (; current != end; current++)
		{
			bitsetType bitset = boost::get(boost::vertex_name, cliqueGraph, *current).contents;
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					if (bitset[i] && bitset[j])
					{
						ADDONEEDGE(&(nautyGraph[0]), i, j, m);
					}
				}
			}
		}
		densenauty(&(nautyGraph[0]), &(lab[0]), &(ptn[0]), &(orbits[0]), &options, &stats, m, n, cannonicalNautyGraph);
	}
	void cliqueTreeAdjacencyMatrix::convertToNautyWithEdge(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<::graph>& nautyGraph, std::vector<::graph>& cannonicalNautyGraph, int v1, int v2)
	{
		int n = nVertices;
		int m = SETWORDSNEEDED(n);
		cannonicalNautyGraph.resize(n * m);
		convertToNautyWithEdge(lab, ptn, orbits, nautyGraph, &(cannonicalNautyGraph[0]), v1, v2);
	}
	void cliqueTreeAdjacencyMatrix::convertToNautyWithEdge(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<::graph>& nautyGraph, ::graph* cannonicalNautyGraph, int v1, int v2)
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
		cliqueTreeAdjacencyMatrix::cliqueTreeGraphType::vertex_iterator current, end;
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
		densenauty(&(nautyGraph[0]), &(lab[0]), &(ptn[0]), &(orbits[0]), &options, &stats, m, n, cannonicalNautyGraph);
	}
#endif
	void cliqueTreeAdjacencyMatrix::formRemovalTree(std::vector<int>& stateCounts, cliqueTreeAdjacencyMatrix& copied, int u, int v, std::unordered_set<bitsetType>& uniqueSubsets, formRemovalTreeTemporaries& temporaries)
	{
		stateCounts.resize(nVertices);

		std::vector<stackEntry>& stack = temporaries.stack;
		stack.clear();
		
		uniqueSubsets.clear();

		std::vector<int>& counts1 = temporaries.counts1;
		std::vector<int>& counts2 = temporaries.counts2;
		counts1.resize(nVertices);
		counts2.resize(nVertices);

		std::vector<boost::default_color_type>& colourVector = temporaries.colourVector;
		colourVector.resize(nVertices);

		std::fill(stateCounts.begin(), stateCounts.end(), 0);
		copied.makeCopy(*this);

		bitsetType startSet;
		startSet[v] = true;

		bitsetType excluded;
		excluded[u] = excluded[v] = true;

		stack.push_back(stackEntry(startSet, excluded, bitsetType(), v, 0));
		int cliqueVertex;
		copied.canRemoveEdge(u, v, counts1, cliqueVertex);
		copied.removeEdgeKnownCliqueVertex(u, v, colourVector, counts2, cliqueVertex);
		stateCounts[0] = 1;
		while(stack.size() > 0)
		{
			bitsetType unionMinimalSeparatorBitset;
			stackEntry& current = stack.back();
			if(current.currentSearchStart == 0)
			{
				copied.canRemoveEdge(u, v, counts2, cliqueVertex);
				bitsetType& countsToBitset = current.countsToBitset;
				countsToBitset.reset();
				for(int currentVertex = 0; currentVertex != nVertices; currentVertex++)
				{
					if(counts2[currentVertex] == 1 && counts1[currentVertex] != 1 && !current.excluded[currentVertex])
					{
						countsToBitset[currentVertex] = true;
					}
				}
			}
			for(int currentVertex = current.currentSearchStart; currentVertex != nVertices; currentVertex++)
			{
				if(current.countsToBitset[currentVertex])
				{
					bitsetType newSubset = current.currentSet;
					newSubset[currentVertex] = true;
					if(uniqueSubsets.find(newSubset) == uniqueSubsets.end())
					{
						stateCounts[newSubset.count() - 1]++;
						uniqueSubsets.insert(newSubset);
						current.currentSearchStart = currentVertex+1;
						stack.push_back(stackEntry(newSubset, current.excluded, bitsetType(), currentVertex, 0));
						copied.tryRemoveEdge(u, currentVertex, temporaries.colourVector, temporaries.counts2);
						goto continueWhileLoop;
					}
				}
			}
			copied.addEdge(u, current.lastVertex, unionMinimalSeparatorBitset, temporaries.vertexSequence, temporaries.edgeSequence, temporaries.addEdges, temporaries.removeEdges, temporaries.unionTemporaries, false);
			stack.pop_back();
continueWhileLoop:
			;
		}
	}
}
