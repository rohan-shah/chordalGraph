#ifndef CLIQUE_TREE_ADJACENCY_MATRIX_HEADER_GUARD
#define CLIQUE_TREE_ADJACENCY_MATRIX_HEADER_GUARD
#include <bitset>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#define MAX_STORAGE_VERTICES 64
#define USE_ADJACENCY_MATRIX_FOR_GRAPH
#ifdef HAS_NAUTY
#include "nauty.h"
#endif
#include <boost/multiprecision/mpfr.hpp>
#include "cliqueTree.h"
#ifndef NDEBUG
#define TRACK_GRAPH
#endif
namespace chordalGraph
{
	class cliqueTreeAdjacencyMatrix
	{
	public:
		typedef moveable_adjacency_matrix<boost::property<boost::vertex_name_t, cliqueVertex>, boost::no_property> cliqueTreeGraphType;
#ifdef USE_ADJACENCY_MATRIX_FOR_GRAPH
		typedef moveable_adjacency_matrix<> graphType;
#else
		typedef moveable_adjacency_list<> graphType;
#endif
		struct unionMinimalSeparatorsTemporaries
		{
			std::vector<boost::default_color_type> colorMap;
			std::vector<cliqueTreeGraphType::edge_descriptor> predecessorEdges;
			std::vector<cliqueTreeGraphType::vertex_descriptor> predecessorVertices;
		};
		struct externalEdge
		{
		public:
			externalEdge(int source, int target)
				:source(source), target(target)
			{}
			int source, target;
		};
	public:
		cliqueTreeAdjacencyMatrix(cliqueTreeAdjacencyMatrix&& other)
			:cliqueGraph(std::move(other.cliqueGraph)), 
#ifdef TRACK_GRAPH
			graph(std::move(other.graph)), 
#endif
			nVertices(other.nVertices), nMaxVertices(other.nMaxVertices),
			verticesToCliqueVertices(std::move(other.verticesToCliqueVertices)), componentIDs(std::move(other.componentIDs)), remainingCliqueTreeVertices(std::move(other.remainingCliqueTreeVertices))
		{}
		cliqueTreeAdjacencyMatrix& operator=(cliqueTreeAdjacencyMatrix&& other);
		cliqueTreeAdjacencyMatrix(int maximumVertices);
		cliqueTreeAdjacencyMatrix(const cliqueTreeAdjacencyMatrix& other)
			:cliqueGraph(other.cliqueGraph), 
#ifdef TRACK_GRAPH
			graph(other.graph), 
#endif
			nVertices(other.nVertices), nMaxVertices(other.nMaxVertices),
			verticesToCliqueVertices(other.verticesToCliqueVertices), componentIDs(other.componentIDs), remainingCliqueTreeVertices(other.remainingCliqueTreeVertices)
		{}
		void makeCopy(const cliqueTreeAdjacencyMatrix& other);
		const cliqueTreeGraphType& getCliqueGraph() const;
		bool canRemoveEdge(int u, int v, std::vector<int>& counts, int& cliqueVertex);
		bool tryRemoveEdge(int u, int v, std::vector<boost::default_color_type>& colourVector, std::vector<int>& countsAfter);
		bool removeEdgeKnownCliqueVertex(int u, int v, std::vector<boost::default_color_type>& colourVector, std::vector<int>& counts2, int cliqueVertex);
#ifdef TRACK_GRAPH
		const graphType& getGraph() const;
#endif
		void swap(cliqueTreeAdjacencyMatrix& other);
		bool tryAddVertexWithEdges(const bitsetType& involvedEdges, unionMinimalSeparatorsTemporaries& temp);
		void addVertex();
		int getNVertices();
		void unionMinimalSeparators(int u, int v, bitsetType& vertices, std::list<cliqueTreeGraphType::vertex_descriptor>& vertexSequence, std::list<externalEdge>& edgeSequence, std::vector<externalEdge>& addEdges, std::vector<externalEdge>& removeEdges, unionMinimalSeparatorsTemporaries& temp);
		void addEdge(int vertexForExtraEdges, int v, bitsetType& unionMinimalSeparatorBitset, std::list<cliqueTreeGraphType::vertex_descriptor>& vertexSequence, std::list<externalEdge>& edgeSequence, std::vector<externalEdge>& addEdges, std::vector<externalEdge>& removeEdges, unionMinimalSeparatorsTemporaries& temp, bool hasPrecomputedUnionMinimalSeparator);
		void check() const;
#ifdef HAS_NAUTY
		void convertToNauty(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<::graph>& nautyGraph, std::vector<::graph>& cannonicalNautyGraph);
		void convertToNautyWithEdge(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<::graph>& nautyGraph, std::vector<::graph>& cannonicalNautyGraph, int v1, int v2);
		void convertToNauty(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<::graph>& nautyGraph, ::graph* cannonicalNautyGraph);
		void convertToNautyWithEdge(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<::graph>& nautyGraph, ::graph* cannonicalNautyGraph, int v1, int v2);
		//void convertToNautyAndCountAutomorphisms(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<graph>& nautyGraph, std::vector<graph>& cannonicalNautyGraph, mpz_class& automorphismCount);
		static void userlevelproc(int* lab, int* ptn, int level, int* orbits, statsblk* stats, int tv, int index, int tcellsize, int numcells, int childcount, int n);
#endif
		int getNVertices() const;
		int maxVertices() const;
		void reset();
		struct stackEntry
		{
		public:
			stackEntry(bitsetType currentSet, bitsetType excluded, bitsetType countsToBitset, int lastVertex, int currentSearchStart)
				: currentSet(currentSet), excluded(excluded), countsToBitset(countsToBitset), lastVertex(lastVertex), currentSearchStart(currentSearchStart)
			{}
			bitsetType currentSet, excluded, countsToBitset;
			int lastVertex, currentSearchStart;
		};
		struct formRemovalTreeTemporaries
		{
			std::vector<stackEntry> stack;
			std::vector<int> counts1, counts2;
			std::vector<boost::default_color_type> colourVector;
			std::list<cliqueTreeGraphType::vertex_descriptor> vertexSequence;
			std::list<externalEdge> edgeSequence;
			std::vector<externalEdge> addEdges;
			std::vector<externalEdge> removeEdges;
			unionMinimalSeparatorsTemporaries unionTemporaries;
		};
		void formRemovalTree(std::vector<int>& outputCounts, cliqueTreeAdjacencyMatrix& copied, int u, int v, std::unordered_set<bitsetType>& uniqueSubsets, formRemovalTreeTemporaries& temporaries);
	private:
		cliqueTreeGraphType cliqueGraph;
#ifdef TRACK_GRAPH
		graphType graph;
#endif
		int nVertices;
		int nMaxVertices;
		//A vector that converts any vertex of the graph to a vertex of the clique tree
		//which contains that vertex in the relevant subset. 
		std::vector<int> verticesToCliqueVertices;

		//Store a list of the connected components
		std::vector<int> componentIDs;
		//clique tree vertices that can still be used
		std::vector<int> remainingCliqueTreeVertices;
	};
}
#endif
