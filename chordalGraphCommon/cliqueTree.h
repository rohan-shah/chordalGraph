#ifndef CLIQUE_TREE_HEADER_GUARD
#define CLIQUE_TREE_HEADER_GUARD
#include <bitset>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#define MAX_STORAGE_VERTICES 64
#define USE_ADJACENCY_MATRIX_FOR_GRAPH
#ifdef HAS_NAUTY
#include "nauty.h"
#endif
#include <boost/multiprecision/mpfr.hpp>
typedef boost::multiprecision::mpz_int mpz_class;
namespace chordalGraph
{
	typedef std::bitset<MAX_STORAGE_VERTICES> bitsetType;
	class cliqueVertex
	{
	public:
		cliqueVertex(bitsetType& inputContents)
			:contents(inputContents)
		{}
		static cliqueVertex createEmpty()
		{
			return cliqueVertex();
		}
		cliqueVertex()
			:contents(0)
		{}
		bitsetType contents;
	};
	class cliqueEdge
	{
	public:
		cliqueEdge(bitsetType contents)
			:contents(contents)
		{}
		cliqueEdge()
			:contents(0)
		{}
		cliqueEdge& operator=(const cliqueEdge& other)
		{
			contents = other.contents;
			return *this;
		}
		bitsetType contents;
	};
}
namespace chordalGraph
{
	template <class VertexProperty = boost::no_property, class EdgeProperty = boost::no_property, class EdgeListS = boost::listS>
	class moveable_adjacency_list : public boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProperty, EdgeProperty, boost::no_property, EdgeListS>
	{
	public:
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProperty, EdgeProperty, boost::no_property, EdgeListS> base;
		typedef typename boost::detail::adj_list_gen<base, boost::vecS, boost::vecS, boost::undirectedS, VertexProperty, EdgeProperty, boost::no_property, EdgeListS>::type detailType;
		typedef moveable_adjacency_list<VertexProperty, EdgeProperty, EdgeListS> type;
		moveable_adjacency_list()
		{}
		moveable_adjacency_list(typename type::base::vertices_size_type num_vertices)
			:base(num_vertices)
		{}
		moveable_adjacency_list(const type& other)
			:base(other)
		{
		}
		moveable_adjacency_list(type&& other)
		{
			detailType::m_vertices.swap(other.detailType::m_vertices);
			detailType::m_edges.swap(other.detailType::m_edges);
		}
		type& operator=(type&& other)
		{
			this->detailType::m_vertices.swap(other.detailType::m_vertices);
			this->detailType::m_edges.swap(other.detailType::m_edges);
			return *this;
		}
	};
	template <typename VertexProperty = boost::no_property, typename EdgeProperty = boost::no_property, typename Allocator = std::allocator<bool> >
	class moveable_adjacency_matrix : public boost::adjacency_matrix<boost::undirectedS, VertexProperty, EdgeProperty, boost::no_property, Allocator>
	{
	public:
		typedef boost::adjacency_matrix<boost::undirectedS, VertexProperty, EdgeProperty, boost::no_property, Allocator> base;
		typedef moveable_adjacency_matrix<VertexProperty, EdgeProperty, Allocator> type;
		moveable_adjacency_matrix()
			:base(0), num_vertices(0)
		{}
		moveable_adjacency_matrix(typename type::base::vertices_size_type maxVertices)
			:base(maxVertices), num_vertices(0)
		{}
		moveable_adjacency_matrix(const type& other)
			:base(other), num_vertices(other.num_vertices)
		{}
		moveable_adjacency_matrix(type&& other)
			: base(0)
		{
			type::base::m_matrix.swap(other.m_matrix);
			type::base::m_vertex_set = other.m_vertex_set;
			type::base::m_vertex_properties.swap(other.m_vertex_properties);
			type::base::m_num_edges = other.m_num_edges;
			num_vertices = other.num_vertices;
		}
		type& operator=(type&& other)
		{
			type::base::m_matrix.swap(other.m_matrix);
			type::base::m_vertex_set = other.m_vertex_set;
			type::base::m_vertex_properties.swap(other.m_vertex_properties);
			type::base::m_num_edges = other.m_num_edges;
			num_vertices = other.num_vertices;
			return *this;
		}
		int num_vertices;
	};
}
namespace boost
{
	template <typename VertexProperty, typename EdgeProperty, typename EdgeListS> typename ::chordalGraph::moveable_adjacency_matrix<VertexProperty, EdgeProperty, EdgeListS>::vertices_size_type num_vertices(const ::chordalGraph::moveable_adjacency_matrix<VertexProperty, EdgeProperty, EdgeListS>& g_)
	{
		return g_.num_vertices;
	}
	template <typename VertexProperty, typename EdgeProperty, typename EdgeListS> typename ::chordalGraph::moveable_adjacency_matrix<VertexProperty, EdgeProperty, EdgeListS>::vertex_descriptor add_vertex(::chordalGraph::moveable_adjacency_matrix<VertexProperty, EdgeProperty, EdgeListS>& g_)
	{
		typename ::chordalGraph::moveable_adjacency_matrix<VertexProperty, EdgeProperty, EdgeListS>::vertex_descriptor ret = g_.num_vertices;
		g_.num_vertices++;
		return ret;
	}
}
namespace chordalGraph
{
	class cliqueTree
	{
	public:
		//typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_name_t, cliqueVertex>, boost::property<boost::edge_name_t, cliqueEdge> > cliqueTreeGraphType;
		//typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graphType;
		typedef moveable_adjacency_list<boost::property<boost::vertex_name_t, cliqueVertex>, boost::property<boost::edge_name_t, cliqueEdge> > cliqueTreeGraphType;
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
			externalEdge(int source, int target, bitsetType contents)
				:source(source), target(target), contents(contents)
			{}
			int source, target;
			bitsetType contents;
		};
	public:
		cliqueTree(cliqueTree&& other)
			:cliqueGraph(std::move(other.cliqueGraph)), 
#ifdef TRACK_GRAPH
			graph(std::move(other.graph)), 
#else
			nVertices(other.nVertices),
#endif
			verticesToCliqueVertices(std::move(other.verticesToCliqueVertices)), componentIDs(std::move(other.componentIDs))
		{}
		cliqueTree(int maximumVertices);
		cliqueTree(const cliqueTree& other)
			:cliqueGraph(other.cliqueGraph), 
#ifdef TRACK_GRAPH
			graph(other.graph), 
#else
			nVertices(other.nVertices),
#endif
			verticesToCliqueVertices(other.verticesToCliqueVertices), componentIDs(other.componentIDs)
		{}
		const cliqueTreeGraphType& getCliqueGraph() const;
#ifdef TRACK_GRAPH
		const graphType& getGraph() const;
#endif
		bool tryAddVertexWithEdges(const bitsetType& involvedEdges, unionMinimalSeparatorsTemporaries& temp);
		void addVertex();
		int getNVertices();
		void unionMinimalSeparators(int u, int v, bitsetType& vertices, std::list<cliqueTreeGraphType::vertex_descriptor>& vertexSequence, std::list<externalEdge>& edgeSequence, std::vector<externalEdge>& addEdges, std::vector<externalEdge>& removeEdges, unionMinimalSeparatorsTemporaries& temp);
		void addEdge(int vertexForExtraEdges, int v, bitsetType& unionMinimalSeparatorBitset, std::list<cliqueTreeGraphType::vertex_descriptor>& vertexSequence, std::list<externalEdge>& edgeSequence, std::vector<externalEdge>& addEdges, std::vector<externalEdge>& removeEdges, unionMinimalSeparatorsTemporaries& temp, bool hasPrecomputedUnionMinimalSeparator);
		void check() const;
#ifdef HAS_NAUTY
		void convertToNauty(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<graph>& nautyGraph, std::vector<graph>& cannonicalNautyGraph);
		//void convertToNautyAndCountAutomorphisms(std::vector<int>& lab, std::vector<int>& ptn, std::vector<int>& orbits, std::vector<graph>& nautyGraph, std::vector<graph>& cannonicalNautyGraph, mpz_class& automorphismCount);
		static void userlevelproc(int* lab, int* ptn, int level, int* orbits, statsblk* stats, int tv, int index, int tcellsize, int numcells, int childcount, int n);
#endif
		int getNVertices() const;
	private:
		cliqueTreeGraphType cliqueGraph;
#ifdef TRACK_GRAPH
		graphType graph;
#else
		int nVertices;
#endif
		//A vector that converts any vertex of the graph to a vertex of the clique tree
		//which contains that vertex in the relevant subset. 
		std::vector<int> verticesToCliqueVertices;

		//Store a list of the connected components
		std::vector<int> componentIDs;
	};
}
#endif
