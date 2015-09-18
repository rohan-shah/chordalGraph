#ifndef CLIQUE_TREE_HEADER_GUARD
#define CLIQUE_TREE_HEADER_GUARD
#include <bitset>
#include <boost/graph/adjacency_list.hpp>
#define MAX_STORAGE_VERTICES 64
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
	class cliqueTree
	{
	public:
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_name_t, cliqueVertex>, boost::property<boost::edge_name_t, cliqueEdge> > cliqueTreeGraphType;
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graphType;
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
			:cliqueGraph(std::move(other.cliqueGraph)), graph(std::move(other.graph)), verticesToCliqueVertices(std::move(other.verticesToCliqueVertices)), componentIDs(std::move(other.componentIDs))
		{}
		cliqueTree(int maximumVertices)
		{
			if(maximumVertices > MAX_STORAGE_VERTICES)
			{
				throw std::runtime_error("Requested number of vertices exceeded amount of storage available");
			}
		}
		const cliqueTreeGraphType& getCliqueGraph() const;
		const graphType& getGraph() const;
		bool tryAddVertexWithEdges(const bitsetType& involvedEdges);
		void addVertex();
		int getNVertices();
		void unionMinimalSeparators(int u, int v, bitsetType& vertices, std::list<cliqueTreeGraphType::vertex_descriptor>& vertexSequence, std::list<externalEdge>& edgeSequence, std::vector<externalEdge>& addEdges, std::vector<externalEdge>& removeEdges);
		void addEdge(int vertexForExtraEdges, int v, bitsetType& unionMinimalSeparatorBitset, std::list<cliqueTreeGraphType::vertex_descriptor>& vertexSequence, std::list<externalEdge>& edgeSequence, std::vector<externalEdge>& addEdges, std::vector<externalEdge>& removeEdges, bool hasPrecomputedUnionMinimalSeparator);
		void check() const;
	private:
		cliqueTreeGraphType cliqueGraph;
		graphType graph;
		//A vector that converts any vertex of the graph to a vertex of the clique tree
		//which contains that vertex in the relevant subset. 
		std::vector<int> verticesToCliqueVertices;

		//Store a list of the connected components
		std::vector<int> componentIDs;
	};
}
#endif