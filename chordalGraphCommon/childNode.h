#ifndef CHILD_NODE_HEADER_GUARD
#define CHILD_NODE_HEADER_GUARD
namespace chordalGraph
{
	template <typename weightType> struct childNode
	{
	public:
		childNode(int parentIndex, bool includesEdge, weightType weight)
			:weight(weight), value(2*parentIndex + includesEdge)
		{}
		int getParentIndex() const
		{
			return value / 2;
		}
		bool includesEdge() const 
		{
			return (value % 2) == 1;
		}
		weightType weight;
		bool operator<(const childNode& other) const
		{
			return value < other.value;
		}
	private:
		int value;
	};
}
#endif
