#ifndef CHILD_NODE_HEADER_GUARD
#define CHILD_NODE_HEADER_GUARD
namespace chordalGraph
{
	struct childNode
	{
	public:
		childNode(int parentIndex, bool includesEdge, double weight)
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
		double weight;
	private:
		int value;
	};
}
#endif
