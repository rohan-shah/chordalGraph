#include "semiDeterministicSampling.h"
#include <boost/range/algorithm.hpp>
#include <boost/random/random_number_generator.hpp>
#include <algorithm>
#include <functional>
namespace chordalGraph
{
	void semiDeterministic(semiDeterministicSamplingArgs& args, std::vector<int>& indices, std::vector<numericType>& inclusionProbabilities, std::vector<numericType>& weights, boost::mt19937& randomSource)
	{
		indices.clear();
		int nUnits = (int)weights.size();
		if((int)args.n > nUnits)
		{
			throw std::runtime_error("Input n was too big");
		}
		else if((int)args.n == nUnits)
		{
			indices.reserve(nUnits);
			for(int i = 0; i < nUnits; i++)
			{
				indices.push_back(i);
			}
			return;
		}
		args.deterministicInclusion.resize(nUnits);
		inclusionProbabilities.resize(nUnits);
		std::fill(args.deterministicInclusion.begin(), args.deterministicInclusion.end(), false);
		//Work out which units are going to be deterministically selected. 
		numericType cumulative;
		bool hasDeterministic = false;
		do
		{
			hasDeterministic = false;
			//Work out sum of weights
			cumulative = 0;
			for(int i = 0; i < nUnits; i++)
			{
				cumulative += weights[i];
			}
			//Make the maxAllowed value slightly smaller. This is because we don't want inclusion probabilities that are numerically close to 1, because in that case the alternating series used to compute the inclusion probabilities becomes unstable. 
			numericType maxAllowed = 0.9999* cumulative / numericType(args.n - indices.size());
			//Any weights that are too big are included with probability 1
			for(int i = 0; i < nUnits; i++)
			{
				if(weights[i] >= maxAllowed)
				{
					args.deterministicInclusion[i] = true;
					indices.push_back(i);
					hasDeterministic = true;
					weights[i] = 0;
					inclusionProbabilities[i] = 1;
				}
			}
		} while(hasDeterministic);
		int deterministicIndices = (int)indices.size();
		for(int i = 0; i < (int)nUnits; i++)
		{
			if(!args.deterministicInclusion[i])
			{
				indices.push_back(i);
				inclusionProbabilities[i] = (double)(args.n - deterministicIndices) / (double)(nUnits - deterministicIndices);
			}
		}
		boost::random_number_generator<boost::mt19937> generator(randomSource);
		std::random_shuffle(indices.begin()+deterministicIndices, indices.end(), generator);
		indices.resize(args.n);
	}
}
