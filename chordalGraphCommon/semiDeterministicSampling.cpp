#include "semiDeterministicSampling.h"
#include <boost/range/algorithm.hpp>
#include <boost/random/random_number_generator.hpp>
#include <algorithm>
#include <functional>
namespace chordalGraph
{
	void semiDeterministic(semiDeterministicSamplingArgs& args)
	{
		args.indices.clear();
		int nUnits = (int)args.weights.size();
		if((int)args.n > nUnits)
		{
			throw std::runtime_error("Input n was too big");
		}
		else if((int)args.n == nUnits)
		{
			args.indices.reserve(nUnits);
			for(int i = 0; i < nUnits; i++)
			{
				args.indices.push_back(i);
			}
			return;
		}
		args.deterministicInclusion.resize(nUnits);
		args.inclusionProbabilities.resize(nUnits);
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
				cumulative += args.weights[i];
			}
			//Make the maxAllowed value slightly smaller. This is because we don't want inclusion probabilities that are numerically close to 1, because in that case the alternating series used to compute the inclusion probabilities becomes unstable. 
			numericType maxAllowed = 0.9999* cumulative / numericType(args.n - args.indices.size());
			//Any weights that are too big are included with probability 1
			for(int i = 0; i < nUnits; i++)
			{
				if(args.weights[i] >= maxAllowed)
				{
					args.deterministicInclusion[i] = true;
					args.indices.push_back(i);
					hasDeterministic = true;
					args.weights[i] = 0;
					args.inclusionProbabilities[i] = 1;
				}
			}
		} while(hasDeterministic);
		int deterministicIndices = (int)args.indices.size();
		for(int i = 0; i < (int)nUnits; i++)
		{
			if(!args.deterministicInclusion[i])
			{
				args.indices.push_back(i);
				args.inclusionProbabilities[i] = (double)(args.n - deterministicIndices) / (double)(nUnits - deterministicIndices);
			}
		}
		boost::random_number_generator<boost::mt19937> generator(args.randomSource);
		std::random_shuffle(args.indices.begin()+deterministicIndices, args.indices.end(), generator);
		args.indices.resize(args.n);
	}
}
