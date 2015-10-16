#include "sampford.h"
#include <boost/random/uniform_real_distribution.hpp>
#include <functional>
namespace chordalGraph
{
	void sampfordBruteForce(sampfordBruteForceArgs& args)
	{
		std::size_t nUnits = args.weights.size();
		if(args.n > nUnits)
		{
			throw std::runtime_error("Input n was too big");
		}
		else if(args.n == nUnits)
		{
			args.indices.clear();
			args.indices.reserve(nUnits);
			for(int i = 0; i < (int)nUnits; i++)
			{
				args.indices.push_back(i);
			}
			return;
		}
		args.accumulated1.resize(nUnits);
		args.accumulated2.resize(nUnits);
		std::fill(args.accumulated1.begin(), args.accumulated1.end(), 0);
		std::fill(args.accumulated2.begin(), args.accumulated2.end(), 0);

		args.inclusionProbabilities.resize(nUnits);
		std::fill(args.inclusionProbabilities.begin(), args.inclusionProbabilities.end(), 0);
		args.deterministicIndices.clear();
		args.deterministicInclusion.resize(nUnits);
		std::fill(args.deterministicInclusion.begin(), args.deterministicInclusion.end(), false);
		double cumulative;
		bool hasDeterministic = false;
		do
		{
			hasDeterministic = false;
			//Work out sum of weights
			cumulative = 0;
			for(int i = 0; i < (int)nUnits; i++)
			{
				cumulative += args.weights[i];
			}
			double maxAllowed = (double)cumulative / (double)(args.n - args.deterministicIndices.size());
			//Any weights that are too big are included with probability 1
			for(int i = 0; i < (int)nUnits; i++)
			{
				if(args.weights[i] >= maxAllowed)
				{
					args.deterministicInclusion[i] = true;
					args.deterministicIndices.push_back(i);
					hasDeterministic = true;
					args.weights[i] = 0;
					args.inclusionProbabilities[i] = 1;
				}
			}
		} while(hasDeterministic);
		int deterministicIndices = (int)args.deterministicIndices.size();
		//Rescale weights so they sum to args.n. These are the weights used for the first sample
		double rescaledCumulative = 0;
		double factor = (double)(args.n - deterministicIndices)/ (double)cumulative;
		for(int i = 0; i < (int)nUnits; i++)
		{
			double prob;
			if(!args.deterministicInclusion[i])
			{
				prob = args.inclusionProbabilities[i] = args.weights[i]*factor;
			}
			else prob = 0;
			rescaledCumulative += prob;
			args.accumulated1[i] = rescaledCumulative;
		}
		if(rescaledCumulative == std::numeric_limits<double>::infinity())
		{
			throw std::runtime_error("Infinity encountered in sampford sampling");
		}

		//Accumulated weights for the other samples
		double cumulativeOther = 0;
		for(int i = 0; i < (int)nUnits; i++)
		{
			double rescaled = args.weights[i]*factor;
			double asProb = rescaled / (1 - rescaled);
			cumulativeOther += asProb;
			args.accumulated2[i] = cumulativeOther;
		}
		if(cumulativeOther == std::numeric_limits<double>::infinity())
		{
			throw std::runtime_error("Infinity encountered in sampford sampling");
		}
		if(args.deterministicIndices.size() >= args.n)
		{
			throw std::runtime_error("Deterministic sample taken");
		}
		boost::random::uniform_real_distribution<> firstSampleDist(0, rescaledCumulative);
getSample:
		args.indices = args.deterministicIndices;
		double firstSample = firstSampleDist(args.randomSource);
		int firstIndex = std::distance(args.accumulated1.begin(), std::upper_bound(args.accumulated1.begin(), args.accumulated1.end(), firstSample, std::less_equal<double>()));
		if(args.deterministicInclusion[firstIndex] || args.deterministicIndices.size() == args.n)
		{
			throw std::runtime_error("Internal error");
		}
		args.indices.push_back(firstIndex);
		{
			boost::random::uniform_real_distribution<> otherSamplesDist(0, cumulativeOther);
			for(int i = 0; i < (int)args.n-deterministicIndices-1; i++)
			{
				double sampledWeight = otherSamplesDist(args.randomSource);
				int index = std::distance(args.accumulated2.begin(), std::upper_bound(args.accumulated2.begin(), args.accumulated2.end(), sampledWeight, std::less_equal<double>()));
				if(args.deterministicInclusion[index])
				{
					throw std::runtime_error("Internal error");
				}
				args.indices.push_back(index);
			}
			std::sort(args.indices.begin(), args.indices.end());
			if(std::unique(args.indices.begin(), args.indices.end()) != args.indices.end())
			{
				goto getSample;
			}
		}
	}
}
