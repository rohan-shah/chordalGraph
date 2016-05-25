#include "sampford.h"
#include <boost/random/uniform_real_distribution.hpp>
#include <functional>
namespace chordalGraph
{
	void sampfordMultinomialRejective(sampfordMultinomialRejectiveArgs& args, std::vector<int>& indices, std::vector<numericType>& inclusionProbabilities, std::vector<numericType>& weights, boost::mt19937& randomSource)
	{
		int nUnits = (int)weights.size();
		if((int)args.n > nUnits)
		{
			throw std::runtime_error("Input n was too big");
		}
		else if((int)args.n == nUnits)
		{
			indices.clear();
			indices.reserve(nUnits);
			for(int i = 0; i < nUnits; i++)
			{
				indices.push_back(i);
			}
			return;
		}
		args.accumulated1.resize(nUnits);
		args.accumulated2.resize(nUnits);
		std::fill(args.accumulated1.begin(), args.accumulated1.end(), 0);
		std::fill(args.accumulated2.begin(), args.accumulated2.end(), 0);

		inclusionProbabilities.resize(nUnits);
		std::fill(inclusionProbabilities.begin(), inclusionProbabilities.end(), 0);
		args.deterministicIndices.clear();
		args.deterministicInclusion.resize(nUnits);
		std::fill(args.deterministicInclusion.begin(), args.deterministicInclusion.end(), false);
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
			numericType maxAllowed = cumulative / numericType(args.n - args.deterministicIndices.size());
			//Any weights that are too big are included with probability 1
			for(int i = 0; i < nUnits; i++)
			{
				if(weights[i] >= maxAllowed)
				{
					args.deterministicInclusion[i] = true;
					args.deterministicIndices.push_back(i);
					hasDeterministic = true;
					weights[i] = 0;
					inclusionProbabilities[i] = 1;
				}
			}
			if(indices.size() > args.n)
			{
				throw std::runtime_error("Internal error");
			}
		} while(hasDeterministic);
		int deterministicIndices = (int)args.deterministicIndices.size();
		//Rescale weights so they sum to args.n. These are the weights used for the first sample
		numericType rescaledCumulative = 0;
		numericType factor = numericType(args.n - deterministicIndices)/ cumulative;
		if(cumulative == 0)
		{
			throw std::runtime_error("Divide by zero encountered");
		}
		for(int i = 0; i < nUnits; i++)
		{
			numericType prob;
			if(!args.deterministicInclusion[i])
			{
				prob = inclusionProbabilities[i] = weights[i]*factor;
			}
			else prob = 0;
			rescaledCumulative += prob;
			args.accumulated1[i] = (double)rescaledCumulative;
		}

		//Accumulated weights for the other samples
		numericType cumulativeOther = 0;
		for(int i = 0; i < nUnits; i++)
		{
			numericType rescaled = weights[i]*factor;
			if(rescaled == 1)
			{
				throw std::runtime_error("Divide by zero encountered");
			}
			numericType asProb = rescaled / (1 - rescaled);
			cumulativeOther += asProb;
			args.accumulated2[i] = (double)cumulativeOther;
		}
		if(args.deterministicIndices.size() >= args.n)
		{
			throw std::runtime_error("Deterministic sample taken");
		}
		boost::random::uniform_real_distribution<> firstSampleDist(0, (double)rescaledCumulative);
getSample:
		indices = args.deterministicIndices;
		double firstSample = firstSampleDist(randomSource);
		int firstIndex = (int)std::distance(args.accumulated1.begin(), std::upper_bound(args.accumulated1.begin(), args.accumulated1.end(), firstSample, std::less_equal<double>()));
		if(args.deterministicInclusion[firstIndex] || args.deterministicIndices.size() == args.n)
		{
			throw std::runtime_error("Internal error");
		}
		indices.push_back(firstIndex);
		{
			boost::random::uniform_real_distribution<> otherSamplesDist(0, (double)cumulativeOther);
			for(int i = 0; i < (int)args.n-deterministicIndices-1; i++)
			{
				double sampledWeight = otherSamplesDist(randomSource);
				int index = (int)std::distance(args.accumulated2.begin(), std::upper_bound(args.accumulated2.begin(), args.accumulated2.end(), sampledWeight, std::less_equal<double>()));
				if(args.deterministicInclusion[index])
				{
					throw std::runtime_error("Internal error");
				}
				indices.push_back(index);
			}
			std::sort(indices.begin(), indices.end());
			if(std::unique(indices.begin(), indices.end()) != indices.end())
			{
				goto getSample;
			}
		}
	}
}
