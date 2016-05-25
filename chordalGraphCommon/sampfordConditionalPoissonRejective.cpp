#include "sampford.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <functional>
namespace chordalGraph
{
	using boost::multiprecision::log;
	using std::log;
	using boost::multiprecision::exp;
	using std::exp;
	void sampfordConditionalPoissonRejective(sampfordConditionalPoissonRejectiveArgs& args, std::vector<int>& indices, std::vector<numericType>& inclusionProbabilities, std::vector<numericType>& weights, boost::mt19937& randomSource)
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
		args.accumulated.resize(nUnits);
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
				args.accumulated[i] = double(cumulative);
			}
			numericType maxAllowed = cumulative / numericType(args.n - indices.size());
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
			if(indices.size() > args.n)
			{
				throw std::runtime_error("Internal error");
			}
		} while(hasDeterministic);
		int deterministicIndices = (int)indices.size();

		//Rescale the weights so that they sum to n
		numericType factor = numericType(args.n - deterministicIndices)/ cumulative;
		if(cumulative == 0)
		{
			throw std::runtime_error("Divide by zero encountered");
		}
		//And also work out the exponential parameters
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i])
			{
				inclusionProbabilities[i] = weights[i] = weights[i]*factor;
			}
		}
		boost::random::uniform_real_distribution<> firstSampleDist(0, (double)cumulative);
beginSample:
		indices.resize(deterministicIndices);
		double firstSample = firstSampleDist(randomSource);
		int firstIndex = (int)std::distance(args.accumulated.begin(), std::upper_bound(args.accumulated.begin(), args.accumulated.end(), firstSample, std::less_equal<double>()));
		if(args.deterministicInclusion[firstIndex] || indices.size() == args.n)
		{
			//In rare cases this can fail, e.g. getting a random number of exactly 0. 
			goto beginSample;
		}

		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i])
			{
				boost::random::bernoulli_distribution<double> currentUnitDist((double)weights[i]);
				if(currentUnitDist(randomSource))
				{
					indices.push_back(i);
				}
			}
			//If we can't reach the target of n units, then start again. 
			if(indices.size() > args.n-1 || indices.size() + (nUnits - i - 1) < args.n-1) goto beginSample;
		}
		if(indices.size() != args.n-1)
		{
			throw std::runtime_error("Internal error");
		}
		if(args.n > 1 && std::find(indices.begin(), indices.end(), firstIndex) != indices.end())
		{
			goto beginSample;
		}
		indices.push_back(firstIndex);
		std::sort(indices.begin(), indices.end());
		if(std::unique(indices.begin(), indices.end()) != indices.end())
		{
			throw std::runtime_error("Internal error");
		}
	}
}
