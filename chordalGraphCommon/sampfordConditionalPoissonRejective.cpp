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
	void sampfordConditionalPoissonRejective(sampfordConditionalPoissonRejectiveArgs& args)
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
				cumulative += args.weights[i];
				args.accumulated[i] = double(cumulative);
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
				args.inclusionProbabilities[i] = args.weights[i] = args.weights[i]*factor;
			}
		}
		boost::random::uniform_real_distribution<> firstSampleDist(0, (double)cumulative);
beginSample:
		args.indices.resize(deterministicIndices);
		double firstSample = firstSampleDist(args.randomSource);
		int firstIndex = (int)std::distance(args.accumulated.begin(), std::upper_bound(args.accumulated.begin(), args.accumulated.end(), firstSample, std::less_equal<double>()));
		if(args.deterministicInclusion[firstIndex] || args.indices.size() == args.n)
		{
			throw std::runtime_error("Internal error");
		}

		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i])
			{
				boost::random::bernoulli_distribution<double> currentUnitDist((double)args.weights[i]);
				if(currentUnitDist(args.randomSource))
				{
					args.indices.push_back(i);
				}
			}
			//If we can't reach the target of n units, then start again. 
			if(args.indices.size() > args.n-1 || args.indices.size() + (nUnits - i - 1) < args.n-1) goto beginSample;
		}
		if(args.indices.size() != args.n-1)
		{
			throw std::runtime_error("Internal error");
		}
		if(args.n > 1 && std::find(args.indices.begin(), args.indices.end(), firstIndex) != args.indices.end())
		{
			goto beginSample;
		}
		args.indices.push_back(firstIndex);
		std::sort(args.indices.begin(), args.indices.end());
		if(std::unique(args.indices.begin(), args.indices.end()) != args.indices.end())
		{
			throw std::runtime_error("Internal error");
		}
	}
}
