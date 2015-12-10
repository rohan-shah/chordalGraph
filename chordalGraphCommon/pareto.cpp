#include "pareto.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <functional>
namespace chordalGraph
{
	using boost::multiprecision::log;
	using std::log;
	using boost::multiprecision::exp;
	using std::exp;
	void paretoSampling(paretoSamplingArgs& args, std::vector<int>& indices, std::vector<numericType>& inclusionProbabilities, std::vector<numericType>& weights, boost::mt19937& randomSource)
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
				}
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
				weights[i] = weights[i]*factor;
			}
		}
		boost::random::uniform_real_distribution<> standardUniform(0, 1);
		//Now compute the pareto statistics
		args.paretoStatistics.clear();
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i])
			{
				double uniform = standardUniform(randomSource);
				paretoSamplingArgs::paretoStatistic newStatistic;
				mpfr_class value = ((uniform * (1 - weights[i]))/(weights[i]*(1-uniform)));
				newStatistic.statistic = value.convert_to<double>();
				newStatistic.order = i;
				args.paretoStatistics.push_back(newStatistic);
			}
		}
		std::sort(args.paretoStatistics.begin(), args.paretoStatistics.end());
		//Select so many smallest values
		for(int i = 0; i < (int)(args.n - deterministicIndices); i++)
		{
			indices.push_back(args.paretoStatistics[i].order);
		}
		
		if(args.calculateInclusionProbabilities)
		{
			throw std::runtime_error("Calculation of Pareto inclusion probabilities not implemented yet");
		}
	}
}
