#include "conditionalPoisson.h"
#include <boost/random/bernoulli_distribution.hpp>
namespace chordalGraph
{
	using boost::multiprecision::log;
	using std::log;
	using boost::multiprecision::exp;
	using std::exp;
	void conditionalPoisson(conditionalPoissonArgs& args)
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

		//Rescale the weights so that they sum to n
		numericType factor = numericType(args.n - deterministicIndices)/ cumulative;
		if(cumulative == 0)
		{
			throw std::runtime_error("Divide by zero encountered");
		}
		//And also work out the exponential parameters
		args.exponentialParameters.resize(nUnits);
		args.expExponentialParameters.resize(nUnits);
		numericType sumExponentialParameters = 0;
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i])
			{
				args.weights[i] = args.weights[i]*factor;
				args.expExponentialParameters[i] = args.weights[i] / (1 - args.weights[i]);
				args.exponentialParameters[i] = log(args.expExponentialParameters[i]);
				sumExponentialParameters += args.exponentialParameters[i];
			}
		}
		//Rescale so the exponential parameters sum no zero
		numericType toSubtract = sumExponentialParameters / nUnits;
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i])
			{
				args.exponentialParameters[i] -= toSubtract;
				args.expExponentialParameters[i] = exp(args.exponentialParameters[i]);
			}
		}
beginSample:
		args.indices.resize(deterministicIndices);
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
			if(args.indices.size() > args.n || args.indices.size() + (nUnits - i - 1) < args.n) goto beginSample;
		}
		if(args.indices.size() != args.n)
		{
			throw std::runtime_error("Internal error");
		}
		//Now compute the inclusion probabilities
		args.inclusionProbabilities.resize(nUnits);
		calculateExpNormalisingConstants(args.expExponentialParameters, args.exponentialParameters, args.expNormalisingConstant, (int)args.n - deterministicIndices, nUnits - deterministicIndices, args.deterministicInclusion);
		numericType expNormalisingConstant = args.expNormalisingConstant(0, args.n - deterministicIndices - 1);
		for(int unitCounter = 0; unitCounter < nUnits; unitCounter++)
		{
			if(args.deterministicInclusion[unitCounter]) continue;
			//First term of the alternating sum
			if((args.n-deterministicIndices)% 2 == 1)
			{
				args.inclusionProbabilities[unitCounter] = 1;
			}
			else
			{
				args.inclusionProbabilities[unitCounter] = -1;
			}
			for(int j = 2; j <= (int)args.n-deterministicIndices; j++)
			{
				if((args.n -deterministicIndices- j) % 2 == 1)
				{
					args.inclusionProbabilities[unitCounter] -= args.expNormalisingConstant(0, j-2) / exp((j - 1)*args.exponentialParameters[unitCounter]);
				}
				else
				{
					args.inclusionProbabilities[unitCounter] += args.expNormalisingConstant(0, j-2) / exp((j - 1)*args.exponentialParameters[unitCounter]);
				}
			}
			args.inclusionProbabilities[unitCounter] *= exp((args.n - deterministicIndices) * args.exponentialParameters[unitCounter]) / expNormalisingConstant;
			if(args.inclusionProbabilities[unitCounter] <= 0)
			{
				std::stringstream ss;
				ss << "Inclusion probability had negative value " << args.inclusionProbabilities[unitCounter] << ", probably because of numerical instability";
				throw std::runtime_error(ss.str().c_str());
			}
		}
	}
	void calculateExpNormalisingConstants(std::vector<numericType>& expExponentialParameters, std::vector<numericType>& exponentialParameters, boost::numeric::ublas::matrix<numericType>& expNormalisingConstant, int n, int nUnits, std::vector<bool>& ignore)
	{
		//We start by computing the normalising constants. First index is k, second is z. All indices in this loop are 1 indexed
		expNormalisingConstant.resize(nUnits,n);
		//This will skip over the *ignored* units (the ones that were deterministically included)
		int k = (int)expExponentialParameters.size();
		for(int unitIndex = nUnits; unitIndex >= 1; unitIndex--)
		{
			while(ignore[k-1]) k--;
			for(int z = 1; z <= std::min(n, nUnits-unitIndex+1); z++)
			{
				if(z == 1)
				{
					numericType sum = 0;
					for(int unitIndex2 = k; unitIndex2 <= (int)expExponentialParameters.size(); unitIndex2++)
					{
						if(!ignore[unitIndex2-1]) sum += expExponentialParameters[unitIndex2-1];
					}
					expNormalisingConstant(unitIndex-1, z-1) = sum;
				}
				else if(z == nUnits - unitIndex + 1)
				{
					numericType sum = 0;
					for(int unitIndex2 = k; unitIndex2 <= (int)expExponentialParameters.size(); unitIndex2++)
					{
						if(!ignore[unitIndex2-1]) sum += exponentialParameters[unitIndex2-1];
					}
					expNormalisingConstant(unitIndex-1, z-1) = exp(sum);
				}
				else
				{
					expNormalisingConstant(unitIndex-1, z-1) = expExponentialParameters[k-1] * expNormalisingConstant(unitIndex, z-2) + expNormalisingConstant(unitIndex, z-1);
				}
			}
			k--;
		}
	}
}
