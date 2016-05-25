#include "conditionalPoisson.h"
#include "numericType.h"
#include <iostream>
namespace chordalGraph
{
	using std::exp;
	using boost::multiprecision::exp;
	int main(int argc, char** argv)
	{
		boost::mt19937 randomSource;
		std::vector<int> indices;
		std::vector<numericType> inclusionProbabilities;
		std::vector<numericType> weights;
		::sampling::conditionalPoissonArgs args;
		randomSource.seed(1);

		args.exponentialParameters.push_back(-2.151);
		args.exponentialParameters.push_back(-1.221);
		args.exponentialParameters.push_back(-0.211);
		args.exponentialParameters.push_back(0.344);
		args.exponentialParameters.push_back(1.285);
		args.exponentialParameters.push_back(1.954);

		args.deterministicInclusion.resize(args.exponentialParameters.size());
		std::fill(args.deterministicInclusion.begin(), args.deterministicInclusion.end(), false);

		for(int i = 0; i < (int)args.exponentialParameters.size(); i++)
		{
			args.expExponentialParameters.push_back(exp(args.exponentialParameters[i]));
		}
		::sampling::calculateExpNormalisingConstants(args.expExponentialParameters, args.exponentialParameters, args.expNormalisingConstant, (int)args.exponentialParameters.size(), (int)args.exponentialParameters.size(), args.deterministicInclusion);
		for(int i = 0; i < 6; i++)
		{
			for(int j = 0; j < 6; j++)
			{
				std::cout << double(args.expNormalisingConstant(i, j)) << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl << std::endl;



		args.n = 3;
		weights.clear();
		weights.push_back(1);
		weights.push_back(1);
		weights.push_back(2);
		weights.push_back(2);
		weights.push_back(2);
		weights.push_back(4);
		const int nSamples = 1000;
		std::vector<int> counts(weights.size());
		for(int i = 0; i < nSamples; i++)
		{
			std::vector<numericType> weightsCopy = weights;
			sampling::conditionalPoisson(args, indices, inclusionProbabilities, weightsCopy, randomSource);
			for(std::vector<int>::iterator j = indices.begin(); j != indices.end(); j++)
			{
				counts[*j]++;
			}
		}
		std::cout << "Inclusion probabilities:" << std::endl;
		for(std::vector<numericType>::iterator j = inclusionProbabilities.begin(); j != inclusionProbabilities.end(); j++)
		{
			std::cout << *j << " ";
		}
		std::cout << std::endl;
		std::cout << "Counts:" << std::endl;
		for(std::vector<int>::iterator j = counts.begin(); j != counts.end(); j++)
		{
			std::cout << *j << " ";
		}
		std::cout << std::endl;
		return 0;
	}
}
int main(int argc, char** argv)
{
	return chordalGraph::main(argc, argv);
}
