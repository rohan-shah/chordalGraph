#include "conditionalPoisson.h"
#include <iostream>
namespace chordalGraph
{
	using std::exp;
	using boost::multiprecision::exp;
	int main(int argc, char** argv)
	{
		boost::mt19937 randomSource;
		conditionalPoissonArgs args(randomSource);
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
		calculateExpNormalisingConstants(args.expExponentialParameters, args.exponentialParameters, args.expNormalisingConstant, (int)args.exponentialParameters.size(), (int)args.exponentialParameters.size(), args.deterministicInclusion);
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
		std::vector<numericType> weights;
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
			args.weights = weights;
			conditionalPoisson(args);
			for(std::vector<int>::iterator j = args.indices.begin(); j != args.indices.end(); j++)
			{
				counts[*j]++;
			}
		}
		std::cout << "Inclusion probabilities:" << std::endl;
		for(std::vector<numericType>::iterator j = args.inclusionProbabilities.begin(); j != args.inclusionProbabilities.end(); j++)
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
