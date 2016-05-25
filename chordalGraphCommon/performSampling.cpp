#include "performSampling.h"
namespace chordalGraph
{
	void performSampling(performSamplingArgs& args, std::vector<int>& indices, std::vector<numericType>& inclusionProbabilities, std::vector<numericType>& weights, boost::mt19937& randomSource)
	{
		if(args.sampling == sampfordSamplingConditionalPoisson)
		{
			args.sampfordConditionalPoissonArgs.n = args.toTake;
			sampling::sampfordConditionalPoissonRejective(args.sampfordConditionalPoissonArgs, indices, inclusionProbabilities, weights, randomSource);
		}
		else if(args.sampling == semiDeterministicSampling)
		{
			args.semiDetArgs.n = args.toTake;
			sampling::semiDeterministic(args.semiDetArgs, indices, inclusionProbabilities, weights, randomSource);
		}
		//The multinomial rejective version of sampford sampling
		else if(args.sampling == sampfordSamplingMultinomial)
		{
			args.sampfordMultinomialArgs.n = args.toTake;
			sampling::sampfordMultinomialRejective(args.sampfordMultinomialArgs, indices, inclusionProbabilities, weights, randomSource);
		}
		else if(args.sampling == conditionalPoissonSampling)
		{
			args.conditionalArgs.n = args.toTake;
			sampling::conditionalPoisson(args.conditionalArgs, indices, inclusionProbabilities, weights, randomSource);
		}
		else if(args.sampling == sampfordSamplingFromParetoNaive)
		{
			args.sampfordFromParetoArgs.n = args.toTake;
			sampling::sampfordFromParetoNaive(args.sampfordFromParetoArgs, indices, inclusionProbabilities, weights, randomSource);
		}
		else
		{
			throw std::runtime_error("This type of sampling is still unsupported");
		}
	}
}
