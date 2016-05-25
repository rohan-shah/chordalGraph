#ifndef PERFORM_SAMPLING_HEADER_GUARD
#define PERFORM_SAMPLING_HEADER_GUARD
#include "sampford.h"
#include "numericType.h"
#include "semiDeterministicSampling.h"
#include "conditionalPoisson.h"
#include "samplingType.h"
namespace chordalGraph
{
	struct performSamplingArgs
	{
		performSamplingArgs()
		{
			conditionalArgs.calculateInclusionProbabilities = true;
		}
		int toTake;
		//Arguments for calling the sampford sampling function
		sampling::sampfordMultinomialRejectiveArgs sampfordMultinomialArgs;
		//Arguments for calling the conditional poisson sampling function
		sampling::conditionalPoissonArgs conditionalArgs;
		//Arguments for semi-deterministic sampling
		sampling::semiDeterministicSamplingArgs semiDetArgs;
		//Arguments for calling the sampford conditional poisson rejective sampling
		sampling::sampfordConditionalPoissonRejectiveArgs sampfordConditionalPoissonArgs;
		//Arguments for performing pareto sampling and pretending it's a sampford sample
		sampling::sampfordFromParetoNaiveArgs sampfordFromParetoArgs;
		samplingType sampling;
	};
	void performSampling(performSamplingArgs& args, std::vector<int>& indices, std::vector<numericType>& inclusionProbabilities, std::vector<numericType>& weights, boost::mt19937& randomSource);
}
#endif
