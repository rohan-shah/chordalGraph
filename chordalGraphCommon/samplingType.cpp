#include "samplingType.h"
#include <stdexcept>
namespace chordalGraph
{
	samplingType toSamplingType(std::string samplingString)
	{
		if(samplingString == "sampfordMultinomial")
		{
			return sampfordSamplingMultinomial;
		}
		else if(samplingString == "sampfordConditionalPoisson")
		{
			return sampfordSamplingConditionalPoisson;
		}
		else if(samplingString == "conditionalPoisson")
		{
			return conditionalPoissonSampling;
		}
		else if(samplingString == "sampfordFromParetoNaive")
		{
			return sampfordSamplingFromParetoNaive;
		}
		else if(samplingString == "pareto")
		{
			return paretoSampling;
		}
		else if(samplingString == "semiDeterministic")
		{
			return semiDeterministicSampling;
		}
		else
		{
			throw std::runtime_error("Sampling type must be one of \"sampfordMultinomial\", \"sampfordConditionalPoisson\", \"conditionalPoisson\", \"sampfordFromParetoNaive\", \"pareto\" or \"semiDeterministic\"");
		}
	}
}
