#ifndef SAMPLING_TYPE_HEADER_GUARD
#define SAMPLING_TYPE_HEADER_GUARD
#include <string>
namespace chordalGraph
{
	enum samplingType
	{
		sampfordSamplingMultinomial, conditionalPoissonSampling, paretoSampling, semiDeterministicSampling, sampfordSamplingConditionalPoisson, sampfordSamplingFromParetoNaive
	};
	samplingType toSamplingType(std::string samplingString);
}
#endif
