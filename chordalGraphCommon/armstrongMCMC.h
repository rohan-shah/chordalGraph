#ifndef ARMSTRONG_MCMC_HEADER_GUARD
#define ARMSTRONG_MCMC_HEADER_GUARD
#include "numericType.h"
#include <boost/random/mersenne_twister.hpp>
#include "customMCMC.h"
namespace chordalGraph
{
	void armstrongMCMC(mcmcArgs& args);
}
#endif
