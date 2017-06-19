#ifndef CUSTOM_SYMMETRIC_POSTERIOR_INFERENCE_RPACKAGE_HEADER_GUARD
#define CUSTOM_SYMMETRIC_POSTERIOR_INFERENCE_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
SEXP customSymmetricPosteriorInference(SEXP outerProductsSum, SEXP delta, SEXP dimension, SEXP dataPoints, SEXP psi, SEXP exactCounts, SEXP burnIn, SEXP runSize, SEXP seed, SEXP uniqueGraphsLimit);
#endif
