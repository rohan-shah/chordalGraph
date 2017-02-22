#ifndef POSTERIOR_INFERENCE_RPACKAGE_HEADER_GUARD
#define POSTERIOR_INFERENCE_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
SEXP posteriorInference(SEXP outerProductsSum, SEXP delta, SEXP dimension, SEXP dataPoints, SEXP psi, SEXP exactCounts, SEXP burnIn, SEXP runSize);
#endif
