#ifndef IS_CHORDAL_ITERATIVE_HEADER_GUARD
#define IS_CHORDAL_ITERATIVE_HEADER_GUARD
#include <Rcpp.h>
SEXP isChordalIterative_graphAM(SEXP graph_sexp);
SEXP isChordalIterative_graphNEL(SEXP graph_sexp);
SEXP isChordalIterative_igraph(SEXP graph_sexp);
#endif