#ifndef NUMERIC_TYPE_HEADER_GUARD
#define NUMERIC_TYPE_HEADER_GUARD
#include <boost/multiprecision/mpfr.hpp>
typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<50> > mpfr_class;
typedef boost::multiprecision::mpz_int mpz_class;
typedef mpfr_class numericType;
#endif
