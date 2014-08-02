#include "FB6.h"

/*!
 *  Null constructor
 */
FB6::FB6() : kappa(1), beta(1), gamma(1) 
{}

/*!
 *  Constructor
 */
FB6::FB6(long double kappa, long double beta, long double gamma) : 
        kappa(kappa), beta(beta), gamma(gamma)
{}

/*!
 *  Copy the elements
 */
FB6 FB6::operator=(const FB6 &source)
{
  if (this != &source) {
    kappa = source.kappa;
    beta = source.beta;
    gamma = source.gamma;
  }
  return *this;
}

std::vector<std::vector<long double> > FB6::generate(int sample_size)
{
}
