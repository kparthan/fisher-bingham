#include "FB6.h"

extern std::vector<long double> XAXIS,YAXIS,ZAXIS;

/*!
 *  Null constructor
 */
FB6::FB6() : kappa(1), beta(1), gamma(1) 
{
  mu = ZAXIS;
  major_axis = XAXIS;
  minor_axis = YAXIS;
}

/*!
 *  Constructor
 */
FB6::FB6(long double kappa, long double beta, long double gamma) : 
        kappa(kappa), beta(beta), gamma(gamma)
{
  mu = ZAXIS;
  major_axis = XAXIS;
  minor_axis = YAXIS;
}

/*!
 *  Constructor
 */
FB6::FB6(std::vector<long double> &mu, std::vector<long double> &major_axis, 
         std::vector<long double> &minor_axis, long double kappa, long double beta,
         long double gamma) : mu(mu), major_axis(major_axis), minor_axis(minor_axis), 
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
    mu = source.mu;
    major_axis = source.major_axis;
    minor_axis = source.minor_axis;
  }
  return *this;
}

std::vector<std::vector<long double> > FB6::generate(int sample_size)
{
}
