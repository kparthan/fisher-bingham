#include "vMC.h"
#include "Support.h"
#include "Normal.h"

/*!
 *  \brief This is a constructor module
 */
vMC::vMC()
{
  mu = Vector(2,0); mu[0] = 1;
  kappa = 1;
}

/*!
 *  \brief constructor function which sets the value of mean and 
 *  kappa of the distribution
 *  \param mu a reference to a vector<long double>
 *  \param kappa a long double
 */
vMC::vMC(Vector &mu, long double kappa) : mu(mu), kappa(kappa)
{}

vMC::vMC(long double kappa) : kappa(kappa)
{
  mu = Vector(2,0); mu[0] = 1;
}

/*!
 *  \brief This function assigns a source vMC distribution.
 *  \param source a reference to a vMC
 */
vMC vMC::operator=(const vMC &source)
{
  if (this != &source) {
    mu = source.mu;
    kappa = source.kappa;
  }
  return *this;
}

/*!
 *  \brief This function returns the mean of the distribution
 *  \return the mean of the distribution
 */
Vector vMC::Mean(void)
{
	return mu;
}

/*!
 *  \brief This function returns the kappa of the distribution
 *  \return the kappa of the distribution
 */
long double vMC::Kappa(void)
{
	return kappa;
}

void vMC::generateCanonical(std::vector<Vector> &canonical_sample, int sample_size)
{
  long double tmp = 1 + (4 * kappa * kappa);
  long double a = 1 + sqrt(tmp);

  tmp = a - sqrt(2*a);
  long double b = tmp / (2*kappa);

  tmp = 1 + (b*b);
  long double r = tmp / (2*b);

  long double u1,u2,u3,z,f,c,check1,check2,angle;
  Vector thetas;
  for (int i=0; i<sample_size; i++) {
    repeat:
    u1 = uniform_random();
    z = cos(PI*u1);
    f = (1+r*z) / (r+z);
    c = kappa * (r-f);
    u2 = uniform_random();
    check1 = (c * (2-c)) - u2;
    if (check1 > 0) {
      goto accept;
    } else if (check1 <= 0) {
      check2 = log(c) - log(u2) + 1 - c;
      if (check2 < 0) goto repeat;
      else if (check2 >= 0) goto accept;
    }
    accept:
    u3 = uniform_random();
    angle = (PI/2) + (sign(u3-0.5)*acos(f));
    thetas.push_back(angle);
  }
  //assert(thetas.size() == sample_size);
  canonical_sample.clear();
  Vector v(2,0);
  for (int i=0; i<sample_size; i++) {
    v[0] = cos(thetas[i]);
    v[1] = sin(thetas[i]);
    canonical_sample.push_back(v);
  }
}

