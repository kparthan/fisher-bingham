#include "FB4.h"
#include "vMF.h"
#include "Normal.h"

/*!
 *  Null constructor
 */
FB4::FB4() : kappa(1), gamma(1) 
{}

/*!
 *  Constructor
 */
FB4::FB4(long double kappa, long double gamma) : kappa(kappa), gamma(gamma)
{}

/*!
 *  Copy the elements
 */
FB4 FB4::operator=(const FB4 &source)
{
  if (this != &source) {
    kappa = source.kappa;
    gamma = source.gamma;
  }
  return *this;
}

/*!
 *  Determines the best envelope function
 */
int FB4::determine_best_case()
{
  // calculate acceptance ratios to determine the best envelope 
  long double a0,a1,a2;
  Normal normal(0,1);
  long double x1 = (kappa - 2*gamma) / sqrt(-2*gamma);
  long double x2 = (kappa + 2*gamma) / sqrt(-2*gamma);
  long double cdf1 = normal.cumulativeDensity(x1);
  long double cdf2 = normal.cumulativeDensity(x2);
  a0 = cdf1 - cdf2;

  // a1/a0
  long double abs_gamma = fabs(gamma);
  long double tmp = ((kappa * kappa) / (4 * abs_gamma)) + (0.5 * log(PI/gamma));
  long double log_tmp = log(kappa) + tmp - log(2*sinhl(kappa));
  long double ratio = exp(log_tmp);
  a1 = ratio * a0;

  // a2/a0
  log_tmp = log(kappa+2*gamma) + tmp + gamma - log(2*sinhl(kappa+2*gamma));
  ratio = exp(log_tmp);
  a2 = ratio * a0;

  if (a0 > a1 && a0 > a2) {
    return 1;
  } else if (a1 > a0 && a1 > a2) {
    return 2;
  } else {
    return 3;
  }
}

std::vector<std::vector<long double> > FB4::generate(int sample_size)
{
  if (gamma < 0) {
    int best = determine_best_case();
    return generate_FB4_minus(best,sample_size);
  } else if (gamma > 0) {
  } else if (gamma == 0) {  // vMF
    std::vector<long double> mu(3,0); mu[2] = 1;
    vMF vmf(mu,kappa);
    return vmf.generate(sample_size);
  }
}

std::vector<std::vector<long double> > FB4::generate_FB4_minus(int best, int sample_size)
{
  vector<long double> u(sample_size,0);
  switch(best) {
    case 1:
      break;

    case 2:
      break;

    case 3:
      break;

    default:
      cout << "Error: best case = " << best << endl;
      exit(1);
  }
}

