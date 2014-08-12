#include "Kent.h"
#include "FB6.h"
#include "Support.h"

extern std::vector<long double> XAXIS,YAXIS,ZAXIS;

/*!
 *  Null constructor
 */
Kent::Kent() : kappa(1), beta(1)
{
  mu = ZAXIS;
  major_axis = XAXIS;
  minor_axis = YAXIS;
}

/*!
 *  Constructor
 */
Kent::Kent(long double kappa, long double beta) : kappa(kappa), beta(beta)
{
  mu = ZAXIS;
  major_axis = XAXIS;
  minor_axis = YAXIS;
}

/*!
 *  Constructor
 */
Kent::Kent(std::vector<long double> &mu, std::vector<long double> &major_axis, 
         std::vector<long double> &minor_axis, long double kappa, long double beta)
         : mu(mu), major_axis(major_axis), minor_axis(minor_axis), kappa(kappa), beta(beta)
{}

/*!
 *  Copy the elements
 */
Kent Kent::operator=(const Kent &source)
{
  if (this != &source) {
    mu = source.mu;
    major_axis = source.major_axis;
    minor_axis = source.minor_axis;
    kappa = source.kappa;
    beta = source.beta;
  }
  return *this;
}

std::vector<std::vector<long double> > Kent::generate(int sample_size)
{
  cout << "\nGenerating from Kent:\n";
  cout << "mean: "; print(cout,mu,3); cout << endl;
  cout << "major: "; print(cout,major_axis,3); cout << endl;
  cout << "minor: "; print(cout,minor_axis,3); cout << endl;
  cout << "Kappa: " << kappa << "; Beta: " << beta 
       << "; sample size: " << sample_size << endl;
  std::vector<std::vector<long double> > canonical_sample = generateCanonical(sample_size);
  matrix<long double> transformation = computeOrthogonalTransformation(mu,major_axis);
  return transform(canonical_sample,transformation);
}

std::vector<std::vector<long double> > Kent::generateCanonical(int sample_size)
{
  FB6 fb6(kappa,beta,0);
  return fb6.generateCanonical(sample_size);
}

/*!
 *  e = 2 beta / kappa
 */
long double Kent::eccentricity()
{
  return 2 * beta / kappa;
}

/*!
 *  Normalization constant
 */
long double Kent::computeLogNormalizationConstant(void)
{
  return computeLogNormalizationConstant(kappa,beta);
}

/*!
 *  Normalization constant
 */
long double Kent::computeLogNormalizationConstant(long double k, long double b)
{
  long double ex = 2 * b/k; 
  long double log_ex = 2 * log(ex);
  long double log_bessel = logModifiedBesselFirstKind(0,k);
  long double log_f0 = lgamma<long double>(0.5) + log_bessel;
  long double log_fj,prev=1,current,diff,delta,series_sum=1;
  int j = 1;
  while (1) {
    log_fj = lgamma<long double>(j+0.5) - lgamma<long double>(j+1);
    log_fj += j * log_ex;
    log_bessel = logModifiedBesselFirstKind(2*j+0.5,k);
    log_fj += log_bessel;
    diff = log_fj - log_f0; 
    current = exp(diff); // < 1
    delta = (prev-current)/prev;
    if (delta < TOLERANCE) {
      break;
    }
    series_sum += current;
    prev = current;
    j++;
  }
  long double log_norm = log(2*PI) + 0.5 * log(2.0/k);
  log_norm += log_f0;
  log_norm += log(series_sum);
  cout << "#j: " << j << "; log_norm: " << log_norm << endl;
  return log_norm;
}

