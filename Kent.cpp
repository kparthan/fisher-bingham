#include "Kent.h"
#include "FB6.h"

extern std::vector<long double> XAXIS,YAXIS,ZAXIS;

/*!
 *  Null constructor
 */
Kent::Kent() : kappa(1), beta(0)
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
  assert(eccentricity() < 1);
}

/*!
 *  Constructor
 */
Kent::Kent(std::vector<long double> &mu, std::vector<long double> &major_axis, 
         std::vector<long double> &minor_axis, long double kappa, long double beta)
         : mu(mu), major_axis(major_axis), minor_axis(minor_axis), kappa(kappa), beta(beta)
{
  assert(eccentricity() < 1);
}

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
  //long double log_bessel = logModifiedBesselFirstKind(0,k);
  long double log_bessel = log(cyl_bessel_i(0,k));
  long double log_f0 = lgamma<long double>(0.5) + log_bessel;
  long double log_fj,prev=1,current,diff,delta,series_sum=1;
  int j = 1;
  while (1) {
    log_fj = lgamma<long double>(j+0.5) - lgamma<long double>(j+1);
    log_fj += j * log_ex;
    //log_bessel = logModifiedBesselFirstKind(2*j+0.5,k);
    log_bessel = log(cyl_bessel_i(2*j+0.5,k));
    log_fj += log_bessel;
    diff = log_fj - log_f0; 
    current = exp(diff); // < 1
    series_sum += current;
    if (current/series_sum <= ZERO) {
      break;
    }
    prev = current;
    j++;
  }
  long double log_norm = log(2*PI) + 0.5 * log(2.0/k);
  log_norm += log_f0;
  log_norm += log(series_sum);
  //cout << "log_norm: " << log_norm << endl;
  return log_norm;
}

/*!
 *  Moment estimation
 */
struct Estimates Kent::computeMomentEstimates(std::vector<std::vector<long double> > &data)
{
  int N = data.size();
  std::vector<long double> sample_mean = computeVectorSum(data);
  matrix<long double> S = computeDispersionMatrix(data);
  return computeMomentEstimates(N,sample_mean,S);
}

struct Estimates Kent::computeMomentEstimates(
  int sample_size, 
  std::vector<long double> &sample_mean,
  matrix<long double> &T
) { 
  struct Estimates moment_estimates;
  matrix<long double> Ht = align_zaxis_with_vector(sample_mean);
  matrix<long double> H = trans(Ht);
  matrix<long double> tmp = prod(Ht,T);
  matrix<long double> B = prod(tmp,H);
  cout << "B: " << B << endl;
  long double ratio = 2 * B(0,1) / (B(0,0) - B(1,1));
  long double psi = 0.5 * atan(ratio);
  cout << "psi (degrees): " << psi * 180/PI << endl;
  matrix<long double> K = identity_matrix<long double>(3,3);
  K(0,0) = cos(psi); K(0,1) = -sin(psi);
  K(1,0) = -K(0,1); K(1,1) = K(0,0);
  cout << "H: " << H << endl;
  cout << "K: " << K << endl;
  matrix<long double> est = prod(H,K);
  cout << "HK: " << est << endl;
  moment_estimates.mean = std::vector<long double>(3,0);
  moment_estimates.major_axis = std::vector<long double>(3,0);
  moment_estimates.minor_axis = std::vector<long double>(3,0);
  for (int i=0; i<3; i++) {
    moment_estimates.mean[i] = est(i,0);
    moment_estimates.major_axis[i] = est(i,1);
    moment_estimates.minor_axis[i] = est(i,2);
  }
  // compute r1:
  long double r1 = norm(sample_mean);
  cout << "r1: " << r1 << endl;

  // compute r2:
  long double t1 = (B(0,0)-B(1,1)) * (B(0,0)-B(1,1));
  long double t2 = 4 * B(0,1) * B(1,0);
  long double r2 = sqrt(t1+t2);
  cout << "r2: " << r2 << endl;

  return moment_estimates;
}

/*struct Estimates Kent::computeMomentEstimates(
  int sample_size, 
  std::vector<long double> &sample_mean,
  matrix<long double> &T
) { 
  struct Estimates moment_estimates;
  std::vector<long double> spherical(3,0);
  cartesian2spherical(sample_mean,spherical);
  //long double theta = spherical[1];
  //long double phi = spherical[2];
  long double theta = 1.485;
  long double phi = 3;
  matrix<long double> H(3,3);
  H(0,0) = cos(theta); H(0,1) = -sin(theta); H(0,2) = 0;
  H(1,0) = sin(theta)*cos(phi); H(1,1) = cos(theta)*cos(phi); H(1,2) = -sin(phi);
  H(2,0) = sin(theta)*sin(phi); H(2,1) = cos(theta)*sin(phi); H(2,2) = cos(phi);
  matrix<long double> Ht = trans(H);
  matrix<long double> tmp = prod(Ht,T);
  matrix<long double> B = prod(tmp,H);
  cout << "B: " << B << endl;
  long double ratio = 2 * B(1,2) / (B(1,1) - B(2,2));
  long double psi = 0.5 * atan(ratio);
  cout << "psi (degrees): " << psi * 180/PI << endl;
  matrix<long double> K = identity_matrix<long double>(3,3);
  K(1,1) = cos(psi); K(1,2) = -sin(psi);
  K(2,1) = -K(1,2); K(2,2) = K(1,1);
  cout << "H: " << H << endl;
  cout << "K: " << K << endl;
  matrix<long double> est = prod(H,K);
  cout << "HK: " << est << endl;
  moment_estimates.mean = std::vector<long double>(3,0);
  moment_estimates.major_axis = std::vector<long double>(3,0);
  moment_estimates.minor_axis = std::vector<long double>(3,0);
  for (int i=0; i<3; i++) {
    moment_estimates.mean[i] = est(i,0);
    moment_estimates.major_axis[i] = est(i,1);
    moment_estimates.minor_axis[i] = est(i,2);
  }
  // compute r1:
  long double r1 = norm(sample_mean);
  cout << "r1: " << r1 << endl;

  // compute r2:
  long double t1 = (B(1,1)-B(2,2)) * (B(1,1)-B(2,2));
  long double t2 = 4 * B(1,2) * B(1,2);
  long double r2 = sqrt(t1+t2);
  cout << "r2: " << r2 << endl;

  return moment_estimates;
}*/

