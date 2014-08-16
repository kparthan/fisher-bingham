#include "Kent.h"
#include "FB6.h"
#include "Optimize.h"

extern Vector XAXIS,YAXIS,ZAXIS;

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
  //assert(eccentricity() < 1);
}

/*!
 *  Constructor
 */
Kent::Kent(Vector &mu, Vector &major_axis, Vector &minor_axis,
          long double kappa, long double beta): mu(mu), 
          major_axis(major_axis), minor_axis(minor_axis), kappa(kappa), beta(beta)
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

std::vector<Vector> Kent::generate(int sample_size)
{
  cout << "\nGenerating from Kent:\n";
  cout << "mean: "; print(cout,mu,3); cout << endl;
  cout << "major: "; print(cout,major_axis,3); cout << endl;
  cout << "minor: "; print(cout,minor_axis,3); cout << endl;
  cout << "Kappa: " << kappa << "; Beta: " << beta 
       << "; sample size: " << sample_size << endl;
  std::vector<Vector> canonical_sample = generateCanonical(sample_size);
  Matrix transformation = computeOrthogonalTransformation(mu,major_axis);
  return transform(canonical_sample,transformation);
}

std::vector<Vector> Kent::generateCanonical(int sample_size)
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
 *  Moment estimation (with +Z as the north pole)
 */
struct Estimates Kent::computeMomentEstimates(std::vector<Vector> &data)
{
  int N = data.size();
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeMomentEstimates(sample_mean,S);
}

/*!
 *  Moment estimation (with +X as the north pole) 
 *  (similar to used in Kent (1982) paper)
 */
struct Estimates Kent::computeMomentEstimates(Vector &sample_mean, Matrix &S)
{ 
  // compute r1:
  long double r1 = norm(sample_mean);

  struct Estimates estimates;
  Vector spherical(3,0);
  cartesian2sphericalPoleXAxis(sample_mean,spherical);
  long double theta = spherical[1];
  long double phi = spherical[2];

  // rotation matrix to align north pole
  Matrix H(3,3);
  H(0,0) = cos(theta); H(0,1) = -sin(theta); H(0,2) = 0;
  H(1,0) = sin(theta)*cos(phi); H(1,1) = cos(theta)*cos(phi); H(1,2) = -sin(phi);
  H(2,0) = sin(theta)*sin(phi); H(2,1) = cos(theta)*sin(phi); H(2,2) = cos(phi);
  Matrix Ht = trans(H);
  Matrix tmp = prod(Ht,S);
  Matrix B = prod(tmp,H);
  long double ratio = 2 * B(1,2) / (B(1,1) - B(2,2));
  long double psi = 0.5 * atan(ratio);
  Matrix K = IdentityMatrix(3,3);
  K(1,1) = cos(psi); K(1,2) = -sin(psi);
  K(2,1) = -K(1,2); K(2,2) = K(1,1);
  Matrix HK = prod(H,K);

  // compute r2:
  long double t1 = (B(1,1)-B(2,2)) * (B(1,1)-B(2,2));
  long double t2 = 4 * B(1,2) * B(1,2);
  long double r2 = sqrt(t1+t2);

  estimates.mean = Vector(3,0);
  Vector axis1(3,0),axis2(3,0);
  for (int i=0; i<3; i++) {
    estimates.mean[i] = HK(i,0);
    axis1[i] = HK(i,1);
    axis2[i] = HK(i,2);
  }
  if (B(1,1) >= B(2,2)) {
    estimates.major_axis = axis1;
    estimates.minor_axis = axis2;
  } else {
    for (int j=0; j<3; j++) axis2[j] *= -1;
    estimates.major_axis = axis2;
    estimates.minor_axis = axis1;
  }

  // estimate kappa, beta
  long double f1 = 1/(2 - 2*r1 - r2);
  long double f2 = 1/(2 - 2*r1 + r2);
  estimates.kappa = f1 + f2;
  estimates.beta = 0.5 * (f1-f2);

  //Vector cross = crossProduct(estimates.major_axis,estimates.minor_axis);
  // submatrix 2 X 2
  /*Matrix B_sub(2,2);
  B_sub(0,0) = B(1,1); B_sub(0,1) = B(1,2);
  B_sub(1,0) = B(2,1); B_sub(1,1) = B(2,2);
  // ... find eigen values of submatrix by solving a simple quadratic equation
  Vector eigen_vals(2,0);
  Matrix eigen_vecs = IdentityMatrix(2,2);
  eigenDecomposition(B_sub,eigen_vals,eigen_vecs);*/
  //cout << "theta: " << theta*180/PI << "; phi: " << phi*180/PI << endl;
  //cout << "B: " << B << endl;
  //cout << "psi (degrees): " << psi * 180/PI << endl;
  //cout << "H: " << H << endl;
  //cout << "K: " << K << endl;
  //cout << "HK: " << HK << endl;
  //cout << "cross: "; print(cout,cross,0); cout << endl;
  //cout << "r1: " << r1 << endl;
  //cout << "r2: " << r2 << endl;

  Optimize opt(estimates.mean,estimates.major_axis,estimates.minor_axis);
  opt.computeMomentEstimates(sample_mean,S);
  
  return estimates;
}

/*struct Estimates Kent::computeMomentEstimates(
  int sample_size, 
  Vector &sample_mean,
  Matrix &T
) { 
  struct Estimates moment_estimates;
  Matrix Ht = align_zaxis_with_vector(sample_mean);
  Matrix H = trans(Ht);
  Matrix tmp = prod(Ht,T);
  Matrix B = prod(tmp,H);
  cout << "B: " << B << endl;
  long double ratio = 2 * B(0,1) / (B(0,0) - B(1,1));
  long double psi = 0.5 * atan(ratio);
  cout << "psi (degrees): " << psi * 180/PI << endl;
  Matrix K = IdentityMatrix(3,3);
  K(0,0) = cos(psi); K(0,1) = -sin(psi);
  K(1,0) = -K(0,1); K(1,1) = K(0,0);
  cout << "Ht: " << Ht << endl;
  cout << "H: " << H << endl;
  cout << "K: " << K << endl;
  Matrix est = prod(H,K);
  cout << "HK: " << est << endl;
  moment_estimates.mean = Vector(3,0);
  moment_estimates.major_axis = Vector(3,0);
  moment_estimates.minor_axis = Vector(3,0);
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
  long double t2 = 4 * B(0,1) * B(0,1);
  long double r2 = sqrt(t1+t2);
  cout << "r2: " << r2 << endl;

  return moment_estimates;
}*/

