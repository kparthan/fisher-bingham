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
  //assert(eccentricity() < 1);
}

Kent::Kent(long double alpha, long double eta, long double psi, long double delta, 
           long double kappa, long double beta) : alpha(alpha), eta(eta), psi(psi), 
           delta(delta), kappa(kappa), beta(beta)
{
  //assert(eccentricity() < 1);
  Vector spherical(3,1);
  mu = Vector(3,0);
  // compute mu
  spherical[1] = alpha; spherical[2] = eta;
  spherical2cartesian(spherical,mu);
  // compute major_axis
  spherical[1] = psi;
  spherical[2] = delta;
  major_axis = Vector(3,0);
  spherical2cartesian(spherical,major_axis);
  // compute minor_axis
  minor_axis = crossProduct(mu,major_axis);
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
    alpha = source.alpha;
    eta = source.eta;
    psi = source.psi;
    delta = source.delta;
    kappa = source.kappa;
    beta = source.beta;
    constants = source.constants;
    computed = source.computed;
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

Kent::Constants Kent::getConstants()
{
  if (computed != SET) {
    computeExpectation();
  }
  return constants;
}

/*!
 *  Normalization constant
 */
long double Kent::computeLogNormalizationConstant(void)
{
  return computeSeriesSum(kappa,beta,0.5);
}

long double Kent::log_dc_dk(void)
{
  return computeSeriesSum(kappa,beta,1.5);
}

long double Kent::log_d2c_dk2(void)
{
  return computeSeriesSum(kappa,beta,2.5);
}

long double Kent::computeSeriesSum(long double k, long double b, long double d)
{
  long double ex = 2 * b/k; 
  long double log_ex = 2 * log(ex);
  //long double log_bessel = logModifiedBesselFirstKind(0,k);
  long double log_bessel_prev,log_bessel_current;
  long double log_f0,log_fj_current,log_fj_prev;
  long double current,log_diff,series_sum=1;
  int j = 0;
  long double gj;

  log_bessel_prev = log(cyl_bessel_i(d,k));
  log_f0 = lgamma<long double>(0.5) + log_bessel_prev;
  log_fj_prev = log_f0;
  while (1) {
    d += 2;
    log_bessel_current = log(cyl_bessel_i(d,k));
    gj = log_bessel_current - log_bessel_prev;
    gj += log_ex;
    gj += (log(j+0.5) - log(j+1));

    log_fj_current = log_fj_prev + gj;
    log_diff = log_fj_current - log_f0; 
    current = exp(log_diff); // < 1
    series_sum += current;
    if (current/series_sum <= TOLERANCE) {
      break;
    }
    j++;
    log_bessel_prev = log_bessel_current;
    log_fj_prev = log_fj_current;
  }
  long double ans = log(2*PI) + 0.5 * log(2.0/k);
  ans += log_f0;
  ans += log(series_sum);
  cout << "j: " << j << endl;
  //cout << "ans: " << ans << endl;
  return ans;
}

long double Kent::log_dc_db()
{
  return computeSeriesSum2(kappa,beta,2.5);
}

long double Kent::log_d2c_dkdb()
{
  return computeSeriesSum2(kappa,beta,3.5);
}

long double Kent::computeSeriesSum2(long double k, long double b, long double d)
{
  long double ex = 2 * b/k; 
  long double log_ex = 2 * log(ex);
  long double log_bessel_prev,log_bessel_current;
  long double log_f1,log_fj_current,log_fj_prev;
  long double current,log_diff,series_sum=1;
  int j = 1;
  long double gj;

  log_bessel_prev = log(cyl_bessel_i(d,k));
  log_f1 = lgamma<long double>(1.5) + log_ex + log_bessel_prev;
  log_fj_prev = log_f1;
  while (1) {
    d += 2;
    log_bessel_current = log(cyl_bessel_i(d,k));
    gj = log_bessel_current - log_bessel_prev;
    gj += log_ex;
    gj += (log(j+0.5) - log(j));

    log_fj_current = log_fj_prev + gj;
    log_diff = log_fj_current - log_f1; 
    current = exp(log_diff); // < 1
    series_sum += current;
    if (current/series_sum <= TOLERANCE) {
      break;
    }
    j++;
    log_bessel_prev = log_bessel_current;
    log_fj_prev = log_fj_current;
  }
  long double ans = log(2*PI) + 0.5 * log(2.0/k) + log(2/b);
  ans += log_f1;
  ans += log(series_sum);
  cout << "j: " << j << endl;
  return ans;
}

long double Kent::log_d2c_db2()
{
  long double b = beta;
  long double k = kappa;
  long double ex = 2 * b/k; 
  long double log_ex = 2 * log(ex);
  long double log_bessel_prev,log_bessel_current;
  long double log_f1,log_fj_current,log_fj_prev;
  long double current,log_diff,series_sum=1;
  long double d = 2.5;
  int j = 1;
  long double gj;

  log_bessel_prev = log(cyl_bessel_i(d,k));
  log_f1 = lgamma<long double>(1.5) + log_ex + log_bessel_prev;
  log_fj_prev = log_f1;
  while (1) {
    d += 2;
    log_bessel_current = log(cyl_bessel_i(d,k));
    gj = log_bessel_current - log_bessel_prev;
    gj += log_ex;
    gj += (2*log(2*j+1) - log(2*j) - log(2*j-1));

    log_fj_current = log_fj_prev + gj;
    log_diff = log_fj_current - log_f1; 
    current = exp(log_diff); // < 1
    series_sum += current;
    if (current/series_sum <= TOLERANCE) {
      break;
    }
    j++;
    log_bessel_prev = log_bessel_current;
    log_fj_prev = log_fj_current;
  }
  long double ans = log(2*PI) + 0.5 * log(2.0/k) + log(2) - 2*log(b);
  ans += log_f1;
  ans += log(series_sum);
  cout << "j: " << j << endl;
  return ans;
}

void Kent::computeConstants()
{
  constants.log_c = computeLogNormalizationConstant();
  constants.log_cb = log_dc_db();
  constants.log_cbb = log_d2c_db2();
  constants.log_ck = log_dc_dk();
  constants.log_ckk = log_d2c_dk2();
  constants.log_ckb = log_d2c_dkdb();

  long double tmp;
  tmp = constants.log_ck - constants.log_c;
  constants.ck_c = exp(tmp);

  tmp = constants.log_ckk - constants.log_c;
  constants.ckk_c = exp(tmp);

  tmp = constants.log_cb - constants.log_c;
  constants.cb_c = exp(tmp);

  tmp = constants.log_cbb - constants.log_c;
  constants.cbb_c = exp(tmp);

  tmp = constants.log_ckb - constants.log_c;
  constants.ckb_c = exp(tmp);

  constants.R = computeOrthogonalTransformation(mu,major_axis);
  constants.Rt = trans(constants.R);
}

/*!
 *  E[x x'] = R E[y y'] R'  where
 *  R is the rotation matrix from standard frame of reference 
 *  to the current orientation of axes
 */
void Kent::computeExpectation()
{
  computeConstants();

  constants.E_x = Vector(3,0);
  for (int i=0; i<3; i++) {
    constants.E_x[i] = mu[i] * constants.ck_c;
  }

  Matrix expectation_std = ZeroMatrix(3,3);
  expectation_std(0,0) = 0.5 * (1 - constants.ckk_c + constants.cb_c);  // lambda_1
  expectation_std(1,1) = 0.5 * (1 - constants.ckk_c - constants.cb_c);  // lambda_2
  expectation_std(2,2) = constants.ckk_c;  // lambda_3

  assert(expectation_std(0,0) > 0);
  assert(expectation_std(1,1) > 0);
  assert(expectation_std(2,2) > 0);
  cout << "E_std: " << expectation_std << endl;

  Matrix tmp1 = prod(constants.R,expectation_std);
  constants.E_xx = prod(tmp1,constants.Rt);

  computed = SET;
}

long double Kent::computeNegativeLogLikelihood(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeNegativeLogLikelihood(sample_mean,S);
}

long double Kent::computeNegativeLogLikelihood(Vector &sample_mean, Matrix &S)
{
  long double c1 = computeDotProduct(sample_mean,mu);
  Vector x = prod(S,major_axis);
  long double mj = computeDotProduct(major_axis,x);
  x = prod(S,minor_axis);
  long double mi = computeDotProduct(minor_axis,x);
  long double c2 = mj - mi;
  long double log_norm = computeLogNormalizationConstant();
  long double ans = log_norm - kappa * c1 - beta * c2;
  return ans;
}

long double Kent::computeLogPriorProbability()
{
}

long double Kent::computeLogFisherInformation()
{
  computeConstants();
  computeExpectation();
  long double log_det_axes = computeLogFisherAxes();
  long double log_det_kb = computeLogFisherScale();
  return log_det_axes + log_det_kb;
}

long double Kent::computeLogFisherAxes()
{
}

long double Kent::computeLogFisherScale()
{
  // E [d^2 L / d k^2]
  long double t1 = constants.ckk_c - (constants.ck_c * constants.ck_c);

  // E [d^2 L / d b^2]
  long double t2 = constants.cbb_c - (constants.cb_c * constants.cb_c);

  // E [d^2 L / dk db]
  long double t3 = constants.ckb_c - (constants.ck_c * constants.cb_c);

  long double det = t1 * t2 - t3;
  return log(det);
}

/*!
 *  Moment estimation (with +Z as the north pole)
 */
struct Estimates Kent::computeMomentEstimates(std::vector<Vector> &data)
{
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

  Optimize opt;
  opt.initialize(estimates.mean,estimates.major_axis,estimates.minor_axis,
                 estimates.kappa,estimates.beta);
  opt.computeMomentEstimates(sample_mean,S,estimates);
  return estimates;
}

/*!
 *  Max LH estimation
 */
struct Estimates Kent::computeMLEstimates(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeMLEstimates(sample_mean,S);
}

struct Estimates Kent::computeMLEstimates(Vector &sample_mean, Matrix &S)
{
  struct Estimates estimates = computeMomentEstimates(sample_mean,S);
  Optimize opt;
  opt.initialize(estimates.mean,estimates.major_axis,estimates.minor_axis,
                 estimates.kappa,estimates.beta);
  opt.computeMLEstimates(sample_mean,S,estimates);
  return estimates;
}

struct Estimates Kent::computeMMLEstimates(Vector &sample_mean, Matrix &S)
{
  struct Estimates estimates = computeMomentEstimates(sample_mean,S);
  Optimize opt;
  opt.initialize(estimates.mean,estimates.major_axis,estimates.minor_axis,
                 estimates.kappa,estimates.beta);
  opt.computeMMLEstimates(sample_mean,S,estimates);
  return estimates;
}

Vector Kent::Mean()
{
  return mu;
}

Vector Kent::MajorAxis()
{
  return major_axis;
}

Vector Kent::MinorAxis()
{
  return minor_axis;
}

long double Kent::Kappa()
{
  return kappa;
}

long double Kent::Beta()
{
  return beta;
}

long double Kent::computeKLDivergence(Kent &other)
{
  struct Constants constants1 = getConstants();
  struct Constants constants2 = other.getConstants();

  long double ans = constants2.log_c - constants1.log_c;
  
  long double kappa2 = other.Kappa();
  Vector mu2 = other.Mean();
  Vector kmu(3,0);
  for (int i=0; i<3; i++) {
    kmu[i] = kappa * mu[i] - kappa2 * mu2[i];
  }
  ans += computeDotProduct(kmu,constants1.E_x);

  long double tmp1,tmp2;
  tmp1 = prod_vMv(major_axis,constants1.E_xx);
  tmp2 = prod_vMv(minor_axis,constants1.E_xx);
  ans += beta * (tmp1 - tmp2);

  long double beta2 = other.Beta();
  Vector mj2 = other.MajorAxis();
  Vector mi2 = other.MinorAxis();
  tmp1 = prod_vMv(mj2,constants1.E_xx);
  tmp2 = prod_vMv(mi2,constants1.E_xx);
  ans -= beta2 * (tmp1 - tmp2);

  return ans;
}

