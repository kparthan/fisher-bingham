#include "Kent.h"
#include "FB6.h"
#include "Optimize.h"
#include "Support.h"

extern Vector XAXIS,YAXIS,ZAXIS;

/*!
 *  Null constructor
 */
Kent::Kent() : kappa(1), beta(0)
{
  mu = ZAXIS;
  major_axis = XAXIS;
  minor_axis = YAXIS;
  computed = UNSET;
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
  computed = UNSET;
}

/*!
 *  Constructor
 */
Kent::Kent(Vector &mu, Vector &major_axis, Vector &minor_axis,
          long double kappa, long double beta): mu(mu), 
          major_axis(major_axis), minor_axis(minor_axis), kappa(kappa), beta(beta)
{
  //assert(eccentricity() < 1);
  Vector spherical(3,0);
  cartesian2spherical(mu,spherical);
  alpha = spherical[1];
  eta = spherical[2];
  Matrix r = align_vector_with_zaxis(mu);
  Vector mj = prod(r,major_axis);
  cartesian2spherical(mj,spherical);
  psi = spherical[2];
  computed = UNSET;
}

Kent::Kent(long double psi, long double alpha, long double eta, 
           long double kappa, long double beta) : psi(psi), alpha(alpha), eta(eta), 
           kappa(kappa), beta(beta)
{
  Matrix r = computeOrthogonalTransformation(psi,alpha,eta);
  mu = Vector(3,0); major_axis = mu; minor_axis = mu;
  for (int i=0; i<3; i++) {
    major_axis[i] = r(i,0);
    minor_axis[i] = r(i,1);
    mu[i] = r(i,2);
  }
  computed = UNSET;
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
    psi = source.psi;
    alpha = source.alpha;
    eta = source.eta;
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
  long double x = computeSeriesSum(kappa,beta,2.5); // log a
  long double y = constants.log_ck; // log b
  long double tmp1 = x - y;
  long double tmp2 = exp(tmp1) + (1/kappa);
  long double ans = y + log(tmp2);
  return ans;
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
  //cout << "j: " << j << endl;
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
  //cout << "j: " << j << endl;
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
  //cout << "j: " << j << endl;
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
  expectation_std(0,0) = 0.5 * (1 - constants.ckk_c + constants.cb_c);  // lambda_2
  constants.lambda2 = expectation_std(0,0);
  expectation_std(1,1) = 0.5 * (1 - constants.ckk_c - constants.cb_c);  // lambda_3
  constants.lambda3 = expectation_std(1,1);
  expectation_std(2,2) = constants.ckk_c;  // lambda_1
  constants.lambda1 = expectation_std(2,2);

  /*assert(expectation_std(0,0) > 0);
  assert(expectation_std(1,1) > 0);
  assert(expectation_std(2,2) > 0);*/
  //cout << "E_std: " << expectation_std << endl;

  Matrix tmp1 = prod(constants.R,expectation_std);
  constants.E_xx = prod(tmp1,constants.Rt);

  computed = SET;
}

long double Kent::log_density(Vector &x)
{
  long double c1 = computeDotProduct(mu,x);

  Matrix xx = outer_prod(x,x);
  long double mj = prod_xMy(major_axis,xx,major_axis);
  long double mi = prod_xMy(minor_axis,xx,minor_axis);
  long double c2 = mj - mi;

  long double log_norm; 
  if (computed == UNSET) {
    log_norm = computeLogNormalizationConstant();
  } else {
    log_norm = constants.log_c;
  }

  long double ans = -log_norm + kappa * c1 + beta * c2;
  return ans;
}

long double Kent::computeNegativeLogLikelihood(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeNegativeLogLikelihood(sample_mean,S,data.size());
}

long double Kent::computeNegativeLogLikelihood(Vector &sample_mean, Matrix &S, long double N)
{
  long double c1 = computeDotProduct(sample_mean,mu);

  long double mj = prod_xMy(major_axis,S,major_axis);
  long double mi = prod_xMy(minor_axis,S,minor_axis);
  long double c2 = mj - mi;

  long double log_norm;
  if (computed == UNSET) {
    log_norm = computeLogNormalizationConstant();
  } else {
    log_norm = constants.log_c;
  }

  long double ans = N * log_norm - kappa * c1 - beta * c2;
  return ans;
}

long double Kent::computeLogPriorProbability()
{
  long double log_prior_axes = computeLogPriorAxes();
  long double log_prior_scale = computeLogPriorScale();
  long double log_joint_prior = log_prior_axes + log_prior_scale;
  assert(!boost::math::isnan(log_joint_prior));
  return log_joint_prior;
}

long double Kent::computeLogPriorAxes()
{
  long double angle = alpha;
  /*while (angle < 0) {
    angle += PI;
  }
  while (angle > PI) {
    angle -= PI;
  }*/

  if (angle < TOLERANCE) angle = TOLERANCE;
  
  long double log_prior = 0;
  log_prior += -log(4) - 3*log(PI);
  log_prior += log(sin(angle));
  return log_prior;
}

long double Kent::computeLogPriorScale()
{
  long double log_prior = 0;
  /*log_prior += log(-log(TOLERANCE)) - log(1-TOLERANCE);
  log_prior += 2 * log(kappa); 
  log_prior -= 2 * log(1+kappa*kappa);
  log_prior += log(8/PI);*/
  /*log_prior += log(4/PI);
  log_prior += 2 * log(kappa);
  log_prior -= 2 * log(1+kappa*kappa);
  log_prior += 2 * log(beta);
  log_prior -= 2 * log(1+beta*beta);
  long double cinv = 0.5 * atan(kappa/2.0);
  long double tmp = kappa + (4/kappa);
  cinv -= 1/tmp;
  log_prior -= log(cinv);*/
  long double ex = (2 * beta) / kappa;
  log_prior -= (log(PI-2) - log(8));
  log_prior += 2 * log(ex);
  log_prior -= 2 * log(1+ex*ex);
  return log_prior;
}

long double Kent::computeLogFisherInformation()
{
  if (computed == UNSET) {
    computeExpectation();
  }
  long double log_det_axes = computeLogFisherAxes();
  long double log_det_kb = computeLogFisherScale();
  return log_det_axes + log_det_kb; 
}

long double Kent::computeLogFisherInformation(long double N)
{
  long double log_fisher = computeLogFisherInformation(); 
  return log_fisher + 5 * log(N);
}

long double Kent::computeLogFisherAxes()
{
  long double sin_psi = sin(psi);
  long double sinsq_psi = sin_psi * sin_psi;
  long double cos_psi = cos(psi);
  long double cossq_psi = cos_psi * cos_psi; 
  long double sin_alpha = sin(alpha);
  long double sinsq_alpha = sin_alpha * sin_alpha;
  long double cos_alpha = cos(alpha);
  long double cossq_alpha = cos_alpha * cos_alpha; 

  constants.fisher_axes = ZeroMatrix(3,3);
  
  long double ans,tmp1,tmp2,tmp3;
  // E [d^2 L / d a^2]
  ans = kappa * constants.ck_c;
  tmp1 = (constants.lambda1 - constants.lambda3) * sinsq_psi;
  tmp2 = (constants.lambda1 - constants.lambda2) * cossq_psi;
  ans += 2 * beta * (tmp1 - tmp2);
  constants.fisher_axes(0,0) = ans;

  // E [d^2 L / d n^2]
  ans = kappa * constants.ck_c * sinsq_alpha;
  tmp1 = cossq_alpha * cossq_psi + sinsq_psi;
  tmp1 *= constants.lambda2;
  tmp2 = cossq_alpha * sinsq_psi + cossq_psi;
  tmp2 *= constants.lambda3;
  tmp3 = tmp1 - tmp2;
  tmp1 = constants.cb_c * cossq_alpha;
  tmp2 = constants.lambda1 * sinsq_alpha * (cossq_psi - sinsq_psi);
  tmp3 += (tmp1 + tmp2);
  ans += 2 * beta * tmp3;
  constants.fisher_axes(1,1) = ans;
  
  // E [d^2 L / d s^2]
  ans = 4 * beta * constants.cb_c;
  constants.fisher_axes(2,2) = ans;

  // E [d^2 L / dadn]
  tmp1 = 1 - (3 * constants.lambda1);
  ans = tmp1 * 2 * beta * sin_alpha * sin_psi * cos_psi;
  constants.fisher_axes(0,1) = ans;

  // E [d^2 L / dads]
  constants.fisher_axes(0,2) = 0;

  // E [d^2 L / dnds]
  ans = constants.fisher_axes(2,2) * cos_alpha;
  constants.fisher_axes(1,2) = ans;

  constants.fisher_axes(1,0) = constants.fisher_axes(0,1);
  constants.fisher_axes(2,0) = constants.fisher_axes(0,2);
  constants.fisher_axes(2,1) = constants.fisher_axes(1,2);

  long double det = determinant(constants.fisher_axes);
  //cout << "det: " << det << endl; //exit(1);
  return log(det);
}

long double Kent::computeLogFisherScale()
{
  // E [d^2 L / d k^2]
  long double t1 = constants.ckk_c - (constants.ck_c * constants.ck_c);

  // E [d^2 L / d b^2]
  long double t2 = constants.cbb_c - (constants.cb_c * constants.cb_c);

  // E [d^2 L / dk db]
  long double t3 = constants.ckb_c - (constants.ck_c * constants.cb_c);

  long double det = t1 * t2 - t3 * t3;
  //cout << "det: " << det << endl;
  return log(det);
}

void Kent::computeAllEstimators(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  computeAllEstimators(sample_mean,S,data.size());
}

void Kent::computeAllEstimators(Vector &sample_mean, Matrix &S, long double N)
{
  string type = "MOMENT";
  struct Estimates moment_est = computeMomentEstimates(sample_mean,S,N);
  print(type,moment_est);
  cout << "msglen: " << computeMessageLength(moment_est,sample_mean,S,N) << endl;
  cout << "KL-divergence: " << computeKLDivergence(moment_est) << endl << endl;

  type = "MLE";
  struct Estimates ml_est = moment_est;
  Optimize opt1(type);
  opt1.initialize(N,ml_est.mean,ml_est.major_axis,ml_est.minor_axis,
                  ml_est.kappa,ml_est.beta);
  opt1.computeEstimates(sample_mean,S,ml_est);
  print(type,ml_est);
  cout << "msglen: " << computeMessageLength(ml_est,sample_mean,S,N) << endl;
  cout << "KL-divergence: " << computeKLDivergence(ml_est) << endl << endl;

  type = "MAP";
  struct Estimates map_est = moment_est;
  Optimize opt2(type);
  opt2.initialize(N,map_est.mean,map_est.major_axis,map_est.minor_axis,
                  map_est.kappa,map_est.beta);
  opt2.computeEstimates(sample_mean,S,map_est);
  print(type,map_est);
  cout << "msglen: " << computeMessageLength(map_est,sample_mean,S,N) << endl;
  cout << "KL-divergence: " << computeKLDivergence(map_est) << endl << endl;

  type = "MML_2";
  struct Estimates mml_est1 = map_est;
  Optimize opt3(type);
  opt3.initialize(N,mml_est1.mean,mml_est1.major_axis,mml_est1.minor_axis,
                  mml_est1.kappa,mml_est1.beta);
  opt3.computeEstimates(sample_mean,S,mml_est1);
  print(type,mml_est1);
  cout << "msglen: " << computeMessageLength(mml_est1,sample_mean,S,N) << endl;
  cout << "KL-divergence: " << computeKLDivergence(mml_est1) << endl << endl;

  type = "MML_5";
  struct Estimates mml_est2 = map_est;
  Optimize opt4(type);
  opt4.initialize(N,mml_est2.mean,mml_est2.major_axis,mml_est2.minor_axis,
                  mml_est2.kappa,mml_est2.beta);
  opt4.computeEstimates(sample_mean,S,mml_est2);
  print(type,mml_est2);
  cout << "msglen: " << computeMessageLength(mml_est2,sample_mean,S,N) << endl;
  cout << "KL-divergence: " << computeKLDivergence(mml_est2) << endl << endl;
}

void Kent::computeAllEstimators(
  std::vector<Vector> &data, 
  std::vector<struct Estimates> &all_estimates
) {
  int N = data.size();
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);

  all_estimates.clear();

  string type = "MOMENT";
  struct Estimates moment_est = computeMomentEstimates(sample_mean,S,N);
  print(type,moment_est);
  all_estimates.push_back(moment_est);

  type = "MLE";
  struct Estimates ml_est = moment_est;
  Optimize opt1(type);
  opt1.initialize(N,ml_est.mean,ml_est.major_axis,ml_est.minor_axis,
                  ml_est.kappa,ml_est.beta);
  opt1.computeEstimates(sample_mean,S,ml_est);
  print(type,ml_est);
  all_estimates.push_back(ml_est);

  type = "MAP";
  struct Estimates map_est = moment_est;
  Optimize opt2(type);
  opt2.initialize(N,map_est.mean,map_est.major_axis,map_est.minor_axis,
                  map_est.kappa,map_est.beta);
  opt2.computeEstimates(sample_mean,S,map_est);
  print(type,map_est);
  all_estimates.push_back(map_est);

  type = "MML_2";
  struct Estimates mml_est = map_est;
  Optimize opt3(type);
  opt3.initialize(N,mml_est.mean,mml_est.major_axis,mml_est.minor_axis,
                  mml_est.kappa,mml_est.beta);
  opt3.computeEstimates(sample_mean,S,mml_est);
  all_estimates.push_back(mml_est);
  print(type,mml_est);

  type = "MML_5";
  struct Estimates mml_est2 = map_est;
  Optimize opt4(type);
  opt4.initialize(N,mml_est2.mean,mml_est2.major_axis,mml_est2.minor_axis,
                  mml_est2.kappa,mml_est2.beta);
  opt4.computeEstimates(sample_mean,S,mml_est2);
  print(type,mml_est2);
  all_estimates.push_back(mml_est2);
}

/*!
 *  Moment estimation (with +Z as the north pole)
 */
struct Estimates Kent::computeMomentEstimates(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeMomentEstimates(sample_mean,S,data.size());
}

/*!
 *  Moment estimation (with +X as the north pole) 
 *  (similar to used in Kent (1982) paper)
 *  (sample_mean1 and S1 are \sum_x and \sum_ x x')
 */
struct Estimates Kent::computeMomentEstimates(Vector &sample_mean1, Matrix &S1, long double N)
{
  Vector sample_mean = sample_mean1;
  Matrix S = S1; 
  for (int i=0; i<3; i++) {
    sample_mean[i] /= N;
    for (int j=0; j<3; j++) {
      S(i,j) /= N;
    }
  }
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

  Optimize opt("MOMENT");
  opt.initialize(N,estimates.mean,estimates.major_axis,estimates.minor_axis,
                 estimates.kappa,estimates.beta);
  opt.computeEstimates(sample_mean1,S1,estimates);
  return estimates;
}

/*!
 *  Max LH estimation
 */
struct Estimates Kent::computeMLEstimates(std::vector<Vector> &data, string type)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeMLEstimates(sample_mean,S,data.size(),type);
}

struct Estimates Kent::computeMLEstimates(Vector &sample_mean, Matrix &S, long double N, string type)
{
  struct Estimates estimates = computeMomentEstimates(sample_mean,S,N);
  Optimize opt(type);
  opt.initialize(N,estimates.mean,estimates.major_axis,estimates.minor_axis,
                 estimates.kappa,estimates.beta);
  opt.computeEstimates(sample_mean,S,estimates);
  return estimates;
}

/*!
 *  MML estimation
 */
struct Estimates Kent::computeMMLEstimates(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeMMLEstimates(sample_mean,S,data.size());
}

struct Estimates Kent::computeMMLEstimates(Vector &sample_mean, Matrix &S, long double N)
{
  string type = "MOMENT";
  struct Estimates estimates = computeMomentEstimates(sample_mean,S,N);
  print(type,estimates);
  long double msglen = computeMessageLength(estimates,sample_mean,S,N);
  cout << "msglen: " << msglen << endl;
  cout << "msglen (bpr): " << msglen/N << endl;

  type = "MAP";
  Optimize opt2(type);
  opt2.initialize(N,estimates.mean,estimates.major_axis,estimates.minor_axis,
                  estimates.kappa,estimates.beta);
  opt2.computeEstimates(sample_mean,S,estimates);
  print(type,estimates);
  msglen = computeMessageLength(estimates,sample_mean,S,N);
  cout << "msglen: " << msglen << endl;
  cout << "msglen (bpr): " << msglen/N << endl;

  /*if (N > 10) {
    //type = "MML_2";
    type = "MML_5";
    Optimize opt(type);
    opt.initialize(N,estimates.mean,estimates.major_axis,estimates.minor_axis,
                   estimates.kappa,estimates.beta);
    opt.computeEstimates(sample_mean,S,estimates);
    print(type,estimates);
    msglen = computeMessageLength(estimates,sample_mean,S,N);
    cout << "msglen: " << msglen << endl;
    cout << "msglen (bpr): " << msglen/N << endl;
  }*/
  return estimates;
}

void Kent::estimateParameters(std::vector<Vector> &data, Vector &weights)
{
  long double Neff;
  Vector sample_mean = computeVectorSum(data,weights,Neff);
  Matrix S = computeDispersionMatrix(data,weights);
  struct Estimates estimates = computeMMLEstimates(sample_mean,S,Neff);
  updateParameters(estimates);
}

void Kent::updateParameters(struct Estimates &estimates)
{
  mu = estimates.mean;
  major_axis = estimates.major_axis;
  minor_axis = estimates.minor_axis;
  kappa = estimates.kappa;
  beta = estimates.beta;
  psi = estimates.psi;
  alpha = estimates.alpha;
  eta = estimates.eta;
  assert(!boost::math::isnan(alpha));
  assert(!boost::math::isnan(psi));
  assert(!boost::math::isnan(eta));
  computeExpectation();
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

// 'other' is the approximate to the true distribution
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

long double Kent::computeKLDivergence(struct Estimates &estimates)
{
  Kent kent_est(estimates.mean,estimates.major_axis,estimates.minor_axis,
                estimates.kappa,estimates.beta);
  return computeKLDivergence(kent_est);
}

long double Kent::computeMessageLength(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeMessageLength(sample_mean,S,data.size());
}

long double Kent::computeMessageLength(Vector &sample_mean, Matrix &S, long double N)
{
  long double log_prior = computeLogPriorProbability();
  long double log_fisher = computeLogFisherInformation(N);
  long double part1 = -6.455 - log_prior + 0.5 * log_fisher;
  long double part2 = computeNegativeLogLikelihood(sample_mean,S,N) + 2.5
                 - 2 * N * log(AOM);
  long double msglen = part1 + part2;
  return msglen/log(2);
}

long double Kent::computeMessageLength(struct Estimates &estimates, 
                                       Vector &sample_mean, Matrix &S, long double N)
{
  Kent kent_est(estimates.mean,estimates.major_axis,estimates.minor_axis,
                estimates.kappa,estimates.beta);
  return kent_est.computeMessageLength(sample_mean,S,N);
}

void Kent::printParameters(ostream &os)
{
  os << "[mu]: "; print(os,mu,3);
  os << "\t[mj]: "; print(os,major_axis,3);
  os << "\t[mi]: "; print(os,minor_axis,3);
  os << "\t[kappa]:" << setw(10) << setprecision(3) << kappa;
  os << "\t[beta]:" << setw(10) << setprecision(3) << beta << endl;
  /*vector<long double> spherical(3,0);
  cartesian2spherical(estimates.mu,spherical);
  spherical[1] *= 180/PI; 
  spherical[2] *= 180/PI; 
  print(os,spherical);*/
}

