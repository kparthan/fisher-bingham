#include "Kent.h"
#include "FB6.h"
#include "Optimize.h"
#include "Support.h"
#include "Bingham.h"
#include "vMF.h"

extern Vector XAXIS,YAXIS,ZAXIS;
extern int VERBOSE,COMPUTE_KLDIV;
extern int ESTIMATION;
extern int PRIOR;

/*!
 *  Null constructor
 */
Kent::Kent() : kappa(1), beta(0)
{
  mu = XAXIS;
  major_axis = YAXIS;
  minor_axis = ZAXIS;
  computed = UNSET;
}

/*!
 *  Constructor
 */
Kent::Kent(double kappa, double beta) : kappa(kappa), beta(beta)
{
  mu = XAXIS;
  major_axis = YAXIS;
  minor_axis = ZAXIS;
  //assert(eccentricity() < 1);
  computed = UNSET;
}

/*!
 *  Constructor
 */
Kent::Kent(Vector &mu, Vector &major_axis, Vector &minor_axis,
          double kappa, double beta): 
          mu(mu), major_axis(major_axis), minor_axis(minor_axis), 
          kappa(kappa), beta(beta)
{
  //assert(eccentricity() < 1);
  Vector spherical(3,0);
  cartesian2spherical(mu,spherical);
  alpha = spherical[1];
  eta = spherical[2];
  Matrix r = align_vector_with_xaxis(alpha,eta);
  Vector mj = prod(r,major_axis);
  cartesian2spherical(mj,spherical);
  psi = spherical[2];
  computed = UNSET;
}

Kent::Kent(
  double psi, double alpha, double eta, double kappa, double beta
) : psi(psi), alpha(alpha), eta(eta), kappa(kappa), beta(beta)
{
  Matrix r = computeOrthogonalTransformation(psi,alpha,eta);
  mu = Vector(3,0); 
  major_axis = Vector(3,0); 
  minor_axis = Vector(3,0); 
  for (int i=0; i<3; i++) {
    mu[i] = r(i,0);
    major_axis[i] = r(i,1);
    minor_axis[i] = r(i,2);
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
  Matrix r1 = rotate_about_yaxis(PI/2);
  Matrix r2 = computeOrthogonalTransformation(mu,major_axis);
  Matrix transformation = prod(r2,r1);
  std::vector<Vector> sample = transform(canonical_sample,transformation);
  //std::vector<Vector> scaled_sample = scale_to_aom(sample);
  //writeToFile("random_sample_uniform.dat",sample,3);
  //writeToFile("random_sample_vmf.dat",sample,3);
  //writeToFile("random_sample_beta.dat",sample,3);
  //writeToFile("random_sample_new2.dat",sample,3);
  //writeToFile("random_sample.dat",sample,3);
  //computeLambertProjection(sample);
  return sample;
}

std::vector<Vector> Kent::generateCanonical(int N)
{
  //FB6 fb6(kappa,beta,0);
  //return fb6.generateCanonical(sample_size);

  // A = -beta * (mj mj' - mi mi')
  Matrix A = ZeroMatrix(3,3);
  A(0,0) = -beta;
  A(1,1) = beta;

  // B = 0.5 * kappa * (I - mu mu')
  Matrix B = ZeroMatrix(3,3);
  B(0,0) = kappa * 0.5;
  B(1,1) = kappa * 0.5;

  Matrix A1 = A + B;
  Bingham bingham(A1);

  std::vector<Vector> canonical_sample;
  Vector x;
  double u,check,exp1,exp2,tmp1,tmp2;
  int num_samples;
  std::vector<Vector> many_samples; 
  int accepted = 0;
  Vector kmu(3,0); kmu[2] = kappa;
  while (accepted != N) {
    num_samples = 2 * (N-accepted);
    many_samples = bingham.generate(num_samples);
    for (int i=0; i<many_samples.size(); i++) {
      u = uniform_random();
      x = many_samples[i];

      tmp1 = computeDotProduct(kmu,x);
      tmp2 = prod_vMv(x,A);
      exp1 = tmp1 - tmp2;
      tmp1 = prod_vMv(x,A1);
      exp2 = kappa - tmp1; 

      check =  exp1 - exp2;
      if (log(u) < check) {
        canonical_sample.push_back(x);
        if (++accepted == N) {
          goto finish;
        }
      }
    } // for() ends ...
  } // while() ends ..
  finish:
  return canonical_sample;
}

/*!
 *  e = 2 beta / kappa
 */
double Kent::eccentricity()
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
double Kent::computeLogNormalizationConstant(void)
{
  return computeSeriesSum(kappa,beta,0.5);
}

double Kent::log_dc_dk(void)
{
  return computeSeriesSum(kappa,beta,1.5);
}

double Kent::log_d2c_dk2(void)
{
  double x = computeSeriesSum(kappa,beta,2.5); // log a
  double y = constants.log_ck; // log b
  double tmp1 = x - y;
  double tmp2 = exp(tmp1) + (1/kappa);
  double ans = y + log(tmp2);
  return ans;
}

double Kent::computeSeriesSum(double k, double b, double d)
{
  double ex = 2 * b/k; 
  double log_ex = 2 * log(ex);
  //double log_bessel = logModifiedBesselFirstKind(0,k);
  double log_bessel_prev,log_bessel_current;
  double log_f0,log_fj_current,log_fj_prev;
  double current,log_diff,series_sum=1;
  int j = 0;
  double gj;

  //log_bessel_prev = log(cyl_bessel_i(d,k));
  log_bessel_prev = computeLogModifiedBesselFirstKind(d,k);
  log_f0 = lgamma<double>(0.5) + log_bessel_prev;
  log_fj_prev = log_f0;
  while (1) {
    d += 2;
    //log_bessel_current = log(cyl_bessel_i(d,k));
    log_bessel_current = computeLogModifiedBesselFirstKind(d,k);
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
  double ans = log(2*PI) + 0.5 * log(2.0/k);
  ans += log_f0;
  ans += log(series_sum);
  //cout << "j: " << j << endl;
  //cout << "ans: " << ans << endl;
  return ans;
}

double Kent::log_dc_db()
{
  return computeSeriesSum2(kappa,beta,2.5);
}

double Kent::log_d2c_dkdb()
{
  return computeSeriesSum2(kappa,beta,3.5);
}

double Kent::computeSeriesSum2(double k, double b, double d)
{
  double ex = 2 * b/k; 
  double log_ex = 2 * log(ex);
  double log_bessel_prev,log_bessel_current;
  double log_f1,log_fj_current,log_fj_prev;
  double current,log_diff,series_sum=1;
  int j = 1;
  double gj;

  //log_bessel_prev = log(cyl_bessel_i(d,k));
  log_bessel_prev = computeLogModifiedBesselFirstKind(d,k);
  log_f1 = lgamma<double>(1.5) + log_ex + log_bessel_prev;
  log_fj_prev = log_f1;
  while (1) {
    d += 2;
    //log_bessel_current = log(cyl_bessel_i(d,k));
    log_bessel_current = computeLogModifiedBesselFirstKind(d,k);
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
  double ans = log(2*PI) + 0.5 * log(2.0/k) + log(2/b);
  ans += log_f1;
  ans += log(series_sum);
  //cout << "j: " << j << endl;
  return ans;
}

double Kent::log_d2c_db2()
{
  double b = beta;
  double k = kappa;
  double ex = 2 * b/k; 
  double log_ex = 2 * log(ex);
  double log_bessel_prev,log_bessel_current;
  double log_f1,log_fj_current,log_fj_prev;
  double current,log_diff,series_sum=1;
  double d = 2.5;
  int j = 1;
  double gj;

  //log_bessel_prev = log(cyl_bessel_i(d,k));
  log_bessel_prev = computeLogModifiedBesselFirstKind(d,k);
  log_f1 = lgamma<double>(1.5) + log_ex + log_bessel_prev;
  log_fj_prev = log_f1;
  while (1) {
    d += 2;
    //log_bessel_current = log(cyl_bessel_i(d,k));
    log_bessel_current = computeLogModifiedBesselFirstKind(d,k);
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
  double ans = log(2*PI) + 0.5 * log(2.0/k) + log(2) - 2*log(b);
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

  double tmp;
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
  expectation_std(0,0) = constants.ckk_c;  // lambda_1
  constants.lambda1 = expectation_std(0,0);
  expectation_std(1,1) = 0.5 * (1 - constants.ckk_c + constants.cb_c);  // lambda_2
  constants.lambda2 = expectation_std(1,1);
  expectation_std(2,2) = 0.5 * (1 - constants.ckk_c - constants.cb_c);  // lambda_3
  constants.lambda3 = expectation_std(2,2);

  /*assert(expectation_std(0,0) > 0);
  assert(expectation_std(1,1) > 0);
  assert(expectation_std(2,2) > 0);*/
  //cout << "E_std: " << expectation_std << endl;

  Matrix tmp1 = prod(constants.R,expectation_std);
  constants.E_xx = prod(tmp1,constants.Rt);

  computed = SET;
}

double Kent::log_density(Vector &x)
{
  double c1 = computeDotProduct(mu,x);

  Matrix xx = outer_prod(x,x);
  double mj = prod_xMy(major_axis,xx,major_axis);
  double mi = prod_xMy(minor_axis,xx,minor_axis);
  double c2 = mj - mi;

  double log_norm; 
  if (computed == UNSET) {
    log_norm = computeLogNormalizationConstant();
  } else {
    log_norm = constants.log_c;
  }

  double ans = -log_norm + kappa * c1 + beta * c2;
  return ans;
}

double Kent::computeNegativeLogLikelihood(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeNegativeLogLikelihood(sample_mean,S,data.size());
}

double Kent::computeNegativeLogLikelihood(Vector &sample_mean, Matrix &S, double N)
{
  double c1 = computeDotProduct(sample_mean,mu);

  double mj = prod_xMy(major_axis,S,major_axis);
  double mi = prod_xMy(minor_axis,S,minor_axis);
  double c2 = mj - mi;

  double log_norm;
  if (computed == UNSET) {
    log_norm = computeLogNormalizationConstant();
  } else {
    log_norm = constants.log_c;
  }

  double ans = N * log_norm - kappa * c1 - beta * c2;
  return ans;
}

double Kent::computeNegativeLogLikelihood(
  struct Estimates &estimates, 
  Vector &sample_mean, 
  Matrix &S, 
  double N
) {
  Kent kent_est(estimates.mean,estimates.major_axis,estimates.minor_axis,
                estimates.kappa,estimates.beta);
  return kent_est.computeNegativeLogLikelihood(sample_mean,S,N);
}

double Kent::computeLogParametersProbability(double Neff)
{
  double log_prior_density = computeLogPriorProbability();
  double log_expected_fisher = computeLogFisherInformation(Neff);
  double logp = -log_prior_density + 0.5 * log_expected_fisher;
  return logp;
}

double Kent::computeLogPriorProbability()
{
  double log_prior_axes = computeLogPriorAxes();
  double log_prior_scale = computeLogPriorScale();
  double log_joint_prior = log_prior_axes + log_prior_scale;
  assert(!boost::math::isnan(log_joint_prior));
  return log_joint_prior;
}

double Kent::computeLogPriorAxes()
{
  double angle = alpha;
  if (angle < TOLERANCE) angle = TOLERANCE;
  if (fabs(angle-PI) < TOLERANCE) angle = PI-TOLERANCE;
  
  double log_prior = 0;
  //if (ESTIMATION == MML) {
    log_prior += (log(sin(angle)));
  //}
  log_prior -= (log(4) + 2*log(PI));
  return log_prior;
}

double Kent::computeLogPriorScale()
{
  double log_prior = 0;

  if (PRIOR == 3) { // ... using vMF (3D) kappa prior
    log_prior += log(8/PI);
    log_prior += log(kappa);
    log_prior -= (2 * log(1+kappa*kappa));
  } else if (PRIOR == 2) {  // ... using vMF (2D) kappa prior
    log_prior += log(2);
    log_prior -= (1.5 * log(1+kappa*kappa));
  }

  // vmf beta prior
  /*log_prior += (2 * log(4/PI));
  log_prior += (2 * log(kappa));
  log_prior -= (2 * log(1+kappa*kappa));

  log_prior += (2 * log(beta));
  log_prior -= (2 * log(1+beta*beta));*/
  /*double tmp = atan(0.5*kappa);
  tmp -= (2*kappa/(kappa*kappa+4));
  log_prior -= log(tmp);*/
  //log_prior -= log(0.39629);

  // beta distribution priors
  /*double ex = 2 * beta/kappa;
  if (ex >= 1) ex = 1 - TOLERANCE;
  beta_distribution<> beta_dist(2,2);
  //beta_distribution<> beta_dist(2,10);
  //beta_distribution<> beta_dist(10,2);
  double f = pdf(beta_dist,ex);
  log_prior += log(f);*/

  // gamma distribution prior
  //boost::gamma_distribution<> gamma_dist(2);
  //double k = 3, t = 2;  // k = shape; t = scale
  //double f = boost::math::gamma_p_derivative(k,beta/t)/t;
  //log_prior += log(f);
  //double log_f = 2 * log(beta) - log(16) - (0.5 * beta);
  //log_prior += log_f;

  // uniform priors
  //double K1 = 1,K2=100;
  //log_prior = log(4) - log(K2+K1) - log(K2-K1);
  //log_prior -= log(log(K2/K1));
  //log_prior += (log(2) - 2*log(kappa));

  return log_prior;
}

double Kent::computeLogFisherInformation_Single(double N)
{
  if (computed == UNSET) {
    computeExpectation();
  }
  double log_det_axes = computeLogFisherAxes(N);
  double log_det_kb = computeLogFisherScale();
  return log_det_axes + log_det_kb; 
}

double Kent::computeLogFisherInformation(double N)
{
  double log_fisher = computeLogFisherInformation_Single(N); 
  return log_fisher + 5 * log(N);
}

double Kent::computeLogFisherAxes(double N)
{
  double sin_psi = sin(psi);
  double sinsq_psi = sin_psi * sin_psi;
  double cos_psi = cos(psi);
  double cossq_psi = cos_psi * cos_psi; 
  double sin_alpha = sin(alpha);
  double sinsq_alpha = sin_alpha * sin_alpha;
  double cos_alpha = cos(alpha);
  double cossq_alpha = cos_alpha * cos_alpha; 

  constants.fisher_axes = ZeroMatrix(3,3);
  
  //double N = 10;
  double ans,tmp1,tmp2,tmp3;
  // E [d^2 L / d a^2]
  ans = kappa * constants.ck_c;
  tmp1 = (constants.lambda1 - constants.lambda3) * sinsq_psi;
  tmp2 = (constants.lambda1 - constants.lambda2) * cossq_psi;
  ans += (2 * beta * (tmp1 - tmp2));
  constants.fisher_axes(0,0) = ans + (3/(PI*PI*N));

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
  ans += (2 * beta * tmp3);
  constants.fisher_axes(1,1) = ans+ (3/(PI*PI*N));
  
  // E [d^2 L / d s^2]
  ans = 4 * beta * constants.cb_c;
  constants.fisher_axes(2,2) = ans+ (3/(PI*PI*N));

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

  double det = determinant(constants.fisher_axes);
  //cout << "det: " << det << endl; //exit(1);
  return log(det);
}

double Kent::computeLogFisherScale()
{
  // E [d^2 L / d k^2]
  double t1 = constants.ckk_c - (constants.ck_c * constants.ck_c);

  // E [d^2 L / d b^2]
  double t2 = constants.cbb_c - (constants.cb_c * constants.cb_c);

  // E [d^2 L / dk db]
  double t3 = constants.ckb_c - (constants.ck_c * constants.cb_c);

  double det = t1 * t2 - t3 * t3;
  //cout << "det: " << det << endl;
  return log(det);
}

void Kent::computeAllEstimators(
  std::vector<Vector> &data, 
  std::vector<struct Estimates> &all_estimates,
  int verbose,
  int compute_kldiv
) {
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);

  computeAllEstimators(
    sample_mean,S,data.size(),all_estimates,verbose,compute_kldiv
  );
}

void Kent::computeAllEstimators(
  Vector &sample_mean, 
  Matrix &S, 
  double N,
  std::vector<struct Estimates> &all_estimates,
  int verbose,
  int compute_kldiv
) {
  double msglen,negloglike,kldiv,min_msg;
  int min_index;

  all_estimates.clear();

  struct Estimates asymptotic_est = computeAsymptoticMomentEstimates(sample_mean,S,N);

  string type = "MOMENT";
  struct Estimates moment_est = asymptotic_est;
  Optimize opt(type);
  opt.initialize(
    N,moment_est.mean,moment_est.major_axis,moment_est.minor_axis,
    moment_est.kappa,moment_est.beta
  );
  opt.computeEstimates(sample_mean,S,moment_est);
  moment_est.msglen = computeMessageLength(moment_est,sample_mean,S,N);
  moment_est.negloglike = computeNegativeLogLikelihood(moment_est,sample_mean,S,N);
  if (compute_kldiv) {
    moment_est.kldiv = computeKLDivergence(moment_est);
  }
  if (verbose) {
    print(type,moment_est);
    cout << fixed << "msglen: " << moment_est.msglen << endl;
    cout << "negloglike: " << moment_est.negloglike << endl;
    cout << "KL-divergence: " << moment_est.kldiv << endl << endl;
  }
  all_estimates.push_back(moment_est);
  min_msg = moment_est.msglen;
  min_index = MOMENT;

  type = "MLE";
  struct Estimates ml_est = moment_est;
  //struct Estimates ml_est = asymptotic_est;
  Optimize opt1(type);
  opt1.initialize(
    N,ml_est.mean,ml_est.major_axis,ml_est.minor_axis,
    ml_est.kappa,ml_est.beta
  );
  opt1.computeEstimates(sample_mean,S,ml_est);
  ml_est.msglen = computeMessageLength(ml_est,sample_mean,S,N);
  ml_est.negloglike = computeNegativeLogLikelihood(ml_est,sample_mean,S,N);
  if (compute_kldiv) {
    ml_est.kldiv = computeKLDivergence(ml_est);
  }
  if (verbose) {
    print(type,ml_est);
    cout << fixed << "msglen: " << ml_est.msglen << endl;
    cout << "negloglike: " << ml_est.negloglike << endl;
    cout << "KL-divergence: " << ml_est.kldiv << endl << endl;
  }
  all_estimates.push_back(ml_est);
  if (ml_est.msglen < min_msg) {
    min_index = MLE;
    min_msg = ml_est.msglen;
  }

  type = "MAP";
  struct Estimates map_est = moment_est;
  Optimize opt2(type);
  opt2.initialize(
    N,map_est.mean,map_est.major_axis,map_est.minor_axis,
    map_est.kappa,map_est.beta
  );
  opt2.computeEstimates(sample_mean,S,map_est);
  map_est.msglen = computeMessageLength(map_est,sample_mean,S,N);
  map_est.negloglike = computeNegativeLogLikelihood(map_est,sample_mean,S,N);
  if (compute_kldiv) {
    map_est.kldiv = computeKLDivergence(map_est);
  }
  if (verbose) {
    print(type,map_est);
    cout << fixed << "msglen: " << map_est.msglen << endl;
    cout << "negloglike: " << map_est.negloglike << endl;
    cout << "KL-divergence: " << map_est.kldiv << endl << endl;
  }
  all_estimates.push_back(map_est);
  if (map_est.msglen < min_msg) {
    min_index = MAP;
    min_msg = map_est.msglen;
  }

  type = "MML";
  //struct Estimates mml_est = asymptotic_est;
  struct Estimates mml_est = all_estimates[min_index];
  Optimize opt3(type);
  opt3.initialize(
    N,mml_est.mean,mml_est.major_axis,mml_est.minor_axis,
    mml_est.kappa,mml_est.beta
  );
  opt3.computeEstimates(sample_mean,S,mml_est);
  mml_est.msglen = computeMessageLength(mml_est,sample_mean,S,N);
  mml_est.negloglike = computeNegativeLogLikelihood(mml_est,sample_mean,S,N);
  if (compute_kldiv) {
    mml_est.kldiv = computeKLDivergence(mml_est);
  }
  if (verbose) {
    cout << endl;
    print(type,mml_est);
    cout << fixed << "msglen: " << mml_est.msglen << endl;
    cout << "negloglike: " << mml_est.negloglike << endl;
    cout << "KL-divergence: " << mml_est.kldiv << endl << endl;
  }
  all_estimates.push_back(mml_est);

  /***************************** MAP Variants *********************************/

  type = "MAP_ECCENTRICITY_TRANSFORM";
  struct Estimates map2_est = moment_est;
  Optimize optmap2(type);
  optmap2.initialize(
    N,map2_est.mean,map2_est.major_axis,map2_est.minor_axis,
    map2_est.kappa,map2_est.beta
  );
  optmap2.computeEstimates(sample_mean,S,map2_est);
  map2_est.msglen = computeMessageLength(map2_est,sample_mean,S,N);
  map2_est.negloglike = computeNegativeLogLikelihood(map2_est,sample_mean,S,N);
  if (compute_kldiv) {
    map2_est.kldiv = computeKLDivergence(map2_est);
  }
  if (verbose) {
    cout << endl;
    print(type,map2_est);
    cout << fixed << "msglen: " << map2_est.msglen << endl;
    cout << "negloglike: " << map2_est.negloglike << endl;
    cout << "KL-divergence: " << map2_est.kldiv << endl << endl;
  }
  all_estimates.push_back(map2_est);

  if (PRIOR == 2) {
    type = "MAP_UNIFORM_TRANSFORM";
    struct Estimates map3_est = moment_est;
    Optimize optmap3(type);
    optmap3.initialize(
      N,map3_est.mean,map3_est.major_axis,map3_est.minor_axis,
      map3_est.kappa,map3_est.beta
    );
    optmap3.computeEstimates(sample_mean,S,map3_est);
    map3_est.msglen = computeMessageLength(map3_est,sample_mean,S,N);
    map3_est.negloglike = computeNegativeLogLikelihood(map3_est,sample_mean,S,N);
    if (compute_kldiv) {
      map3_est.kldiv = computeKLDivergence(map3_est);
    }
    if (verbose) {
      cout << endl;
      print(type,map3_est);
      cout << fixed << "msglen: " << map3_est.msglen << endl;
      cout << "negloglike: " << map3_est.negloglike << endl;
      cout << "KL-divergence: " << map3_est.kldiv << endl << endl;
    }
    all_estimates.push_back(map3_est);
  }
}

/*!
 *  Moment estimation (with +X as the north pole)
 */
struct Estimates Kent::computeAsymptoticMomentEstimates(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeAsymptoticMomentEstimates(sample_mean,S,data.size());
}

/*!
 *  Moment estimation (with +X as the north pole) 
 *  (similar to used in Kent (1982) paper)
 *  (sample_mean1 and S1 are \sum_x and \sum_ x x')
 */
struct Estimates Kent::computeAsymptoticMomentEstimates(
  Vector &sample_mean1, Matrix &S1, double N
) {
  Vector sample_mean = sample_mean1;
  Matrix S = S1; 
  for (int i=0; i<3; i++) {
    sample_mean[i] /= N;
    for (int j=0; j<3; j++) {
      S(i,j) /= N;
    }
  }
  // compute r1:
  double r1 = norm(sample_mean);

  struct Estimates estimates;
  Vector spherical(3,0);
  cartesian2spherical(sample_mean,spherical);
  double theta = spherical[1];
  double phi = spherical[2];

  // rotation matrix to align north pole
  Matrix Ht = align_vector_with_xaxis(theta,phi);
  Matrix H = trans(Ht);
  Matrix tmp = prod(Ht,S);
  Matrix B = prod(tmp,H);
  double ratio = 2 * B(1,2) / (B(1,1) - B(2,2));
  double psi = 0.5 * atan(ratio);
  Matrix K = IdentityMatrix(3,3);
  K(1,1) = cos(psi); K(1,2) = -sin(psi);
  K(2,1) = -K(1,2); K(2,2) = K(1,1);
  Matrix HK = prod(H,K);

  /*Matrix BL(2,2);
  BL(0,0) = B(1,1); BL(0,1) = B(1,2);
  BL(1,0) = B(2,1); BL(1,1) = B(2,2);
  Vector eigen_values(2,0);
  Matrix eigen_vectors = IdentityMatrix(2,2);
  eigenDecomposition(BL,eigen_values,eigen_vectors);
  cout << "eigen values: "; print(cout,eigen_values,3); cout << endl;*/

  // compute r2:
  //cout << "B: " << B << endl;
  double t1 = (B(1,1)-B(2,2)) * (B(1,1)-B(2,2));
  double t2 = 4 * B(1,2) * B(1,2);
  double r2 = sqrt(t1+t2);
  //cout << "r2: " << r2 << endl;

  estimates.eig_max = 0.5 * (B(1,1) + B(2,2) + r2);
  //cout << "eig_max: " << eig_max << endl;
  //cout << "angle: " << acos(sqrt(1-eig_max)) * 180/PI << endl;

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

  double s,a,e;
  computeOrthogonalTransformation(estimates.mean,estimates.major_axis,s,a,e);
  /*cout << "before: \n";
  cout << "mean_est: "; print(cout,estimates.mean,3); cout << endl;
  cout << "major_est: "; print(cout,estimates.major_axis,3); cout << endl;
  cout << "minor_est: "; print(cout,estimates.minor_axis,3); cout << endl;
  cout << "psi: " << s*180/PI << "; alpha: " << a*180/PI << "; eta: " << e*180/PI << endl;*/
  if (s >= PI) s -= PI;
  assert(s >= 0 && s < PI);
  Matrix r = computeOrthogonalTransformation(s,a,e);
  estimates.psi = s; estimates.alpha = a; estimates.eta = e;
  estimates.mean = prod(r,XAXIS);
  estimates.major_axis = prod(r,YAXIS);
  estimates.minor_axis = prod(r,ZAXIS);
  /*cout << "after: \n";
  cout << "mean_est: "; print(cout,estimates.mean,3); cout << endl;
  cout << "major_est: "; print(cout,estimates.major_axis,3); cout << endl;
  cout << "minor_est: "; print(cout,estimates.minor_axis,3); cout << endl;
  cout << "psi: " << s*180/PI << "; alpha: " << a*180/PI << "; eta: " << e*180/PI << endl;*/

  // estimate kappa, beta
  double f1 = 1/(2 - 2*r1 - r2);
  double f2 = 1/(2 - 2*r1 + r2);
  estimates.kappa = f1 + f2;
  estimates.beta = 0.5 * (f1-f2);

  //cout << "(asymptotic) kappa: " << estimates.kappa << endl;
  //cout << "(asymptotic) beta: " << estimates.beta << endl;

  /*Optimize opt("MOMENT");
  opt.initialize(N,estimates.mean,estimates.major_axis,estimates.minor_axis,
                 estimates.kappa,estimates.beta);
  opt.computeEstimates(sample_mean1,S1,estimates);*/
  return estimates;
}

struct Estimates Kent::computeMomentEstimates(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeMomentEstimates(sample_mean,S,data.size());
}

struct Estimates Kent::computeMomentEstimates(Vector &sample_mean, Matrix &S, double N)
{
  struct Estimates estimates = computeAsymptoticMomentEstimates(sample_mean,S,N);
  string type = "MOMENT";
  Optimize opt(type);
  opt.initialize(N,estimates.mean,estimates.major_axis,estimates.minor_axis,
                 estimates.kappa,estimates.beta);
  opt.computeEstimates(sample_mean,S,estimates);
  return estimates;
}

/*!
 *  Max LH estimation
 */
struct Estimates Kent::computeMLEstimates(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeMLEstimates(sample_mean,S,data.size());
}

struct Estimates Kent::computeMLEstimates(Vector &sample_mean, Matrix &S, double N)
{
  struct Estimates estimates = computeAsymptoticMomentEstimates(sample_mean,S,N);
  string type = "MLE";
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

struct Estimates Kent::computeMMLEstimates(Vector &sample_mean, Matrix &S, double N)
{
  struct Estimates asymptotic_est = computeAsymptoticMomentEstimates(sample_mean,S,N);

  string type = "MOMENT";
  struct Estimates moment_est = asymptotic_est;
  Optimize opt1(type);
  opt1.initialize(N,moment_est.mean,moment_est.major_axis,moment_est.minor_axis,
                 moment_est.kappa,moment_est.beta);
  opt1.computeEstimates(sample_mean,S,moment_est);
  print(type,moment_est);
  double msglen = computeMessageLength(moment_est,sample_mean,S,N);
  cout << "msglen (bpr): " << msglen/N << endl;

  /*type = "MAP";
  struct Estimates map_est = moment_est;
  Optimize opt2(type);
  opt2.initialize(N,map_est.mean,map_est.major_axis,map_est.minor_axis,
                  map_est.kappa,map_est.beta);
  opt2.computeEstimates(sample_mean,S,map_est);
  print(type,map_est);
  msglen = computeMessageLength(map_est,sample_mean,S,N);
  cout << "msglen (bpr): " << msglen/N << endl;*/

  //struct Estimates mml_est = map_est;
  struct Estimates mml_est = moment_est;
  //if (N > 10) {
    type = "MML";
    Optimize opt3(type);
    opt3.initialize(N,mml_est.mean,mml_est.major_axis,mml_est.minor_axis,
                   mml_est.kappa,mml_est.beta);
    opt3.computeEstimates(sample_mean,S,mml_est);
    print(type,mml_est);
    msglen = computeMessageLength(mml_est,sample_mean,S,N);
    cout << "msglen (bpr): " << msglen/N << endl;
  //}
  return mml_est;
}

void Kent::estimateParameters(std::vector<Vector> &data, Vector &weights)
{
  double Neff;
  Vector sample_mean = computeVectorSum(data,weights,Neff);
  Matrix S = computeDispersionMatrix(data,weights);

  struct Estimates estimates;

  struct Estimates asymptotic_est = computeAsymptoticMomentEstimates(sample_mean,S,Neff);

  string type = "MOMENT";
  struct Estimates moment_est = asymptotic_est;
  Optimize opt(type);
  opt.initialize(Neff,moment_est.mean,moment_est.major_axis,moment_est.minor_axis,
                 moment_est.kappa,moment_est.beta);
  opt.computeEstimates(sample_mean,S,moment_est);
  moment_est.msglen = computeMessageLength(moment_est,sample_mean,S,Neff);
  //print(type,moment_est);

  switch(ESTIMATION) {
    case MOMENT:
    {
      estimates = moment_est;
      break;
    }

    case MLE:
    {
      type = "MLE";
      struct Estimates ml_est = moment_est;
      Optimize opt1(type);
      opt1.initialize(Neff,ml_est.mean,ml_est.major_axis,ml_est.minor_axis,
                      ml_est.kappa,ml_est.beta);
      opt1.computeEstimates(sample_mean,S,ml_est);
      estimates = ml_est;
      break;
    }

    case MAP:
    {
      type = "MAP";
      struct Estimates map_est = moment_est;
      Optimize opt2(type);
      opt2.initialize(Neff,map_est.mean,map_est.major_axis,map_est.minor_axis,
                      map_est.kappa,map_est.beta);
      opt2.computeEstimates(sample_mean,S,map_est);
      estimates = map_est;
      break;
    }

    case MML:
    {
      type = "MAP";
      struct Estimates map_est = moment_est;
      Optimize opt2(type);
      opt2.initialize(Neff,map_est.mean,map_est.major_axis,map_est.minor_axis,
                      map_est.kappa,map_est.beta);
      opt2.computeEstimates(sample_mean,S,map_est);
      map_est.msglen = computeMessageLength(map_est,sample_mean,S,Neff);

      type = "MML";
      struct Estimates mml_est; 
      if (moment_est.msglen < map_est.msglen) mml_est = moment_est;
      else mml_est = map_est;
      Optimize opt3(type);
      opt3.initialize(Neff,mml_est.mean,mml_est.major_axis,mml_est.minor_axis,
                      mml_est.kappa,mml_est.beta);
      opt3.computeEstimates(sample_mean,S,mml_est);
      estimates = mml_est;
      //print(type,mml_est);
      break;
    }
  } // switch()
  
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

double Kent::Kappa()
{
  return kappa;
}

double Kent::Beta()
{
  return beta;
}

// 'other' is the approximate to the true distribution
double Kent::computeKLDivergence(Kent &other)
{
  struct Constants constants1 = getConstants();
  struct Constants constants2 = other.getConstants();

  double ans = constants2.log_c - constants1.log_c;
  
  double kappa2 = other.Kappa();
  Vector mu2 = other.Mean();
  Vector kmu(3,0);
  for (int i=0; i<3; i++) {
    kmu[i] = kappa * mu[i] - kappa2 * mu2[i];
  }
  ans += computeDotProduct(kmu,constants1.E_x);

  double tmp1,tmp2;
  tmp1 = prod_vMv(major_axis,constants1.E_xx);
  tmp2 = prod_vMv(minor_axis,constants1.E_xx);
  ans += beta * (tmp1 - tmp2);

  double beta2 = other.Beta();
  Vector mj2 = other.MajorAxis();
  Vector mi2 = other.MinorAxis();
  tmp1 = prod_vMv(mj2,constants1.E_xx);
  tmp2 = prod_vMv(mi2,constants1.E_xx);
  ans -= beta2 * (tmp1 - tmp2);

  return ans/log(2);
}

double Kent::computeKLDivergence(struct Estimates &estimates)
{
  Kent kent_est(estimates.mean,estimates.major_axis,estimates.minor_axis,
                estimates.kappa,estimates.beta);
  return computeKLDivergence(kent_est);
}

double Kent::computeMessageLength(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeMessageLength(sample_mean,S,data.size());
}

double Kent::computeMessageLength(Vector &sample_mean, Matrix &S, double N)
{
  double log_prior = computeLogPriorProbability();
  double log_fisher = computeLogFisherInformation(N);
  double part1 = -6.455 - log_prior + 0.5 * log_fisher;
  double part2 = computeNegativeLogLikelihood(sample_mean,S,N) + 2.5
                 - 2 * N * log(AOM);
  double msglen = part1 + part2;
  return msglen/log(2);
}

double Kent::computeMessageLength(struct Estimates &estimates, 
                                  Vector &sample_mean, Matrix &S, double N)
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
  os << "\t[kappa]: " << fixed << setprecision(3) << kappa << "\t\t";
  os << "\t[beta]: " << fixed << setprecision(3) << beta << endl;
  /*vector<double> spherical(3,0);
  cartesian2spherical(estimates.mu,spherical);
  spherical[1] *= 180/PI; 
  spherical[2] *= 180/PI; 
  print(os,spherical);*/
}

double Kent::computeTestStatistic_vMF(std::vector<Vector> &x)
{
  int N = x.size();

  Vector sample_mean = computeNormalizedVectorSum(x);
  double rbar = norm(sample_mean);
  //cout << "sample mean: "; print(cout,sample_mean,3); cout << endl;

  Matrix S = computeNormalizedDispersionMatrix(x);
  
  Matrix Ht = align_vector_with_xaxis(sample_mean);
  Matrix H = trans(Ht);
  //cout << "H: " << H << endl;
  Matrix tmp = prod(Ht,S);
  Matrix B = prod(tmp,H);

  double t1 = (B(1,1)-B(2,2)) * (B(1,1)-B(2,2));
  double t2 = 4 * B(1,2) * B(1,2);
  double w = t1+t2;  // r2^2

  // Assumption: x comes from vMF
  // Estimate vMF kappa
  vMF vmf;
  std::vector<struct Estimates_vMF> all_estimates;
  vmf.computeAllEstimators(x,all_estimates);
  double k,t;
  chi_squared chisq(2);
  double alpha = 0.05,pvalue;
  for (int i=0; i<all_estimates.size(); i++) {
    k = all_estimates[i].kappa;
    t = computeTestStatistic(k,w,rbar,N);
    pvalue = compute_pvalue(t,chisq);
    cout << scientific << "t: " << t << "; pvalue: " << pvalue << endl;
  }
  return t;
}

double Kent::computeConfidenceRegion(std::vector<Vector> &x)
{
  std::vector<struct Estimates> all_estimates;
  computeAllEstimators(x,all_estimates,1,1);

  double m,l2,l3,s2,s3,area;
  for (int i=0; i<all_estimates.size(); i++) {
    Kent kent(all_estimates[i].mean,all_estimates[i].major_axis,all_estimates[i].minor_axis,
              all_estimates[i].kappa,all_estimates[i].beta);
    struct Constants consts = kent.getConstants();
    m = consts.ck_c;
    l2 = consts.lambda2;
    l3 = consts.lambda3;
    s2 = sqrt(l2);
    s3 = sqrt(l3);
    area = s2 * s3 / (m * m);
    area *= (PI / x.size());
    cout << "area: " << scientific << area << endl;
  }
}

