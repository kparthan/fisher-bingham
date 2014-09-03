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
  Vector spherical(3,0);
  cartesian2spherical(mu,spherical);
  alpha = spherical[1];
  eta = spherical[2];
  Matrix r = align_vector_with_zaxis(mu);
  Vector mj = prod(r,major_axis);
  cartesian2spherical(mj,spherical);
  psi = spherical[2];
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
  constants.kappa_E_x = Vector(3,0);
  for (int i=0; i<3; i++) {
    constants.E_x[i] = mu[i] * constants.ck_c;
    constants.kappa_E_x[i] = kappa * constants.E_x[i];
  }

  Matrix expectation_std = ZeroMatrix(3,3);
  expectation_std(0,0) = 0.5 * (1 - constants.ckk_c + constants.cb_c);  // lambda_1
  expectation_std(1,1) = 0.5 * (1 - constants.ckk_c - constants.cb_c);  // lambda_2
  expectation_std(2,2) = constants.ckk_c;  // lambda_3

  assert(expectation_std(0,0) > 0);
  assert(expectation_std(1,1) > 0);
  assert(expectation_std(2,2) > 0);
  //cout << "E_std: " << expectation_std << endl;

  Matrix tmp1 = prod(constants.R,expectation_std);
  constants.E_xx = prod(tmp1,constants.Rt);
  constants.beta_E_xx = ZeroMatrix(3,3);
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      constants.beta_E_xx(i,j) = constants.E_xx(i,j) * beta;
    }
  }

  computed = SET;
}

long double Kent::computeNegativeLogLikelihood(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeNegativeLogLikelihood(sample_mean,S,data.size());
}

long double Kent::computeNegativeLogLikelihood(Vector &sample_mean, Matrix &S, int N)
{
  long double c1 = computeDotProduct(sample_mean,mu);

  long double mj = prod_xMy(major_axis,S,major_axis);
  long double mi = prod_xMy(minor_axis,S,minor_axis);
  long double c2 = mj - mi;

  long double log_norm = computeLogNormalizationConstant();

  long double ans = N * log_norm - kappa * c1 - beta * c2;
  return ans;
}

long double Kent::computeLogPriorProbability()
{
  long double log_prior_axes = computeLogPriorAxes();
  long double log_prior_scale = computeLogPriorScale();
  long double log_joint_prior = log_prior_axes + log_prior_scale;
  return log_joint_prior;
}

long double Kent::computeLogPriorAxes()
{
  long double log_prior = 0;
  log_prior += log(4) - 4*log(PI);
  log_prior += log(sin(alpha));
  return log_prior;
}

long double Kent::computeLogPriorScale()
{
  long double log_prior = 0;
  log_prior += 2 * log(kappa);
  log_prior -= 2 * log(1+kappa*kappa);
  return log_prior;
}

long double Kent::computeLogFisherInformation()
{
  computeExpectation();
  long double log_det_axes = computeLogFisherAxes();
  long double log_det_kb = computeLogFisherScale();
  return log_det_axes + log_det_kb; 
}

long double Kent::computeLogFisherInformation(int N)
{
  long double log_fisher = computeLogFisherInformation(); 
  return log_fisher + 5 * log(N);
}

long double Kent::computeLogFisherAxes()
{
  computeFirstOrderDifferentials();
  computeSecondOrderDifferentials();
  computeFisherMatrixAxes();
  long double det = determinant(df.fisher_axes);
  //cout << "det: " << det << endl; //exit(1);
  return log(det);
}

void Kent::computeFirstOrderDifferentials()
{
  std::vector<Vector > first_order(3);
  Vector d1(3,0);

  // d1_mu[0]: dmu_da
  d1[0] = tc.cos_alpha * tc.cos_eta;  // d1 y11_da
  d1[1] = tc.cos_alpha * tc.sin_eta;  // d1 y12_da
  d1[2] = -tc.sin_alpha;              // d1 y13_da
  first_order[0] = d1;

  // d1_mu[1]: dmu_dn
  d1[0] = -tc.sin_alpha * tc.sin_eta;  // d1 y11_dn
  d1[1] = tc.sin_alpha * tc.cos_eta;   // d1 y12_dn
  d1[2] = 0;                           // d1 y13_dn
  first_order[1] = d1;

  // d1_mu[2]: dmu_ds = 0
  d1[0] = 0;
  d1[1] = 0;
  d1[2] = 0;
  first_order[2] = d1;

  df.d1_mu = first_order;

  computeDeltaDifferentials();
  // d1_mj[0]: dmj_da
  d1[0] = -tc.sin_psi * tc.sin_delta * df.ddel_da;
  d1[1] = tc.sin_psi * tc.cos_delta * df.ddel_da;
  d1[2] = 0; 
  first_order[0] = d1;

  // d1_mj[1]: dmj_dn
  d1[0] = -tc.sin_psi * tc.sin_delta;
  d1[1] = tc.sin_psi * tc.cos_delta;
  d1[2] = 0; 
  first_order[1] = d1;

  // d1_mj[2]: dmj_ds
  d1[0] = tc.cos_psi * tc.cos_delta - tc.sin_psi * tc.sin_delta * df.ddel_ds;
  d1[1] = tc.cos_psi * tc.sin_delta + tc.sin_psi * tc.cos_delta * df.ddel_ds;
  d1[2] = -tc.sin_psi;
  first_order[2] = d1;

  df.d1_mj = first_order;

  // d1_mi[0]: dmi_da
  // d1_mi[1]: dmi_dn
  // d1_mi[2]: dmi_ds
  for (int i=0; i<3; i++) {
    first_order[i] = computeFirstOrderDifferentialsMinorAxis(i);
  }
  df.d1_mi = first_order;
}

Vector Kent::computeFirstOrderDifferentialsMinorAxis(int i)
{
  Vector cross1 = crossProduct(df.d1_mu[i],major_axis);
  Vector cross2 = crossProduct(mu,df.d1_mj[i]);
  Vector ans(3,0);
  for (int i=0; i<3; i++) {
    ans[i] = cross1[i] + cross2[i];
  }
  return ans;
}

void Kent::computeDeltaDifferentials()
{
  long double tmp1,tmp2;

  // d(delta)_d(alpha)
  tmp1 = tc.sin_d_n * tc.tan_psi * tc.sin_alpha * tc.sin_alpha;
  df.ddel_da = -1/tmp1;

  // d(delta)_d(psi)
  tmp1 = tc.sin_d_n * tc.tan_alpha * tc.sin_psi * tc.sin_psi;
  df.ddel_ds = -1/tmp1;

  // d2(delta)_d(alpha)2
  tmp1 = tc.tan_psi * tc.sin_alpha * tc.sin_alpha * tc.tan_alpha;
  tmp2 = 2/tmp1;
  tmp1 = df.ddel_da * df.ddel_da * tc.cos_d_n;
  df.d2del_da2 = (tmp2 - tmp1) / tc.sin_d_n;

  // d2(delta)_d(psi)2
  tmp1 = tc.tan_alpha * tc.sin_psi * tc.sin_psi * tc.tan_psi;
  tmp2 = 2/tmp1;
  tmp1 = df.ddel_ds * df.ddel_ds * tc.cos_d_n;
  df.d2del_ds2 = (tmp2 - tmp1) / tc.sin_d_n;

  // d2(delta)_d(alpha)d(psi)
  tmp1 = tc.sin_alpha * tc.sin_alpha * tc.sin_psi * tc.sin_psi;
  tmp2 = 1/tmp1;
  tmp1 = df.ddel_da * df.ddel_ds * tc.cos_d_n;
  df.d2del_dads = (tmp2 - tmp1) / tc.sin_d_n;
}

void Kent::computeSecondOrderDifferentials()
{
  std::vector<Vector > second_order(6);
  Vector d2(3,0);

  // d2_mu[0]: d2mu_da2
  d2[0] = -tc.sin_alpha * tc.cos_eta;
  d2[1] = -tc.sin_alpha * tc.sin_eta;
  d2[2] = -tc.cos_alpha;
  second_order[0] = d2;

  // d2_mu[1]: d2mu_dn2
  d2[2] = 0;
  second_order[1] = d2;

  // d2_mu[2]: d2mu_ds2
  d2[0] = 0;
  d2[1] = 0;
  second_order[2] = d2;

  // d2_mu[3]: d2mu_dadn
  d2[0] = -tc.cos_alpha * tc.sin_eta;
  d2[1] = tc.cos_alpha * tc.cos_eta;
  second_order[3] = d2;

  // d2_mu[4]: d2mu_dads
  d2[0] = 0;
  d2[1] = 0;
  second_order[4] = d2;

  // d2_mu[5]: d2mu_dnds
  second_order[5] = d2;

  df.d2_mu = second_order;

  // d2_mj[0]: d2mj_da2
  d2[0] = -tc.sin_psi * (tc.sin_delta * df.d2del_da2 + tc.cos_delta * df.ddel_da * df.ddel_da);
  d2[1] = tc.sin_psi * (tc.cos_delta * df.d2del_da2 - tc.sin_delta * df.ddel_da * df.ddel_da);
  d2[2] = 0;
  second_order[0] = d2;

  // d2_mj[1]: d2mj_dn2
  d2[0] = -tc.sin_psi * tc.cos_delta;
  d2[1] = -tc.sin_psi * tc.sin_delta;
  d2[2] = 0; 
  second_order[1] = d2;

  // d2_mj[2]: d2mj_ds2
  d2[0] = -tc.cos_delta * tc.sin_psi - 2 * tc.cos_psi * tc.sin_delta * df.ddel_ds
          - tc.sin_psi * (tc.sin_delta * df.d2del_ds2 + tc.cos_delta * df.ddel_ds * df.ddel_ds);
  d2[1] = -tc.sin_delta * tc.sin_psi + 2 * tc.cos_delta * tc.cos_psi * df.ddel_ds
          + tc.sin_psi * (tc.cos_delta * df.d2del_ds2 - tc.sin_delta * df.ddel_ds * df.ddel_ds);
  d2[2] = -tc.cos_psi;
  second_order[2] = d2;

  // d2_mj[3]: d2mj_dadn
  d2[0] = -df.ddel_da * major_axis[0];
  d2[1] = -df.ddel_da * major_axis[1];
  d2[2] = 0;
  second_order[3] = d2;

  // d2_mj[4]: d2mj_dads
  d2[0] = -tc.sin_psi * (tc.sin_delta * df.d2del_dads + tc.cos_delta * df.ddel_da * df.ddel_ds)
          - tc.cos_psi * tc.sin_delta * df.ddel_da;
  d2[1] = tc.sin_psi * (tc.cos_delta * df.d2del_dads - tc.sin_delta * df.ddel_da * df.ddel_ds)
          + tc.cos_psi * tc.cos_delta * df.ddel_da;
  d2[2] = 0; 
  second_order[4] = d2;

  // d2_mj[5]: d2mj_dnds
  d2[0] = -tc.cos_psi * tc.sin_delta - tc.sin_psi * tc.cos_delta * df.ddel_ds;
  d2[1] = tc.cos_psi * tc.cos_delta - tc.sin_psi * tc.sin_delta * df.ddel_ds;
  d2[2] = 0; 
  second_order[5] = d2;

  df.d2_mj = second_order;

  // d2_mi[0]: d2mi_da2
  second_order[0] = computeSecondOrderDifferentialsMinorAxis(0,0,0,0,0,0);
  // d2_mi[1]: d2mi_dn2
  second_order[1] = computeSecondOrderDifferentialsMinorAxis(1,1,1,1,1,1);
  // d2_mi[2]: d2mi_ds2
  second_order[2] = computeSecondOrderDifferentialsMinorAxis(2,2,2,2,2,2);
  // d2_mi[3]: d2mi_dadn
  second_order[3] = computeSecondOrderDifferentialsMinorAxis(3,1,0,3,0,1);
  // d2_mi[4]: d2mi_dads
  second_order[4] = computeSecondOrderDifferentialsMinorAxis(4,2,0,4,0,2);
  // d2_mi[5]: d2mi_dnds
  second_order[5] = computeSecondOrderDifferentialsMinorAxis(5,2,1,5,1,2);

  df.d2_mi = second_order;
}

Vector Kent::computeSecondOrderDifferentialsMinorAxis(
  int i0, int i1, int i2, int i3, int i4, int i5
) {
  Vector c1 = crossProduct(df.d2_mu[i0],major_axis);
  Vector c2 = crossProduct(df.d1_mu[i1],df.d1_mj[i2]);
  Vector c3 = crossProduct(mu,df.d2_mj[i3]);
  Vector c4 = crossProduct(df.d1_mu[i4],df.d1_mj[i5]);

  Vector ans(3,0);
  for (int i=0; i<3; i++) {
    ans[i] = c1[i] + c2[i] + c3[i] + c4[i];
  }
  return ans;
}

void Kent::computeFisherMatrixAxes()
{
  df.fisher_axes = ZeroMatrix(3,3);
  // E[d^2 L / da^2]
  df.fisher_axes(0,0) = computeExpectationLikelihood(0,0,0);
  // E[d^2 L / dadn]
  df.fisher_axes(0,1) = computeExpectationLikelihood(3,0,1);
  // E[d^2 L / dads]
  df.fisher_axes(0,2) = computeExpectationLikelihood(4,0,2);
  // E[d^2 L / dn^2]
  df.fisher_axes(1,1) = computeExpectationLikelihood(1,1,1);
  // E[d^2 L / dnds]
  df.fisher_axes(1,2) = computeExpectationLikelihood(5,1,2);
  // E[d^2 L / ds^2]
  df.fisher_axes(2,2) = computeExpectationLikelihood(2,2,2);
  df.fisher_axes(1,0) = df.fisher_axes(0,1);
  df.fisher_axes(2,0) = df.fisher_axes(0,2);
  df.fisher_axes(2,1) = df.fisher_axes(1,2);
  //cout << "FisherAxes: " << df.fisher_axes << endl;
}

long double Kent::computeExpectationLikelihood(int i1, int i2, int i3)
{
  long double t1 = computeDotProduct(constants.kappa_E_x,df.d2_mu[i1]);

  long double t2 = prod_xMy(major_axis,constants.beta_E_xx,df.d2_mj[i1]);
  t2 += prod_xMy(df.d1_mj[i2],constants.beta_E_xx,df.d1_mj[i3]);
  t2 *= 2;

  long double t3 = prod_xMy(minor_axis,constants.beta_E_xx,df.d2_mi[i1]);
  t3 += prod_xMy(df.d1_mi[i2],constants.beta_E_xx,df.d1_mi[i3]);
  t3 *= 2;

  return -t1-t2+t3;
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

void Kent::computeAllEstimators(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  computeAllEstimators(sample_mean,S,data.size());
}

void Kent::computeAllEstimators(Vector &sample_mean, Matrix &S, int N)
{
  string type = "MOMENT";
  struct Estimates moment_est = computeMomentEstimates(sample_mean,S,N);
  print(type,moment_est);

  /*type = "MLE_UNCONSTRAINED";
  struct Estimates ml_est1 = moment_est;
  Optimize opt1(type);
  opt1.initialize(N,ml_est1.mean,ml_est1.major_axis,ml_est1.minor_axis,
                 ml_est1.kappa,ml_est1.beta);
  opt1.computeEstimates(sample_mean,S,ml_est1);
  print(type,ml_est1);*/

  type = "MML_SCALE";
  struct Estimates mml_est1 = moment_est;
  Optimize opt2(type);
  opt2.initialize(N,mml_est1.mean,mml_est1.major_axis,mml_est1.minor_axis,
                 mml_est1.kappa,mml_est1.beta);
  opt2.computeEstimates(sample_mean,S,mml_est1);
  print(type,mml_est1);

  /*type = "MML";
  struct Estimates mml_est2 = moment_est;
  Optimize opt3(type);
  opt3.initialize(N,mml_est2.mean,mml_est2.major_axis,mml_est2.minor_axis,
                 mml_est2.kappa,mml_est2.beta);
  opt3.computeEstimates(sample_mean,S,mml_est2);
  print(type,mml_est2);*/
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
struct Estimates Kent::computeMomentEstimates(Vector &sample_mean1, Matrix &S1, int N)
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

struct Estimates Kent::computeMLEstimates(Vector &sample_mean, Matrix &S, int N, string type)
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

struct Estimates Kent::computeMMLEstimates(Vector &sample_mean, Matrix &S, int N)
{
  struct Estimates estimates = computeMomentEstimates(sample_mean,S,N);
  Optimize opt("MML");
  opt.initialize(N,estimates.mean,estimates.major_axis,estimates.minor_axis,
                 estimates.kappa,estimates.beta);
  opt.computeEstimates(sample_mean,S,estimates);
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

