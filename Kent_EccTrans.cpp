#include "Kent_EccTrans.h"
#include "Support.h"

extern Vector XAXIS,YAXIS,ZAXIS;
extern int VERBOSE,COMPUTE_KLDIV;
extern int ESTIMATION;
extern int PRIOR;

/*!
 *  Null constructor
 */
Kent_EccTrans::Kent_EccTrans() : kappa(1), ecc(0)
{
  mu = XAXIS;
  major_axis = YAXIS;
  minor_axis = ZAXIS;
  computed = UNSET;
}

/*!
 *  Constructor
 */
Kent_EccTrans::Kent_EccTrans(double kappa, double ecc) : kappa(kappa), ecc(ecc)
{
  mu = XAXIS;
  major_axis = YAXIS;
  minor_axis = ZAXIS;
  computed = UNSET;
}

Kent_EccTrans::Kent_EccTrans(
  double psi, double alpha, double eta, double kappa, double ecc
) : psi(psi), alpha(alpha), eta(eta), kappa(kappa), ecc(ecc)
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
Kent_EccTrans Kent_EccTrans::operator=(const Kent_EccTrans &source)
{
  if (this != &source) {
    mu = source.mu;
    major_axis = source.major_axis;
    minor_axis = source.minor_axis;
    psi = source.psi;
    alpha = source.alpha;
    eta = source.eta;
    kappa = source.kappa;
    ecc = source.ecc;
    constants = source.constants;
    computed = source.computed;
  }
  return *this;
}

/*!
 *  Normalization constant
 */
double Kent_EccTrans::computeLogNormalizationConstant(void)
{
  double beta = 0.5 * kappa * ecc;
  constants.log_c = computeSeriesSum(kappa,beta,0.5);
  computed = SET;
  return constants.log_c;
}

double Kent_EccTrans::computeSeriesSum(double k, double b, double d)
{
  double ex = ecc; 
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

double Kent_EccTrans::log_density(Vector &x)
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

  double beta = 0.5 * kappa * ecc;
  double ans = -log_norm + kappa * c1 + beta * c2;
  return ans;
}

double Kent_EccTrans::computeNegativeLogLikelihood(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeNegativeLogLikelihood(sample_mean,S,data.size());
}

double Kent_EccTrans::computeNegativeLogLikelihood(Vector &sample_mean, Matrix &S, double N)
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

  double beta = 0.5 * kappa * ecc;
  double ans = N * log_norm - kappa * c1 - beta * c2;
  return ans;
}

double Kent_EccTrans::computeLogPriorProbability()
{
  double log_prior_axes = computeLogPriorAxes();
  double log_prior_scale = computeLogPriorScale();
  double log_factor = log(kappa) - log(2);  // Jacobian
  double log_joint_prior = log_prior_axes + log_prior_scale + log_factor;
  assert(!boost::math::isnan(log_joint_prior));
  return log_joint_prior;
}

double Kent_EccTrans::computeLogPriorAxes()
{
  double angle = alpha;
  if (angle < TOLERANCE) angle = TOLERANCE;
  if (fabs(angle-PI) < TOLERANCE) angle = PI-TOLERANCE;
  
  double log_prior = 0;
  log_prior += (log(sin(angle)));
  log_prior -= (log(4) + 2*log(PI));
  return log_prior;
}

double Kent_EccTrans::computeLogPriorScale()
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

  return log_prior;
}

Vector Kent_EccTrans::Mean()
{
  return mu;
}

Vector Kent_EccTrans::MajorAxis()
{
  return major_axis;
}

Vector Kent_EccTrans::MinorAxis()
{
  return minor_axis;
}

double Kent_EccTrans::Kappa()
{
  return kappa;
}

double Kent_EccTrans::Beta()
{
  return 0.5 * kappa * ecc;
}

double Kent_EccTrans::eccentricity()
{
  return ecc;
}

void Kent_EccTrans::printParameters(ostream &os)
{
  os << "[mu]: "; print(os,mu,3);
  os << "\t[mj]: "; print(os,major_axis,3);
  os << "\t[mi]: "; print(os,minor_axis,3);
  os << "\t[kappa]: " << fixed << setprecision(3) << kappa << "\t\t";
  os << "\t[beta]: " << fixed << setprecision(3) << 0.5 * kappa * ecc << "\t\t";
  os << "\t[ecc]: " << fixed << setprecision(3) << ecc << endl;
}

