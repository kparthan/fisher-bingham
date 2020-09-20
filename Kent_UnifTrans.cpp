#include "Kent_UnifTrans.h"
#include "Support.h"

extern Vector XAXIS,YAXIS,ZAXIS;
extern int VERBOSE,COMPUTE_KLDIV;
extern int ESTIMATION;

Kent_UnifTrans::Kent_UnifTrans(
  double z1, double z2, double z3, double z4, double z5
) : z1(z1), z2(z2), z3(z3), z4(z4), z5(z5)
{
  double psi = PI * z1;
  double alpha = acos(1 - 2 * z2);
  double eta = 2 * PI * z3;

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
Kent_UnifTrans Kent_UnifTrans::operator=(const Kent_UnifTrans &source)
{
  if (this != &source) {
    mu = source.mu;
    major_axis = source.major_axis;
    minor_axis = source.minor_axis;
    z1 = source.z1;
    z2 = source.z2;
    z3 = source.z3;
    z4 = source.z4;
    z5 = source.z5;
    constants = source.constants;
    computed = source.computed;
  }
  return *this;
}

/*!
 *  Normalization constant
 */
double Kent_UnifTrans::computeLogNormalizationConstant(void)
{
  double tmp = acos(1 - z4);
  double kappa = tan(tmp);
  double beta = 0.5 * kappa * z5;

  constants.log_c = computeSeriesSum(kappa,beta,0.5);
  computed = SET;
  return constants.log_c;
}

double Kent_UnifTrans::computeSeriesSum(double k, double b, double d)
{
  double ex = 2 * b / k; 
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

double Kent_UnifTrans::log_density(Vector &x)
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

  double tmp = acos(1 - z4);
  double kappa = tan(tmp);
  double beta = 0.5 * kappa * z5;
  double ans = -log_norm + kappa * c1 + beta * c2;
  return ans;
}

double Kent_UnifTrans::computeNegativeLogLikelihood(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  Matrix S = computeDispersionMatrix(data);
  return computeNegativeLogLikelihood(sample_mean,S,data.size());
}

double Kent_UnifTrans::computeNegativeLogLikelihood(Vector &sample_mean, Matrix &S, double N)
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

  double tmp = acos(1 - z4);
  double kappa = tan(tmp);
  double beta = 0.5 * kappa * z5;
  double ans = N * log_norm - kappa * c1 - beta * c2;
  return ans;
}

Vector Kent_UnifTrans::Mean()
{
  return mu;
}

Vector Kent_UnifTrans::MajorAxis()
{
  return major_axis;
}

Vector Kent_UnifTrans::MinorAxis()
{
  return minor_axis;
}

double Kent_UnifTrans::Kappa()
{
  double tmp = acos(1 - z4);
  double kappa = tan(tmp);
  return kappa;
}

double Kent_UnifTrans::Beta()
{
  double tmp = acos(1 - z4);
  double kappa = tan(tmp);
  double beta = 0.5 * kappa * z5;
  return beta;
}

double Kent_UnifTrans::eccentricity()
{
  double tmp = acos(1 - z4);
  double kappa = tan(tmp);
  double beta = 0.5 * kappa * z5;
  return 2 * beta / kappa;
}

