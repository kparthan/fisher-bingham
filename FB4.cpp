#include "FB4.h"
#include "vMF.h"
#include "Normal.h"
#include "Support.h"

extern std::vector<long double> XAXIS,YAXIS,ZAXIS;

/*!
 *  Null constructor
 */
FB4::FB4() : kappa(1), gamma(1) 
{
  mu = ZAXIS;
  major_axis = XAXIS;
  minor_axis = YAXIS;
}

/*!
 *  Constructor
 */
FB4::FB4(long double kappa, long double gamma) : kappa(kappa), gamma(gamma)
{
  mu = ZAXIS;
  major_axis = XAXIS;
  minor_axis = YAXIS;
}

/*!
 *  Constructor
 */
FB4::FB4(std::vector<long double> &mu, std::vector<long double> &major_axis, 
         std::vector<long double> &minor_axis, long double kappa, long double gamma) : 
         mu(mu), major_axis(major_axis), minor_axis(minor_axis), kappa(kappa), gamma(gamma)
{}

/*!
 *  Copy the elements
 */
FB4 FB4::operator=(const FB4 &source)
{
  if (this != &source) {
    mu = source.mu;
    major_axis = source.major_axis;
    minor_axis = source.minor_axis;
    kappa = source.kappa;
    gamma = source.gamma;
  }
  return *this;
}

/*!
 *  Normalization constant
 */
long double FB4::computeNormalizationConstant(void)
{
  long double c;
  if (gamma < 0) {
    Normal normal(0,1);
    long double x1 = (kappa - 2*gamma) / sqrt(-2*gamma);
    long double x2 = (kappa + 2*gamma) / sqrt(-2*gamma);
    long double cdf1 = normal.cumulativeDensity(x1);
    long double cdf2 = normal.cumulativeDensity(x2);
    long double a0 = cdf1 - cdf2;
    long double c1 = sqrt(-PI/gamma);
    long double tmp = (-kappa*kappa)/(4*gamma);
    c = c1 * exp(tmp) * a0;
  } else if (gamma > 0) {
    long double tmp = 2*sqrt(gamma);
    long double tau1 = (kappa+2*gamma)/tmp;
    long double tau2 = (kappa-2*gamma)/tmp;
    long double int1 = computeIntegral(tau1);
    long double int2 = computeIntegral(fabs(tau2));
    tmp = exp(kappa);
    long double c1 = tmp * int1 - (sign(tau2) * int2 / tmp);
    c = c1 * exp(gamma) / sqrt(gamma);
  } else if (gamma == 0) {  // vMF
    long double c1 = exp(kappa);
    c = (c1 - (1/c1)) / kappa;
  }
  return c;
}

std::vector<std::vector<long double> > FB4::generate(int sample_size)
{
  cout << "\nGenerating from FB4:\n";
  cout << "mean: "; print(cout,mu,3); cout << endl;
  cout << "major: "; print(cout,major_axis,3); cout << endl;
  cout << "minor: "; print(cout,minor_axis,3); cout << endl;
  cout << "Kappa: " << kappa << "; Gamma: " << gamma 
       << "; sample size: " << sample_size << endl;
  std::vector<std::vector<long double> > canonical_sample = generateCanonical(sample_size);
  matrix<long double> transformation = computeOrthogonalTransformation(mu,major_axis);
  return transform(canonical_sample,transformation);
}

std::vector<std::vector<long double> > FB4::generateCanonical(int sample_size)
{
  if (gamma != 0) {
    std::vector<long double> u = generate_u(sample_size);
    // step 2 (of FB6_0)
    std::vector<long double> phi = generate_spherical_coordinates(u);
    // step 3 (of FB6_0)
    return generate_cartesian_coordinates(u,phi);
  } else if (gamma == 0) { // vMF
    std::vector<long double> mu(3,0); mu[2] = 1;
    vMF vmf(mu,kappa);
    return vmf.generate(sample_size);
  }
}

std::vector<long double> FB4::generate_u(int sample_size)
{
  std::vector<long double> u;
  if (gamma < 0) {  // FB4-
    int best = determine_best_case();
    u = generate_FB4_minus(best,sample_size);
  } else if (gamma > 0) { // FB4+
    u = generate_FB4_plus(sample_size);
  }  
  return u;
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

/*!
 *  FB4- algorithm
 */
std::vector<long double> FB4::generate_FB4_minus(int best, int sample_size)
{
  std::vector<long double> u(sample_size,0);
  if (best == 1) {
    // step 0
    long double sigma_inv = sqrt(-2*gamma);
    long double mu1 = (kappa+2*gamma)/sigma_inv;
    long double mu2 = (kappa-2*gamma)/sigma_inv;
    long double q1 = 1/sigma_inv;
    long double q2 = -kappa/(2*gamma);
    Normal normal(0,1);
    long double t,u_accept;
    std::vector<long double> tmp;
    for (int i=0; i<sample_size; i++) {
      // step 1
      tmp = normal.generate(1);
      t = tmp[0];
      // step 2
      if (t >= mu1 && t <= mu2) {
        u_accept = q1 * t + q2;
        u[i] = u_accept;
      } else {
        i--;
      }
    } // for loop ...
  } else if (best == 2 || best == 3) {
    long double r,q1,q2,s1,s2,u1,tmp1,tmp2;
    // step 0
    if (best == 2) {
      r = kappa;
    } else if (best == 3) {
      r = kappa + 2 * gamma;
    }
    q1 = exp(r);
    q2 = 1/q1;
    for (int i=0; i<sample_size; i++) {
      // step 1
      s1 = rand() / (long double)RAND_MAX;
      tmp1 = q1 * s1 + q2 * (1-s1);
      u1 = log(tmp1)/r;
      // step 2
      s2 = rand() / (long double)RAND_MAX;
      tmp1 = exp(gamma*u1*u1);
      tmp2 = exp(gamma*(1-u1)*(1-u1));
      if (s2 <= tmp1 && s2 <= tmp2) {
        u[i] = u1;
      } else {
        i--;
      }
    } // for loop ...
  } else {
    cout << "Error: best case = " << best << endl;
    exit(1);
  } 
  return u;
}

/*!
 *  FB4+ algorithm
 */
std::vector<long double> FB4::generate_FB4_plus(int sample_size)
{
  std::vector<long double> u(sample_size,0);

  long double r1,r2,m1,m2,n1,n2,lambda,num,denom,p1,s1,s2,s3,r,u1,q1,q2,tmp;
  // step 0
  r1 = kappa+gamma; r2 = kappa-gamma;
  m1 = exp(r1); m2 = 1/m1;
  n1 = exp(r2), n2 = 1/n1;
  lambda = 1+exp(-2*gamma);
  num = r2 * sinhl(r1);
  denom = num + r1 * sinhl(r2);
  p1 = num/denom;
  for(int i=0; i<sample_size; i++) {
    // step 1
    s1 = rand()/(long double)RAND_MAX;
    if (s1 <= p1) {
      r = r1; q1 = m1; q2 = m2;
    } else if (s1 > p1) {
      r = r2; q1 = n1; q2 = n2;
    }
    // step 2
    s2 = rand()/(long double)RAND_MAX;
    tmp = q1 * s2 + q2 * (1-s2);
    u1 = log(tmp)/r;
    // step 3
    s3 = rand()/(long double)RAND_MAX;
    num = lambda * exp(gamma*u1*u1);
    denom = 2*coshl(gamma*u1);
    tmp = num/denom;
    if (s3 <= tmp) {
      u[i] = u1;
    } else {
      i--;
    }
  }
  return u;
}

/*!
 *  step 2 (of FB6_0) (actually 2* for FB4)
 */
std::vector<long double> FB4::generate_spherical_coordinates(std::vector<long double> &u)
{
  std::vector<long double> phi(u.size(),0);
  long double s;
  for (int i=0; i<u.size(); i++) {
    s = rand()/(long double)RAND_MAX;
    phi[i] = 2 * PI * s;
  }
  return phi;
}

/*!
 *  step 3 (of FB6_0) : transformation to (unit) Cartesian coordinates
 */
std::vector<std::vector<long double> > FB4::generate_cartesian_coordinates(
  std::vector<long double> &u, 
  std::vector<long double> &phi
) {
  std::vector<std::vector<long double> > coordinates(u.size());
  std::vector<long double> x(3,0);
  long double tmp;
  for (int i=0; i<u.size(); i++) {
    tmp = sqrt(1-u[i]*u[i]);
    x[0] = tmp * cos(phi[i]);
    x[1] = tmp * sin(phi[i]);
    x[2] = u[i];
    coordinates[i] = x;
  }
  return coordinates;
}

