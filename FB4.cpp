#include "FB4.h"
#include "vMF.h"
#include "Normal.h"

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
    kappa = source.kappa;
    gamma = source.gamma;
    mu = source.mu;
    major_axis = source.major_axis;
    minor_axis = source.minor_axis;
  }
  return *this;
}

std::vector<std::vector<long double> > FB4::generate(int sample_size)
{
  std::vector<long double> u;
  if (gamma < 0) {
    int best = determine_best_case();
    u = generate_FB4_minus(best,sample_size);
  } else if (gamma > 0) {
  } else if (gamma == 0) {  // vMF
    std::vector<long double> mu(3,0); mu[2] = 1;
    vMF vmf(mu,kappa);
    return vmf.generate(sample_size);
  }
  // step 2 (of FB6_0)
  std::vector<long double> phi = generate_spherical_coordinates(u);
  // step 3 (of FB6_0)
  return generate_cartesian_coordinates(u,phi);
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
std::vector<std::vector<long double> > 
FB4::generate_cartesian_coordinates(
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

