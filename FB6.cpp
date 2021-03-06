#include "FB6.h"
#include "FB4.h"
#include "vMC.h"
#include "Support.h"

extern Vector XAXIS,YAXIS,ZAXIS;

/*!
 *  Null constructor
 */
FB6::FB6() : kappa(1), beta(1), gamma(1) 
{
  mu = ZAXIS;
  major_axis = XAXIS;
  minor_axis = YAXIS;
}

/*!
 *  Constructor
 */
FB6::FB6(double kappa, double beta, double gamma) : 
        kappa(kappa), beta(beta), gamma(gamma)
{
  mu = ZAXIS;
  major_axis = XAXIS;
  minor_axis = YAXIS;
}

/*!
 *  Constructor
 */
FB6::FB6(Vector &mu, Vector &major_axis, Vector &minor_axis,
         double kappa, double beta, double gamma) :
         mu(mu), major_axis(major_axis), minor_axis(minor_axis), 
         kappa(kappa), beta(beta), gamma(gamma)
{}

/*!
 *  Copy the elements
 */
FB6 FB6::operator=(const FB6 &source)
{
  if (this != &source) {
    mu = source.mu;
    major_axis = source.major_axis;
    minor_axis = source.minor_axis;
    kappa = source.kappa;
    beta = source.beta;
    gamma = source.gamma;
  }
  return *this;
}

std::vector<Vector> FB6::generate(int sample_size)
{
  cout << "\nGenerating from FB6:\n";
  cout << "mean: "; print(cout,mu,3); cout << endl;
  cout << "major: "; print(cout,major_axis,3); cout << endl;
  cout << "minor: "; print(cout,minor_axis,3); cout << endl;
  cout << "Kappa: " << kappa << "; Beta: " << beta << "; Gamma: " << gamma 
       << "; sample size: " << sample_size << endl;
  std::vector<Vector> canonical_sample = generateCanonical(sample_size);
  Matrix transformation = computeOrthogonalTransformation(mu,major_axis);
  return transform(canonical_sample,transformation);
}

std::vector<Vector> FB6::generateCanonical(int sample_size)
{
  if (beta != 0) {  // general case
    double psi1 = gamma - beta;
    double psi2 = gamma + beta;
    FB4 g1(kappa,psi1);
    FB4 g2(kappa,psi2);
    double c1 = g1.computeNormalizationConstant();
    double c2 = g2.computeNormalizationConstant();
    double num,denom;
    Vector u(sample_size,0);
    // step 1: FB6_0
    // step 0
    denom = c1 + exp(-2*beta) * c2;
    double p2 = c1 / denom;
    double s1,s2,eta,tmp;
    Vector u1;
    for (int i=0; i<sample_size; i++) {
      // step 1
      s1 = uniform_random();
      if (s1 <= p2) {
        eta = psi1;
      } else if (s1 > p2) {
        eta = psi2;
      }
      // step 2
      FB4 fb4(kappa,eta);
      u1 = fb4.generate_u(1);
      // step 3
      s2 = uniform_random();
      tmp = beta * (1-u1[0]*u1[0]);
      num = cyl_bessel_i(0,tmp);
      denom = coshl(tmp);
      if (s2 <= num/denom) {
        u[i] = u1[0];
      } else {
        i--;
      }
    } // for loop ends ...
    // step 2: FB6_0
    Vector phi(sample_size,0);
    double p,psi,angle;
    int delta;
    Vector vmf2dmean(2,0); vmf2dmean[0] = 1; 
    std::vector<Vector> random_sample;
    Vector x(2,0);
    for (int i=0; i<sample_size; i++) {
      // generate delta from Bernoulli distribution
      p = uniform_random();
      if (p <= 0.5) {
        delta = 0;
      } else {
        delta = 1;
      }
      // generate psi from von Mises circular
      tmp = beta * (1-u[i]*u[i]);
      vMC vmc(vmf2dmean,tmp);
      vmc.generateCanonical(random_sample,1);
      //x = random_sample[0];
      /*angle = atan(fabs(x[0])/fabs(x[1]));
      if (x[0] < 0 && x[1] > 0) {
        psi = angle;
      } else if (x[0] < 0 && x[1] < 0) {
        psi = PI - angle;
      } else if (x[0] > 0 && x[1] < 0) {
        psi = PI + angle;
      } else {
        psi = 2 * PI - angle;
      }*/
      x[0] = random_sample[0][1];
      x[1] = -random_sample[0][0];
      angle = atan(fabs(x[1])/fabs(x[0]));
      if (x[0] == 0) {
        if (x[1] < 0) psi = 1.5 * PI;
        else if (x[1] >= 0) psi = PI/2;
      } else if (x[1] == 0) {
        if (x[0] < 0) psi = PI;
        else if (x[0] >= 0) psi = 0;
      } else if (x[0] < 0 && x[1] > 0) {
        psi = PI - angle;
      } else if (x[0] < 0 && x[1] < 0) {
        psi = PI + angle;
      } else if (x[0] > 0 && x[1] < 0) {
        psi = 2*PI - angle;
      } else {
        psi = angle;
      }
      phi[i] = delta*PI + psi/2;
    } // for loop ends ...
    // step 3: FB6_0
    return generate_cartesian_coordinates(u,phi);
  } else if (beta == 0) { // FB6 becomes FB4
    FB4 fb4(kappa,gamma);
    return fb4.generateCanonical(sample_size);
  }
}

/*!
 *  step 3 (of FB6_0) : transformation to (unit) Cartesian coordinates
 */
std::vector<Vector> FB6::generate_cartesian_coordinates(Vector &u, Vector &phi)
{
  std::vector<Vector> coordinates(u.size());
  Vector x(3,0);
  double tmp;
  for (int i=0; i<u.size(); i++) {
    tmp = sqrt(1-u[i]*u[i]);
    x[0] = tmp * cos(phi[i]);
    x[1] = tmp * sin(phi[i]);
    x[2] = u[i];
    coordinates[i] = x;
  }
  return coordinates;
}

