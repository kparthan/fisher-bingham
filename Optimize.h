#ifndef OPTIMIZE
#define OPTIMIZE

#include <dlib/optimization.h>
#include "Kent.h"

using namespace dlib;

typedef dlib::matrix<double,0,1> column_vector;

class MomentObjectiveFunction
{
  private:
    long double C1,C2;

  public:
    MomentObjectiveFunction(Vector &m0, Vector &m1, Vector &m2,
                            Vector &sample_mean, Matrix &S) {
      C1 = computeDotProduct(sample_mean,m0);
      Vector x = prod(S,m1);
      long double mj = computeDotProduct(m1,x);
      x = prod(S,m2);
      long double mi = computeDotProduct(m2,x);
      C2 = mj - mi;
    }

    /*!
     *  minimize function: log c(k,b) - k * a - b * delta
     */
    double operator() (const column_vector& x) const {
      const double k = x(0);
      const double b = x(1);

      Kent kent(k,b);
      long double log_norm = kent.computeLogNormalizationConstant(k,b);
      double fval = log_norm - k * C1 - b * C2;
      return fval;
    }
}; 

class MaximumLikelihoodObjectiveFunction
{
  private:
    Vector sample_mean;

    Matrix S;

  public:
    MaximumLikelihoodObjectiveFunction(Vector &sample_mean, Matrix &S) :
                                       sample_mean(sample_mean), S(S)
    {}

    /*!
     *  minimize function: log c(k,b) - k  (m0' x) - b (mj' xx' mj -  mi xx' mi)
     *  x: sample mean
     *  xx' : dispersion matrix (S)
     *  k,b,m0,mj,mi are parameters
     */
    double operator() (const column_vector& x) const {
      const double alpha = x(0);
      const double eta = x(1);
      const double psi = x(2);
      const double k = x(3);
      const double b = x(4);

      Vector m0(3,0),m1(3,0),m2,spherical(3,0);
      spherical[0] = 1;
      // compute m0
      spherical[1] = alpha; spherical[2] = eta;
      spherical2cartesian(spherical,m0);
      // compute m1
      spherical[1] = psi;
      // cos(delta-eta) = -cot(alpha) cot(psi)
      double tmp = -1/(tan(alpha) * tan(psi));
      const double delta = eta + acos(tmp);
      spherical[2] = delta;
      spherical2cartesian(spherical,m1);
      // compute m2: m0 X m1
      m2 = crossProduct(m0,m1);

      Kent kent(m0,m1,m2,k,b);
      Vector sample_mean1 = sample_mean; Matrix S1 = S;
      double fval = kent.computeNegativeLogLikelihood(sample_mean1,S1);
      return fval;
    }
};

class Optimize
{
  private:
    Vector mean,major,minor; 

    double alpha,eta,psi,delta;

    double kappa,beta;

  public:
    Optimize();

    void initialize(Vector &, Vector &, Vector &, long double, long double);

    void computeMomentEstimates(Vector &, Matrix &, struct Estimates &);

    void computeMLEstimates(Vector &, Matrix &, struct Estimates &);

    column_vector minimize(Vector &, Matrix &, int, int);
};

#endif

