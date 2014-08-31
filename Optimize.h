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
                            Vector &sample_mean, Matrix &S, int sample_size) {
      C1 = computeDotProduct(sample_mean,m0) / sample_size;

      long double mj = prod_xMy(m1,S,m1);
      long double mi = prod_xMy(m2,S,m2);
      C2 = (mj - mi) / sample_size;
    }

    /*!
     *  minimize function: log c(k,b) - k * a - b * delta
     */
    double operator() (const column_vector& x) const {
      const double k = x(0);
      const double b = x(1);

      Kent kent(k,b);
      long double log_norm = kent.computeLogNormalizationConstant();
      double fval = log_norm - k * C1 - b * C2;
      return fval;
    }
}; 

class MaximumLikelihoodObjectiveFunction
{
  private:
    Vector sample_mean;

    Matrix S;

    int N;

    double psi_init,delta_init;

  public:
    MaximumLikelihoodObjectiveFunction(Vector &sample_mean, Matrix &S, int sample_size,
                                       double psi_init, double delta_init) : 
                                       sample_mean(sample_mean), S(S), N(sample_size), 
                                       psi_init(psi_init), delta_init(delta_init)
    {}

    /*!
     *  minimize function: N * log c(k,b) - k  (m0' x) - b (mj' xx' mj -  mi xx' mi)
     *  \sum_x: sample mean
     *  \sum_xx' : dispersion matrix (S)
     *  k,b,m0,mj,mi are parameters
     */
    double operator() (const column_vector& x) const {
      double alpha = x(0);
      double eta = x(1);
      double psi = x(2);
      double k = x(3);
      double b = x(4);

      Vector m0(3,0),m1(3,0),m2,spherical(3,0);
      spherical[0] = 1;
      // compute m0
      spherical[1] = alpha; spherical[2] = eta;
      spherical2cartesian(spherical,m0);
      // compute m1
      // cos(delta-eta) = -cot(alpha) cot(psi)
      double tmp = -1/(tan(alpha) * tan(psi));
      double delta;
      if (fabs(tmp) > 1) {
        cout << "yes\n";
        psi = psi_init;
        delta = delta_init;
      } else {
        double acos_tmp = acos(tmp);
        delta = eta + acos_tmp;
        if (!(delta >= eta && delta <= PI+eta)) {
          cout << "yes2\n";
          delta = eta - acos_tmp;
          assert(delta >= eta && delta <= PI+eta);
        }
      }
      spherical[1] = psi;
      spherical[2] = delta;
      spherical2cartesian(spherical,m1);
      // compute m2: m0 X m1
      m2 = crossProduct(m0,m1);

      Kent kent(m0,m1,m2,k,b);
      Vector sample_mean1 = sample_mean; Matrix S1 = S;
      double fval = kent.computeNegativeLogLikelihood(sample_mean1,S1,N);
      assert(!boost::math::isnan(fval));
      return fval;
    }
};

class MMLObjectiveFunction
{
  private:
    Vector sample_mean;

    Matrix S;

    int N;

    double psi_init,delta_init;

  public:
    MMLObjectiveFunction(Vector &sample_mean, Matrix &S, int sample_size, 
                         double psi_init, double delta_init) :
                         sample_mean(sample_mean), S(S), N(sample_size),
                         psi_init(psi_init), delta_init(delta_init)
    {}

    /*!
     *  minimize function: log c(k,b) - k  (m0' x) - b (mj' xx' mj -  mi xx' mi)
     *  x: sample mean
     *  xx' : dispersion matrix (S)
     *  k,b,m0,mj,mi are parameters
     */
    double operator() (const column_vector& x) const {
      double alpha = x(0);
      double eta = x(1);
      double psi = x(2);
      double k = x(3);
      double b = x(4);

      double tmp = -1/(tan(alpha) * tan(psi));
      double delta;
      if (fabs(tmp) > 1) {
        psi = psi_init;
        delta = delta_init;
      } else {
        delta = eta + acos(tmp);
      }
      Kent kent(alpha,eta,psi,delta,k,b);
      long double log_prior = kent.computeLogPriorProbability();
      long double log_fisher = kent.computeLogFisherInformation();

    }
};

class Optimize
{
  private:
    int N;

    Vector mean,major,minor; 

    double alpha,eta,psi,delta;

    double kappa,beta;

  public:
    Optimize();

    void initialize(int, Vector &, Vector &, Vector &, long double, long double);

    void finalize(column_vector &, struct Estimates &);

    void computeMomentEstimates(Vector &, Matrix &, struct Estimates &);

    void computeMLEstimates(Vector &, Matrix &, struct Estimates &);

    void computeMMLEstimates(Vector &, Matrix &, struct Estimates &);

    column_vector minimize(Vector &, Matrix &, int, int);
};

#endif

