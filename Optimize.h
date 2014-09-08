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

// MLE 
class MaximumLikelihoodObjectiveFunction
{
  private:
    Vector sample_mean;

    Matrix S;

    int N;

  public:
    MaximumLikelihoodObjectiveFunction(
      Vector &sample_mean, Matrix &S, int sample_size, column_vector &sp
    ) : sample_mean(sample_mean), S(S), N(sample_size)
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

      Kent kent(psi,alpha,eta,k,b);
      Vector sample_mean1 = sample_mean; Matrix S1 = S;
      double fval = kent.computeNegativeLogLikelihood(sample_mean1,S1,N)
                    - 2 * N * log(AOM);
      assert(!boost::math::isnan(fval));
      return fval;
    }
};


class MMLObjectiveFunctionScale
{
  private:
    Vector sample_mean;

    Matrix S;

    int N;

    double psi,alpha,eta,kappa_init,beta_init;

    double k2;

    double log_prior_axes,log_fisher_axes;

  public:
    MMLObjectiveFunctionScale(
      double psi, double alpha, double eta, double k, double b,
      Vector &sample_mean, Matrix &S, int sample_size
    ) : psi(psi), alpha(alpha), eta(eta), kappa_init(k), beta_init(b),
        sample_mean(sample_mean), S(S), N(sample_size)
    {
      k2 = 0.08019;
      Kent kent(psi,alpha,eta,kappa_init,beta_init);
      kent.computeExpectation();
      log_prior_axes = kent.computeLogPriorAxes();
      log_fisher_axes = kent.computeLogFisherAxes();
    }

    /*!
     *  minimize function: (d/2) log(kd) - log h + 0.5 log (det(fisher)) +
     *          N * log c(k,b) - k  (m0' x) - b (mj' xx' mj -  mi xx' mi) + (d/2)
     *  d: 2
     *  x: sample mean
     *  xx' : dispersion matrix (S)
     *  k,b,m0,mj,mi are parameters
     */
    double operator() (const column_vector& x) const {
      double k = x(0);
      double b = x(1);

      Kent kent(psi,alpha,eta,k,b);
      kent.computeExpectation();
      long double log_prior_scale = kent.computeLogPriorScale();
      long double log_fisher_scale = kent.computeLogFisherScale();
      long double log_fisher = log_fisher_axes + log_fisher_scale + 5 * log(N);

      Vector sample_mean1 = sample_mean; Matrix S1 = S;
      double part1 = log(k2) - log_prior_axes - log_prior_scale + 0.5 * log_fisher;
      double part2 = kent.computeNegativeLogLikelihood(sample_mean1,S1,N) + 1
                     - 2 * N * log(AOM);
      double fval = part1 + part2;
      assert(!boost::math::isnan(fval));
      return fval;
    }
};

//class MMLObjectiveFunction
//{
//  private:
//    Vector sample_mean;
//
//    Matrix S;
//
//    int N;
//
//    double alpha_init,eta_init,psi_init,delta_init;
//
//    double kd;
//
//  public:
//    MMLObjectiveFunction(
//      Vector &sample_mean, Matrix &S, int sample_size, column_vector &sp
//    ) : sample_mean(sample_mean), S(S), N(sample_size)
//    {
//      alpha_init = sp(0);
//      eta_init = sp(1);
//      psi_init = sp(2);
//      double tmp = -1/(tan(alpha_init) * tan(psi_init));
//      delta_init = eta_init + acos(tmp);
//      kd = 1;
//    }
//
//    /*!
//     *  minimize function: (d/2) log(kd) - log h + 0.5 log (det(fisher)) +
//     *          N * log c(k,b) - k  (m0' x) - b (mj' xx' mj -  mi xx' mi) + (d/2)
//     *  d: 5
//     *  x: sample mean
//     *  xx' : dispersion matrix (S)
//     *  k,b,m0,mj,mi are parameters
//     */
//    double operator() (const column_vector& x) const {
//      double alpha = x(0);
//      double eta = x(1);
//      double psi = x(2);
//      double k = x(3);
//      double b = x(4);
//
//      Vector m0(3,0),m1(3,0),m2,spherical(3,0);
//      spherical[0] = 1;
//      // compute m0
//      spherical[1] = alpha; spherical[2] = eta;
//      spherical2cartesian(spherical,m0);
//      // compute m1
//      // cos(delta-eta) = -cot(alpha) cot(psi)
//      double tmp = -1/(tan(alpha) * tan(psi));
//      double delta;
//      if (fabs(tmp)-1 < TOLERANCE) {
//        cout << "here\n";
//        if (tmp < 0) tmp = -1;
//        else if (tmp > 0) tmp = 1;
//      }
//      if (fabs(tmp) > 1) {
//        cout << "yes, tmp: " << tmp << endl;;
//        alpha = alpha_init;
//        eta = eta_init;
//        psi = psi_init;
//        delta = delta_init;
//      } else {
//        double acos_tmp = acos(tmp);
//        delta = eta + acos_tmp;
//        if (!(delta >= eta && delta <= PI+eta)) {
//          cout << "yes2\n";
//          delta = eta - acos_tmp;
//          assert(delta >= eta && delta <= PI+eta);
//        }
//      }
//      spherical[1] = psi;
//      spherical[2] = delta;
//      spherical2cartesian(spherical,m1);
//      // compute m2: m0 X m1
//      m2 = crossProduct(m0,m1);
//
//      Kent kent(alpha,eta,psi,delta,k,b);
//      long double log_prior = kent.computeLogPriorProbability();
//      long double log_fisher = kent.computeLogFisherInformation(N);
//      Vector sample_mean1 = sample_mean; Matrix S1 = S;
//      double part1 = 2.5 * log(kd) - log_prior + 0.5 * log_fisher;
//      double part2 = kent.computeNegativeLogLikelihood(sample_mean1,S1,N) + 2.5
//                     - 2 * N * log(AOM);
//      double fval = part1 + part2;
//      assert(!boost::math::isnan(fval));
//      return fval;
//    }
//};
//
class Optimize
{
  private:
    int estimation;

    int N;

    Vector mean,major,minor; 

    long double psi,alpha,eta;

    double kappa,beta;

  public:
    Optimize(string);

    void initialize(int, Vector &, Vector &, Vector &, long double, long double);

    void computeEstimates(Vector &, Matrix &, struct Estimates &);

    void finalize(column_vector &, struct Estimates &);

    column_vector minimize(Vector &, Matrix &, int);
};

#endif

