#ifndef OPTIMIZE
#define OPTIMIZE

#include <nlopt.hpp>
#include "Kent.h"
#include "Support.h"

// MOMENT
class MomentObjectiveFunction
{
  private:
    long double C1,C2;

  public:
    MomentObjectiveFunction(Vector &m0, Vector &m1, Vector &m2,
                            Vector &sample_mean, Matrix &S, double sample_size) {
      C1 = computeDotProduct(sample_mean,m0) / sample_size;

      long double mj = prod_xMy(m1,S,m1);
      long double mi = prod_xMy(m2,S,m2);
      C2 = (mj - mi) / sample_size;
    }

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        return (*reinterpret_cast<MomentObjectiveFunction*>(data))(x, grad); 
    }

    /*!
     *  minimize function: log c(k,b) - k * a - b * delta
     */
    double operator() (const std::vector<double> &x, std::vector<double> &grad)
    {
      double k = x[0];
      double b = x[1];

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

    double N;

  public:
    MaximumLikelihoodObjectiveFunction(
      Vector &sample_mean, Matrix &S, double sample_size
    ) : sample_mean(sample_mean), S(S), N(sample_size)
    {}

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        return (*reinterpret_cast<MaximumLikelihoodObjectiveFunction*>(data))(x, grad); 
    }

    /*!
     *  minimize function: N * log c(k,b) - k  (m0' x) - b (mj' xx' mj -  mi xx' mi)
     *  \sum_x: sample mean
     *  \sum_xx' : dispersion matrix (S)
     *  k,b,m0,mj,mi are parameters
     */
    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      double alpha = x[0];
      double eta = x[1];
      double psi = x[2];
      double k = x[3];
      double b = x[4];

      Kent kent(psi,alpha,eta,k,b);
      double fval = kent.computeNegativeLogLikelihood(sample_mean,S,N)
                    - 2 * N * log(AOM);
      return fval;
    }
};

// MAP 
class MAPObjectiveFunction
{
  private:
    Vector sample_mean;

    Matrix S;

    double N;

  public:
    MAPObjectiveFunction(
      Vector &sample_mean, Matrix &S, double sample_size
    ) : sample_mean(sample_mean), S(S), N(sample_size)
    {}

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        return (*reinterpret_cast<MAPObjectiveFunction*>(data))(x, grad); 
    }

    /*!
     *  minimize function: N * log c(k,b) - k  (m0' x) - b (mj' xx' mj -  mi xx' mi)
     *  \sum_x: sample mean
     *  \sum_xx' : dispersion matrix (S)
     *  k,b,m0,mj,mi are parameters
     */
    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      double alpha = x[0];
      double eta = x[1];
      double psi = x[2];
      double k = x[3];
      double b = x[4];

      Kent kent(psi,alpha,eta,k,b);
      long double log_prior = kent.computeLogPriorProbability();// - log(sin(alpha));
      double fval = -log_prior + kent.computeNegativeLogLikelihood(sample_mean,S,N)
                    - 2 * N * log(AOM);
      return fval;
    }
};

// MML Scale
class MMLObjectiveFunctionScale
{
  private:
    Vector sample_mean;

    Matrix S;

    double N;

    double psi,alpha,eta,kappa_init,beta_init;

    double k2,const_lattk;

    double log_prior_axes,log_fisher_axes;

  public:
    MMLObjectiveFunctionScale(
      double psi, double alpha, double eta, double k, double b,
      Vector &sample_mean, Matrix &S, double sample_size
    ) : psi(psi), alpha(alpha), eta(eta), kappa_init(k), beta_init(b),
        sample_mean(sample_mean), S(S), N(sample_size)
    {
      /*k2 = 0.08019;
      Kent kent(psi,alpha,eta,kappa_init,beta_init);
      kent.computeExpectation();
      log_prior_axes = kent.computeLogPriorAxes();
      log_fisher_axes = kent.computeLogFisherAxes();*/
      const_lattk = -6.455;
    }

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        return (*reinterpret_cast<MMLObjectiveFunctionScale*>(data))(x, grad); 
    }

    /*!
     *  minimize function: (d/2) log(kd) - log h + 0.5 log (det(fisher)) +
     *          N * log c(k,b) - k  (m0' x) - b (mj' xx' mj -  mi xx' mi) + (d/2)
     *  d: 2
     *  x: sample mean
     *  xx' : dispersion matrix (S)
     *  k,b,m0,mj,mi are parameters
     */
    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      double k = fabs(x[0]);
      double b = fabs(x[1]);

      Kent kent(psi,alpha,eta,k,b);
      long double log_prior = kent.computeLogPriorProbability();
      long double log_fisher = kent.computeLogFisherInformation(N);
      double part1 = const_lattk - log_prior + 0.5 * log_fisher;
      double part2 = kent.computeNegativeLogLikelihood(sample_mean,S,N) + 2.5
                     - 2 * N * log(AOM);
      double fval = part1 + part2;
      //assert(!boost::math::isnan(fval));
      return fval;
    }
};

// MML
class MMLObjectiveFunction
{
  private:
    Vector sample_mean;

    Matrix S;

    double N;

    double const_lattk;

  public:
    MMLObjectiveFunction(
      Vector &sample_mean, Matrix &S, double sample_size
    ) : sample_mean(sample_mean), S(S), N(sample_size)
    {
      //k5 = 0.0756226;
      const_lattk = -6.455;
    }

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        for (int i=0; i<x.size(); i++) {
          assert(!boost::math::isnan(x[i]));
        }
        return (*reinterpret_cast<MMLObjectiveFunction*>(data))(x, grad); 
    }

    /*!
     *  minimize function: (d/2) log(kd) - log h + 0.5 log (det(fisher)) +
     *          N * log c(k,b) - k  (m0' x) - b (mj' xx' mj -  mi xx' mi) + (d/2)
     *  d: 5
     *  x: sample mean
     *  xx' : dispersion matrix (S)
     *  k,b,m0,mj,mi are parameters
     */
    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      double alpha = x[0];
      double eta = x[1];
      double psi = x[2];
      double k = fabs(x[3]);
      double b = fabs(x[4]);

      Kent kent(psi,alpha,eta,k,b);
      long double log_prior = kent.computeLogPriorProbability();
      long double log_fisher = kent.computeLogFisherInformation(N);
      double part1 = const_lattk - log_prior + 0.5 * log_fisher;
      double part2 = kent.computeNegativeLogLikelihood(sample_mean,S,N) + 2.5
                     - 2 * N * log(AOM);
      double fval = part1 + part2;
      //assert(!boost::math::isnan(fval));
      return fval;
    }
};

class Optimize
{
  private:
    int estimation;

    double N;

    Vector mean,major,minor; 

    long double psi,alpha,eta;

    double kappa,beta;

  public:
    Optimize(string);

    void initialize(double, Vector &, Vector &, Vector &, long double, long double);

    void computeEstimates(Vector &, Matrix &, struct Estimates &);

    void finalize(std::vector<double> &, struct Estimates &);

    void validate_scale(long double &, long double &);

    std::vector<double> minimize(Vector &, Matrix &, int);
};

#endif

