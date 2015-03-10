#ifndef OPTIMIZE
#define OPTIMIZE

#include <nlopt.hpp>
#include "Kent.h"
#include "Kent_EccTrans.h"
#include "Kent_UnifTrans.h"
#include "Support.h"

extern int MOMENT_FAIL,MLE_FAIL,MAP_FAIL,MML_FAIL;
extern bool FAIL_STATUS;

// MOMENT
class MomentObjectiveFunction
{
  private:
    double N;

    double C1,C2;

  public:
    MomentObjectiveFunction(Vector &m0, Vector &m1, Vector &m2,
                            Vector &sample_mean, Matrix &S, double sample_size) {
      N = sample_size;
      C1 = computeDotProduct(sample_mean,m0) / sample_size;

      double mj = prod_xMy(m1,S,m1);
      double mi = prod_xMy(m2,S,m2);
      C2 = (mj - mi) / sample_size;
    }

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        if(boost::math::isnan(x[0]) || boost::math::isnan(x[1])) { // return this sub-optimal state
          if (FAIL_STATUS == 0) {
            MOMENT_FAIL++;
            FAIL_STATUS = 1;
          }
          return -HUGE_VAL;
        } // if() ends ...
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
      double log_norm = kent.computeLogNormalizationConstant();
      double fval = log_norm - k * C1 - b * C2
                    - 2 * N * log(AOM);
      //cout << "k: " << k << "\tb: " << b << "\tfval: " << fval << endl;
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
        for (int i=0; i<x.size(); i++) {
          if(boost::math::isnan(x[i])) {  // return this sub-optimal state 
            if (FAIL_STATUS == 0) {
              MLE_FAIL++;
              FAIL_STATUS = 1;
            }
            return -HUGE_VAL;
          } // if() ends ...
        } // for() ends ...
        return (*reinterpret_cast<MaximumLikelihoodObjectiveFunction*>(data))(x, grad); 
    }

    /*!
     *  minimize function: N * log c(k,b) - k  (m0' x) - b (mj' xx' mj -  mi xx' mi)
     *  \sum_x: sample mean
     *  \sum_xx' : dispersion matrix (S)
     *  k,b,m0,mj,mi are parameters
     */
    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      double psi = x[0];
      double alpha = x[1];
      double eta = x[2];
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
        for (int i=0; i<x.size(); i++) {
          if(boost::math::isnan(x[i])) {  // return this sub-optimal state 
            if (FAIL_STATUS == 0) {
              MAP_FAIL++;
              FAIL_STATUS = 1;
            }
            return -HUGE_VAL;
          } // if() ends ...
        } // for() ends ...
        return (*reinterpret_cast<MAPObjectiveFunction*>(data))(x, grad); 
    }

    /*!
     *  minimize function: N * log c(k,b) - k  (m0' x) - b (mj' xx' mj -  mi xx' mi)
     *  \sum_x: sample mean
     *  \sum_xx' : dispersion matrix (S)
     *  k,b,m0,mj,mi are parameters
     */
    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      double psi = x[0];
      double alpha = x[1];
      double eta = x[2];
      double k = x[3];
      double b = x[4];

      Kent kent(psi,alpha,eta,k,b);
      double log_prior = kent.computeLogPriorProbability();// - log(sin(alpha));
      double fval = -log_prior + kent.computeNegativeLogLikelihood(sample_mean,S,N)
                    - 2 * N * log(AOM);
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
          //assert(!boost::math::isnan(x[i]));
          if(boost::math::isnan(x[i])) {  // return this sub-optimal state 
            if (FAIL_STATUS == 0) {
              FAIL_STATUS = 1;
            }
            return 0;
          } // if() ends ...
        } // for() ends ...
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
      double psi = x[0];
      double alpha = x[1];
      double eta = x[2];
      double k = x[3];
      double b = x[4];

      Kent kent(psi,alpha,eta,k,b);
      double log_prior = kent.computeLogPriorProbability();
      double log_fisher = kent.computeLogFisherInformation(N);
      double part1 = const_lattk - log_prior + 0.5 * log_fisher;
      if (part1 < 0) {
        part1 = 0; cout << "Part 1 is negative ...";
        //MML_FAIL = 1;
      }
      double part2 = kent.computeNegativeLogLikelihood(sample_mean,S,N) + 2.5
                     - 2 * N * log(AOM);
      double fval = part1 + part2;
      //assert(!boost::math::isnan(fval));
      return fval;
    }
};

// MAP_ECCENTRICITY_TRANSFORM 
class MAPObjectiveFunction_EccTrans
{
  private:
    Vector sample_mean;

    Matrix S;

    double N;

  public:
    MAPObjectiveFunction_EccTrans(
      Vector &sample_mean, Matrix &S, double sample_size
    ) : sample_mean(sample_mean), S(S), N(sample_size)
    {}

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        for (int i=0; i<x.size(); i++) {
          if(boost::math::isnan(x[i])) {  // return this sub-optimal state 
            if (FAIL_STATUS == 0) {
              MAP_FAIL++;
              FAIL_STATUS = 1;
            }
            return -HUGE_VAL;
          } // if() ends ...
        } // for() ends ...
        return (*reinterpret_cast<MAPObjectiveFunction_EccTrans*>(data))(x, grad); 
    }

    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      double psi = x[0];
      double alpha = x[1];
      double eta = x[2];
      double k = x[3];
      double e = x[4];

      Kent_EccTrans kent(psi,alpha,eta,k,e);
      double log_prior = kent.computeLogPriorProbability();
      //log_prior = 0;
      double fval = -log_prior + kent.computeNegativeLogLikelihood(sample_mean,S,N)
                    - 2 * N * log(AOM);
      return fval;
    }
};

// MAP_UNIFORM_TRANSFORM 
class MAPObjectiveFunction_UnifTrans
{
  private:
    Vector sample_mean;

    Matrix S;

    double N;

  public:
    MAPObjectiveFunction_UnifTrans(
      Vector &sample_mean, Matrix &S, double sample_size
    ) : sample_mean(sample_mean), S(S), N(sample_size)
    {}

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        for (int i=0; i<x.size(); i++) {
          if(boost::math::isnan(x[i])) {  // return this sub-optimal state 
            if (FAIL_STATUS == 0) {
              MAP_FAIL++;
              FAIL_STATUS = 1;
            }
            return -HUGE_VAL;
          } // if() ends ...
        } // for() ends ...
        return (*reinterpret_cast<MAPObjectiveFunction_UnifTrans*>(data))(x, grad); 
    }

    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      double z1 = x[0];
      double z2 = x[1];
      double z3 = x[2];
      double z4 = x[3];
      double z5 = x[4];

      Kent_UnifTrans kent(z1,z2,z3,z4,z5);
      double fval = kent.computeNegativeLogLikelihood(sample_mean,S,N)
                    - 2 * N * log(AOM);
      return fval;
    }
};

class Optimize
{
  private:
    int estimation;

    double N;

    Vector mean,major,minor; 

    double psi,alpha,eta;

    double kappa,beta;

  public:
    Optimize(string);

    void initialize(double, Vector &, Vector &, Vector &, double, double);

    void computeEstimates(Vector &, Matrix &, struct Estimates &);

    void finalize(std::vector<double> &, struct Estimates &);

    void validate_scale(double &, double &);

    std::vector<double> minimize(Vector &, Matrix &, int);
};

#endif

