#ifndef OPTIMIZE2
#define OPTIMIZE2

#include <nlopt.hpp>
#include "vMF.h"
#include "Support.h"

extern int MLE_FAIL,MAP_FAIL,MML2_FAIL;
extern bool FAIL_STATUS;

// MLE 
class MaximumLikelihoodObjectiveFunction_vMF
{
  private:
    double N,R;

  public:
    MaximumLikelihoodObjectiveFunction_vMF(double sample_size, double magnitude) :
                                           N(sample_size), R(magnitude)
    {}

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        if(boost::math::isnan(x[0])) {  // return this sub-optimal state 
          if (FAIL_STATUS == 0) {
            MLE_FAIL++;
            FAIL_STATUS = 1;
          }
          return -HUGE_VAL;
        } // for() ends ...
        return (*reinterpret_cast<MaximumLikelihoodObjectiveFunction_vMF*>(data))(x, grad); 
    }

    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      double k = x[0];

      vMF vmf(k);
      double log_norm = vmf.getLogNormalizationConstant();
      double fval = (N * log_norm) + (k * R) + (2 * N * log(AOM));
      return -fval;
    }
};

// MAP 
class MAPObjectiveFunction_vMF
{
  private:
    double N,R;

    Vector mean;

  public:
    MAPObjectiveFunction_vMF(double N, double R, Vector &mean) : N(N), R(R), mean(mean)
    {}

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        if(boost::math::isnan(x[0])) {  // return this sub-optimal state 
          if (FAIL_STATUS == 0) {
            MAP_FAIL++;
            FAIL_STATUS = 1;
          }
          return -HUGE_VAL;
        } // if() ends ...
        return (*reinterpret_cast<MAPObjectiveFunction_vMF*>(data))(x, grad); 
    }

    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      double k = x[0];

      vMF vmf(mean,k);
      double log_prior = vmf.computeLogPriorProbability(); 
      double fval = -log_prior + vmf.computeNegativeLogLikelihood(R,N)
                    - 2 * N * log(AOM);
      return fval;
    }
};

// MML
class MMLObjectiveFunction_vMF
{
  private:
    double N,R;

    Vector mean;

    double const_lattk;

  public:
    MMLObjectiveFunction_vMF(double N, double R, Vector &mean) : N(N), R(R), mean(mean)
    {
      const_lattk = -3.816;
    }

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        if(boost::math::isnan(x[0])) {  // return this sub-optimal state 
          if (FAIL_STATUS == 0) {
            MML2_FAIL++;
            FAIL_STATUS = 1;
          }
          return 0;
        } // if() ends ...
        return (*reinterpret_cast<MMLObjectiveFunction_vMF*>(data))(x, grad); 
    }

    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      double k = x[0];

      vMF vmf(mean,k);
      double log_prior = vmf.computeLogPriorProbability();
      double log_fisher = vmf.computeLogFisherInformation(N);
      double part1 = const_lattk - log_prior + 0.5 * log_fisher;
      double part2 = vmf.computeNegativeLogLikelihood(R,N) + 1.5
                     - 2 * N * log(AOM);
      double fval = part1 + part2;
      //assert(!boost::math::isnan(fval));
      return fval;
    }
};

class Optimize2
{
  private:
    int estimation;

    double N,R;

    Vector mean;

    double kappa;

  public:
    Optimize2(string);

    void initialize(double, double, Vector &, double);

    void computeEstimates(struct Estimates_vMF &);
};

#endif

