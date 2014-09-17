#ifndef KENT_H
#define KENT_H

#include "Header.h"

class Kent  // FB5
{
  friend class Test;
  friend class Experiments;

  private:
    Vector mu,major_axis,minor_axis;

    long double psi,alpha,eta;
    
    long double kappa,beta; // gamma = 0

    struct Constants {
      long double log_c,log_cb,log_ck,log_ckk,log_ckb,log_cbb;
      long double ck_c,ckk_c,cb_c,cbb_c,ckb_c;
      long double lambda1,lambda2,lambda3;
      Vector E_x,kappa_E_x;
      Matrix E_xx,beta_E_xx;
      Matrix R,Rt;  // R: standard -> current orientation
      Matrix fisher_axes;
    } constants;

    int computed;

  public:
    Kent();

    Kent(long double, long double);

    Kent(Vector &, Vector &, Vector &, long double, long double);

    Kent(long double, long double, long double, long double, long double);
 
    Kent operator=(const Kent &);

    std::vector<Vector> generate(int);

    std::vector<Vector> generateCanonical(int);

    long double eccentricity();

    struct Constants getConstants();

    long double computeLogNormalizationConstant();

    long double log_dc_dk();

    long double log_d2c_dk2();

    long double computeSeriesSum(long double, long double, long double);

    long double log_dc_db();

    long double log_d2c_dkdb();

    long double computeSeriesSum2(long double, long double, long double);

    long double log_d2c_db2();

    void computeConstants();

    void computeExpectation();

    long double computeLogFisherAxes();

    long double computeLogFisherScale();

    long double log_density(Vector &);

    long double computeNegativeLogLikelihood(std::vector<Vector> &);

    long double computeNegativeLogLikelihood(Vector &, Matrix &, long double);

    long double computeLogPriorProbability();

    long double computeLogPriorAxes();

    long double computeLogPriorScale();

    long double computeLogFisherInformation();

    long double computeLogFisherInformation(long double);

    void computeAllEstimators(std::vector<Vector> &);

    void computeAllEstimators(Vector &, Matrix &, long double);

    void computeAllEstimators(std::vector<Vector > &, std::vector<struct Estimates> &);

    struct Estimates computeMomentEstimates(std::vector<Vector> &);

    struct Estimates computeMomentEstimates(Vector &, Matrix &, long double);

    struct Estimates computeMLEstimates(std::vector<Vector> &, string);

    struct Estimates computeMLEstimates(Vector &, Matrix &, long double, string);

    struct Estimates computeMMLEstimates(std::vector<Vector> &);

    struct Estimates computeMMLEstimates(Vector &, Matrix &, long double);

    void estimateParameters(std::vector<Vector > &, Vector &);

    void updateParameters(struct Estimates &);

    Vector Mean();

    Vector MajorAxis();

    Vector MinorAxis();

    long double Kappa();

    long double Beta();

    long double computeKLDivergence(Kent &);

    long double computeKLDivergence(struct Estimates &);

    long double computeMessageLength(Vector &, Matrix &, long double);

    long double computeMessageLength(struct Estimates &, Vector &, Matrix &, long double);

    void printParameters(ostream &);
};

#endif

