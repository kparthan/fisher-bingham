#ifndef KENT_H
#define KENT_H

#include "Header.h"
#include "Support.h"

class Kent  // FB5
{
  friend class Test;

  private:
    Vector mu,major_axis,minor_axis;

    long double alpha,eta,psi,delta;
    
    long double kappa,beta; // gamma = 0

    struct Constants {
      long double log_c,log_cb,log_ck,log_ckk,log_ckb,log_cbb;
      long double ck_c,ckk_c,cb_c,cbb_c,ckb_c;
      Vector E_x;
      Matrix E_xx;
      Matrix R,Rt;  // R: standard -> current orientation
    } constants;

    int computed;

  public:
    Kent();

    Kent(long double, long double);

    Kent(Vector &, Vector &, Vector &, long double, long double);

    Kent(long double, long double, long double, long double, long double, long double);
 
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

    long double computeNegativeLogLikelihood(std::vector<Vector> &);

    long double computeNegativeLogLikelihood(Vector &, Matrix &);

    long double computeLogPriorProbability();

    long double computeLogFisherInformation();

    struct Estimates computeMomentEstimates(std::vector<Vector> &);

    struct Estimates computeMomentEstimates(Vector &, Matrix &);

    struct Estimates computeMLEstimates(std::vector<Vector> &);

    struct Estimates computeMLEstimates(Vector &, Matrix &);

    struct Estimates computeMMLEstimates(Vector &, Matrix &);

    Vector Mean();

    Vector MajorAxis();

    Vector MinorAxis();

    long double Kappa();

    long double Beta();

    long double computeKLDivergence(Kent &);
};

#endif

