#ifndef KENT_H
#define KENT_H

#include "Header.h"

class Kent  // FB5
{
  friend class Test;
  friend class Experiments;

  private:
    Vector mu,major_axis,minor_axis;

    double psi,alpha,eta;
    
    double kappa,beta; // gamma = 0

    struct Constants {
      double log_c,log_cb,log_ck,log_ckk,log_ckb,log_cbb;
      double ck_c,ckk_c,cb_c,cbb_c,ckb_c;
      double lambda1,lambda2,lambda3;
      Vector E_x;
      Matrix E_xx;
      Matrix R,Rt;  // R: standard -> current orientation
      Matrix fisher_axes;
    } constants;

    int computed;

  public:
    Kent();

    Kent(double, double);

    Kent(Vector &, Vector &, Vector &, double, double);

    Kent(double, double, double, double, double);
 
    Kent operator=(const Kent &);

    std::vector<Vector> generate(int);

    std::vector<Vector> generateCanonical(int);

    double eccentricity();

    struct Constants getConstants();

    double computeLogNormalizationConstant();

    double log_dc_dk();

    double log_d2c_dk2();

    double computeSeriesSum(double, double, double);

    double log_dc_db();

    double log_d2c_dkdb();

    double computeSeriesSum2(double, double, double);

    double log_d2c_db2();

    void computeConstants();

    void computeExpectation();

    double computeLogFisherAxes();

    double computeLogFisherScale();

    double log_density(Vector &);

    double computeNegativeLogLikelihood(std::vector<Vector> &);

    double computeNegativeLogLikelihood(Vector &, Matrix &, double);

    double computeNegativeLogLikelihood(struct Estimates &, Vector &, Matrix &, double);

    double computeLogParametersProbability(double);

    double computeLogPriorProbability();

    double computeLogPriorAxes();

    double computeLogPriorScale();

    double computeLogFisherInformation();

    double computeLogFisherInformation(double);

    void computeAllEstimators(
      std::vector<Vector> &, std::vector<struct Estimates> &, int, int
    );

    void computeAllEstimators(
      Vector &, Matrix &, double, std::vector<struct Estimates> &, int, int
    );

    struct Estimates computeAsymptoticMomentEstimates(std::vector<Vector> &);

    struct Estimates computeAsymptoticMomentEstimates(Vector &, Matrix &, double);

    struct Estimates computeMomentEstimates(std::vector<Vector> &);

    struct Estimates computeMomentEstimates(Vector &, Matrix &, double);

    struct Estimates computeMLEstimates(std::vector<Vector> &);

    struct Estimates computeMLEstimates(Vector &, Matrix &, double);

    struct Estimates computeMMLEstimates(std::vector<Vector> &);

    struct Estimates computeMMLEstimates(Vector &, Matrix &, double);

    void estimateParameters(std::vector<Vector> &, Vector &);

    void updateParameters(struct Estimates &);

    Vector Mean();

    Vector MajorAxis();

    Vector MinorAxis();

    double Kappa();

    double Beta();

    double computeKLDivergence(Kent &);

    double computeKLDivergence(struct Estimates &);

    double computeMessageLength(std::vector<Vector> &);

    double computeMessageLength(Vector &, Matrix &, double);

    double computeMessageLength(struct Estimates &, Vector &, Matrix &, double);

    void printParameters(ostream &);

    double computeTestStatistic_vMF(std::vector<Vector> &);

    double computeConfidenceRegion(std::vector<Vector> &);
};

#endif

