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
      Vector E_x,kappa_E_x;
      Matrix E_xx,beta_E_xx;
      Matrix R,Rt;  // R: standard -> current orientation
    } constants;

    struct TrignometryConstants {
      long double cos_alpha,sin_alpha,tan_alpha;
      long double cos_eta,sin_eta;
      long double cos_psi,sin_psi,tan_psi;
      long double cos_delta,sin_delta;
      long double d_n,cos_d_n,sin_d_n;
    } tc;

    //  d1_mu (3 X 3): <dmu_da> <dmu_dn> <dmu_ds>
    //  d2_mu (6 X 3): [0] <d2mu_da2> 
    //                 [1] <d2mu_dn2> 
    //                 [2] <d2mu_ds2> 
    //                 [3] <d2mu_dadn> 
    //                 [4] <d2mu_dads> 
    //                 [5] <d2mu_dnds>
    struct Differentials {
      std::vector<Vector > d1_mu,d1_mj,d1_mi;
      std::vector<Vector > d2_mu,d2_mj,d2_mi;
      long double ddel_da,ddel_ds,d2del_da2,d2del_ds2,d2del_dads;
      Matrix fisher_axes;
    } df;

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

    long double computeExpectationLikelihood(int, int, int);

    void computeFirstOrderDifferentials();

    Vector computeFirstOrderDifferentialsMinorAxis(int);

    void computeDeltaDifferentials();

    void computeSecondOrderDifferentials();

    Vector computeSecondOrderDifferentialsMinorAxis(int, int, int, int, int, int);

    void computeFisherMatrixAxes();

    long double computeLogFisherScale();

    long double computeNegativeLogLikelihood(std::vector<Vector> &);

    long double computeNegativeLogLikelihood(Vector &, Matrix &, int);

    long double computeLogPriorProbability();

    long double computeLogPriorAxes();

    long double computeLogPriorScale();

    long double computeLogFisherInformation();

    long double computeLogFisherInformation(int);

    void computeAllEstimators(std::vector<Vector> &);

    void computeAllEstimators(Vector &, Matrix &, int);

    struct Estimates computeMomentEstimates(std::vector<Vector> &);

    struct Estimates computeMomentEstimates(Vector &, Matrix &, int);

    struct Estimates computeMLEstimates(std::vector<Vector> &, string);

    struct Estimates computeMLEstimates(Vector &, Matrix &, int, string);

    struct Estimates computeMMLEstimates(std::vector<Vector> &);

    struct Estimates computeMMLEstimates(Vector &, Matrix &, int);

    Vector Mean();

    Vector MajorAxis();

    Vector MinorAxis();

    long double Kappa();

    long double Beta();

    long double computeKLDivergence(Kent &);
};

#endif

