#ifndef KENT_ECC_TRANS_H
#define KENT_ECC_TRANS_H

#include "Header.h"

class Kent_EccTrans  
{
  private:
    Vector mu,major_axis,minor_axis;

    double psi,alpha,eta;
    
    double kappa,ecc; // gamma = 0

    struct Constants {
      double log_c;
    } constants;

    int computed;

  public:
    Kent_EccTrans();

    Kent_EccTrans(double, double);

    Kent_EccTrans(Vector &, Vector &, Vector &, double, double);

    Kent_EccTrans(double, double, double, double, double);
 
    Kent_EccTrans operator=(const Kent_EccTrans &);

    double computeLogNormalizationConstant();

    double computeSeriesSum(double, double, double);

    double log_density(Vector &);

    double computeNegativeLogLikelihood(std::vector<Vector> &);

    double computeNegativeLogLikelihood(Vector &, Matrix &, double);

    double computeLogPriorProbability();

    double computeLogPriorAxes();

    double computeLogPriorScale();

    Vector Mean();

    Vector MajorAxis();

    Vector MinorAxis();

    double Kappa();

    double Beta();

    double eccentricity();

    void printParameters(ostream &);
};

#endif

