#ifndef KENT_UNIF_TRANS_H
#define KENT_UNIF_TRANS_H

#include "Header.h"

class Kent_UnifTrans  
{
  private:
    Vector mu,major_axis,minor_axis;

    double z1,z2,z3,z4,z5;
    
    struct Constants {
      double log_c;
    } constants;

    int computed;

  public:
    Kent_UnifTrans(double, double, double, double, double);
 
    Kent_UnifTrans operator=(const Kent_UnifTrans &);

    double computeLogNormalizationConstant();

    double computeSeriesSum(double, double, double);

    double log_density(Vector &);

    double computeNegativeLogLikelihood(std::vector<Vector> &);

    double computeNegativeLogLikelihood(Vector &, Matrix &, double);

    Vector Mean();

    Vector MajorAxis();

    Vector MinorAxis();

    double Kappa();

    double Beta();

    double eccentricity();
};

#endif

