#ifndef KENT_H
#define KENT_H

#include "Header.h"

class Kent  // FB5
{
  private:
    std::vector<long double> mu,major_axis,minor_axis;

    long double kappa,beta; // gamma = 0

  public:
    Kent();

    Kent(long double, long double);

    Kent(std::vector<long double> &, std::vector<long double> &, std::vector<long double> &, 
        long double, long double);
 
    Kent operator=(const Kent &);

    std::vector<std::vector<long double> > generate(int);

    std::vector<std::vector<long double> > generateCanonical(int);

    long double eccentricity();

    long double computeLogNormalizationConstant();

    long double computeLogNormalizationConstant(long double, long double);
};

#endif

