#ifndef FB4_H
#define FB4_H

#include "Header.h"

class FB4
{
  private:
    std::vector<long double> mu,major_axis,minor_axis;

    long double kappa,gamma;

  public:
    FB4();

    FB4(long double, long double);

    FB4(std::vector<long double> &, std::vector<long double> &, std::vector<long double> &, 
        long double, long double);
 
    FB4 operator=(const FB4 &);

    long double computeNormalizationConstant(void);

    std::vector<std::vector<long double> > generate(int);

    std::vector<std::vector<long double> > generateCanonical(int);

    std::vector<long double> generate_u(int);

    int determine_best_case();

    std::vector<long double> generate_FB4_minus(int, int);

    std::vector<long double> generate_FB4_plus(int);

    std::vector<long double> generate_spherical_coordinates(std::vector<long double> &);

    std::vector<std::vector<long double> > generate_cartesian_coordinates(
      std::vector<long double> &, std::vector<long double> &
    );

};

#endif

