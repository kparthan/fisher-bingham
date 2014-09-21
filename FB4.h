#ifndef FB4_H
#define FB4_H

#include "Header.h"

class FB4
{
  private:
    Vector mu,major_axis,minor_axis;

    long double kappa,gamma;

  public:
    FB4();

    FB4(long double, long double);

    FB4(Vector &, Vector &, Vector &, long double, long double);
 
    FB4 operator=(const FB4 &);

    long double computeNormalizationConstant(void);

    std::vector<Vector> generate(int);

    std::vector<Vector> generateCanonical(int);

    Vector generate_u(int);

    int determine_best_case();

    Vector generate_FB4_minus(int, int);

    Vector generate_FB4_plus(int);

    Vector generate_spherical_coordinates(Vector &);

    std::vector<Vector> generate_cartesian_coordinates(Vector &, Vector &);

};

#endif

