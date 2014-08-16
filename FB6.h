#ifndef FB6_H
#define FB6_H

#include "Header.h"

class FB6
{
  private:
    Vector mu,major_axis,minor_axis;

    long double kappa,beta,gamma;

  public:
    FB6();

    FB6(long double, long double, long double);

    FB6(Vector &, Vector &, Vector &, long double, long double, long double);
 
    FB6 operator=(const FB6 &);

    std::vector<Vector> generate(int);

    std::vector<Vector> generateCanonical(int);

    std::vector<Vector> generate_cartesian_coordinates(Vector &, Vector &);
};

#endif

