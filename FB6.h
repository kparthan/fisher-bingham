#ifndef FB6_H
#define FB6_H

#include "Header.h"

class FB6
{
  private:
    Vector mu,major_axis,minor_axis;

    double kappa,beta,gamma;

  public:
    FB6();

    FB6(double, double, double);

    FB6(Vector &, Vector &, Vector &, double, double, double);
 
    FB6 operator=(const FB6 &);

    std::vector<Vector> generate(int);

    std::vector<Vector> generateCanonical(int);

    std::vector<Vector> generate_cartesian_coordinates(Vector &, Vector &);
};

#endif

