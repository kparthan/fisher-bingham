#ifndef FB6_H
#define FB6_H

#include "Header.h"

class FB6
{
  private:
    long double kappa,beta,gamma;

  public:
    FB6();

    FB6(long double, long double, long double);

    FB6 operator=(const FB6 &);

    std::vector<std::vector<long double> > generate(int);
};

#endif

