#ifndef FB4_H
#define FB4_H

#include "Header.h"

class FB4
{
  private:
    long double kappa,gamma;

  protected:
    int determine_best_case();

    std::vector<std::vector<long double> > generate_FB4_minus(int, int);

  public:
    FB4();

    FB4(long double, long double);

    FB4 operator=(const FB4 &);

    std::vector<std::vector<long double> > generate(int);
};

#endif

