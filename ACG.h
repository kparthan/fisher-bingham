#ifndef ACG_H
#define ACG_H

#include "Header.h"

/*!
 *  ACG: Angular Central Gaussian distribution
 *  pdf(x) = |W|^(1/2) (x' W x)^(-D/2)
 *  x is a D X 1 matrix (vector)
 *  W is a +ve definite symmetric matrix
 */
class ACG
{
  private:
    int D;

    Matrix W;

  public:
    ACG();

    ACG(Matrix &);

    ACG operator=(const ACG &);

    void printParameters();

    std::vector<Vector> generate(int);
};

#endif

