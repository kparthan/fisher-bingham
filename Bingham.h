#ifndef BINGHAM_H
#define BINGHAM_H

#include "Header.h"

/*!
 *  pdf(x) = constant * exp(- x' A x)
 */
class Bingham
{
  private:
    int D;

    Matrix A;

  public:
    Bingham();

    Bingham(Matrix &);

    Bingham operator=(const Bingham &);

    void printParameters();

    std::vector<Vector> generate(int);
};

#endif

