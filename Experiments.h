#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "Header.h"

class Experiments
{
  private:
    int iterations;

  public:
    Experiments(int);

    void plotBias(long double, long double);
};

#endif

