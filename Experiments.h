#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "Header.h"

class Experiments
{
  private:
    int iterations;

  public:
    Experiments(int);

    void simulate();

    void computeMeasures(long double, long double, std::vector<Vector> &, std::vector<Vector> &, int);

    Vector computeEstimateMedians(ostream &, std::vector<Vector> &);

    Vector computeEstimateMeans(ostream &, std::vector<Vector> &);

    void computeBias(ostream &, long double, std::vector<Vector> &);

    void computeVariance(ostream &, long double, std::vector<Vector> &);

    void computeMeanAbsoluteError(ostream &, long double, std::vector<Vector> &);

    void computeMeanSquaredError(ostream &, long double, std::vector<Vector> &);

    void computeWinsRatio(bool, const char *, string, int, string);
};

#endif

