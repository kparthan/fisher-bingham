#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "Header.h"

class Experiments
{
  private:
    int iterations;

  public:
    Experiments(int);

    void simulate(double, double);

    void computeMeasures(double, double, std::vector<Vector> &, std::vector<Vector> &, int);

    Vector computeEstimateMedians(ostream &, std::vector<Vector> &);

    Vector computeEstimateMeans(ostream &, std::vector<Vector> &);

    void computeBias(ostream &, double, std::vector<Vector> &);

    void computeVariance(ostream &, double, std::vector<Vector> &);

    void computeMeanAbsoluteError(ostream &, double, std::vector<Vector> &);

    void computeMeanSquaredError(ostream &, double, std::vector<Vector> &);

    void computeWinsRatio(bool, const char *, string, int, string);
};

#endif

