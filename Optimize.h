#ifndef OPTIMIZE
#define OPTIMIZE

#include <dlib/optimization.h>
#include "Header.h"

using namespace dlib;

typedef dlib::matrix<double,0,1> column_vector;

class Optimize
{
  private:
    Vector mu_est,major_est,minor_est; 

    long double kappa,beta;

    //long double C1,C2;

  public:
    Optimize();

    Optimize(Vector &, Vector &, Vector &);

    void computeMomentEstimates(Vector &, Matrix &);

    //double momentObjectiveFunction(const column_vector &);

    void minimize(void);
};

#endif

