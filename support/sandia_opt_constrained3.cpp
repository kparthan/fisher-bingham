#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#if defined(SGI) || defined(RS6K)
#define WANT_MATH
#else
#define WANT_STREAM
#define WANT_MATH
#endif

#ifdef HAVE_STD
#include <cstdio>
#include <cstdlib>
#include <cmath>
#else
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#endif

#include "newmatap.h"

#include "NLP.h"
#include "NLF.h"
#include "BoundConstraint.h"
#include "NonLinearInequality.h"
#include "CompoundConstraint.h"
#include "OptNIPS.h"
#include "OptQNIPS.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

using namespace OPTPP;

void init_hs65(int ndim, ColumnVector& x)
{
  if (ndim != 3)
  {
    exit (1);
  }
  double factor = 0.0;
  x(1) = -5.0  - (factor - 1)*8.6505  ;
  x(2) =  5.0  + (factor - 1)*1.3495  ;
  x(3) =  0.0  - (factor - 1)*4.6204  ;
}

void hs65(int n, const ColumnVector& x, double& fx, int& result)
{ // Hock and Schittkowski's Problem 65 (the objective fcn)
  double f1, f2, f3, x1, x2, x3;

  if (n != 3) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  f1 = x1 - x2;
  f2 = x1 + x2 - 10.0;
  f3 = x3 - 5.0;

  //if (mode & NLPFunction) {
    fx  = f1*f1+ (f2*f2)*(1.0/9.0) +f3*f3;
    result = NLPFunction;
  //}
  /*if (mode & NLPGradient) {
    g(1) =  2*f1 + (2.0/9.0)*f2;
    g(2) = -2*f1 + (2.0/9.0)*f2;
    g(3) =  2*f3;
    result = NLPGradient;
  }
  if (mode & NLPHessian) {
    Hx(1,1) =  2 + (2.0/9.0);

    Hx(2,1) = -2 + (2.0/9.0);
    Hx(2,2) =  2 + (2.0/9.0);

    Hx(3,1) = 0.0;
    Hx(3,2) =  0.0;
    Hx(3,3) =  2.0;
    result = NLPHessian;
  }*/
}

void ineq_hs65(int mode, int ndim, const ColumnVector& x, ColumnVector& cx, Matrix& cgx, OptppArray<SymmetricMatrix>& cHx, int& result)
{ // Hock and Schittkowski's Problem 65 
  double f1, f2, f3, x1, x2, x3;
  SymmetricMatrix Htmp(ndim);

  if (ndim != 3)
     exit(1);

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  f1 = x1;
  f2 = x2;
  f3 = x3;

  if (mode & NLPFunction) {
    cx(1)  = 48 - f1*f1 - f2*f2 - f3*f3;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    cgx(1,1) = -2*x1;
    cgx(2,1) = -2*x2;
    cgx(3,1) = -2*x3;
    result = NLPGradient;
  }
  if (mode & NLPHessian) {
    Htmp(1,1) = -2;
    Htmp(1,2) = 0.0;
    Htmp(1,3) = 0.0;
    Htmp(2,1) = 0.0;
    Htmp(2,2) = -2;
    Htmp(2,3) = 0.0;
    Htmp(3,1) = 0.0;
    Htmp(3,2) = 0.0;
    Htmp(3,3) = -2;

    cHx[0] = Htmp;
    result = NLPHessian;
  }
}

int main ()
{
  int ndim = 3;
  ColumnVector lower(ndim), upper(ndim); 

  // Here is one way to assign values to a ColumnVector.

  lower << -4.5 << -4.5 << -5.0;
  upper <<  4.5 <<  4.5 <<  5.0 ;

  // simple boundary constraints ...
  Constraint c1 = new BoundConstraint(ndim, lower, upper);

  // non-linear constraint
  NLP* chs65 = new NLP(new NLF2(ndim, 1, ineq_hs65, init_hs65));
  Constraint nleqn = new NonLinearInequality(chs65);

  CompoundConstraint* constraints = new CompoundConstraint(nleqn, c1);

  FDNLF1 nips(ndim, hs65, init_hs65, constraints);
  OptQNIPS objfcn(&nips);
  objfcn.setOutputFile("example2.out", 0);
  objfcn.setFcnTol(1.0e-06);
  objfcn.setMaxIter(150);
  objfcn.setMeritFcn(ArgaezTapia);
  objfcn.optimize();

  objfcn.printStatus("Solution from nips");
  objfcn.cleanup();

}

