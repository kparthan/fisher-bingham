#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <fstream>
#ifdef HAVE_STD
#include <cstdio>
#else
#include <stdio.h>
#endif

#include "OptLBFGS.h"
#include "NLF.h"

//#include "tstfcn.h"
using namespace std;
using NEWMAT::ColumnVector;

using namespace OPTPP;
void update_model(int, int, ColumnVector) {}

void init_rosen(int ndim, ColumnVector& x)
{
  if (ndim != 2)
    exit (1);
  x(1) = 1.2;
  x(2) =  1.0;
}

void rosen(int ndim, const ColumnVector& x, double& fx, int& result)
{
  double f1, f2, x1, x2;

  if (ndim != 2)
    exit (1);

  x1 = x(1);
  x2 = x(2);
  f1 = (x2 - x1 * x1);
  f2 = 1. - x1;
  
  fx  = 100.* f1*f1 + f2*f2;
  result = NLPFunction;
}

int main ()
{
  int n = 2;
  
  static char *status_file = {"tstLBFGS.out"};

  //  Create a Nonlinear problem object

  FDNLF1 nlp(n,rosen,init_rosen);
  
  //  Build a LBFGS object and optimize 

  OptLBFGS objfcn(&nlp);   
  objfcn.setUpdateModel(update_model);
  if (!objfcn.setOutputFile(status_file, 0))
    cerr << "main: output file open failed" << endl;
  objfcn.setGradTol(1.e-8);
  objfcn.setMaxBacktrackIter(10);
  objfcn.setPrintFinalX(true);
  objfcn.optimize();
    
  objfcn.printStatus("Solution from LBFGS: More and Thuente's linesearch");

//#ifdef REG_TEST
  ColumnVector x_sol = nlp.getXc();
  double f_sol = nlp.getF();

  cout << "sol: (" << x_sol(1) << ", " << x_sol(2) << ")" << endl;
  cout << "f_sol: " << f_sol << endl;

  ostream* optout = objfcn.getOutputFile();
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
                                                                 <=
                                                                 1.e-2))
    *optout << "LBFGS 1 PASSED" << endl;
  else
    *optout << "LBFGS 1 FAILED" << endl;
//#endif

  objfcn.cleanup();

}

