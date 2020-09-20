#include <fstream>

#include "OptNewton.h"
#include "NLF.h"
//#include "tstfcn.h"

using NEWMAT::ColumnVector;
using namespace OPTPP;

void update_model(int, int, ColumnVector) {}

void init_rosen(int ndim, ColumnVector& x)
{
  if (ndim != 2)
    exit (1);
  x(1) = -1.2;
  x(2) =  1.0;
}

void rosen2(int ndim, const ColumnVector& x, double& fx, int& result)
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
  
  static char *status_file = {"tstnewton.out"};

//----------------------------------------------------------------------------
// 1. Newton with trust regions
//----------------------------------------------------------------------------

  //  Create a Nonlinear problem object

  FDNLF1 nlp(n,rosen2,init_rosen);
  
  //  Build a Newton object and optimize 

  OptNewton objfcn(&nlp, update_model);
  if (!objfcn.setOutputFile(status_file, 0))
    cerr << "main: output file open failed" << endl;
  objfcn.setTRSize(1.0e2);
  objfcn.optimize();
  objfcn.printStatus("Solution from newton: trust regions");

#ifdef REG_TEST
  ColumnVector x_sol = nlp.getXc();
  double f_sol = nlp.getF();
  ostream* optout = objfcn.getOutputFile();
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
                                                                 <=
                                                                 1.e-2))
    *optout << "Newton 1 PASSED" << endl;
  else
    *optout << "Newton 1 FAILED" << endl;
#endif

  objfcn.cleanup();      
    
//----------------------------------------------------------------------------
// 2. Newton with More and Thuente's line search
//----------------------------------------------------------------------------
//
//  FDNLF1 nlp2(n,rosen2,init_rosen);
//  
//  OptNewton objfcn2(&nlp2,update_model);   
//  objfcn2.setSearchStrategy(LineSearch);
//  objfcn2.setOutputFile(status_file, 1);
//  objfcn2.optimize();
//  objfcn2.printStatus("Solution from newton: More and Thuente linesearch");
//
//#ifdef REG_TEST
//  x_sol = nlp2.getXc();
//  f_sol = nlp2.getF();
//  optout = objfcn2.getOutputFile();
//  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
//                                                                 <=
//                                                                 1.e-2))
//    *optout << "Newton 2 PASSED" << endl;
//  else
//    *optout << "Newton 2 FAILED" << endl;
//#endif
//
//  objfcn2.cleanup();     
//    
////----------------------------------------------------------------------------
//// 3. Newton with backtracking line search
////----------------------------------------------------------------------------
//
//  FDNLF1 nlp3(n,rosen2,init_rosen);
//  nlp3.setIsExpensive(true);
//  
//  OptNewton objfcn3(&nlp3,update_model);   
//  objfcn3.setOutputFile(status_file, 1);
//  objfcn3.setSearchStrategy(LineSearch);
//  objfcn3.optimize();
//  objfcn3.printStatus("Solution from newton: backtracking linesearch");
//
//#ifdef REG_TEST
//  x_sol = nlp3.getXc();
//  f_sol = nlp3.getF();
//  optout = objfcn3.getOutputFile();
//  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
//                                                                 <=
//                                                                 1.e-2))
//    *optout << "Newton 3 PASSED" << endl;
//  else
//    *optout << "Newton 3 FAILED" << endl;
//#endif
//
//  objfcn3.cleanup();     
}
