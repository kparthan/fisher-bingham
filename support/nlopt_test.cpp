#include <iostream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <nlopt.hpp>

using namespace std;

typedef struct {
    double a, b;
} my_constraint_data;

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    /*if (!grad.empty()) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }*/
    return sqrt(x[1]);
}

double myvconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
    //my_constraint_data *d = (my_constraint_data*) data;
    double a = d->a, b = d->b;
    /*if (!grad.empty()) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }*/
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}

int main(int argc, char **argv)
{
  cout << "Testing NLopt ...\n";

  // C++ style ...
  //nlopt::opt opt(nlopt::LD_MMA, 2);
  nlopt::opt opt(nlopt::LN_COBYLA, 2);
  //nlopt::opt opt(nlopt::GN_ISRES, 2);
  //opt::opt opt(nlopt::GN_ORIG_DIRECT, 2);

  std::vector<double> lb(2);
  lb[0] = -HUGE_VAL; lb[1] = 1e-6;
  opt.set_lower_bounds(lb);
  std::vector<double> ub(2);
  ub[0] = HUGE_VAL; ub[1] = HUGE_VAL;
  opt.set_upper_bounds(ub);

  opt.set_min_objective(myvfunc, NULL);

  my_constraint_data data[2] = { {2,0}, {-1,1} };
  opt.add_inequality_constraint(myvconstraint, &data[0], 1e-8);
  opt.add_inequality_constraint(myvconstraint, &data[1], 1e-8);

  opt.set_xtol_rel(1e-4);

  std::vector<double> x(2);
  x[0] = 1.234; x[1] = 5.678;
  double minf;
  nlopt::result result = opt.optimize(x, minf);

  cout << "sol: (" << x[0] << ", " << x[1] << ")\n";
  cout << "minf: " << minf << endl;

  return 0;
}

/*double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}

double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
    my_constraint_data *d = (my_constraint_data *) data;
    double a = d->a, b = d->b;
    if (grad) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
 }*/

  // C style ...
//  double lb[2] = { -HUGE_VAL, 0 }; /* lower bounds */
//  nlopt_opt opt;
//
//  opt = nlopt_create(NLOPT_LD_MMA, 2); /* algorithm and dimensionality */
//  nlopt_set_lower_bounds(opt, lb);
//  nlopt_set_min_objective(opt, myfunc, NULL);
//
//  my_constraint_data data[2] = { {2,0}, {-1,1} };
//  nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-8);
//  nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-8);
//
//  nlopt_set_xtol_rel(opt, 1e-4);
//
//  double x[2] = { 1.234, 5.678 };  /* some initial guess */
//  double minf; /* the minimum objective value, upon return */
//
//  if (nlopt_optimize(opt, x, &minf) < 0) {
//      printf("nlopt failed!\n");
//  }
//  else {
//      printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
//  }
//
//  nlopt_destroy(opt);
