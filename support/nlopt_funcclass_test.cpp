#include <iostream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <nlopt.hpp>

using namespace std;

typedef struct {
    double a, b;
} my_constraint_data;

class Function
{
  private:
    //my_constraint_data *d;
    //void *d;

  public:
    //Function() {}

    //Function(void *data) {
      //d = reinterpret_cast<my_constraint_data*>(data);
      //d = data;
    //}

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        return (*reinterpret_cast<Function*>(data))(x, grad); 
    }

    double operator()(const std::vector<double> &x, std::vector<double> &grad)
    {
      return sqrt(x[1]);
    }
};

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    /*if (!grad.empty()) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }*/
    return sqrt(x[1]);
};

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

  //opt.set_min_objective(myvfunc, NULL);
  Function function;
  opt.set_min_objective(Function::wrap, &function);

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

