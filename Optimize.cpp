#include "Optimize.h"
#include "Support.h"
#include "Kent.h"

long double C1,C2;

Optimize::Optimize()
{}

Optimize::Optimize(Vector &mu_est, Vector &major_est, Vector &minor_est) : 
          mu_est(mu_est), major_est(major_est), minor_est(minor_est)
{}

void Optimize::computeMomentEstimates(Vector &sample_mean, Matrix &S)
{
  C1 = computeDotProduct(sample_mean,mu_est);
  Vector x = prod(S,major_est);
  long double mj = computeDotProduct(major_est,x);
  x = prod(S,minor_est);
  long double mi = computeDotProduct(minor_est,x);
  C2 = mj - mi;
  // minimize function: log c(k,b) - k * a - b * delta
  minimize();
}

/*!
 *  minimize function: log c(k,b) - k * a - b * delta
 */
double momentObjectiveFunction(const column_vector &x)
{
  const double k = x(0);
  const double b = x(1);

  Kent kent(k,b);
  long double log_norm = kent.computeLogNormalizationConstant();
  double fval = log_norm - k * C1 - b * C2;
  cout << setprecision(6);
  //cout << "log_norm: " << log_norm << "; C1: " << C1 << "; C2: " << C2 << endl;
  cout << "x: (" << x(0) << "," << x(1) << "); fx: " << fval << endl;
  return fval;
}

void Optimize::minimize(void)
{
  column_vector starting_point(2);
  starting_point = 40,9;
  find_min_using_approximate_derivatives(bfgs_search_strategy(),
                                         objective_delta_stop_strategy(1e-10),
                                         momentObjectiveFunction,starting_point,-100);
  cout << "solution:\n" << starting_point << endl;
}

