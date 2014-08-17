#include "Optimize.h"

Optimize::Optimize()
{}

void Optimize::initialize(Vector &m0, Vector &m1, Vector &m2, long double k, long double b)
{
  mean = m0;
  major = m1;
  minor = m2;
  kappa = k;
  beta = b;
}

void Optimize::computeMomentEstimates(Vector &sample_mean, Matrix &S, struct Estimates &estimates)
{
  column_vector theta = minimize(sample_mean,S,MOMENT,2);
  estimates.kappa = theta(0);
  estimates.beta = theta(1);
}

void Optimize::computeMLEstimates(Vector &sample_mean, Matrix &S, struct Estimates &estimates)
{
  // initialize
  computeMomentEstimates(sample_mean,S,estimates);
  Vector spherical(3,0); 
  cartesian2spherical(estimates.mean,spherical);
  alpha = spherical[1];
  eta = spherical[2];
  cartesian2spherical(estimates.major_axis,spherical);
  psi = spherical[1];
  delta = spherical[2];
  kappa = estimates.kappa;
  beta = estimates.beta;

  // minimize
  column_vector theta = minimize(sample_mean,S,MLE,5);
  alpha = theta(0);
  eta = theta(1);
  spherical[1] = alpha; spherical[2] = eta;
  spherical2cartesian(spherical,estimates.mean);
  psi = theta(2);
  double tmp = -1/(tan(alpha) * tan(psi));
  double delta = eta + acos(tmp);
  spherical[1] = psi; spherical[2] = delta;
  spherical2cartesian(spherical,estimates.major_axis);
  estimates.minor_axis = crossProduct(estimates.major_axis,estimates.minor_axis);
  estimates.kappa = theta(3);
  estimates.beta = theta(4);

}

column_vector Optimize::minimize(Vector &sample_mean, Matrix &S, int estimation, int num_params)
{
  column_vector starting_point(num_params);
  switch(estimation) {
    case MOMENT:
    {
      starting_point = kappa,beta; 
      find_min_using_approximate_derivatives(
        bfgs_search_strategy(),
        objective_delta_stop_strategy(1e-10),
        MomentObjectiveFunction(mean,major,minor,sample_mean,S),
        starting_point,
        -100
      );
      break;
    }

    case MLE:
    {
      starting_point = alpha,eta,psi,kappa,beta; 
      find_min_using_approximate_derivatives(
        bfgs_search_strategy(),
        objective_delta_stop_strategy(1e-10),
        MaximumLikelihoodObjectiveFunction(sample_mean,S),
        starting_point,
        -100
      );
      break;
    }
  }
  /*find_min_box_constrained(bfgs_search_strategy(),  
                           objective_delta_stop_strategy(1e-9),  
                           momentObjectiveFunction, derivative(momentObjectiveFunction), 
                           starting_point, 
                           uniform_matrix<double>(2,1,0.0),
                           uniform_matrix<double>(2,1,2000.0));*/
  /*find_min_bobyqa(momentObjectiveFunction,
                  starting_point, 
                  5,    // number of interpolation points
                  uniform_matrix<double>(2,1,-1e100),  // lower bound constraint
                  uniform_matrix<double>(2,1, 2000),   // upper bound constraint
                  10,    // initial trust region radius
                  1e-6,  // stopping trust region radius
                  100    // max number of objective function evaluations
  );*/
  cout << "solution:\n" << starting_point << endl;
  return starting_point;
}

