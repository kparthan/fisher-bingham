#include "Optimize.h"

Optimize::Optimize(string type)
{
  if (type.compare("MOMENT") == 0) {
    estimation = MOMENT;
  } else if (type.compare("MLE") == 0) {
    estimation = MLE;
  } else if (type.compare("MML_SCALE") == 0) {
    estimation = MML_SCALE;
  } else if (type.compare("MML") == 0) {
    estimation = MML;
  }
}

void Optimize::initialize(int sample_size, Vector &m0, Vector &m1, Vector &m2, long double k, long double b)
{
  N = sample_size;
  mean = m0;
  major = m1;
  minor = m2;
  kappa = k;
  beta = b;
}

void Optimize::computeEstimates(Vector &sample_mean, Matrix &S, struct Estimates &estimates)
{
  switch(estimation) {
    case MOMENT:
    {
      column_vector theta = minimize(sample_mean,S,2);
      estimates.kappa = theta(0);
      estimates.beta = theta(1);
      break;
    }

    case MLE:
    {
      computeOrthogonalTransformation(mean,major,psi,alpha,eta);
      column_vector theta = minimize(sample_mean,S,5);
      finalize(theta,estimates);
      break;
    }

    case MML_SCALE:
    {
      computeOrthogonalTransformation(mean,major,psi,alpha,eta);
      column_vector theta = minimize(sample_mean,S,2);
      estimates.kappa = theta(0);
      estimates.beta = theta(1);
      break;
    }

    /*case MML:
    {
      column_vector theta = minimize(sample_mean,S,5);
      finalize(theta,estimates);
      break;
    }*/
  }
}

void Optimize::finalize(column_vector &theta, struct Estimates &estimates)
{
  double alpha_f = theta(0);
  double eta_f = theta(1);
  double psi_f = theta(2);
  estimates.kappa = theta(3);
  estimates.beta = theta(4);
  Kent kent(psi_f,alpha_f,eta_f,estimates.kappa,estimates.beta);
  estimates.mean = kent.Mean();
  estimates.major_axis = kent.MajorAxis();
  estimates.minor_axis = kent.MinorAxis();
}

column_vector Optimize::minimize(Vector &sample_mean, Matrix &S, int num_params)
{
  column_vector starting_point(num_params);
  switch(estimation) {
    case MOMENT:
    {
      starting_point = kappa,beta; 
      find_min_using_approximate_derivatives(
        bfgs_search_strategy(),
        objective_delta_stop_strategy(1e-10),
        MomentObjectiveFunction(mean,major,minor,sample_mean,S,N),
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
        MaximumLikelihoodObjectiveFunction(sample_mean,S,N,starting_point),
        starting_point,
        -100
      );
      break;
    }

    case MML_SCALE:
    {
      starting_point = kappa,beta; 
      find_min_using_approximate_derivatives(
        bfgs_search_strategy(),
        objective_delta_stop_strategy(1e-10),
        MMLObjectiveFunctionScale(psi,alpha,eta,kappa,beta,sample_mean,S,N),
        starting_point,
        -100
      );
      break;
    }

    /*case MML:
    {
      starting_point = alpha,eta,psi,kappa,beta; 
      find_min_using_approximate_derivatives(
        bfgs_search_strategy(),
        objective_delta_stop_strategy(1e-10),
        MMLObjectiveFunction(sample_mean,S,N,starting_point),
        starting_point,
        -100
      );
      break;
    }*/
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

