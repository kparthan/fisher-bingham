#include "Optimize.h"

Optimize::Optimize(string type)
{
  if (type.compare("MOMENT") == 0) {
    estimation = MOMENT;
  } else if (type.compare("MLE_UNCONSTRAINED") == 0) {
    estimation = MLE_UNCONSTRAINED;
  } else if (type.compare("MLE_CONSTRAINED") == 0) {
    estimation = MLE_CONSTRAINED;
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
  Vector spherical(3,0);
  cartesian2spherical(m0,spherical);
  alpha = spherical[1]; eta = spherical[2];
  cartesian2spherical(m1,spherical);
  psi = spherical[1]; delta = spherical[2];
  c1 = 1000;  
  c2 = c1;
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

    case MLE_UNCONSTRAINED:
    {
      column_vector theta = minimize(sample_mean,S,5);
      finalize(theta,estimates);
      break;
    }

    case MLE_CONSTRAINED:
    {
      column_vector theta = minimize(sample_mean,S,7);
      finalize(theta,estimates);
      break;
    }

    case MML_SCALE:
    {
      column_vector theta = minimize(sample_mean,S,2);
      estimates.kappa = theta(0);
      estimates.beta = theta(1);
      break;
    }

    case MML:
    {
      column_vector theta = minimize(sample_mean,S,5);
      finalize(theta,estimates);
      break;
    }
  }
}

void Optimize::finalize(column_vector &theta, struct Estimates &estimates)
{
  double alpha_f = theta(0);
  double eta_f = theta(1);
  Vector spherical(3,1);
  spherical[1] = alpha_f; spherical[2] = eta_f;
  spherical2cartesian(spherical,estimates.mean);
  double psi_f = theta(2);
  double tmp = -1/(tan(alpha_f) * tan(psi_f));
  double delta_f;
  if (fabs(tmp) > 1) { 
    cout << "yesf, tmp: " << tmp << endl;;
    psi_f = psi;
    delta_f = delta;
  } else {
    double acos_tmp = acos(tmp);
    delta_f = eta_f + acos_tmp;
    if (!(delta_f >= eta_f && delta_f <= PI+eta_f)) {
      cout << "yes2\n";
      delta_f = eta_f - acos_tmp;
      assert(delta_f >= eta_f && delta_f <= PI+eta_f);
    }
  }
  assert(!boost::math::isnan(delta_f));
  spherical[1] = psi_f; spherical[2] = delta_f;
  spherical2cartesian(spherical,estimates.major_axis);
  estimates.minor_axis = crossProduct(estimates.mean,estimates.major_axis);
  estimates.kappa = theta(3);
  estimates.beta = theta(4);
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

    case MLE_UNCONSTRAINED:
    {
      starting_point = alpha,eta,psi,kappa,beta; 
      find_min_using_approximate_derivatives(
        bfgs_search_strategy(),
        objective_delta_stop_strategy(1e-10),
        MaximumLikelihoodObjectiveFunctionUnconstrained(sample_mean,S,N,starting_point),
        starting_point,
        -100
      );
      break;
    }

    case MLE_CONSTRAINED:
    {
      starting_point = alpha,eta,psi,kappa,beta,c1,c2; 
      column_vector min_values(num_params);
      column_vector max_values(num_params);
      min_values = -1e10,-1e10,-1e10,0,0,1e-10,1e-10;
      max_values = uniform_matrix<double>(num_params,1,1e10);
      find_min_box_constrained(
        bfgs_search_strategy(),  
        objective_delta_stop_strategy(1e-9),  
        MaximumLikelihoodObjectiveFunctionConstrained(sample_mean,S,N,starting_point), 
        derivative(MaximumLikelihoodObjectiveFunctionConstrained(sample_mean,S,N,starting_point)), 
        starting_point,min_values,max_values 
      );
      break;
    }

    case MML_SCALE:
    {
      starting_point = kappa,beta; 
      find_min_using_approximate_derivatives(
        bfgs_search_strategy(),
        objective_delta_stop_strategy(1e-10),
        MMLObjectiveFunctionScale(alpha,eta,psi,delta,sample_mean,S,N),
        starting_point,
        -100
      );
      break;
    }

    case MML:
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

