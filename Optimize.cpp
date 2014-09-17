#include "Optimize.h"

Optimize::Optimize(string type)
{
  if (type.compare("MOMENT") == 0) {
    estimation = MOMENT;
  } else if (type.compare("MLE") == 0) {
    estimation = MLE;
  } else if (type.compare("MAP") == 0) {
    estimation = MAP;
  } else if (type.compare("MML_SCALE") == 0) {
    estimation = MML_SCALE;
  } else if (type.compare("MML") == 0) {
    estimation = MML;
  }
}

void Optimize::initialize(double sample_size, Vector &m0, Vector &m1, Vector &m2, long double k, long double b)
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
      std::vector<double> theta = minimize(sample_mean,S,2);
      estimates.kappa = theta(0);
      estimates.beta = theta(1);
      break;
    }

    case MLE:
    {
      computeOrthogonalTransformation(mean,major,psi,alpha,eta);
      std::vector<double> theta = minimize(sample_mean,S,5);
      finalize(theta,estimates);
      break;
    }

    case MAP:
    {
      computeOrthogonalTransformation(mean,major,psi,alpha,eta);
      std::vector<double> theta = minimize(sample_mean,S,5);
      finalize(theta,estimates);
      break;
    }

    case MML_SCALE:
    {
      computeOrthogonalTransformation(mean,major,psi,alpha,eta);
      std::vector<double> theta = minimize(sample_mean,S,2);
      estimates.kappa = theta(0);
      estimates.beta = theta(1);
      break;
    }

    case MML:
    {
      computeOrthogonalTransformation(mean,major,psi,alpha,eta);
      std::vector<double> theta = minimize(sample_mean,S,5);
      finalize(theta,estimates);
      break;
    }
  }
}

void Optimize::finalize(std::vector<double> &theta, struct Estimates &estimates)
{
  estimates.alpha = theta(0);
  estimates.eta = theta(1);
  estimates.psi = theta(2);
  estimates.kappa = theta(3);
  estimates.beta = theta(4);
  Kent kent(estimates.psi,estimates.alpha,estimates.eta,estimates.kappa,estimates.beta);
  estimates.mean = kent.Mean();
  estimates.major_axis = kent.MajorAxis();
  estimates.minor_axis = kent.MinorAxis();
}

std::vector<double> Optimize::minimize(Vector &sample_mean, Matrix &S, int num_params)
{
  std::vector<double> starting_point(num_params);
  switch(estimation) {
    case MOMENT:
    {
      nlopt::opt opt(nlopt::LN_COBYLA, 2);
      std::vector<double> lb(2);
      lb[0] = TOLERANCE; lb[1] = TOLERANCE;
      opt.set_lower_bounds(lb);
      std::vector<double> ub(2);
      ub[0] = HUGE_VAL; ub[1] = HUGE_VAL;
      opt.set_upper_bounds(ub);

      //opt.set_min_objective(myvfunc, NULL);
      MomentObjectiveFunction moment(mean,major,minor,sample_mean,S,N);
      opt.set_min_objective(MomentObjectiveFunction::wrap, &moment);

      my_constraint_data data[2] = { {2,0}, {-1,1} };
      opt.add_inequality_constraint(myvconstraint, &data[0], 1e-8);
      opt.add_inequality_constraint(myvconstraint, &data[1], 1e-8);

      opt.set_xtol_rel(1e-4);

      std::vector<double> x(2);
      x[0] = 1.234; x[1] = 5.678;
      double minf;
      nlopt::result result = opt.optimize(x, minf);
      break;
    }

    /*case MLE:
    {
      starting_point = alpha,eta,psi,kappa,beta; 
      find_min_using_approximate_derivatives(
        bfgs_search_strategy(),
        objective_delta_stop_strategy(1e-10),
        MaximumLikelihoodObjectiveFunction(sample_mean,S,N),
        starting_point,
        -100
      );
      break;
    }

    case MAP:
    {
      starting_point = alpha,eta,psi,kappa,beta; 
      find_min_using_approximate_derivatives(
        bfgs_search_strategy(),
        objective_delta_stop_strategy(1e-10),
        MAPObjectiveFunction(sample_mean,S,N),
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

    case MML:
    {
      starting_point = alpha,eta,psi,kappa,beta; 
      find_min_using_approximate_derivatives(
        bfgs_search_strategy(),
        objective_delta_stop_strategy(1e-10),
        MMLObjectiveFunction(sample_mean,S,N),
        starting_point,
        -100
      );
      break;
    }*/

    default:
      break;
  }
  //cout << "solution:\n" << starting_point << endl;
  return starting_point;
}

