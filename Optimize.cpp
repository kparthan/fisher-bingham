#include "Optimize.h"

Optimize::Optimize(string type)
{
  if (type.compare("MOMENT") == 0) {
    estimation = MOMENT;
  } else if (type.compare("MLE") == 0) {
    estimation = MLE;
  } else if (type.compare("MAP") == 0) {
    estimation = MAP;
  } else if (type.compare("MML_2") == 0) {
    estimation = MML_2;
  } else if (type.compare("MML_5") == 0) {
    estimation = MML_5;
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
      estimates.kappa = theta[0];
      estimates.beta = theta[1];
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

    case MML_2:
    {
      computeOrthogonalTransformation(mean,major,psi,alpha,eta);
      std::vector<double> theta = minimize(sample_mean,S,2);
      estimates.kappa = theta[0];
      estimates.beta = theta[1];
      break;
    }

    case MML_5:
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
  estimates.alpha = theta[0];
  estimates.eta = theta[1];
  estimates.psi = theta[2];
  estimates.kappa = theta[3];
  estimates.beta = theta[4];
  Kent kent(estimates.psi,estimates.alpha,estimates.eta,estimates.kappa,estimates.beta);
  estimates.mean = kent.Mean();
  estimates.major_axis = kent.MajorAxis();
  estimates.minor_axis = kent.MinorAxis();
}

std::vector<double> Optimize::minimize(Vector &sample_mean, Matrix &S, int num_params)
{
  std::vector<double> x(num_params);
  //nlopt::opt opt(nlopt::LN_NELDERMEAD, num_params);
  nlopt::opt opt(nlopt::LN_COBYLA, num_params);
  std::vector<double> lb(num_params,TOLERANCE);
  std::vector<double> ub(num_params,HUGE_VAL);

  // GN_ISRES and GN_ORIG_DIRECT work with finite bounds
  //nlopt::opt opt(nlopt::GN_ISRES, num_params);
  //nlopt::opt opt(nlopt::GN_ORIG_DIRECT, num_params);
  // std::vector<double> ub(num_params,1000);
  double minf;

  switch(estimation) {
    case MOMENT:
    {
      opt.set_lower_bounds(lb);
      opt.set_upper_bounds(ub);

      MomentObjectiveFunction moment(mean,major,minor,sample_mean,S,N);
      opt.set_min_objective(MomentObjectiveFunction::wrap, &moment);
      opt.add_inequality_constraint(Constraint2, NULL, TOLERANCE);
      opt.set_xtol_rel(TOLERANCE);

      x[0] = kappa; x[1] = beta;
      nlopt::result result = opt.optimize(x, minf);
      break;
    }

    case MLE:
    {
      opt.set_lower_bounds(lb);
      ub[0] = PI; ub[1] = 2*PI; ub[2] = 2*PI;
      opt.set_upper_bounds(ub);

      MaximumLikelihoodObjectiveFunction mle(sample_mean,S,N);
      opt.set_min_objective(MaximumLikelihoodObjectiveFunction::wrap, &mle);
      opt.add_inequality_constraint(Constraint5, NULL, TOLERANCE);
      opt.set_xtol_rel(1e-4);

      x[0] = alpha; x[1] = eta; x[2] = psi; x[3] = kappa; x[4] = beta;
      nlopt::result result = opt.optimize(x, minf);
      break;
    }

    case MAP:
    {
      opt.set_lower_bounds(lb);
      ub[0] = PI; ub[1] = 2*PI; ub[2] = 2*PI;
      opt.set_upper_bounds(ub);

      MAPObjectiveFunction map(sample_mean,S,N);
      opt.set_min_objective(MAPObjectiveFunction::wrap, &map);
      opt.add_inequality_constraint(Constraint5, NULL, TOLERANCE);
      opt.set_xtol_rel(1e-4);

      x[0] = alpha; x[1] = eta; x[2] = psi; x[3] = kappa; x[4] = beta;
      nlopt::result result = opt.optimize(x, minf);
      break;
    }

    case MML_2:
    {
      opt.set_lower_bounds(lb);
      opt.set_upper_bounds(ub);

      MMLObjectiveFunctionScale mml(psi,alpha,eta,kappa,beta,sample_mean,S,N);
      opt.set_min_objective(MMLObjectiveFunctionScale::wrap, &mml);
      opt.add_inequality_constraint(Constraint2, NULL, TOLERANCE);
      opt.set_xtol_rel(1e-4);

      x[0] = kappa; x[1] = beta;
      nlopt::result result = opt.optimize(x, minf);
      break;
      break;
    }

    case MML_5:
    {
      opt.set_lower_bounds(lb);
      ub[0] = PI; ub[1] = 2*PI; ub[2] = 2*PI;
      opt.set_upper_bounds(ub);

      MMLObjectiveFunction mml(sample_mean,S,N);
      opt.set_min_objective(MMLObjectiveFunction::wrap, &mml);
      opt.add_inequality_constraint(Constraint5, NULL, TOLERANCE);
      opt.set_xtol_rel(1e-4);

      x[0] = alpha; x[1] = eta; x[2] = psi; x[3] = kappa; x[4] = beta;
      nlopt::result result = opt.optimize(x, minf);
      break;
    }

    default:
      break;
  }
  //cout << "solution: (" << x[0] << ", " << x[1] << ")\n";
  cout << "minf: " << minf << endl;
  return x;
}

