#include "Optimize2.h"

extern int CONSTRAIN_KAPPA;
extern double MAX_KAPPA;

Optimize2::Optimize2(string type)
{
  if (type.compare("MLE") == 0) {
    estimation = MLE;
  } else if (type.compare("MAP") == 0) {
    estimation = MAP;
  } else if (type.compare("MML_2") == 0) {
    estimation = MML_2;
  }
}

void Optimize2::initialize(double sample_size, double magnitude, Vector &mu, double k)
{
  N = sample_size;
  R = magnitude;
  mean = mu;
  kappa = k;
  if (CONSTRAIN_KAPPA == SET) {
    if (kappa >= MAX_KAPPA) {
      kappa = MAX_KAPPA - TOLERANCE;
    }
  }
}

void Optimize2::computeEstimates(struct Estimates_vMF &estimates)
{
  std::vector<double> x(1);
  nlopt::opt opt(nlopt::LN_COBYLA, 1);

  std::vector<double> lb(1,TOLERANCE);
  std::vector<double> ub(1,HUGE_VAL);
  //std::vector<double> ub(1,MAX_KAPPA);

  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  opt.set_xtol_rel(1e-4);
  x[0] = kappa;

  double minf;
  FAIL_STATUS = 0;

  switch(estimation) {
    case MLE:
    {
      MaximumLikelihoodObjectiveFunction_vMF mle(N,R);
      opt.set_min_objective(MaximumLikelihoodObjectiveFunction_vMF::wrap, &mle);
      nlopt::result result = opt.optimize(x, minf);
      break;
    }

    case MAP:
    {
      MAPObjectiveFunction_vMF map(N,R,mean);
      opt.set_min_objective(MAPObjectiveFunction_vMF::wrap, &map);
      nlopt::result result = opt.optimize(x, minf);
      break;
    }

    case MML_2:
    {
      MMLObjectiveFunction_vMF mml2(N,R,mean);
      opt.set_min_objective(MMLObjectiveFunction_vMF::wrap, &mml2);
      nlopt::result result = opt.optimize(x, minf);
      break;
    }

    default:
      break;
  }
  if (boost::math::isnan(minf)) {
    x[0] = kappa;
  }
  //cout << "solution: (" << x[0] << ")\n";
  cout << "minf: " << minf << endl;
  estimates.kappa = x[0];
}

