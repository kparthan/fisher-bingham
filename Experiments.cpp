#include "Support.h"
#include "Experiments.h"
#include "Kent.h"

extern Vector XAXIS,YAXIS,ZAXIS;

Experiments::Experiments(int iterations) : iterations(iterations)
{}

void Experiments::plotBias(long double kappa, long double beta)
{
  std::vector<int> sample_sizes;
  //sample_sizes.push_back(5);
  //sample_sizes.push_back(10);
  //sample_sizes.push_back(20);
  //sample_sizes.push_back(30);
  sample_sizes.push_back(50);
  sample_sizes.push_back(100);
  sample_sizes.push_back(250);
  sample_sizes.push_back(500);
  sample_sizes.push_back(1000);

  Kent kent(XAXIS,YAXIS,ZAXIS,kappa,beta);
  long double kappa_est,beta_est,diffk,diffb;
  std::vector<struct Estimates> all_estimates;

  string kappa_str = boost::lexical_cast<string>(kappa);
  string beta_str = boost::lexical_cast<string>(beta);
  string tmp = "k_" + kappa_str + "_b_" + beta_str;

  string folder = "./experiments/bias_tests/" + tmp + "/";
  string medians_file = folder + "median_kappa";
  ofstream mediansk(medians_file.c_str(),ios::app);
  string variance_file = folder + "variance_kappa";
  ofstream variancek(variance_file.c_str(),ios::app);
  medians_file = folder + "median_beta";
  ofstream mediansb(medians_file.c_str(),ios::app);
  variance_file = folder + "variance_beta";
  ofstream varianceb(variance_file.c_str(),ios::app);

  for (int i=0; i<sample_sizes.size(); i++) {
    string outputk = folder + "n_" + boost::lexical_cast<string>(sample_sizes[i])
                     + "_k_" + kappa_str;
    string outputb = folder + "n_" + boost::lexical_cast<string>(sample_sizes[i])
                     + "_b_" + beta_str;
    ofstream logk(outputk.c_str());
    ofstream logb(outputb.c_str());
    std::vector<std::vector<long double> > kappa_est_all(NUM_METHODS),errors_kappa(NUM_METHODS);
    std::vector<std::vector<long double> > beta_est_all(NUM_METHODS),errors_beta(NUM_METHODS);
    for (int iter=1; iter<=iterations; iter++) {
      std::vector<Vector > data = kent.generate(sample_sizes[i]);
      Kent kent_est;
      kent_est.computeAllEstimators(data,all_estimates);
      logk << fixed << setw(10) << sample_sizes[i] << "\t";
      logb << fixed << setw(10) << sample_sizes[i] << "\t";
      for (int j=0; j<NUM_METHODS; j++) {
        kappa_est = all_estimates[j].kappa;
        kappa_est_all[j].push_back(kappa_est);
        diffk = fabs(kappa_est - kappa);
        errors_kappa[j].push_back(diffk);
        logk << scientific << kappa_est << "\t";
        beta_est = all_estimates[j].beta;
        beta_est_all[j].push_back(beta_est);
        diffb = fabs(beta_est - beta);
        errors_beta[j].push_back(diffb);
        logb << scientific << beta_est << "\t";
      } // j loop ends ...
      logk << endl; logb << endl;
    } // iter loop ends ..
    logk.close(); logb.close();
    mediansk << fixed << setw(10) << sample_sizes[i] << "\t";
    variancek << fixed << setw(10) << sample_sizes[i] << "\t";
    mediansb << fixed << setw(10) << sample_sizes[i] << "\t";
    varianceb << fixed << setw(10) << sample_sizes[i] << "\t";
    for (int j=0; j<NUM_METHODS; j++) {
      long double median_kappa = computeMean(kappa_est_all[j]);
      long double var_kappa = computeVariance(errors_kappa[j]);
      mediansk << scientific << median_kappa << "\t";
      variancek << scientific << var_kappa << "\t";
      long double median_beta = computeMean(beta_est_all[j]);
      long double var_beta = computeVariance(errors_beta[j]);
      mediansb << scientific << median_beta << "\t";
      varianceb << scientific << var_beta << "\t";
    } // j loop ends ...
    mediansk << endl; mediansb << endl;
    variancek << endl;  varianceb << endl;
  } // i loop ends ...
  mediansk.close(); mediansb.close();
  variancek.close();  varianceb.close();
}

