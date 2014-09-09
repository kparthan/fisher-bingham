#include "Support.h"
#include "Experiments.h"
#include "Kent.h"

Experiments::Experiments(int iterations) : iterations(iterations)
{}

void Experiments::plotBias(long double kappa, long double beta)
{
  vector<int> sample_sizes;
  sample_sizes.push_back(5);
  sample_sizes.push_back(10);
  sample_sizes.push_back(20);
  sample_sizes.push_back(30);
  sample_sizes.push_back(50);
  sample_sizes.push_back(100);
  sample_sizes.push_back(500);
  sample_sizes.push_back(1000);

  Kent kent(kappa,beta);
  long double kappa_est,beta_est,diffk,diffb;
  std::vector<Kent::Estimates> all_estimates;

  string kappa_str = boost::lexical_cast<string>(kappa);
  string beta_str = boost::lexical_cast<string>(beta);
  string tmp = "k_" + kappa_str + "_b_" + beta_str;

  string folder = "./experiments/bias_tests/" + tmp + "/";
  string medians_file = folder + "median_kappas";
  ofstream medians(medians_file.c_str());
  string bias_file = folder + "bias";
  ofstream bias(bias_file.c_str());
  string variance_file = folder + "variance";
  ofstream variance(variance_file.c_str());

  for (int i=0; i<sample_sizes.size(); i++) {
    string output_kappa = folder + "_n_" + boost::lexical_cast<string>(sample_sizes[i])
                          + "_k_" + kappa_str;
    string output_beta = folder + "_n_" + boost::lexical_cast<string>(sample_sizes[i])
                         + "_k_" + beta_str;
    ofstream log_kappa(output_kappa.c_str());
    ofstream log_beta(output_beta.c_str());
    std::vector<std::vector<long double> > kappa_est_all(NUM_METHODS),errors_kappa(NUM_METHODS);
    std::vector<std::vector<long double> > beta_est_all(NUM_METHODS),errors_beta(NUM_METHODS);
    for (int iter=1; iter<=iterations; iter++) {
      vector<vector<long double> > data = kent.generate(sample_sizes[i]);
      Kent kent_est;
      kent_est.computeAllEstimators(data,all_estimates);
      log_kappa << fixed << setw(10) << sample_sizes[i] << "\t";
      for (int j=0; j<NUM_METHODS; j++) {
        
      } // j loop ends ...
      log_kappa << endl;
    } // iter loop ends ..
    log_file.close();
    medians << fixed << setw(10) << sample_sizes[i] << "\t";
    bias << fixed << setw(10) << sample_sizes[i] << "\t";
    variance << fixed << setw(10) << sample_sizes[i] << "\t";
    for (int j=0; j<NUM_METHODS; j++) {
      long double median_kappa = computeMedian(kappa_est_all[j]);
      long double mean_error = computeMean(errors[j]);
      long double var = computeVariance(errors[j]);
      medians << scientific << median_kappa << "\t";
      bias << scientific << mean_error << "\t";
      variance << scientific << var << "\t";
    }
    medians << endl;
    bias << endl;
    variance << endl;
  } // i loop ends ...
  medians.close();
  bias.close(); 
  variance.close();
}
