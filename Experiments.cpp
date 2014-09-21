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
  sample_sizes.push_back(10);
  /*sample_sizes.push_back(20);
  sample_sizes.push_back(30);
  sample_sizes.push_back(50);
  sample_sizes.push_back(100);
  sample_sizes.push_back(250);
  sample_sizes.push_back(500);
  sample_sizes.push_back(1000);*/

  Kent kent(XAXIS,YAXIS,ZAXIS,kappa,beta);
  long double kappa_est,beta_est,diffk,diffb;
  std::vector<struct Estimates> all_estimates;

  string kappa_str = boost::lexical_cast<string>(kappa);
  string beta_str = boost::lexical_cast<string>(beta);
  string tmp = "k_" + kappa_str + "_b_" + beta_str;

  string folder = "./experiments/bias_tests/" + tmp + "/";
  string bias_file = folder + "bias_sq_kappa";
  ofstream biask(bias_file.c_str());
  string variance_file = folder + "variance_kappa";
  ofstream variancek(variance_file.c_str());
  string error_file = folder + "mean_abs_error_kappa";
  ofstream abs_error_k(error_file.c_str());
  error_file = folder + "mean_sqd_error_kappa";
  ofstream sqd_error_k(error_file.c_str());
  string medians_file = folder + "medians_kappa";
  ofstream mediansk(medians_file.c_str());
  string means_file = folder + "means_kappa";
  ofstream meansk(means_file.c_str());

  bias_file = folder + "bias_sq_beta";
  ofstream biasb(bias_file.c_str());
  variance_file = folder + "variance_beta";
  ofstream varianceb(variance_file.c_str());
  error_file = folder + "mean_abs_error_beta";
  ofstream abs_error_b(error_file.c_str());
  error_file = folder + "mean_sqd_error_beta";
  ofstream sqd_error_b(error_file.c_str());
  medians_file = folder + "medians_beta";
  ofstream mediansb(medians_file.c_str());
  means_file = folder + "means_beta";
  ofstream meansb(means_file.c_str());

  for (int i=0; i<sample_sizes.size(); i++) { // for each sample size ...
    string outputk = folder + "n_" + boost::lexical_cast<string>(sample_sizes[i])
                     + "_k_" + kappa_str;
    string outputb = folder + "n_" + boost::lexical_cast<string>(sample_sizes[i])
                     + "_b_" + beta_str;
    ofstream logk(outputk.c_str());
    ofstream logb(outputb.c_str());
    Vector tmp(NUM_METHODS,0);
    std::vector<Vector> kappa_est_all(iterations,tmp),beta_est_all(iterations,tmp);
    Vector avg_kappa_est(NUM_METHODS,0),avg_beta_est(NUM_METHODS,0);
    for (int iter=0; iter<iterations; iter++) {  // for each iteration ...
      repeat:
      std::vector<Vector> data = kent.generate(sample_sizes[i]);
      Kent kent_est;
      kent_est.computeAllEstimators(data,all_estimates);
      long double beta_est_mml = all_estimates[MML_5].beta;
      if (all_estimates[MML_5].beta <= 1e-5) {
        cout << "*** IGNORING ITERATION ***\n";
        goto repeat;
      } else {  // all good with optimization ...
        logk << fixed << setw(10) << sample_sizes[i] << "\t";
        logb << fixed << setw(10) << sample_sizes[i] << "\t";
        for (int j=0; j<NUM_METHODS; j++) { // for each method ...
          // beta_est
          beta_est = all_estimates[j].beta;
          beta_est_all[iter][j] = beta_est;
          logb << scientific << beta_est << "\t";
          avg_beta_est[j] += beta_est;
          // kappa_est
          kappa_est = all_estimates[j].kappa;
          kappa_est_all[iter][j] = kappa_est;
          logk << scientific << kappa_est << "\t";
          avg_kappa_est[j] += kappa_est;
        } // j loop ends ...
        logk << endl; logb << endl;
      } // if() ends ...
    } // iter loop ends ..
    logk.close(); logb.close();
    // compute errors
    for (int j=0; j<NUM_METHODS; j++) {
      avg_kappa_est[j] /= iterations;
      avg_beta_est[j] /= iterations;
      //cout << "avg_kappa: " << avg_kappa_est[j] << endl;
      //cout << "avg_beta: " << avg_beta_est[j] << endl;
    }
    biask << fixed << setw(10) << sample_sizes[i] << "\t";
    variancek << fixed << setw(10) << sample_sizes[i] << "\t";
    abs_error_k << fixed << setw(10) << sample_sizes[i] << "\t";
    sqd_error_k << fixed << setw(10) << sample_sizes[i] << "\t";
    mediansk << fixed << setw(10) << sample_sizes[i] << "\t";
    meansk << fixed << setw(10) << sample_sizes[i] << "\t";
    biasb << fixed << setw(10) << sample_sizes[i] << "\t";
    varianceb << fixed << setw(10) << sample_sizes[i] << "\t";
    abs_error_b << fixed << setw(10) << sample_sizes[i] << "\t";
    sqd_error_b << fixed << setw(10) << sample_sizes[i] << "\t";
    mediansb << fixed << setw(10) << sample_sizes[i] << "\t";
    meansb << fixed << setw(10) << sample_sizes[i] << "\t";
    for (int j=0; j<NUM_METHODS; j++) {
      long double diff_kappa = kappa - avg_kappa_est[j];
      long double bias_sq_kappa = (diff_kappa * diff_kappa)/(long double)iterations;
      long double variance_kappa=0,mean_sqd_error_kappa=0,mean_abs_error_kappa=0;
      long double diff_beta = beta - avg_beta_est[j];
      long double bias_sq_beta = (diff_beta * diff_beta)/(long double)iterations;
      long double variance_beta=0,mean_sqd_error_beta=0,mean_abs_error_beta=0;
      for (int iter=0; iter<iterations; iter++) {
        // compute error terms for kappa estimation
        diff_kappa = avg_kappa_est[j] - kappa_est_all[iter][j];
        variance_kappa += (diff_kappa * diff_kappa);
        diff_kappa = fabs(kappa - kappa_est_all[iter][j]);
        mean_sqd_error_kappa += (diff_kappa * diff_kappa);
        mean_abs_error_kappa += diff_kappa;
        // compute error terms for beta estimation
        diff_beta = avg_beta_est[j] - beta_est_all[iter][j];
        variance_beta += (diff_beta * diff_beta);
        diff_beta = fabs(beta - beta_est_all[iter][j]);
        mean_sqd_error_beta += (diff_beta * diff_beta);
        mean_abs_error_beta += diff_beta;
      } // iter loop ends ...
      variance_kappa /= iterations;
      mean_sqd_error_kappa /= iterations;
      mean_abs_error_kappa /= iterations;
      variance_beta /= iterations;
      mean_sqd_error_beta /= iterations;
      mean_abs_error_beta /= iterations;

      biask << scientific << bias_sq_kappa << "\t";
      variancek << scientific << variance_kappa << "\t";
      abs_error_k << scientific << mean_abs_error_kappa << "\t";
      sqd_error_k << scientific << mean_sqd_error_kappa << "\t";
      biasb << scientific << bias_sq_beta << "\t";
      varianceb << scientific << variance_beta << "\t";
      abs_error_b << scientific << mean_abs_error_beta << "\t";
      sqd_error_b << scientific << mean_sqd_error_beta << "\t";
    } // j loop ends ...
    biask << endl; biasb << endl;
    variancek << endl; varianceb << endl;
    abs_error_k << endl; abs_error_b << endl;
    sqd_error_k << endl; sqd_error_b << endl;
    // compute medians/means of kappa estimates ...
    Vector medians_kappa = computeMedians(kappa_est_all);
    Vector medians_beta = computeMedians(beta_est_all);
    Vector means_kappa = computeMeans(kappa_est_all);
    Vector means_beta = computeMeans(beta_est_all);
    for (int j=0; j<NUM_METHODS; j++) {
      mediansk << scientific << medians_kappa[j] << "\t";
      mediansb << scientific << medians_beta[j] << "\t";
      meansk << scientific << means_kappa[j] << "\t";
      meansb << scientific << means_beta[j] << "\t";
    } // j loop ends ...
    mediansk << endl; mediansb << endl;
    meansk << endl; meansb << endl;
  } // i loop ends ...
  biask.close(); biasb.close();
  variancek.close(); varianceb.close();
  abs_error_k.close(); abs_error_b.close();
  sqd_error_k.close(); sqd_error_b.close();
  mediansk.close(); mediansb.close();
  meansk.close(); meansb.close();
}

