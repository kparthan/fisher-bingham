#include "Support.h"
#include "Experiments.h"
#include "Kent.h"

extern Vector XAXIS,YAXIS,ZAXIS;

Experiments::Experiments(int iterations) : iterations(iterations)
{}

void Experiments::simulate(long double kappa, long double beta)
{
  std::vector<int> sample_sizes;
  //sample_sizes.push_back(5);
  sample_sizes.push_back(100);
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
  string folder = "./experiments/single_kent/" + tmp + "/";

  for (int i=0; i<sample_sizes.size(); i++) { // for each sample size ...
    string size_str = boost::lexical_cast<string>(sample_sizes[i]);
    string kappas_file = folder + "n_" + size_str + "_kappas";
    string betas_file = folder + "n_" + size_str + "_betas";
    string negloglkhd_file = folder + "n_" + size_str + "_negloglikelihood";
    string kldvg_file = folder + "n_" + size_str + "_kldiv";
    string msglens_file = folder + "n_" + size_str + "_msglens";
    ofstream logk(kappas_file.c_str());
    ofstream logb(betas_file.c_str());
    ofstream logneg(negloglkhd_file.c_str());
    ofstream logkldiv(kldvg_file.c_str());
    ofstream logmsg(msglens_file.c_str());

    Vector emptyvec(NUM_METHODS,0);
    std::vector<Vector> kappa_est_all(iterations,emptyvec),beta_est_all(iterations,emptyvec);
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
        logneg << fixed << setw(10) << sample_sizes[i] << "\t";
        logkldiv << fixed << setw(10) << sample_sizes[i] << "\t";
        logmsg << fixed << setw(10) << sample_sizes[i] << "\t";
        long double actual_negloglkhd = kent.computeNegativeLogLikelihood(data);
        long double actual_msglen = kent.computeMessageLength(data);
        logneg << scientific << actual_negloglkhd << "\t";
        logmsg << scientific << actual_msglen << "\t";
        for (int j=0; j<NUM_METHODS; j++) { // for each method ...
          Kent fit(all_estimates[j].mean,all_estimates[j].major_axis,
          all_estimates[j].minor_axis,all_estimates[j].kappa,all_estimates[j].beta);
          long double negloglkhd = fit.computeNegativeLogLikelihood(data);
          long double msglen = fit.computeMessageLength(data);
          long double kldiv = kent.computeKLDivergence(fit);
          logneg << scientific << negloglkhd << "\t";
          logmsg << scientific << msglen << "\t";
          logkldiv << scientific << kldiv << "\t";
          // beta_est
          beta_est = all_estimates[j].beta;
          beta_est_all[iter][j] = beta_est;
          logb << scientific << beta_est << "\t";
          // kappa_est
          kappa_est = all_estimates[j].kappa;
          kappa_est_all[iter][j] = kappa_est;
          logk << scientific << kappa_est << "\t";
        } // j loop ends ...
        logk << endl; logb << endl; logneg << endl; logmsg << endl; logkldiv << endl;
      } // if() ends ...
    } // iter loop ends ..
    logk.close(); logb.close(); logneg.close(); logmsg.close(); logkldiv.close();
    computeMeasures(kappa,beta,kappa_est_all,beta_est_all,sample_sizes[i]);
  }
}

void 
Experiments::computeMeasures(
  long double kappa,
  long double beta,
  std::vector<Vector> &kappa_est_all,
  std::vector<Vector> &beta_est_all,
  int N
) {
  long double kappa_est,beta_est,diffk,diffb;

  string kappa_str = boost::lexical_cast<string>(kappa);
  string beta_str = boost::lexical_cast<string>(beta);
  string tmp = "k_" + kappa_str + "_b_" + beta_str;
  string folder = "./experiments/single_kent/" + tmp + "/";

  string medians_file = folder + "medians_kappa";
  ofstream mediansk(medians_file.c_str(),ios::app);
  mediansk << fixed << setw(10) << N << "\t";
  Vector medians_kappa = computeEstimateMedians(mediansk,kappa_est_all);
  mediansk.close();

  medians_file = folder + "medians_beta";
  ofstream mediansb(medians_file.c_str(),ios::app);
  mediansb << fixed << setw(10) << N << "\t";
  Vector medians_beta = computeEstimateMedians(mediansb,beta_est_all);
  mediansb.close();

  string means_file = folder + "means_kappa";
  ofstream meansk(means_file.c_str(),ios::app);
  meansk << fixed << setw(10) << N << "\t";
  Vector means_kappa = computeEstimateMeans(meansk,kappa_est_all);
  meansk.close();

  means_file = folder + "means_beta";
  ofstream meansb(means_file.c_str(),ios::app);
  meansb << fixed << setw(10) << N << "\t";
  Vector means_beta = computeEstimateMeans(meansb,beta_est_all);
  meansb.close();

  string bias_file = folder + "bias_sq_kappa";
  ofstream biask(bias_file.c_str(),ios::app);
  biask << fixed << setw(10) << N << "\t";
  computeBias(biask,kappa,kappa_est_all);
  biask.close();

  bias_file = folder + "bias_sq_beta";
  ofstream biasb(bias_file.c_str(),ios::app);
  biasb << fixed << setw(10) << N << "\t";
  computeBias(biasb,beta,beta_est_all);
  biasb.close();

  string variance_file = folder + "variance_kappa";
  ofstream variancek(variance_file.c_str(),ios::app);
  variancek << fixed << setw(10) << N << "\t";
  computeVariance(variancek,kappa,kappa_est_all);
  variancek.close();

  variance_file = folder + "variance_beta";
  ofstream varianceb(variance_file.c_str(),ios::app);
  variancek << fixed << setw(10) << N << "\t";
  computeVariance(varianceb,beta,beta_est_all);
  varianceb.close();

  string error_file = folder + "mean_abs_error_kappa";
  ofstream abs_error_k(error_file.c_str(),ios::app);
  computeMeanAbsoluteError(abs_error_k,kappa,kappa_est_all);
  abs_error_k.close();

  error_file = folder + "mean_abs_error_beta";
  ofstream abs_error_b(error_file.c_str(),ios::app);
  computeMeanAbsoluteError(abs_error_b,beta,beta_est_all);
  abs_error_b.close();

  error_file = folder + "mean_sqd_error_kappa";
  ofstream sqd_error_k(error_file.c_str(),ios::app);
  computeMeanSquaredError(sqd_error_k,kappa,kappa_est_all);
  sqd_error_k.close();

  error_file = folder + "mean_sqd_error_beta";
  ofstream sqd_error_b(error_file.c_str(),ios::app);
  computeMeanSquaredError(sqd_error_b,beta,beta_est_all);
  sqd_error_b.close();
}

Vector Experiments::computeEstimateMedians(ostream &out, std::vector<Vector> &p_est_all)
{
  Vector medians = computeMedians(p_est_all);
  for (int i=0; i<NUM_METHODS; i++) {
    out << scientific << medians[i] << "\t";
  }
  out << endl;
  return medians;
}

Vector Experiments::computeEstimateMeans(ostream &out, std::vector<Vector> &p_est_all)
{
  Vector means = computeMeans(p_est_all);
  for (int i=0; i<NUM_METHODS; i++) {
    out << scientific << means[i] << "\t";
  }
  out << endl;
  return means;
}

void Experiments::computeBias(ostream &out, long double p, std::vector<Vector> &p_est_all)
{
  Vector avg_p_est = computeMeans(p_est_all);

  Vector bias_sq(NUM_METHODS,0);
  for (int j=0; j<NUM_METHODS; j++) {
    bias_sq[j] = (p - avg_p_est[j]) * (p - avg_p_est[j]);
    out << scientific << bias_sq[j] << "\t";
  }
  out << endl;
}

void Experiments::computeVariance(ostream &out, long double p, std::vector<Vector> &p_est_all)
{
  Vector avg_p_est = computeMeans(p_est_all);
  int num_elements = p_est_all.size();

  Vector variance(NUM_METHODS,0);
  for (int i=0; i<num_elements; i++) {
    for (int j=0; j<NUM_METHODS; j++) {
      variance[j] += (avg_p_est[j] - p_est_all[i][j]) * (avg_p_est[j] - p_est_all[i][j]);
    }
  }
  for (int j=0; j<NUM_METHODS; j++) {
    variance[j] /= num_elements;
    out << scientific << variance[j] << "\t";
  }
  out << endl;
}

void Experiments::computeMeanAbsoluteError(ostream &out, long double p, std::vector<Vector> &p_est_all)
{
  int num_elements = p_est_all.size();

  Vector error(NUM_METHODS,0);
  for (int i=0; i<num_elements; i++) {
    for (int j=0; j<NUM_METHODS; j++) {
      error[j] += fabs(p - p_est_all[i][j]);
    }
  }
  for (int j=0; j<NUM_METHODS; j++) {
    error[j] /= num_elements;
    out << scientific << error[j] << "\t";
  }
  out << endl;
}

void Experiments::computeMeanSquaredError(ostream &out, long double p, std::vector<Vector> &p_est_all)
{
  int num_elements = p_est_all.size();

  Vector error(NUM_METHODS,0);
  for (int i=0; i<num_elements; i++) {
    for (int j=0; j<NUM_METHODS; j++) {
      error[j] += (p - p_est_all[i][j]) * (p - p_est_all[i][j]);
    }
  }
  for (int j=0; j<NUM_METHODS; j++) {
    error[j] /= num_elements;
    out << scientific << error[j] << "\t";
  }
  out << endl;
}

