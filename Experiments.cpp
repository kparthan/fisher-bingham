#include "Support.h"
#include "Experiments.h"
#include "Kent.h"

extern Vector XAXIS,YAXIS,ZAXIS;
extern int MOMENT_FAIL,MLE_FAIL,MAP_FAIL,MML_FAIL;
extern struct stat st;

Experiments::Experiments(int iterations) : iterations(iterations)
{}

void Experiments::simulate()
{
  int N = 10;
  double kappa,beta,eccentricity;

  string n_str = "N_" + boost::lexical_cast<string>(N) + "_uniform_prior";
  //string n_str = "N_" + boost::lexical_cast<string>(N) + "_vmf_prior";
  //string n_str = "N_" + boost::lexical_cast<string>(N) + "_beta_prior";
  string parent_dir = "./experiments/single_kent/" + n_str + "/";
  string current_dir,kappa_str,eccentricity_str;
  string kappas,betas,negloglike,kldivs,msglens;
  std::vector<Vector> random_sample;
  std::vector<struct Estimates> all_estimates;
  //string data_file = "random_sample_uniform.dat";
  //string data_file = "random_sample_vmf.dat";
  //string data_file = "random_sample_beta.dat";
  string data_file = "random_sample.dat";

  kappa = 50;
  while (kappa <= 100) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    kappa_str = ssk.str();
    eccentricity = 0.1;
    while (eccentricity <= 0.95) {
      beta = 0.5 * kappa * eccentricity;
      ostringstream sse;
      sse << fixed << setprecision(1);
      sse << eccentricity;
      eccentricity_str = sse.str();
      current_dir = parent_dir + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
      if (stat(current_dir.c_str(), &st) == -1) {
          mkdir(current_dir.c_str(), 0700);
      }
      kappas = current_dir + "kappas";
      betas = current_dir + "betas";
      negloglike = current_dir + "negloglike";
      kldivs = current_dir + "kldivs";
      msglens = current_dir + "msglens";

      ofstream fk(kappas.c_str());
      ofstream fb(betas.c_str());
      ofstream fnlh(negloglike.c_str());
      ofstream fkl(kldivs.c_str());
      ofstream fmsg(msglens.c_str());
      
      cout << "kappa: " << kappa << "; beta: " << beta << "; e: " << eccentricity << endl;
      //Kent kent(ZAXIS,XAXIS,YAXIS,kappa,beta);
      double psi = 60; psi *= PI/180;
      double alpha = 60; alpha *= PI/180;
      double eta = 70; eta *= PI/180;
      Kent kent(psi,alpha,eta,kappa,beta);

      for (int i=0; i<iterations; i++) {
        cout << "Iteration: " << i+1 << endl;
        kent.generate(N);
        random_sample = load_matrix(data_file);
        //random_sample = kent.generate(N);
        kent.computeAllEstimators(random_sample,all_estimates,0,1);
        for (int j=0; j<all_estimates.size(); j++) {
          fk << scientific << all_estimates[j].kappa << "\t";
          fb << scientific << all_estimates[j].beta << "\t";
          fnlh << scientific << all_estimates[j].negloglike << "\t";
          fkl << scientific << all_estimates[j].kldiv << "\t";
          fmsg << scientific << all_estimates[j].msglen << "\t";
        } // for j ()
        fk << endl;
        fb << endl;
        fnlh << endl;
        fkl << endl;
        fmsg << endl;
      } // for i ()

      eccentricity += 0.1;
      fk.close();
      fb.close();
      fnlh.close();
      fkl.close();
      fmsg.close();
    } // eccentricity
    kappa += 10;
  } // kappa
}

void Experiments::computeMeasures(
  double kappa,
  double beta,
  std::vector<Vector> &kappa_est_all,
  std::vector<Vector> &beta_est_all,
  int N
) {
  double kappa_est,beta_est,diffk,diffb;

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
  varianceb << fixed << setw(10) << N << "\t";
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

void Experiments::computeBias(ostream &out, double p, std::vector<Vector> &p_est_all)
{
  Vector avg_p_est = computeMeans(p_est_all);

  Vector bias_sq(NUM_METHODS,0);
  for (int j=0; j<NUM_METHODS; j++) {
    bias_sq[j] = (p - avg_p_est[j]) * (p - avg_p_est[j]);
    out << scientific << bias_sq[j] << "\t";
  }
  out << endl;
}

void Experiments::computeVariance(ostream &out, double p, std::vector<Vector> &p_est_all)
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

void Experiments::computeMeanAbsoluteError(
  ostream &out, double p, std::vector<Vector> &p_est_all
) {
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

void Experiments::computeMeanSquaredError(
  ostream &out, double p, std::vector<Vector> &p_est_all
) {
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

void Experiments::computeWinsRatio(
  bool ignore_first, 
  const char *measure, 
  string input_file, 
  int N,
  string folder
) {
  ifstream input(input_file.c_str());
  // read values from the file
  string line;
  int start;
  if (ignore_first == 0) {
    start = 1;
  } else if (ignore_first == 1) {
    start = 2;
  }
  Vector numbers;
  std::vector<Vector> table;
  while (getline(input,line)) {
    boost::char_separator<char> sep(" \t");
    boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
    BOOST_FOREACH (const string& t, tokens) {
      istringstream iss(t);
      double x;
      iss >> x;
      numbers.push_back(x);
    }
    Vector values_to_compare;
    for (int i=start; i<numbers.size(); i++) {
      values_to_compare.push_back(numbers[i]);
    }
    table.push_back(values_to_compare);
    numbers.clear();
  }
  input.close();

  string wins_file = folder + "wins";
  ofstream out(wins_file.c_str(),ios::app);
  out << setw(10) << N << "\t";
  out << setw(20) << measure << "\t";
  std::vector<int> wins(NUM_METHODS,0);
  for (int i=0; i<table.size(); i++) {
    int winner = 0;
    double min = table[i][0];
    for (int j=1; j<NUM_METHODS; j++) {
      if (table[i][j] <= min) {
        min = table[i][j];
        winner = j;
      }
    } // j loop ends ...
    wins[winner]++;
  } // i loop ends ...
  out << "[";
  for (int i=0; i<wins.size()-1; i++) {
    out << wins[i] << " : ";
  }
  out << wins[wins.size()-1] << "]\n";
  out.close();
}

