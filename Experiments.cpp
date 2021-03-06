#include "Support.h"
#include "Experiments.h"
#include "Kent.h"

extern Vector XAXIS,YAXIS,ZAXIS;
extern int ESTIMATION,CRITERION;
extern int INFER_COMPONENTS;
extern string EM_LOG_FOLDER;
extern struct stat st;
extern int MML_FAIL;
extern int PRIOR;
extern int NUM_THREADS;
extern int ENABLE_DATA_PARALLELISM;

Experiments::Experiments(int iterations) : iterations(iterations)
{}

void Experiments::fisher_uncertainty()
{
  int N = 10;

  double kappa = 100;
  for (double ecc=0.1; ecc<0.95; ecc+=0.1) {
    double beta = 0.5 * kappa * ecc;
    Kent kent(ZAXIS,XAXIS,YAXIS,kappa,beta);
    kent.computeExpectation();
    //double log_det_fisher = kent.computeLogFisherScale();
    double log_det_fisher = kent.computeLogFisherInformation(N);
    double log_vkb = -2.5*logLatticeConstant(5) - 0.5 * log_det_fisher;
    cout << "log_vkb: " << log_vkb << endl;
  }
}

void Experiments::simulate()
{
  int N = 20;
  double kappa,beta,eccentricity;

  //string n_str = "N_" + boost::lexical_cast<string>(N) + "_uniform_prior/";
  //string n_str = "N_" + boost::lexical_cast<string>(N) + "_vmf_prior/";
  //string n_str = "N_" + boost::lexical_cast<string>(N) + "_beta_prior/";
  //string n_str = "N_" + boost::lexical_cast<string>(N) + "_new_prior/";
  //string n_str = "N_" + boost::lexical_cast<string>(N) + "_new2_prior/";
  string n_str = "N_" + boost::lexical_cast<string>(N) 
                 + "_prior" + boost::lexical_cast<string>(PRIOR) + "/";
  string parent_dir = "experiments/single_kent/estimates/" + n_str + "/";
  check_and_create_directory(parent_dir);

  string current_dir,kappa_str,eccentricity_str;
  string psi_est,alpha_est,eta_est,kappas,betas,ecc_file;
  string negloglike,kldivs,msglens,chisq_stat,pvalues_file;
  std::vector<Vector> random_sample;
  std::vector<struct Estimates> all_estimates;
  Vector statistics,pvalues;

  //string data_file = "random_sample_uniform.dat";
  //string data_file = "random_sample_vmf.dat";
  //string data_file = "random_sample_beta.dat";
  //string data_file = "random_sample.dat";
  //string data_file = "random_sample_new2.dat";

  double INIT_KAPPA = 1;
  double MAX_KAPPA = 100;
  double KAPPA_INCREMENT = 10;
  double ecc;

  kappa = INIT_KAPPA;
  while (kappa <= MAX_KAPPA + 1) {
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
      check_and_create_directory(current_dir);
      psi_est = current_dir + "psi_est";
      alpha_est = current_dir + "alpha_est";
      eta_est = current_dir + "eta_est";
      kappas = current_dir + "kappa_est";
      betas = current_dir + "beta_est";
      ecc_file = current_dir + "ecc_est";
      negloglike = current_dir + "negloglike";
      kldivs = current_dir + "kldivs";
      msglens = current_dir + "msglens";
      chisq_stat = current_dir + "chisq_stat";
      pvalues_file = current_dir + "pvalues";

      ofstream fpsi(psi_est.c_str(),ios::app);
      ofstream falpha(alpha_est.c_str(),ios::app);
      ofstream feta(eta_est.c_str(),ios::app);
      ofstream fk(kappas.c_str(),ios::app);
      ofstream fb(betas.c_str(),ios::app);
      ofstream fecc(ecc_file.c_str(),ios::app);
      ofstream fnlh(negloglike.c_str(),ios::app);
      ofstream fkl(kldivs.c_str(),ios::app);
      ofstream fmsg(msglens.c_str(),ios::app);
      ofstream fchi(chisq_stat.c_str(),ios::app);
      ofstream fpval(pvalues_file.c_str(),ios::app);
      
      cout << "kappa: " << kappa << "; beta: " << beta << "; e: " << eccentricity << endl;
      //Kent kent(ZAXIS,XAXIS,YAXIS,kappa,beta);
      double psi = 45; psi *= PI/180;
      double alpha = 90; alpha *= PI/180;
      double eta = 90; eta *= PI/180;
      Kent kent(psi,alpha,eta,kappa,beta);

      for (int i=0; i<iterations; i++) {
        repeat:
        cout << "Iteration: " << i+1 << endl;
        random_sample = kent.generate(N);
        kent.computeAllEstimators(random_sample,all_estimates,0,1);
        Vector msglens(all_estimates.size(),0);
        for (int j=0; j<all_estimates.size(); j++) {
          if (all_estimates[j].kldiv < 0 || all_estimates[j].msglen < 0) goto repeat;
          msglens[j] = all_estimates[j].msglen;
        }
        int min_index = minimumIndex(msglens);
        if (min_index != MML) {  // ignore iteration
          goto repeat;
        }
        chisquare_hypothesis_testing(random_sample,all_estimates,statistics,pvalues);
        for (int j=0; j<all_estimates.size(); j++) {
          ecc = 2 * all_estimates[j].beta / all_estimates[j].kappa;
          fpsi << scientific << all_estimates[j].psi << "\t";
          falpha << scientific << all_estimates[j].alpha << "\t";
          feta << scientific << all_estimates[j].eta << "\t";
          fecc << scientific << ecc << "\t";
          fk << scientific << all_estimates[j].kappa << "\t";
          fb << scientific << all_estimates[j].beta << "\t";
          fnlh << scientific << all_estimates[j].negloglike << "\t";
          fkl << scientific << all_estimates[j].kldiv << "\t";
          fmsg << scientific << all_estimates[j].msglen << "\t";
          fchi << scientific << statistics[j] << "\t";
          fpval << scientific << pvalues[j] << "\t";
        } // for j ()
        fpsi << endl; falpha << endl; feta << endl;
        fecc << endl; fk << endl; fb << endl;
        fnlh << endl; fkl << endl; fmsg << endl;
        fchi << endl; fpval << endl;
      } // for i ()

      eccentricity += 0.4;
      fecc.close();
      fk.close();
      fb.close();
      fnlh.close();
      fkl.close();
      fmsg.close();
      fchi.close();
      fpval.close();
    } // eccentricity
    //kappa += KAPPA_INCREMENT;
    kappa *= KAPPA_INCREMENT;
  } // kappa
}

void Experiments::checkFolders(string &exp_folder)
{
  string data_folder = exp_folder + "data/";
  check_and_create_directory(data_folder);

  string search = exp_folder + "proposed_search/";
  check_and_create_directory(search);
  search = exp_folder + "traditional_search/";
  check_and_create_directory(search);

  string criterion;
  string logs_folder = exp_folder + "proposed_search/logs/";
  check_and_create_directory(logs_folder);
  criterion = "aic/";
  create_sub_folders(logs_folder,criterion);
  criterion = "bic/";
  create_sub_folders(logs_folder,criterion);
  criterion = "icl/";
  create_sub_folders(logs_folder,criterion);

  string mixtures_folder = exp_folder + "proposed_search/mixtures/";
  check_and_create_directory(mixtures_folder);
  criterion = "aic/";
  create_sub_folders(mixtures_folder,criterion);
  criterion = "bic/";
  create_sub_folders(mixtures_folder,criterion);
  criterion = "icl/";
  create_sub_folders(mixtures_folder,criterion);

  mixtures_folder = exp_folder + "traditional_search/mixtures/";
  check_and_create_directory(mixtures_folder);
  criterion = "aic/";
  create_sub_folders(mixtures_folder,criterion);
  criterion = "bic/";
  create_sub_folders(mixtures_folder,criterion);
  criterion = "icl/";
  create_sub_folders(mixtures_folder,criterion);

  string results_folder = exp_folder + "proposed_search/results/";
  check_and_create_directory(results_folder);
  criterion = "aic/";
  string results = results_folder + criterion;
  check_and_create_directory(results);
  criterion = "bic/";
  results = results_folder + criterion;
  check_and_create_directory(results);
  criterion = "icl/";
  results = results_folder + criterion;
  check_and_create_directory(results);

  results_folder = exp_folder + "traditional_search/results/";
  check_and_create_directory(results_folder);
  criterion = "aic/";
  results = results_folder + criterion;
  check_and_create_directory(results);
  criterion = "bic/";
  results = results_folder + criterion;
  check_and_create_directory(results);
  criterion = "icl/";
  results = results_folder + criterion;
  check_and_create_directory(results);
}

void Experiments::create_sub_folders(string &folder, string &criterion)
{
  string tmpc,tmp;

  tmpc = folder + criterion;
  check_and_create_directory(tmpc);
  tmp = tmpc + "moment/";
  check_and_create_directory(tmp);
  tmp = tmpc + "mle/";
  check_and_create_directory(tmp);
  tmp = tmpc + "map/";
  check_and_create_directory(tmp);
  tmp = folder + "mml/";
  check_and_create_directory(tmp);
}

void Experiments::exp1()
{
  int K = 2;
  int num_mixtures = 50;
  int N = 100;

  //generate_data_exp1(N,K,num_mixtures);

  infer_components_exp1(K,num_mixtures);
}

void Experiments::generate_data_exp1(int N, int K, int num_mixtures)
{
  string common = "./experiments/infer_components/exp1/";
  check_and_create_directory(common);

  string K_str = boost::lexical_cast<string>(K);

  string exp_folder = common + "K_" + K_str + "/";
  check_and_create_directory(exp_folder);
  checkFolders(exp_folder);

  string data_folder = exp_folder + "data/";
  string original_mix_folder = exp_folder + "original/";
  check_and_create_directory(original_mix_folder);

  for (int index=1; index<=num_mixtures; index++) {
    string index_str = boost::lexical_cast<string>(index);

    Vector weights = generateFromSimplex(K);
    std::vector<Kent> components = generateRandomComponents(K);
    Mixture original(K,components,weights);
    string mix_file = original_mix_folder + "mixture_" + index_str;
    original.printParameters(mix_file);

    std::vector<Vector> random_sample = original.generate(N,0);
    string data_file = data_folder + "mixture_" + index_str + ".dat";
    writeToFile(data_file,random_sample);
  } // for()
}

void Experiments::infer_components_exp1(int K, int num_mixtures)
{
  INFER_COMPONENTS = SET;

  string common = "./experiments/infer_components/exp1/";
  string K_str = boost::lexical_cast<string>(K);

  EM_LOG_FOLDER = "./infer/logs/kent/K_" + K_str + "/";
  check_and_create_directory(EM_LOG_FOLDER);

  string exp_folder = common + "K_" + K_str + "/";
  string data_folder = exp_folder + "data/";
  string original_mix_folder = exp_folder + "original/";

  int large_N = 100000;

  for (int index=1; index<=num_mixtures; index++) {
    cout << "Mixture: " << index << " ...\n";
    string index_str = boost::lexical_cast<string>(index);

    string original_mix_file = original_mix_folder + "mixture_" + index_str;
    Mixture original;
    string mix_file = original_mix_folder + "mixture_" + index_str;
    original.load(mix_file);

    std::vector<Vector> large_sample = original.generate(large_N,0);
    string data_file = data_folder + "mixture_" + index_str + ".dat";
    std::vector<Vector> random_sample = load_data_table(data_file);

    /* using traditional search */
    cout << "Traditional search ...\n";
    traditional_search(index_str,original,random_sample,large_sample,exp_folder);

    /* using the proposed search */
    cout << "Proposed search ...\n";
    proposed_search(index_str,original,random_sample,large_sample,exp_folder);
  } // for()
}

void Experiments::exp2()
{
  int K = 5;
  int num_mixtures = 50;
  int N = 50;

  generate_data_exp2(N,K,num_mixtures);

  infer_components_exp2(N,num_mixtures);
}

void Experiments::generate_data_exp2(int N, int K, int num_mixtures)
{
  string common = "./experiments/infer_components/exp2/";
  check_and_create_directory(common);

  string N_str = boost::lexical_cast<string>(N);

  string exp_folder = common + "N_" + N_str + "/";
  check_and_create_directory(exp_folder);
  checkFolders(exp_folder);

  string data_folder = exp_folder + "data/";
  string original_mix_folder = exp_folder + "original/";
  check_and_create_directory(original_mix_folder);

  for (int index=1; index<=num_mixtures; index++) {
    string index_str = boost::lexical_cast<string>(index);

    Vector weights = generateFromSimplex(K);
    std::vector<Kent> components = generateRandomComponents(K);
    Mixture original(K,components,weights);
    string mix_file = original_mix_folder + "mixture_" + index_str;
    original.printParameters(mix_file);

    std::vector<Vector> random_sample = original.generate(N,0);
    string data_file = data_folder + "mixture_" + index_str + ".dat";
    writeToFile(data_file,random_sample);
  } // for()
}

void Experiments::infer_components_exp2(int N, int num_mixtures)
{
  INFER_COMPONENTS = SET;

  string common = "./experiments/infer_components/exp2/";
  string N_str = boost::lexical_cast<string>(N);

  EM_LOG_FOLDER = "./infer/logs/kent/N_" + N_str + "/";
  check_and_create_directory(EM_LOG_FOLDER);

  string exp_folder = common + "N_" + N_str + "/";
  string data_folder = exp_folder + "data/";
  string original_mix_folder = exp_folder + "original/";

  int large_N = 100000;

  for (int index=1; index<=num_mixtures; index++) {
    cout << "Mixture: " << index << " ...\n";
    string index_str = boost::lexical_cast<string>(index);

    string original_mix_file = original_mix_folder + "mixture_" + index_str;
    Mixture original;
    string mix_file = original_mix_folder + "mixture_" + index_str;
    original.load(mix_file);

    std::vector<Vector> large_sample = original.generate(large_N,0);
    string data_file = data_folder + "mixture_" + index_str + ".dat";
    std::vector<Vector> random_sample = load_data_table(data_file);

    /* using traditional search */
    cout << "Traditional search ...\n";
    traditional_search(index_str,original,random_sample,large_sample,exp_folder);

    /* using the proposed search */
    cout << "Proposed search ...\n";
    proposed_search(index_str,original,random_sample,large_sample,exp_folder);
  } // for()
}

void Experiments::traditional_search(
  string &index_str,
  Mixture &original, 
  std::vector<Vector> &random_sample,
  std::vector<Vector> &large_sample,
  string &exp_folder
) {
  ESTIMATION = MOMENT; // Moment estimation ...
  string est_type_file = "moment";
  string est_type_folder = est_type_file + "/";
  string results_file = exp_folder + "traditional_search/results/aic/" + est_type_file;
  ofstream aic_out_mom(results_file.c_str(),ios::app);
  results_file = exp_folder + "traditional_search/results/bic/" + est_type_file;
  ofstream bic_out_mom(results_file.c_str(),ios::app);
  results_file = exp_folder + "traditional_search/results/icl/" + est_type_file;
  ofstream icl_out_mom(results_file.c_str(),ios::app);
  traditional_search(
    index_str,
    original,random_sample,large_sample,
    exp_folder,est_type_folder,
    aic_out_mom,bic_out_mom,icl_out_mom
  );
  aic_out_mom.close(); bic_out_mom.close(); icl_out_mom.close();

  ESTIMATION = MLE; // Max LH estimation ...
  est_type_file = "mle";
  est_type_folder = est_type_file + "/";
  results_file = exp_folder + "traditional_search/results/aic/" + est_type_file;
  ofstream aic_out_mle(results_file.c_str(),ios::app);
  results_file = exp_folder + "traditional_search/results/bic/" + est_type_file;
  ofstream bic_out_mle(results_file.c_str(),ios::app);
  results_file = exp_folder + "traditional_search/results/icl/" + est_type_file;
  ofstream icl_out_mle(results_file.c_str(),ios::app);
  traditional_search(
    index_str,
    original,random_sample,large_sample,
    exp_folder,est_type_folder,
    aic_out_mle,bic_out_mle,icl_out_mle
  );
  aic_out_mle.close(); bic_out_mle.close(); icl_out_mle.close();

  ESTIMATION = MAP; // MAP estimation ...
  est_type_file = "map";
  est_type_folder = est_type_file + "/";
  results_file = exp_folder + "traditional_search/results/aic/" + est_type_file;
  ofstream aic_out_map(results_file.c_str(),ios::app);
  results_file = exp_folder + "traditional_search/results/bic/" + est_type_file;
  ofstream bic_out_map(results_file.c_str(),ios::app);
  results_file = exp_folder + "traditional_search/results/icl/" + est_type_file;
  ofstream icl_out_map(results_file.c_str(),ios::app);
  traditional_search(
    index_str,
    original,random_sample,large_sample,
    exp_folder,est_type_folder,
    aic_out_map,bic_out_map,icl_out_map
  );
  aic_out_map.close(); bic_out_map.close(); icl_out_map.close();

  ESTIMATION = MML; // MML estimation ...
  est_type_file = "mml";
  est_type_folder = est_type_file + "/";
  results_file = exp_folder + "traditional_search/results/" + est_type_file;
  ofstream mml_out(results_file.c_str(),ios::app);
  traditional_search_mml(
    index_str,
    original,random_sample,large_sample,
    exp_folder,est_type_folder,
    mml_out
  );
  mml_out.close();
}

void Experiments::traditional_search(
  string &index_str,
  Mixture &original, 
  std::vector<Vector> &random_sample,
  std::vector<Vector> &large_sample,
  string &exp_folder,
  string &est_type_folder,
  ostream &aic_out,
  ostream &bic_out,
  ostream &icl_out
) {
  int K_MAX = 10; 

  Vector data_weights(random_sample.size(),1);
  Mixture mixture(1,random_sample,data_weights);
  mixture.estimateParameters();

  double aic_best = mixture.computeAIC();
  double bic_best = mixture.computeBIC();
  double icl_best = mixture.computeICL();
  Mixture aic_best_mix = mixture;
  Mixture bic_best_mix = mixture;
  Mixture icl_best_mix = mixture;

  for (int k=2; k<=K_MAX; k++) {
    cout << "k: " << k << endl;
    Mixture mixture(k,random_sample,data_weights);
    mixture.estimateParameters();
    double aic = mixture.computeAIC();
    if (aic < aic_best) {
      aic_best = aic;
      aic_best_mix = mixture;
    }
    double bic = mixture.computeBIC();
    if (bic < bic_best) {
      bic_best = bic;
      bic_best_mix = mixture;
    }
    double icl = mixture.computeICL();
    if (icl < icl_best) {
      icl_best = icl;
      icl_best_mix = mixture;
    }
  } // for()

  string mixtures_folder = exp_folder + "traditional_search/mixtures/";
  /* save inferred mixtures */
  string inferred_mix_file = mixtures_folder + "aic/" + est_type_folder + "mixture_" + index_str;
  aic_best_mix.printParameters(inferred_mix_file);
  aic_out << fixed << setw(10) << aic_best_mix.getNumberOfComponents() << "\t\t";
  aic_out << fixed << scientific << original.computeKLDivergence(aic_best_mix,large_sample) << "\t\t";
  aic_out << fixed << scientific << aic_best_mix.computeMinimumMessageLength() << endl;

  inferred_mix_file = mixtures_folder + "bic/" + est_type_folder + "mixture_" + index_str;  
  bic_best_mix.printParameters(inferred_mix_file);
  bic_out << fixed << setw(10) << bic_best_mix.getNumberOfComponents() << "\t\t";
  bic_out << fixed << scientific << original.computeKLDivergence(bic_best_mix,large_sample) << "\t\t";
  bic_out << fixed << scientific << bic_best_mix.computeMinimumMessageLength() << endl;

  inferred_mix_file = mixtures_folder + "icl/" + est_type_folder + "mixture_" + index_str;  
  icl_best_mix.printParameters(inferred_mix_file);
  icl_out << fixed << setw(10) << icl_best_mix.getNumberOfComponents() << "\t\t";
  icl_out << fixed << scientific << original.computeKLDivergence(icl_best_mix,large_sample) << "\t\t";
  icl_out << fixed << scientific << icl_best_mix.computeMinimumMessageLength() << endl;
}

void Experiments::traditional_search_mml(
  string &index_str,
  Mixture &original, 
  std::vector<Vector> &random_sample,
  std::vector<Vector> &large_sample,
  string &exp_folder,
  string &est_type_folder,
  ostream &mml_out
) {
  int K_MAX = 10; 

  Vector data_weights(random_sample.size(),1);
  Mixture mixture(1,random_sample,data_weights);
  mixture.estimateParameters();

  double msglen_best = mixture.getMinimumMessageLength();
  Mixture mml_best_mix = mixture;

  for (int k=2; k<=K_MAX; k++) {
    cout << "k: " << k << endl;
    Mixture mixture(k,random_sample,data_weights);
    mixture.estimateParameters();
    double msglen = mixture.getMinimumMessageLength();
    if (msglen < msglen_best) {
      msglen_best = msglen;
      mml_best_mix = mixture;
    }
  } // for()

  string mixtures_folder = exp_folder + "traditional_search/mixtures/";
  /* save inferred mixtures */
  string inferred_mix_file = mixtures_folder + est_type_folder + "mixture_" + index_str;
  mml_best_mix.printParameters(inferred_mix_file);
  mml_out << fixed << setw(10) << mml_best_mix.getNumberOfComponents() << "\t\t";
  mml_out << fixed << scientific << original.computeKLDivergence(mml_best_mix,large_sample) << "\t\t";
  mml_out << fixed << scientific << mml_best_mix.getMinimumMessageLength() << endl;
}

void Experiments::proposed_search(
  string &index_str,
  Mixture &original, 
  std::vector<Vector> &random_sample,
  std::vector<Vector> &large_sample,
  string &exp_folder
) {
  string logs_folder = exp_folder + "proposed_search/logs/";
  string mixtures_folder = exp_folder + "proposed_search/mixtures/";
  string results_folder = exp_folder + "proposed_search/results/";
  string criterion;

  criterion = "aic/";
  CRITERION = AIC;

  ESTIMATION = MOMENT;
  proposed_search(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  ESTIMATION = MLE;
  proposed_search(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  ESTIMATION = MAP;
  proposed_search(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  criterion = "bic/";
  CRITERION = BIC;

  ESTIMATION = MOMENT;
  proposed_search(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  ESTIMATION = MLE;
  proposed_search(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  ESTIMATION = MAP;
  proposed_search(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  criterion = "icl/";
  CRITERION = ICL;

  ESTIMATION = MOMENT;
  proposed_search(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  ESTIMATION = MLE;
  proposed_search(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  ESTIMATION = MAP;
  proposed_search(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  criterion = "mml/";
  CRITERION = MMLC;

  ESTIMATION = MML;
  proposed_search(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );
}

void Experiments::proposed_search(
  string &index_str,
  Mixture &original, 
  std::vector<Vector> &random_sample,
  std::vector<Vector> &large_sample,
  string &logs_folder,
  string &mixtures_folder,
  string &results_folder,
  string &criterion
) {
  string log_file,mixture_file,results_file,iter_str;
  string logs_folder_specific,mixtures_folder_specific;
  Mixture inferred;
  double kldiv;
  
  switch(ESTIMATION) {
    case MOMENT:
      logs_folder_specific = logs_folder + criterion + "moment/";
      mixtures_folder_specific = mixtures_folder + criterion + "moment/";
      results_file = results_folder + criterion + "moment";
      break;

    case MLE:
      logs_folder_specific = logs_folder + criterion + "mle/";
      mixtures_folder_specific = mixtures_folder + criterion + "mle/";
      results_file = results_folder + criterion + "mle";
      break;

    case MAP:
      logs_folder_specific = logs_folder + criterion + "map/";
      mixtures_folder_specific = mixtures_folder + criterion + "map/";
      results_file = results_folder + criterion + "map";
      break;

    case MML:
      logs_folder_specific = logs_folder + "mml/";
      mixtures_folder_specific = mixtures_folder + "mml/";
      results_file = results_folder + "mml";
      break;
  } // switch() ends ...

  // for empirical KL-divergence ...
  ofstream out(results_file.c_str(),ios::app);

  log_file = logs_folder_specific + "mixture_" + index_str + ".log";
  mixture_file = mixtures_folder_specific + "mixture_" + index_str;

  inferred = inferComponents(random_sample,log_file);
  inferred.printParameters(mixture_file);
  kldiv = original.computeKLDivergence(inferred,large_sample);

  out << setw(10) << inferred.getNumberOfComponents() << "\t\t";
  out << scientific << kldiv << "\t\t";
  out << scientific << inferred.getMinimumMessageLength() << endl;
  out.close();
}

void Experiments::exp3()
{
  struct Parameters parameters;
  parameters.profiles_dir = "./data/profiles-b/";
  NUM_THREADS = 2;

  std::vector<Vector> data;
  gatherData(parameters,data);
  cout << "data.size(): " << data.size() << endl;

  string exp_folder = "./experiments/infer_components/exp3/";
  check_and_create_directory(exp_folder);
  checkFolders(exp_folder);

  traditional_search(data,exp_folder);
}

void Experiments::traditional_search(
  std::vector<Vector> &random_sample,
  string &exp_folder
) {
  string est_type_file,est_type_folder,results_file;
/*
  ESTIMATION = MOMENT; // Moment estimation ...
  est_type_file = "moment";
  est_type_folder = est_type_file + "/";
  results_file = exp_folder + "traditional_search/results/aic/" + est_type_file;
  ofstream aic_out_mom(results_file.c_str(),ios::app);
  results_file = exp_folder + "traditional_search/results/bic/" + est_type_file;
  ofstream bic_out_mom(results_file.c_str(),ios::app);
  results_file = exp_folder + "traditional_search/results/icl/" + est_type_file;
  ofstream icl_out_mom(results_file.c_str(),ios::app);
  traditional_search(
    random_sample,
    exp_folder,est_type_folder,
    aic_out_mom,bic_out_mom,icl_out_mom
  );
  aic_out_mom.close(); bic_out_mom.close(); icl_out_mom.close();
*/
/*
  ESTIMATION = MLE; // Max LH estimation ...
  est_type_file = "mle";
  est_type_folder = est_type_file + "/";
  results_file = exp_folder + "traditional_search/results/aic/" + est_type_file;
  ofstream aic_out_mle(results_file.c_str(),ios::app);
  results_file = exp_folder + "traditional_search/results/bic/" + est_type_file;
  ofstream bic_out_mle(results_file.c_str(),ios::app);
  results_file = exp_folder + "traditional_search/results/icl/" + est_type_file;
  ofstream icl_out_mle(results_file.c_str(),ios::app);
  traditional_search(
    random_sample,
    exp_folder,est_type_folder,
    aic_out_mle,bic_out_mle,icl_out_mle
  );
  aic_out_mle.close(); bic_out_mle.close(); icl_out_mle.close();
*/

/*
  ESTIMATION = MAP; // MAP estimation ...
  est_type_file = "map";
  est_type_folder = est_type_file + "/";
  results_file = exp_folder + "traditional_search/results/aic/" + est_type_file;
  ofstream aic_out_map(results_file.c_str(),ios::app);
  results_file = exp_folder + "traditional_search/results/bic/" + est_type_file;
  ofstream bic_out_map(results_file.c_str(),ios::app);
  results_file = exp_folder + "traditional_search/results/icl/" + est_type_file;
  ofstream icl_out_map(results_file.c_str(),ios::app);
  traditional_search(
    random_sample,
    exp_folder,est_type_folder,
    aic_out_map,bic_out_map,icl_out_map
  );
  aic_out_map.close(); bic_out_map.close(); icl_out_map.close();
*/


  ESTIMATION = MML; // MML estimation ...
  est_type_file = "mml";
  est_type_folder = est_type_file + "/";
  results_file = exp_folder + "traditional_search/results/" + est_type_file;
  ofstream mml_out(results_file.c_str(),ios::app);
  traditional_search_mml(
    random_sample,
    exp_folder,est_type_folder,
    mml_out
  );
  mml_out.close();

}

void Experiments::traditional_search(
  std::vector<Vector> &random_sample,
  string &exp_folder,
  string &est_type_folder,
  ostream &aic_out,
  ostream &bic_out,
  ostream &icl_out
) {
  string mixtures_folder = exp_folder + "traditional_search/mixtures/";

  string all_inferred_mix = mixtures_folder + "aic/" + est_type_folder + "aic_log";
  ofstream aic_log(all_inferred_mix.c_str());
  all_inferred_mix = mixtures_folder + "bic/" + est_type_folder + "bic_log";
  ofstream bic_log(all_inferred_mix.c_str());
  all_inferred_mix = mixtures_folder + "icl/" + est_type_folder + "icl_log";
  ofstream icl_log(all_inferred_mix.c_str());

  int K_MAX = 50; 

  Vector data_weights(random_sample.size(),1);
  Mixture mixture(1,random_sample,data_weights);
  mixture.estimateParameters();
  double msg = mixture.computeMinimumMessageLength();
  double aic_best = mixture.computeAIC();
  double bic_best = mixture.computeBIC();
  double icl_best = mixture.computeICL();
  aic_log << fixed << setw(10) << 1 << "\t";
  aic_log << fixed << scientific << aic_best << "\t";
  aic_log << fixed << scientific << bic_best << "\t";
  aic_log << fixed << scientific << icl_best << "\t";
  aic_log << fixed << scientific << mixture.first_part() << "\t";
  aic_log << fixed << scientific << mixture.second_part() << "\t";
  aic_log << fixed << scientific << msg << endl;

  bic_log << fixed << setw(10) << 1 << "\t";
  bic_log << fixed << scientific << aic_best << "\t";
  bic_log << fixed << scientific << bic_best << "\t";
  bic_log << fixed << scientific << mixture.first_part() << "\t";
  bic_log << fixed << scientific << mixture.second_part() << "\t";
  bic_log << fixed << scientific << msg << endl;

  icl_log << fixed << setw(10) << 1 << "\t";
  icl_log << fixed << scientific << aic_best << "\t\t";
  icl_log << fixed << scientific << bic_best << "\t\t";
  icl_log << fixed << scientific << icl_best << "\t\t";
  icl_log << fixed << scientific << mixture.first_part() << "\t";
  icl_log << fixed << scientific << mixture.second_part() << "\t";
  icl_log << fixed << scientific << msg << endl;

  Mixture aic_best_mix = mixture;
  Mixture bic_best_mix = mixture;
  Mixture icl_best_mix = mixture;

  int NUM_TRIALS = 5;
  for (int k=2; k<=K_MAX; k++) {
    cout << "k: " << k << endl;
    Mixture aic_best_mix_iter,bic_best_mix_iter,icl_best_mix_iter;
    double aic,bic,icl,aic_best_iter,bic_best_iter,icl_best_iter;
    for (int j=0; j<NUM_TRIALS; j++) {
      Mixture mixture(k,random_sample,data_weights);
      mixture.estimateParameters();
      aic = mixture.computeAIC();
      bic = mixture.computeBIC();
      icl = mixture.computeICL();
      if (j == 0) {
        aic_best_iter = aic;
        aic_best_mix_iter = mixture;
        bic_best_iter = bic;
        bic_best_mix_iter = mixture;
        icl_best_iter = icl;
        icl_best_mix_iter = mixture;
      } else {
        if (aic < aic_best_iter) {
          aic_best_iter = aic;
          aic_best_mix_iter = mixture;
        }
        if (bic < bic_best_iter) {
          bic_best_iter = bic;
          bic_best_mix_iter = mixture;
        }
        if (icl < icl_best_iter) {
          icl_best_iter = icl;
          icl_best_mix_iter = mixture;
        }
      } // if (j == 0)
    } // for (j)

    aic_log << fixed << setw(10) << k << "\t";
    aic_log << fixed << scientific << aic_best_iter << "\t";
    aic_log << fixed << scientific << aic_best_mix_iter.computeBIC() << "\t";
    aic_log << fixed << scientific << aic_best_mix_iter.computeICL() << "\t";
    msg = aic_best_mix_iter.computeMinimumMessageLength();
    aic_log << fixed << scientific << aic_best_mix_iter.first_part() << "\t";
    aic_log << fixed << scientific << aic_best_mix_iter.second_part() << "\t";
    aic_log << fixed << scientific << msg << endl;
    if (aic_best_iter < aic_best) {
      aic_best = aic_best_iter;
      aic_best_mix = aic_best_mix_iter;
    }
    bic_log << fixed << setw(10) << k << "\t";
    bic_log << fixed << scientific << bic_best_mix_iter.computeAIC() << "\t";
    bic_log << fixed << scientific << bic_best_iter << "\t";
    bic_log << fixed << scientific << bic_best_mix_iter.computeICL() << "\t";
    msg = bic_best_mix_iter.computeMinimumMessageLength();
    bic_log << fixed << scientific << bic_best_mix_iter.first_part() << "\t";
    bic_log << fixed << scientific << bic_best_mix_iter.second_part() << "\t";
    bic_log << fixed << scientific << msg << endl;
    if (bic_best_iter < bic_best) {
      bic_best = bic_best_iter;
      bic_best_mix = bic_best_mix_iter;
    }
    icl_log << fixed << setw(10) << k << "\t";
    icl_log << fixed << scientific << icl_best_mix_iter.computeAIC() << "\t";
    icl_log << fixed << scientific << icl_best_mix_iter.computeBIC() << "\t";
    icl_log << fixed << scientific << icl_best_iter << "\t";
    msg = icl_best_mix_iter.computeMinimumMessageLength();
    icl_log << fixed << scientific << icl_best_mix_iter.first_part() << "\t";
    icl_log << fixed << scientific << icl_best_mix_iter.second_part() << "\t";
    icl_log << fixed << scientific << msg << endl;
    if (icl_best_iter < icl_best) {
      icl_best = icl_best_iter;
      icl_best_mix = icl_best_mix_iter;
    }
  } // for (k)
  aic_log.close();
  bic_log.close();
  icl_log.close();

  /* save inferred mixtures */
  string inferred_mix_file = mixtures_folder + "aic/" + est_type_folder + "mixture";
  aic_best_mix.printParameters(inferred_mix_file);
  aic_out << fixed << setw(10) << aic_best_mix.getNumberOfComponents() << "\t\t";
  aic_out << fixed << scientific << aic_best << "\t\t";
  aic_out << fixed << scientific << aic_best_mix.computeMinimumMessageLength() << endl;

  inferred_mix_file = mixtures_folder + "bic/" + est_type_folder + "mixture";  
  bic_best_mix.printParameters(inferred_mix_file);
  bic_out << fixed << setw(10) << bic_best_mix.getNumberOfComponents() << "\t\t";
  bic_out << fixed << scientific << bic_best << "\t\t";
  bic_out << fixed << scientific << bic_best_mix.computeMinimumMessageLength() << endl;

  inferred_mix_file = mixtures_folder + "icl/" + est_type_folder + "mixture";  
  icl_best_mix.printParameters(inferred_mix_file);
  icl_out << fixed << setw(10) << icl_best_mix.getNumberOfComponents() << "\t\t";
  icl_out << fixed << scientific << icl_best << "\t\t";
  icl_out << fixed << scientific << icl_best_mix.computeMinimumMessageLength() << endl;
}

void Experiments::traditional_search_mml(
  std::vector<Vector> &random_sample,
  string &exp_folder,
  string &est_type_folder,
  ostream &mml_out
) {
  string mixtures_folder = exp_folder + "traditional_search/mixtures/";

  string all_inferred_mix = mixtures_folder + est_type_folder + "mml_log";
  ofstream mml_log(all_inferred_mix.c_str());

  int K_MAX = 50; 

  Vector data_weights(random_sample.size(),1);
  Mixture mixture(1,random_sample,data_weights);
  mixture.estimateParameters();

  double msglen_best = mixture.getMinimumMessageLength();
  mml_log << fixed << setw(10) << 1 << "\t";
  mml_log << fixed << scientific << mixture.computeAIC() << "\t";
  mml_log << fixed << scientific << mixture.computeBIC() << "\t";
  mml_log << fixed << scientific << mixture.computeICL() << "\t";
  mml_log << fixed << scientific << mixture.first_part() << "\t";
  mml_log << fixed << scientific << mixture.second_part() << "\t";
  mml_log << fixed << scientific << msglen_best << endl;
  Mixture mml_best_mix = mixture;

  int NUM_TRIALS = 5;
  for (int k=38; k<=K_MAX; k++) {
    cout << "k: " << k << endl;
    Mixture mml_best_mix_iter;
    double msglen,msglen_best_iter;
    for (int j=0; j<NUM_TRIALS; j++) {
      cout << "\ttrial #" << j+1 << endl;
      Mixture mixture(k,random_sample,data_weights);
      mixture.estimateParameters();
      msglen = mixture.getMinimumMessageLength();
      if (j == 0) {
        msglen_best_iter = msglen;
        mml_best_mix_iter = mixture;
      } else {
        if (msglen < msglen_best_iter) {
          msglen_best_iter = msglen;
          mml_best_mix_iter = mixture;
        }
      } // if (j == 0)
    } // for (j)
    mml_log << fixed << setw(10) << k << "\t\t";
    mml_log << fixed << scientific << mml_best_mix_iter.computeAIC() << "\t";
    mml_log << fixed << scientific << mml_best_mix_iter.computeBIC() << "\t";
    mml_log << fixed << scientific << mml_best_mix_iter.computeICL() << "\t";
    mml_log << fixed << scientific << mml_best_mix_iter.first_part() << "\t";
    mml_log << fixed << scientific << mml_best_mix_iter.second_part() << "\t";
    mml_log << fixed << scientific << msglen_best_iter << endl;
    if (msglen_best_iter < msglen_best) {
      msglen_best = msglen_best_iter;
      mml_best_mix = mml_best_mix_iter;
    }
  } // for()
  mml_log.close();

  /* save inferred mixtures */
  string inferred_mix_file = mixtures_folder + est_type_folder + "mixture";
  mml_best_mix.printParameters(inferred_mix_file);
  mml_out << fixed << setw(10) << mml_best_mix.getNumberOfComponents() << "\t\t";
  mml_out << fixed << scientific << msglen_best << "\t\t";
  mml_out << fixed << scientific << mml_best_mix.getMinimumMessageLength() << endl;
}

void Experiments::exp4()
{
  //string exp_folder = "./experiments/infer_components/exp4/";
  //check_and_create_directory(exp_folder);
  //string exp_folder = "./";
  //checkFolders(exp_folder);

  int N = 10000;
  std::vector<Vector> sampled_data = generate_data_exp4(N);
  //string data_file = "constrained_N_10000.dat";
  string data_file = "./data/arun_N_10000.dat";
  writeToFile(data_file,sampled_data);

  /*for (int N=1000; N<=10000; N+=1000) {
    generate_data_exp4(exp_folder,N);
  }*/

//  int N = 1000;
//  std::vector<Vector> data = generate_data_exp4(N);
//  std::vector<Vector> combined_data = data;
//  double res = 1;
//  do {
//    cout << "N: " << N << endl;
//    if (N > 1000) {
//      std::vector<Vector> data2 = generate_data_exp4(1000);
//      /* merge data */
//      for (int i=0; i<data2.size(); i++) {
//        combined_data.push_back(data2[i]);
//      } // for (i)
//    } // if (N2 != 0)
//
//    /* save data */
//    std::vector<std::vector<int> > sampled_bins = updateBins(combined_data,res);
//    string N_str = boost::lexical_cast<string>(N);
//    string data_folder = exp_folder + "data/";
//    check_and_create_directory(data_folder);
//    string sampled_bins_file = data_folder + "N_" + N_str + "_sampled_bins.dat";
//    writeToFile(sampled_bins_file,sampled_bins);
//    string data_file = data_folder + "N_" + N_str + ".dat";
//    writeToFile(data_file,combined_data);
//
//    N += 1000;
//  } while (N <= 50000);

  //infer_components_exp4(exp_folder,N);
}

std::vector<Vector> Experiments::generate_data_exp4(int N)
{
  struct Parameters parameters;
  //parameters.profiles_dir = "./data/profiles-b/";
  //parameters.profiles_dir = "./data/profiles/";
  //parameters.profile_file = "constrained_bins.dat";
  parameters.profile_file = "./data/arun.dat";

  std::vector<Vector> data;
  gatherData(parameters,data);
  cout << "data.size(): " << data.size() << endl;
  double res = 1;
  std::vector<std::vector<int> > true_bins = updateBins(data,res);
  //string true_bins_file = "true_bins.dat";  // integers
  //writeToFile(true_bins_file,true_bins);

  int num_rows = true_bins.size();
  int num_cols = true_bins[0].size();
  int num_bins = num_rows * num_cols;
  cout << "num_bins: " << num_bins << endl;

  Vector emptyvec(num_cols,0);
  std::vector<Vector> prob_bins(num_rows,emptyvec);
  Vector elements(num_bins,0);
  int count = 0;
  for (int i=0; i<num_rows; i++) {
    for (int j=0; j<num_cols; j++) {
      prob_bins[i][j] = true_bins[i][j] / (double) data.size();
      elements[count++] = prob_bins[i][j];
    } // for (j)
  } // for (i)
  //string prob_bins_file = "prob_bins.dat";  // fractional values
  //writeToFile(prob_bins_file,prob_bins);

  std::vector<int> sorted_index;
  Vector sorted_elements = sort(elements,sorted_index);
  Vector cumsum(num_bins,0);
  cumsum[0] = sorted_elements[0];
  for (int i=1; i<num_bins; i++) {
    cumsum[i] = cumsum[i-1] + sorted_elements[i];
    //cout << sorted_index[i] << "\t\t" << cumsum[i] << endl;
  }

  std::vector<Vector> random_sample;
  for (int i=0; i<N; i++) {
    double cdf = uniform_random();
    int bin;
    for (int j=0; j<num_bins; j++) {
      if (cdf <= cumsum[j]) {
        bin = sorted_index[j];
        break;
      } // if ()
    } // for (j)
    int row = bin / num_cols;
    double theta = (row + uniform_random()) * res;  // in degrees
    //double theta = row * res;
    int col = bin % num_cols;
    double phi = (col + uniform_random()) * res;   // in degrees`
    //double phi = col * res;
    Vector spherical(3,0),cartesian(3,0);
    spherical[0] = 1;
    spherical[1] = theta * PI/180;
    spherical[2] = phi * PI/180;
    spherical2cartesian(spherical,cartesian);
    random_sample.push_back(cartesian);
  } // for (i)

  return random_sample;
}

void Experiments::infer_components_exp4(string &exp_folder, int N)
{
  string N_str = boost::lexical_cast<string>(N);
  string data_file = exp_folder + "data/N_" + N_str + ".dat";
  std::vector<Vector> data = load_data_table(data_file);

  traditional_search(data,exp_folder);
}

// mix_example
void Experiments::exp5()
{
/*
  double psi,alpha,eta,kappa,beta,ecc;
  std::vector<Kent> components;

  // (1)
  psi = 0; alpha = 60; eta = 45; 
  psi *= PI/180; alpha *= PI/180; eta *= PI/180;
  kappa = 100; ecc = 0.1;
  beta = 0.5 * kappa * ecc;
  Kent kent1(psi,alpha,eta,kappa,beta);
  components.push_back(kent1);

  // (2)
  psi = 150; alpha = 45; eta = 30; 
  psi *= PI/180; alpha *= PI/180; eta *= PI/180;
  kappa = 100; ecc = 0.5;
  beta = 0.5 * kappa * ecc;
  Kent kent2(psi,alpha,eta,kappa,beta);
  components.push_back(kent2);

  // (3)
  psi = 30; alpha = 45; eta = 60; 
  psi *= PI/180; alpha *= PI/180; eta *= PI/180;
  kappa = 100; ecc = 0.9;
  beta = 0.5 * kappa * ecc;
  Kent kent3(psi,alpha,eta,kappa,beta);
  components.push_back(kent3);

  double w = 1.0/3.0;
  Vector weights(3,w);
  int N = 1000;

  Mixture mix(3,components,weights);
  ofstream file("./simulation/simulated_mixture");
  for (int i=0; i<3; i++) {
    file << fixed << setw(10) << setprecision(5) << weights[i];
    file << "\t";
    components[i].printParameters(file);
  }
  file.close();
  mix.generate(N,1);
*/

  // mml infer components (traditional search)
  string exp_folder = "./experiments/infer_components/exp5/";
  string data_file = exp_folder + "mix_example/random_sample_mix.dat";
  std::vector<Vector> data = load_data_table(data_file);

  string all_inferred_mix = exp_folder + "mml_log";
  ofstream mml_log(all_inferred_mix.c_str());

  int K_MAX = 15; 

  Vector data_weights(data.size(),1);
  Mixture mixture(1,data,data_weights);
  mixture.estimateParameters();

  double msglen_best = mixture.getMinimumMessageLength();
  mml_log << fixed << setw(10) << 1 << "\t\t";
  mml_log << fixed << scientific << mixture.first_part() << "\t\t";
  mml_log << fixed << scientific << mixture.second_part() << "\t\t";
  mml_log << fixed << scientific << msglen_best << endl;
  Mixture mml_best_mix = mixture;

  for (int k=2; k<=K_MAX; k++) {
    cout << "k: " << k << endl;
    Mixture mixture(k,data,data_weights);
    mixture.estimateParameters();
    double msglen = mixture.getMinimumMessageLength();
    mml_log << fixed << setw(10) << k << "\t\t";
    mml_log << fixed << scientific << mixture.first_part() << "\t\t";
    mml_log << fixed << scientific << mixture.second_part() << "\t\t";
    mml_log << fixed << scientific << msglen << endl;
    if (msglen < msglen_best) {
      msglen_best = msglen;
      mml_best_mix = mixture;
    }
  } // for()
  mml_log.close();

  /* save inferred mixtures */
  string inferred_mix_file = exp_folder + "best_mml_mixture";
  mml_best_mix.printParameters(inferred_mix_file);
}

void Experiments::exp6()
{
  struct Parameters parameters;
  parameters.profiles_dir = "./data/profiles/";
  //parameters.profile_file = "random_sample_mix.dat";

  ENABLE_DATA_PARALLELISM = SET;
  NUM_THREADS = 4;
  std::vector<Vector> data;
  gatherData(parameters,data);
  cout << "data.size(): " << data.size() << endl;

  string log_file = "infer_vmf_all_jk.log";
  inferComponents_vMF(data,log_file);
}


