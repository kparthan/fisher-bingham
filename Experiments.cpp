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

Experiments::Experiments(int iterations) : iterations(iterations)
{}

void Experiments::simulate()
{
  int N = 45;
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
  string psi_est,alpha_est,eta_est,kappas,betas,ecc_file,negloglike,kldivs,msglens;
  std::vector<Vector> random_sample;
  std::vector<struct Estimates> all_estimates;
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

      ofstream fpsi(psi_est.c_str(),ios::app);
      ofstream falpha(alpha_est.c_str(),ios::app);
      ofstream feta(eta_est.c_str(),ios::app);
      ofstream fk(kappas.c_str(),ios::app);
      ofstream fb(betas.c_str(),ios::app);
      ofstream fecc(ecc_file.c_str(),ios::app);
      ofstream fnlh(negloglike.c_str(),ios::app);
      ofstream fkl(kldivs.c_str(),ios::app);
      ofstream fmsg(msglens.c_str(),ios::app);
      
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
        } // for j ()
        fpsi << endl; falpha << endl; feta << endl;
        fecc << endl; fk << endl; fb << endl;
        fnlh << endl;
        fkl << endl;
        fmsg << endl;
      } // for i ()

      eccentricity += 0.4;
      fecc.close();
      fk.close();
      fb.close();
      fnlh.close();
      fkl.close();
      fmsg.close();
    } // eccentricity
    //kappa += KAPPA_INCREMENT;
    kappa *= KAPPA_INCREMENT;
  } // kappa
}

void Experiments::checkFolders(string &exp_folder)
{
  string data_folder = exp_folder + "data/";
  check_and_create_directory(data_folder);

  string criterion;
  string logs_folder = exp_folder + "logs/";
  check_and_create_directory(logs_folder);
  criterion = "bic/";
  create_sub_folders(logs_folder,criterion);
  criterion = "icl/";
  create_sub_folders(logs_folder,criterion);

  string mixtures_folder = exp_folder + "mixtures/";
  check_and_create_directory(mixtures_folder);
  criterion = "bic/";
  create_sub_folders(mixtures_folder,criterion);
  criterion = "icl/";
  create_sub_folders(mixtures_folder,criterion);

  string results_folder = exp_folder + "results/";
  check_and_create_directory(results_folder);
  criterion = "bic/";
  string results = results_folder + criterion;
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

void Experiments::infer_components_exp1()
{
  INFER_COMPONENTS = SET;

  double kappa,ecc;

  kappa = 100;
  //while(kappa <= 101) {
    ecc = 0.6;
    while(ecc < 0.95) {
      infer_components_exp1(kappa,ecc);
      ecc += 0.1;
    } // ecc()
    //kappa += 10;
  //} // kappa()
}

void Experiments::infer_components_exp1(double kappa, double ecc)
{
  iterations = 50;
  int N = 100;

  string common = "./experiments/infer_components/exp1/";
  string parent_folder,exp_folder;

  double alpha,alpha_rad,sin_alpha,cos_alpha;
  double beta;
  Vector mean,major,minor;

  for (alpha=10; alpha<=50; alpha+=10) {
    //alpha = 50; // in degrees
    alpha_rad = alpha * PI / 180;

    mean = ZAXIS;
    major = XAXIS;
    minor = YAXIS;
    //kappa = 10; ecc = 0.9;
    beta = 0.5 * kappa * ecc;
    Kent kent1(mean,major,minor,kappa,beta);

    sin_alpha = sin(alpha_rad);
    cos_alpha = cos(alpha_rad);
    mean[0] = sin_alpha; mean[1] = 0; mean[2] = cos_alpha;
    major[0] = cos_alpha; major[1] = 0; major[2] = -sin_alpha;
    //kappa = 10; ecc = 0.9;
    beta = 0.5 * kappa * ecc;
    Kent kent2(mean,major,minor,kappa,beta);

    Vector weights(2,0.5);

    std::vector<Kent> components;
    components.push_back(kent1);
    components.push_back(kent2);
    Mixture original(2,components,weights);

    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();

    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << ecc;
    string eccentricity_str = sse.str();

    std::ostringstream ss;
    ss << fixed << setprecision(0);
    ss << alpha;
    string alpha_str = ss.str();

    parent_folder = common + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
    check_and_create_directory(parent_folder);
    exp_folder = parent_folder + "alpha_" + alpha_str + "/" ;
    check_and_create_directory(exp_folder);

    generateData(original,exp_folder,N);
    inferMixtures_exp1(original,exp_folder);
  } // alpha()
}

void Experiments::generateData(
  Mixture &original, 
  string &exp_folder,
  int sample_size
) {
  string data_folder = exp_folder + "data/";
  check_and_create_directory(data_folder);

  string iter_str,data_file;
  std::vector<Vector> data;
  for (int iter=1; iter<=iterations; iter++) {
    iter_str = boost::lexical_cast<string>(iter);
    data_file =  data_folder + "mixture_iter_" + iter_str + ".dat";
    data = original.generate(sample_size,0);
    writeToFile(data_file,data);
  }

  // generate large data for computing empirical KL-divergence later ...  
  int num_samples = 100000;
  std::vector<Vector> random_sample = original.generate(num_samples,0);
  data_file = data_folder + "random_sample.dat";
  writeToFile(data_file,random_sample);
}

void Experiments::inferMixtures_exp1(
  Mixture &original, 
  string &exp_folder
) {
  checkFolders(exp_folder);

  string data_folder = exp_folder + "data/";
  string logs_folder = exp_folder + "logs/";
  string mixtures_folder = exp_folder + "mixtures/";
  string results_folder = exp_folder + "results/";
  string criterion;

  string big_data_file = data_folder + "random_sample.dat";
  std::vector<Vector> large_data = load_data_table(big_data_file);

  criterion = "bic/";
  CRITERION = BIC;

  ESTIMATION = MOMENT;
  inferMixtures_exp1(
    original,large_data,data_folder,logs_folder,mixtures_folder,results_folder,criterion
  );

  ESTIMATION = MLE;
  inferMixtures_exp1(
    original,large_data,data_folder,logs_folder,mixtures_folder,results_folder,criterion
  );

  ESTIMATION = MAP;
  inferMixtures_exp1(
    original,large_data,data_folder,logs_folder,mixtures_folder,results_folder,criterion
  );

  /*criterion = "icl/";
  CRITERION = ICL;

  ESTIMATION = MOMENT;
  inferMixtures_exp1(
    original,large_data,data_folder,logs_folder,mixtures_folder,results_folder,criterion
  );

  ESTIMATION = MLE;
  inferMixtures_exp1(
    original,large_data,data_folder,logs_folder,mixtures_folder,results_folder,criterion
  );

  ESTIMATION = MAP;
  inferMixtures_exp1(
    original,large_data,data_folder,logs_folder,mixtures_folder,results_folder,criterion
  );*/

  criterion = "mml/";
  CRITERION = MMLC;
  ESTIMATION = MML;
  inferMixtures_exp1(
    original,large_data,data_folder,logs_folder,mixtures_folder,results_folder,criterion
  );
}

void Experiments::inferMixtures_exp1(
  Mixture &original, 
  std::vector<Vector> &large_data,
  string &data_folder,
  string &logs_folder,
  string &mixtures_folder,
  string &results_folder,
  string &criterion
) {
  string data_file,log_file,mixture_file,results_file,iter_str,big_data_file;
  string logs_folder_specific,mixtures_folder_specific;
  std::vector<Vector> data;
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
  //big_data_file = data_folder + "random_sample.dat";
  //std::vector<Vector> large_data = load_data_table(big_data_file);

  ofstream out(results_file.c_str());
  for (int iter=1; iter<=iterations; iter++) {
    iter_str = boost::lexical_cast<string>(iter);

    data_file = data_folder + "mixture_iter_" + iter_str + ".dat";
    cout << "data_file: " << data_file << endl;
    data = load_data_table(data_file);

    log_file = logs_folder_specific + "mixture_iter_" + iter_str + ".log";
    mixture_file = mixtures_folder_specific + "mixture_iter_" + iter_str;

    inferred = inferComponents(data,log_file);
    inferred.printParameters(mixture_file);
    kldiv = original.computeKLDivergence(inferred,large_data);

    out << setw(10) << inferred.getNumberOfComponents() << "\t\t";
    out << scientific << kldiv << "\t\t";
    out << scientific << inferred.getMinimumMessageLength() << endl;
  } // iter() loop ...
  out.close();
}

void Experiments::infer_components_exp2()
{
  INFER_COMPONENTS = SET;
  int K = 6;
  int N = 1000;
  int large_N = 100000;
  int num_mixtures = 50;
  iterations = 1;

  string common = "./experiments/infer_components/exp2/";
  string K_str = boost::lexical_cast<string>(K);
  EM_LOG_FOLDER = "./infer/logs/kent/K_" + K_str + "/";
  check_and_create_directory(EM_LOG_FOLDER);

  string exp_folder = common + "K_" + K_str + "/";
  check_and_create_directory(exp_folder);
  checkFolders(exp_folder);

  string data_folder = exp_folder + "data/";
  string original_mix_folder = exp_folder + "original/";
  check_and_create_directory(original_mix_folder);

  string mix_file,index_str,data_file;
  std::vector<Vector> random_sample,large_sample;
  for (int index=1; index<=num_mixtures; index++) {
    index_str = boost::lexical_cast<string>(index);

    Vector weights = generateFromSimplex(K);
    std::vector<Kent> components = generateRandomComponents(K);
    Mixture original(K,components,weights);
    mix_file = original_mix_folder + "mixture_" + index_str;
    original.printParameters(mix_file);

    large_sample = original.generate(large_N,0);
    random_sample = original.generate(N,0);
    data_file = data_folder + "mixture_" + index_str + ".dat";
    writeToFile(data_file,random_sample);

    inferMixtures_exp2(index_str,original,random_sample,large_sample,exp_folder);
  } // for()
}

void Experiments::inferMixtures_exp2(
  string &index_str,
  Mixture &original, 
  std::vector<Vector> &random_sample,
  std::vector<Vector> &large_sample,
  string &exp_folder
) {
  string logs_folder = exp_folder + "logs/";
  string mixtures_folder = exp_folder + "mixtures/";
  string results_folder = exp_folder + "results/";
  string criterion;

  criterion = "bic/";
  CRITERION = BIC;

  ESTIMATION = MOMENT;
  inferMixtures_exp2(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  ESTIMATION = MLE;
  inferMixtures_exp2(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  ESTIMATION = MAP;
  inferMixtures_exp2(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  criterion = "icl/";
  CRITERION = ICL;

  ESTIMATION = MOMENT;
  inferMixtures_exp2(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  ESTIMATION = MLE;
  inferMixtures_exp2(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  ESTIMATION = MAP;
  inferMixtures_exp2(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );

  criterion = "mml/";
  CRITERION = MMLC;

  ESTIMATION = MML;
  inferMixtures_exp2(
    index_str,
    original,random_sample,large_sample,
    logs_folder,mixtures_folder,results_folder,
    criterion
  );
}

void Experiments::inferMixtures_exp2(
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

void Experiments::infer_components_exp3()
{
  INFER_COMPONENTS = SET;

  double kappa1,kappa2,ecc;

  kappa1 = 10; kappa2 = 100;
  //while(kappa <= 101) {
    ecc = 0.5;
    while(ecc < 0.95) {
      infer_components_exp3(kappa1,kappa2,ecc);
      ecc += 0.1;
    } // ecc()
    //kappa += 10;
  //} // kappa()
}

void Experiments::infer_components_exp3(double kappa1, double kappa2, double ecc)
{
  iterations = 25;
  int N = 100;

  string common = "./experiments/infer_components/exp3/";
  string parent_folder,exp_folder;

  double alpha,alpha_rad,sin_alpha,cos_alpha;
  double beta1,beta2;
  Vector mean,major,minor;

  //for (alpha=10; alpha<=50; alpha+=10) {
    alpha = 50; // in degrees
    alpha_rad = alpha * PI / 180;

    mean = ZAXIS;
    major = XAXIS;
    minor = YAXIS;
    beta1 = 0.5 * kappa1 * ecc;
    Kent kent1(mean,major,minor,kappa1,beta1);

    sin_alpha = sin(alpha_rad);
    cos_alpha = cos(alpha_rad);
    mean[0] = sin_alpha; mean[1] = 0; mean[2] = cos_alpha;
    major[0] = cos_alpha; major[1] = 0; major[2] = -sin_alpha;
    beta2 = 0.5 * kappa2 * ecc;
    Kent kent2(mean,major,minor,kappa2,beta2);

    Vector weights(2,0);
    weights[0] = 0.8; weights[1] = 0.2;

    std::vector<Kent> components;
    components.push_back(kent1);
    components.push_back(kent2);
    Mixture original(2,components,weights);

    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << ecc;
    string eccentricity_str = sse.str();

    std::ostringstream ss;
    ss << fixed << setprecision(0);
    ss << alpha;
    string alpha_str = ss.str();

    parent_folder = common + "e_" + eccentricity_str + "/";
    check_and_create_directory(parent_folder);
    exp_folder = parent_folder + "alpha_" + alpha_str + "/" ;
    check_and_create_directory(exp_folder);

    generateData(original,exp_folder,N);
    inferMixtures_exp1(original,exp_folder);
  //} // alpha()
}

