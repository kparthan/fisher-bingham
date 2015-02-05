#include "Support.h"
#include "Experiments.h"
#include "Kent.h"

extern Vector XAXIS,YAXIS,ZAXIS;
extern int MOMENT_FAIL,MLE_FAIL,MAP_FAIL,MML_FAIL,ESTIMATION;
extern struct stat st;

Experiments::Experiments(int iterations) : iterations(iterations)
{}

void Experiments::simulate()
{
  int N = 100;
  double kappa,beta,eccentricity;

  //string n_str = "N_" + boost::lexical_cast<string>(N) + "_uniform_prior/";
  //string n_str = "N_" + boost::lexical_cast<string>(N) + "_vmf_prior/";
  //string n_str = "N_" + boost::lexical_cast<string>(N) + "_beta_prior/";
  string n_str = "N_" + boost::lexical_cast<string>(N) + "_new_prior/";
  string parent_dir = "experiments/single_kent/" + n_str + "/";
  string current_dir,kappa_str,eccentricity_str;
  string kappas,betas,negloglike,kldivs,msglens;
  std::vector<Vector> random_sample;
  std::vector<struct Estimates> all_estimates;
  //string data_file = "random_sample_uniform.dat";
  //string data_file = "random_sample_vmf.dat";
  //string data_file = "random_sample_beta.dat";
  //string data_file = "random_sample.dat";
  string data_file = "random_sample_new.dat";

  double INIT_KAPPA = 10;
  double MAX_KAPPA = 100;
  double KAPPA_INCREMENT = 10;

  kappa = INIT_KAPPA;
  while (kappa <= MAX_KAPPA) {
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
        random_sample = load_data_table(data_file);
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
    kappa += KAPPA_INCREMENT;
  } // kappa
}

void Experiments::checkFolders(string &exp_folder)
{
  string data_folder = exp_folder + "data/";
  if (stat(data_folder.c_str(), &st) == -1) {
    cout << "Error: data folder not present ...\n";
    exit(1);
  }

  string logs_folder = exp_folder + "logs/";
  if (stat(logs_folder.c_str(), &st) == -1) {
    mkdir(logs_folder.c_str(), 0700);
  }
  string tmp = logs_folder + "moment/";
  if (stat(tmp.c_str(), &st) == -1) {
    mkdir(tmp.c_str(), 0700);
  }
  tmp = logs_folder + "mle/";
  if (stat(tmp.c_str(), &st) == -1) {
    mkdir(tmp.c_str(), 0700);
  }
  tmp = logs_folder + "map/";
  if (stat(tmp.c_str(), &st) == -1) {
    mkdir(tmp.c_str(), 0700);
  }
  tmp = logs_folder + "mml/";
  if (stat(tmp.c_str(), &st) == -1) {
    mkdir(tmp.c_str(), 0700);
  }

  string mixtures_folder = exp_folder + "mixtures/";
  if (stat(mixtures_folder.c_str(), &st) == -1) {
    mkdir(mixtures_folder.c_str(), 0700);
  }
  tmp = mixtures_folder + "moment/";
  if (stat(tmp.c_str(), &st) == -1) {
    mkdir(tmp.c_str(), 0700);
  }
  tmp = mixtures_folder + "mle/";
  if (stat(tmp.c_str(), &st) == -1) {
    mkdir(tmp.c_str(), 0700);
  }
  tmp = mixtures_folder + "map/";
  if (stat(tmp.c_str(), &st) == -1) {
    mkdir(tmp.c_str(), 0700);
  }
  tmp = mixtures_folder + "mml/";
  if (stat(tmp.c_str(), &st) == -1) {
    mkdir(tmp.c_str(), 0700);
  }

  string results_folder = exp_folder + "results/";
  if (stat(results_folder.c_str(), &st) == -1) {
    mkdir(results_folder.c_str(), 0700);
  }
}

void Experiments::infer_components_exp1()
{
  iterations = 50;
  int N = 100;

  string parent_folder = "./experiments/infer_components/exp1/";
  string exp_folder;

  double alpha,alpha_rad,sin_alpha,cos_alpha;
  double kappa,beta,ecc;
  Vector mean,major,minor;

  //for (double alpha=5; alpha<=25; alpha+=5) {
    alpha = 10; // in degrees
    alpha_rad = alpha * PI / 180;

    mean = ZAXIS;
    major = XAXIS;
    minor = YAXIS;
    kappa = 10; ecc = 0.1;
    beta = 0.5 * kappa * ecc;
    Kent kent1(mean,major,minor,kappa,beta);

    sin_alpha = sin(alpha_rad);
    cos_alpha = cos(alpha_rad);
    mean[0] = sin_alpha; mean[1] = 0; mean[2] = cos_alpha;
    major[0] = cos_alpha; major[1] = 0; major[2] = -sin_alpha;
    kappa = 10; ecc = 0.1;
    beta = 0.5 * kappa * ecc;
    Kent kent2(mean,major,minor,kappa,beta);

    Vector weights(2,0.5);

    std::vector<Kent> components;
    components.push_back(kent1);
    components.push_back(kent2);
    Mixture original(2,components,weights);

    std::ostringstream ss;
    ss << fixed << setprecision(0);
    ss << alpha;
    string alpha_str = "alpha_" + ss.str();
    exp_folder = parent_folder + alpha_str + "/" ;
    if (stat(exp_folder.c_str(), &st) == -1) {
      mkdir(exp_folder.c_str(), 0700);
    }

    //generateData(original,exp_folder,N);
    inferMixtures(original,exp_folder);
  //}
}

void Experiments::generateData(
  Mixture &original, 
  string &exp_folder,
  int sample_size
) {
  string data_folder = exp_folder + "data/";
  if (stat(data_folder.c_str(), &st) == -1) {
      mkdir(data_folder.c_str(), 0700);
  }
  
  string iter_str,data_file;
  std::vector<Vector> data;

  for (int iter=1; iter<=iterations; iter++) {
    iter_str = boost::lexical_cast<string>(iter);
    data_file =  data_folder + "mixture_iter_" + iter_str + ".dat";
    data = original.generate(sample_size,0);
    writeToFile(data_file,data);
  }
}

void Experiments::inferMixtures(
  Mixture &original, 
  string &exp_folder
) {
  checkFolders(exp_folder);

  string data_folder = exp_folder + "data/";
  string logs_folder = exp_folder + "logs/";
  string mixtures_folder = exp_folder + "mixtures/";
  string results_folder = exp_folder + "results/";

  string data_file,log_file,mixture_file,iter_str;
  std::vector<Vector> data;
  Mixture inferred;
  
  int num_success = 0;
  double avg_number,variance;

  for (int iter=1; iter<=iterations; iter++) {
    iter_str = boost::lexical_cast<string>(iter);

    data_file = data_folder + "mixture_iter_" + iter_str + ".dat";
    cout << "data_file: " << data_file << endl;
    data = load_data_table(data_file);

    for (int j=0; j<NUM_METHODS; j++) {
      ESTIMATION = j;
      switch(ESTIMATION) {
        case MOMENT:
          log_file = logs_folder + "moment/mixture_iter_" + iter_str + ".log";
          mixture_file = mixtures_folder + "moment/mixture_iter_" + iter_str;
          break;

        case MLE:
          log_file = logs_folder + "mle/mixture_iter_" + iter_str + ".log";
          mixture_file = mixtures_folder + "mle/mixture_iter_" + iter_str;
          break;

        case MAP:
          log_file = logs_folder + "map/mixture_iter_" + iter_str + ".log";
          mixture_file = mixtures_folder + "map/mixture_iter_" + iter_str;
          break;

        case MML:
          log_file = logs_folder + "mml/mixture_iter_" + iter_str + ".log";
          mixture_file = mixtures_folder + "mml/mixture_iter_" + iter_str;
          break;

      } // switch() ends ...
      if (ESTIMATION != MML) {
        inferred = inferComponents_ML(data,log_file);
      } else {
        inferred = inferComponents(data,log_file);
      }
      inferred.printParameters(mixture_file);
    }
  } // iter() loop ...
  /*avg_number = computeMean(inferred);
  variance = computeVariance(inferred);
  summary << "\nsuccess rate: " //<< setprecision(2) 
          << num_success * 100 / (double)(iterations) << " %\n";
  summary << "Avg: " << avg_number << endl;
  summary << "Variance: " << variance << endl;
  summary.close();*/
}

