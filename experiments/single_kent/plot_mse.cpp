#include <iostream>
#include <fstream>
#include <sstream>
#include <vector> 
#include <cstdlib>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost::program_options;
using namespace boost::filesystem;

typedef std::vector<double> Vector;

#define NUM_METHODS 4 
#define MOMENT 0 
#define MLE 1
#define MAP 2
#define MML 3 

std::vector<Vector> load_matrix(string &file_name, int D)
{
  std::vector<Vector> sample;
  ifstream file(file_name.c_str());
  string line;
  Vector numbers(D,0);
  int i;
  while (getline(file,line)) {
    boost::char_separator<char> sep(" \t");
    boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
    i = 0;
    BOOST_FOREACH(const string &t, tokens) {
      istringstream iss(t);
      double x;
      iss >> x;
      numbers[i++] = x;
    }
    sample.push_back(numbers);
  }
  file.close();
  return sample;
}

void computeMeanSquaredError(
  ostream &out, double p, std::vector<Vector> &p_est_all
) {
  int num_elements = p_est_all.size();

  Vector error(p_est_all[0].size(),0);
  for (int i=0; i<num_elements; i++) {
    for (int j=0; j<p_est_all[0].size(); j++) {
      error[j] += (p - p_est_all[i][j]) * (p - p_est_all[i][j]);
    }
  }
  for (int j=0; j<p_est_all[0].size(); j++) {
    error[j] /= num_elements;
    out << fixed << scientific << setprecision(6) << error[j] << "\t";
  }
  out << endl;
}

struct Parameters
{
  int N;
  double kappa,eccentricity;
  int quantity;
};

// all k fixed e, MSE(k)
void tabulate_eccentricity_mse_kappas(
  string &data_file, string &n_str, string &eccentricity_str
) {
  std::vector<Vector> all_kappas;
  string kappas_file;
  ofstream out(data_file.c_str());
  double kappa = 10;
  while (kappa <= 100) {
    out << fixed << setw(10) << setprecision(0) << kappa << "\t";
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
    string kappas_file = current_dir + "kappas";
    all_kappas = load_matrix(kappas_file,NUM_METHODS);
    computeMeanSquaredError(out,kappa,all_kappas);
    kappa += 10;
  } // while()
  out.close();
}

// all k fixed e, MSE(b)
void tabulate_eccentricity_mse_betas(
  string &data_file, string &n_str, double eccentricity, string &eccentricity_str
) {
  std::vector<Vector> all_betas;
  string betas_file;
  ofstream out(data_file.c_str());
  double kappa = 10,beta;
  while (kappa <= 100) {
    beta = 0.5 * kappa * eccentricity;
    out << fixed << setw(10) << setprecision(0) << kappa << "\t";
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
    string betas_file = current_dir + "betas";
    all_betas = load_matrix(betas_file,NUM_METHODS);
    computeMeanSquaredError(out,beta,all_betas);
    kappa += 10;
  } // while()
  out.close();
}

// all e fixed k, MSE(k)
void tabulate_kappa_mse_kappas(
  string &data_file, string &n_str, double kappa, string &kappa_str
) {
  std::vector<Vector> all_kappas;
  string kappas_file;
  ofstream out(data_file.c_str());
  double eccentricity = 0.1;
  while (eccentricity < 0.95) {
    out << fixed << setw(10) << setprecision(1) << eccentricity << "\t";
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << eccentricity;
    string eccentricity_str = sse.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
    kappas_file = current_dir + "kappas";
    all_kappas = load_matrix(kappas_file,NUM_METHODS);
    computeMeanSquaredError(out,kappa,all_kappas);
    eccentricity += 0.1;
  } // while()
  out.close();
}

// all e fixed k, MSE(b)
void tabulate_kappa_mse_betas(
  string &data_file, string &n_str, double kappa, string &kappa_str
) {
  std::vector<Vector> all_betas;
  string betas_file;
  ofstream out(data_file.c_str());
  double eccentricity = 0.1,beta;
  while (eccentricity < 0.95) {
    beta = 0.5 * kappa * eccentricity;
    out << fixed << setw(10) << setprecision(1) << eccentricity << "\t";
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << eccentricity;
    string eccentricity_str = sse.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
    betas_file = current_dir + "betas";
    all_betas = load_matrix(betas_file,NUM_METHODS);
    computeMeanSquaredError(out,beta,all_betas);
    eccentricity += 0.1;
  } // while()
  out.close();
}

void plot_script_mse_kappas(
  string &data_file, string &n_str, string &type_str, string &xlabel, string &ylabel
) {
  string script_file = n_str + type_str + "_mse_kappas.p";
  string plot_file = n_str + type_str + "_mse_kappas.eps";
  ofstream out(script_file.c_str());
  out << "set terminal postscript eps enhanced color\n\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "set autoscale\n";	
  out << "set xtics 0.1\n";
  out << "set ytic auto\n";
  out << "set style data linespoints\n\n";
  out << "set xtics font \"Times-Roman, 20\"\n";
  out << "set ytics font \"Times-Roman, 20\"\n";
  out << "set xlabel font \"Times-Roman, 25\"\n";
  out << "set ylabel font \"Times-Roman, 25\"\n\n";
  out << "set xlabel \"" << xlabel << "\\n\"\n";
  out << "set ylabel \"" << ylabel << "\\n\"\n";
  //out << "set xr[10:]\n\n";
  out << "plot \"" << data_file << "\" using 1:2 title \"MOMENT\" lc rgb \"red\", \\\n"
      << "\"\" using 1:3 title \"MLE\" lc rgb \"blue\", \\\n"
      << "\"\" using 1:4 title \"MAP\" lc rgb \"dark-green\", \\\n"
      << "\"\" using 1:5 title \"MML\" lc rgb \"black\"\n";
  out.close();
  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}

/*void plot_script_mse_betas(
  string &data_file, string &n_str, string &type_str, string &xlabel, string &ylabel
) {
  string script_file = n_str + type_str + "_mse_betas.p";
  string plot_file = n_str + type_str + "_mse_betas.eps";
  ofstream out(script_file.c_str());
  out << "set terminal postscript eps enhanced color\n\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "set autoscale\n";	
  out << "set xtics 0.1\n";
  out << "set ytic auto\n";
  out << "set style data linespoints\n\n";
  out << "set xtics font \"Times-Roman, 20\"\n";
  out << "set ytics font \"Times-Roman, 20\"\n";
  out << "set xlabel font \"Times-Roman, 25\"\n";
  out << "set ylabel font \"Times-Roman, 25\"\n\n";
  out << "set xlabel \"" << xlabel << "\"\n";
  out << "set ylabel \"" << ylabel << "\"\n";
  //out << "set xr[10:]\n\n";
  out << "plot \"" << data_file << "\" using 1:2 title \"MOMENT\" lc rgb \"red\", \\\n"
      << "\"\" using 1:3 title \"MLE\" lc rgb \"blue\", \\\n"
      << "\"\" using 1:4 title \"MAP\" lc rgb \"dark-green\", \\\n"
      << "\"\" using 1:5 title \"MML\" lc rgb \"black\"\n";
  out.close();
  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}*/

struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;
  string quantity;

  options_description desc("Allowed options");
  desc.add_options()
       ("n",value<int>(&parameters.N),"N")
       ("kappa",value<double>(&parameters.kappa),"kappa")
       ("eccentricity",value<double>(&parameters.eccentricity),"eccentricity")
       ("all",value<string>(&quantity),"all (kappa/eccentricity) ?")
  ;
  variables_map vm;
  store(command_line_parser(argc,argv).options(desc).run(),vm);
  notify(vm);

  if (vm.count("all")) {
    if (quantity.compare("k") == 0) { // all k fixed e
      parameters.quantity = 1;
    } else if (quantity.compare("e") == 0) { // all e fixed k
      parameters.quantity = 2;
    }
  } 

  return parameters;
}

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);
  double kappa = parameters.kappa;
  double eccentricity = parameters.eccentricity;

  string data_file,script_file;
  //string n_str = "./N_" + boost::lexical_cast<string>(parameters.N) + "_uniform_prior/";
  string n_str = "./N_" + boost::lexical_cast<string>(parameters.N) + "_vmf_prior/";
  //string n_str = "./N_" + boost::lexical_cast<string>(parameters.N) + "/";

  if (parameters.quantity == 1) { // all k fixed e
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << eccentricity;
    string eccentricity_str = sse.str();
    string type_str = "e_" + eccentricity_str;

    data_file = n_str + "e_" + eccentricity_str + "_mse_kappas.dat";
    tabulate_eccentricity_mse_kappas(data_file,n_str,eccentricity_str);
    string xlabel = "Kappa";
    string ylabel = "MSE (Kappa)";
    plot_script_mse_kappas(data_file,n_str,type_str,xlabel,ylabel);

    /*data_file = n_str + "e_" + eccentricity_str + "_mse_betas.dat";
    tabulate_eccentricity_mse_betas(data_file,n_str,eccentricity,eccentricity_str);
    ylabel = "MSE (Beta)";
    plot_script_mse_betas(data_file,n_str,type_str,xlabel,ylabel);*/
  } else if (parameters.quantity == 2) { // all e fixed k
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string type_str = "k_" + kappa_str;

    data_file = n_str + "k_" + kappa_str + "_mse_kappas.dat";
    tabulate_kappa_mse_kappas(data_file,n_str,kappa,kappa_str);
    string xlabel = "eccentricity";
    string ylabel = "MSE (Kappa)";
    plot_script_mse_kappas(data_file,n_str,type_str,xlabel,ylabel);

    /*data_file = n_str + "k_" + kappa_str + "_mse_betas.dat";
    tabulate_kappa_mse_betas(data_file,n_str,kappa,kappa_str);
    ylabel = "MSE (Beta)";
    plot_script_mse_betas(data_file,n_str,type_str,xlabel,ylabel);*/
  }

  return 0;
}

