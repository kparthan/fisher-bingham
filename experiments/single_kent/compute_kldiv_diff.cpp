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

std::vector<Vector> load_table(string &file_name, int D)
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

struct Parameters
{
  int N;
  double kappa,eccentricity;
  int quantity;
};

struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;
  string quantity;

  options_description desc("Allowed options");
  desc.add_options()
       ("n",value<int>(&parameters.N),"N")
       ("kappa",value<double>(&parameters.kappa),"kappa")
       ("eccentricity",value<double>(&parameters.eccentricity),"eccentricity")
       ("all",value<string>(&quantity),"all (k/e) ?")
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

void process_kldivs(
  std::vector<Vector> &kldivs, 
  Vector &kldivs_wins, 
  Vector &kldivs_losses
) {
  double diff;
  kldivs_wins = Vector(3,0);
  kldivs_losses = Vector(3,0);
  for (int i=0; i<kldivs.size(); i++) {
    for (int j=0; j<NUM_METHODS-1; j++) {
      diff = kldivs[i][j] - kldivs[i][MML];
      if (diff < 0) {  // MML loses ...
        kldivs_wins[j] -= diff;
      } else {  // MML wins ...
        kldivs_losses[j] -= diff;
      } // if ()
    } // j ()
  } // i ()
}

// Usage: ./compute_kldiv_diff --n 100 --kappa 10 --all e
// Usage: ./compute_kldiv_diff --n 100 --eccentricity 0.1 --all k
int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);
  double kappa = parameters.kappa;
  double eccentricity = parameters.eccentricity;
  //string n_str = "N_" + boost::lexical_cast<string>(parameters.N) + "_new2_prior/";
  string n_str = "N_" + boost::lexical_cast<string>(parameters.N) + "_new32_prior/";

  if (parameters.quantity == 1) { // all k fixed e
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << eccentricity;
    string eccentricity_str = sse.str();
    Vector kldivs_wins,kldivs_losses;
    string output = n_str + "kldivs_diff_e_" + eccentricity_str + ".dat";
    ofstream out(output.c_str());
    Vector kappas(3,0); kappas[0] = 10; kappas[1] = 50; kappas[2] = 100;
    //while (kappa <= 100) {
    for (int i=0; i<kappas.size(); i++) {
      kappa = kappas[i];
      ostringstream ssk;
      ssk << fixed << setprecision(0);
      ssk << kappa;
      string kappa_str = ssk.str();
      string current_dir = n_str + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
      string kldivs_file = current_dir + "kldivs";
      std::vector<Vector> kldivs_table = load_table(kldivs_file,NUM_METHODS);
      process_kldivs(kldivs_table,kldivs_wins,kldivs_losses);
      out << fixed << setw(10) << kappa << "\t\t";
      for (int i=0; i<3; i++) {
        out << fixed << scientific << setprecision(6) << kldivs_wins[i] << "\t";
      }
      for (int i=0; i<3; i++) {
        out << fixed << scientific << setprecision(6) << kldivs_losses[i] << "\t";
      }
      out << endl;
    } // for ()
    out.close();
  } else if (parameters.quantity == 2) { // all e fixed k
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    Vector kldivs_wins,kldivs_losses;
    string output = n_str + "kldivs_diff_k_" + kappa_str + ".dat";
    ofstream out(output.c_str());
    eccentricity = 0.1;
    while (eccentricity <= 0.95) {
      ostringstream sse;
      sse << fixed << setprecision(1);
      sse << eccentricity;
      string eccentricity_str = sse.str();
      string current_dir = n_str + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
      string kldivs_file = current_dir + "kldivs";
      std::vector<Vector> kldivs_table = load_table(kldivs_file,NUM_METHODS);
      process_kldivs(kldivs_table,kldivs_wins,kldivs_losses);
      out << fixed << setw(10) << setprecision(1) << eccentricity << "\t\t";
      for (int i=0; i<3; i++) {
        out << fixed << scientific << setprecision(6) << kldivs_wins[i] << "\t";
      }
      for (int i=0; i<3; i++) {
        out << fixed << scientific << setprecision(6) << kldivs_losses[i] << "\t";
      }
      out << endl;
      eccentricity += 0.4;
    }
    out.close();
  }
}

