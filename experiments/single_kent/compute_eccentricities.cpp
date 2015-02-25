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

std::vector<Vector> compute_eccentricities(
  std::vector<Vector> &kappas,
  std::vector<Vector> &betas
) {
  Vector emptyvec(NUM_METHODS,0);
  std::vector<Vector> ecc(kappas.size(),emptyvec);
  double k,b;

  for (int i=0; i<kappas.size(); i++) {
    for (int j=0; j<NUM_METHODS; j++) {
      k = kappas[i][j];
      b = betas[i][j];
      ecc[i][j] = 2 * b / k;
    }
  }

  return ecc;
}

void writeToFile(string &file_name, std::vector<Vector> &v)
{
  ofstream file(file_name.c_str());
  for (int i=0; i<v.size(); i++) {
    for (int j=0; j<v[i].size(); j++) {
      file << fixed << setprecision(6) << scientific << v[i][j] << "\t";
    }
    file << endl;
  }
  file.close(); 
}

void compute_eccentricities(string &n_str) 
{
  double kappa,beta,eccentricity;
  string kappa_str,eccentricity_str,current_dir,kappas_file,betas_file,ecc_file;
  std::vector<Vector> kappas_table,betas_table,ecc_table;

  kappa = 10;
  while (kappa <= 100) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    kappa_str = ssk.str();
    eccentricity = 0.1;
    while (eccentricity < 0.95) {
      ostringstream sse;
      sse << fixed << setprecision(1);
      sse << eccentricity;
      eccentricity_str = sse.str();
      current_dir = n_str + "k_" + kappa_str + "_e_" + eccentricity_str + "/";

      kappas_file = current_dir + "kappas";
      kappas_table = load_table(kappas_file,NUM_METHODS);
      betas_file = current_dir + "betas";
      betas_table = load_table(betas_file,NUM_METHODS);

      ecc_table = compute_eccentricities(kappas_table,betas_table);
      ecc_file = current_dir + "ecc";
      writeToFile(ecc_file,ecc_table);

      eccentricity += 0.1;
    } // while (eccentricity)
    kappa += 10;
  } // while (kappa)
}

int main(int argc, char **argv)
{
  int N;
  string n_str;

  N = 100;
  n_str = "N_" + boost::lexical_cast<string>(N) + "_prior1/";
  compute_eccentricities(n_str);
  n_str = "N_" + boost::lexical_cast<string>(N) + "_prior2/";
  compute_eccentricities(n_str);
  n_str = "N_" + boost::lexical_cast<string>(N) + "_prior3/";
  compute_eccentricities(n_str);
  n_str = "N_" + boost::lexical_cast<string>(N) + "_prior4/";
  compute_eccentricities(n_str);

  N = 1000;
  n_str = "N_" + boost::lexical_cast<string>(N) + "_prior1/";
  compute_eccentricities(n_str);
  n_str = "N_" + boost::lexical_cast<string>(N) + "_prior2/";
  compute_eccentricities(n_str);
  n_str = "N_" + boost::lexical_cast<string>(N) + "_prior3/";
  compute_eccentricities(n_str);
  n_str = "N_" + boost::lexical_cast<string>(N) + "_prior4/";
  compute_eccentricities(n_str);
}

