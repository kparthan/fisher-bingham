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

int minimumIndex(int ignore_map, Vector &values)
{
  int min_index = 0;
  double min_val = values[0];
  if (ignore_map) {
    for (int i=1; i<values.size(); i++) { 
      if (i != MAP) {
        if (values[i] <= min_val) {
          min_index = i;
          min_val = values[i];
        } // if()
      } 
    } // for()
  } else { // don't ignore_map
    for (int i=1; i<values.size(); i++) { 
      if (values[i] <= min_val) {
        min_index = i;
        min_val = values[i];
      } // if()
    } // for()
  } // if(ignore_map)
  return min_index;
}

void computeWins(
  int ignore_map, ostream &out, std::vector<Vector> &values
) {
  int num_elements = values.size();
  std::vector<int> wins(values[0].size(),0);

  int min_index;
  for (int i=0; i<num_elements; i++) {
    min_index = minimumIndex(ignore_map,values[i]);
    wins[min_index]++;
  }
  double percent_wins;
  for (int j=0; j<wins.size(); j++) {
    percent_wins = wins[j] * 100.0 / num_elements;
    out << fixed << setw(10) << setprecision(2) << percent_wins;
  }
  out << endl;
}

// all k fixed e
void tabulate_eccentricity_kldivs(string &data_file, string &n_str, string &eccentricity_str)
{
  std::vector<Vector> kldivs;
  string kldivs_file;
  ofstream out(data_file.c_str());
  double kappa = 10;
  while (kappa <= 100) {
    out << fixed << setw(10) << setprecision(0) << kappa;
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
    string kldivs_file = current_dir + "kldivs";
    kldivs = load_matrix(kldivs_file,NUM_METHODS);
    computeWins(0,out,kldivs);
    kappa += 10;
  } // while()
  out.close();
}

// all e fixed k
void tabulate_kappa_kldivs(string &data_file, string &n_str, string &kappa_str)
{
  std::vector<Vector> kldivs;
  string kldivs_file;
  ofstream out(data_file.c_str());
  double eccentricity = 0.1;
  while (eccentricity < 0.95) {
    out << fixed << setw(10) << setprecision(1) << eccentricity;
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << eccentricity;
    string eccentricity_str = sse.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
    kldivs_file = current_dir + "kldivs";
    kldivs = load_matrix(kldivs_file,NUM_METHODS);
    computeWins(0,out,kldivs);
    eccentricity += 0.1;
  } // while()
  out.close();
}

void plot_script_kldivs(string &data_file, string &n_str, string &type_str)
{
  string script_file = n_str + type_str + "_kldivs.p";
  string plot_file = n_str + type_str + "_kldivs.eps";
  ofstream out(script_file.c_str());
  out << "set terminal postscript eps enhanced color\n\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "set key invert reverse left top\n\n";
  out << "set grid y\n";
  out << "set style data histograms\n";
  out << "set style histogram rowstacked\n";
  out << "set boxwidth 0.5\n";
  out << "set style fill solid 1.0 border -1\n";
  out << "set ytics 10 nomirror\n";
  out << "set yrange [:100]\n";
  out << "set ylabel \"\% of wins\"\n";
  out << "set ytics 10\n\n"; 
  out << "plot \"" << data_file << "\" using 2 t \"MOM\", \\\n"
      << "\"\" using 3 t \"MLE\", \\\n"
      << "\"\" using 4 t \"MAP\", \\\n"
      << "\"\" using 5:xtic(1) t \"MML\"";
  out.close();
  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
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
  //string n_str = "./N_" + boost::lexical_cast<string>(parameters.N) + "_vmf_prior/";
  string n_str = "./N_" + boost::lexical_cast<string>(parameters.N) + "_beta_prior/";
  //string n_str = "N_" + boost::lexical_cast<string>(parameters.N) + "_kappa_until_50/";
  //string n_str = "./N_" + boost::lexical_cast<string>(parameters.N) + "/";

  if (parameters.quantity == 1) { // all k fixed e
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << eccentricity;
    string eccentricity_str = sse.str();
    data_file = n_str + "e_" + eccentricity_str + "_kldivs.dat";
    //script_file = n_str + "e_" + eccentricity_str + "_kldivs.p";
    tabulate_eccentricity_kldivs(data_file,n_str,eccentricity_str);
    string type_str = "e_" + eccentricity_str;
    plot_script_kldivs(data_file,n_str,type_str);
  } else if (parameters.quantity == 2) { // all e fixed k
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    data_file = n_str + "k_" + kappa_str + "_kldivs.dat";
    tabulate_kappa_kldivs(data_file,n_str,kappa_str);
    string type_str = "k_" + kappa_str;
    plot_script_kldivs(data_file,n_str,type_str);
  }

  return 0;
}

