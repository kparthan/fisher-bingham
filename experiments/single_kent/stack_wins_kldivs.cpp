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

struct Parameters
{
  double kappa,eccentricity;
  int quantity;
};

struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;
  string quantity;

  options_description desc("Allowed options");
  desc.add_options()
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

  string parent_folder = "./N_100/";

  if (parameters.quantity == 1) { // all k fixed e

  } else if (parameters.quantity == 2) { // all e fixed k
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();

    eccentricity = 0.1;
    while (eccentricity < 0.95) {
      beta = 0.5 * kappa * eccentricity;
      ostringstream sse;
      sse << fixed << setprecision(1);
      sse << eccentricity;
      eccentricity_str = sse.str();
      current_dir = parent_dir + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
  }

  return 0;
}

