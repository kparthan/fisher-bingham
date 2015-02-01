#include <iostream>
#include <cstdlib>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost::program_options;

double scale_to_aom(double x, double aom)
{
  cout << "x: " << x << endl;

  int aom_inv = 1.0/aom;
  cout << "aom_inv: " << aom_inv << endl;

  int tmp = x * aom_inv;
  cout << "tmp: " << tmp << endl;

  double y = tmp * aom;
  cout << "y: " << y << endl;
  return y;
}

struct Parameters
{
  double x,aom;
};

struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;

  options_description desc("Allowed options");
  desc.add_options()
       ("x",value<double>(&parameters.x),"value to be scaled")
       ("aom",value<double>(&parameters.aom),"value to be scaled")
  ;
  variables_map vm;
  store(command_line_parser(argc,argv).options(desc).run(),vm);
  notify(vm);

  if (!vm.count("aom")) {
    parameters.aom = 0.001; // default AOM
  }

  return parameters;
}

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  double y = scale_to_aom(parameters.x,parameters.aom);

  cout << "x: " << parameters.x << "; scaled: " << y << endl; 

  return 0;
}

