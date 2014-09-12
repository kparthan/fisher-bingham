#include "Support.h"

extern Vector XAXIS,YAXIS,ZAXIS;

int main(int argc, char **argv)
{
  srand(time(NULL));

  XAXIS = Vector(3,0); XAXIS[0] = 1;
  YAXIS = Vector(3,0); YAXIS[1] = 1;
  ZAXIS = Vector(3,0); ZAXIS[2] = 1;

  struct Parameters parameters = parseCommandLineInput(argc,argv);

  if (parameters.test == SET) {
    TestFunctions();
  }

  if (parameters.experiments == SET) {
    RunExperiments(parameters.iterations);
  }

  if (parameters.read_profiles == SET && parameters.simulation == UNSET) {
    computeEstimators(parameters);
  } 

  if (parameters.simulation == SET) {
    //simulateMixtureModel(parameters);
  }

  return 0;
}

