#include "Support.h"

extern std::vector<long double> XAXIS,YAXIS,ZAXIS;

int main(int argc, char **argv)
{
  srand(time(NULL));

  XAXIS = std::vector<long double>(3,0); XAXIS[0] = 1;
  YAXIS = std::vector<long double>(3,0); YAXIS[1] = 1;
  ZAXIS = std::vector<long double>(3,0); ZAXIS[2] = 1;

  struct Parameters parameters = parseCommandLineInput(argc,argv);

  if (parameters.test == SET) {
    TestFunctions();
  }

  return 0;
}

