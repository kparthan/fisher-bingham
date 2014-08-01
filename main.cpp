#include "Support.h"

int main(int argc, char **argv)
{
  srand(time(NULL));

  struct Parameters parameters = parseCommandLineInput(argc,argv);

  if (parameters.test == SET) {
    Test();
  }

  return 0;
}

