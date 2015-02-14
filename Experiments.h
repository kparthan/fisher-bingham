#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "Header.h"
#include "Mixture.h"

class Experiments
{
  private:
    int iterations;

  public:
    Experiments(int);

    void simulate();

    void checkFolders(string &);

    void create_sub_folders(string &, string &);

    void infer_components_exp1();
    void infer_components_exp1(double, double);

    void generateData(Mixture &, string &, int);

    void inferMixtures(Mixture &, string &);

    void inferMixtures(
      Mixture &, string &, string &, string &, string &, string &
    );
};

#endif

