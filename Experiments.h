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

    void fisher_uncertainty();

    void simulate();

    void checkFolders(string &);

    void create_sub_folders(string &, string &);

    void infer_components_exp1();
    void infer_components_exp1(double, double);

    void generateData(Mixture &, string &, int);

    void inferMixtures_exp1(Mixture &, string &);
    void inferMixtures_exp1(
      Mixture &, std::vector<Vector> &, string &, string &, string &, string &, string &
    );

    void infer_components_exp2();
    void inferMixtures_exp2(
      string &, Mixture &, std::vector<Vector> &, std::vector<Vector> &, string &
    );
    void inferMixtures_exp2(
      string &, 
      Mixture &, std::vector<Vector> &, std::vector<Vector> &, 
      string &, string &, string &, string &
    );

    void infer_components_exp3();
    void infer_components_exp3(double, double, double);
};

#endif

