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

    void exp1();
    void generate_data_exp1(int, int, int);
    void infer_components_exp1(int, int);

    void traditional_search(
      string &,
      Mixture &, std::vector<Vector> &, std::vector<Vector> &,
      string &, string &,
      ostream &, ostream &, ostream &
    );

    void proposed_search(
      string &, Mixture &, std::vector<Vector> &, std::vector<Vector> &, string &
    );
    void proposed_search(
      string &, 
      Mixture &, std::vector<Vector> &, std::vector<Vector> &, 
      string &, string &, string &, string &
    );
};

#endif

