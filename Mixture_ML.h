#ifndef MIXTURE_ML_H
#define MIXTURE_ML_H

#include "Kent.h"

class Mixture_ML
{
  private:
    //! ID
    int id;

    //! Sample size
    int N;

    //! Number of components
    int K;

    //! List of components
    std::vector<Kent> components;
    
    //! Sample (x_i) -- Cartesian coordinates
    //! (on a sphere of unit radius)
    std::vector<Vector> data;

    //! Data weights
    Vector data_weights;

    //! Responsibility matrix (K X N)
    std::vector<Vector> responsibility;

    //! Effective sample size for each component (n_k)
    Vector sample_size;

    //! Weights of the components (a_k)
    Vector weights;

    //! Optimal encoding length
    double negloglike;

  public:
    //! Null constructor
    Mixture_ML();

    //! Constructor
    Mixture_ML(int, std::vector<Kent> &, Vector &);

    //! Constructor
    Mixture_ML(int, std::vector<Vector> &, Vector &);

    //! Constructor
    Mixture_ML(int, std::vector<Kent> &, Vector &, Vector &, 
               std::vector<Vector> &, std::vector<Vector> &, Vector &);

    //! Overloading = operator
    Mixture_ML operator=(const Mixture_ML &);

    //! Overloading == operator
    bool operator==(const Mixture_ML &);

    //! Prepare log file
    string getLogFile();

    //! Gets the list of weights
    Vector getWeights();

    //! Gets the list of components
    std::vector<Kent> getComponents();

    //! Returns number of components
    int getNumberOfComponents();

    //! Gets the responsibility matrix
    std::vector<Vector> getResponsibilityMatrix();

    //! Gets the sample size
    Vector getSampleSize();

    //! Initialize parameters
    void initialize();

    //! Updates the effective sample size
    void updateEffectiveSampleSize();

    void updateWeights_ML();

    //! Update components
    void updateComponents();

    //! Update the responsibility matrix
    void updateResponsibilityMatrix();

    //! Computes the responsibility matrix
    void computeResponsibilityMatrix(std::vector<Vector> &, string &);
                                          
    //! Probability of a datum
    double log_probability(Vector &);

    //! Computes the negative log likelihood
    double computeNegativeLogLikelihood(std::vector<Vector> &);

    //! Gets the minimum message length
    double getNegativeLogLikelihood();

    //! Estimate mixture parameters
    double estimateParameters();

    //! EM loop
    void EM();

    //! Prints the model parameters
    void printParameters(ostream &, int, double);

    //! Prints the model parameters
    void printParameters(ostream &, int);

    //! Prints the model parameters
    void printParameters(ostream &);

    //! Loads the mixture file
    void load(string &);

    //! Loads the mixture file with the corresponding data
    void load(string &, std::vector<Vector> &, Vector &);

    //! Randomly choose a component
    int randomComponent();

    //! Saves the data generated from a component
    void saveComponentData(int, std::vector<Vector> &);

    //! Generate random data from the distribution using mixture proportions
    std::vector<Vector> generate(int, bool);

    //! Splits a component
    Mixture_ML split(int, ostream &);

    //! Deltes a component
    Mixture_ML kill(int, ostream &);

    //! Joins two  components
    Mixture_ML join(int, int, ostream &);

    //! Generate heat maps (for d=3)
    void generateHeatmapData(double);

    //! Get the nearest component
    int getNearestComponent(int);

    //! Computes the approx KL divergence between two mixtures
    double computeKLDivergence(Mixture_ML &);

    double computeKLDivergence(Mixture_ML &, std::vector<Vector> &);
};

#endif

