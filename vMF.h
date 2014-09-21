#ifndef VMF_H
#define VMF_H

#include "Header.h"

class vMF
{
  private:
    //! Dimensionality
    int D;

    //! (unit) mean of the distribution
		Vector mu;

    //! Concentration parameter 
		long double kappa;

    //! Normalization constant
    long double cd,log_cd;

    //! kappa * mu
    Vector kmu;

  protected:
    //! Computes the vMF constants
    void updateConstants();

    //! Computes the normalization constant
    long double computeLogNormalizationConstant();

  public:
		//! Constructor
		vMF();

		//! Constructor that sets value of parameters
		vMF(Vector &, long double);

    //! Updates the model parameters
    void updateParameters();

    //! Assignment of an existing vMF distribution
    vMF operator=(const vMF &);

		//! Gets the mean 
		Vector mean();

    //! Gets the Kappa 
    long double Kappa(); 

    //! Gets the normalization constant
    long double getNormalizationConstant();
    long double getLogNormalizationConstant();

    //! Gets the dimensionality of the data
    int getDimensionality();

		//! Function value
		long double density(Vector &);

    //! Computes the log of probability density
		long double log_density(Vector &);

    //! Computes the negative log likelihood of a sample
    long double negativeLogLikelihood(Vector &);

    //! Computes the negative log likelihood of a sample
    long double negativeLogLikelihood(std::vector<Vector> &);

    //! Prints the model parameters
    void printParameters(ostream &);

    //! Generate random sample
    std::vector<Vector> generate(int);

    //! Generate a random canonical sample
    void generateCanonical(std::vector<Vector> &, int);

};

#endif

