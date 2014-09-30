#ifndef VMF_H
#define VMF_H

#include "Header.h"

class vMF
{
  private:
    //! (unit) mean of the distribution
		Vector mu,kmu;

    long double theta,phi;

    //! Concentration parameter 
		long double kappa;

    //! Normalization constant
    long double log_cd;

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

		vMF(long double);

    //! Updates the model parameters
    void updateParameters();

    //! Assignment of an existing vMF distribution
    vMF operator=(const vMF &);

		//! Gets the mean 
		Vector Mean();

    //! Gets the Kappa 
    long double Kappa(); 

    //! Gets the normalization constant
    long double getLogNormalizationConstant();

		//! Function value
		long double density(Vector &);

    //! Computes the log of probability density
		long double log_density(Vector &);

    //! Computes the negative log likelihood of a sample
    long double computeNegativeLogLikelihood(Vector &);

    //! Computes the negative log likelihood of a sample
    long double computeNegativeLogLikelihood(std::vector<Vector> &);

    long double computeNegativeLogLikelihood(long double, long double);

    //! Prints the model parameters
    void printParameters(ostream &);

    //! Generate random sample
    std::vector<Vector> generate(int, int D = 3);

    //! Generate a random canonical sample
    void generateCanonical(std::vector<Vector> &, int, int D = 3);

    long double computeLogPriorProbability();

    long double computeLogPriorMean();

    long double computeLogPriorScale();

    long double computeLogFisherInformation();

    long double computeLogFisherInformation(long double);

    void computeAllEstimators(std::vector<Vector> &);

    void computeAllEstimators(std::vector<Vector> &, std::vector<struct Estimates_vMF> &);

    void estimateMean(struct Estimates_vMF &, std::vector<Vector> &, Vector &);

    void estimateMLApproxKappa(struct Estimates_vMF &);

    struct Estimates_vMF computeMMLEstimates(std::vector<Vector> &);

    struct Estimates_vMF computeMMLEstimates(struct Estimates_vMF &);

    void estimateParameters(std::vector<Vector> &, Vector &);

    void updateParameters(struct Estimates_vMF &);

    long double computeMessageLength(std::vector<Vector> &);

    long double computeMessageLength(long double, long double);

    long double computeMessageLength(struct Estimates_vMF &);

    long double computeKLDivergence(vMF &);

    long double computeKLDivergence(struct Estimates_vMF &);
};

#endif

