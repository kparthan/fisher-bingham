#ifndef VMF_H
#define VMF_H

#include "Header.h"

class vMF
{
  private:
    //! (unit) mean of the distribution
		Vector mu,kmu;

    double theta,phi;

    //! Concentration parameter 
		double kappa;

    //! Normalization constant
    double log_cd;

  protected:
    //! Computes the vMF constants
    void updateConstants();

    //! Computes the normalization constant
    double computeLogNormalizationConstant();

  public:
		//! Constructor
		vMF();

		//! Constructor that sets value of parameters
		vMF(Vector &, double);

		vMF(double);

    //! Updates the model parameters
    void updateParameters();

    //! Assignment of an existing vMF distribution
    vMF operator=(const vMF &);

		//! Gets the mean 
		Vector Mean();

    //! Gets the Kappa 
    double Kappa(); 

    //! Gets the normalization constant
    double getLogNormalizationConstant();

		//! Function value
		double density(Vector &);

    //! Computes the log of probability density
		double log_density(Vector &);

    //! Computes the negative log likelihood of a sample
    double computeNegativeLogLikelihood(Vector &);

    //! Computes the negative log likelihood of a sample
    double computeNegativeLogLikelihood(std::vector<Vector> &);

    double computeNegativeLogLikelihood(double, double);

    //! Prints the model parameters
    void printParameters(ostream &);

    //! Generate random sample
    std::vector<Vector> generate(int, int D = 3);

    //! Generate a random canonical sample
    void generateCanonical(std::vector<Vector> &, int, int D = 3);

    double computeLogPriorProbability();

    double computeLogPriorMean();

    double computeLogPriorScale();

    double computeLogFisherInformation();

    double computeLogFisherInformation(double);

    void computeAllEstimators(std::vector<Vector> &);

    void computeAllEstimators(std::vector<Vector> &, std::vector<struct Estimates_vMF> &);

    void estimateMean(struct Estimates_vMF &, std::vector<Vector> &, Vector &);

    void estimateMLApproxKappa(struct Estimates_vMF &);

    struct Estimates_vMF computeMMLEstimates(std::vector<Vector> &);

    struct Estimates_vMF computeMMLEstimates(struct Estimates_vMF &);

    void estimateParameters(std::vector<Vector> &, Vector &);

    void updateParameters(struct Estimates_vMF &);

    double computeMessageLength(std::vector<Vector> &);

    double computeMessageLength(double, double);

    double computeMessageLength(struct Estimates_vMF &);

    double computeKLDivergence(vMF &);

    double computeKLDivergence(struct Estimates_vMF &);
};

#endif

