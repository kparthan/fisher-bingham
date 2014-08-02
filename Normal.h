#ifndef NORMAL_H
#define NORMAL_H

#include "Header.h"

class Normal
{
  private:
    //! Mean of the distribution
		long double mu;

    //! Standard deviation of the distribution
		long double sigma;

  public:
		//! Constructor
		Normal() ;

		//! Constructor that sets value of parameters
		Normal(long double, long double);

    //! Assignment of an existing Normal distribution
    Normal operator=(const Normal &);

		//! Gets the mean 
		const long double mean();

    //! Gets the standard deviation
    const long double standardDeviation(); 

		//! Function value
		long double density(long double);

    long double cumulativeDensity(long double);

    //! Computes the negative log likelihood of a sample
    long double negativeLogLikelihood(long double);

    //! Computes the negative log likelihood of a sample
    long double negativeLogLikelihood(std::vector<long double> &);

    //! Prints the model parameters
    void printParameters(ostream &);

    //! Generate random sample
    std::vector<long double> generate(int);

};

#endif

