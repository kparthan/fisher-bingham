#ifndef VMC_H
#define VMC_H

#include "Header.h"

class vMC
{
  private:
    //! (unit) mean of the distribution
		Vector mu;

    //! Concentration parameter 
		long double kappa;

  public:
		//! Constructor
		vMC();

		//! Constructor that sets value of parameters
		vMC(Vector &, long double);

		vMC(long double);

    //! Assignment of an existing vMC distribution
    vMC operator=(const vMC &);

		//! Gets the mean 
		Vector Mean();

    //! Gets the Kappa 
    long double Kappa(); 

    //! Generate a random canonical sample
    void generateCanonical(std::vector<Vector> &, int);
};

#endif

