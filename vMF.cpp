#include "vMF.h"
#include "Support.h"
#include "Normal.h"

/*!
 *  \brief This is a constructor module
 */
vMF::vMF()
{
  D = 3;
  mu = std::vector<long double>(3,0); mu[2] = 1;
  kappa = 1;
  updateConstants();
}

/*!
 *  \brief constructor function which sets the value of mean and 
 *  kappa of the distribution
 *  \param mu a reference to a vector<long double>
 *  \param kappa a long double
 */
vMF::vMF(std::vector<long double> &mu, long double kappa) : mu(mu), kappa(kappa)
{
  D = mu.size();
  updateConstants();
}

/*!
 *  \brief This function computes the constants associated with the vMF.
 */
void vMF::updateConstants()
{
  kmu = std::vector<long double>(D,0);
  for (int i=0; i<D; i++) {
    kmu[i] = kappa * mu[i];
  }
  log_cd = computeLogNormalizationConstant();
  cd = exp(log_cd);
}

/*!
 *  \brief This function computes the normalization constant of the distribution.
 *  \return the normalization constant
 */
long double vMF::computeLogNormalizationConstant()
{
  if (kappa < ZERO) {
    long double log_area = computeLogSurfaceAreaSphere(D);
    return log_area;
  } else {
    long double log_bessel = logModifiedBesselFirstKind(D/2.0-1,kappa);
    if (log_bessel >= INFINITY) {
      my_float d2_1 = (D / 2.0) - 1;
      my_float bessel = cyl_bessel_i(d2_1,kappa);
      my_float large_log_bessel = log(bessel);
      log_bessel = large_log_bessel.convert_to<long double>();
      cout << "Normalization constant error ...\n";
      cout << scientific << "bessel: " << bessel;
      cout << scientific << "; log_bessel: " << log_bessel << endl;
      cout << "D: " << D << "; kappa: " << kappa << endl;
    }
    if (D == 2) {
      long double log_cd = -log(2*PI) - log_bessel;
    } else {
      long double log_tmp = (D/2.0) * log (kappa/(2*PI));
      long double log_cd = log_tmp - log(kappa) - log_bessel;
    }
    return log_cd;
  }
}

/*!
 *  \brief This function assigns a source vMF distribution.
 *  \param source a reference to a vMF
 */
vMF vMF::operator=(const vMF &source)
{
  if (this != &source) {
    D = source.D;
    mu = source.mu;
    kmu = source.kmu;
    kappa = source.kappa;
    cd = source.cd;
    log_cd = source.log_cd;
  }
  return *this;
}

/*!
 *  \brief This function returns the mean of the distribution
 *  \return the mean of the distribution
 */
std::vector<long double> vMF::mean(void)
{
	return mu;
}

/*!
 *  \brief This function returns the kappa of the distribution
 *  \return the kappa of the distribution
 */
long double vMF::Kappa(void)
{
	return kappa;
}

/*!
 *  \brief This function returns the normalization constant of the distribution
 *  \return the normalization constant of the distribution
 */
long double vMF::getNormalizationConstant()
{
  return cd;
}

long double vMF::getLogNormalizationConstant()
{
  return log_cd;
}

/*!
 *  \brief This function returns the dimensionality of the data.
 */
int vMF::getDimensionality()
{
  return D;
}

/*!
 *  \brief This function computes the value of the distribution at a given x
 *  \param x a reference to a vector<long double>
 *  \return density of the function given x
 */
long double vMF::density(std::vector<long double> &x)
{
  long double value = log_density(x);
  return exp(value);
}

/*!
 *  \brief This function computes the value of the distribution at a given x
 *  \param x a reference to a vector<long double>
 *  \return density of the function given x
 */
long double vMF::log_density(std::vector<long double> &x)
{
  long double expnt = computeDotProduct(kmu,x);
  long double log_density = log_cd + expnt;// + (D-1) * log(AOM);
  return log_density;
}

/*!
 *  \brief This function computes the negative log likelihood of given datum.
 *  \param x a reference to a vector<long double>
 *  \return the negative log likelihood (base e)
 */
long double vMF::negativeLogLikelihood(std::vector<long double> &x)
{
  return -log_density(x);
}

/*!
 *  \brief This function computes the negative log likelihood of given data.
 *  \param sample a reference to a vector<vector<long double> >
 *  \return the negative log likelihood (base e)
 */
long double vMF::negativeLogLikelihood(std::vector<std::vector<long double> > &sample)
{
  long double value = 0;
  int N = sample.size();
  value -= N * log_cd;
  std::vector<long double> sum(D,0);
  for (int i=0; i<N; i++) {
    for (int j=0; j<D; j++) {
      sum[j] += sample[i][j];
    }
  }
  value -= computeDotProduct(kmu,sum);
  return value;
}

/*!
 *  \brief This function prints the parameters of the model.
 */
void vMF::printParameters(ostream &os)
{
  os << "[mu]: ";
  print(os,mu,3);
  os << "\t[kappa]:" << setw(10) << setprecision(3) << kappa << endl;
}

/*!
 *  \brief This function generates a random sample from a canonical vMF 
 *  based on Wood's algorithm. Steps involved are as follows:
 *  Step 0:
 *          Calculate b = [-2K + {4 K^2 + (D-1)^2}^1/2] / (D-1)
 *               put x0 = (1 - b) / (1 + b)
 *               and c  = K x0 + (D-1) log(1 - x0^2)
 *  Step 1:
 *          Generate Z ~ Beta{(D-1)/2,(D-1)/2}
 *                   U ~ Uniform(0,1)
 *     and calculate W = {1-(1+b)Z} / {1-(1-b)Z}
 *  Step 2:
 *          If KW + (D-1) log (1-x0 W) - c < log(U), then go to Step 1
 *  Step 3:
 *          Generate a uniform (D-1)-dimensional unit vector V, and
 *          return X^T = ((1-W^2)^1/2 V^T,W)
 *  Then X has the vMF distribution with mean direction (0,...,0,1)^T and Kappa = K
 *  \param canonical_sample a reference to a vector<vector<long double> >
 *  \param sample_size an integer
 *  \return the random list of points
 */
void vMF::generateCanonical(std::vector<std::vector<long double> > &canonical_sample, int sample_size)
{
  canonical_sample.clear();
  int count = 0;
  beta_distribution<> beta((D-1)/2.0,(D-1)/2.0);
  Normal normal(0,1);
  std::vector<long double> V(D-1,0);  // (D-1)-unit vector
  std::vector<long double> random_vmf(D,0);

  // step 0
  long double tmp1,tmp2,p,Z,U,W,check;
  tmp1 = (4 * kappa * kappa) + (D-1) * (D-1);
  long double b = (-2 * kappa + sqrt(tmp1)) / (long double) (D-1);
  long double x0 =  (1 - b) / (1 + b);
  long double c = (kappa * x0) + ((D-1) * log(1 - x0*x0));

  while (count < sample_size) {
    // step 1
    p = rand() / (long double) RAND_MAX;
    //cout << "p: " << p << endl;
    Z = quantile(beta,p);
    U = rand() / (long double) RAND_MAX;
    //cout << "U: " << U << endl;
    tmp1 = 1 - ((1+b) * Z);
    tmp2 = 1 - ((1-b) * Z);
    W = tmp1 / tmp2;
    //cout << "W: " << W << endl;

    // step 2
    check = kappa * W + ((D-1) * log(1 - x0*W)) - c;
    if (check >= log(U)) {
      // step 3
      std::vector<long double> random_normal = normal.generate(D-1);
      normalize(random_normal,V);
      //print(cout,V);
      tmp1 = sqrt(1-W*W);
      for (int i=0; i<D-1; i++) {
        random_vmf[i] = tmp1 * V[i];
      }
      random_vmf[D-1] = W;
      canonical_sample.push_back(random_vmf);
      count++;
    }
  }
}

/*!
 *  \brief This function generates a random sample from vMF distribution.
 *  \param sample_size an integer
 *  \return the random list of points
 */
std::vector<std::vector<long double> > vMF::generate(int sample_size)
{
  cout << "\nGenerating from vMF with mean: ";
  if (D == 3) {
    std::vector<long double> spherical(3,0);
    cartesian2spherical(mu,spherical);
    spherical[1] *= 180 / PI;
    spherical[2] *= 180 / PI;
    print(cout,spherical,3);
    cout << "; Kappa: " << kappa << "; sample size = " << sample_size << endl;
  } else {
    cout << "[D,K] = [" << D << "," << kappa 
         << "] and sample size = " << sample_size;
  }
  if (sample_size != 0) {
    std::vector<std::vector<long double> > canonical_sample;
    generateCanonical(canonical_sample,sample_size);
    if (fabs(mu[D-1] - 1) <= TOLERANCE) {  // check if mu is Z-axis
      return canonical_sample;
    } else {
      matrix<long double> transformation = align_zaxis_with_vector(mu);
      return transform(canonical_sample,transformation);
    }
  } else if (sample_size == 0) {
    return std::vector<std::vector<long double> >(); 
  }
}

