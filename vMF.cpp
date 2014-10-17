#include "vMF.h"
#include "Support.h"
#include "Normal.h"
#include "Optimize2.h"

/*!
 *  \brief This is a constructor module
 */
vMF::vMF()
{
  mu = Vector(3,0); mu[2] = 1;
  kappa = 1;
  updateConstants();
}

/*!
 *  \brief constructor function which sets the value of mean and 
 *  kappa of the distribution
 *  \param mu a reference to a vector<long double>
 *  \param kappa a long double
 */
vMF::vMF(Vector &mu, long double kappa) : mu(mu), kappa(kappa)
{
  updateConstants();
}

vMF::vMF(long double kappa) : kappa(kappa)
{
  mu = Vector(3,0); mu[2] = 1;
  updateConstants();
}

/*!
 *  \brief This function computes the constants associated with the vMF.
 */
void vMF::updateConstants()
{
  kmu = Vector(3,0);
  for (int i=0; i<3; i++) {
    kmu[i] = kappa * mu[i];
  }
  log_cd = computeLogNormalizationConstant();
  Vector spherical(3,0);
  cartesian2spherical(mu,spherical);
  theta = spherical[1];
  phi = spherical[2];
}

/*!
 *  \brief This function computes the normalization constant of the distribution.
 *  \return the normalization constant
 */
long double vMF::computeLogNormalizationConstant()
{
  if (kappa < ZERO) {
    long double log_area = computeLogSurfaceAreaSphere(3);
    return log_area;
  } else {
    log_cd = log(kappa) - log(2*PI) - kappa;
    long double tmp = 1 - exp(-2*kappa);
    log_cd -= log(tmp);
    /*log_cd = log(kappa) - log(2*PI);
    long double tmp = exp(kappa) - exp(-kappa);
    log_cd -= log(tmp);*/
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
    mu = source.mu;
    kmu = source.kmu;
    theta = source.theta;
    phi = source.phi;
    kappa = source.kappa;
    log_cd = source.log_cd;
  }
  return *this;
}

/*!
 *  \brief This function returns the mean of the distribution
 *  \return the mean of the distribution
 */
Vector vMF::Mean(void)
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

long double vMF::getLogNormalizationConstant()
{
  return log_cd;
}

/*!
 *  \brief This function computes the value of the distribution at a given x
 *  \param x a reference to a vector<long double>
 *  \return density of the function given x
 */
long double vMF::density(Vector &x)
{
  long double value = log_density(x);
  return exp(value);
}

/*!
 *  \brief This function computes the value of the distribution at a given x
 *  \param x a reference to a vector<long double>
 *  \return density of the function given x
 */
long double vMF::log_density(Vector &x)
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
long double vMF::computeNegativeLogLikelihood(Vector &x)
{
  return -log_density(x);
}

/*!
 *  \brief This function computes the negative log likelihood of given data.
 *  \param sample a reference to a vector<vector<long double> >
 *  \return the negative log likelihood (base e)
 */
long double vMF::computeNegativeLogLikelihood(std::vector<Vector> &sample)
{
  long double value = 0;
  int N = sample.size();
  value -= N * log_cd;
  Vector sum = computeVectorSum(sample);
  value -= computeDotProduct(kmu,sum);
  return value;
}

long double vMF::computeNegativeLogLikelihood(long double R, long double N)
{
  long double value = 0;
  value -= N * log_cd;
  value -= kappa * R;
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
void vMF::generateCanonical(std::vector<Vector> &canonical_sample, int sample_size)
{
  int D = 3;
  canonical_sample.clear();
  int count = 0;
  beta_distribution<> beta((D-1)/2.0,(D-1)/2.0);
  Normal normal(0,1);
  Vector V(D-1,0);  // (D-1)-unit vector
  Vector random_vmf(D,0);

  // step 0
  long double tmp1,tmp2,p,Z,U,W,check;
  tmp1 = (4 * kappa * kappa) + (D-1) * (D-1);
  long double b = (-2 * kappa + sqrt(tmp1)) / (long double) (D-1);
  long double x0 =  (1 - b) / (1 + b);
  long double c = (kappa * x0) + ((D-1) * log(1 - x0*x0));

  while (count < sample_size) {
    // step 1
    p = uniform_random();
    Z = quantile(beta,p);
    U = uniform_random();
    //cout << "U: " << U << endl;
    tmp1 = 1 - ((1+b) * Z);
    tmp2 = 1 - ((1-b) * Z);
    W = tmp1 / tmp2;
    //cout << "W: " << W << endl;

    // step 2
    check = kappa * W + ((D-1) * log(1 - x0*W)) - c;
    if (check >= log(U)) {
      // step 3
      Vector random_normal = normal.generate(D-1);
      normalize(random_normal,V);
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
std::vector<Vector> vMF::generate(int sample_size)
{
  cout << "\nGenerating from vMF with mean: ";
  Vector spherical(3,0);
  cartesian2spherical(mu,spherical);
  spherical[1] *= 180 / PI;
  spherical[2] *= 180 / PI;
  print(cout,spherical,3);
  cout << "; Kappa: " << kappa << "; sample size = " << sample_size << endl;

  if (sample_size != 0) {
    std::vector<Vector> canonical_sample;
    generateCanonical(canonical_sample,sample_size);
    if (fabs(mu[2] - 1) <= TOLERANCE) {  // check if mu is Z-axis
      return canonical_sample;
    } else {
      Matrix transformation = align_zaxis_with_vector(mu);
      return transform(canonical_sample,transformation);
    }
  } else if (sample_size == 0) {
    return std::vector<Vector>(); 
  }
}

long double vMF::computeLogPriorProbability()
{
  long double log_prior_mean = computeLogPriorMean();
  long double log_prior_scale = computeLogPriorScale();
  long double log_joint_prior = log_prior_mean + log_prior_scale;
  assert(!boost::math::isnan(log_joint_prior));
  return log_joint_prior;
}

long double vMF::computeLogPriorMean()
{
  long double angle = theta;
  if (angle < TOLERANCE) angle = TOLERANCE;
  
  long double log_prior = 0;
  log_prior = log(sin(angle));
  log_prior -= log(4*PI);
  return log_prior;
}

long double vMF::computeLogPriorScale()
{
  long double log_prior = 0;
  log_prior = log(4/PI);
  log_prior += 2 * log(kappa);
  log_prior -= 2 * log(1+kappa*kappa);
  return log_prior;
}

/*!
 *  \brief This function computes the expected Fisher information using the
 *  component parameters.
 *  \return the expected Fisher value
 */
long double vMF::computeLogFisherInformation()
{
  long double log_fisher = 0;
  if (theta < TOLERANCE) {
    log_fisher += 2 * log(sin(TOLERANCE));
  } else {
    log_fisher += 2 * log(sin(theta));
  }
  log_fisher += log(kappa);
  long double kappa_inv = 1 / kappa;
  // A_3(k)
  long double a3k = (1 / tanh(kappa)) - kappa_inv;
  // A_3(k) -- derivative
  long double a3k_der = (kappa_inv * kappa_inv) - (1 /(sinh(kappa)*sinh(kappa)));
  log_fisher += 2 * log(fabs(a3k));
  log_fisher += log(fabs(a3k_der));
  assert(log_fisher < INFINITY);
  return log_fisher;
}

long double vMF::computeLogFisherInformation(long double N)
{
  long double log_fisher = computeLogFisherInformation();
  log_fisher += 3 * log(N);
  return log_fisher;
}

void vMF::computeAllEstimators(std::vector<Vector> &data)
{
  std::vector<struct Estimates_vMF> all_estimates;
  computeAllEstimators(data,all_estimates);
}

void vMF::computeAllEstimators(
  std::vector<Vector> &data,
  std::vector<struct Estimates_vMF> &all_estimates
) {
  int N = data.size();

  all_estimates.clear();

  Vector weights(N,1.0);
  struct Estimates_vMF mlapprox_est;
  estimateMean(mlapprox_est,data,weights);

  // ML_APPROX
  estimateMLApproxKappa(mlapprox_est);
  all_estimates.push_back(mlapprox_est);
  long double msglen = computeMessageLength(mlapprox_est);
  cout << "msglen: " << msglen << endl;
  cout << "KL-divergence: " << computeKLDivergence(mlapprox_est) << endl << endl;

  // MLE
  string type = "MLE";
  struct Estimates_vMF ml_est = mlapprox_est;
  Optimize2 opt_mle(type);
  opt_mle.initialize(N,ml_est.R,ml_est.mean,ml_est.kappa);
  opt_mle.computeEstimates(ml_est);
  print(type,ml_est);
  msglen = computeMessageLength(ml_est);
  cout << "msglen: " << msglen << endl;
  cout << "KL-divergence: " << computeKLDivergence(ml_est) << endl << endl;
  all_estimates.push_back(ml_est);

  // MAP
  type = "MAP";
  struct Estimates_vMF map_est = mlapprox_est;
  Optimize2 opt_map(type);
  opt_map.initialize(N,map_est.R,map_est.mean,map_est.kappa);
  opt_map.computeEstimates(map_est);
  print(type,map_est);
  msglen = computeMessageLength(map_est);
  cout << "msglen: " << msglen << endl;
  cout << "KL-divergence: " << computeKLDivergence(map_est) << endl << endl;
  all_estimates.push_back(map_est);

  // MML
  type = "MML_2";
  struct Estimates_vMF mml_est = mlapprox_est;
  //struct Estimates_vMF mml_est = map_est;
  Optimize2 opt_mml(type);
  opt_mml.initialize(N,mml_est.R,mml_est.mean,mml_est.kappa);
  opt_mml.computeEstimates(mml_est);
  print(type,mml_est);
  msglen = computeMessageLength(mml_est);
  cout << "msglen: " << msglen << endl;
  cout << "KL-divergence: " << computeKLDivergence(mml_est) << endl << endl;
  all_estimates.push_back(mml_est);
}

void vMF::estimateMean(
  struct Estimates_vMF &estimates, 
  std::vector<Vector> &data, 
  Vector &weights
) {
  Vector resultant = computeVectorSum(data,weights,estimates.Neff);
  estimates.mean = Vector(3,0);
  estimates.R = normalize(resultant,estimates.mean); // norm of resultant
  estimates.Rbar = estimates.R / estimates.Neff;
  if (estimates.Rbar >= 1) {
    assert((estimates.Rbar - 1) <= 0.0001);
    estimates.Rbar = 1 - TOLERANCE;
    estimates.R = estimates.Neff * estimates.Rbar;
  } 

  Vector spherical(3,0);
  cartesian2spherical(estimates.mean,spherical);
  estimates.theta = spherical[1];
  estimates.phi = spherical[2];
}

void vMF::estimateMLApproxKappa(struct Estimates_vMF &estimates)
{
  long double rbar = estimates.Rbar;
  long double num = rbar * (3 - (rbar * rbar));
  long double denom = 1 - (rbar * rbar);

  estimates.kappa = num / denom;
  cout << "Kappa (ML approx): " << estimates.kappa << endl;
}

struct Estimates_vMF vMF::computeMMLEstimates(std::vector<Vector> &data)
{
  Vector weights(data.size(),1);
  struct Estimates_vMF mlapprox_est;
  estimateMean(mlapprox_est,data,weights);
  estimateMLApproxKappa(mlapprox_est);
  return computeMMLEstimates(mlapprox_est);
}

struct Estimates_vMF vMF::computeMMLEstimates(struct Estimates_vMF &mlapprox_est)
{
  string type;
  long double msglen;

  /*type = "MAP";
  struct Estimates_vMF map_est = mlapprox_est;
  Optimize2 opt_map(type);
  opt_map.initialize(map_est.Neff,map_est.R,map_est.mean,map_est.kappa);
  opt_map.computeEstimates(map_est);
  print(type,map_est);
  msglen = computeMessageLength(map_est);
  cout << "msglen (bpr): " << msglen/map_est.Neff << endl;*/

  type = "MML_2";
  struct Estimates_vMF mml_est = mlapprox_est;
  Optimize2 opt_mml(type);
  opt_mml.initialize(mml_est.Neff,mml_est.R,mml_est.mean,mml_est.kappa);
  opt_mml.computeEstimates(mml_est);
  print(type,mml_est);
  msglen = computeMessageLength(mml_est);
  cout << "msglen (bpr): " << msglen/mml_est.Neff << endl;

  return mml_est;
}

void vMF::estimateParameters(std::vector<Vector> &data, Vector &weights)
{
  struct Estimates_vMF mlapprox_est;
  estimateMean(mlapprox_est,data,weights);
  estimateMLApproxKappa(mlapprox_est);

  struct Estimates_vMF mml_est = computeMMLEstimates(mlapprox_est);
  updateParameters(mml_est);
}

void vMF::updateParameters(struct Estimates_vMF &estimates)
{
  mu = estimates.mean;
  kappa = estimates.kappa;
  updateConstants();
}

long double vMF::computeMessageLength(std::vector<Vector> &data)
{
  Vector sample_mean = computeVectorSum(data);
  long double R = computeDotProduct(sample_mean,mu);
  return computeMessageLength(R,data.size());
}

long double vMF::computeMessageLength(long double R, long double N)
{
  long double log_prior = computeLogPriorProbability();
  long double log_fisher = computeLogFisherInformation(N);
  long double part1 = -3.816 - log_prior + 0.5 * log_fisher;
  long double part2 = computeNegativeLogLikelihood(R,N) + 1.5
                      - 2 * N * log(AOM);
  long double msglen = part1 + part2;
  return msglen/log(2);
}

long double vMF::computeMessageLength(struct Estimates_vMF &estimates)
{
  vMF vmf_est(estimates.mean,estimates.kappa);
  return vmf_est.computeMessageLength(estimates.R,estimates.Neff);
}

// 'other' is the approximate to the true distribution
long double vMF::computeKLDivergence(vMF &other)
{
  long double log_cd2 = other.getLogNormalizationConstant();

  long double ans = log_cd - log_cd2; 
  
  long double kappa2 = other.Kappa();
  Vector mu2 = other.Mean();
  long double dp = computeDotProduct(mu,mu2);
  long double diff = kappa - kappa2 * dp;
  long double ak1 = (1/tanh(kappa)) - (1/kappa);

  ans += (ak1 * diff);

  return ans;
}

long double vMF::computeKLDivergence(struct Estimates_vMF &estimates)
{
  vMF vmf_est(estimates.mean,estimates.kappa);
  return computeKLDivergence(vmf_est);
}

