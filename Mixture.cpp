#include "Mixture.h"
#include "Support.h"
#include "Optimize.h"

extern int MIXTURE_ID;
extern int MIXTURE_SIMULATION;
extern int INFER_COMPONENTS;
extern int ENABLE_DATA_PARALLELISM;
extern int NUM_THREADS;
extern int ESTIMATION,CRITERION;
extern double IMPROVEMENT_RATE;
extern int SPLITTING;
extern int IGNORE_SPLIT;
extern double MIN_N;
extern int MSGLEN_FAIL;

/*!
 *  \brief Null constructor module
 */
Mixture::Mixture()
{
  id = MIXTURE_ID++;
}

/*!
 *  \brief This is a constructor function which instantiates a Mixture
 *  \param K an integer
 *  \param components a reference to a std::vector<Kent>
 *  \param weights a reference to a Vector 
 */
Mixture::Mixture(int K, std::vector<Kent> &components, Vector &weights):
                 K(K), components(components), weights(weights)
{
  assert(components.size() == K);
  assert(weights.size() == K);
  id = MIXTURE_ID++;
  minimum_msglen = 0;
  negloglike = 0;
  aic = 0; bic = 0; icl = 0;
}

/*!
 *  \brief This is a constructor function.
 *  \param K an integer
 *  \param data a reference to a std::vector<Vector>
 *  \param data_weights a reference to a Vector
 */
Mixture::Mixture(int K, std::vector<Vector> &data, Vector &data_weights) : 
                 K(K), data(data), data_weights(data_weights)
{
  id = MIXTURE_ID++;
  N = data.size();
  assert(data_weights.size() == N);
  minimum_msglen = 0;
  negloglike = 0;
  aic = 0; bic = 0; icl = 0;
}

/*!
 *  \brief This is a constructor function.
 *  \param K an integer
 *  \param components a reference to a std::vector<Kent>
 *  \param weights a reference to a Vector
 *  \param sample_size a reference to a Vector
 *  \param responsibility a reference to a std::vector<Vector>
 *  \param data a reference to a std::vector<Vector>
 *  \param data_weights a reference to a Vector
 */
Mixture::Mixture(
  int K, 
  std::vector<Kent> &components, 
  Vector &weights,
  Vector &sample_size, 
  std::vector<Vector> &responsibility,
  std::vector<Vector> &data,
  Vector &data_weights
) : K(K), components(components), weights(weights), sample_size(sample_size),
    responsibility(responsibility), data(data), data_weights(data_weights)
{
  id = MIXTURE_ID++;
  assert(components.size() == K);
  assert(weights.size() == K);
  assert(sample_size.size() == K);
  assert(responsibility.size() == K);
  N = data.size();
  assert(data_weights.size() == N);
  minimum_msglen = 0;
  negloglike = 0;
  aic = 0; bic = 0; icl = 0;
}

/*!
 *  \brief This function assigns a source Mixture distribution.
 *  \param source a reference to a Mixture
 */
Mixture Mixture::operator=(const Mixture &source)
{
  if (this != &source) {
    id = source.id;
    N = source.N;
    K = source.K;
    components = source.components;
    data = source.data;
    data_weights = source.data_weights;
    responsibility = source.responsibility;
    sample_size = source.sample_size;
    weights = source.weights;
    msglens = source.msglens;
    null_msglen = source.null_msglen;
    minimum_msglen = source.minimum_msglen;
    part1 = source.part1;
    part2 = source.part2;
    Ik = source.Ik;
    Iw = source.Iw;
    It = source.It;
    sum_It = source.sum_It;
    Il = source.Il;
    kd_term = source.kd_term;
    negloglike = source.negloglike;
    aic = source.aic;
    bic = source.bic;
    icl = source.icl;
  }
  return *this;
}

/*!
 *  \brief This function checks whether the two Mixture objects are the same.
 *  \param other a reference to a Mixture
 *  \return whether they are the same object or not
 */
bool Mixture::operator==(const Mixture &other)
{
  if (id == other.id) {
    return 1;
  } else {
    return 0;
  }
}

/*!
 *  \brief This function returns the list of all weights.
 *  \return the list of weights
 */
Vector Mixture::getWeights()
{
  return weights;
}

/*!
 *  \brief This function returns the list of components.
 *  \return the components
 */
std::vector<Kent> Mixture::getComponents()
{
  return components;
}

/*!
 *  \brief Gets the number of components
 */
int Mixture::getNumberOfComponents()
{
  return components.size();
}

/*!
 *  \brief This function returns the responsibility matrix.
 */
std::vector<Vector> Mixture::getResponsibilityMatrix()
{
  return responsibility;
}

/*!
 *  \brief This function returns the sample size of the mixture.
 *  \return the sample size
 */
Vector Mixture::getSampleSize()
{
  return sample_size;
}

/*!
 *  \brief This function initializes the parameters of the model.
 */
void Mixture::initialize()
{
  N = data.size();
  //cout << "Sample size: " << N << endl;

  // initialize responsibility matrix
  //srand(time(NULL));
  Vector tmp(N,0);
  responsibility = std::vector<Vector>(K,tmp);
  /*for (int i=0; i<K; i++) {
    responsibility.push_back(tmp);
  }*/

  #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) 
  for (int i=0; i<N; i++) {
    int index = rand() % K;
    responsibility[index][i] = 1;
  }
  sample_size = Vector(K,0);
  updateEffectiveSampleSize();
  weights = Vector(K,0);
  if (ESTIMATION == MML) {
    updateWeights();
  } else {
    updateWeights_ML();
  }

  // initialize parameters of each component
  components = std::vector<Kent>(K);
  updateComponents();
}

void Mixture::initialize_children_1()
{
  N = data.size();

  // initialize responsibility matrix
  Vector tmp(N,0);
  responsibility = std::vector<Vector>(K,tmp);

  #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) 
  for (int i=0; i<N; i++) {
    responsibility[0][i] = uniform_random();
    responsibility[1][i] = 1 - responsibility[0][i];
  }
  sample_size = Vector(K,0);
  updateEffectiveSampleSize();
  weights = Vector(K,0);
  if (ESTIMATION == MML) {
    updateWeights();
  } else {
    updateWeights_ML();
  }

  // initialize parameters of each component
  components = std::vector<Kent>(K);
  updateComponents();
}

void Mixture::initialize_children_2()
{
  N = data.size();

  // initialize responsibility matrix
  Vector tmp(N,0);
  responsibility = std::vector<Vector>(K,tmp);

  #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) 
  for (int i=0; i<N; i++) {
    responsibility[0][i] = uniform_random();
    responsibility[1][i] = 1 - responsibility[0][i];
  }
  sample_size = Vector(K,0);
  updateEffectiveSampleSize();
  weights = Vector(K,0);
  if (ESTIMATION == MML) {
    updateWeights();
  } else {
    updateWeights_ML();
  }

  // initialize parameters of each component
  double Neff;
  Vector sample_mean = computeVectorSum(data,data_weights,Neff);
  Matrix S = computeDispersionMatrix(data,data_weights);
  Kent parent;

  struct Estimates asymptotic_est = parent.computeAsymptoticMomentEstimates(sample_mean,S,Neff);
  string type = "MOMENT";
  struct Estimates moment_est = asymptotic_est;
  Optimize opt(type);
  opt.initialize(Neff,moment_est.mean,moment_est.major_axis,moment_est.minor_axis,
                 moment_est.kappa,moment_est.beta);
  opt.computeEstimates(sample_mean,S,moment_est);

  Vector mu,rotation_axis,spherical(3,0);
  double theta,psi,alpha,eta,kappa,beta;

  rotation_axis = moment_est.minor_axis;
  kappa = moment_est.kappa;
  beta = moment_est.beta;
  //beta = 0;

  theta = 0.5 * PI * uniform_random();
  Matrix R = rotate_about_arbitrary_axis(rotation_axis,theta);
  mu = prod(R,moment_est.mean);
  cartesian2spherical(mu,spherical);
  alpha = spherical[1]; eta = spherical[2];
  psi = uniform_random() * 2 * PI;
  Kent child1(psi,alpha,eta,kappa,beta);

  theta = 0.5 * PI * uniform_random();
  R = rotate_about_arbitrary_axis(rotation_axis,-theta);
  mu = prod(R,moment_est.mean);
  cartesian2spherical(mu,spherical);
  alpha = spherical[1]; eta = spherical[2];
  psi = uniform_random() * 2 * PI;
  Kent child2(psi,alpha,eta,kappa,beta);

  components = std::vector<Kent>(K);
  components[0] = child1; components[1] = child2;
  //updateComponents();
}

void Mixture::initialize_children_3()
{
  N = data.size();

  // initialize responsibility matrix
  Vector tmp(N,0);
  responsibility = std::vector<Vector>(K,tmp);

  #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) 
  for (int i=0; i<N; i++) {
    responsibility[0][i] = uniform_random();
    responsibility[1][i] = 1 - responsibility[0][i];
  }
  sample_size = Vector(K,0);
  updateEffectiveSampleSize();
  weights = Vector(K,0);
  if (ESTIMATION == MML) {
    updateWeights();
  } else {
    updateWeights_ML();
  }

  // initialize parameters of each component
  double Neff;
  Vector sample_mean = computeVectorSum(data,data_weights,Neff);
  Matrix S = computeDispersionMatrix(data,data_weights);
  Kent parent;

  struct Estimates asymptotic_est = parent.computeAsymptoticMomentEstimates(sample_mean,S,Neff);
  string type = "MOMENT";
  struct Estimates moment_est = asymptotic_est;
  Optimize opt(type);
  opt.initialize(Neff,moment_est.mean,moment_est.major_axis,moment_est.minor_axis,
                 moment_est.kappa,moment_est.beta);
  opt.computeEstimates(sample_mean,S,moment_est);

  Vector mu,mj,mi,spherical(3,0);
  double theta,psi,alpha,eta,kappa,beta;

  mi = moment_est.minor_axis;
  kappa = moment_est.kappa;
  beta = moment_est.beta;

  //double ecc = 2 * beta / kappa;
  //theta = 0.5 * PI * ecc;
  double span = uniform_random();
  theta = 0.5 * PI * span;
  Matrix R = rotate_about_arbitrary_axis(mi,theta);
  mu = prod(R,moment_est.mean);
  mj = crossProduct(mi,mu);
  Kent child1(mu,mj,mi,kappa,beta);

  theta = 0.5 * PI * span;
  R = rotate_about_arbitrary_axis(mi,-theta);
  mu = prod(R,moment_est.mean);
  mj = crossProduct(mi,mu);
  Kent child2(mu,mj,mi,kappa,beta);

  components = std::vector<Kent>(K);
  components[0] = child1; components[1] = child2;
  //updateComponents();
}

/*!
 *  \brief This function updates the effective sample size of each component.
 */
void Mixture::updateEffectiveSampleSize()
{
  for (int i=0; i<K; i++) {
    double count = 0;
    #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) reduction(+:count)
    for (int j=0; j<N; j++) {
      count += responsibility[i][j];
    }
    sample_size[i] = count;
  }
}

/*!
 *  \brief This function is used to update the weights of the components.
 */
void Mixture::updateWeights()
{
  double normalization_constant = N + (K/2.0);
  for (int i=0; i<K; i++) {
    weights[i] = (sample_size[i] + 0.5) / normalization_constant;
  }
}

void Mixture::updateWeights_ML()
{
  double normalization_constant = N;
  for (int i=0; i<K; i++) {
    weights[i] = sample_size[i] / normalization_constant;
  }
}

/*!
 *  \brief This function is used to update the components.
 */
void Mixture::updateComponents()
{
  Vector comp_data_wts(N,0);
  for (int i=0; i<K; i++) {
    #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) 
    for (int j=0; j<N; j++) {
      comp_data_wts[j] = responsibility[i][j] * data_weights[j];
    }
    components[i].estimateParameters(data,comp_data_wts);
  }
}

/*!
 *  \brief This function updates the terms in the responsibility matrix.
 */
void Mixture::updateResponsibilityMatrix()
{
  #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) //private(j)
  for (int i=0; i<N; i++) {
    Vector log_densities(K,0);
    for (int j=0; j<K; j++) {
      log_densities[j] = components[j].log_density(data[i]);
    }
    int max_index = maximumIndex(log_densities);
    double max_log_density = log_densities[max_index];
    for (int j=0; j<K; j++) {
      log_densities[j] -= max_log_density; 
    }
    double px = 0;
    Vector probabilities(K,0);
    for (int j=0; j<K; j++) {
      probabilities[j] = weights[j] * exp(log_densities[j]);
      px += probabilities[j];
    }
    for (int j=0; j<K; j++) {
      responsibility[j][i] = probabilities[j] / px;
      assert(!boost::math::isnan(responsibility[j][i]));
      /*if(boost::math::isnan(responsibility[j][i])) {
        cout << "error: resp (pj,px): " << probabilities[j] << "; " << px << endl;
        cout << "pj: "; print(cout,probabilities); cout << endl;
        cout << "log_den: "; print(cout,log_densities); cout << endl;
        cout << "weights: "; print(cout,weights); cout << endl;
        exit(1);
      }*/
    }
  }
}

/*!
 *  \brief This function updates the terms in the responsibility matrix.
 */
void Mixture::computeResponsibilityMatrix(std::vector<Vector> &sample,
                                          string &output_file)
{
  int sample_size = sample.size();
  Vector tmp(sample_size,0);
  std::vector<Vector> resp(K,tmp);
  #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) //private(j)
  for (int i=0; i<sample_size; i++) {
    Vector log_densities(K,0);
    for (int j=0; j<K; j++) {
      log_densities[j] = components[j].log_density(sample[i]);
    }
    int max_index = maximumIndex(log_densities);
    double max_log_density = log_densities[max_index];
    for (int j=0; j<K; j++) {
      log_densities[j] -= max_log_density; 
    }
    double px = 0;
    Vector probabilities(K,0);
    for (int j=0; j<K; j++) {
      probabilities[j] = weights[j] * exp(log_densities[j]);
      px += probabilities[j];
    }
    for (int j=0; j<K; j++) {
      resp[j][i] = probabilities[j] / px;
    }
  }
  ofstream out(output_file.c_str());
  for (int i=0; i<sample_size; i++) {
    for (int j=0; j<K; j++) {
      out << fixed << setw(10) << setprecision(5) << resp[j][i];
    }
    out << endl;
  }
  out << "Cumulative memberships:\n";
  double comp_sum;
  for (int j=0; j<K; j++) {
    comp_sum = computeSum(resp[j]);
    out << "Component " << j+1 << ": " << comp_sum << endl;
  }
  out.close();
  out.close();
}

/*!
 *
 */
double Mixture::log_probability(Vector &x)
{
  Vector log_densities(K,0);
  for (int j=0; j<K; j++) {
    log_densities[j] = components[j].log_density(x);
    assert(!boost::math::isnan(log_densities[j]));
  }
  int max_index = maximumIndex(log_densities);
  double max_log_density = log_densities[max_index];
  for (int j=0; j<K; j++) {
    log_densities[j] -= max_log_density;
  }
  double density = 0;
  for (int j=0; j<K; j++) {
    density += weights[j] * exp(log_densities[j]);
  }
  return max_log_density + log(density);
}

/*!
 *  \brief This function computes the negative log likelihood of a data
 *  sample.
 *  \param a reference to a std::vector<array<double,2> >
 *  \return the negative log likelihood (base e)
 */
double Mixture::computeNegativeLogLikelihood(std::vector<Vector> &sample)
{
  double value=0,log_density;
  #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) private(log_density) reduction(-:value)
  for (int i=0; i<sample.size(); i++) {
    log_density = log_probability(sample[i]);
    if(boost::math::isnan(log_density)) {
      writeToFile("resp",responsibility,3); 
    }
    assert(!boost::math::isnan(log_density));
    value -= log_density;
  }
  negloglike = value;
  return negloglike;
}

double Mixture::computeNegativeLogLikelihood(int verbose)
{
  //return computeNegativeLogLikelihood(data);
  double neglog = computeNegativeLogLikelihood(data);
  return neglog - (2 * N * log(AOM));
}

/*!
 *  \brief This function computes the minimum message length using the current
 *  model parameters.
 *  \return the minimum message length
 */
double Mixture::computeMinimumMessageLength(int verbose /* default = 1 (print) */)
{
  MSGLEN_FAIL = 0;
  part1 = 0;
  part2 = 0;
  minimum_msglen = 0;

  /****************** PART 1 *********************/

  // encode the number of components
  // assume uniform priors
  //double Ik = log(MAX_COMPONENTS);
  Ik = K;
  //Ik = log(MAX_COMPONENTS) / log(2);
  //cout << "Ik: " << Ik << endl;

  // enocde the weights
  Iw = ((K-1)/2.0) * log(N);
  Iw -= boost::math::lgamma<double>(K); // log(K-1)!
  for (int i=0; i<K; i++) {
    Iw -= 0.5 * log(weights[i]);
  }
  Iw /= log(2);
  //cout << "Iw: " << Iw << endl;
  //assert(Iw >= 0);

  // encode the parameters of the components
  It.clear();
  sum_It = 0;
  double logp;
  for (int i=0; i<K; i++) {
    logp = components[i].computeLogParametersProbability(sample_size[i]);
    logp /= log(2);
    It.push_back(logp);
    sum_It += logp;
  }
  //cout << "It: " << sum_It << endl;
  /*if (It <= 0) { cout << It << endl;}
  fflush(stdout);
  assert(It > 0);*/

  // the constant term
  //int D = data[0].size();
  int num_free_params = (5 * K) + (K - 1);
  double log_lattice_constant = logLatticeConstant(num_free_params);
  kd_term = 0.5 * num_free_params * log_lattice_constant;
  kd_term /= log(2);

  part1 = Ik + Iw + sum_It + kd_term;

  /****************** PART 2 *********************/

  // encode the likelihood of the sample
  double Il_partial = computeNegativeLogLikelihood(data);
  Il = Il_partial - (2 * N * log(AOM));
  Il /= log(2);
  //cout << "Il: " << Il << endl;
  //assert(Il > 0);
  //if (Il < 0 || boost::math::isnan(Il)) {
  if (boost::math::isnan(Il)) {
    cout << "isnan(Il)\n"; sleep(5);
    minimum_msglen = LARGE_NUMBER;
    MSGLEN_FAIL = 1;
    return minimum_msglen;
  }

  double constant = 0.5 * num_free_params;
  constant /= log(2);
  part2 = Il + constant;

  minimum_msglen = part1 + part2;

  if (verbose == 1) {
    cout << "Ik: " << Ik << endl;
    cout << "Iw: " << Iw << endl;
    cout << "It: " << sum_It << endl;
    cout << "Il: " << Il << endl;
  }

  return minimum_msglen;
}

void Mixture::printIndividualMsgLengths(ostream &log_file)
{
  log_file << "\t\tIk: " << Ik << endl;
  log_file << "\t\tIw: " << Iw << endl;
  log_file << "\t\tIt: " << sum_It << " "; print(log_file,It,3); log_file << endl;
  log_file << "\t\tlatt: " << kd_term << endl;
  log_file << "\t\tIl: " << Il << endl;
  log_file << "\t\tpart1 (Ik+Iw+It+latt): " << part1 << " + " 
           << "part2 (Il+d/(2*log(2))): " << part2 << " = "
           << part1 + part2 << " bits." << endl << endl;
}

/*!
 *  \brief Prepares the appropriate log file
 */
string Mixture::getLogFile()
{
  string file_name;
  if (INFER_COMPONENTS == UNSET) {
    if (MIXTURE_SIMULATION == UNSET) {
      file_name = "./mixture/logs/";
    } else if (MIXTURE_SIMULATION == SET) {
      file_name = "./simulation/logs/";
    }
  } else if (INFER_COMPONENTS == SET) {
    file_name = "./infer/logs/kent/";
    file_name += "m_" + boost::lexical_cast<string>(id) + "_";
  }
  file_name += boost::lexical_cast<string>(K) + ".log";
  return file_name;
}

/*!
 *  \brief This function is used to estimate the model parameters by running
 *  an EM algorithm.
 *  \return the stable message length
 */
double Mixture::estimateParameters()
{
  if (SPLITTING == 1) {
    //initialize_children_1();
    //initialize_children_2();
    initialize_children_3();
  } else {
    initialize();
  }

  //initialize();

  EM();

  return minimum_msglen;
}

/*!
 *  \brief This function runs the EM method.
 */
void Mixture::EM()
{
  /* prepare log file */
  string log_file = getLogFile();
  ofstream log(log_file.c_str());

  computeNullModelMessageLength();
  //cout << "null_msglen: " << null_msglen << endl;

  printParameters(log,0,0);

  if (ESTIMATION == MML) {
    EM(
      log,
      &Mixture::updateWeights,
      &Mixture::computeMinimumMessageLength
    );
  } else {
    EM(
      log,
      &Mixture::updateWeights_ML,
      &Mixture::computeNegativeLogLikelihood
    );
  }

  switch(CRITERION) {
    case AIC:
      aic = computeAIC();
      break;

    case BIC:
      bic = computeBIC();
      break;

    case ICL:
      icl = computeICL();
      break;

    case MMLC:
      break;
  }

  log.close();
}

void Mixture::EM(
  ostream &log,
  void (Mixture::*update_weights)(),
  double (Mixture::*objective_function)(int)
) {
  double prev=0,current;
  int iter = 1;
  double impr_rate = 0.00001;
  int MIN_ITER = 5;

  while (1) {
    // Expectation (E-step)
    updateResponsibilityMatrix();
    updateEffectiveSampleSize();
    //if (SPLITTING == 1) {
      for (int i=0; i<K; i++) {
        if (sample_size[i] < MIN_N) {
          current = (this->*objective_function)(0);
          cout << "stopping \n";
          goto stop;
        }
      }
    //}
    // Maximization (M-step)
    (this->*update_weights)();

    updateComponents();
    current = (this->*objective_function)(0);
    printParameters(log,iter,current);
    if (iter != 1) {
      //assert(current > 0);
      //if ((iter > 3 && (prev - current) <= impr_rate * prev) ||
      //      (iter > 1 && current > prev) || current <= 0 || MSGLEN_FAIL == 1) {
      if (
          (iter > MIN_ITER && current >= prev) ||
          MSGLEN_FAIL == 1 ||
          (iter > MIN_ITER && (fabs(prev - current) <= impr_rate * fabs(prev)))
         ) {
        stop:
        current = computeMinimumMessageLength();
        log << "\nSample size: " << N << endl;
        log << "Kent encoding rate: " << current << " bits.";
        log << "\t(" << current/N << " bits/point)" << endl;
        log << "Null model encoding: " << null_msglen << " bits.";
        log << "\t(" << null_msglen/N << " bits/point)" << endl;
        break;
      }
    }
    prev = current;
    iter++;
  }
}

/*!
 *  \brief This function computes the null model message length.
 *  \return the null model message length
 */
double Mixture::computeNullModelMessageLength()
{
  // compute logarithm of surface area of nd-sphere
  double log_area = log(4*PI);
  null_msglen = N * (log_area - (2*log(AOM)));
  null_msglen /= log(2);
  return null_msglen;
}

/*!
 *  \brief This function returns the minimum message length of this mixture
 *  model.
 */
double Mixture::getMinimumMessageLength()
{
  return minimum_msglen;
}

/*!
 *  \brief Returns the first part of the msg.
 */
double Mixture::first_part()
{
  return part1;
}

/*!
 *  \brief Returns the second part of the msg.
 */
double Mixture::second_part()
{
  return part2;
}

double Mixture::getNegativeLogLikelihood()
{
  return negloglike;
}

double Mixture::getAIC()
{
  return aic;
}

double Mixture::getBIC()
{
  return bic;
}

double Mixture::getICL()
{
  return icl;
}

/*!
 *  \brief This function prints the parameters to a log file.
 *  \param os a reference to a ostream
 *  \param iter an integer
 *  \param msglen a double
 */
void Mixture::printParameters(ostream &os, int iter, double value)
{
  os << "Iteration #: " << iter << endl;
  for (int k=0; k<K; k++) {
    os << "\t" << fixed << setw(5) << "[" << k+1 << "]";
    os << "\t" << fixed << setw(10) << setprecision(3) << sample_size[k];
    os << "\t" << fixed << setw(10) << setprecision(5) << weights[k];
    os << "\t";
    components[k].printParameters(os);
  }
  if (ESTIMATION == MML) {
    os << "\t\t\tmsglen: " << value << " bits." << endl;
  } else {
    os << "\t\t\tnegloglike: " << value/log(2) << " bits." << endl;
  }
}

/*!
 *  \brief This function prints the parameters to a log file.
 *  \param os a reference to a ostream
 */
void Mixture::printParameters(ostream &os, int num_tabs)
{
  string tabs = "\t";
  if (num_tabs == 2) {
    tabs += "\t";
  }
  for (int k=0; k<K; k++) {
    os << tabs << "[";// << fixed << setw(5) << "[" << k+1 << "]";
    os << fixed << setw(2) << k+1;
    os << "]";
    os << "\t" << fixed << setw(10) << setprecision(3) << sample_size[k];
    os << "\t" << fixed << setw(10) << setprecision(5) << weights[k];
    os << "\t";
    components[k].printParameters(os);
  }
  os << tabs << "ID: " << id << endl;
  os << tabs << "Kent encoding: " << minimum_msglen << " bits. "
     << "(" << minimum_msglen/N << " bits/point)" << endl;

  switch(CRITERION) {
    case AIC:
      aic = computeAIC();
      os << tabs << "AIC: " << aic << endl;
      break;

    case BIC:
      bic = computeBIC();
      os << tabs << "BIC: " << bic << endl;
      break;

    case ICL:
      icl = computeICL();
      os << tabs << "ICL: " << icl << endl;
      break;

    case MMLC:
      break;
  }
  os << endl;
}

void Mixture::printParameters(string &file)
{
  ofstream out(file.c_str());
  printParameters(out);
  out.close();
}

/*!
 *  \brief Outputs the mixture weights and component parameters to a file.
 */
void Mixture::printParameters(ostream &os)
{
  for (int k=0; k<K; k++) {
    os << "\t" << fixed << setw(10) << setprecision(5) << weights[k];
    os << "\t";
    components[k].printParameters(os);
  }
}

/*!
 *  \brief This function is used to read the mixture details to aid in
 *  visualization.
 *  \param file_name a reference to a string
 *  \param D an integer
 */
void Mixture::load(string &file_name)
{
  sample_size.clear();
  weights.clear();
  components.clear();
  K = 0;
  ifstream file(file_name.c_str());
  string line;
  Vector numbers;
  Vector unit_mean(3,0),mean(3,0);
  Vector unit_mj(3,0),mj(3,0),unit_mi(3,0),mi(3,0);
  double sum_weights = 0;
  while (getline(file,line)) {
    K++;
    boost::char_separator<char> sep("mujikapbet,:()[] \t");
    boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
    BOOST_FOREACH (const string& t, tokens) {
      istringstream iss(t);
      double x;
      iss >> x;
      numbers.push_back(x);
    }
    weights.push_back(numbers[0]);
    sum_weights += numbers[0];
    for (int i=1; i<=3; i++) {
      mean[i-1] = numbers[i];
      mj[i-1] = numbers[i+3];
      mi[i-1] = numbers[i+6];
    }
    double kappa = numbers[10];
    double beta = numbers[11];
    double ex = 2 * beta / kappa;
    if (fabs(ex - 1) <= TOLERANCE) {
      beta -= TOLERANCE;
    }
    normalize(mean,unit_mean);
    normalize(mj,unit_mj); normalize(mi,unit_mi);
    Kent kent(unit_mean,unit_mj,unit_mi,kappa,beta);
    components.push_back(kent);
    numbers.clear();
  }
  file.close();
  for (int i=0; i<K; i++) {
    weights[i] /= sum_weights;
  }
}

/*!
 *  \brief This function is used to read the mixture details 
 *  corresponding to the given data.
 *  \param file_name a reference to a string
 *  \param D an integer
 *  \param d a reference to a std::vector<Vector>
 *  \param dw a reference to a Vector
 */
void Mixture::load(string &file_name, std::vector<Vector> &d, Vector &dw)
{
  load(file_name);
  data = d;
  N = data.size();
  data_weights = dw;
  Vector tmp(N,0);
  responsibility = std::vector<Vector>(K,tmp);
  /*for (int i=0; i<K; i++) {
    responsibility.push_back(tmp);
  }*/
  updateResponsibilityMatrix();
  sample_size = Vector(K,0);
  updateEffectiveSampleSize();
  updateComponents();
  minimum_msglen = computeMinimumMessageLength();
}

/*!
 *  \brief This function is used to randomly choose a component.
 *  \return the component index
 */
int Mixture::randomComponent()
{
  double random = uniform_random();
  //cout << random << endl;
  double previous = 0;
  for (int i=0; i<weights.size(); i++) {
    if (random <= weights[i] + previous) {
      return i;
    }
    previous += weights[i];
  }
}

/*!
 *  \brief This function saves the data generated from a component to a file.
 *  \param index an integer
 *  \param data a reference to a std::vector<Vector>
 */
void Mixture::saveComponentData(int index, std::vector<Vector> &data)
{
  string data_file = "./visualize/sampled_data/comp";
  data_file += boost::lexical_cast<string>(index+1) + ".dat";
  ofstream file(data_file.c_str());

  Vector projection(2,0);
  double theta,phi,rho,z1,z2;
  string transformed_file = "./visualize/sampled_data/transformed_comp" 
                            + boost::lexical_cast<string>(index+1) + ".dat";
  ofstream transformed_comp(transformed_file.c_str());

  for (int j=0; j<data.size(); j++) {
    for (int k=0; k<3; k++) {
      file << fixed << setw(10) << setprecision(3) << data[j][k];
    }
    file << endl;
    computeLambertProjection(data[j],projection);
    transformed_comp << fixed << setw(10) << setprecision(3) << projection[0];
    transformed_comp << fixed << setw(10) << setprecision(3) << projection[1] << endl;
  }
  file.close();
  transformed_comp.close();
}

/*!
 *  \brief This function is used to randomly sample from the mixture
 *  distribution.
 *  \param num_samples an integer
 *  \param save_data a boolean variable
 *  \return the random sample
 */
std::vector<Vector> Mixture::generate(int num_samples, bool save_data) 
{
  sample_size = Vector(K,0);
  for (int i=0; i<num_samples; i++) {
    // randomly choose a component
    int k = randomComponent();
    sample_size[k]++;
  }
  /*ofstream fw("sample_size");
  for (int i=0; i<sample_size.size(); i++) {
    fw << sample_size[i] << endl;
  }
  fw.close();*/

  std::vector<std::vector<Vector> > random_data;
  std::vector<Vector> sample;
  for (int i=0; i<K; i++) {
    std::vector<Vector> x = components[i].generate((int)sample_size[i]);
    random_data.push_back(x);
    for (int j=0; j<random_data[i].size(); j++) {
      sample.push_back(random_data[i][j]);
    }
  } // for i

  if (save_data) {
    writeToFile("random_sample.dat",sample);
    string comp_density_file,transformed_file;
    string mix_density_file = "./visualize/sampled_data/mixture_density.dat";
    ofstream mix(mix_density_file.c_str());
    double comp_density,mix_density;
    for (int i=0; i<K; i++) {
      saveComponentData(i,random_data[i]);
      comp_density_file = "./visualize/sampled_data/comp" 
                          + boost::lexical_cast<string>(i+1) + "_density.dat";
      ofstream comp(comp_density_file.c_str());
      for (int j=0; j<random_data[i].size(); j++) {
        comp_density = exp(components[i].log_density(random_data[i][j]));
        mix_density = exp(log_probability(random_data[i][j]));
        for (int k=0; k<random_data[i][j].size(); k++) {
          comp << fixed << setw(10) << setprecision(3) << random_data[i][j][k];
          mix << fixed << setw(10) << setprecision(3) << random_data[i][j][k];
        } // k
        comp << "\t\t" << scientific << comp_density << endl;
        mix <<  "\t\t" << scientific << mix_density << endl;
      } // j
      comp.close();
    } // i
    mix.close();
    generateHeatmapData(1);
  } // if()
  return sample;
}

/*!
 *  \brief This function splits a component into two.
 *  \return c an integer
 *  \param log a reference to a ostream
 *  \return the modified Mixture
 */
Mixture Mixture::split(int c, ostream &log)
{
  SPLITTING = 1;
  log << "\tSPLIT component " << c + 1 << " ... " << endl;

  int num_children = 2; 
  Mixture m(num_children,data,responsibility[c]);
  m.estimateParameters();
  log << "\t\tChildren:\n";
  m.printParameters(log,2); // print the child mixture

  // adjust weights
  Vector weights_c = m.getWeights();
  weights_c[0] *= weights[c];
  weights_c[1] *= weights[c];

  // adjust responsibility matrix
  std::vector<Vector> responsibility_c = m.getResponsibilityMatrix();
  for (int i=0; i<2; i++) {
    #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) 
    for (int j=0; j<N; j++) {
      responsibility_c[i][j] *= responsibility[c][j];
    }
  }

  // adjust effective sample size
  Vector sample_size_c(2,0);
  for (int i=0; i<2; i++) {
    double sum = 0;
    #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) reduction(+:sum) 
    for (int j=0; j<N; j++) {
      sum += responsibility_c[i][j];
    }
    sample_size_c[i] = sum;
    if (sample_size_c[i] < MIN_N) {
      IGNORE_SPLIT = 1;
    }
  }

  // child components
  std::vector<Kent> components_c = m.getComponents();

  // merge with the remaining components
  int K_m = K + 1;
  std::vector<Vector> responsibility_m(K_m);
  Vector weights_m(K_m,0),sample_size_m(K_m,0);
  std::vector<Kent> components_m(K_m);
  int index = 0;
  for (int i=0; i<K; i++) {
    if (i != c) {
      weights_m[index] = weights[i];
      sample_size_m[index] = sample_size[i];
      responsibility_m[index] = responsibility[i];
      components_m[index] = components[i];
      index++;
    } else if (i == c) {
      for (int j=0; j<2; j++) {
        weights_m[index] = weights_c[j];
        sample_size_m[index] = sample_size_c[j];
        responsibility_m[index] = responsibility_c[j];
        components_m[index] = components_c[j];
        index++;
      }
    }
  }

  Vector data_weights_m(N,1);
  Mixture merged(K_m,components_m,weights_m,sample_size_m,responsibility_m,data,data_weights_m);
  log << "\t\tBefore adjustment ...\n";
  merged.computeMinimumMessageLength();
  merged.printParameters(log,2);
  merged.EM();
  log << "\t\tAfter adjustment ...\n";
  merged.computeMinimumMessageLength();
  merged.printParameters(log,2);
  merged.printIndividualMsgLengths(log);
  SPLITTING = 0;
  return merged;
}

/*!
 *  \brief This function deletes a component.
 *  \return c an integer
 *  \param log a reference to a ostream
 *  \return the modified Mixture
 */
Mixture Mixture::kill(int c, ostream &log)
{
  log << "\tKILL component " << c + 1 << " ... " << endl;

  int K_m = K - 1;
  // adjust weights
  Vector weights_m(K_m,0);
  double residual_sum = 0;
  for (int i=0; i<K; i++) {
    if (i != c) {
      residual_sum += weights[i];
    }
  }
  double wt;
  int index = 0;
  for (int i=0; i<K; i++) {
    if (i != c) {
      weights_m[index++] = weights[i] / residual_sum;
    }
  }

  // adjust responsibility matrix
  Vector residual_sums(N,0);
  for (int i=0; i<N; i++) {
    for (int j=0; j<K; j++) {
      if (j != c) {
        residual_sums[i] += responsibility[j][i];
      }
    }
    //if (residual_sums[i] < TOLERANCE) residual_sums[i] = TOLERANCE;
  }
  Vector resp(N,0);
  std::vector<Vector> responsibility_m(K_m,resp);
  index = 0;
  for (int i=0; i<K; i++) {
    if (i != c) {
      #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) //private(residual_sum) 
      for (int j=0; j<N; j++) {
        //residual_sum = 1 - responsibility[c][j];
        if (residual_sums[j] <= 0.06) {
          responsibility_m[index][j] = 1.0 / K_m;
        } else {
          responsibility_m[index][j] = responsibility[i][j] / residual_sums[j];
        }
        assert(responsibility_m[index][j] >= 0 && responsibility_m[index][j] <= 1);
      }
      index++;
    }
  }

  // adjust effective sample size
  Vector sample_size_m(K_m,0);
  for (int i=0; i<K-1; i++) {
    double sum = 0;
    #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) reduction(+:sum) 
    for (int j=0; j<N; j++) {
      sum += responsibility_m[i][j];
    }
    sample_size_m[i] = sum;
  }

  // child components
  std::vector<Kent> components_m(K_m);
  index = 0;
  for (int i=0; i<K; i++) {
    if (i != c) {
      components_m[index++] = components[i];
    }
  }

  log << "\t\tResidual:\n";
  Vector data_weights_m(N,1);
  Mixture modified(K_m,components_m,weights_m,sample_size_m,responsibility_m,data,data_weights_m);
  log << "\t\tBefore adjustment ...\n";
  modified.computeMinimumMessageLength();
  modified.printParameters(log,2);
  modified.EM();
  log << "\t\tAfter adjustment ...\n";
  //modified.computeMinimumMessageLength();
  modified.printParameters(log,2);
  modified.printIndividualMsgLengths(log);
  return modified;
}

/*!
 *  \brief This function joins two components.
 *  \return c1 an integer
 *  \return c2 an integer
 *  \param log a reference to a ostream
 *  \return the modified Mixture
 */
Mixture Mixture::join(int c1, int c2, ostream &log)
{
  log << "\tJOIN components " << c1+1 << " and " << c2+1 << " ... " << endl;

  int K_m = K - 1;
  // adjust weights
  Vector weights_m(K_m,0);
  int index = 0;
  for (int i=0; i<K; i++) {
    if (i != c1 && i != c2) {
      weights_m[index++] = weights[i];
    }
  }
  weights_m[index] = weights[c1] + weights[c2];

  // adjust responsibility matrix
  std::vector<Vector> responsibility_m(K_m);
  index = 0;
  for (int i=0; i<K; i++) {
    if (i != c1 && i != c2) {
      responsibility_m[index++] = responsibility[i];
    }
  }
  Vector resp(N,0);
  #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) 
  for (int i=0; i<N; i++) {
    resp[i] = responsibility[c1][i] + responsibility[c2][i];
  }
  responsibility_m[index] = resp;

  // adjust effective sample size 
  Vector sample_size_m(K_m,0);
  index = 0;
  for (int i=0; i<K; i++) {
    if (i != c1 && i != c2) {
      sample_size_m[index++] = sample_size[i];
    }
  }
  sample_size_m[index] = sample_size[c1] + sample_size[c2];

  // child components
  std::vector<Kent> components_m(K_m);
  index = 0;
  for (int i=0; i<K; i++) {
    if (i != c1 && i != c2) {
      components_m[index++] = components[i];
    }
  }
  Mixture joined(1,data,resp);
  joined.estimateParameters();
  log << "\t\tResultant join:\n";
  joined.printParameters(log,2); // print the joined pair mixture

  std::vector<Kent> joined_comp = joined.getComponents();
  components_m[index++] = joined_comp[0];
  Vector data_weights_m(N,1);
  Mixture modified(K-1,components_m,weights_m,sample_size_m,responsibility_m,data,data_weights_m);
  log << "\t\tBefore adjustment ...\n";
  modified.computeMinimumMessageLength();
  modified.printParameters(log,2);
  modified.EM();
  log << "\t\tAfter adjustment ...\n";
  modified.printParameters(log,2);
  return modified;
}

/*!
 *  \brief This function generates data to visualize the 2D/3D heat maps.
 *  \param res a double
 */
void Mixture::generateHeatmapData(double res)
{
  string data_fbins2D = "./visualize/sampled_data/prob_bins2D.dat";
  string data_fbins3D = "./visualize/sampled_data/prob_bins3D.dat";
  ofstream fbins2D(data_fbins2D.c_str());
  ofstream fbins3D(data_fbins3D.c_str());
  Vector x(3,1);
  Vector point(3,0);
  for (double theta=0; theta<180; theta+=res) {
    x[1] = theta * PI/180;
    for (double phi=0; phi<360; phi+=res) {
      x[2] = phi * PI/180;
      spherical2cartesian(x,point);
      double pr = exp(log_probability(point));
      // 2D bins
      fbins2D << fixed << setw(10) << setprecision(4) << floor(pr * 100);
      // 3D bins
      for (int k=0; k<3; k++) {
        fbins3D << fixed << setw(10) << setprecision(4) << point[k];
      }
      fbins3D << fixed << setw(10) << setprecision(4) << pr << endl;
    }
    fbins2D << endl;
  }
  fbins2D.close();
  fbins3D.close();
}

/*!
 *  \brief This function computes the nearest component to a given component.
 *  \param c an integer
 *  \return the index of the closest component
 */
int Mixture::getNearestComponent(int c)
{
  double current,dist = LARGE_NUMBER;
  int nearest;

  for (int i=0; i<K; i++) {
    if (i != c) {
      current = components[c].computeKLDivergence(components[i]);
      if (current < dist) {
        dist = current;
        nearest = i;
      }
    }
  }
  return nearest;
}

double Mixture::computeKLDivergence(Mixture &other)
{
  return computeKLDivergence(other,data);
}

// Monte Carlo
double Mixture::computeKLDivergence(Mixture &other, std::vector<Vector> &sample)
{
  double kldiv = 0,log_fx,log_gx;
  for (int i=0; i<sample.size(); i++) {
    log_fx = log_probability(sample[i]);
    log_gx = other.log_probability(sample[i]);
    kldiv += (log_fx - log_gx);
  }
  return kldiv/(log(2) * sample.size());
}

/*!
 *  \brief This function computes the Akaike information criteria (AIC)
 *  \return the AIC value (natural log -- nits)
 */
double Mixture::computeAIC()
{
  int k = 6 * K - 1;
  negloglike = computeNegativeLogLikelihood(data);
  return compute_aic(k,N,negloglike);
}

/*!
 *  \brief This function computes the Bayesian information criteria (AIC)
 *  \return the BIC value (natural log -- nits)
 */
double Mixture::computeBIC()
{
  int k = 6 * K - 1;
  negloglike = computeNegativeLogLikelihood(data);
  return compute_bic(k,N,negloglike);
}

std::vector<std::vector<int> > Mixture::compute_cluster_indicators()
{
  std::vector<int> emptyvec(N,0);
  std::vector<std::vector<int> > z(K,emptyvec);
  std::vector<Vector> flipped_resp = flip(responsibility);
  int max_index;
  for (int i=0; i<N; i++) {
    max_index = maximumIndex(flipped_resp[i]);
    z[max_index][i] = 1;
  }
  return z;
}

double Mixture::computeICL()
{
  negloglike = computeNegativeLogLikelihood(data);
  std::vector<std::vector<int> > indicators = compute_cluster_indicators();
  double ec = 0,term;
  for (int i=0; i<K; i++) {
    for (int j=0; j<N; j++) {
      term = indicators[i][j] * log(responsibility[i][j]);
      ec -= term;
    }
  }
  return (negloglike + ec) / log(2);
}

