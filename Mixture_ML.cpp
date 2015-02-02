#include "Mixture_ML.h"
#include "Support.h"

extern int MIXTURE_ID;
extern int MIXTURE_SIMULATION;
extern int INFER_COMPONENTS;
extern int ENABLE_DATA_PARALLELISM;
extern int NUM_THREADS;
extern int ESTIMATION;
extern double IMPROVEMENT_RATE;
extern int SPLITTING;
extern int IGNORE_SPLIT;
extern double MIN_N;
extern int MSGLEN_FAIL;

/*!
 *  \brief Null constructor module
 */
Mixture_ML::Mixture_ML()
{
  id = MIXTURE_ID++;
}

/*!
 *  \brief This is a constructor function which instantiates a Mixture_ML
 *  \param K an integer
 *  \param components a reference to a std::vector<Kent>
 *  \param weights a reference to a Vector 
 */
Mixture_ML::Mixture_ML(int K, std::vector<Kent> &components, Vector &weights):
                 K(K), components(components), weights(weights)
{
  assert(components.size() == K);
  assert(weights.size() == K);
  id = MIXTURE_ID++;
  negloglike = 0;
}

/*!
 *  \brief This is a constructor function.
 *  \param K an integer
 *  \param data a reference to a std::vector<Vector>
 *  \param data_weights a reference to a Vector
 */
Mixture_ML::Mixture_ML(int K, std::vector<Vector> &data, Vector &data_weights) : 
                 K(K), data(data), data_weights(data_weights)
{
  id = MIXTURE_ID++;
  N = data.size();
  assert(data_weights.size() == N);
  negloglike = 0;
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
Mixture_ML::Mixture_ML(
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
  negloglike = 0;
}

/*!
 *  \brief This function assigns a source Mixture_ML distribution.
 *  \param source a reference to a Mixture_ML
 */
Mixture_ML Mixture_ML::operator=(const Mixture_ML &source)
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
  }
  return *this;
}

/*!
 *  \brief This function checks whether the two Mixture_ML objects are the same.
 *  \param other a reference to a Mixture_ML
 *  \return whether they are the same object or not
 */
bool Mixture_ML::operator==(const Mixture_ML &other)
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
Vector Mixture_ML::getWeights()
{
  return weights;
}

/*!
 *  \brief This function returns the list of components.
 *  \return the components
 */
std::vector<Kent> Mixture_ML::getComponents()
{
  return components;
}

/*!
 *  \brief Gets the number of components
 */
int Mixture_ML::getNumberOfComponents()
{
  return components.size();
}

/*!
 *  \brief This function returns the responsibility matrix.
 */
std::vector<Vector> Mixture_ML::getResponsibilityMatrix()
{
  return responsibility;
}

/*!
 *  \brief This function returns the sample size of the mixture.
 *  \return the sample size
 */
Vector Mixture_ML::getSampleSize()
{
  return sample_size;
}

/*!
 *  \brief This function initializes the parameters of the model.
 */
void Mixture_ML::initialize()
{
  N = data.size();
  cout << "Sample size: " << N << endl;

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
  updateWeights_ML();

  // initialize parameters of each component
  components = std::vector<Kent>(K);
  updateComponents();
}

/*!
 *  \brief This function updates the effective sample size of each component.
 */
void Mixture_ML::updateEffectiveSampleSize()
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

void Mixture_ML::updateWeights_ML()
{
  double normalization_constant = N;
  for (int i=0; i<K; i++) {
    weights[i] = sample_size[i] / normalization_constant;
  }
}

/*!
 *  \brief This function is used to update the components.
 */
void Mixture_ML::updateComponents()
{
  Vector comp_data_wts(N,0);
  for (int i=0; i<K; i++) {
    #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) 
    for (int j=0; j<N; j++) {
      comp_data_wts[j] = responsibility[i][j] * data_weights[j];
    }
    components[i].estimateParameters(data,comp_data_wts);
    //components[i].updateParameters();
  }
}

/*!
 *  \brief This function updates the terms in the responsibility matrix.
 */
void Mixture_ML::updateResponsibilityMatrix()
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
void Mixture_ML::computeResponsibilityMatrix(std::vector<Vector> &sample,
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
double Mixture_ML::log_probability(Vector &x)
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
double Mixture_ML::computeNegativeLogLikelihood(std::vector<Vector> &sample)
{
  double value = 0,log_density;
  #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) private(log_density) reduction(-:value)
  for (int i=0; i<sample.size(); i++) {
    log_density = log_probability(sample[i]);
    if(boost::math::isnan(log_density)) {
      writeToFile("resp",responsibility,3); 
    }
    assert(!boost::math::isnan(log_density));
    value -= log_density;
  }
  return value;
}

/*!
 *  \brief Prepares the appropriate log file
 */
string Mixture_ML::getLogFile()
{
  string file_name;
  if (INFER_COMPONENTS == UNSET) {
    if (MIXTURE_SIMULATION == UNSET) {
      file_name = "./mixture/logs/";
    } else if (MIXTURE_SIMULATION == SET) {
      file_name = "./simulation/logs/";
    }
  } else if (INFER_COMPONENTS == SET) {
    file_name = "./infer/logs/";
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
double Mixture_ML::estimateParameters()
{
  /*if (SPLITTING == 1) {
    initialize_children_1();
  } else {
    initialize();
  }*/

  initialize();

  EM();

  return negloglike;
}

/*!
 *  \brief This function runs the EM method.
 */
void Mixture_ML::EM()
{
  /* prepare log file */
  string log_file = getLogFile();
  ofstream log(log_file.c_str());

  double prev=0,current;
  int iter = 1;
  printParameters(log,0,0);

  double impr_rate = 0.00001;
  /* EM loop */
  while (1) {
    // Expectation (E-step)
    updateResponsibilityMatrix();
    updateEffectiveSampleSize();
    //if (SPLITTING == 1) {
      for (int i=0; i<K; i++) {
        if (sample_size[i] < MIN_N) {
          current = computeNegativeLogLikelihood(data);
          cout << "stopping 2\n";
          goto stop2;
        }
      }
    //}
    // Maximization (M-step)
    updateWeights_ML();
    updateComponents();
    //current = negativeLogLikelihood(data);
    current = computeNegativeLogLikelihood(data);
    printParameters(log,iter,current);
    if (iter != 1) {
      //assert(current > 0);
      // because EM has to consistently produce lower 
      // -ve likelihood values otherwise something wrong!
      if ((iter > 3 && (prev - current) <= impr_rate * prev) ||
            (iter > 1 && current > prev) || current <= 0) {
        stop2:
        log << "\nSample size: " << N << endl;
        log << "Kent negloglikelihood: " << current/(std::log(2)*N) << " bits/point" << endl;
        break;
      }
    }
    prev = current;
    iter++;
  }
  log.close();
}


/*!
 *  \brief This function returns the minimum message length of this mixture
 *  model.
 */
double Mixture_ML::getNegativeLogLikelihood()
{
  return negloglike;
}

/*!
 *  \brief This function prints the parameters to a log file.
 *  \param os a reference to a ostream
 *  \param iter an integer
 *  \param msglen a double
 */
void Mixture_ML::printParameters(ostream &os, int iter, double neglog)
{
  os << "Iteration #: " << iter << endl;
  for (int k=0; k<K; k++) {
    os << "\t" << fixed << setw(5) << "[" << k+1 << "]";
    os << "\t" << fixed << setw(10) << setprecision(3) << sample_size[k];
    os << "\t" << fixed << setw(10) << setprecision(5) << weights[k];
    os << "\t";
    components[k].printParameters(os);
  }
  os << "\t\t\tneglog: " << neglog/log(2) << " bits." << endl;
}

/*!
 *  \brief This function prints the parameters to a log file.
 *  \param os a reference to a ostream
 */
void Mixture_ML::printParameters(ostream &os, int num_tabs)
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
  os << tabs << "Kent encoding: " << negloglike/log(2) << " bits. "
     << "(" << negloglike/(N*log(2)) << " bits/point)" << endl << endl;
}

/*!
 *  \brief Outputs the mixture weights and component parameters to a file.
 */
void Mixture_ML::printParameters(ostream &os)
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
void Mixture_ML::load(string &file_name)
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
void Mixture_ML::load(string &file_name, std::vector<Vector> &d, Vector &dw)
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
  negloglike = computeNegativeLogLikelihood(data);
}

/*!
 *  \brief This function is used to randomly choose a component.
 *  \return the component index
 */
int Mixture_ML::randomComponent()
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
void Mixture_ML::saveComponentData(int index, std::vector<Vector> &data)
{
  string data_file = "./visualize/sampled_data/comp";
  data_file += boost::lexical_cast<string>(index+1) + ".dat";
  //components[index].printParameters(cout);
  ofstream file(data_file.c_str());
  for (int j=0; j<data.size(); j++) {
    for (int k=0; k<3; k++) {
      file << fixed << setw(10) << setprecision(3) << data[j][k];
    }
    file << endl;
  }
  file.close();
}

/*!
 *  \brief This function is used to randomly sample from the mixture
 *  distribution.
 *  \param num_samples an integer
 *  \param save_data a boolean variable
 *  \return the random sample
 */
std::vector<Vector> Mixture_ML::generate(int num_samples, bool save_data) 
{
  sample_size = Vector(K,0);
  for (int i=0; i<num_samples; i++) {
    // randomly choose a component
    int k = randomComponent();
    sample_size[k]++;
  }
  ofstream fw("sample_size");
  for (int i=0; i<sample_size.size(); i++) {
    fw << sample_size[i] << endl;
  }
  fw.close();

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
    string comp_density_file;
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
  } // if()
  return sample;
}

/*!
 *  \brief This function splits a component into two.
 *  \return c an integer
 *  \param log a reference to a ostream
 *  \return the modified Mixture_ML
 */
Mixture_ML Mixture_ML::split(int c, ostream &log)
{
  SPLITTING = 1;
  log << "\tSPLIT component " << c + 1 << " ... " << endl;

  int num_children = 2; 
  Mixture_ML m(num_children,data,responsibility[c]);
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
    /*if (sample_size_c[i] < MIN_N) {
      IGNORE_SPLIT = 1;
    }*/
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
  Mixture_ML merged(K_m,components_m,weights_m,sample_size_m,responsibility_m,data,data_weights_m);
  log << "\t\tBefore adjustment ...\n";
  merged.computeNegativeLogLikelihood(data);
  merged.printParameters(log,2);
  merged.EM();
  log << "\t\tAfter adjustment ...\n";
  merged.computeNegativeLogLikelihood(data);
  merged.printParameters(log,2);
  SPLITTING = 0;
  return merged;
}

/*!
 *  \brief This function deletes a component.
 *  \return c an integer
 *  \param log a reference to a ostream
 *  \return the modified Mixture_ML
 */
Mixture_ML Mixture_ML::kill(int c, ostream &log)
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
  Mixture_ML modified(K_m,components_m,weights_m,sample_size_m,responsibility_m,data,data_weights_m);
  log << "\t\tBefore adjustment ...\n";
  modified.computeNegativeLogLikelihood(data);
  modified.printParameters(log,2);
  modified.EM();
  log << "\t\tAfter adjustment ...\n";
  modified.printParameters(log,2);
  return modified;
}

/*!
 *  \brief This function joins two components.
 *  \return c1 an integer
 *  \return c2 an integer
 *  \param log a reference to a ostream
 *  \return the modified Mixture_ML
 */
Mixture_ML Mixture_ML::join(int c1, int c2, ostream &log)
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
  Mixture_ML joined(1,data,resp);
  joined.estimateParameters();
  log << "\t\tResultant join:\n";
  joined.printParameters(log,2); // print the joined pair mixture

  std::vector<Kent> joined_comp = joined.getComponents();
  components_m[index++] = joined_comp[0];
  Vector data_weights_m(N,1);
  Mixture_ML modified(K-1,components_m,weights_m,sample_size_m,responsibility_m,data,data_weights_m);
  log << "\t\tBefore adjustment ...\n";
  modified.computeNegativeLogLikelihood(data);
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
void Mixture_ML::generateHeatmapData(double res)
{
  string data_fbins2D = "./visualize/prob_bins2D.dat";
  string data_fbins3D = "./visualize/prob_bins3D.dat";
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
int Mixture_ML::getNearestComponent(int c)
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

double Mixture_ML::computeKLDivergence(Mixture_ML &other)
{
  return computeKLDivergence(other,data);
}

// Monte Carlo
double Mixture_ML::computeKLDivergence(Mixture_ML &other, std::vector<Vector> &sample)
{
  double kldiv = 0,fx,log_fx,log_gx;
  for (int i=0; i<sample.size(); i++) {
    log_fx = log_probability(sample[i]);
    log_gx = other.log_probability(sample[i]);
    kldiv += (log_fx - log_gx);
  }
  return kldiv/(log(2) * sample.size());
}

