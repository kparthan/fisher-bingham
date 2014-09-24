#include "vMF.h"
#include "Support.h"
#include "Normal.h"

/*!
 *  \brief This is a constructor module
 */
vMF::vMF()
{
  D = 3;
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
  D = mu.size();
  updateConstants();
}

/*!
 *  \brief This function computes the constants associated with the vMF.
 */
void vMF::updateConstants()
{
  kmu = Vector(D,0);
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
    long double log_bessel = log(cyl_bessel_i(D/2.0-1,kappa));
    /*if (log_bessel >= INFINITY) {
      my_float d2_1 = (D / 2.0) - 1;
      my_float bessel = cyl_bessel_i(d2_1,kappa);
      my_float large_log_bessel = log(bessel);
      log_bessel = large_log_bessel.convert_to<long double>();
      cout << "Normalization constant error ...\n";
      cout << scientific << "bessel: " << bessel;
      cout << scientific << "; log_bessel: " << log_bessel << endl;
      cout << "D: " << D << "; kappa: " << kappa << endl;
    }*/
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
Vector vMF::mean(void)
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
long double vMF::negativeLogLikelihood(Vector &x)
{
  return -log_density(x);
}

/*!
 *  \brief This function computes the negative log likelihood of given data.
 *  \param sample a reference to a vector<vector<long double> >
 *  \return the negative log likelihood (base e)
 */
long double vMF::negativeLogLikelihood(std::vector<Vector> &sample)
{
  long double value = 0;
  int N = sample.size();
  value -= N * log_cd;
  Vector sum(D,0);
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
void vMF::generateCanonical(std::vector<Vector> &canonical_sample, int sample_size)
{
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
      Vector random_normal = normal.generate(D-1);
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
std::vector<Vector> vMF::generate(int sample_size)
{
  cout << "\nGenerating from vMF with mean: ";
  if (D == 3) {
    Vector spherical(3,0);
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
    std::vector<Vector> canonical_sample;
    generateCanonical(canonical_sample,sample_size);
    if (fabs(mu[D-1] - 1) <= TOLERANCE) {  // check if mu is Z-axis
      return canonical_sample;
    } else {
      Matrix transformation = align_zaxis_with_vector(mu);
      return transform(canonical_sample,transformation);
    }
  } else if (sample_size == 0) {
    return std::vector<Vector>(); 
  }
}

void vMF::computeAllEstimators(
  std::vector<Vector> &data,
  std::vector<struct Estimates_vMF> &all_estimates
) {
  all_estimates.clear();

  Vector weights(data.size(),1.0);

  struct Estimates_vMF mlapprox_est;
  estimateMean(mlapprox_est,data,weights);

  // ML_APPROX
  estimateMLApproxKappa(mlapprox_est);
  all_estimates.push_back(mlapprox_est);

  // MLE
  struct Estimates_vMF ml_est = mlapprox_est;
  estimateMLKappa(ml_est);
  all_estimates.push_back(ml_est);

  // MAP
  struct Estimates_vMF map_est = mlapprox_est;
  estimateMAPKappa(map_est);
  all_estimates.push_back(map_est);

  // MML
  struct Estimates_vMF mml_est = mlapprox_est;
  estimateMMLKappa(mml_est);
  all_estimates.push_back(map_est);
}

void vMF::estimateMean(
  struct Estimates_vMF &estimates, 
  std::vector<Vector> &data, 
  Vector &weights
) {
  Vector resultant(3,0);  // resultant direction
  for (int i=0; i<data.size(); i++) {
    for (int j=0; j<3; j++) {
      resultant[j] += data[i][j];
    }
  }

  Vector mean(3,0);
  estimates.R = normalize(resultant,mean); // norm of resultant
  estimates.Rbar = estimates.R / data.size();
}

void VonMises::estimateMLApproxKappa(struct Estimates_vMF &estimates)
{
  long double rbar = estimates.Rbar;
  long double num = rbar * (3 - (rbar * rbar));
  long double denom = 1 - (rbar * rbar);

  estimates.kappa = num / denom;

  cout << "Kappa (ML approx): " << estimates.kappa << endl;
}

void VonMises::estimateMLKappa(struct Estimates_vMF &estimates)
{
  tryRootFindingMethods(estimates.kappa,
                        &VonMises::computeFunctionValue_ML,
                        &VonMises::computeFunctionAndDerivative_ML);
}

long double VonMises::computeFunctionValue_ML(long double x)
{
  long double A = computeRatioBessel(D,x);                   // A
  long double Ader = computeDerivativeOfRatioBessel(D,x,A);  // A'

  return A - estimates.Rbar;          // A - Rbar
}

void VonMises::computeFunctionAndDerivative_ML(
  long double x, 
  long double &Fx, 
  long double &Fx_der
) {
  long double A = computeRatioBessel(D,x);                   // A
  long double Ader = computeDerivativeOfRatioBessel(D,x,A);  // A'

  Fx = A - estimates.Rbar;  // A - Rbar 
  Fx_der = Ader;            // A'         
}

void VonMises::tryRootFindingMethods(
  long double &kappa_est,
  long double (VonMises::*function_value)(long double),
  void (VonMises::*function_and_derivative)(long double, long double &, long double &)
) {
  array<long double,2> interval;
  interval[0] = -1; interval[1] = -1;

  long double initial = kappa_est;

  // try newton raphson
  bool status = 0;
  status = NewtonRaphson(interval,kappa_est,function_and_derivative);

  if (!status) {  // newton raphson failed ...
    // check whether interval is valid
    if (interval[0] > 0 && interval[1] > 0) {
      status = IntervalBisection(interval,kappa_est,function_value); // unsuitable interval ...
    }
    if (!status) {
      // determine an interval to run the root bisection method
      status = determineInterval(interval,function_value,initial);
      if (!status) {
        cout << "No interval found ...\n";
        goto stop_searching;
      } else {
        status = IntervalBisection(interval,kappa_est,function_value); // try to find a root ...
      }
    } 
  }

  stop_searching:
  if (!status) {  // every attempt to find a root failed at this point ...
    kappa_est = initial;
  }
}

bool VonMises::NewtonRaphson(
  array<long double,2> &interval, 
  long double &kappa_est,
  void (VonMises::*function_and_derivative)(long double, long double &, long double &)
) {
  long double prev = estimates.kappa_ml_approx; // initial value
  long double current;
  int NUM_ITERATIONS = 20;
  long double Fx,Fx_der;
  long double new_range,current_range;

  for (int iter=1; iter<=NUM_ITERATIONS; iter++) {
    // get the function and derivative value
    (this->*function_and_derivative)(prev,Fx,Fx_der);

    //update(iter,prev,current,Fx,Fx_der,kappa_est,interval);
    if (fabs(Fx) <= ZERO) {
      assert(prev > 0);
      kappa_est = prev;
      cout << "Iteration " << iter << " Fx(" << prev << ") ~ 0" << endl;
      return 1;
    } else if (fabs(Fx_der) > TOLERANCE) {
      current = prev - (Fx/Fx_der);
      cout << "Iteration " << iter << ": [" << prev << ", " 
           << current <<  ", " << Fx << ", " << Fx_der << "]" << endl;
      if (fabs(current - prev) > TOLERANCE) {
        prev = fabs(current);
        // check if prev is HUGE ...
        //if (prev >= DEFAULT_MAX_KAPPA) return 0;
      } else {
        cout << "No significant change in Kappa ..." << endl;
        cout << "current: " << current << endl;
        if (current < 0) {
          kappa_est = prev;
        } else {
          kappa_est = current;
        }
        return 1;
      }
    } else if (fabs(Fx_der) <= ZERO) {    // Fx_der ~ 0
      cout << "Iteration " << iter << ": [" << prev << ", " 
           << current <<  ", " << Fx << ", " << Fx_der << "]" << endl;
      cout << "Derivative is zero ..." << endl;
      //kappa_est = prev;
      return 0;
    }
  } // iter loop ends ...

  cout << "Newton Raphson failed to converge ...\n";
  return 0;
}


/*!
 *  \brief This function computes the MML estimate using the Interval
 *  Bisection method.
 *  \param interval a reference to an array<long double,2>
 *  \param kappa_est a reference to a long double
 *  \param function_value a pointer to a member function that computes the
 *  function value to minimize
 *  \return the implementation status : convergence(success) or not
 */
bool VonMises::IntervalBisection(
  array<long double,2> &interval, 
  long double &kappa_est,
  long double (VonMises::*function_value)(long double)
) {
  long double a = interval[0];
  long double fa = (this->*function_value)(a);
  if (fabs(fa) <= ZERO) {
    kappa_est = a;
    return 1;
  }
  long double b = interval[1];
  long double fb = (this->*function_value)(b);
  if (fabs(fb) <= ZERO) {
    kappa_est = b;
    return 1;
  }
  // assert Fx has opposite signs
  if ((sign(fa) * sign(fb)) != -1) {
    cout << "a: " << a << " fa: " << fa;
    cout << "\tb: " << b << " fb: " << fb << endl;
    cout << "Unsuitable interval ...\n";
    return 0;
  }

  long double mid = (a+b)/ 2;
  int iter = 1;
  while(1) {
    long double fmid = (this->*function_value)(mid);
    cout << "Iteration " << iter++ << ": [" << a << ", " << b << ", " << mid << ", "
         << fa << ", " << fb << ", " << fmid << "]\n";
    if (fabs(fmid) <= ZERO) {
      kappa_est = mid;
      return 1;
    } else if (sign(fa) != sign(fmid)) {
      b = mid;
    } else {
      a = mid;
    }
    long double previous = mid;
    mid = (a + b) / 2;
    if (fabs(mid - previous) < TOLERANCE) {
      kappa_est = mid;
      return 1;
    }
    fa = (this->*function_value)(a);
    fb = (this->*function_value)(b);
  }
  return 0;   // should not happen ...
}

/*!
 *  \brief This function determines a suitable interval in which the root exists.
 *  \param interval a reference to an array<long double,2>
 *  \param function_value a pointer to a member function that computes the
 *  function value to minimize
 */
bool VonMises::determineInterval(
  array<long double,2> &interval,
  long double (VonMises::*function_value)(long double),
  long double kappa_ml_approx
) {
  // try to find a value lower than Kappa ML approx ...
  long double b = kappa_ml_approx;
  interval[1] = b;
  long double fb = (this->*function_value)(b);
  int signb = sign(fb);
  if (signb == 0) {
    interval[0] = b;
    cout << "Interval: (" << interval[0] << ", " << interval[1] << ")\n";
    return 1;
  }

  // compute appropriate increment/decrement based on the value of kappa_ml_approx
  long double CHANGE = 0.01 * kappa_ml_approx;

  cout << "Determining interval ...\n";
  long double decrement = CHANGE;
  long double a = b - decrement;
  long double fa;
  int signa;
  while (a > AOM) {
    fa = (this->*function_value)(a);
    signa = sign(fa);
    if (!(signa == 0 || signa == 1 || signa == -1)) {
      cout << "Error in calculating the function value ...\n";
      exit(1);
    }
    if (signa == 0 || signa * signb == -1) {
      interval[0] = a;
      cout << "Interval: (" << interval[0] << ", " << interval[1] << ")\n";
      return 1;
    } 
    a -= decrement;
  }

  return 0;
}

void VonMises::estimateMMLKappa(struct Estimates_vMF &estimates)
{
  long double prev = estimates.kappa;
  long double kappa_est,current,num,denom;
  long double Fx,Fx_der;
  bool status = 0;

  for (int i=0; i<20; i++) {
    computeFunctionAndDerivative_MML(prev,Fx,Fx_der);
    if (fabs(Fx) <= TOLERANCE) {
      cout << "Iteration " << i+1 << " Fx(" << prev << ") ~ 0" << endl;
      status = 1;
      break;
    }

    current = prev - (Fx/Fx_der);
    if (current < 0) {
      status = 0;
      break;
    }
    prev = current;
  }

  estimates.kappa = prev;
  if (!status) {
    array<long double,2> interval;
    interval[0] = -1; interval[1] = -1;
    // determine an interval to run the root bisection method
    status = determineInterval(interval,&VonMises::computeFunctionValue_MML);
    if (!status) {
      cout << "No interval found ...\n";
    } else {
      status = IntervalBisection(interval,estimates.kappa,&VonMises::computeFunctionValue_MML); // try to find a root ...
    }
  } 

  //validate_kappa(kappa_est,estimates.kappa_mml_complete);

  cout << "Kappa (MML Complete): " << estimates.kappa << endl;
}

void 
VonMises::computeFunctionAndDerivative_MML(
  long double x, 
  long double &Fx, 
  long double &Fx_der
) {
  long double A = computeRatioBessel(D,x);                   // A
  long double Ader = computeDerivativeOfRatioBessel(D,x,A);  // A'
  long double Ader_A = Ader / A;                             // A'/A
  long double A2der_Ader = A2der_Over_Ader(D,x,A,Ader_A);    // A''/A'

  Fx = computeFirstDerivativeOfMsglen(x,A,Ader_A,A2der_Ader);  // dI/dk
  Fx_der = computeSecondDerivativeOfMsglen(x,A,Ader,Ader_A,A2der_Ader);// d^2 I/dk^2
  if (boost::math::isnan(Fx) || boost::math::isnan(Fx_der)) {
    cout << "Error in computeFunctionAndDerivative_MML Fx: ";
    cout << "x: " << x << endl;
    cout << "A: " << A << "; Ader: " << Ader << "; Ader_A: " << Ader_A;
    cout << "; A2der_Ader: " << A2der_Ader << endl;
    cout << Fx << "; Fx_der: " << Fx_der << endl;
    //exit(1);
  }
}

/*!
 *  \brief This function computes the first derivative (wrt kappa) of the
 *  message length expression.
 *  \param k a reference to a long double
 *  \param Ad a reference to a long double
 *  \param Ader_A a reference to a long double
 *  \param A2der_Ader a reference to a long double
 *  \return the first derivative value
 */
long double 
VonMises::computeFirstDerivativeOfMsglen(
  long double &k, 
  long double &Ad, 
  long double &Ader_A, 
  long double &A2der_Ader
) {
  long double ans = -(D-1) / (2*k);
  ans += (D+1) * k / (1 + k*k);
  ans += ((D-1) /2.0) * Ader_A;
  ans += A2der_Ader;
  ans += estimates.Neff * Ad;
  ans -= estimates.R;
  return ans;
}

/*!
 *  \brief This function computes the first derivative (wrt kappa) of the
 *  message length expression.
 *  \param k a reference to a long double
 *  \param Ad a reference to a long double
 *  \param Ader a reference to a long double
 *  \param Ader_A a reference to a long double
 *  \param A2der_Ader a reference to a long double
 *  \return the first derivative value
 */
long double VonMises::computeSecondDerivativeOfMsglen(
  long double &k, long double &Ad, long double &Ader, long double &Ader_A, long double &A2der_Ader
) {
  long double ksq = k * k;
  long double ans = (D-1) / (2 * ksq);
  ans += ((D+1) * (1-ksq)) / ((1+ksq) * (1+ksq));
  ans += ((D-1)/2.0) * computeDerivativeOf_Ader_A(D,k,Ad,Ader);
  ans += computeDerivativeOf_A2der_Ader(D,k,Ader,Ader_A,A2der_Ader);
  ans += estimates.Neff * Ader;
  return ans;
}


