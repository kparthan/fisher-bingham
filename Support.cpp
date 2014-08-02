#include "Support.h"
#include "Normal.h"

int CONSTRAIN_KAPPA;

////////////////////// GENERAL PURPOSE FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

/*!
 *  \brief This function checks to see if valid arguments are given to the 
 *  command line output.
 *  \param argc an integer
 *  \param argv an std::vector of strings
 *  \return the parameters of the model
 */
struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;
  string constrain;

  bool noargs = 1;

  cout << "Checking command-line input ..." << endl;
  options_description desc("Allowed options");
  desc.add_options()
       ("help","produce help component")
       ("test","run some test cases")
  ;
  variables_map vm;
  store(command_line_parser(argc,argv).options(desc).run(),vm);
  notify(vm);

  if (vm.count("help")) {
    Usage(argv[0],desc);
  }

  if (vm.count("test")) {
    parameters.test = SET;
  } else {
    parameters.test = UNSET;
  }

  return parameters;
}

/*!
 *  \brief This module prints the acceptable input format to the program
 *  \param exe a reference to a const char
 *  \param desc a reference to a options_description object
 */
void Usage(const char *exe, options_description &desc)
{
  cout << "Usage: " << exe << " [options]" << endl;
  cout << desc << endl;
  exit(1);
}

/*!
 *  \brief This module checks whether the input file exists or not.
 *  \param file_name a reference to a string
 *  \return true or false depending on whether the file exists or not.
 */
bool checkFile(string &file_name)
{
  ifstream file(file_name.c_str());
  return file;
}

/*!
 *  \brief This module prints the elements of a std::vector<std::vector<long double><> > to a file
 *  \param v a reference to std::vector<std::vector<long double> >
 *  \param file_name a pointer to a const char
 */
void writeToFile(std::vector<std::vector<long double> > &v, const char *file_name)
{
  ofstream file(file_name);
  for (int i=0; i<v.size(); i++) {
    file << "(";
    for (int j=0; j<v[i].size()-1; j++) {
      file << v[i][j] << ", ";
    }
    file << v[i][v[i].size()-1] << ")" << endl;
  }
  file.close(); 
}

/*!
 *  \brief This module extracts the file name from the path
 *  \param file a reference to a string
 *  \return the extracted portion of the file name
 */
string extractName(string &file)
{
  unsigned pos1 = file.find_last_of("/");
  unsigned pos2 = file.find(".");
  int length = pos2 - pos1 - 1;
  string sub = file.substr(pos1+1,length);
  return sub;
}

/*!
 *  \brief This function prints the elements of an std::vector<long double>.
 *  \param os a reference to a ostream
 *  \param v a reference to a std::vector<long double>
 */
void print(ostream &os, std::vector<long double> &v, int precision)
{
  if (precision == 0) {
    if (v.size() == 1) {
      os << scientific << "(" << v[0] << ")";
    } else if (v.size() > 1) {
      os << scientific << "(" << v[0] << ", ";
      for (int i=1; i<v.size()-1; i++) {
        os << scientific << v[i] << ", ";
      }
      os << scientific << v[v.size()-1] << ")\t";
    } else {
      os << "No elements in v ...";
    }
  } else if (precision != 0) { // scientific notation
    if (v.size() == 1) {
      os << fixed << setprecision(3) << "(" << v[0] << ")";
    } else if (v.size() > 1) {
      os << fixed << setprecision(3) << "(" << v[0] << ", ";
      for (int i=1; i<v.size()-1; i++) {
        os << fixed << setprecision(3) << v[i] << ", ";
      }
      os << fixed << setprecision(3) << v[v.size()-1] << ")\t";
    } else {
      os << "No elements in v ...";
    }
  }
}

////////////////////// MATH FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

/*!
 *  \brief This module returns the sign of a number.
 *  \param number a long double
 *  \return the sign
 */
int sign(long double number)
{
  if (fabs(number) <= ZERO) {
    return 0;
  } else if (number > 0) {
    return 1;
  } else {
    return -1;
  }
}

/*!
 *  \brief This function computes the exponent a^x
 *  \param a a long double
 *  \param x a long double
 *  \return the exponent value
 */
long double exponent(long double a, long double x)
{
  assert(a > 0);
  long double tmp = x * log(a);
  return exp(tmp);
}

/*!
 *  \brief Normalizes a std::vector<long double>
 *  \param x a reference to a std::vector<long double>
 *  \param unit a reference to a std::vector<long double>
 *  \return the norm of the std::vector<long double>
 */
long double normalize(std::vector<long double> &x, std::vector<long double> &unit)
{
  long double normsq = 0;
  for (int i=0; i<x.size(); i++) {
    normsq += x[i] * x[i];
  }
  long double norm = sqrt(normsq);
  for (int i=0; i<x.size(); i++) {
    unit[i] = x[i] / norm;
  }
  return norm;
}

/*!
 *  \brief This function converts the cartesian coordinates into spherical.
 *  \param cartesian a reference to a std::vector<long double> 
 *  \param spherical a reference to a std::vector<long double> 
 */
void cartesian2spherical(std::vector<long double> &cartesian, std::vector<long double> &spherical)
{
  std::vector<long double> unit(3,0);
  long double r = normalize(cartesian,unit);

  long double x = unit[0];
  long double y = unit[1];
  long double z = unit[2];

  // theta \in [0,PI]: angle with Z-axis
  long double theta = acos(z);

  // phi \in[0,2 PI]: angle with positive X-axis
  long double ratio = x/sin(theta);
  if (ratio > 1) {
    ratio = 1;
  } else if (ratio < -1) {
    ratio = -1;
  }
  long double angle = acos(ratio);
  long double phi = 0;
  if (x == 0 && y == 0) {
    phi = 0;
  } else if (x == 0) {
    if (y > 0) {
      phi = angle;
    } else {
      phi = 2 * PI - angle;
    }
  } else if (y >= 0) {
    phi = angle;
  } else if (y < 0) {
    phi = 2 * PI - angle;
  }

  spherical[0] = r;
  spherical[1] = theta;
  spherical[2] = phi;
}

/*!
 *  \brief This function converts the spherical coordinates into cartesian.
 *  \param spherical a reference to a std::vector<long double> 
 *  \param cartesian a reference to a std::vector<long double> 
 */
void spherical2cartesian(std::vector<long double> &spherical, std::vector<long double> &cartesian)
{
  cartesian[0] = spherical[0] * sin(spherical[1]) * cos(spherical[2]);
  cartesian[1] = spherical[0] * sin(spherical[1]) * sin(spherical[2]);
  cartesian[2] = spherical[0] * cos(spherical[1]);
}

/*!
 *  \brief This funciton computes the dot product between two std::vector<long double>s.
 *  \param v1 a reference to a std::vector<long double>
 *  \param v2 a reference to a std::vector<long double>
 *  \return the dot product
 */
long double computeDotProduct(std::vector<long double> &v1, std::vector<long double> &v2) 
{
  assert(v1.size() == v2.size());
  long double dot_product = 0;
  for (int i=0; i<v1.size(); i++) {
    dot_product += v1[i] * v2[i];
  }
  return dot_product;
}

/*!
 *  \brief This function computes the surface area of nd-sphere
 *  Surface area = Gamma(d/2+1)/d \pi^(d/2)
 *  \param d an integer
 *  \return the surface area
 */
long double computeLogSurfaceAreaSphere(int d)
{
  long double log_num = log(d) + ((d/2.0) * log(PI));
  long double log_denom = boost::math::lgamma<long double>(d/2.0+1);
  return (log_num - log_denom);
}

/*!
 *  \brief This function computes the log of modified Bessel function value
 *  (used for numerical stability reasons)
 */
long double logModifiedBesselFirstKind(long double alpha, long double x)
{
  if (!(alpha >= 0 && fabs(x) >= 0)) {
    cout << "Error logModifiedBesselFirstKind: (alpha,x) = (" << alpha << "," << x << ")\n";
    //exit(1);
  }
  long double t;
  if (alpha == 0) {
    t = 1;
  } else if (fabs(x) <= TOLERANCE) {
    t = 0;
    return 0;
  } 
  long double x2_4 = x * x * 0.25;
  long double R = 1.0,             // 0-th term
         I = 1.0,             // m=0 sum
         m = 1.0;             // next, m=1
  // m! = m*(m-1)!, and
  // Gamma(m+alpha+1) = (m+alpha)*Gamma(m+alpha), 
  // because Gamma(x+1) = x*Gamma(x), for all x > 0.
  do { 
    //long double tmp = log(x2_4) - log(m) - log(alpha+m);
    long double tmp = x2_4 / (m*(alpha+m));  // the m-th term
    R *= tmp;  // the m-th term
    I += R;                     // total
    if (R >= INFINITY || I >= INFINITY) {
      return INFINITY;
    }
    m += 1.0;
    //cout << "m: " << m << "; tmp: " << tmp << "; R: " << R << "; I: " << I << endl;
  } while( R >= I * TOLERANCE);
  long double log_mod_bessel = log(I) + (alpha * log(x/2.0)) - boost::math::lgamma<long double>(alpha+1);
  /*if (log_mod_bessel >= INFINITY) {
    log_mod_bessel = approximate_bessel(alpha,x);
  }*/
  return log_mod_bessel;
}

/*!
 *  \brief This function computes the orthogonal transformation matrix
 *  to align a vector x = (0,...,0,1)^T with another vector y
 *  Source: http://math.stackexchange.com/questions/598750/finding-the-rotation-matrix-in-n-dimensions
 *  \param y a reference to a vector<long double>
 *  \return the transformation matrix
 */
matrix<long double> computeOrthogonalTransformation(std::vector<long double> &y)
{
  int D = y.size();

  std::vector<long double> x(D,0);
  x[D-1] = 1;
  std::vector<long double> u = x;
  long double uy = y[D-1];
  std::vector<long double> vn(D,0),v(D,0);
  for (int i=0; i<D; i++) {
    vn[i] = y[i] - uy * x[i];
  }
  normalize(vn,v);
  matrix<long double> I = identity_matrix<long double>(D,D);
  boost_vector U(D),V(D);
  for (int i=0; i<D; i++) {
    U[i] = u[i];
    V[i] = v[i];
  }

  matrix<long double> UUt = outer_prod(U,U);
  matrix<long double> VVt = outer_prod(V,V);
  matrix<long double> tmp1 = UUt + VVt;
  matrix<long double> tmp2 = I - tmp1;
  matrix<long double> UV(D,2);
  for (int i=0; i<D; i++) {
    UV(i,0) = U(i);
    UV(i,1) = V(i);
  }
  matrix<long double> UVt = trans(UV);
  long double cost = 0;
  for (int i=0; i<D; i++) {
    cost += x[i] * y[i];
  }
  long double sint = sqrt(1-cost*cost);
  matrix<long double> R(2,2);
  R(0,0) = cost;
  R(0,1) = -sint;
  R(1,0) = sint;
  R(1,1) = cost;
  matrix<long double> tmp3 = prod(UV,R);
  matrix<long double> tmp4 = prod(tmp3,UVt);
  matrix<long double> Q = tmp2 + tmp4;
  return Q;
}

/*
 *  \brief Transformation of x using T
 *  \param x a reference to a vector<vector<long double> >
 *  \param T a reference to a Matrix<long double>
 *  \return the transformed vector list
 */
std::vector<std::vector<long double> > transform(std::vector<std::vector<long double> > &x, matrix<long double> &T)
{
  int N = x.size();
  int D = x[0].size();
  matrix<long double> X(D,N);
  for (int i=0; i<N; i++) {
    for (int j=0; j<D; j++) {
      X(j,i) = x[i][j];
    }
  }
  matrix<long double> Y = prod(T,X);
  std::vector<std::vector<long double> > y;
  std::vector<long double> tmp(D,0);
  for (int i=0; i<N; i++) {
    for (int j=0; j<D; j++) {
      tmp[j] = Y(j,i);
    }
    y.push_back(tmp);
  }
  return y;
}

/*!
 *  Matrix inverse C++ Boost::ublas
 */
bool invertMatrix(const matrix<long double> &input, matrix<long double> &inverse)
{
  typedef permutation_matrix<std::size_t> pmatrix;

  // create a working copy of the input
  matrix<long double> A(input);

  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

  // perform LU-factorization
  int res = lu_factorize(A, pm);
  if (res != 0)
    return false;

  // create identity matrix of "inverse"
  inverse.assign(identity_matrix<long double> (A.size1()));

  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);

  return true;
}

/*!
 *  Eigen decomposition
 */
void eigenDecomposition(
  matrix<long double> m, 
  boost_vector &eigen_values,
  matrix<long double> &eigen_vectors
) {
  // check if m is symmetric
  int num_rows = m.size1();
  int num_cols = m.size2();
  if (num_rows != num_cols) {
    cout << "Error: rows: " << num_rows << " != cols: " << num_cols << endl;
    exit(1);
  }
  for (int i=0; i<num_rows; i++) {
    for (int j=0; j<num_cols; j++) {
      if (fabs(m(i,j)-m(j,i)) >= ZERO) {
        cout << "Error: Matrix is not symmetric ...\n";
        cout << "m: " << m << endl;
        cout << "m(" << i << "," << j << ") != m(" << j << "," << i << ")\n";
        exit(1);
      }
    }
  }

  // matrix is now validated ...
  int MAX_ITERATIONS = 100;
  for (int i=0; i < MAX_ITERATIONS; i++) {
    //find the largest off-diagonal 
    int max_row = 0, max_col = 1;
    int cur_row, cur_col;
    long double max_val = m(max_row,max_col);
    for (cur_row = 0; cur_row < num_rows-1; ++cur_row) {
      for (cur_col = cur_row + 1; cur_col < num_cols; ++cur_col) {
        if (fabs(m(cur_row,cur_col)) > max_val) {
          max_row = cur_row;
          max_col = cur_col;
          max_val = fabs(m(cur_row,cur_col));
        }
      }
    }

    if (max_val <= ZERO) {
      break; //finished
    }

    jacobiRotateMatrix(m,eigen_vectors, max_row, max_col);
  }

  for (int i = 0; i < num_cols; i++) {
    eigen_values(i) = m(i,i);
  }

  cout << "eigen_values: " << eigen_values << endl;
  cout << "eigen_vectors: " << eigen_vectors << endl;
}

void jacobiRotateMatrix(
  matrix<long double> &m,
  matrix<long double> &eigen_vectors, 
  int max_row, 
  int max_col
) {
  long double diff = m(max_col,max_col) - m(max_row,max_row);
  long double phi, t, c, s, tau, temp;
  int i;
  
  phi = diff / (2.0 * m(max_row,max_col));
  t = 1.0 / (std::fabs(phi) + std::sqrt((phi*phi) + 1.0));
  if(phi < 0){ t = -t; }

  c = 1.0 / std::sqrt(t*t + 1.0);
  s = t*c;
  tau = s/(1.0 + c);

  temp = m(max_row,max_col);
  m(max_row,max_col) = 0;
  m(max_row,max_row) = m(max_row,max_row) - (t*temp);
  m(max_col,max_col) = m(max_col,max_col) + (t*temp);
  for(i = 0; i < max_row; i++){ // Case i < max_row
    temp = m(i,max_row);
    m(i,max_row) = temp - (s*(m(i,max_col) + (tau*temp)));
    m(i,max_col) = m(i,max_col) + (s*(temp - (tau*m(i,max_col))));
  }
  for(i = max_row + 1; i < max_col; i++){ // Case max_row < i < max_col
    temp = m(max_row,i);
    m(max_row,i) = temp - (s*(m(i,max_col) + (tau*m(max_row,i))));
    m(i,max_col) = m(i,max_col) + (s*(temp - (tau*m(i,max_col))));
  }
  for(i = max_col + 1; i < 3; i++){ // Case i > max_col
    temp = m(max_row,i);
    m(max_row,i) = temp - (s*(m(max_col,i) + (tau*temp)));
    m(max_col,i) = m(max_col,i) + (s*(temp - (tau*m(max_col,i))));
  }

  for (i = 0; i < 3; i++) { // update the transformation matrix
    temp = eigen_vectors(i,max_row);
    eigen_vectors(i,max_row) = temp
      - (s*(eigen_vectors(i,max_col) + (tau*temp)));
    eigen_vectors(i,max_col) = eigen_vectors(i,max_col)
      + (s*(temp - (tau*eigen_vectors(i,max_col))));
  }
  return;
}


void rhs(const state_type &x, state_type &dxdt, const double t)
{
    dxdt[0] = 1 - (2 * t * x[0]); 
}

void track(const state_type &x, const double t)
{
    cout << t << "\t" << x[0] << endl;
}

void integrate_function(double limit)
{
  state_type x (1,0.0);
  //integrate(rhs,x,0.0,limit,0.1,track);
  integrate(rhs,x,0.0,limit,0.1);
  cout << "ans: " << x[0] << endl;
}

////////////////////// TESTING FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

void Test(void)
{
  cout << "Testing matrices ...\n";

  //cout << "Matrix declaration ...\n";
  matrix<long double> m1 (3, 3);
  for (int i = 0; i < m1.size1(); ++ i) {
    for (int j = 0; j < m1.size2(); ++ j) {
      m1 (i, j) = 3 * i + j;
    }
  }
  cout << "m1: " << m1 << endl;
  cout << "2 * m1: " << 2*m1 << endl;
  cout << "m1/2: " << m1/2 << endl;

  //cout << "Matrix transpose ...\n";
  matrix<long double> m2;
  m2 = trans(m1);
  cout << "m1' = m2: " << m2 << endl;

  //cout << "Matrix inverse ...\n";
  matrix<long double> inverse(3,3);
  m1(0,0) = 1;
  invertMatrix(m1,inverse);
  cout << "m1: " << m1 << endl;
  cout << "inv(m1): " << inverse << endl;

  //cout << "Identity matrix ...\n";
  identity_matrix<long double> id(3,3);
  cout << "id: " << id << endl;
  matrix<long double> add = id + m1;
  cout << "id + m1: " << add << endl;

  // matrix row
  matrix_row<matrix<long double> > mr(m1,0);
  cout << "mr: " << mr << endl;

  // boost vector
  boost_vector v(3);
  for (int i=0; i<3; i++) {
    v[i] = i + 3;
  }
  cout << "v: " << v << endl;
  cout << "2 * v: " << 2 * v << endl;
  cout << "v/2: " << v/2 << endl;

  // adding matrix row and vector
  boost_vector v1 = mr + v;
  cout << "v1: " << v1 << endl;

  // multiplication
  matrix<long double> m3 = prod(m1,m2);
  cout << "m1 * m2 = m3: " << m3 << endl;

  // multiplying matrices and vectors
  boost_vector mv = prod(v1,m1);
  cout << "v1 * m1 = mv: " << mv << endl;
  mv = prod(m1,v1);
  cout << "m1 * v1 = mv: " << mv << endl;

  long double v1_T_v1 = inner_prod(v1,v1);
  cout << "v1' * v1 = : " << v1_T_v1 << endl;

  matrix<long double> m4 = outer_prod(v1,v1);
  cout << "v1 * v1' = m4: " << m4 << endl;

  // eigen values & vectors
  matrix<long double> symm = (m1 + trans(m1))/2;
  cout << "symmetric matrix: " << symm << endl;
  boost_vector eigen_values(3,0);
  matrix<long double> eigen_vectors = identity_matrix<long double>(3,3);
  eigenDecomposition(symm,eigen_values,eigen_vectors);

  symm = identity_matrix<long double>(3,3);
  cout << "symmetric matrix: " << symm << endl;
  eigen_vectors = identity_matrix<long double>(3,3);
  eigenDecomposition(symm,eigen_values,eigen_vectors);

  // integration
  integrate_function(10);

  // Normal cdf
  Normal normal(0,1);
  long double cdf = normal.cumulativeDensity(2);
  cout << "cdf: " << cdf << endl;

  long double x = sqrt(PI/2);
  cdf = normal.cumulativeDensity(x);
  cout << "cdf: " << cdf << endl;
  cout << "2*cdf-1: " << 2*cdf-1 << endl;
}

