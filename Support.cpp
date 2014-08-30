#include "Support.h"
#include "Test.h"

Vector XAXIS,YAXIS,ZAXIS;

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
 *  \brief This module prints the elements of a std::vector<Vector<> > to a file
 *  \param v a reference to std::vector<Vector >
 *  \param file_name a pointer to a const char
 */
void writeToFile(const char *file_name, std::vector<Vector > &v, int precision)
{
  ofstream file(file_name);
  for (int i=0; i<v.size(); i++) {
    for (int j=0; j<v[i].size(); j++) {
      file << fixed << setw(10) << setprecision(precision) << v[i][j];
    }
    file << endl;
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
 *  \brief This function prints the elements of an Vector.
 *  \param os a reference to a ostream
 *  \param v a reference to a Vector
 */
void print(ostream &os, Vector &v, int precision)
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
 *  \brief Normalizes a Vector
 *  \param x a reference to a Vector
 *  \param unit a reference to a Vector
 *  \return the norm of the Vector
 */
long double normalize(Vector &x, Vector &unit)
{
  long double l2norm = norm(x);
  for (int i=0; i<x.size(); i++) {
    unit[i] = x[i] / l2norm;
  }
  return l2norm;
}

long double norm(Vector &v)
{
  long double normsq = 0;
  for (int i=0; i<v.size(); i++) {
    normsq += v[i] * v[i];
  }
  return sqrt(normsq);
}

/*!
 *  \brief This function converts the cartesian coordinates into spherical.
 *  (theta with +Z and phi with +X)
 *  \param cartesian a reference to a Vector 
 *  \param spherical a reference to a Vector 
 */
void cartesian2spherical(Vector &cartesian, Vector &spherical)
{
  Vector unit(3,0);
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
 *  \brief This function converts the cartesian coordinates into spherical.
 *  (theta with +X and phi with +Y)
 *  \param cartesian a reference to a Vector 
 *  \param spherical a reference to a Vector 
 */
void cartesian2sphericalPoleXAxis(Vector &cartesian, Vector &spherical)
{
  Vector unit(3,0);
  long double r = normalize(cartesian,unit);

  long double x = unit[0];
  long double y = unit[1];
  long double z = unit[2];

  // theta \in [0,PI]: angle with X-axis
  long double theta = acos(x);

  // phi \in[0,2 PI]: angle with positive Y-axis
  long double ratio = y/sin(theta);
  if (ratio > 1) {
    ratio = 1;
  } else if (ratio < -1) {
    ratio = -1;
  }
  long double angle = acos(ratio);
  long double phi = 0;
  if (y == 0 && z == 0) {
    phi = 0;
  } else if (y == 0) {
    if (z > 0) {
      phi = angle;
    } else {
      phi = 2 * PI - angle;
    }
  } else if (z >= 0) {
    phi = angle;
  } else if (z < 0) {
    phi = 2 * PI - angle;
  }

  spherical[0] = r;
  spherical[1] = theta;
  spherical[2] = phi;
}

/*!
 *  \brief This function converts the spherical coordinates into cartesian.
 *  \param spherical a reference to a Vector 
 *  \param cartesian a reference to a Vector 
 */
void spherical2cartesian(Vector &spherical, Vector &cartesian)
{
  cartesian[0] = spherical[0] * sin(spherical[1]) * cos(spherical[2]);
  cartesian[1] = spherical[0] * sin(spherical[1]) * sin(spherical[2]);
  cartesian[2] = spherical[0] * cos(spherical[1]);
}

/*!
 *  \brief This funciton computes the dot product between two Vectors.
 *  \param v1 a reference to a Vector
 *  \param v2 a reference to a Vector
 *  \return the dot product
 */
long double computeDotProduct(Vector &v1, Vector &v2) 
{
  assert(v1.size() == v2.size());
  long double dot_product = 0;
  for (int i=0; i<v1.size(); i++) {
    dot_product += v1[i] * v2[i];
  }
  return dot_product;
}

Vector crossProduct(Vector &v1, Vector &v2) 
{
  Vector ans(3,0);
  ans[0] = v1[1] * v2[2] - v1[2] * v2[1];
  ans[1] = v1[2] * v2[0] - v1[0] * v2[2];
  ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return ans;
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
    long double tmp = x2_4 / (m*(alpha+m));  // the m-th term
    R *= tmp;  // the m-th term
    I += R;                     // total
    if (R >= INFINITY || I >= INFINITY) {
      return INFINITY;
    }
    m += 1.0;
    //cout << "m: " << m << "; tmp: " << tmp << "; R: " << R << "; I: " << I << endl;
  } while( R >= I * ZERO);
  long double log_mod_bessel = log(I) + (alpha * log(x/2.0)) - lgamma<long double>(alpha+1);
  return log_mod_bessel;
}

/*!
 *  Find the roots of a quadratic: ax^2 + bx + c = 0
 */
void solveQuadratic(
  Vector &roots, 
  long double a,
  long double b, 
  long double c
) {
  long double D = sqrt(b*b-4*a*c);
  roots[0] = (-b + D) / (2*a);
  roots[1] = (-b - D) / (2*a);
}

////////////////////// GEOMETRY FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

std::vector<Vector > load_matrix(string &file_name)
{
  std::vector<Vector > sample;
  ifstream file(file_name.c_str());
  string line;
  Vector numbers(3,0),unit_vector(3,0);
  int i;
  while (getline(file,line)) {
    boost::char_separator<char> sep(" \t");
    boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
    i = 0;
    BOOST_FOREACH(const string &t, tokens) {
      istringstream iss(t);
      long double x;
      iss >> x;
      numbers[i++] = x;
    }
    normalize(numbers,unit_vector);
    sample.push_back(unit_vector);
  }
  file.close();
  return sample;
}

/*!
 *  v1 and v2 are considered to be a column std::vectors
 *  output: v1 * v2' (the outer product matrix)
 */
Matrix outer_prod(Vector &v1, Vector &v2)
{
  assert(v1.size() == v2.size());
  int m = v1.size();
  Matrix ans(m,m);
  for (int i=0; i<m; i++) {
    for (int j=0; j<m; j++) {
      ans(i,j) = v1[i] * v2[j];
    }
  }
  return ans;
}

/*!
 *  v is considered to be a column std::vector
 *  output: m * v (a row std::vector)
 */
Vector prod(Matrix &m, Vector &v)
{
  assert(m.size2() == v.size());
  Vector ans(m.size1(),0);
  for (int i=0; i<m.size1(); i++) {
    for (int j=0; j<m.size2(); j++) {
      ans[i] += m(i,j) * v[j];
    }
  }
  return ans;
}

/*!
 *  v is considered to be a column std::vector
 *  output: v' * m (a row std::vector)
 */
Vector prod(Vector &v, Matrix &m)
{
  assert(m.size1() == v.size());
  Vector ans(m.size2(),0);
  for (int i=0; i<m.size2(); i++) {
    for (int j=0; j<m.size1(); j++) {
      ans[i] += v[j] * m(j,i);
    }
  }
  return ans;
}

/*!
 *  v is considered to be a column std::vector
 *  output: v' M v
 */
long double prod_vMv(Vector &v, Matrix &M)
{
  Vector vM = prod(v,M);
  return computeDotProduct(vM,v);
}

/*!
 *  x,y are considered to be a column std::vectors
 *  output: x' M y
 */
long double prod_xMy(Vector &x, Matrix &M, Vector &y)
{
  Vector xM = prod(x,M);
  return computeDotProduct(xM,y);
}

/*!
 *  determinant of 3 X 3 matrix
 */
long double determinant(Matrix &m)
{
  long double det = 0,subdet;
  subdet = m(1,1) * m(2,2) - m(1,2) * m(2,1);
  det += m(0,0) * subdet;

  subdet = m(1,0) * m(2,2) - m(1,2) * m(2,0);
  det -= m(0,1) * subdet;

  subdet = m(1,0) * m(2,1) - m(1,1) * m(2,0);
  det += m(0,2) * subdet;

  return det;
}

/*!
 *  Computes \sum x (x is a vector)
 */
Vector computeVectorSum(std::vector<Vector > &sample) 
{
  int d = sample[0].size();
  Vector sum(d,0);
  for (int i=0; i<sample.size(); i++) {
    for (int j=0; j<d; j++) {
      sum[j] += sample[i][j];
    }
  }
  return sum;
}

/*!
 *  Computes \sum x / N (x is a vector)
 */
Vector computeNormalizedVectorSum(std::vector<Vector > &sample) 
{
  Vector sum = computeVectorSum(sample);
  for (int j=0; j<sum.size(); j++) {
    sum[j] /= sample.size();
  }
  return sum;
}

/*!
 *  Computes \sum x * x' (x is a vector)
 */
Matrix computeDispersionMatrix(std::vector<Vector > &sample)
{
  int d = sample[0].size();
  Matrix dispersion = ZeroMatrix(d,d);
  for (int i=0; i<sample.size(); i++) {
    dispersion += outer_prod(sample[i],sample[i]);
  }
  return dispersion;
}

/*!
 *  Computes \sum x * x' / N (x is a vector)
 */
Matrix computeNormalizedDispersionMatrix(std::vector<Vector > &sample)
{
  Matrix dispersion = computeDispersionMatrix(sample);
  return dispersion/sample.size();
}

/*!
 *  Matrix to rotate FROM the standard frame of reference.
 */
Matrix computeOrthogonalTransformation(
  Vector &mean,
  Vector &major_axis
) {
  Matrix r1 = align_zaxis_with_vector(mean);
  Matrix inverse(3,3);
  invertMatrix(r1,inverse);
  // rotate major axis onto XY-plane
  Vector major1 = prod(inverse,major_axis);
  Matrix r2 = align_xaxis_with_major_axis(major1);
  Matrix r = prod(r1,r2);
  return r;
}

Matrix align_xaxis_with_major_axis(Vector &major_axis)
{
  Vector spherical(3,0);
  cartesian2spherical(major_axis,spherical);
  long double theta = spherical[1]; // theta = PI/2
  long double phi = spherical[2];

  Matrix r = IdentityMatrix(3,3);
  r(0,0) = cos(phi);
  r(0,1) = -sin(phi);
  r(1,0) = -r(0,1); // sin(phi)
  r(1,1) = r(0,0); // cos(phi)
  return r;
}

Matrix align_zaxis_with_vector(Vector &y)
{
  Vector spherical(3,0);
  cartesian2spherical(y,spherical);
  long double theta = spherical[1];
  long double phi = spherical[2];

  Matrix r1 = IdentityMatrix(3,3);
  r1(0,0) = cos(theta);
  r1(0,2) = sin(theta);
  r1(2,0) = -r1(0,2); //-sin(theta)
  r1(2,2) = r1(0,0); //cos(theta)

  Matrix r2 = IdentityMatrix(3,3);
  r2(0,0) = cos(phi);
  r2(0,1) = -sin(phi);
  r2(1,0) = -r2(0,1); //sin(phi)
  r2(1,1) = r2(0,0);  // cos(phi)

  Matrix r = prod(r2,r1);
  return r;
}

void generateRandomOrthogonalVectors(
  Vector &mean,
  Vector &major_axis,
  Vector &minor_axis
) {
  long double phi = rand()*2*PI/(long double)RAND_MAX;
  Vector spherical(3,1),major1(3,0);
  spherical[1] = PI/2;
  spherical[2] = phi;
  spherical2cartesian(spherical,major1); // major axis
  Vector mu1 = ZAXIS;
  //Vector minor1 = crossProduct(mu1,major1);

  long double theta = rand()*PI/(long double)RAND_MAX;
  phi = rand()*2*PI/(long double)RAND_MAX;
  spherical[1] = theta;
  spherical[2] = phi;
  mean = Vector(3,0);
  spherical2cartesian(spherical,mean); 

  Matrix r = align_zaxis_with_vector(mean);
  major_axis = prod(r,major1);
  //minor_axis = prod(r,minor1);
  minor_axis = crossProduct(mean,major_axis);
}

/*
 *  \brief Transformation of x using T
 *  \param x a reference to a vector<vector<long double> >
 *  \param T a reference to a Matrix<long double>
 *  \return the transformed vector list
 */
std::vector<Vector > transform(
  std::vector<Vector > &x, 
  Matrix &T
) {
  std::vector<Vector > y(x.size());
  for (int i=0; i<x.size(); i++) {
    y[i] = prod(T,x[i]);
  }
  return y;
}

/*!
 *  Matrix inverse C++ Boost::ublas
 */
bool invertMatrix(const Matrix &input, Matrix &inverse)
{
  typedef permutation_matrix<std::size_t> pmatrix;

  // create a working copy of the input
  Matrix A(input);

  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

  // perform LU-factorization
  int res = lu_factorize(A, pm);
  if (res != 0)
    return false;

  // create identity matrix of "inverse"
  inverse.assign(IdentityMatrix (A.size1()));

  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);

  return true;
}

/*!
 *  Eigen decomposition
 *  Inputs:
 *                m -- a symmetric matrix
 *    eigen_vectors -- an identity matrix
 *  Outputs:
 *    eigen_vectors -- each column is a unit eigen vector
 */
void eigenDecomposition(
  Matrix m, 
  Vector &eigen_values,
  Matrix &eigen_vectors
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

    jacobiRotateMatrix(m,eigen_vectors,max_row,max_col);
  }

  for (int i = 0; i < num_cols; i++) {
    eigen_values[i] = m(i,i);
  }

  cout << "eigen_values: "; print(cout,eigen_values,0); cout << endl;
  cout << "eigen_vectors: " << eigen_vectors << endl;
}

void jacobiRotateMatrix(
  Matrix &m,
  Matrix &eigen_vectors, 
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
  for(i = max_col + 1; i < m.size2(); i++){ // Case i > max_col
    temp = m(max_row,i);
    m(max_row,i) = temp - (s*(m(max_col,i) + (tau*temp)));
    m(max_col,i) = m(max_col,i) + (s*(temp - (tau*m(max_col,i))));
  }

  for (i = 0; i < eigen_vectors.size1(); i++) { // update the transformation matrix
    temp = eigen_vectors(i,max_row);
    eigen_vectors(i,max_row) = temp
      - (s*(eigen_vectors(i,max_col) + (tau*temp)));
    eigen_vectors(i,max_col) = eigen_vectors(i,max_col)
      + (s*(temp - (tau*eigen_vectors(i,max_col))));
  }
  return;
}

/*
 *  Dawson's integral required to compute FB4 normalization constant
 */
long double computeDawsonsIntegral(double limit)
{
  //state_type x (1,0.0);
  std::vector<double> x (1,0.0);
  //integrate(rhs,x,0.0,limit,0.1,track);
  integrate(rhs,x,0.0,limit,0.1);
  //cout << "ans: " << x[0] << endl;
  return x[0];
}

void rhs(const std::vector<double> &x, std::vector<double> &dxdt, const double t)
{
    dxdt[0] = 1 - (2 * t * x[0]); 
}

void track(const std::vector<double> &x, const double t)
{
    cout << t << "\t" << x[0] << endl;
}

////////////////////// TESTING FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

void TestFunctions(void)
{
  Test test;

  //test.matrixFunctions();

  //test.productMatrixVector();

  //test.dispersionMatrix();

  //test.numericalIntegration();

  //test.normalDistributionFunctions();

  //test.orthogonalTransformations();

  //test.orthogonalTransformations2();

  //test.randomSampleGeneration();

  //test.normalization_constant();

  //test.optimization();

  //test.moment_estimation();

  test.ml_estimation();

  //test.expectation();

  //test.kl_divergence();

  //test.fisher();
}

