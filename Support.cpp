#include "Support.h"

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
 *  \brief This module prints the elements of a std::vector<std::vector<> > to a file
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
 *  \brief This function prints the elements of an std::vector.
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
 *  \brief Normalizes a std::vector
 *  \param x a reference to a std::vector<long double>
 *  \param unit a reference to a std::vector<long double>
 *  \return the norm of the std::vector
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
 *  \brief This funciton computes the dot product between two std::vectors.
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
 *  Matrix inverse C++ Boost::ublas
 */
bool invertMatrix(const matrix<long double>& input, matrix<long double>& inverse)
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
  boost::numeric::ublas::vector<long double> v(3);
  for (int i=0; i<3; i++) {
    v[i] = i + 3;
  }
  cout << "v: " << v << endl;
  cout << "2 * v: " << 2 * v << endl;
  cout << "v/2: " << v/2 << endl;

  // adding matrix row and vector
  boost::numeric::ublas::vector<long double> v1 = mr + v;
  cout << "v1: " << v1 << endl;

  // multiplication
  matrix<long double> m3 = prod(m1,m2);
  cout << "m1 * m2 = m3: " << m3 << endl;

  // multiplying matrices and vectors
  boost::numeric::ublas::vector<long double> mv = prod(v1,m1);
  cout << "v1 * m1 = mv: " << mv << endl;
  mv = prod(m1,v1);
  cout << "m1 * v1 = mv: " << mv << endl;

  long double v1_T_v1 = inner_prod(v1,v1);
  cout << "v1' * v1 = : " << v1_T_v1 << endl;

  matrix<long double> m4 = outer_prod(v1,v1);
  cout << "v1 * v1' = m4: " << m4 << endl;

  // eigen values & vectors
}

