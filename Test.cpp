#include "Test.h"
#include "Support.h"
#include "Normal.h"
#include "vMF.h"
#include "FB4.h"

extern std::vector<long double> XAXIS,YAXIS,ZAXIS;

void Test::matrixFunctions()
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
}

void Test::numericalIntegration(void)
{
  // integration
  integrate_function(10);
}

void Test::normalDistributionFunctions(void)
{
  // Normal cdf
  Normal normal(0,1);
  long double cdf = normal.cumulativeDensity(2);
  cout << "cdf: " << cdf << endl;

  long double x = sqrt(PI/2);
  cdf = normal.cumulativeDensity(x);
  cout << "cdf: " << cdf << endl;
  cout << "2*cdf-1: " << 2*cdf-1 << endl;
}

void Test::orthogonalTransformations(void)
{
  std::vector<long double> spherical(3,0),cartesian(3,0);
  spherical[0] = 1;
  spherical[1] = PI/3;  // 60 degrees
  spherical[2] = 250 * PI/180;  // 250 degrees
  spherical2cartesian(spherical,cartesian);
  cout << "cartesian: "; print(cout,cartesian,3); cout << endl;

  boost_vector vec(3); convert2boostvector(cartesian,vec);
  matrix<long double> r1 = align_zaxis_with_vector(cartesian);
  cout << "r1: " << r1 << endl;
  matrix<long double> inverse(3,3);
  invertMatrix(r1,inverse);
  cout << "inverse: " << inverse << endl;
  matrix<long double> check = prod(r1,inverse);
  cout << "check: " << check << endl;

  boost_vector zaxis(3); convert2boostvector(ZAXIS,zaxis);
  boost_vector ans1 = prod(r1,zaxis);
  cout << "ans1: " << ans1 << endl;
  boost_vector ans2 = prod(inverse,vec);
  cout << "ans2: " << ans2 << endl;
}

void Test::randomSampleGeneration(void)
{
  // vMF generation
  std::vector<long double> mu = ZAXIS;
  vMF vmf(ZAXIS,100);
  std::vector<std::vector<long double> > random_sample = vmf.generate(1000);
  writeToFile("./visualize/vmf.dat",random_sample,3);

  // FB4 generation
  FB4 fb4(100,-10);
  random_sample = fb4.generate(1000);
  writeToFile("./visualize/fb4.dat",random_sample,3);
}

