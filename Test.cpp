#include "Test.h"
#include "Support.h"
#include "Normal.h"
#include "vMF.h"
#include "FB4.h"
#include "FB6.h"
#include "Kent.h"

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

void Test::productMatrixVector(void)
{
  matrix<long double> m1(4,3);
  for (int i=0; i<4; i++) {
    for (int j=0; j<3; j++) {
      m1(i,j) = i + 1 + j/2.0;
    }
  }
  std::vector<long double> v(3,0);
  v[0] = 1; v[1] = -2; v[2] = 3;
  cout << "m1: " << m1 << endl;
  cout << "v: "; print(cout,v,3); cout << endl;
  std::vector<long double> ans = prod(m1,v);
  cout << "m1*v: "; print(cout,ans,3); cout << endl;

  matrix<long double> m2(3,4);
  for (int i=0; i<3; i++) {
    for (int j=0; j<4; j++) {
      m2(i,j) = i + 1 + j/2.0;
    }
  }
  cout << "m2: " << m2 << endl;
  cout << "v: "; print(cout,v,3); cout << endl;
  ans = prod(v,m2);
  cout << "m2*v: "; print(cout,ans,3); cout << endl;
}

void Test::dispersionMatrix(void)
{
  string file_name = "./visualize/kent.dat";
  std::vector<std::vector<long double> > sample = load_matrix(file_name);
  std::vector<long double> unit_mean = computeVectorSum(sample);
  cout << "unit_mean: "; print(cout,unit_mean,3); cout << endl;
  matrix<long double> m = computeDispersionMatrix(sample);
  cout << "dispersion: " << m << endl;

  matrix<long double> eigen_vectors = identity_matrix<long double>(3,3);
  boost_vector eigen_values(3);
  eigenDecomposition(m,eigen_values,eigen_vectors);

  m(0,0) = 0.341; m(0,1) = -0.221; m(0,2) = 0.408;
  m(1,0) = -0.221; m(1,1) = 0.153; m(1,2) = -0.272;
  m(2,0) = 0.408; m(2,1) = -0.272; m(2,2) = 0.506;
  eigen_vectors = identity_matrix<long double>(3,3);
  eigenDecomposition(m,eigen_values,eigen_vectors);
}

void Test::numericalIntegration(void)
{
  // integration
  long double val = computeIntegral(10);
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

void Test::orthogonalTransformations2(void)
{
  std::vector<long double> m0,m1,m2;
  generateRandomOrthogonalVectors(m0,m1,m2);
  cout << "mean: "; print(cout,m0,3); cout << endl;
  cout << "major: "; print(cout,m1,3); cout << endl;
  cout << "minor: "; print(cout,m2,3); cout << endl;

  matrix<long double> r = computeOrthogonalTransformation(m0,m1);
  cout << "r: " << r << endl;
  std::vector<long double> ztransform = prod(r,ZAXIS);
  cout << "ztransform: "; print(cout,ztransform,3); cout << endl; // = mu0
  std::vector<long double> xtransform = prod(r,XAXIS);
  cout << "xtransform: "; print(cout,xtransform,3); cout << endl; // = mu1
  std::vector<long double> ytransform = prod(r,YAXIS);
  cout << "ytransform: "; print(cout,ytransform,3); cout << endl; // = mu2
}

void Test::randomSampleGeneration(void)
{
  std::vector<std::vector<long double> > random_sample;
  std::vector<long double> m0,m1,m2;
  generateRandomOrthogonalVectors(m0,m1,m2);

  // FB4 generation (gamma < 0)
  FB4 fb4_1(m0,m1,m2,100,-10);
  random_sample = fb4_1.generate(1000);
  writeToFile("./visualize/fb4_1.dat",random_sample,3);

  // FB4 generation (gamma < 0)
  FB4 fb4_2(m0,m1,m2,100,10);
  random_sample = fb4_2.generate(1000);
  writeToFile("./visualize/fb4_2.dat",random_sample,3);

  // vMF generation
  vMF vmf(m0,100);
  random_sample = vmf.generate(1000);
  writeToFile("./visualize/vmf.dat",random_sample,3);

  // vMF (2D)
  std::vector<long double> mean(2,0); mean[0] = 1;
  vMF vmf2(mean,10);
  vmf2.generateCanonical(random_sample,100);
  writeToFile("./visualize/vmf2.dat",random_sample,3);

  // FB6 generation 
  FB6 fb6(m0,m1,m2,100,15,-10);
  random_sample = fb6.generate(1000);
  writeToFile("./visualize/fb6.dat",random_sample,3);

  // FB6 generation 
  FB6 kent(m0,m1,m2,1000,475,0);
  random_sample = kent.generate(1000);
  writeToFile("./visualize/kent.dat",random_sample,3);

  // FB6 generation 
  /*FB6 kent1(m0,m1,m2,1,60,-10);
  random_sample = kent1.generate(1000);
  writeToFile("./visualize/kent1.dat",random_sample,3);*/
}

void Test::normalization_constant(void)
{
  cout << "ZERO: " << ZERO << endl;
  std::vector<long double> m0 = ZAXIS;
  std::vector<long double> m1 = XAXIS;
  std::vector<long double> m2 = YAXIS;
  matrix<long double> a1 = outer_prod(m1,m1);
  cout << "a1: " << a1 << endl;
  matrix<long double> a2 = outer_prod(m2,m2);
  cout << "a2: " << a2 << endl;
  matrix<long double> A = a1 - a2;
  cout << "A: " << A << endl;
  //Kent kent(100,30);
  Kent kent(100,47.5);
  kent.computeLogNormalizationConstant();
}

