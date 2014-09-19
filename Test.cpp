#include "Test.h"
#include "Support.h"
#include "Normal.h"
#include "vMF.h"
#include "FB4.h"
#include "FB6.h"
#include "Kent.h"

extern Vector XAXIS,YAXIS,ZAXIS;

void Test::matrixFunctions()
{
  cout << "Testing matrices ...\n";

  //cout << "Matrix declaration ...\n";
  Matrix m1 (3, 3);
  for (int i = 0; i < m1.size1(); ++ i) {
    for (int j = 0; j < m1.size2(); ++ j) {
      m1 (i, j) = 3 * i + j;
    }
  }
  cout << "m1: " << m1 << endl;
  cout << "2 * m1: " << 2*m1 << endl;
  cout << "m1/2: " << m1/2 << endl;

  //cout << "Matrix transpose ...\n";
  Matrix m2;
  m2 = trans(m1);
  cout << "m1' = m2: " << m2 << endl;

  //cout << "Matrix inverse ...\n";
  Matrix inverse(3,3);
  m1(0,0) = 1;
  invertMatrix(m1,inverse);
  cout << "m1: " << m1 << endl;
  cout << "inv(m1): " << inverse << endl;

  //cout << "Identity matrix ...\n";
  IdentityMatrix id(3,3);
  cout << "id: " << id << endl;
  Matrix add = id + m1;
  cout << "id + m1: " << add << endl;

  // matrix row
  matrix_row<Matrix > mr(m1,0);
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
  Matrix m3 = prod(m1,m2);
  cout << "m1 * m2 = m3: " << m3 << endl;

  // multiplying matrices and vectors
  boost::numeric::ublas::vector<long double> mv = prod(v1,m1);
  cout << "v1 * m1 = mv: " << mv << endl;
  mv = prod(m1,v1);
  cout << "m1 * v1 = mv: " << mv << endl;

  long double v1_T_v1 = inner_prod(v1,v1);
  cout << "v1' * v1 = : " << v1_T_v1 << endl;

  Matrix m4 = outer_prod(v1,v1);
  cout << "v1 * v1' = m4: " << m4 << endl;

  // eigen values & vectors
  Matrix symm = (m1 + trans(m1))/2;
  cout << "symmetric matrix: " << symm << endl;
  Vector eigen_values(3,0);
  Matrix eigen_vectors = IdentityMatrix(3,3);
  eigenDecomposition(symm,eigen_values,eigen_vectors);

  symm = IdentityMatrix(3,3);
  cout << "symmetric matrix: " << symm << endl;
  eigen_vectors = IdentityMatrix(3,3);
  eigenDecomposition(symm,eigen_values,eigen_vectors);

  symm(0,0) = 2.294628e+01; symm(0,1) = -2.162988e+01; symm(0,2) = -1.516247e+01;
  symm(1,0) = -2.162988e+01; symm(1,1) = -2.794379e+01; symm(1,2) = 2.987167e+01;
  symm(2,0) = -1.516247e+01; symm(2,1) = 2.987167e+01; symm(2,2) = 4.997508e+00;
  cout << "symmetric matrix: " << symm << endl;
  eigen_vectors = IdentityMatrix(3,3);
  eigenDecomposition(symm,eigen_values,eigen_vectors);

  symm(0,0) = 16.8974; symm(0,1) = -20.5575; symm(0,2) = 11.4795;
  symm(1,0) = -20.5575; symm(1,1) = -4.5362; symm(1,2) = -38.3720;
  symm(2,0) = 11.4795; symm(2,1) = -38.3720; symm(2,2) = -12.3612;
  cout << "symmetric matrix: " << symm << endl;
  eigen_vectors = IdentityMatrix(3,3);
  eigenDecomposition(symm,eigen_values,eigen_vectors);
}

void Test::productMatrixVector(void)
{
  Matrix m1(4,3);
  for (int i=0; i<4; i++) {
    for (int j=0; j<3; j++) {
      m1(i,j) = i + 1 + j/2.0;
    }
  }
  Vector v(3,0);
  v[0] = 1; v[1] = -2; v[2] = 3;
  cout << "m1: " << m1 << endl;
  cout << "v: "; print(cout,v,3); cout << endl;
  Vector ans = prod(m1,v);
  cout << "m1*v: "; print(cout,ans,3); cout << endl;

  Matrix m2(3,4);
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
  std::vector<Vector > sample = load_matrix(file_name);
  Vector unit_mean = computeVectorSum(sample);
  cout << "unit_mean: "; print(cout,unit_mean,3); cout << endl;
  Matrix m = computeDispersionMatrix(sample);
  cout << "dispersion: " << m << endl;

  Matrix eigen_vectors = IdentityMatrix(3,3);
  Vector eigen_values(3);
  eigenDecomposition(m,eigen_values,eigen_vectors);

  m(0,0) = 0.341; m(0,1) = -0.221; m(0,2) = 0.408;
  m(1,0) = -0.221; m(1,1) = 0.153; m(1,2) = -0.272;
  m(2,0) = 0.408; m(2,1) = -0.272; m(2,2) = 0.506;
  eigen_vectors = IdentityMatrix(3,3);
  eigenDecomposition(m,eigen_values,eigen_vectors);
}

void Test::numericalIntegration(void)
{
  // integration
  long double val = computeDawsonsIntegral(10);
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
  Vector spherical(3,0),cartesian(3,0);
  spherical[0] = 1;
  spherical[1] = PI/3;  // 60 degrees
  spherical[2] = 250 * PI/180;  // 250 degrees
  spherical2cartesian(spherical,cartesian);
  cout << "cartesian: "; print(cout,cartesian,3); cout << endl;

  boost::numeric::ublas::vector<long double> vec(3);
  for (int i=0; i<3; i++) vec(i) = cartesian[i];
  Matrix r1 = align_zaxis_with_vector(cartesian);
  cout << "r1: " << r1 << endl;
  Matrix inverse(3,3);
  invertMatrix(r1,inverse);
  cout << "inverse: " << inverse << endl;
  Matrix check = prod(r1,inverse);
  cout << "check: " << check << endl;

  boost::numeric::ublas::vector<long double> zaxis(3);
  for (int i=0; i<3; i++) zaxis(i) = ZAXIS[i];
  boost::numeric::ublas::vector<long double> ans1 = prod(r1,zaxis);
  cout << "ans1: " << ans1 << endl;
  boost::numeric::ublas::vector<long double> ans2 = prod(inverse,vec);
  cout << "ans2: " << ans2 << endl;
}

void Test::orthogonalTransformations2(void)
{
  Vector m0,m1,m2;
  generateRandomOrthogonalVectors(m0,m1,m2);
  cout << "mean: "; print(cout,m0,3); cout << endl;
  cout << "major: "; print(cout,m1,3); cout << endl;
  cout << "minor: "; print(cout,m2,3); cout << endl;

  Matrix r = computeOrthogonalTransformation(m0,m1);
  cout << "r: " << r << endl;
  Vector ztransform = prod(r,ZAXIS);
  cout << "ztransform: "; print(cout,ztransform,3); cout << endl; // = mu0
  Vector xtransform = prod(r,XAXIS);
  cout << "xtransform: "; print(cout,xtransform,3); cout << endl; // = mu1
  Vector ytransform = prod(r,YAXIS);
  cout << "ytransform: "; print(cout,ytransform,3); cout << endl; // = mu2
}

void Test::randomSampleGeneration(void)
{
  std::vector<Vector > random_sample;
  Vector m0,m1,m2;
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
  Vector mean(2,0); mean[0] = 1;
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
  long double kappa = 100;
  long double beta = 14.5;
  Vector m0 = ZAXIS;
  Vector m1 = XAXIS;
  Vector m2 = YAXIS;
  Vector kmu(3,0);
  for (int i=0; i<3; i++) kmu[i] = kappa * m0[i];
  cout << "kmu: "; print(cout,kmu,0); cout << endl;
  Matrix a1 = outer_prod(m1,m1);
  cout << "a1: " << a1 << endl;
  Matrix a2 = outer_prod(m2,m2);
  cout << "a2: " << a2 << endl;
  Matrix A = beta * (a1 - a2);
  cout << "A: " << A << endl;
  //Kent kent(100,30);
  Kent kent(kappa,beta);
  Kent::Constants constants = kent.getConstants();
  cout << "log_norm: " << constants.log_c << endl;
  cout << "dc_db: " << constants.log_cb << endl;
  cout << "dc_dk: " << constants.log_ck << endl;
  cout << "d2c_dk2: " << constants.log_ckk << endl;
  cout << "d2c_db2: " << constants.log_cbb << endl;
  cout << "d2c_dkdb: " << constants.log_ckb << endl;

  /*std::vector<Vector > random_sample;
  generateRandomOrthogonalVectors(m0,m1,m2);
  cout << "m0: "; print(cout,m0,0); cout << endl;
  cout << "m1: "; print(cout,m1,0); cout << endl;
  cout << "m2: "; print(cout,m2,0); cout << endl;
  for (int i=0; i<3; i++) kmu[i] = kappa * m0[i];
  cout << "kmu: "; print(cout,kmu,0); cout << endl;
  a1 = outer_prod(m1,m1);
  cout << "a1: " << a1 << endl;
  a2 = outer_prod(m2,m2);
  cout << "a2: " << a2 << endl;
  A = beta * (a1 - a2);
  cout << "A: " << A << endl;
  Kent kent1(m0,m1,m2,kappa,beta);
  log_norm = kent1.computeLogNormalizationConstant();
  cout << "log_norm: " << log_norm << endl;
  cout << "dc_db: " << kent1.log_dc_db() << endl;
  cout << "dc_dk: " << kent1.log_dc_dk() << endl;
  cout << "d2c_dk2: " << kent1.log_d2c_dk2() << endl;
  cout << "d2c_db2: " << kent1.log_d2c_db2() << endl;
  cout << "d2c_dkdb: " << kent1.log_d2c_dkdb() << endl;*/
}

void Test::optimization(void)
{
  Kent kent;
  struct Estimates estimates;
  Vector spherical(3,0),sample_mean(3,0);

  // Kent example from paper
  cout << "\nExample from paper:\n";
  kent= Kent(100,20);
  int N = 34;
  sample_mean[0] = 0.083; 
  sample_mean[1] = -0.959; 
  sample_mean[2] = 0.131;
  cartesian2spherical(sample_mean,spherical);
  cout << "m0: "; print(cout,sample_mean,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  Matrix S(3,3);
  S(0,0) = 0.045; S(0,1) = -0.075; S(0,2) = 0.014;
  S(1,0) = -0.075; S(1,1) = 0.921; S(1,2) = -0.122;
  S(2,0) = 0.014; S(2,1) = -0.122; S(2,2) = 0.034;
  for (int i=0; i<3; i++) {
    sample_mean[i] *= N;
    for (int j=0; j<3; j++) {
      S(i,j) *= N;
    }
  }
  estimates = kent.computeMomentEstimates(sample_mean,S,N);
}

void Test::moment_estimation(void)
{
  Vector spherical(3,0);
  struct Estimates estimates;
  std::vector<Vector> random_sample;
  Vector m0,m1,m2;
  long double kappa = 100;
  long double beta = 47.5;

  generateRandomOrthogonalVectors(m0,m1,m2);
  cartesian2spherical(m0,spherical);
  cout << "m0: "; print(cout,m0,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cartesian2spherical(m1,spherical);
  cout << "m1: "; print(cout,m1,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cartesian2spherical(m2,spherical);
  cout << "m2: "; print(cout,m2,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";

  Kent kent(m0,m1,m2,kappa,beta);
  random_sample = kent.generate(100);
  estimates = kent.computeMomentEstimates(random_sample);

  cartesian2spherical(estimates.mean,spherical);
  cout << "m0_est: "; print(cout,estimates.mean,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cartesian2spherical(estimates.major_axis,spherical);
  cout << "m1_est: "; print(cout,estimates.major_axis,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cartesian2spherical(estimates.minor_axis,spherical);
  cout << "m2_est: "; print(cout,estimates.minor_axis,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cout << "kappa_est: " << estimates.kappa << "; beta_est: " << estimates.beta << endl;

  // Kent example from paper
  cout << "\nExample from paper:\n";
  kent = Kent(100,20);
  Vector sample_mean(3,0);
  sample_mean[0] = 0.083; sample_mean[1] = -0.959; sample_mean[2] = 0.131;
  //sample_mean[0] = -0.959; sample_mean[1] = 0.131; sample_mean[2] = 0.083;
  cartesian2spherical(sample_mean,spherical);
  cout << "m0: "; print(cout,sample_mean,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  Matrix S(3,3);
  S(0,0) = 0.045; S(0,1) = -0.075; S(0,2) = 0.014;
  S(1,0) = -0.075; S(1,1) = 0.921; S(1,2) = -0.122;
  S(2,0) = 0.014; S(2,1) = -0.122; S(2,2) = 0.034;
  /*S(0,0) = 0.921; S(0,1) = -0.122; S(0,2) = -0.075;
  S(1,0) = -0.122; S(1,1) = 0.034; S(1,2) = 0.014;
  S(2,0) = -0.075; S(2,1) = 0.014; S(2,2) = 0.045;*/
  int N = 34;
  for (int i=0; i<3; i++) {
    sample_mean[i] *= N;
    for (int j=0; j<3; j++) {
      S(i,j) *= N;
    }
  }
  estimates = kent.computeMomentEstimates(sample_mean,S,N);
  cartesian2spherical(estimates.mean,spherical);
  cout << "m0_est: "; print(cout,estimates.mean,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cartesian2spherical(estimates.major_axis,spherical);
  cout << "m1_est: "; print(cout,estimates.major_axis,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cartesian2spherical(estimates.minor_axis,spherical);
  cout << "m2_est: "; print(cout,estimates.minor_axis,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cout << "kappa_est: " << estimates.kappa << "; beta_est: " << estimates.beta << endl;
}

void Test::ml_estimation(void)
{
  Vector spherical(3,0);
  struct Estimates estimates;
  std::vector<Vector> random_sample;
  Vector m0,m1,m2;
  long double kappa = 100;
  long double beta = 40;

  generateRandomOrthogonalVectors(m0,m1,m2);
  cartesian2spherical(m0,spherical);
  cout << "m0: "; print(cout,m0,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cartesian2spherical(m1,spherical);
  cout << "m1: "; print(cout,m1,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cartesian2spherical(m2,spherical);
  cout << "m2: "; print(cout,m2,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";

  Kent kent(m0,m1,m2,kappa,beta);
  random_sample = kent.generate(100);
  kent.computeAllEstimators(random_sample);

  // Kent example from paper
  cout << "\nExample from paper:\n";
  kent = Kent(100,20);
  Vector sample_mean(3,0);
  sample_mean[0] = 0.083; sample_mean[1] = -0.959; sample_mean[2] = 0.131;
  cartesian2spherical(sample_mean,spherical);
  cout << "m0: "; print(cout,sample_mean,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  Matrix S(3,3);
  S(0,0) = 0.045; S(0,1) = -0.075; S(0,2) = 0.014;
  S(1,0) = -0.075; S(1,1) = 0.921; S(1,2) = -0.122;
  S(2,0) = 0.014; S(2,1) = -0.122; S(2,2) = 0.034;
  int N = 34;
  for (int i=0; i<3; i++) {
    sample_mean[i] *= N;
    for (int j=0; j<3; j++) {
      S(i,j) *= N;
    }
  }
  kent.computeAllEstimators(sample_mean,S,N);
}

void Test::expectation()
{
  Vector m0,m1,m2;
  long double kappa = 100;
  long double beta = 47.5;

  generateRandomOrthogonalVectors(m0,m1,m2);
  Kent kent(m0,m1,m2,kappa,beta);
  kent.computeExpectation();
  Kent::Constants constants = kent.getConstants();
  cout << "E_x: "; print(cout,constants.E_x,0); cout << endl;
  cout << "E_xx: " << constants.E_xx << endl;
}

void Test::kl_divergence()
{
  Vector m0,m1,m2;
  long double kappa = 100;
  long double beta = 47.5;
  generateRandomOrthogonalVectors(m0,m1,m2);
  Kent kent1(m0,m1,m2,kappa,beta);

  kappa = 200; beta = 60;
  Kent kent2(m0,m1,m2,kappa,beta);
  cout << "KL-Div: " << kent1.computeKLDivergence(kent2) << endl;

  generateRandomOrthogonalVectors(m0,m1,m2);
  kappa = 100; beta = 47.5;
  Kent kent3(m0,m1,m2,kappa,beta);
  cout << "KL-Div: " << kent1.computeKLDivergence(kent3) << endl;
}

void Test::fisher()
{
  Vector m0,m1,m2;
  long double kappa = 100;
  long double beta = 30;
  long double psi,alpha,eta;

  generateRandomOrthogonalVectors(m0,m1,m2);
  Vector spherical(3,0);

  computeOrthogonalTransformation(m0,m1,psi,alpha,eta);
  cartesian2spherical(m0,spherical);
  cout << "m0: "; print(cout,spherical,0); cout << endl;

  cartesian2spherical(m1,spherical);
  cout << "m1: "; print(cout,spherical,0); cout << endl;

  Kent kent(psi,alpha,eta,kappa,beta);
  kent.computeExpectation();
  Kent::Constants constants = kent.getConstants();
  long double log_det_fkb = kent.computeLogFisherScale();
  cout << "log(det(f_kb)): " << log_det_fkb << endl;
  cout << "det(f_kb): " << exp(log_det_fkb) << endl;
  long double log_det_faxes = kent.computeLogFisherAxes();
  cout << "log(det(f_axes)): " << log_det_faxes << endl;
  cout << "det(f_axes): " << exp(log_det_faxes) << endl;
  cout << "log(det(fisher)): " << log_det_fkb+log_det_faxes << endl;
  cout << "log(fisher): " << kent.computeLogFisherInformation() << endl;
}

void Test::mml_estimation(void)
{
  Vector spherical(3,0);
  struct Estimates estimates;
  std::vector<Vector> random_sample;
  Vector m0,m1,m2;
  long double kappa = 100;
  long double beta = 47.5;
  int sample_size = 10;

  generateRandomOrthogonalVectors(m0,m1,m2);
  cartesian2spherical(m0,spherical);
  cout << "m0: "; print(cout,m0,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cartesian2spherical(m1,spherical);
  cout << "m1: "; print(cout,m1,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cartesian2spherical(m2,spherical);
  cout << "m2: "; print(cout,m2,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";

  Kent kent(m0,m1,m2,kappa,beta);
  random_sample = kent.generate(sample_size);
  //writeToFile("random_sample.dat",random_sample,3);
  kent.computeAllEstimators(random_sample);
}

