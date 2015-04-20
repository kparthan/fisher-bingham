#include "Test.h"
#include "Support.h"
#include "Normal.h"
#include "vMC.h"
#include "vMF.h"
#include "FB4.h"
#include "FB6.h"
#include "Kent.h"
#include "Kent_EccTrans.h"
#include "Kent_UnifTrans.h"
#include "MultivariateNormal.h"
#include "ACG.h"
#include "Bingham.h"

extern Vector XAXIS,YAXIS,ZAXIS;
extern int ENABLE_DATA_PARALLELISM;
extern int NUM_THREADS;
extern int ESTIMATION,CRITERION;
extern int PRIOR;

void Test::load_data()
{
  string data_file = "random_sample.dat";
  std::vector<Vector> data = load_data_table(data_file);
  print(cout,data[10],3); cout << endl;
  string output = "copy.dat";
  writeToFile(output,data);
}

void Test::bessel()
{
  double d,k,log_bessel;
  
  d = 3;
  k = 730;

  log_bessel = computeLogModifiedBesselFirstKind(d,k);
  cout << "log I(" << d << "," << k << "): " << log_bessel << endl;

  log_bessel = log(cyl_bessel_i(d,k));
  cout << "Boost log I(" << d << "," << k << "): " << log_bessel << endl;
}

void Test::testing_cartesian2sphericalPoleXAxis()
{
  double theta,phi;
  Vector x(3,0),spherical(3,0);

  // in degrees
  theta = 90;
  phi = 100;

  // in radians
  theta *= PI/180;
  phi *= PI/180;

  x[0] = cos(theta);
  x[1] = sin(theta) * cos(phi);
  x[2] = sin(theta) * sin(phi);

  cartesian2spherical(x,spherical);

  cout << "Spherical coordinates: ";
  cout << "r: " << spherical[0] << "\n";
  cout << "theta: " << spherical[1] * 180/PI << "\n";
  cout << "phi: " << spherical[2] * 180/PI << "\n";
}

void Test::parallel_sum_computation(void)
{
  int N = 100;
  Vector ans;
  std::vector<Vector> sample;
  Vector weights(N,0);
  Vector spherical(3,1),x(3,0);
  double Neff;

  for (int i=0; i<N; i++) {
    spherical[1] = PI * uniform_random();
    spherical[2] = (2 * PI) * uniform_random();
    spherical2cartesian(spherical,x);
    sample.push_back(x);
    weights[i] = uniform_random();
  }
  // without parallelization ...
  ans = computeVectorSum(sample);
  cout << "sum: "; print(cout,ans,3); cout << endl;
  ans = computeVectorSum(sample,weights,Neff);
  cout << "sum: "; print(cout,ans,3); cout << endl;
  cout << "Neff: " << Neff << endl;

  // with parallelization ...
  ENABLE_DATA_PARALLELISM = SET;
  NUM_THREADS = 42;
  ans = computeVectorSum(sample);
  cout << "sum(parallel): "; print(cout,ans,3); cout << endl;
  ans = computeVectorSum(sample,weights,Neff);
  cout << "sum: "; print(cout,ans,3); cout << endl;
  cout << "Neff: " << Neff << endl;
}

void Test::uniform_number_generation()
{
  for (int i=0; i<10; i++) {
    cout << uniform_random() << endl;
  }
}

void Test::arbitrary_rotation()
{
  double psi,alpha,eta,theta;
  Vector mu,major_axis,minor_axis;
  Matrix R,Rx,Ry,Rz;

  psi = 60;
  alpha = 90;
  eta = 45;

  psi *= PI/180;
  alpha *= PI/180;
  eta *= PI/180;

  Matrix r = computeOrthogonalTransformation(psi,alpha,eta);
  mu = Vector(3,0); 
  major_axis = Vector(3,0); 
  minor_axis = Vector(3,0); 
  for (int i=0; i<3; i++) {
    mu[i] = r(i,0);
    major_axis[i] = r(i,1);
    minor_axis[i] = r(i,2);
  }
  //cout << "r: " << r << endl;

  theta = 50;
  theta *= PI/180;

  Rx = rotate_about_xaxis(theta);
  cout << "Rx: " << Rx << endl;
  R = rotate_about_arbitrary_axis(XAXIS,theta);
  cout << "R: " << R << endl;

  Ry = rotate_about_yaxis(theta);
  cout << "Ry: " << Ry << endl;
  R = rotate_about_arbitrary_axis(YAXIS,theta);
  cout << "R: " << R << endl;

  Rz = rotate_about_zaxis(theta);
  cout << "Rz: " << Rz << endl;
  R = rotate_about_arbitrary_axis(ZAXIS,theta);
  cout << "R: " << R << endl;
}

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
  boost::numeric::ublas::vector<double> v(3);
  for (int i=0; i<3; i++) {
    v[i] = i + 3;
  }
  cout << "v: " << v << endl;
  cout << "2 * v: " << 2 * v << endl;
  cout << "v/2: " << v/2 << endl;

  // adding matrix row and vector
  boost::numeric::ublas::vector<double> v1 = mr + v;
  cout << "v1: " << v1 << endl;

  // multiplication
  Matrix m3 = prod(m1,m2);
  cout << "m1 * m2 = m3: " << m3 << endl;

  // multiplying matrices and vectors
  boost::numeric::ublas::vector<double> mv = prod(v1,m1);
  cout << "v1 * m1 = mv: " << mv << endl;
  mv = prod(m1,v1);
  cout << "m1 * v1 = mv: " << mv << endl;

  double v1_T_v1 = inner_prod(v1,v1);
  cout << "v1' * v1 = : " << v1_T_v1 << endl;

  Matrix m4 = outer_prod(v1,v1);
  cout << "v1 * v1' = m4: " << m4 << endl;

  // eigen values & vectors
  Matrix symm = (m1 + trans(m1))/2;
  cout << "symmetric matrix: " << symm << endl;
  Vector eigen_values(3,0);
  Matrix eigen_vectors = IdentityMatrix(3,3);
  eigenDecomposition(symm,eigen_values,eigen_vectors);
  Matrix diag = IdentityMatrix(3,3);
  for (int i=0; i<3; i++) diag(i,i) = eigen_values[i];
  Matrix tmp = prod(eigen_vectors,diag);
  Matrix trans_eigenvec = trans(eigen_vectors);
  Matrix check = prod(tmp,trans_eigenvec);
  cout << "check (V * diag * V'): " << check << endl << endl;
  
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
  diag = IdentityMatrix(3,3);
  for (int i=0; i<3; i++) diag(i,i) = eigen_values[i];
  tmp = prod(eigen_vectors,diag);
  trans_eigenvec = trans(eigen_vectors);
  check = prod(tmp,trans_eigenvec);
  cout << "check (V * diag * V'): " << check << endl << endl;

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
  string file_name = "./visualize/sampled_data/kent.dat";
  std::vector<Vector> sample = load_data_table(file_name);
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
  double val = computeDawsonsIntegral(10);
}

void Test::normalDistributionFunctions(void)
{
  // Normal cdf
  Normal normal(0,1);
  double cdf = normal.cumulativeDensity(2);
  cout << "cdf: " << cdf << endl;

  double x = sqrt(PI/2);
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

  boost::numeric::ublas::vector<double> vec(3);
  for (int i=0; i<3; i++) vec(i) = cartesian[i];
  Matrix r1 = align_xaxis_with_vector(cartesian);
  cout << "r1: " << r1 << endl;
  Matrix inverse(3,3);
  invertMatrix(r1,inverse);
  cout << "inverse: " << inverse << endl;
  Matrix check = prod(r1,inverse);
  cout << "check: " << check << endl;

  boost::numeric::ublas::vector<double> xaxis(3);
  for (int i=0; i<3; i++) xaxis(i) = XAXIS[i];
  boost::numeric::ublas::vector<double> ans1 = prod(r1,xaxis);
  cout << "ans1: " << ans1 << endl;
  boost::numeric::ublas::vector<double> ans2 = prod(inverse,vec);
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
  Vector xtransform = prod(r,XAXIS);
  cout << "xtransform: "; print(cout,xtransform,3); cout << endl; // = mu0
  Vector ytransform = prod(r,YAXIS);
  cout << "ytransform: "; print(cout,ytransform,3); cout << endl; // = mu1
  Vector ztransform = prod(r,ZAXIS);
  cout << "ztransform: "; print(cout,ztransform,3); cout << endl; // = mu2
}

void Test::randomSampleGeneration(void)
{
  std::vector<Vector> random_sample;
  Vector m0,m1,m2;
  generateRandomOrthogonalVectors(m0,m1,m2);

  // FB4 generation (gamma < 0)
  FB4 fb4_1(m0,m1,m2,100,-10);
  random_sample = fb4_1.generate(1000);
  writeToFile("./visualize/sampled_data/fb4_1.dat",random_sample,3);

  // FB4 generation (gamma < 0)
  FB4 fb4_2(m0,m1,m2,100,10);
  random_sample = fb4_2.generate(1000);
  writeToFile("./visualize/sampled_data/fb4_2.dat",random_sample,3);

  // vMF generation
  vMF vmf(XAXIS,1000);
  random_sample = vmf.generate(10000);
  writeToFile("./visualize/sampled_data/vmf.dat",random_sample,3);

  // vMF (2D)
  Vector mean(2,0); mean[0] = 1;
  vMC vmc(mean,10);
  vmc.generateCanonical(random_sample,100);
  writeToFile("./visualize/sampled_data/vmc.dat",random_sample,3);

  // FB6 generation 
  FB6 fb6(m0,m1,m2,100,15,-10);
  random_sample = fb6.generate(1000);
  writeToFile("./visualize/sampled_data/fb6.dat",random_sample,3);

  // FB6 generation 
  FB6 kent(m0,m1,m2,1000,475,0);
  random_sample = kent.generate(1000);
  writeToFile("./visualize/sampled_data/kent.dat",random_sample,3);

  // FB6 generation 
  /*FB6 kent1(m0,m1,m2,1,60,-10);
  random_sample = kent1.generate(1000);
  writeToFile("./visualize/sampled_data/kent1.dat",random_sample,3);*/
}

void Test::multivariate_normal(void)
{
  MultivariateNormal mvnrm;
  //mvnorm.printParameters();

  Vector mean;
  Matrix cov;
  int D;
  int N = 10000;
  std::vector<Vector> random_sample;

  // 2 D
  D = 2;
  mean = Vector(D,0);
  //cov = generateRandomCovarianceMatrix(D);
  cov = IdentityMatrix(2,2); cov(0,0) = 1; cov(1,1) = 10;
  MultivariateNormal mvnorm2d(mean,cov);
  mvnorm2d.printParameters();
  random_sample = mvnorm2d.generate(N);
  writeToFile("./visualize/sampled_data/mvnorm2d.dat",random_sample,3);

  // 3 D
  D = 3;
  mean = Vector(D,0);
  //cov = generateRandomCovarianceMatrix(D);
  cov = IdentityMatrix(3,3); cov(0,0) = 1; cov(1,1) = 10; cov(2,2) = 100;
  MultivariateNormal mvnorm3d(mean,cov);
  mvnorm3d.printParameters();
  random_sample = mvnorm3d.generate(N);
  writeToFile("./visualize/sampled_data/mvnorm3d.dat",random_sample,3);
}

void Test::acg(void)
{
  Vector mean;
  Matrix cov;
  int D;
  int N = 10000;
  std::vector<Vector> random_sample;

  // 3 D
  D = 3;
  mean = Vector(D,0);
  //cov = generateRandomCovarianceMatrix(D);
  cov = IdentityMatrix(3,3); cov(0,0) = 1; cov(1,1) = 10; cov(2,2) = 1000;
  Matrix W(D,D);
  invertMatrix(cov,W);
  ACG acg(W);
  acg.printParameters();
  random_sample = acg.generate(N);
  writeToFile("./visualize/sampled_data/acg.dat",random_sample,3);
}

void Test::bingham(void)
{
  int D = 3;
  int N = 10000;
  double beta = 47.5;
  std::vector<Vector> random_sample;

  Vector mean,major,minor;
  generateRandomOrthogonalVectors(mean,major,minor);
  Matrix a1 = outer_prod(major,major);
  Matrix a2 = outer_prod(minor,minor);
  Matrix tmp = a2 - a1;
  Matrix A = (beta * tmp);

  /*Vector eigen_values(3,0);
  Matrix eigen_vectors = IdentityMatrix(3,3);
  eigenDecomposition(A,eigen_values,eigen_vectors);
  Vector sorted = sort(eigen_values);
  cout << "eig: "; print(cout,eigen_values,3); cout << endl;
  cout << "sorted: "; print(cout,sorted,3); cout << endl;*/
  Bingham bingham(A);
  bingham.printParameters();
  random_sample = bingham.generate(N);
  writeToFile("./visualize/sampled_data/bingham.dat",random_sample,3);
}

void Test::kent_bingham_generation(void)
{
  double ecc,kappa,beta;
  std::vector<Vector> random_sample;
  Vector m0,m1,m2;
  //generateRandomOrthogonalVectors(m0,m1,m2);

  m0 = ZAXIS; m1 = XAXIS; m2 = YAXIS;
  int N = 10000;
  kappa = 1000;
  ecc = 0.99;
  beta = 0.5 * kappa * ecc;
  Kent kent(m0,m1,m2,kappa,beta);
  random_sample = kent.generate(N);
  writeToFile("./visualize/sampled_data/kent.dat",random_sample,3);
}

void Test::normalization_constant(void)
{
  cout << "ZERO: " << ZERO << endl;
  double kappa = 100;
  double beta = 14.5;
  Vector m0 = XAXIS;
  Vector m1 = YAXIS;
  Vector m2 = ZAXIS;
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

  /*std::vector<Vector> random_sample;
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
  double kappa = 100;
  //double beta = 45;
  double beta = 47.5;

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
  random_sample = kent.generate(1000);
  writeToFile("./visualize/sampled_data/kent.dat",random_sample,3);
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

  // Kent example from paper
  cout << "\nReading Whin Sill data ...\n";
  string file_name = "./support/R_codes/whin_sill.txt";
  std::vector<Vector> whin_sill = load_data_table(file_name);
  estimates = kent.computeMomentEstimates(whin_sill);

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
  std::vector<struct Estimates> all_estimates;
  std::vector<Vector> random_sample;
  Vector m0,m1,m2;
  double kappa = 100;
  double beta = 40;

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
  kent.computeAllEstimators(random_sample,all_estimates,1,1);

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
  kent.computeAllEstimators(sample_mean,S,N,all_estimates,1,0);
}

void Test::expectation()
{
  Vector m0,m1,m2;
  double kappa = 100;
  double beta = 47.5;

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
  double kappa = 100;
  double beta = 47.5;
  generateRandomOrthogonalVectors(m0,m1,m2);
  Kent kent1(m0,m1,m2,kappa,beta);

  kappa = 200; beta = 60;
  Kent kent2(m0,m1,m2,kappa,beta);
  cout << "KL-Div: " << kent1.computeKLDivergence(kent2) << endl;

  //generateRandomOrthogonalVectors(m0,m1,m2);
  kappa = 100; beta = 47.5;
  Kent kent3(m0,m1,m2,kappa,beta);
  cout << "KL-Div: " << kent1.computeKLDivergence(kent3) << endl;
}

void Test::fisher()
{
  Vector m0,m1,m2;
  double kappa = 100;
  double beta = 30;
  double psi,alpha,eta;

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
  double log_det_fkb = kent.computeLogFisherScale();
  cout << "log(det(f_kb)): " << log_det_fkb << endl;
  cout << "det(f_kb): " << exp(log_det_fkb) << endl;
  double log_det_faxes = kent.computeLogFisherAxes();
  cout << "log(det(f_axes)): " << log_det_faxes << endl;
  cout << "det(f_axes): " << exp(log_det_faxes) << endl;
  cout << "log(det(fisher)): " << log_det_fkb+log_det_faxes << endl;
  cout << "log(fisher): " << kent.computeLogFisherInformation() << endl;
}

void Test::fisher2()
{
  double kappa,beta,ecc;
  double log_det_fkb,log_det_faxes,log_det_fisher;
  double log_latt,log_prior_axes,log_prior_scale,log_prior;
  double ip1,ip2,ip;
  double Neff = 1;

  kappa = 1;
  ecc = TOLERANCE;
  cout << "kappa: " << kappa << endl << endl;
  while (ecc < 0.95) {
    beta = 0.5 * kappa * ecc;
    cout << "(ecc,beta): (" << ecc << ", " << beta << ")\n";

    Kent kent(ZAXIS,XAXIS,YAXIS,kappa,beta);
    kent.computeExpectation();
    Kent::Constants constants = kent.getConstants();

    log_latt = 0;
    log_prior_scale = kent.computeLogPriorScale();
    log_det_fkb = kent.computeLogFisherScale();
    //log_latt = logLatticeConstant(2);
    //cout << "log(latt_2): " << log_latt << endl;
    cout << "log(prior_scale): " << log_prior_scale << endl;
    cout << "log(det(f_kb)): " << log_det_fkb << endl;
    ip1 = -log_prior_scale + 0.5 * log_det_fkb + log_latt;
    cout << "ip1: " << ip1 << endl;

    //log_latt = logLatticeConstant(3);
    //cout << "log(latt_3): " << log_latt << endl;
    log_prior_axes = kent.computeLogPriorAxes();
    log_det_faxes = kent.computeLogFisherAxes();
    cout << "log(prior_axes): " << log_prior_axes << endl;
    cout << "log(det(f_axes)): " << log_det_faxes << endl;
    ip2 = -log_prior_axes + 0.5 * log_det_faxes + 1.5 * log_latt;
    cout << "ip2: " << ip2 << endl;

    log_latt = logLatticeConstant(5);
    cout << "log(latt_5): " << log_latt << endl;
    log_prior = kent.computeLogPriorProbability();
    log_det_fisher = kent.computeLogFisherInformation(1);
    cout << "log(prior): " << log_prior << endl;
    cout << "log(det(fisher)): " << log_det_fisher << endl;
    ip = -log_prior + 0.5 * log_det_fisher + 2.5 * log_latt;
    cout << "ip: " << ip << endl;

    Neff = 25;
    ip += 2.5 * log(Neff);
    cout << "ip: " << ip;

    cout << "\n\n";

    ecc += 0.005;
  }
}

void Test::mml_estimation(void)
{
  Vector spherical(3,0);
  std::vector<struct Estimates> all_estimates;
  std::vector<Vector> random_sample;
  Vector m0,m1,m2;
  double kappa = 100;
  double beta;
  int sample_size = 10;
  string data_file = "random_sample.dat";

  beta = 45;
  //beta = 2;

  //generateRandomOrthogonalVectors(m0,m1,m2);
  m0 = ZAXIS;
  m1 = XAXIS;
  m2 = YAXIS;
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
  writeToFile("random_sample.dat",random_sample,3);
  random_sample = load_data_table(data_file);
  kent.computeAllEstimators(random_sample,all_estimates,1,1);
}

void Test::mml_estimation2(void)
{
  double psi,alpha,eta;
  std::vector<struct Estimates> all_estimates;
  std::vector<Vector> random_sample;
  double kappa,beta,ecc;
  int sample_size = 10;
  string data_file = "random_sample.dat";

  kappa = 10;
  ecc = 0.5;
  beta = 0.5 * kappa * ecc;

  // in degrees
  psi = 90;
  alpha = 90;
  eta = 90;

  psi *= PI/180;
  alpha *= PI/180;
  eta *= PI/180;
  Kent kent(psi,alpha,eta,kappa,beta);
  //Kent kent(ZAXIS,XAXIS,YAXIS,kappa,beta);
  random_sample = kent.generate(sample_size);
  writeToFile("random_sample.dat",random_sample,3);

  //data_file = "./support/R_codes/whin_sill.txt";

  random_sample = load_data_table(data_file);
  kent.computeAllEstimators(random_sample,all_estimates,1,1);
  //kent.computeAllEstimators(random_sample,all_estimates,1,0);

  Vector statistics,pvalues;
  chisquare_hypothesis_testing(random_sample,all_estimates,statistics,pvalues);
}

void Test::plot_posterior_density(void)
{
  double psi,alpha,eta;
  string data_file = "./visualize/sampled_data/random_sample_ex1.dat";
  std::vector<Vector> random_sample = load_data_table(data_file);

  double kappa_inc = 0.5;
  double kappa_max = 30;
  double ecc_inc = 0.02;
  double ecc_max = 0.99;

/*
  string posterior_file = "./visualize/sampled_data/prior3d_posterior1.dat";
  ofstream out1(posterior_file.c_str());
  posterior_file = "./visualize/sampled_data/prior3d_posterior2.dat";
  ofstream out2(posterior_file.c_str());
  for (double k=kappa_inc; k<=kappa_max; k+=kappa_inc) {
    for (double e=ecc_inc; e<=ecc_max; e+=ecc_inc) {
      //psi = 63.873; alpha = 95.926; eta = 91.623;
      psi = 118.632; alpha = 85.533; eta = 87.221;
      psi *= PI/180; alpha *= PI/180; eta *= PI/180;
      double b = 0.5 * k * e;
      Kent kent1(psi,alpha,eta,k,b);
      double log_prior = kent1.computeLogPriorProbability();
      double fval = -log_prior + kent1.computeNegativeLogLikelihood(random_sample);
      double posterior = exp(-fval);
      out1 << fixed << scientific << setprecision(6) << k << "\t";
      out1 << fixed << scientific << setprecision(6) << b << "\t";
      out1 << fixed << scientific << setprecision(6) << posterior << "\t";
      out1 << endl;
      
      // eccentricity transform ...
      //psi = 63.869; alpha = 95.907; eta = 91.537;
      psi = 118.634; alpha = 85.547; eta = 87.213;
      psi *= PI/180; alpha *= PI/180; eta *= PI/180;
      Kent_EccTrans kent2(psi,alpha,eta,k,e);
      log_prior = kent2.computeLogPriorProbability();
      fval = -log_prior + kent2.computeNegativeLogLikelihood(random_sample);
      posterior = exp(-fval);
      out2 << fixed << scientific << setprecision(6) << k << "\t";
      out2 << fixed << scientific << setprecision(6) << e << "\t";
      out2 << fixed << scientific << setprecision(6) << posterior << "\t";
      out2 << endl;
    }
  }
  out1.close();
  out2.close();
*/

  PRIOR = 2;
  string posterior_file = "./visualize/sampled_data/prior2d_posterior1.dat";
  ofstream out1(posterior_file.c_str());
  posterior_file = "./visualize/sampled_data/prior2d_posterior2.dat";
  ofstream out2(posterior_file.c_str());
  for (double k=kappa_inc; k<=kappa_max; k+=kappa_inc) {
    for (double e=ecc_inc; e<=ecc_max; e+=ecc_inc) {
      psi = 118.629; alpha = 85.544; eta = 87.222;
      psi *= PI/180; alpha *= PI/180; eta *= PI/180;
      double b = 0.5 * k * e;
      Kent kent1(psi,alpha,eta,k,b);
      double log_prior = kent1.computeLogPriorProbability();
      double fval = -log_prior + kent1.computeNegativeLogLikelihood(random_sample);
      double posterior = exp(-fval);
      out1 << fixed << scientific << setprecision(6) << k << "\t";
      out1 << fixed << scientific << setprecision(6) << b << "\t";
      out1 << fixed << scientific << setprecision(6) << posterior << "\t";
      out1 << endl;
      
      // eccentricity transform ...
      psi = 118.626; alpha = 85.541; eta = 87.204;
      psi *= PI/180; alpha *= PI/180; eta *= PI/180;
      Kent_EccTrans kent2(psi,alpha,eta,k,e);
      log_prior = kent2.computeLogPriorProbability();
      fval = -log_prior + kent2.computeNegativeLogLikelihood(random_sample);
      posterior = exp(-fval);
      out2 << fixed << scientific << setprecision(6) << k << "\t";
      out2 << fixed << scientific << setprecision(6) << e << "\t";
      out2 << fixed << scientific << setprecision(6) << posterior << "\t";
      out2 << endl;
    } // e
    out1 << endl;
    out2 << endl;
  } // k
  out1.close();
  out2.close();

  posterior_file = "./visualize/sampled_data/prior2d_posterior3.dat";
  ofstream out3(posterior_file.c_str());
  double z4_inc = 0.02;
  double z4_max = 0.99;
  double z5_inc = 0.02;
  double z5_max = 0.99;
  // uniform transform ...
  psi = 118.644; alpha = 85.538; eta = 87.207;
  psi *= PI/180; alpha *= PI/180; eta *= PI/180;
  double z1 = psi / PI; 
  double z2 = 0.5 * (1 - cos(alpha)); 
  double z3 = eta / (2 * PI); 
  for (double z4=z4_inc; z4<=z4_max; z4+=z4_inc) {
    //if (z4>=0.75) z4_inc = 0.001;
    for (double z5=z5_inc; z5<=z5_max; z5+=z5_inc) {
      Kent_UnifTrans kent3(z1,z2,z3,z4,z5);
      double fval = kent3.computeNegativeLogLikelihood(random_sample);
      double posterior = exp(-fval);
      out3 << fixed << scientific << setprecision(6) << z4 << "\t";
      out3 << fixed << scientific << setprecision(6) << z5 << "\t";
      out3 << fixed << scientific << setprecision(6) << posterior << "\t";
      out3 << endl;
    } // z4
    out3 << endl;
  } // z5
  out3.close();

}

void Test::vmf_all_estimation()
{
  Vector spherical(3,1);
  spherical[1] = uniform_random() * PI;
  spherical[2] = uniform_random() * 2 * PI;

  Vector mean(3,0);
  spherical2cartesian(spherical,mean);
  double kappa = 100;
  int sample_size = 10;

  vMF vmf(mean,kappa);
  std::vector<Vector> random_sample = vmf.generate(sample_size);
  writeToFile("random_sample.dat",random_sample,3);
  vmf.computeAllEstimators(random_sample);
}

void Test::chi_square()
{
  int df = 500499;
  chi_squared chisq(df);
  double alpha = 0.05;
  double x = quantile(chisq,1-alpha);
  cout << "quantile: " << x << endl;

  double p = 1 - cdf(chisq,x);
  cout << "pvalue: " << p << endl;
}

void Test::hypothesis_testing()
{
  int N;
  double kappa,beta;
  std::vector<Vector> random_sample;
  Vector m0,m1,m2;
  generateRandomOrthogonalVectors(m0,m1,m2);

  N = 20;

  // Generating data from Kent
  kappa = 100;
  beta = 40;
  Kent kent(m0,m1,m2,kappa,beta);
  random_sample = kent.generate(N);
  writeToFile("./visualize/sampled_data/kent.dat",random_sample,3);

  /*string file_name = "./support/R_codes/whin_sill.txt";
  random_sample = load_data_table(file_name);*/
  kent.computeTestStatistic_vMF(random_sample);

  // Generating data from vMF
  /*Vector spherical(3,1);
  spherical[1] = PI * uniform_random();
  spherical[2] = 2 * PI * uniform_random();
  spherical2cartesian(spherical,m0);
  vMF vmf(m0,kappa);
  random_sample = vmf.generate(N);
  writeToFile("random_sample.dat",random_sample,3);
  kent.computeTestStatistic_vMF(random_sample);*/
}

void Test::hypothesis_testing2()
{
  string data_file = "./support/R_codes/whin_sill.txt";
  std::vector<Vector> data = load_data_table(data_file);

  cout << "vMF estimates ...\n";
  std::vector<struct Estimates_vMF> all_vmf_estimates;
  vMF vmf;
  vmf.computeAllEstimators(data,all_vmf_estimates);

  cout << "Kent estimates ...\n";
  std::vector<struct Estimates> all_kent_estimates;
  Kent kent;
  kent.computeAllEstimators(data,all_kent_estimates,1,0);

  cout << "Hypothesis testing ...\n";
  vMF vmf_ml(all_vmf_estimates[MLE].mean,all_vmf_estimates[MLE].kappa);
  double vmf_negloglike = vmf_ml.computeNegativeLogLikelihood(data);
  cout << "vMF (MLE) negloglike: " << vmf_negloglike << endl;

  Kent kent_ml(
    all_kent_estimates[MLE].psi,
    all_kent_estimates[MLE].alpha,
    all_kent_estimates[MLE].eta,
    all_kent_estimates[MLE].kappa,
    all_kent_estimates[MLE].beta
  );
  double kent_ml_negloglike = kent_ml.computeNegativeLogLikelihood(data);
  cout << "Kent (MLE) negloglike: " << kent_ml_negloglike << endl;

  // null: vMF
  double log_ratio_statistic = 2 * (vmf_negloglike - kent_ml_negloglike); // -2 log(ratio)
  cout << "log_ratio_statistic (null:vMF): " << log_ratio_statistic << endl;

  Vector statistics,pvalues;
  chisquare_hypothesis_testing(data,all_kent_estimates,statistics,pvalues);
/*
  // null: Kent(MOMENT)
  Kent kent_mom(
    all_kent_estimates[MOMENT].psi,
    all_kent_estimates[MOMENT].alpha,
    all_kent_estimates[MOMENT].eta,
    all_kent_estimates[MOMENT].kappa,
    all_kent_estimates[MOMENT].beta
  );
  double kent_mom_negloglike = kent_mom.computeNegativeLogLikelihood(data);
  cout << "Kent (MOMENT) negloglike: " << kent_mom_negloglike << endl;
  log_ratio_statistic = 2 * (kent_mom_negloglike - kent_ml_negloglike);
  cout << "log_ratio_statistic (null: Kent[MOMENT]): " << log_ratio_statistic << endl;

  // null: Kent(MAP)
  Kent kent_map(
    all_kent_estimates[MAP].psi,
    all_kent_estimates[MAP].alpha,
    all_kent_estimates[MAP].eta,
    all_kent_estimates[MAP].kappa,
    all_kent_estimates[MAP].beta
  );
  double kent_map_negloglike = kent_map.computeNegativeLogLikelihood(data);
  cout << "Kent (MAP) negloglike: " << kent_map_negloglike << endl;
  log_ratio_statistic = 2 * (kent_map_negloglike - kent_ml_negloglike);
  cout << "log_ratio_statistic (null: Kent[MAP]): " << log_ratio_statistic << endl;

  // null: Kent(MML)
  Kent kent_mml(
    all_kent_estimates[MML].psi,
    all_kent_estimates[MML].alpha,
    all_kent_estimates[MML].eta,
    all_kent_estimates[MML].kappa,
    all_kent_estimates[MML].beta
  );
  double kent_mml_negloglike = kent_mml.computeNegativeLogLikelihood(data);
  cout << "Kent (MML) negloglike: " << kent_mml_negloglike << endl;
  log_ratio_statistic = 2 * (kent_mml_negloglike - kent_ml_negloglike);
  cout << "log_ratio_statistic (null: Kent[MML]): " << log_ratio_statistic << endl;
*/
}

void Test::confidence_region()
{
  int N;
  double kappa,beta;
  std::vector<Vector> random_sample;
  Vector m0,m1,m2;
  generateRandomOrthogonalVectors(m0,m1,m2);

  N = 100;

  // Generating data from Kent
  kappa = 100;
  beta = 47.5;
  Kent kent(m0,m1,m2,kappa,beta);
  random_sample = kent.generate(N);
  writeToFile("./visualize/sampled_data/kent.dat",random_sample,3);

  /*string file_name = "./support/R_codes/whin_sill.txt";
  random_sample = load_data_table(file_name);*/
  kent.computeConfidenceRegion(random_sample);
}

void Test::infer_mixture()
{
  int N;
  double alpha,alpha_rad,sin_alpha,cos_alpha;
  double kappa,beta,ecc;
  Vector mean,major,minor;

  N = 10000;
  kappa = 100;
  ecc = 0.5;
  alpha = 60;

  alpha_rad = alpha * PI / 180;

  mean = ZAXIS;
  major = XAXIS;
  minor = YAXIS;
  beta = 0.5 * kappa * ecc;
  Kent kent1(mean,major,minor,kappa,beta);

  sin_alpha = sin(alpha_rad);
  cos_alpha = cos(alpha_rad);
  mean[0] = sin_alpha; mean[1] = 0; mean[2] = cos_alpha;
  major[0] = cos_alpha; major[1] = 0; major[2] = -sin_alpha;
  beta = 0.5 * kappa * ecc;
  Kent kent2(mean,major,minor,kappa,beta);

  Vector weights(2,0.5);

  std::vector<Kent> components;
  components.push_back(kent1);
  components.push_back(kent2);
  Mixture original(2,components,weights);

  //std::vector<Vector> data = original.generate(N,1);
  std::vector<Vector> large_data = original.generate(100000,0);
  string data_file = "random_sample.dat";
  std::vector<Vector> data = load_data_table(data_file);

  string log_file = "infer_mml.log";
  ESTIMATION = MML; CRITERION = MMLC;
  Mixture inferred = inferComponents(data,log_file);
  double kldiv_mml = original.computeKLDivergence(inferred,large_data);

  log_file = "infer_map.log";
  ESTIMATION = MAP; CRITERION = BIC;
  inferred = inferComponents(data,log_file);
  double kldiv_map = original.computeKLDivergence(inferred,large_data);

  cout << "kldiv_mml: " << kldiv_mml << endl;
  cout << "kldiv_map: " << kldiv_map << endl;
}

void Test::infer_mixture_vmf()
{
  int N;
  double kappa,alpha,alpha_rad,sin_alpha,cos_alpha;
  Vector mean;

  N = 10000;
  kappa = 100;
  alpha = 60;

  alpha_rad = alpha * PI / 180;

  mean = ZAXIS;
  vMF vmf1(mean,kappa);

  sin_alpha = sin(alpha_rad);
  cos_alpha = cos(alpha_rad);
  mean[0] = sin_alpha; mean[1] = 0; mean[2] = cos_alpha;
  vMF vmf2(mean,kappa);

  Vector weights(2,0.5);

  std::vector<vMF> components;
  components.push_back(vmf1);
  components.push_back(vmf2);
  Mixture_vMF original(2,components,weights);

  std::vector<Vector> data;
  data = original.generate(N,1);
  std::vector<Vector> large_data = original.generate(100000,0);
  string data_file = "random_sample.dat";
  data = load_data_table(data_file);

  string log_file = "infer_mml.log";
  ESTIMATION = MML; CRITERION = MMLC;
  Mixture_vMF inferred = inferComponents_vMF(data,log_file);
  double kldiv_mml = original.computeKLDivergence(inferred,large_data);

  log_file = "infer_map.log";
  ESTIMATION = MAP; CRITERION = BIC;
  inferred = inferComponents_vMF(data,log_file);
  double kldiv_map = original.computeKLDivergence(inferred,large_data);

  cout << "kldiv_mml: " << kldiv_mml << endl;
  cout << "kldiv_map: " << kldiv_map << endl;
}



