#include "ACG.h"
#include "Support.h"
#include "MultivariateNormal.h"

ACG::ACG()
{
  D = 3;
  W = IdentityMatrix(3,3);
}

ACG::ACG(Matrix &W) : W(W)
{
  D = W.size1();
  //assert(D == W.size2());
}

ACG ACG::operator=(const ACG &source)
{
  if (this != &source) {
    D = source.D;
    W = source.W;
  }
  return *this;
}

void ACG::printParameters()
{
  cout << "ACG(W): " << W << endl;
}

std::vector<Vector> ACG::generate(int N)
{
  // compute inverse of W
  Matrix cov(D,D);
  invertMatrix(W,cov);
  Vector mean(D,0);

  MultivariateNormal mvnorm(mean,cov);
  std::vector<Vector> y = mvnorm.generate(N);
  Vector x(D,0);
  std::vector<Vector> random_sample;
  for (int i=0; i<N; i++) {
    normalize(y[i],x);
    random_sample.push_back(x);
  }
  return random_sample;
}

