#include "Bingham.h"
#include "Support.h"
#include "ACG.h"

Bingham::Bingham()
{
  D = 3;
  A = IdentityMatrix(3,3);
}

Bingham::Bingham(Matrix &A) : A(A)
{
  D = A.size1();
}

Bingham Bingham::operator=(const Bingham &source)
{
  if (this != &source) {
    D = source.D;
    A = source.A;
  }
  return *this;
}

void Bingham::printParameters()
{
  cout << "Bingham(A): " << A << endl;
}

std::vector<Vector> Bingham::generate(int N)
{
  long double b = 1;
  Matrix id = IdentityMatrix(D,D);
  Matrix scaled = (2.0/b) * A;
  Matrix W = id + scaled;
  ACG acg(W);

  int num_samples = 0;
  long double logm = -0.5 * (D - b);
  logm += 0.5 * D * log(D/b);
  std::vector<Vector> random_sample,tmp;
  Vector x;
  long double u,check,prod1,prod2;
  for (int i=0; i<N; i++) {
    repeat:
    u = uniform_random();
    tmp = acg.generate(1);
    x = tmp[0];
    prod1 = prod_vMv(x,A);
    prod2 = prod_vMv(x,W);
    check = -prod1 - logm + ((0.5*D)*log(prod2));
    if (log(u) < check) {
      random_sample.push_back(x);
      cout << i+1 << endl;
    } else {
      goto repeat;
    }
  }
  return random_sample;
}

