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

std::vector<Vector> Bingham::generateCanonical(Matrix &L, int N)
{
  long double b = 1;
  Matrix id = IdentityMatrix(D,D);
  Matrix scaled = (2.0/b) * L;
  Matrix W = id + scaled;
  ACG acg(W);

  long double logm = -0.5 * (D - b);
  logm += 0.5 * D * log(D/b);

  std::vector<Vector> canonical_sample;
  Vector x;
  long double u,check,prod1,prod2;
  int num_samples;
  std::vector<Vector> many_samples; 
  int accepted = 0;
  while (accepted != N) {
    num_samples = 2 * (N-accepted);
    many_samples = acg.generate(num_samples);
    for (int i=0; i<many_samples.size(); i++) {
      u = uniform_random();
      x = many_samples[i];
      prod1 = prod_vMv(x,L);
      prod2 = prod_vMv(x,W);
      check = -prod1 - logm + ((0.5*D)*log(prod2));
      if (log(u) < check) {
        canonical_sample.push_back(x);
        if (++accepted == N) {
          goto finish;
        }
      }
    } // for() ends ...
  } // while() ends ..
  finish:
  return canonical_sample;
}

std::vector<Vector> Bingham::generate(int N)
{
  Vector eigen_values(D,0);
  Matrix eigen_vectors = IdentityMatrix(D,D);
  eigenDecomposition(A,eigen_values,eigen_vectors);
  Vector sorted = sort(eigen_values); // sorted in increasing order
  Vector lambdas(D,0);
  for (int i=0; i<D-1; i++) {
    lambdas[i] = sorted[D-1-i] - sorted[0];
  }
  Matrix L = ZeroMatrix(D,D);
  for (int i=0; i<D-1; i++) {
    L(i,i) = lambdas[i];
  }
  
  std::vector<Vector> canonical_sample = generateCanonical(L,N);
  std::vector<Vector> random_sample(N);
  Matrix Vt = trans(eigen_vectors);
  for (int i=0; i<N; i++) {
    random_sample[i] = prod(canonical_sample[i],Vt);
  }
  return random_sample;
}

