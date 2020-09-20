#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;

long double testing_optional(long double x = 3.14)
{
  cout << "x^2: " << x * x << endl;
}

long double testing_optional2(int N, vector<long double> v = vector<long double>(N,1))
{
  long double sum = 0;
  for (int i=0; i<v.size(); i++) {
    sum += v[i];
  }
  cout << "sum: " << sum << endl;
}

int main(int argc, char **argv)
{
  testing_optional(5);
  testing_optional();
  testing_optional(15);
  cout << endl;
  vector<long double> v(10,6.5);
  testing_optional2(10,v);
}

