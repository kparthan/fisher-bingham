#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double norm(vector<double> &v)
{
  double normsq = 0;
  for (int i=0; i<v.size(); i++) {
    normsq += v[i] * v[i];
  }
  return sqrt(normsq);
}

double normalize(vector<double> &x, vector<double> &unit)
{
  double l2norm = norm(x);
  for (int i=0; i<x.size(); i++) {
    unit[i] = x[i] / l2norm;
  }
  return l2norm;
}

void cartesian2spherical(vector<double> &cartesian, vector<double> &spherical)
{
  vector<double> unit(3,0);
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
      phi = 2 * M_PI - angle;
    }
  } else if (y >= 0) {
    phi = angle;
  } else if (y < 0) {
    phi = 2 * M_PI - angle;
  }

  spherical[0] = r;
  spherical[1] = theta;
  spherical[2] = phi;
}

void spherical2cartesian(vector<double> &spherical, vector<double> &cartesian)
{
  cartesian[0] = spherical[0] * sin(spherical[1]) * cos(spherical[2]);
  cartesian[1] = spherical[0] * sin(spherical[1]) * sin(spherical[2]);
  cartesian[2] = spherical[0] * cos(spherical[1]);
}

int main(int argc, char **argv)
{
  vector<double> x(3,0);
  x[0] = 0.364;
  x[1] = -0.339;
  x[2] = 0.868;
  vector<double> spherical(3,0);
  cartesian2spherical(x,spherical);
  cout << "theta: " << spherical[1] * 180/M_PI << endl;
  cout << "phi: " << spherical[2] * 180/M_PI << endl;
  return 0;
}

