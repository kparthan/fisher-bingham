#ifndef SUPPORT_H
#define SUPPORT_H

#include "Header.h"

struct Parameters
{
  int test;                 // flag to test some modules
};

// general functions
struct Parameters parseCommandLineInput (int, char **); 
void Usage (const char *, options_description &);
bool checkFile(string &);
void writeToFile(const char *, std::vector<std::vector<long double> > &, int);
string extractName(string &);
void print(ostream &, std::vector<long double> &, int);
void convert2boostvector(std::vector<long double> &stlvec, boost_vector &vec);

int sign(long double);
long double exponent(long double, long double);
long double normalize(std::vector<long double> &, std::vector<long double> &);
void cartesian2spherical(std::vector<long double> &, std::vector<long double> &);
void spherical2cartesian(std::vector<long double> &, std::vector<long double> &);
long double computeDotProduct(std::vector<long double> &, std::vector<long double> &);
long double computeLogSurfaceAreaSphere(int);
long double logModifiedBesselFirstKind(long double, long double);

matrix<long double> align_zaxis_with_vector2(std::vector<long double> &);
matrix<long double> move_axes(std::vector<long double> &);
matrix<long double> computeOrthogonalTransformation(std::vector<long double> &, std::vector<long double> &);
matrix<long double> align_xaxis_with_major_axis(std::vector<long double> &);
matrix<long double> align_zaxis_with_vector(std::vector<long double> &);
void generateRandomOrthogonalVectors(std::vector<long double> &, std::vector<long double> &, std::vector<long double> &);
std::vector<std::vector<long double> > transform(std::vector<std::vector<long double> > &, matrix<long double> &);
bool invertMatrix(const matrix<long double> &, matrix<long double> &);
void eigenDecomposition(matrix<long double>, boost::numeric::ublas::vector<long double> &, matrix<long double> &);
void jacobiRotateMatrix(matrix<long double> &, matrix<long double> &, int, int);
void integrate_function(double);
void track(const state_type &, const double);
void rhs(const state_type &, state_type &, const double);

void TestFunctions(void);

#endif

