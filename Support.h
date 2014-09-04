#ifndef SUPPORT_H
#define SUPPORT_H

#include "Header.h"

struct Parameters
{
  int test;                 // flag to test some modules
};

struct Estimates
{
  Vector mean,major_axis,minor_axis;
  long double kappa,beta;
};

// general functions
struct Parameters parseCommandLineInput (int, char **); 
void Usage (const char *, options_description &);
bool checkFile(string &);
void writeToFile(const char *, std::vector<Vector > &, int);
string extractName(string &);
void print(ostream &, Vector &, int);
void print(string &, struct Estimates &);

int sign(long double);
long double normalize(Vector &, Vector &);
long double norm(Vector &);
void cartesian2spherical(Vector &, Vector &);
void cartesian2sphericalPoleXAxis(Vector &, Vector &);
void spherical2cartesian(Vector &, Vector &);
long double computeDotProduct(Vector &, Vector &);
Vector crossProduct(Vector &, Vector &); 
long double computeLogSurfaceAreaSphere(int);
long double logModifiedBesselFirstKind(long double, long double);
void solveQuadratic(Vector &, long double, long double, long double);

std::vector<Vector> load_matrix(string &);
Matrix outer_prod(Vector &, Vector &);
Vector prod(Matrix &, Vector &);
Vector prod(Vector &, Matrix &);
long double prod_vMv(Vector &, Matrix &);
long double prod_xMy(Vector &, Matrix &, Vector &);
long double determinant(Matrix &);
Vector computeVectorSum(std::vector<Vector > &);
Vector computeNormalizedVectorSum(std::vector<Vector > &);
Matrix computeDispersionMatrix(std::vector<Vector > &);
Matrix computeNormalizedDispersionMatrix(std::vector<Vector > &);
Matrix rotate_about_yaxis(long double);
Matrix rotate_about_zaxis(long double);
Matrix computeOrthogonalTransformation(Vector &, Vector &);
Matrix computeOrthogonalTransformation(long double, long double, long double);
void computeOrthogonalTransformation(Vector &, Vector &, long double &, long double &, long double &);
Matrix align_zaxis_with_vector(Vector &);
Matrix align_vector_with_zaxis(Vector &);
void generateRandomOrthogonalVectors(Vector &, Vector &, Vector &);
std::vector<Vector> transform(std::vector<Vector > &, Matrix &);
bool invertMatrix(const Matrix &, Matrix &);
void eigenDecomposition(Matrix, Vector &, Matrix &);
void jacobiRotateMatrix(Matrix &, Matrix &, int, int);
long double computeDawsonsIntegral(double);
void track(const std::vector<double> &, const double);
void rhs(const std::vector<double> &, std::vector<double> &, const double);

void TestFunctions(void);

#endif

