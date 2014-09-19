#ifndef SUPPORT_H
#define SUPPORT_H

#include "Mixture.h"

struct Parameters
{
  int test;                 // flag to test some modules
  int experiments;          // flag to run some experiments 
  int iterations;           // number of iterations
  string profile_file;      // path to a single profile
  string profiles_dir;      // path to the directory containing the profiles
  int heat_map;             // flag to generate heat map images
  long double res;          // resolution used in heat map images
  int read_profiles;        // flag to read profile(s)
  int mixture_model;        // flag to model a mixture
  int fit_num_components;   // # of components in the mixture model
  int infer_num_components; // flag to infer # of components
  int max_components;       // max components to infer
  string infer_log;         // log file
  int continue_inference;   // flag to continue inference from some state
  int simulation;           // flag to run mixture model simulation
  int load_mixture;         // flag to read mixture from a file
  int simulated_components; // # of components to be simulated
  string mixture_file;      // file containing the mixture information
  int sample_size;          // sample size to be generated from the simulated mixture
  int num_threads;          // flag to enable multithreading
  long double max_kappa;    // max value of kappa allowed
  int start_from;           // starting value of number of components
                            // during inference
  int estimate_all;         // estimate using all methods
  int compute_responsibility_matrix;  // flag
};

struct Estimates
{
  long double psi,alpha,eta;
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
Vector computeVectorSum(std::vector<Vector > &, Vector &, long double &);
Vector computeNormalizedVectorSum(std::vector<Vector > &);
Matrix computeDispersionMatrix(std::vector<Vector > &);
Matrix computeDispersionMatrix(std::vector<Vector > &, Vector &);
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
double Constraint2(const std::vector<double> &, std::vector<double> &, void *);
double Constraint5(const std::vector<double> &, std::vector<double> &, void *);

long double computeConstantTerm(int);
std::vector<std::vector<int> > updateBins(std::vector<Vector > &, long double);
void outputBins(std::vector<std::vector<int> > &, long double);
void computeEstimators(struct Parameters &);
bool gatherData(struct Parameters &, std::vector<Vector > &);
void modelOneComponent(struct Parameters &, std::vector<Vector > &);
void modelMixture(struct Parameters &, std::vector<Vector > &);
void simulateMixtureModel(struct Parameters &);
Vector generateFromSimplex(int);
std::vector<Kent> generateRandomComponents(int);
Vector generateRandomKappas(int);
Vector generateRandomBetas(Vector &);
Mixture inferComponents(Mixture &, int, ostream &);
void updateInference(Mixture &, Mixture &, ostream &, int);

void TestFunctions(void);
void RunExperiments(int);

Vector sort(Vector &);
void quicksort(Vector &, std::vector<int> &, int, int);
int partition(Vector &, std::vector<int> &, int, int);
long double computeMedian(Vector &);
long double computeMean(Vector &);
long double computeVariance(Vector &);
int maximumIndex(Vector &);

#endif

