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
void writeToFile(vector<vector<long double> > &, const char *);
string extractName(string &);
void print(ostream &, vector<long double> &, int);
int sign(long double);
long double exponent(long double, long double);
long double normalize(vector<long double> &, vector<long double> &);
void cartesian2spherical(vector<long double> &, vector<long double> &);
void spherical2cartesian(vector<long double> &, vector<long double> &);
long double computeDotProduct(vector<long double> &, vector<long double> &);

void Test(void);

#endif

