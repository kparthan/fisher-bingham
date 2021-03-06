#include <gsl/gsl_multimin.h>
#include <iostream>

using namespace std;

/* Paraboloid centered on (p[0],p[1]), with  
   scale factors (p[2],p[3]) and minimum p[4] */

double
my_f (const gsl_vector *v, void *params)
{
  double x, y;
  double *p = (double *)params;
  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
 
  return p[2] * (x - p[0]) * (x - p[0]) +
           p[3] * (y - p[1]) * (y - p[1]) + p[4]; 
}

/* The gradient of f, df = (df/dx, df/dy). */
void 
my_df (const gsl_vector *v, void *params, 
       gsl_vector *df)
{
  double x, y;
  double *p = (double *)params;
  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
 
  gsl_vector_set(df, 0, 2.0 * p[2] * (x - p[0]));
  gsl_vector_set(df, 1, 2.0 * p[3] * (y - p[1]));
}

/* Compute both f and df together. */
void 
my_fdf (const gsl_vector *x, void *params, 
        double *f, gsl_vector *df) 
{
  *f = my_f(x, params); 
  my_df(x, params, df);
}

int 
main(void)
{
  double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};

  const gsl_multimin_fminimizer_type *T = 
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  /* Starting point */
  x = gsl_vector_alloc (2);
  gsl_vector_set (x, 0, 5.0);
  gsl_vector_set (x, 1, 7.0);

  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, 1.0);

  /* Initialize method and iterate */
  minex_func.n = 2;
  minex_func.f = my_f;
  minex_func.params = par;

  s = gsl_multimin_fminimizer_alloc (T, 2);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status) 
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-2);

      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }

      printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
              iter,
              gsl_vector_get (s->x, 0), 
              gsl_vector_get (s->x, 1), 
              s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < 100);
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return status;
}

/*double momentObjectiveFunction(const gsl_vector *x, void *params)
{
  double k = gsl_vector_get(x,0);
  double b = gsl_vector_get(x,1);

  double *p = (double *)params;

  Kent kent(k,b);
  long double log_norm = kent.computeLogNormalizationConstant();
  double fval = log_norm - k * C1 - b * C2;

  cout << "x: (" << k << "," << b << "); fx: " << fval << endl;
  return fval;
}

void Optimize::minimize(void)
{
  double constants[2] = {C1,C2};

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  // Starting point 
  x = gsl_vector_alloc (2);
  gsl_vector_set (x, 0, 25.0);
  gsl_vector_set (x, 1, 15.0);

  // Set initial step sizes to 1 
  ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, 1.0);

  // Initialize method and iterate 
  minex_func.n = 2;
  minex_func.f = momentObjectiveFunction;
  minex_func.params = constants;

  s = gsl_multimin_fminimizer_alloc (T, 2);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    
    if (status) 
      break;

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, 1e-2);

    if (status == GSL_SUCCESS) {
      printf ("converged to minimum at\n");
    }

    printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
            iter,
            gsl_vector_get (s->x, 0), 
            gsl_vector_get (s->x, 1), 
            s->fval, size);
  }
  while (status == GSL_CONTINUE && iter < 100);
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  //return status;
}
*/
