#ifndef TEST_H
#define TEST_H

#include "Header.h"

class Test
{
  private:

  public:
    void load_data();

    void bessel();

    void testing_cartesian2sphericalPoleXAxis();

    void parallel_sum_computation(void);

    void uniform_number_generation();

    void arbitrary_rotation();

    void matrixFunctions(void);

    void productMatrixVector(void);

    void dispersionMatrix(void);

    void numericalIntegration(void);

    void normalDistributionFunctions(void);

    void orthogonalTransformations(void);

    void orthogonalTransformations2(void);

    void randomSampleGeneration(void);

    void multivariate_normal(void);

    void acg(void);

    void bingham(void);

    void kent_bingham_generation(void);

    void normalization_constant(void);

    void optimization(void);

    void moment_estimation(void);

    void ml_estimation(void);

    void expectation();

    void kl_divergence();

    void fisher();

    void fisher2();

    void mml_estimation(void);

    void mml_estimation2(void);

    void plot_posterior_density();

    void vmf_all_estimation();

    void chi_square();

    void hypothesis_testing();

    void hypothesis_testing2();

    void confidence_region();

    void infer_mixture();

    void infer_mixture_vmf();

    void contours();

    void gsl_numerical_integration();

    void gsl_numerical_kent_density_integration();

    void gsl_monte_carlo_integration();

    void gsl_monte_carlo_kent_density_integration();

    void testing_sample_empirical_distribution();
};

#endif

