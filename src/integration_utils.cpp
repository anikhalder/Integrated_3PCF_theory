#include <integration_utils.h>
#include <gsl/gsl_integration.h>
#include <constants.h>

// ######################################################################################

// integration utility functions

void qag_1D_integration(double (*func)(double, void*), void *args, const double &lower_limit, const double &upper_limit, const size_t &workspace_size,
                        double &result, double &error)
{
    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(workspace_size);

    gsl_function integrand;
    integrand.function = func;
    integrand.params = args;

    gsl_integration_qag(&integrand, lower_limit, upper_limit, 0, error_m06, workspace_size, 6, work_ptr, &result, &error);

    gsl_integration_workspace_free(work_ptr);
}

void qag_1D_integration_abs_rel(double (*func)(double, void*), void *args, const double &lower_limit, const double &upper_limit, const size_t &workspace_size,
                                double &result, double &error)
{
    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(workspace_size);

    gsl_function integrand;
    integrand.function = func;
    integrand.params = args;

    //gsl_integration_qag(&integrand, lower_limit, upper_limit, error_m04, error_m04, workspace_size, 6, work_ptr, &result, &error);
    gsl_integration_qag(&integrand, lower_limit, upper_limit, error_m03, error_m03, workspace_size, 6, work_ptr, &result, &error);
    //gsl_integration_qag(&integrand, lower_limit, upper_limit, error_m04, 0, workspace_size, 6, work_ptr, &result, &error);

    gsl_integration_workspace_free(work_ptr);
}

void qagiu_1D_integration(double (*func)(double, void *), void *args, const double &lower_limit, const size_t &workspace_size, double &result, double &error)
{
    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(workspace_size);

    gsl_function integrand;
    integrand.function = func;
    integrand.params = args;

    gsl_integration_qagiu(&integrand, lower_limit, 0, error_m06, workspace_size, work_ptr, &result, &error);

    gsl_integration_workspace_free(work_ptr);
}

void qagil_1D_integration(double (*func)(double, void *), void *args, const double &upper_limit, const size_t &workspace_size, double &result, double &error)
{
    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(workspace_size);

    gsl_function integrand;
    integrand.function = func;
    integrand.params = args;

    gsl_integration_qagil(&integrand, upper_limit, 0, error_m06, workspace_size, work_ptr, &result, &error);

    gsl_integration_workspace_free(work_ptr);
}

void monte_carlo_plain_integration(gsl_monte_function *G, std::vector<double> &lower_limits, std::vector<double> &upper_limits,
                                   const size_t &dim, const size_t &calls, const gsl_rng_type *T, double &result, double &error)
{
    gsl_rng *r = gsl_rng_alloc(T);

    gsl_monte_plain_state *s = gsl_monte_plain_alloc (dim);
    gsl_monte_plain_integrate (G, lower_limits.data(), upper_limits.data(), dim, calls, r, s, &result, &error);
    gsl_monte_plain_free (s);

    gsl_rng_free (r);
}

void monte_carlo_miser_integration(gsl_monte_function *G, std::vector<double> &lower_limits, std::vector<double> &upper_limits,
                                   const size_t &dim, const size_t &calls, const gsl_rng_type *T, double &result, double &error)
{
    gsl_rng *r = gsl_rng_alloc(T);

    gsl_monte_miser_state *s = gsl_monte_miser_alloc (dim);
    gsl_monte_miser_integrate (G, lower_limits.data(), upper_limits.data(), dim, calls, r, s, &result, &error);
    gsl_monte_miser_free (s);

    gsl_rng_free (r);
}

void monte_carlo_vegas_integration(gsl_monte_function *G, std::vector<double> &lower_limits, std::vector<double> &upper_limits,
                                   const size_t &dim, const size_t &calls, const gsl_rng_type *T, double &result, double &error)
{
    gsl_rng *r = gsl_rng_alloc(T);
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);

    size_t calls_warmup;
    if (calls >= calls_1e4)
        calls_warmup = calls_1e4;
    else
        calls_warmup = calls_1e3;

    gsl_monte_vegas_integrate (G, lower_limits.data(), upper_limits.data(), dim, calls_warmup, r, s, &result, &error); // to warm up the grid

    gsl_monte_vegas_integrate (G, lower_limits.data(), upper_limits.data(), dim, calls, r, s, &result, &error);
//    do
//    {
//        gsl_monte_vegas_integrate (G, lower_limits.data(), upper_limits.data(), dim, calls, r, s, &result, &error);
//    }
//    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s);
    gsl_rng_free (r);
}

void hcubature_integration(integrand func, void *args, std::vector<double> &lower_limits, std::vector<double> &upper_limits, const size_t &dim,
                           const size_t &max_evals, double &result, double &error)
{
    // The implementation here is only meant for a NON-vector valued integrand i.e. fdim = 1 --> see documentation at https://github.com/stevengj/cubature

    unsigned ndim = dim;
    unsigned fdim = 1;

    //hcubature(fdim, func, args, ndim, lower_limits.data(), upper_limits.data(), max_evals, 0, error_m08, ERROR_INDIVIDUAL, &result, &error);
    hcubature(fdim, func, args, ndim, lower_limits.data(), upper_limits.data(), max_evals, 0, error_m04, ERROR_INDIVIDUAL, &result, &error);
}

void hcubature_v_integration(integrand_v func, void *args, std::vector<double> &lower_limits, std::vector<double> &upper_limits, const size_t &dim,
                             const size_t &max_evals, double &result, double &error)
{
    // The implementation here is only meant for a NON-vector valued integrand i.e. fdim = 1 --> see documentation at https://github.com/stevengj/cubature

    unsigned ndim = dim;
    unsigned fdim = 1;

    hcubature_v(fdim, func, args, ndim, lower_limits.data(), upper_limits.data(), max_evals, 0, error_m08, ERROR_INDIVIDUAL, &result, &error);
}

//void cuba_cuhre_integration(integrand_t func, void *args, std::vector<double> &lower_limits, std::vector<double> &upper_limits, const size_t &dim,
//                            const size_t &max_evals, double &result, double &error)
//{
//    // The implementation here is only meant for a NON-vector valued integrand

//    const int ndim = dim;
//    const int ncomp = 1; // same as fdim in hcubature

//    double prob;
//    int nregions, neval, fail;

//    Cuhre(ndim, ncomp, func, args, 1, error_m03, 0, 0, 0, max_evals, 9, nullptr, nullptr,  &nregions, &neval, &fail, &result, &error, &prob);
//}
