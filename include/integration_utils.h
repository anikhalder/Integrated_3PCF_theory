#ifndef INTEGRATION_UTILS_H
#define INTEGRATION_UTILS_H

#include <ClassEngine.hh>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <cubature.h>
//#include <cuba.h>

// ######################################################################################

// integration utility functions

void qag_1D_integration(double (*func)(double, void*), void *args, const double &lower_limit, const double &upper_limit, const size_t &workspace_size,
                        double &result, double &error);

void qag_1D_integration_abs_rel(double (*func)(double, void*), void *args, const double &lower_limit, const double &upper_limit, const size_t &workspace_size,
                                double &result, double &error);

void qagiu_1D_integration(double (*func)(double, void*), void *args, const double &lower_limit, const size_t &workspace_size, double &result, double &error);

void qagil_1D_integration(double (*func)(double, void*), void *args, const double &upper_limit, const size_t &workspace_size, double &result, double &error);

void monte_carlo_plain_integration(gsl_monte_function *G, std::vector<double> &lower_limits, std::vector<double> &upper_limits,
                                   const size_t &dim, const size_t &calls, const gsl_rng_type *T, double &result, double &error);

void monte_carlo_miser_integration(gsl_monte_function *G, std::vector<double> &lower_limits, std::vector<double> &upper_limits,
                                   const size_t &dim, const size_t &calls, const gsl_rng_type *T, double &result, double &error);

void monte_carlo_vegas_integration(gsl_monte_function *G, std::vector<double> &lower_limits, std::vector<double> &upper_limits,
                                   const size_t &dim, const size_t &calls, const gsl_rng_type *T, double &result, double &error);

void hcubature_integration(integrand func, void *args, std::vector<double> &lower_limits, std::vector<double> &upper_limits, const size_t &dim,
                           const size_t &max_evals, double &result, double &error);

void hcubature_v_integration(integrand_v func, void *args, std::vector<double> &lower_limits, std::vector<double> &upper_limits, const size_t &dim,
                             const size_t &max_evals, double &result, double &error);

//void cuba_cuhre_integration(integrand_t func, void *args, std::vector<double> &lower_limits, std::vector<double> &upper_limits, const size_t &dim,
//                            const size_t &max_evals, double &result, double &error);

#endif // INTEGRATION_UTILS_H
