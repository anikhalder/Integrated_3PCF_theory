#ifndef POWER_SPECTRA_2D_H
#define POWER_SPECTRA_2D_H

#include <ClassEngine.hh>
#include <cosmology_utils.h>
#include <gsl/gsl_rng.h>

// ######################################################################################

// 2D power spectrum (integrated along the line-of-sight)

double evaluate_P2D_z_integrand(const std::string &key, const double &l, ClassEngine *class_obj, bool use_pk_nl, const double &z,
                                projection_kernel *q1, projection_kernel *q2);

struct params_P2D_z_integrand { const std::string &key; double l; ClassEngine *class_obj; bool use_pk_nl; projection_kernel *q1; projection_kernel *q2;};

// Gaussian-quadrature

double P2D_z_qag_integrand(double z, void *params);

double P2D_z_qag(const std::string &key, const double &l, ClassEngine *class_obj, const bool &use_pk_nl,
                 projection_kernel *q1, projection_kernel *q2, const double &z_lower, const double &z_upper);

// Monte-Carlo

double P2D_z_mc_integrand(double *k, size_t dim, void *params);

double P2D_z_mc(const std::string &key, const double &l, ClassEngine *class_obj, const bool &use_pk_nl,
                projection_kernel *q1, projection_kernel *q2, std::vector<double> &z_lower, std::vector<double> &z_upper,
                const gsl_rng_type *T, const std::string &mc_integration_type);

// h-cubature

int P2D_z_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

double P2D_z_hcubature(const std::string &key, const double &l, ClassEngine *class_obj, const bool &use_pk_nl,
                       projection_kernel *q1, projection_kernel *q2, std::vector<double> &z_lower, std::vector<double> &z_upper);

#endif // POWER_SPECTRA_2D_H
