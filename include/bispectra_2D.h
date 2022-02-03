#ifndef BISPECTRA_2D_H
#define BISPECTRA_2D_H

#include <ClassEngine.hh>
#include <cosmology_utils.h>
#include <gsl/gsl_rng.h>

// ######################################################################################

// 2D bispectrum (integrated along the line-of-sight)

double evaluate_B2D_z_integrand(const std::string &key, const double &l1, const double &l2, const double &l3, ClassEngine *class_obj, bool use_pk_nl, const double &z,
                                const double &q_1, const double &q_2, const double &q_3);

struct params_B2D_z_integrand { const std::string &key; double l1; double l2; double l3; ClassEngine *class_obj; bool use_pk_nl; projection_kernel *q1; projection_kernel *q2; projection_kernel *q3;};

// Gaussian-quadrature

double B2D_z_qag_integrand(double z, void *params);

double B2D_z_qag(const std::string &key, const double &l1, const double &l2, const double &l3, ClassEngine *class_obj, const bool &use_pk_nl,
                 projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, const double &z_lower, const double &z_upper);

// Monte-Carlo

double B2D_z_mc_integrand(double *k, size_t dim, void *params);

double B2D_z_mc(const std::string &key, const double &l1, const double &l2, const double &l3, ClassEngine *class_obj, const bool &use_pk_nl,
                projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> &z_lower, std::vector<double> &z_upper,
                const gsl_rng_type *T, const std::string &mc_integration_type);

// h-cubature

int B2D_z_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

double B2D_z_hcubature(const std::string &key, const double &l1, const double &l2, const double &l3, ClassEngine *class_obj, const bool &use_pk_nl,
                       projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> &z_lower, std::vector<double> &z_upper);

#endif // BISPECTRA_2D_H
