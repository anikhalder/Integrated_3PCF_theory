#ifndef INTEGRATED_BISPECTRA_2D_H
#define INTEGRATED_BISPECTRA_2D_H

#include <ClassEngine.hh>
#include <cosmology_utils.h>
#include <interpolation_methods.h>
#include <gsl/gsl_rng.h>
#include <constants.h>

// ######################################################################################

// Integrated 3-point function area pre-factors

// Gaussian-quadrature

struct params_Adelta_theta_integrand {double theta_T;};

double Adelta_theta_qag_integrand(double theta, void *params);

double Adelta_qag(double theta_T); // this is identical to cosmology_utils::calculate_patch_area()

struct params_A2pt_phi_integrand {double theta; double phi_alpha; double alpha; double theta_T;};

double A2pt_phi_qag_integrand(double phi, void *params);

struct params_A2pt_theta_phi_integrand {double phi_alpha; double alpha; double theta_T;};

double A2pt_theta_qag_integrand(double theta, void *params);

double A2pt_qag(double alpha, double theta_T);

// Monte-Carlo

double A2pt_theta_phi_mc_integrand (double *k, size_t dim, void *params);

double A2pt_mc(double alpha, double theta_T, const gsl_rng_type *T, const std::string &mc_integration_type);

// h-cubature

int A2pt_theta_phi_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

double A2pt_hcubature(double alpha, double theta_T);

// ##############################

// testing with phi_alpha integration; this is redundant and time consuming ---> does not affect the results obtain only from w_alpha_area_theta_phi_integration()

struct params_A2pt_phi_alpha_theta_phi_integrand {double alpha; double theta_T;};

double A2pt_phi_alpha_qag_integrand(double phi_alpha, void *params);

double A2pt_angle_averaged_qag(double alpha, double theta_T);

// ##############################

// testing with bin averaging over A2pt(alpha) integration

struct params_A2pt_alpha_theta_phi_integrand {double phi_alpha; double theta_T;};

double A2pt_alpha_qag_integrand(double alpha, void *params);

double A2pt_bin_averaged_qag(double alpha_min, double alpha_max, double theta_T);

// ######################################################################################
// ######################################################################################
// ######################################################################################

struct struct_iB2D_W_FS { double (*W_meandelta_FS)(const double &l, const double &theta_T); const double &theta_T_meandelta;
                        double (*W_2pt_FS)(const double &l, const double &theta_T); const double &theta_T_2pt;};

double W_products(const std::string &key, const double &l_1, const double &phi_1, const double &l_2, const double &phi_2, const double &l, const double &phi_l,
                  const struct_iB2D_W_FS &info_iB2D_W_FS);

struct params_iB2D_phi_1_phi_2_integrand { const std::string &key; const double &l; const double &phi_l; const struct_iB2D_W_FS &info_iB2D_W_FS; const double &l_1;  const double &l_2;};

double evaluate_iB2D_phi_1_phi_2_integrand(const std::string &key, const double &l, const double &phi_l, const struct_iB2D_W_FS &info_iB2D_W_FS,
                                         const double &l_1, const double &l_2, const double &phi_1, const double &phi_2);

// h-cubature

int iB2D_phi_1_phi_2_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

void iB2D_phi_1_phi_2_hcubature(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, const double &l_1, const double &l_2,
                              double &result, double &error, size_t max_evals = calls_1e5);

// ######################################################################################

struct params_iB2D_l_1_l_2_phi_1_phi_2_integrand { const std::string &key; const double &l; const double &phi_l; const double &z; 
                                                   const struct_iB2D_W_FS &info_iB2D_W_FS; ClassEngine *class_obj; const bool &use_pk_nl;};

double evaluate_iB2D_l_1_l_2_phi_1_phi_2_integrand(const std::string &key, const double &l, const double &phi_l, const double &z, 
                                                   const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, bool use_pk_nl, 
                                                   const double &l_1, const double &l_2, const double &phi_1, const double &phi_2);

// Monte-Carlo

double iB2D_l_1_l_2_phi_1_phi_2_mc_integrand(double *k, size_t dim, void *params);

void iB2D_l_1_l_2_phi_1_phi_2_mc(const std::string &key, const double &l, const double &z, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                   std::vector<double> lower_limits, std::vector<double> upper_limits,
                                   const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls);

void iB2D_mc_4_dim(const std::string &key, const double &l, const double &z, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                    std::vector<double> lower_limits, std::vector<double> upper_limits,
                    const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls);

double iB2D_trapz_z(Linear_interp_1D *iB2D_z_interp, double &z_lower_limit, double &z_upper_limit, ClassEngine *class_obj, 
                    projection_kernel *q1, projection_kernel *q2, projection_kernel *q3);

// ######################################################################################

struct params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand { const std::string &key; const double &l; const double &phi_l; const struct_iB2D_W_FS &info_iB2D_W_FS; ClassEngine *class_obj;
                                                     const bool &use_pk_nl; projection_kernel *q1; projection_kernel *q2; projection_kernel *q3;};

struct params_iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_integrand { const std::string &key; const double &l; const struct_iB2D_W_FS &info_iB2D_W_FS; ClassEngine *class_obj;
                                                           const bool &use_pk_nl; projection_kernel *q1; projection_kernel *q2; projection_kernel *q3;};                                                     

double evaluate_iB2D_z_l_1_l_2_phi_1_phi_2_integrand(const std::string &key, const double &l, const double &phi_l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj,
                                                     bool use_pk_nl, const double &z, const double &l_1, const double &l_2, const double &phi_1, const double &phi_2,
                                                     projection_kernel *q1, projection_kernel *q2, projection_kernel *q3);

// Monte-Carlo

double iB2D_z_l_1_l_2_phi_1_phi_2_mc_integrand(double *k, size_t dim, void *params);

void iB2D_z_l_1_l_2_phi_1_phi_2_mc(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                   projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                   const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls);

void iB2D_mc(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
           projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
           const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls);

// Monte-Carlo angle averaged

double iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc_integrand(double *k, size_t dim, void *params);

void iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                        projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                        const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls);

void iB2D_mc_angle_averaged(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                          projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                          const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls);

// h-cubature

int iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

void iB2D_z_l_1_l_2_phi_1_phi_2_hcubature(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                          projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                          double &result, double &error, size_t max_evals);

void iB2D_hcubature(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                  projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                  double &result, double &error, size_t max_evals);

// h-cubature with 4 dimensions

int iB2D_z_l_1_l_2_phi_1_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

void iB2D_z_l_1_l_2_phi_1_hcubature(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                          projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                          double &result, double &error, size_t max_evals);

void iB2D_hcubature_4dim(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                       projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                       double &result, double& error, size_t max_evals);

// h-cubature vectorized

int iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_v_integrand(unsigned ndim, unsigned npts, const double *k, void *params, unsigned fdim, double *value);

void iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_v(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                            projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                            double &result, double &error, size_t max_evals);

double iB2D_hcubature_v(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                      projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                      double &result, double &error, size_t max_evals);


// h-cubature angle averaged

int iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

void iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_hcubature(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                                projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                                double &result, double &error, size_t max_evals);

void iB2D_hcubature_angle_averaged(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                 projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                 double &result, double &error, size_t max_evals);

#endif // INTEGRATED_BISPECTRA_2D_H
