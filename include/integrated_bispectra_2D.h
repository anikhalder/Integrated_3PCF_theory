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

// ######################################################################################
// ######################################################################################
// ######################################################################################

struct struct_iB_W_FS { double (*W_meandelta_FS)(const double &l, const double &theta_T); const double &theta_T_meandelta;
                        double (*W_2pt_FS)(const double &l, const double &theta_T); const double &theta_T_2pt;};

double W_products(const std::string &key, const double &l_1, const double &phi_1, const double &l_2, const double &phi_2, const double &l, const double &phi_l,
                  const struct_iB_W_FS &info_iB_W_FS);

struct params_iB_phi_1_phi_2_integrand { const std::string &key; const double &l; const double &phi_l; const struct_iB_W_FS &info_iB_W_FS; const double &l_1;  const double &l_2;};

double evaluate_iB_phi_1_phi_2_integrand(const std::string &key, const double &l, const double &phi_l, const struct_iB_W_FS &info_iB_W_FS,
                                         const double &l_1, const double &l_2, const double &phi_1, const double &phi_2);

// h-cubature

int iB_phi_1_phi_2_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

void iB_phi_1_phi_2_hcubature(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, const double &l_1, const double &l_2,
                              double &result, double &error, size_t max_evals = calls_1e5);

// ######################################################################################

struct params_iB_los_l_1_l_2_phi_1_phi_2_integrand { const std::string &key; const double &l; const double &phi_l; const struct_iB_W_FS &info_iB_W_FS; ClassEngine *class_obj;
                                                     const bool &use_pk_nl; projection_kernel *q1; projection_kernel *q2; projection_kernel *q3;};

double evaluate_iB_los_l_1_l_2_phi_1_phi_2_integrand(const std::string &key, const double &l, const double &phi_l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj,
                                                     bool use_pk_nl, const double &z, const double &l_1, const double &l_2, const double &phi_1, const double &phi_2,
                                                     projection_kernel *q1, projection_kernel *q2, projection_kernel *q3);

struct params_iB_phi_l_los_l_1_l_2_phi_1_phi_2_integrand { const std::string &key; const double &l; const struct_iB_W_FS &info_iB_W_FS; ClassEngine *class_obj;
                                                           const bool &use_pk_nl; projection_kernel *q1; projection_kernel *q2; projection_kernel *q3;};

// Monte-Carlo

double iB_los_l_1_l_2_phi_1_phi_2_mc_integrand(double *k, size_t dim, void *params);

void iB_los_l_1_l_2_phi_1_phi_2_mc(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                   projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                   const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls);

void iB_mc(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
           projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
           const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls);

// Monte-Carlo angle averaged

double iB_phi_l_los_l_1_l_2_phi_1_phi_2_mc_integrand(double *k, size_t dim, void *params);

void iB_phi_l_los_l_1_l_2_phi_1_phi_2_mc(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                        projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                        const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls);

void iB_mc_angle_averaged(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                          projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                          const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls);

// Monte-Carlo CIGAR

struct params_iB_los_l_1_l_2_phi_1_phi_2_mc_cigar_integrand { const std::string &key; const double &l; const double &phi_l; const struct_iB_W_FS &info_iB_W_FS; ClassEngine *class_obj;
                                                              const bool &use_pk_nl; projection_kernel *q1; projection_kernel *q2; projection_kernel *q3; 
                                                              std::vector<double> &lower_limits; std::vector<double> &upper_limits;};

double iB_phi_l_los_l_1_l_2_phi_1_phi_2_mc_cigar_integrand(std::vector<double> k, void *params);

void iB_los_l_1_l_2_phi_1_phi_2_mc_cigar(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                         projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                         double &result, double &error);

void iB_mc_cigar(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                 projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                 double &result, double &error);

// h-cubature

int iB_los_l_1_l_2_phi_1_phi_2_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

void iB_los_l_1_l_2_phi_1_phi_2_hcubature(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                          projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                          double &result, double &error, size_t max_evals);

void iB_hcubature(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                  projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                  double &result, double &error, size_t max_evals);

// h-cubature with 4 dimensions

int iB_los_l_1_l_2_phi_1_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

void iB_los_l_1_l_2_phi_1_hcubature(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                          projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                          double &result, double &error, size_t max_evals);

void iB_hcubature_4dim(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                       projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                       double &result, double& error, size_t max_evals);

// h-cubature vectorized

int iB_los_l_1_l_2_phi_1_phi_2_hcubature_v_integrand(unsigned ndim, unsigned npts, const double *k, void *params, unsigned fdim, double *value);

void iB_los_l_1_l_2_phi_1_phi_2_hcubature_v(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                            projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                            double &result, double &error, size_t max_evals);

double iB_hcubature_v(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                      projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                      double &result, double &error, size_t max_evals);


// h-cubature angle averaged

int iB_phi_l_los_l_1_l_2_phi_1_phi_2_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

void iB_phi_l_los_l_1_l_2_phi_1_phi_2_hcubature(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                                projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                                double &result, double &error, size_t max_evals);

void iB_hcubature_angle_averaged(const std::string &key, const double &l, const struct_iB_W_FS &info_iB_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                 projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                 double &result, double &error, size_t max_evals);

#endif // INTEGRATED_BISPECTRA_2D_H
