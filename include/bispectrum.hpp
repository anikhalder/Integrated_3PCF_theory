#ifndef BISPECTRUM_HPP
#define BISPECTRUM_HPP

#include <ClassEngine.hh>
#include <cosmology_utils.h>

bool is_triangle_closed(const double &k_1, const double &k_2, const double &k_3);

// ######################################################################################

// bispectrum prescription

double B(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl);

// ######################################################################################

// tree-level SPT matter bispectrum

double B_tree(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl);

// ######################################################################################

// components for tracer bispectrum

double B_P(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl);

double B_PP(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl);

double B_S2PP(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl);

double B_PaPb(const double &k_a, const double &k_b, const double &z, ClassEngine *class_obj, bool use_pk_nl);

double B_S2PaPb(const double &k_a, const double &k_b, const double &k_c, const double &z, ClassEngine *class_obj, bool use_pk_nl);

// ######################################################################################

// 1-loop SPT matter bispectrum

struct params_B_1_loop_integrand { const std::vector<double> &k_1_vec; const std::vector<double> &k_2_vec; const std::vector<double> &k_3_vec;
                                   const double &z; ClassEngine *class_obj; const bool &use_pk_nl;};

double evaluate_B_1_loop_integrand(const std::vector<double> &q_1_vec, const std::vector<double> &k_a_vec, const std::vector<double> &k_b_vec,
                                   const std::vector<double> &k_c_vec, const double &z, ClassEngine *class_obj, bool use_pk_nl);

// h-cubature

int B_1_loop_hcubature_integrand(unsigned ndim, const double *q, void *params, unsigned fdim, double *value);

double B_1_loop_hcubature(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl);

// ######################################################################################

// matter bispectrum fitting formulae

double n(const double &k, const double &z, ClassEngine *class_obj);

double Q3(const double &n);

double F2_eff(const double &k_1, const double &k_2, const double &k_3,
              const double &a_1, const double &a_2, const double &b_1, const double &b_2, const double &c_1, const double &c_2);

// ######################################################################################

// Scoccimarro and Couchman matter bispectrum fitting formula

double a_SC(const double &n, const double &q, const double &sigma8_z);

double b_SC(const double &n, const double &q);

double c_SC(const double &n, const double &q);

double B_SC(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl=true);

// ######################################################################################

// Gil-Marin matter bispectrum fitting formula

double a_GM(const double &n, const double &q, const double &sigma8_z);

double b_GM(const double &n, const double &q);

double c_GM(const double &n, const double &q);

double B_GM(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl=true);

double B_GM_v2(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl=true);

// ######################################################################################

// bihalofit fitting formula

double B_bihalofit(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl);

double B_bihalofit_baryon_ratio(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj); // bihalofit bispectrum ratio with to without baryons

// ######################################################################################

// response function squeezed matter bispectrum

double B_squeezed_RF(const double &k_h, const double &k_m, const double &k_s, const double &z, ClassEngine *class_obj, bool use_pk_nl=true);

// ######################################################################################
// ######################################################################################
// ######################################################################################

// primordial non-Gaussianity

// Local type

double B_primordial_local(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj);

double B_primordial_equilateral(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj);

double B_primordial_orthogonal(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj);

#endif // BISPECTRUM_HPP
