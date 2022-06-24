#ifndef HALO_UTILS_HPP
#define HALO_UTILS_HPP

#include <ClassEngine.hh>

double R_L(const double &M, ClassEngine *class_obj); // comoving Lagrangian radius [Mpc] of a halo of mass M [M_sun]

double halo_b1_PBS(const double &M, const double &z, ClassEngine *class_obj);

double halo_b2_PBS(const double &M, const double &z, ClassEngine *class_obj);

double halo_b2_Lazeyras(const double &M, const double &z, ClassEngine *class_obj);

double halo_bs2_coevolution(const double &M, const double &z, ClassEngine *class_obj);

double halo_b1_nu_Tinker2010(const double &nu);

double halo_b1_Tinker2010(const double &M, const double &z, ClassEngine *class_obj);

double f_nu_PS(const double &nu);

double dn_dM_PS(const double &M, const double &z, ClassEngine *class_obj);

double f_nu_Tinker2010(const double &nu, const double &z);

double g_sigma_Tinker2010(const double &sigma_R_z, const double &z);

double f_sigma_Tinker2008(const double &sigma_R_z, const double &z);

double dn_dM_Tinker2008(const double &M, const double &z, ClassEngine *class_obj);

// h-cubature integration routine to compute the total number of halos in a given redshift range, mass bin

struct params_N_h_z_M_integrand { ClassEngine *class_obj;};

int N_h_z_M_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

double N_h_hcubature(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max);

// HOD functions

double N_gal_expected_Zacharegkas2020(const double &M, const double &z);

struct params_HOD_M_integrand { const double &z; ClassEngine *class_obj;};

double n_z_gal_M_qag_integrand(double &M, void *params);

double n_z_gal_qag(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max); // galaxy number density from HOD

struct params_HOD_z_integrand { ClassEngine *class_obj; const double &M_min; const double &M_max;};

double N_gal_z_qag_integrand(double &z, void *params);

double N_gal_qag(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max); // Number of galaxies in given redshift range

double b1_avg_z_gal_M_qag_integrand(double M, void *params);

double b1_avg_z_gal_qag(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max); // average galaxy bias from HOD

struct params_b1_avg_gal_z_M_integrand {ClassEngine *class_obj; Linear_interp_1D *n_of_z; const double &M_min; const double &M_max;};

int b1_avg_gal_z_M_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

double b1_avg_gal(ClassEngine *class_obj, Linear_interp_1D *n_of_z, const double &z_min, const double &z_max, const double &M_min, const double &M_max);

// bias and halo mass function integral constraints

struct params_integral_constraints_integrand { const double &z; };

double b1_integral_constraint_nu_qag_integrand(double nu, void *params);

double b1_integral_constraint_qag(const double &z);

double g_sigma_integral_constraint_sigma_qag_integrand(double sigma, void *params);

double g_sigma_integral_constraint_qag(const double &z);

double f_nu_integral_constraint_nu_qag_integrand(double nu, void *params);

double f_nu_integral_constraint_qag(const double &z);
 
#endif // HALO_UTILS_HPP
