#ifndef HALO_UTILS_HPP
#define HALO_UTILS_HPP

#include <ClassEngine.hh>

double R_L(const double &M, ClassEngine *class_obj); // comoving Lagrangian radius [Mpc] of a halo of mass M [M_sun]

double halo_b1_PBS(const double &M, const double &z, ClassEngine *class_obj);

double halo_b2_PBS(const double &M, const double &z, ClassEngine *class_obj);

double halo_b2_Lazeyras(const double &M, const double &z, ClassEngine *class_obj);

double halo_bs2_coevolution(const double &M, const double &z, ClassEngine *class_obj);

double halo_b1_Tinker2010(const double &M, const double &z, ClassEngine *class_obj);

double f_nu_PS(const double &nu);

double dn_dM_PS(const double &M, const double &z, ClassEngine *class_obj);

double g_sigma_Tinker2010(const double &sigma_R_z, const double &z);

double f_sigma_Tinker2008(const double &sigma_R_z, const double &z);

double dn_dM_Tinker2008(const double &M, const double &z, ClassEngine *class_obj);

// h-cubature integration routine to compute the total number of halos in a give redshift, mass bin

struct params_N_h_z_M_integrand { ClassEngine *class_obj;};

int N_h_z_M_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value);

double N_h_hcubature(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max);

#endif // HALO_UTILS_HPP
