#ifndef COSMOLOGY_UTILS_H
#define COSMOLOGY_UTILS_H

#include <ClassEngine.hh>
#include <interpolation_methods.h>

// ######################################################################################

// table utility functions

std::vector<double> read_1_column_table(const std::string &file_name);

std::vector<std::vector<double>> read_2_column_table(const std::string &file_name);

std::vector<std::vector<double>> read_3_column_table(const std::string &file_name);

std::vector<std::vector<double>> read_6_column_table(const std::string &file_name);

std::vector<std::vector<double>> read_7_column_table(const std::string &file_name);

std::vector<std::vector<double>> read_11_column_table(const std::string &file_name);

void normalise_nofz(std::vector<std::vector<double>> &matrix);

// ######################################################################################

size_t num_correlations(const size_t &n, const size_t &k);

// ######################################################################################

// 2D window functions

double spherical_cap_sqdeg_2_radius(const double &patch_area_sq_degrees);

double spherical_cap_radius_2_sqradians(const double &theta_T);

double W2D_TH_RS(const double &length, const double &theta_T);

double W2D_TH_RS_unnormalised(const double &length, const double &theta_T);

double W2D_TH_FS(const double &l, const double &theta_T);

double W2D_U_RS(const double &length, const double &theta_T);

double W2D_U_FS(const double &l, const double &theta_T);

// ######################################################################################

// 3D window functions

double W3D_TH_FS(const double &k, const double &R);

double W3D_prime_TH_FS(const double &k, const double &R);

// ######################################################################################

// wavevector arithmetic utility functions

double l_ApB(const double &l_A, const double &phi_A, const double &l_B, const double &phi_B);

double phi_ApB(const double &l_A, const double &phi_A, const double &l_B, const double &phi_B);

double l_ApBpC(const double &l_A, const double &phi_A, const double &l_B, const double &phi_B, const double &l_C, const double &phi_C);

double phi_ApBpC(const double &l_A, const double &phi_A, const double &l_B, const double &phi_B, const double &l_C, const double &phi_C);

// ######################################################################################

// cosmology utility functions

struct params_proper_time_integrand {ClassEngine *class_obj;};

double proper_time_integrand(double zp, void *params);

double proper_time(const double &z, ClassEngine *class_obj);

struct params_sigma_squared_R_z_integrand {double R; double z; ClassEngine *class_obj;};

double sigma_squared_R_z_integrand(double k, void *params);

double sigma_squared_R_z(const double &R, const double &z, ClassEngine *class_obj);

double sigma_squared_prime_R_z_integrand(double k, void *params);

double sigma_squared_prime_R_z(const double &R, const double &z, ClassEngine *class_obj);

double pk_lin_manual(const double &k, const double &z, ClassEngine *class_obj, Linear_interp_1D& Tk_z0_d_tot);

// ######################################################################################

// 2D projection kernels

double q_m(const double &z, ClassEngine *class_obj, Linear_interp_1D *n_m_of_z);

struct params_W_k_zs_distribution_integrand { double z; ClassEngine *class_obj; Linear_interp_1D *n_source_of_z;};

double W_k_zs_distribution_integrand(double zs, void *params);

double q_k_zs_distribution(const double &z, ClassEngine *class_obj, Linear_interp_1D *n_source_of_z, const double &z_max);

double q_k_zs_fixed(const double &z, ClassEngine *class_obj, const double &zs);

double cmb_lensing_convergence_kernel(const double &z, ClassEngine *class_obj);

struct params_q_h_integrand { double z; ClassEngine *class_obj;};

double q_h_b1_integrand(double M, void *params);

double q_h_b1_unnormalised(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max);

double q_h_b2_integrand(double M, void *params);

double q_h_b2_unnormalised(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max);

double q_h_bs2_integrand(double M, void *params);

double q_h_bs2_unnormalised(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max);

double q_h_without_b_integrand(double M, void *params);

double q_h_without_b_unnormalised(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max);

double n_h_z(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max);

class projection_kernel
{
public:
    virtual double evaluate(const double &z) = 0;

    virtual ~projection_kernel() = 0;
};

class projection_kernel_q_m : public projection_kernel
{
    ClassEngine *m_class_obj;
    Linear_interp_1D *m_n_m_of_z;

public:
    projection_kernel_q_m(ClassEngine *class_obj, Linear_interp_1D *n_m_of_z);

    double evaluate(const double &z);

    ~projection_kernel_q_m();
};

class projection_kernel_q_k_zs_fixed : public projection_kernel
{
    ClassEngine *m_class_obj;
    double m_zs;

    Linear_interp_1D m_q_k_zs_fixed_z_array;

public:
    projection_kernel_q_k_zs_fixed(ClassEngine *class_obj, const double &zs);

    double evaluate(const double &z);

    ~projection_kernel_q_k_zs_fixed();
};

class projection_kernel_q_k_zs_distribution : public projection_kernel
{
    ClassEngine *m_class_obj;
    Linear_interp_1D *m_n_source_of_z;
    double m_z_max;

    Linear_interp_1D m_q_k_zs_distribution_z_array;

public:
    projection_kernel_q_k_zs_distribution(ClassEngine *class_obj, Linear_interp_1D *n_source_of_z, const double &z_max);

    double evaluate(const double &z);

    ~projection_kernel_q_k_zs_distribution();
};

class projection_kernel_q_h : public projection_kernel
{
    ClassEngine *m_class_obj;
    double m_z_min;
    double m_z_max;
    double m_M_min;
    double m_M_max;
    double m_N_h;
    double m_n_h;

    Linear_interp_1D m_q_h_without_b_z_interp_array;
    Linear_interp_1D m_n_h_z_interp_array;

public:
    projection_kernel_q_h(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max);

    double evaluate(const double &z);
    double get_N_h(); // Number of halos
    double get_n_h(); // number density of halos [Mpc^-3]
    double get_n_h_z(const double &z); // number density of halos [Mpc^-3] at a given redshift

    ~projection_kernel_q_h();
};

class projection_kernel_q_h_b1 : public projection_kernel
{
    ClassEngine *m_class_obj;
    double m_z_min;
    double m_z_max;
    double m_M_min;
    double m_M_max;
    double m_N_h;
    double m_n_h;

    Linear_interp_1D m_q_h_b1_z_interp_array;

public:
    projection_kernel_q_h_b1(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max);

    double evaluate(const double &z);
    double get_N_h(); // Number of halos
    double get_n_h(); // number density of halos [Mpc^-3]

    ~projection_kernel_q_h_b1();
};

class projection_kernel_q_h_b2 : public projection_kernel
{
    ClassEngine *m_class_obj;
    double m_z_min;
    double m_z_max;
    double m_M_min;
    double m_M_max;
    double m_N_h;
    double m_n_h;

    Linear_interp_1D m_q_h_b2_z_interp_array;

public:
    projection_kernel_q_h_b2(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max);

    double evaluate(const double &z);
    double get_N_h(); // Number of halos
    double get_n_h(); // number density of halos [Mpc^-3]

    ~projection_kernel_q_h_b2();
};

class projection_kernel_q_h_bs2 : public projection_kernel
{
    ClassEngine *m_class_obj;
    double m_z_min;
    double m_z_max;
    double m_M_min;
    double m_M_max;
    double m_N_h;
    double m_n_h;

    Linear_interp_1D m_q_h_bs2_z_interp_array;

public:
    projection_kernel_q_h_bs2(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max);

    double evaluate(const double &z);
    double get_N_h(); // Number of halos
    double get_n_h(); // number density of halos [Mpc^-3]

    ~projection_kernel_q_h_bs2();
};

double T17_box_correction(const double &k, const double &z, ClassEngine *class_obj);

double T17_shell_correction(const double &k, ClassEngine *class_obj);

double T17_resolution_correction(const double &ell);

double evaluate_Pk_shell_correction_integrand(const double &k_parallel, const double &k_perpendicular, const double &z, const double &delta_r, ClassEngine *class_obj);

struct params_Pk_shell_correction_integrand { double k_perpendicular; double z; double delta_r; ClassEngine *class_obj;};

double Pk_shell_correction_qag_integrand(double k_parallel, void *params);

double Pk_shell_correction_qag(const double &k_perpendicular, const double &z, const double &delta_r, ClassEngine *class_obj);

// ######################################################################################

// PT kernels

double F2_EdS(const double &k_1, const double &k_2, const double &k_3);

double F2_EdS_angular(const double &k_1, const double &k_2, const double &cos_phi_12); /* F2 kernel in EdS cosmology */

double kF2_EdS(const double &k_1, const double &k_2, const double &k_3);

double kF2_EdS_angular(const double &k_1, const double &k_2, const double &cos_phi_12);

double S2(const double &k_1, const double &k_2, const double &k_3); /* tidal shear bias kernel*/

double S2_angular(const double &cos_phi_12);

// ######################################################################################

// 3D wavevector (q) arithmetic utility functions

std::vector<double> spherical_to_cartesian_3D(const std::vector<double> &spherical3D_vec);

std::vector<double> cartesian_to_spherical_3D(const std::vector<double> &cartesian3D_vec);

// the first 3 components of the q vectors should be the Cartesian coordinates and the 4th component should be its length i.e. [q_x,q_y,q_z,q]

std::vector<double> q_mA_vec(const std::vector<double> &q_A_vec);

std::vector<double> q_mA_perturbed_vec(const std::vector<double> &q_A_vec);

std::vector<double> q_ApB_vec(const std::vector<double> &q_A_vec, const std::vector<double> &q_B_vec);

double cos_phi_AB(const std::vector<double> &q_A_vec, const std::vector<double> &q_B_vec);

bool is_vec_alright(const std::vector<double> &q_A_vec);

// ######################################################################################

double alpha_EdS(const double &k_1, const double &k_2, const double &cos_phi_12);

double beta_EdS(const double &k_1, const double &k_2, const double &cos_phi_12);

double F2(const double &q_1, const double &q_2, const double &cos_phi_12);

double G2(const double &q_1, const double &q_2, const double &cos_phi_12);

double F2(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec);

double G2(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec);

double F3(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec);

double G3(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec);

double F4(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec, const std::vector<double> &q_4_vec);

double G4(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec, const std::vector<double> &q_4_vec);

// symmetrized kernels

double F2_sym(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec);

double G2_sym(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec);

double F3_sym(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec);

double F3_sym(const double &q_1, const double &q_2, const double &q_3, const double &mu_12, const double &mu_13, const double &mu_23);

double G3_sym(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec);

double F4_sym(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec, const std::vector<double> &q_4_vec);

double F4_sym(const double &q_1, const double &q_2, const double &q_3, const double &q_4,
              const double &mu_12, const double &mu_13, const double &mu_14, const double &mu_23, const double &mu_24, const double &mu_34);

double G4_sym(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec, const std::vector<double> &q_4_vec);

#endif // COSMOLOGY_UTILS_H
