#ifndef REAL_SPACE_2D_H
#define REAL_SPACE_2D_H

#include <vector>
#include <string>

void C_ells_spherical_sky(const std::vector<double> &l_array, const std::vector<double> &P, std::vector<int> &ells, std::vector<double> &C_ells);

double G_ell_2_x_p(const int &l, const double &x);

double G_ell_2_x_m(const int &l, const double &x);

double xi_theta(const double &theta, const std::vector<int> &ells, const std::vector<double> &C_ells);

double xip_theta(const double &theta, const std::vector<int> &ells, const std::vector<double> &C_ells);

double xim_theta(const double &theta, const std::vector<int> &ells, const std::vector<double> &C_ells);

std::vector<double> xi_theta_array(std::vector<double> &alpha_angles_arcmins, std::vector<double> &l_array, std::vector<double> &P,
                                   const std::string &signal_type, int thread_count);

double P_ell_x_bin_averaged(const int &l, const double &x_min, const double &x_max, const double *P_l_x_min_array, const double *P_l_x_max_array);

double G_ell_2_x_p_bin_averaged(const int &l, const double &x_min, const double &x_max,
                                const double *P_l_x_min_array, const double* ddx_P_l_x_min_array, const double *P_l_x_max_array, const double* ddx_P_l_x_max_array);

double G_ell_2_x_m_bin_averaged(const int &l, const double &x_min, const double &x_max,
                                const double *P_l_x_min_array, const double* ddx_P_l_x_min_array, const double *P_l_x_max_array, const double* ddx_P_l_x_max_array);

double xi_theta_bin_averaged(const double &theta_min, const double &theta_max, const std::vector<int> &ells, const std::vector<double> &C_ells);

double xip_theta_bin_averaged(const double &theta_min, const double &theta_max, const std::vector<int> &ells, const std::vector<double> &C_ells);

double xim_theta_bin_averaged(const double &theta_min, const double &theta_max, const std::vector<int> &ells, const std::vector<double> &C_ells);

std::vector<double> xi_theta_array_bin_averaged(const std::vector<double> &alpha_min_angles_arcmins, const std::vector<double> &alpha_max_angles_arcmins,
                                                const std::vector<double> &l_array, const std::vector<double> &P, const std::string &signal_type, int thread_count);

#endif // REAL_SPACE_2D_H

