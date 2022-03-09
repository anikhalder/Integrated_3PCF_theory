#include <cosmology_utils.h>
#include <integration_utils.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <constants.h>
#include <halo_utils.hpp>

// ######################################################################################

// table utility functions

std::vector<double> read_1_column_table(const std::string &file_name)
{
    std::vector<double> matrix;

    std::ifstream file( file_name );

    double c;

    for( std::string line; getline( file, line ); )
    {
        std::istringstream iss(line);

        iss >> c;

        matrix.push_back(c);
    }

    return matrix;
}

std::vector<std::vector<double>> read_2_column_table(const std::string &file_name)
{
    // where each column is a vector of interest 
    std::vector<std::vector<double>> matrix(2);

    std::ifstream file( file_name );

    double c[2];

    for( std::string line; getline( file, line ); )
    {
        std::istringstream iss(line);

        iss >> c[0] >> c[1];

        matrix.at(0).push_back(c[0]);
        matrix.at(1).push_back(c[1]);
    }

    return matrix;
}

std::vector<std::vector<double> > read_n_column_table(const std::string &file_name, const size_t &n)
{
    // where each column is a vector of interest 
    std::vector<std::vector<double>> matrix(n);

    std::ifstream file( file_name );

    double c[n];

    for( std::string line; getline( file, line ); )
    {
        std::istringstream iss(line);

        for (size_t i=0; i<n; i++)
        {
            iss >> c[i];
            matrix.at(i).push_back(c[i]);
        }            
    }

    return matrix;
}

void normalise_nofz(std::vector<std::vector<double>> &matrix)
{
    double sum = 0, dz = 0;

    for (size_t idx = 0; idx < matrix.at(0).size()-1; idx++)
    {
        dz = matrix.at(0).at(idx+1) - matrix.at(0).at(idx);
        sum += matrix.at(1).at(idx)*dz;
    }

    for (size_t idx = 0; idx < matrix.at(0).size(); idx++)
         matrix.at(1).at(idx) =  matrix.at(1).at(idx) / sum;
}

// ######################################################################################

size_t num_correlations(const size_t &n, const size_t &r)
{
    if (n==1)
        return 1;

    // nCr with replacement
    return tgamma(n+r-1+1)/tgamma(n-1+1)/tgamma(r+1); // Note: tgamma(x+1) = factorial(x)
}

// ######################################################################################

// 2D window functions

double spherical_cap_sqdeg_2_radius(const double &patch_area_sq_degrees)
{
    return acos(1-patch_area_sq_degrees*M_PI/(2.0*180*180));
}

double spherical_cap_radius_2_sqradians(const double &theta_T)
{
    // patch is basically a spherical cap centered around the North Pole on a sphere of unit radius
    double r = 1.0;
    return 2*M_PI*r*r*(1-cos(theta_T));
}

double W2D_TH_RS(const double &length, const double &theta_T)
{
    // normalised 2D tophat (real space)
    if (length <= theta_T)
        return 1/(M_PI*theta_T*theta_T);
    else
        return 0;
}

double W2D_TH_RS_unnormalised(const double &length, const double &theta_T)
{
    // unnormalised 2D tophat (real space)
    if (length <= theta_T)
        return 1;
    else
        return 0;
}

double W2D_TH_FS(const double &l, const double &theta_T)
{
    // normalised 2D tophat (Fourier space)
    if (l == 0.0)
        return 1;
    else
    {
        double x = l*theta_T;
        return j1(x)/(x)*2; //*M_PI*theta_T*theta_T
    }
}

double W2D_U_RS(const double &length, const double &theta_T)
{
    // see eqn (10) of https://www.aanda.org/articles/aa/pdf/2005/40/aa3531-05.pdf
    double x = length*length;
    double y = 1/(2*theta_T*theta_T);
    return y/M_PI * (1-x*y)* exp(-x*y);
}


double W2D_U_FS(const double &l, const double &theta_T)
{
    // see the line below eqn (11) of https://www.aanda.org/articles/aa/pdf/2005/40/aa3531-05.pdf
    double x = l*theta_T;
    return x*x/2 * exp(-x*x/2);
}

// ######################################################################################

// 3D window functions

double W3D_TH_FS(const double &k, const double &R)
{
    // normalised 3D tophat (Fourier space)
    // see eqn (A.114) of Dodelson & Schmidt (2020)
    // k in [1/Mpc] and R in [Mpc] --> both are in comoving coordinates
    double x = k*R;
    return -3./(x*x*x)*(-x*cos(x) + sin(x));
}

double W3D_prime_TH_FS(const double &k, const double &R)
{
    // derivative of normalised 3D tophat wrt R (Fourier space)
    // d W(kR) / d R = (d W(kR) / d kR)*(d kR / d R) = k * (d W(kR) / d kR)
    // k in [1/Mpc] and R in [Mpc] --> both are in comoving coordinates
    double x = k*R;
    return k*(-9./(x*x*x*x)*(-x*cos(x) + sin(x)) + 3.*sin(x)/(x*x));
}

// ######################################################################################

// 2D wavevector (l) arithmetic utility functions

double l_ApB(const double &l_A, const double &phi_A, const double &l_B, const double &phi_B)
{
    double val = sqrt(l_A*l_A+l_B*l_B+2*l_A*l_B*cos(phi_A-phi_B));

    if (isnan(val))
        return 0;

    return val;
}

double phi_ApB(const double &l_A, const double &phi_A, const double &l_B, const double &phi_B)
{
    double x = l_A*cos(phi_A) + l_B*cos(phi_B);
    double y = l_A*sin(phi_A) + l_B*sin(phi_B);

    double phi = atan2(y,x);
    if (phi < 0)
        phi += 2*M_PI;

    return phi;
}

double l_ApBpC(const double &l_A, const double &phi_A, const double &l_B, const double &phi_B, const double &l_C, const double &phi_C)
{
    double val = sqrt(l_A*l_A+l_B*l_B+l_C*l_C
                +2*l_A*l_B*cos(phi_A-phi_B)
                +2*l_A*l_C*cos(phi_A-phi_C)
                +2*l_B*l_C*cos(phi_B-phi_C));

    if (isnan(val))
        return 0;

    return val;
}

double phi_ApBpC(const double &l_A, const double &phi_A, const double &l_B, const double &phi_B, const double &l_C, const double &phi_C)
{
    double l_sum = l_ApB(l_A,phi_A,l_B,phi_B);
    double phi_sum = phi_ApB(l_A,phi_A,l_B,phi_B);

    return phi_ApB(l_sum,phi_sum,l_C,phi_C);
}

// ######################################################################################

// cosmology utility functions

/* Notation:
 * chi or x --> comoving distance (denoted also as 'w' in literature)
 * eta --> conformal time (denoted also as 'tau' in literature) */

double proper_time_integrand(double zp, void *params)
{
    params_proper_time_integrand *p = static_cast<params_proper_time_integrand *>(params);

    return 1./(1.+zp)/p->class_obj->get_H_z(zp);
}

double proper_time(const double &z, ClassEngine *class_obj)
{
    // t = \int_{0}^{t} dt' = \int_{+inf}^{z} dz' a(z') d\eta / dz' = \int_{z}^{+inf} H(z') / (1+z')   where the last equality holds since d\eta / dz' = -H(z)

    params_proper_time_integrand args = {class_obj};

    double age = 0, error = 0;

    qagiu_1D_integration(&proper_time_integrand, static_cast<void *>(&args), z, 5*calls_1e3, age, error);

    return age; // in [Mpc]
}

double sigma_squared_R_z_integrand(double k, void *params)
{
    params_sigma_squared_R_z_integrand *p = static_cast<params_sigma_squared_R_z_integrand *>(params);

    double R = p->R;
    double z = p->z;
    double P_L = p->class_obj->pk_lin(k,z);

    return k*k*P_L*W3D_TH_FS(k,R)*W3D_TH_FS(k,R);
}

double sigma_squared_R_z(const double &R, const double &z,  ClassEngine *class_obj)
{
    // see eqn (12.4) of Dodelson & Schmidt (2020)
    // sigma^2(R,z) = 1 / (2 pi^2) \int dk k^2 P_L(k,z) W(kR)^2

    params_sigma_squared_R_z_integrand args = {R, z, class_obj};

    double sigma_squared = 0, error = 0;

    qag_1D_integration_abs_rel(&sigma_squared_R_z_integrand, static_cast<void *>(&args), 1e-5, class_obj->get_k_max_pk(), 5*calls_1e3, sigma_squared, error);

    return sigma_squared / (2.*M_PI*M_PI); // dimensionless --> square root of this quantity gives sigma(R,z) i.e. similar to what's returned e.g. by CLASS
}

double sigma_squared_prime_R_z_integrand(double k, void *params)
{
    params_sigma_squared_R_z_integrand *p = static_cast<params_sigma_squared_R_z_integrand *>(params);

    double R = p->R;
    double z = p->z;
    double P_L = p->class_obj->pk_lin(k,z);

    double W = W3D_TH_FS(k,R);
    double W_prime = W3D_prime_TH_FS(k,R);

    return k*k*P_L*2.*W*W_prime;

//    if (W >= 0)
//        return k*k*P_L*2.*abs(W)*W_prime;
//    else
//        return -k*k*P_L*2.*abs(W)*W_prime;
}

double sigma_squared_prime_R_z(const double &R, const double &z,  ClassEngine *class_obj)
{
    // d sigma^2(R,z) / d R = 1 / (2 pi^2) \int dk k^2 P_L(k,z) (d W(kR)^2 / d R)
    //                      = 1 / (2 pi^2) \int dk k^2 P_L(k,z) 2 W(kR) (d W(kR) / d R)

    params_sigma_squared_R_z_integrand args = {R, z, class_obj};

    double sigma_squared_prime = 0, error = 0;

    qag_1D_integration_abs_rel(&sigma_squared_prime_R_z_integrand, static_cast<void *>(&args), 1e-5, class_obj->get_k_max_pk(), 5*calls_1e3, sigma_squared_prime, error);

    return sigma_squared_prime / (2.*M_PI*M_PI); // in [1/Mpc] --> currenly, gives -ve of the result returned by CLASS; TODO: FIX THIS!!!
}

double pk_lin_manual(const double &k, const double &z, ClassEngine *class_obj, Linear_interp_1D& Tk_z0_d_tot)
{
    // TODO: THIS IS INCORRECT!!! FIX THIS!!!
    // Units: k [1/Mpc] and P(k) [Mpc^3]
    // reproduce the same linear power spectrum as given out by class (for a non-running scalar spectral index)
    return 2.*M_PI*M_PI / (k*k*k) *  pow(class_obj->get_D_plus_z(z)*Tk_z0_d_tot.interp(k),2) *
            class_obj->get_A_s() * pow(k/class_obj->get_k_pivot(), class_obj->get_n_s()-1);
}

// ######################################################################################

// 2D projection kernels

double q_m(const double &z, ClassEngine *class_obj, Linear_interp_1D *n_m_of_z)
{
    /* q_m(x) = n_m(z) dz/dx
     * => q_m(x(z)) = n_m(z) H(z) */

    return n_m_of_z->interp(z)*class_obj->get_H_z(z);
}

double W_k_zs_distribution_integrand(double zs, void *params)
{
    /* dz'/H(z') g(x(z')) (x(z')-x(z))/x(z')
     * = dz'/H(z') n_source(z') H(z') (x(z')-x(z))/x(z') , where g(x') = n_source(z') H(z')
     * = dz' n_source(z') (x(z')-x(z))/x(z')
     * Here, denote z' as zs */

    params_W_k_zs_distribution_integrand *p = static_cast<params_W_k_zs_distribution_integrand *>(params);

    if (p->z > zs)
        return 0;

    double chi_z = p->class_obj->get_chi_z(p->z);
    double chi_zs = p->class_obj->get_chi_z(zs);

    if (p->delta_photoz == 0.0)
        return p->n_source_of_z->interp(zs) * (chi_zs - chi_z) / chi_zs;
    else
        return p->n_source_of_z->interp(zs + p->delta_photoz) * (chi_zs - chi_z) / chi_zs;
}

double q_k_zs_distribution(const double &z, ClassEngine *class_obj, Linear_interp_1D *n_source_of_z, const double &z_max)
{

    /* equations (10.42) and (10.41) combined of https://edoc.ub.uni-muenchen.de/23401/1/Friedrich_Oliver.pdf
     * or equations (6.21) and (6.19) combined of https://arxiv.org/pdf/astro-ph/9912508.pdf
     *
     * q_k(x) = 3/2 H0^2 Omega0_m x/a(x) W_k(x) , where W_k(x) = \int_x^x_max dx'g(x') (x'-x)/x'
     * => q_k(x) = 3/2 H0^2 Omega0_m x/a(x) \int_x^x_max dx'g(x') (x'-x)/x'
     * => q_k(x(z)) = 3/2 H0^2 Omega0_m (1+z) x(z) \int_z^z_max dz'/H(z') g(x(z')) (x(z')-x(z))/x(z')
     * Here, denote z' as zs */

    //parameters in integrand
    params_W_k_zs_distribution_integrand args = {z, class_obj, n_source_of_z, 0.0};

    double W_k_z = 0, error = 0;

    qag_1D_integration_abs_rel(&W_k_zs_distribution_integrand, static_cast<void *>(&args), z, z_max, 5*calls_1e3, W_k_z, error);

    double chi_z = class_obj->get_chi_z(z);
    double H_0 = class_obj->get_H_z(0);

    return 3./2. * H_0 * H_0 * class_obj->get_Omega0_m() * (1.+z) * chi_z * W_k_z;
}

double q_k_zs_distribution_f_IA_NLA_delta_photoz(const double &z, ClassEngine *class_obj, Linear_interp_1D *n_source_of_z, const double &z_max, 
                                                 const double &delta_photoz, const double &A_IA_0_NLA, const double &alpha_IA_0_NLA)
{

    /* equations (10.42) and (10.41) combined of https://edoc.ub.uni-muenchen.de/23401/1/Friedrich_Oliver.pdf
     * or equations (6.21) and (6.19) combined of https://arxiv.org/pdf/astro-ph/9912508.pdf
     *
     * q_k(x) = 3/2 H0^2 Omega0_m x/a(x) W_k(x) , where W_k(x) = \int_x^x_max dx'g(x') (x'-x)/x'
     * => q_k(x) = 3/2 H0^2 Omega0_m x/a(x) \int_x^x_max dx'g(x') (x'-x)/x'
     * => q_k(x(z)) = 3/2 H0^2 Omega0_m (1+z) x(z) \int_z^z_max dz'/H(z') g(x(z')) (x(z')-x(z))/x(z')
     * Here, denote z' as zs */

    //parameters in integrand
    params_W_k_zs_distribution_integrand args = {z, class_obj, n_source_of_z, delta_photoz};

    double W_k_z = 0, error = 0;

    qag_1D_integration_abs_rel(&W_k_zs_distribution_integrand, static_cast<void *>(&args), z, z_max, 5*calls_1e3, W_k_z, error);

    double chi_z = class_obj->get_chi_z(z);
    double H_0 = class_obj->get_H_z(0);

    double q_k_z = 3./2. * H_0 * H_0 * class_obj->get_Omega0_m() * (1.+z) * chi_z * W_k_z;

    if (A_IA_0_NLA == 0.0)
        return q_k_z;
    else
    {
        double f_IA_NLA_z = - A_IA_0_NLA*pow((1.+z)/1.62,alpha_IA_0_NLA)*0.0134; //*class_obj->get_D_plus_z(0)/class_obj->get_D_plus_z(z);
        return  q_k_z + f_IA_NLA_z*n_source_of_z->interp(z+delta_photoz)*class_obj->get_H_z(z);
    }
}

double q_k_zs_fixed(const double &z, ClassEngine *class_obj, const double &zs)
{
    /* q_k(x) = 3/2 H0^2 Omega0_m x/a(x) (x'-x)/x'
     * => q_k(x(z)) = 3/2 H0^2 Omega0_m (1+z) x(z) (x(z')-x(z))/x(z')
     * Here, denote z' as zs */

    if (z > zs)
        return 0;

    double chi_z = class_obj->get_chi_z(z);
    double chi_zs = class_obj->get_chi_z(zs);
    double H_0 = class_obj->get_H_z(0);

    return 3./2. * H_0 * H_0 * class_obj->get_Omega0_m() * (1.+z) * chi_z * (chi_zs - chi_z) / chi_zs;
}

double cmb_lensing_convergence_kernel(const double &z, ClassEngine *class_obj)
{
    // eqn (23) of https://arxiv.org/pdf/1511.05534.pdf

    if (z > z_cmb)
        return 0;

    double chi_z = class_obj->get_chi_z(z);
    double chi_z_cmb = class_obj->get_chi_z(z_cmb);
    double H_0 = class_obj->get_H_z(0);

    return 3./2. * H_0 * H_0 * class_obj->get_Omega0_m() * (1.+z) * chi_z * (chi_z_cmb - chi_z) / chi_z_cmb; // / class_obj->get_Hz(z)  (should this H(z) be here?)
}

double q_h_b1_integrand(double M, void *params)
{
    params_q_h_integrand *p = static_cast<params_q_h_integrand *>(params);

    // Use PBS b1 and PS HMF
    //double b1 = halo_b1_PBS(M, p->z, p->class_obj);
    //return dn_dM_PS(M, p->z, p->class_obj)*b1;

    // Use Tinker b1 fitting function and Tinker HMF
    double b1 = halo_b1_Tinker2010(M, p->z, p->class_obj);
    return dn_dM_Tinker2008(M, p->z, p->class_obj)*b1;
}

double q_h_b1_unnormalised(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max)
{
    // eqns (11) and (15) of https://iopscience.iop.org/article/10.3847/1538-4357/aa943d/pdf

    // M_min and M_mass must be in [M_sun] units
    //parameters in integrand
    params_q_h_integrand args = {z, class_obj};

    double W_h_z = 0, error = 0;

    qag_1D_integration_abs_rel(&q_h_b1_integrand, static_cast<void *>(&args), M_min, M_max, 5*calls_1e3, W_h_z, error);

    double chi_z = class_obj->get_chi_z(z);

    return 4.*M_PI*chi_z*chi_z*W_h_z; // this is an unnormalised quantity --> divide by the total number of halos in bin, N_h, to normalise
}

double q_h_b2_integrand(double M, void *params)
{
    params_q_h_integrand *p = static_cast<params_q_h_integrand *>(params);

    // Use PBS b2 and PS HMF
    //double b2 = halo_b2_PBS(M, p->z, p->class_obj);
    //return dn_dM_PS(M, p->z, p->class_obj)*b2;

    // Use Lazeyras b2 fitting function (with Tinker b1) and Tinker HMF
    double b2 = halo_b2_Lazeyras(M, p->z, p->class_obj);
    return dn_dM_Tinker2008(M, p->z, p->class_obj)*b2;
}

double q_h_b2_unnormalised(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max)
{
    // eqns (11) and (15) of https://iopscience.iop.org/article/10.3847/1538-4357/aa943d/pdf

    // M_min and M_mass must be in [M_sun] units
    //parameters in integrand
    params_q_h_integrand args = {z, class_obj};

    double W_h_z = 0, error = 0;

    qag_1D_integration_abs_rel(&q_h_b2_integrand, static_cast<void *>(&args), M_min, M_max, 5*calls_1e3, W_h_z, error);

    double chi_z = class_obj->get_chi_z(z);

    return 4.*M_PI*chi_z*chi_z*W_h_z; // this is an unnormalised quantity --> divide by the total number of halos in bin, N_h, to normalise
}

double q_h_bs2_integrand(double M, void *params)
{
    params_q_h_integrand *p = static_cast<params_q_h_integrand *>(params);

    // Use bs2 coevolution relation (with Tinker b1) and Tinker HMF
    double bs2 = halo_bs2_coevolution(M, p->z, p->class_obj);
    return dn_dM_Tinker2008(M, p->z, p->class_obj)*bs2;
}

double q_h_bs2_unnormalised(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max)
{
    // eqns (11) and (15) of https://iopscience.iop.org/article/10.3847/1538-4357/aa943d/pdf

    // M_min and M_mass must be in [M_sun] units
    //parameters in integrand
    params_q_h_integrand args = {z, class_obj};

    double W_h_z = 0, error = 0;

    qag_1D_integration_abs_rel(&q_h_bs2_integrand, static_cast<void *>(&args), M_min, M_max, 5*calls_1e3, W_h_z, error);

    double chi_z = class_obj->get_chi_z(z);

    return 4.*M_PI*chi_z*chi_z*W_h_z; // this is an unnormalised quantity --> divide by the total number of halos in bin, N_h, to normalise
}

double q_h_without_b_integrand(double M, void *params)
{
    params_q_h_integrand *p = static_cast<params_q_h_integrand *>(params);

    return dn_dM_Tinker2008(M, p->z, p->class_obj);
}

double q_h_without_b_unnormalised(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max)
{
    // eqns (11) and (15) of https://iopscience.iop.org/article/10.3847/1538-4357/aa943d/pdf

    // M_min and M_mass must be in [M_sun] units
    //parameters in integrand
    params_q_h_integrand args = {z, class_obj};

    double W_h_z = 0, error = 0;

    qag_1D_integration_abs_rel(&q_h_without_b_integrand, static_cast<void *>(&args), M_min, M_max, 5*calls_1e3, W_h_z, error);

    double chi_z = class_obj->get_chi_z(z);

    return 4.*M_PI*chi_z*chi_z*W_h_z; // this is an unnormalised quantity --> divide by the total number of halos in bin, N_h, to normalise
}

double n_h_z(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max)
{
    // M_min and M_mass must be in [M_sun] units
    //parameters in integrand
    params_q_h_integrand args = {z, class_obj};

    double W_h_z = 0, error = 0;

    qag_1D_integration_abs_rel(&q_h_without_b_integrand, static_cast<void *>(&args), M_min, M_max, 5*calls_1e3, W_h_z, error);

    return W_h_z;
}


projection_kernel::~projection_kernel()
{

}

projection_kernel_q_m::projection_kernel_q_m(ClassEngine *class_obj, Linear_interp_1D *n_m_of_z) : m_class_obj(class_obj), m_n_m_of_z(n_m_of_z)
{

}

double projection_kernel_q_m::evaluate(const double &z)
{
    return q_m(z, m_class_obj, m_n_m_of_z);
}

projection_kernel_q_m::~projection_kernel_q_m()
{

}

projection_kernel_q_k_zs_fixed::projection_kernel_q_k_zs_fixed(ClassEngine *class_obj, const double &zs) : m_class_obj(class_obj), m_zs(zs)
{
    std::vector<double> m_z_array;
    std::vector<double> m_q_k_z_array;

    double z=0;
    while (z <= m_zs)
    {
      m_z_array.push_back(z);
      m_q_k_z_array.push_back(q_k_zs_fixed(z, m_class_obj, m_zs));

      z += delta_z_step;
  //    z += 0.000025;
    }

    m_q_k_zs_fixed_z_array = Linear_interp_1D(m_z_array, m_q_k_z_array);
}

double projection_kernel_q_k_zs_fixed::evaluate(const double &z)
{
    return q_k_zs_fixed(z, m_class_obj, m_zs);
    //return m_q_k_zs_fixed_z_array.interp(z);
}

projection_kernel_q_k_zs_fixed::~projection_kernel_q_k_zs_fixed()
{

}

projection_kernel_q_k_zs_distribution::projection_kernel_q_k_zs_distribution(ClassEngine *class_obj, Linear_interp_1D *n_source_of_z, const double &z_max) :
    m_class_obj(class_obj), m_n_source_of_z(n_source_of_z), m_z_max(z_max), m_delta_photoz(0.0), m_A_IA_0_NLA(0.0), m_alpha_IA_0_NLA(0.0)
{
    std::vector<double> m_z_array;
    std::vector<double> m_q_k_z_array;

    double z=0;
    while (z <= m_z_max)
    {
      m_z_array.push_back(z);
      m_q_k_z_array.push_back(q_k_zs_distribution(z, m_class_obj, m_n_source_of_z, m_z_max));

      //z += 0.05;
      z += delta_z_step;
      //z += 0.001;
    }

    m_q_k_zs_distribution_z_array = Linear_interp_1D(m_z_array, m_q_k_z_array);
}

projection_kernel_q_k_zs_distribution::projection_kernel_q_k_zs_distribution(ClassEngine *class_obj, Linear_interp_1D *n_source_of_z, const double &z_max, const double &delta_photoz, const double &A_IA_0_NLA, const double &alpha_IA_0_NLA) :
    m_class_obj(class_obj), m_n_source_of_z(n_source_of_z), m_z_max(z_max), m_delta_photoz(delta_photoz), m_A_IA_0_NLA(A_IA_0_NLA), m_alpha_IA_0_NLA(alpha_IA_0_NLA)
{
    std::vector<double> m_z_array;
    std::vector<double> m_q_k_z_array;

    double z=0;
    while (z <= m_z_max)
    {
      m_z_array.push_back(z);
      m_q_k_z_array.push_back(q_k_zs_distribution_f_IA_NLA_delta_photoz(z, m_class_obj, m_n_source_of_z, m_z_max, m_delta_photoz, m_A_IA_0_NLA, m_alpha_IA_0_NLA));

      z += delta_z_step;
    }

    m_q_k_zs_distribution_z_array = Linear_interp_1D(m_z_array, m_q_k_z_array);
}

double projection_kernel_q_k_zs_distribution::evaluate(const double &z)
{
    //return q_k_zs_distribution(z, m_class_obj, m_n_source_of_z, m_z_max);
    return m_q_k_zs_distribution_z_array.interp(z);
}

projection_kernel_q_k_zs_distribution::~projection_kernel_q_k_zs_distribution()
{

}

projection_kernel_q_h::projection_kernel_q_h(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max) :
    m_class_obj(class_obj), m_z_min(z_min), m_z_max(z_max), m_M_min(M_min), m_M_max(M_max)
{
    std::vector<double> m_z_array;
    std::vector<double> m_q_h_without_b_z_array;
    std::vector<double> m_n_h_z_array;

    m_N_h = N_h_hcubature(class_obj, m_z_min, m_z_max, m_M_min, m_M_max);

    double chi_z_max = class_obj->get_chi_z(z_max);
    double chi_z_min = class_obj->get_chi_z(z_min);

    double V = 4.*M_PI/3.*( chi_z_max*chi_z_max*chi_z_max - chi_z_min*chi_z_min*chi_z_min );
    m_n_h = m_N_h / V ; // [Mpc^-3]

    double z = m_z_min;
    while (z <= m_z_max)
    {
      m_z_array.push_back(z);
      m_q_h_without_b_z_array.push_back(q_h_without_b_unnormalised(z, m_class_obj, m_M_min, m_M_max)/m_N_h);
      m_n_h_z_array.push_back(n_h_z(z, m_class_obj, m_M_min, m_M_max));

      z += delta_z_step;
    }

    m_q_h_without_b_z_interp_array = Linear_interp_1D(m_z_array, m_q_h_without_b_z_array);
    m_n_h_z_interp_array = Linear_interp_1D(m_z_array, m_n_h_z_array);
}

double projection_kernel_q_h::evaluate(const double &z)
{
    //return q_h_without_b_unnormalised(z, m_class_obj, m_M_min, m_M_max)/m_N_h;

    if (z < m_z_min || z > m_z_max)
        return 0;
    else
        return m_q_h_without_b_z_interp_array.interp(z);
}

double projection_kernel_q_h::get_N_h()
{
    return m_N_h;
}

double projection_kernel_q_h::get_n_h()
{
    return m_n_h;
}

double projection_kernel_q_h::get_n_h_z(const double &z)
{
    //return n_h_z(z, m_class_obj, m_M_min, m_M_max);

    if (z < m_z_min || z > m_z_max)
        return 0;
    else
        return m_n_h_z_interp_array.interp(z);
}

projection_kernel_q_h::~projection_kernel_q_h()
{

}

projection_kernel_q_h_b1::projection_kernel_q_h_b1(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max) :
    m_class_obj(class_obj), m_z_min(z_min), m_z_max(z_max), m_M_min(M_min), m_M_max(M_max)
{
    std::vector<double> m_z_array;
    std::vector<double> m_q_h_b1_z_array;

    m_N_h = N_h_hcubature(class_obj, m_z_min, m_z_max, m_M_min, m_M_max);

    double chi_z_max = class_obj->get_chi_z(z_max);
    double chi_z_min = class_obj->get_chi_z(z_min);

    double V = 4.*M_PI/3.*( chi_z_max*chi_z_max*chi_z_max - chi_z_min*chi_z_min*chi_z_min );
    m_n_h = m_N_h / V ; // [Mpc^-3]

    double z = m_z_min;
    while (z <= m_z_max)
    {
      m_z_array.push_back(z);
      m_q_h_b1_z_array.push_back(q_h_b1_unnormalised(z, m_class_obj, m_M_min, m_M_max)/m_N_h);

      z += delta_z_step;
    }

    m_q_h_b1_z_interp_array = Linear_interp_1D(m_z_array, m_q_h_b1_z_array);
}

double projection_kernel_q_h_b1::evaluate(const double &z)
{
    //return q_h_b1_unnormalised(z, m_class_obj, m_M_min, m_M_max)/m_N_h;

    if (z < m_z_min || z > m_z_max)
        return 0;
    else
        return m_q_h_b1_z_interp_array.interp(z);
}

double projection_kernel_q_h_b1::get_N_h()
{
    return m_N_h;
}

double projection_kernel_q_h_b1::get_n_h()
{
    return m_n_h;
}

projection_kernel_q_h_b1::~projection_kernel_q_h_b1()
{

}

projection_kernel_q_h_b2::projection_kernel_q_h_b2(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max) :
    m_class_obj(class_obj), m_z_min(z_min), m_z_max(z_max), m_M_min(M_min), m_M_max(M_max)
{
    std::vector<double> m_z_array;
    std::vector<double> m_q_h_b2_z_array;

    m_N_h = N_h_hcubature(class_obj, m_z_min, m_z_max, m_M_min, m_M_max);

    double chi_z_max = class_obj->get_chi_z(z_max);
    double chi_z_min = class_obj->get_chi_z(z_min);

    double V = 4.*M_PI/3.*( chi_z_max*chi_z_max*chi_z_max - chi_z_min*chi_z_min*chi_z_min );
    m_n_h = m_N_h / V ; // [Mpc^-3]

    double z = m_z_min;
    while (z <= m_z_max)
    {
      m_z_array.push_back(z);
      m_q_h_b2_z_array.push_back(q_h_b2_unnormalised(z, m_class_obj, m_M_min, m_M_max)/m_N_h);

      z += delta_z_step;
    }

    m_q_h_b2_z_interp_array = Linear_interp_1D(m_z_array, m_q_h_b2_z_array);
}

double projection_kernel_q_h_b2::evaluate(const double &z)
{
    //return q_h_b2_unnormalised(z, m_class_obj, m_M_min, m_M_max)/m_N_h;

    if (z < m_z_min || z > m_z_max)
        return 0;
    else
        return m_q_h_b2_z_interp_array.interp(z);
}

double projection_kernel_q_h_b2::get_N_h()
{
    return m_N_h;
}

double projection_kernel_q_h_b2::get_n_h()
{
    return m_n_h;
}

projection_kernel_q_h_b2::~projection_kernel_q_h_b2()
{

}

projection_kernel_q_h_bs2::projection_kernel_q_h_bs2(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max) :
    m_class_obj(class_obj), m_z_min(z_min), m_z_max(z_max), m_M_min(M_min), m_M_max(M_max)
{
    std::vector<double> m_z_array;
    std::vector<double> m_q_h_bs2_z_array;

    m_N_h = N_h_hcubature(class_obj, m_z_min, m_z_max, m_M_min, m_M_max);

    double chi_z_max = class_obj->get_chi_z(z_max);
    double chi_z_min = class_obj->get_chi_z(z_min);

    double V = 4.*M_PI/3.*( chi_z_max*chi_z_max*chi_z_max - chi_z_min*chi_z_min*chi_z_min );
    m_n_h = m_N_h / V ; // [Mpc^-3]

    double z = m_z_min;
    while (z <= m_z_max)
    {
      m_z_array.push_back(z);
      m_q_h_bs2_z_array.push_back(q_h_bs2_unnormalised(z, m_class_obj, m_M_min, m_M_max)/m_N_h);

      z += delta_z_step;
    }

    m_q_h_bs2_z_interp_array = Linear_interp_1D(m_z_array, m_q_h_bs2_z_array);
}

double projection_kernel_q_h_bs2::evaluate(const double &z)
{
    //return q_h_bs2_unnormalised(z, m_class_obj, m_M_min, m_M_max)/m_N_h;

    if (z < m_z_min || z > m_z_max)
        return 0;
    else
        return m_q_h_bs2_z_interp_array.interp(z);
}

double projection_kernel_q_h_bs2::get_N_h()
{
    return m_N_h;
}

double projection_kernel_q_h_bs2::get_n_h()
{
    return m_n_h;
}

projection_kernel_q_h_bs2::~projection_kernel_q_h_bs2()
{

}


double T17_box_correction(const double &k, const double &z, ClassEngine *class_obj)
{
    // finite box size effect in T17

    double h = class_obj->get_h();
    double chi_z_Mpc_h = class_obj->get_chi_z(z)*h; // Mpc / h

    for (int i = 1; i <= 14 ; i++)
        if ( (chi_z_Mpc_h < (i*delta_r_T17)) && (k < 2*M_PI/(i*delta_r_T17)*h) )
            return 0;

    return 1;
}

double T17_shell_correction(const double &k, ClassEngine *class_obj)
{
    // finite lens-shell-thickness effect in T17

    double h = class_obj->get_h();
    return sqrt(pow(1.+c1_T17*pow(k/h,-a1_T17),a1_T17) / pow(1.+c2_T17*pow(k/h,-a2_T17),a3_T17));
}

double T17_resolution_correction(const double &ell)
{
    // finite angular resolution healpy map effect in T17

    return 1. / ( 1. + ell*ell/(ell_res_T17*ell_res_T17) );
}

double evaluate_Pk_shell_correction_integrand(const double &k_parallel, const double &k_perpendicular, const double &z, const double &delta_r, ClassEngine *class_obj)
{
    // finite lens-shell-thickness effect in T17

    double k = sqrt(k_perpendicular*k_perpendicular + k_parallel*k_parallel);
    double x = k_parallel*delta_r/2.0;
    double sinc_x = sin(x)/x;
    return class_obj->pk_nl(k,z)*sinc_x*sinc_x;
}

double Pk_shell_correction_qag_integrand(double k_parallel, void *params)
{
    params_Pk_shell_correction_integrand *p = static_cast<params_Pk_shell_correction_integrand *>(params);

    return evaluate_Pk_shell_correction_integrand(k_parallel, p->k_perpendicular, p->z, p->delta_r, p->class_obj);
}

double Pk_shell_correction_qag(const double &k_perpendicular, const double &z, const double &delta_r, ClassEngine *class_obj)
{
    double result = 0, error = 0;

    //parameters in integrand
    params_Pk_shell_correction_integrand args = {k_perpendicular, z, delta_r, class_obj};

    double k_parallel_lower = 0.0; // [1/Mpc]
    double k_parallel_upper = class_obj->get_k_max_pk(); // [1/Mpc]

    qag_1D_integration(&Pk_shell_correction_qag_integrand, static_cast<void *>(&args), k_parallel_lower, k_parallel_upper, calls_1e3, result, error);

    return result*delta_r/M_PI;
}

// ######################################################################################

// PT kernels

double F2_EdS(const double &k_1, const double &k_2, const double &k_3)
{
    if(k_1 == 0 || k_2 == 0)
        return 0;

    double cos_phi_12 = 0.5*(k_3*k_3-k_1*k_1-k_2*k_2)/(k_1*k_2); // law of cosines (with a negative sign to get the complementary angle outside the triangle)
    return 5./7. + 0.5*cos_phi_12*(k_1/k_2 + k_2/k_1) + 2./7.*cos_phi_12*cos_phi_12;
}

double F2_EdS_angular(const double &k_1, const double &k_2, const double &cos_phi_12)
{
    if(k_1 == 0 || k_2 == 0)
        return 0;

    //return 0.5*((1+k_1/k_2*cos_phi_12)+(1+k_2/k_1*cos_phi_12)) + 2./7.*(cos_phi_12*cos_phi_12-1);
    return 5./7. + 0.5*cos_phi_12*(k_1/k_2 + k_2/k_1) + 2./7.*cos_phi_12*cos_phi_12;
}

double kF2_EdS(const double &k_1, const double &k_2, const double &k_3)
{
    if(k_1 == 0 || k_2 == 0)
        return 0;

    // multiply k_1*k_2 throughout the F2_EdS kernel expression
    double cos_phi_12 = 0.5*(k_3*k_3-k_1*k_1-k_2*k_2)/(k_1*k_2);
    return 5./7.*k_1*k_2 + 0.5*cos_phi_12*(k_1*k_1 + k_2*k_2) + 2./7.*k_1*k_2*cos_phi_12*cos_phi_12;
}

double kF2_EdS_angular(const double &k_1, const double &k_2, const double &cos_phi_12)
{
    if(k_1 == 0 || k_2 == 0)
        return 0;

    // multiply k_1*k_2 throughout the F2_EdS kernel expression
    //return 0.5*((k_1*k_2+k_1*k_1*cos_phi_12)+(k_1*k_2+k_2*k_2*cos_phi_12)) + 2./7.*k_1*k_2*(cos_phi_12*cos_phi_12-1);
    return 5./7.*k_1*k_2 + 0.5*cos_phi_12*(k_1*k_1 + k_2*k_2) + 2./7.*k_1*k_2*cos_phi_12*cos_phi_12;
}

double S2(const double &k_1, const double &k_2, const double &k_3)
{
    if(k_1 == 0 || k_2 == 0)
        return 0;

    double cos_phi_12 = 0.5*(k_3*k_3-k_1*k_1-k_2*k_2)/(k_1*k_2); // law of cosines
    return cos_phi_12*cos_phi_12 - 1./3.;
}

double S2_angular(const double &cos_phi_12)
{
    return cos_phi_12*cos_phi_12 - 1./3.;
}

// ######################################################################################

// 3D wavevector (q) arithmetic utility functions

std::vector<double> spherical_to_cartesian_3D(const std::vector<double> &spherical3D_vec)
{
    return std::vector<double> {spherical3D_vec.at(0)*sin(spherical3D_vec.at(1))*cos(spherical3D_vec.at(2)),
                                spherical3D_vec.at(0)*sin(spherical3D_vec.at(1))*sin(spherical3D_vec.at(2)),
                                spherical3D_vec.at(0)*cos(spherical3D_vec.at(1)),
                                spherical3D_vec.at(3)};
}

std::vector<double> cartesian_to_spherical_3D(const std::vector<double> &cartesian3D_vec)
{
    double x = cartesian3D_vec.at(0);
    double y = cartesian3D_vec.at(1);
    double z = cartesian3D_vec.at(2);

    double r = cartesian3D_vec.at(3);

    double theta = atan2(sqrt(x*x+y*y), z);

    double phi = atan2(y,x);
    if ( phi < 0)
        phi = 2*M_PI + phi;

    return std::vector<double> {r, theta, phi, r};
}

// the first 3 components of the q_vec should be the Cartesian coordinates and the 4th component should be the vector's length i.e. [q_x,q_y,q_z,q]

// m stands for minus (-) and p stands for plus (+)

std::vector<double> q_mA_vec(const std::vector<double> &q_A_vec)
{
    return std::vector<double> {-q_A_vec.at(0), -q_A_vec.at(1), -q_A_vec.at(2), q_A_vec.at(3)};
}

std::vector<double> q_mA_perturbed_vec(const std::vector<double> &q_A_vec)
{
    double pert_factor = 0.99999;

    double x = -q_A_vec.at(0)*pert_factor;
    double y = -q_A_vec.at(1)*pert_factor;
    double z = -q_A_vec.at(2)*pert_factor;

    double r = sqrt(x*x + y*y + z*z);
    return std::vector<double> {x, y, z, r};
}

std::vector<double> q_ApB_vec(const std::vector<double> &q_A_vec, const std::vector<double> &q_B_vec)
{
    double q_ApB_x = q_A_vec.at(0) + q_B_vec.at(0);
    double q_ApB_y = q_A_vec.at(1) + q_B_vec.at(1);
    double q_ApB_z = q_A_vec.at(2) + q_B_vec.at(2);
    double q_ApB = sqrt(q_ApB_x*q_ApB_x + q_ApB_y*q_ApB_y + q_ApB_z*q_ApB_z);
    return std::vector<double> {q_ApB_x, q_ApB_y, q_ApB_z, q_ApB};
}

double cos_phi_AB(const std::vector<double> &q_A_vec, const std::vector<double> &q_B_vec)
{
    return (q_A_vec.at(0)*q_B_vec.at(0) + q_A_vec.at(1)*q_B_vec.at(1) + q_A_vec.at(2)*q_B_vec.at(2) ) / (q_A_vec.at(3)*q_B_vec.at(3));
}

bool is_vec_alright(const std::vector<double> &q_A_vec)
{
    if (isnan(q_A_vec.at(0)) || isnan(q_A_vec.at(1)) || isnan(q_A_vec.at(2)))
        return false;

    else
        return true;
}

// ######################################################################################

double alpha_EdS(const double &k_1, const double &k_2, const double &cos_phi_12)
{
    if(k_1 == 0)
        return 0;

    return 1.0 + k_2/k_1*cos_phi_12;
}

double beta_EdS(const double &k_1, const double &k_2, const double &cos_phi_12)
{
    if(k_1 == 0 || k_2 == 0)
        return 0;

    return 0.5*cos_phi_12*(k_1/k_2 + k_2/k_1 + 2*cos_phi_12);
}

double F2(const double &q_1, const double &q_2, const double &cos_phi_12)
{
    double alpha_1_2 = alpha_EdS(q_1,q_2,cos_phi_12);
    double beta_1_2 = beta_EdS(q_1,q_2,cos_phi_12);
    return 5/7.0*alpha_1_2 + 2/7.0*beta_1_2;
}

double G2(const double &q_1, const double &q_2, const double &cos_phi_12)
{
    double alpha_1_2 = alpha_EdS(q_1,q_2,cos_phi_12);
    double beta_1_2 = beta_EdS(q_1,q_2,cos_phi_12);
    return 3/7.0*alpha_1_2 + 4/7.0*beta_1_2;
}

double F2(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec)
{
    double cos_phi_1_2 = cos_phi_AB(q_1_vec, q_2_vec);

    double alpha_1_2 = alpha_EdS(q_1_vec.at(3),q_2_vec.at(3),cos_phi_1_2);
    double beta_1_2 = beta_EdS(q_1_vec.at(3),q_2_vec.at(3),cos_phi_1_2);
    return 5/7.0*alpha_1_2 + 2/7.0*beta_1_2;
}

double G2(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec)
{
    double cos_phi_1_2 = cos_phi_AB(q_1_vec, q_2_vec);

    double alpha_1_2 = alpha_EdS(q_1_vec.at(3),q_2_vec.at(3),cos_phi_1_2);
    double beta_1_2 = beta_EdS(q_1_vec.at(3),q_2_vec.at(3),cos_phi_1_2);
    return 3/7.0*alpha_1_2 + 4/7.0*beta_1_2;
}

double F3(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec)
{
    std::vector<double> q_1p2_vec = q_ApB_vec(q_1_vec, q_2_vec);
    std::vector<double> q_2p3_vec = q_ApB_vec(q_2_vec, q_3_vec);

    double cos_phi_1_2p3 = cos_phi_AB(q_1_vec, q_2p3_vec);
    double alpha_1_2p3 = alpha_EdS(q_1_vec.at(3),q_2p3_vec.at(3),cos_phi_1_2p3);
    double beta_1_2p3 = beta_EdS(q_1_vec.at(3),q_2p3_vec.at(3),cos_phi_1_2p3);

    double cos_phi_1p2_3 = cos_phi_AB(q_1p2_vec, q_3_vec);
    double alpha_1p2_3 = alpha_EdS(q_1p2_vec.at(3),q_3_vec.at(3),cos_phi_1p2_3);
    double beta_1p2_3 = beta_EdS(q_1p2_vec.at(3),q_3_vec.at(3),cos_phi_1p2_3);

    double F2_2_3 = F2(q_2_vec,q_3_vec);
    double G2_1_2 = G2(q_1_vec,q_2_vec);
    double G2_2_3 = G2(q_2_vec,q_3_vec);

    return 7/18.0*alpha_1_2p3*F2_2_3 + 2/18.0*beta_1_2p3*G2_2_3 + 7/18.0*alpha_1p2_3*G2_1_2 + 2/18.0*beta_1p2_3*G2_1_2;
}

double G3(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec)
{
    std::vector<double> q_2p3_vec = q_ApB_vec(q_2_vec, q_3_vec);
    std::vector<double> q_1p2_vec = q_ApB_vec(q_1_vec, q_2_vec);

    double cos_phi_1_2p3 = cos_phi_AB(q_1_vec, q_2p3_vec);
    double alpha_1_2p3 = alpha_EdS(q_1_vec.at(3),q_2p3_vec.at(3),cos_phi_1_2p3);
    double beta_1_2p3 = beta_EdS(q_1_vec.at(3),q_2p3_vec.at(3),cos_phi_1_2p3);

    double cos_phi_1p2_3 = cos_phi_AB(q_1p2_vec, q_3_vec);
    double alpha_1p2_3 = alpha_EdS(q_1p2_vec.at(3),q_3_vec.at(3),cos_phi_1p2_3);
    double beta_1p2_3 = beta_EdS(q_1p2_vec.at(3),q_3_vec.at(3),cos_phi_1p2_3);

    double F2_2_3 = F2(q_2_vec,q_3_vec);
    double G2_1_2 = G2(q_1_vec,q_2_vec);
    double G2_2_3 = G2(q_2_vec,q_3_vec);

    return 3/18.0*alpha_1_2p3*F2_2_3 + 6/18.0*beta_1_2p3*G2_2_3 + 3/18.0*alpha_1p2_3*G2_1_2 + 6/18.0*beta_1p2_3*G2_1_2;
}

double F4(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec, const std::vector<double> &q_4_vec)
{
    std::vector<double> q_2p3p4_vec = q_ApB_vec(q_ApB_vec(q_2_vec, q_3_vec),q_4_vec);
    std::vector<double> q_1p2_vec = q_ApB_vec(q_1_vec, q_2_vec);
    std::vector<double> q_3p4_vec = q_ApB_vec(q_3_vec, q_4_vec);
    std::vector<double> q_1p2p3_vec = q_ApB_vec(q_ApB_vec(q_1_vec, q_2_vec),q_3_vec);

    double cos_phi_1_2p3p4 = cos_phi_AB(q_1_vec, q_2p3p4_vec);
    double alpha_1_2p3p4 = alpha_EdS(q_1_vec.at(3),q_2p3p4_vec.at(3),cos_phi_1_2p3p4);
    double beta_1_2p3p4 = beta_EdS(q_1_vec.at(3),q_2p3p4_vec.at(3),cos_phi_1_2p3p4);

    double cos_phi_1p2_3p4 = cos_phi_AB(q_1p2_vec, q_3p4_vec);
    double alpha_1p2_3p4 = alpha_EdS(q_1p2_vec.at(3),q_3p4_vec.at(3),cos_phi_1p2_3p4);
    double beta_1p2_3p4 = beta_EdS(q_1p2_vec.at(3),q_3p4_vec.at(3),cos_phi_1p2_3p4);

    double cos_phi_1p2p3_4 = cos_phi_AB(q_1p2p3_vec, q_4_vec);
    double alpha_1p2p3_4 = alpha_EdS(q_1p2p3_vec.at(3),q_4_vec.at(3),cos_phi_1p2p3_4);
    double beta_1p2p3_4 = beta_EdS(q_1p2p3_vec.at(3),q_4_vec.at(3),cos_phi_1p2p3_4);

    double F3_2_3_4 = F3(q_2_vec,q_3_vec,q_4_vec);
    double G3_2_3_4 = G3(q_2_vec,q_3_vec,q_4_vec);

    double G2_1_2 = G2(q_1_vec,q_2_vec);
    double F2_3_4 = F2(q_3_vec,q_4_vec);
    double G2_3_4 = G2(q_3_vec,q_4_vec);

    double G3_1_2_3 = G3(q_1_vec,q_2_vec,q_3_vec);

    return 9/33.0*alpha_1_2p3p4*F3_2_3_4 + 2/33.0*beta_1_2p3p4*G3_2_3_4 + 9/33.0*alpha_1p2_3p4*G2_1_2*F2_3_4 + 2/33.0*beta_1p2_3p4*G2_1_2*G2_3_4
            + 9/33.0*alpha_1p2p3_4*G3_1_2_3 + 2/33.0*beta_1p2p3_4*G3_1_2_3;
}

double G4(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec, const std::vector<double> &q_4_vec)
{
    std::vector<double> q_2p3p4_vec = q_ApB_vec(q_ApB_vec(q_2_vec, q_3_vec),q_4_vec);
    std::vector<double> q_1p2_vec = q_ApB_vec(q_1_vec, q_2_vec);
    std::vector<double> q_3p4_vec = q_ApB_vec(q_3_vec, q_4_vec);
    std::vector<double> q_1p2p3_vec = q_ApB_vec(q_ApB_vec(q_1_vec, q_2_vec),q_3_vec);

    double cos_phi_1_2p3p4 = cos_phi_AB(q_1_vec, q_2p3p4_vec);
    double alpha_1_2p3p4 = alpha_EdS(q_1_vec.at(3),q_2p3p4_vec.at(3),cos_phi_1_2p3p4);
    double beta_1_2p3p4 = beta_EdS(q_1_vec.at(3),q_2p3p4_vec.at(3),cos_phi_1_2p3p4);

    double cos_phi_1p2_3p4 = cos_phi_AB(q_1p2_vec, q_3p4_vec);
    double alpha_1p2_3p4 = alpha_EdS(q_1p2_vec.at(3),q_3p4_vec.at(3),cos_phi_1p2_3p4);
    double beta_1p2_3p4 = beta_EdS(q_1p2_vec.at(3),q_3p4_vec.at(3),cos_phi_1p2_3p4);

    double cos_phi_1p2p3_4 = cos_phi_AB(q_1p2p3_vec, q_4_vec);
    double alpha_1p2p3_4 = alpha_EdS(q_1p2p3_vec.at(3),q_4_vec.at(3),cos_phi_1p2p3_4);
    double beta_1p2p3_4 = beta_EdS(q_1p2p3_vec.at(3),q_4_vec.at(3),cos_phi_1p2p3_4);

    double F3_2_3_4 = F3(q_2_vec,q_3_vec,q_4_vec);
    double G3_2_3_4 = G3(q_2_vec,q_3_vec,q_4_vec);

    double G2_1_2 = G2(q_1_vec,q_2_vec);
    double F2_3_4 = F2(q_3_vec,q_4_vec);
    double G2_3_4 = G2(q_3_vec,q_4_vec);

    double G3_1_2_3 = G3(q_1_vec,q_2_vec,q_3_vec);

    return 3/33.0*alpha_1_2p3p4*F3_2_3_4 + 8/33.0*beta_1_2p3p4*G3_2_3_4 + 3/33.0*alpha_1p2_3p4*G2_1_2*F2_3_4 + 8/33.0*beta_1p2_3p4*G2_1_2*G2_3_4
            + 3/33.0*alpha_1p2p3_4*G3_1_2_3 + 8/33.0*beta_1p2p3_4*G3_1_2_3;
}

// symmetrized kernels

double F2_sym(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec)
{
    if(q_1_vec.at(3) == 0 || q_2_vec.at(3) == 0)
        return 0;

    return 1/2.0*(F2(q_1_vec,q_2_vec)+F2(q_2_vec,q_1_vec));
}

double G2_sym(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec)
{
    if(q_1_vec.at(3) == 0 || q_2_vec.at(3) == 0)
        return 0;

    return 1/2.0*(G2(q_1_vec,q_2_vec)+G2(q_2_vec,q_1_vec));
}

double F3_sym(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec)
{
    if(q_1_vec.at(3) == 0 || q_2_vec.at(3) == 0 || q_3_vec.at(3) == 0)
        return 0;

    return 1/6.0*(F3(q_1_vec,q_2_vec,q_3_vec)+F3(q_2_vec,q_3_vec,q_1_vec)+F3(q_3_vec,q_1_vec,q_2_vec)+
                  F3(q_1_vec,q_3_vec,q_2_vec)+F3(q_2_vec,q_1_vec,q_3_vec)+F3(q_3_vec,q_2_vec,q_1_vec));
}

double F3_sym(const double &q_1, const double &q_2, const double &q_3, const double &mu_12, const double &mu_13, const double &mu_23)
{
    if(q_1 == 0 || q_3 == 0 || q_3 == 0)
        return 0;

    std::vector<double> q_1_vec = {q_1, 0, 0, q_1};
    std::vector<double> q_2_vec = {q_2*mu_12, q_2*sqrt(1 - mu_12*mu_12), 0, q_2};
    double q_3_x = q_3*mu_13;
    double q_3_y = (q_2*q_3*mu_23 - q_3_x*q_2_vec.at(0)) / q_2_vec.at(1);
    double q_3_z = sqrt(q_3*q_3 - q_3_x*q_3_x - q_3_y*q_3_y);
    std::vector<double> q_3_vec = {q_3_x, q_3_y, q_3_z, q_3};

    return F3_sym(q_1_vec,q_2_vec,q_3_vec);
}

double G3_sym(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec)
{
    if(q_1_vec.at(3) == 0 || q_2_vec.at(3) == 0 || q_3_vec.at(3) == 0)
        return 0;

    return 1/6.0*(G3(q_1_vec,q_2_vec,q_3_vec)+G3(q_2_vec,q_3_vec,q_1_vec)+G3(q_3_vec,q_1_vec,q_2_vec)+
                  G3(q_1_vec,q_3_vec,q_2_vec)+G3(q_2_vec,q_1_vec,q_3_vec)+G3(q_3_vec,q_2_vec,q_1_vec));
}

double F4_sym(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec, const std::vector<double> &q_4_vec)
{
    if(q_1_vec.at(3) == 0 || q_2_vec.at(3) == 0 || q_3_vec.at(3) == 0 || q_4_vec.at(3) == 0)
        return 0;

    return 1/24.0*(F4(q_1_vec,q_2_vec,q_3_vec,q_4_vec)+F4(q_2_vec,q_3_vec,q_4_vec,q_1_vec)+F4(q_3_vec,q_4_vec,q_1_vec,q_2_vec)+F4(q_4_vec,q_1_vec,q_2_vec,q_3_vec)+
                   F4(q_1_vec,q_3_vec,q_4_vec,q_2_vec)+F4(q_2_vec,q_4_vec,q_1_vec,q_3_vec)+F4(q_3_vec,q_1_vec,q_2_vec,q_4_vec)+F4(q_4_vec,q_2_vec,q_3_vec,q_1_vec)+
                   F4(q_1_vec,q_4_vec,q_2_vec,q_3_vec)+F4(q_2_vec,q_1_vec,q_3_vec,q_4_vec)+F4(q_3_vec,q_2_vec,q_4_vec,q_1_vec)+F4(q_4_vec,q_3_vec,q_1_vec,q_2_vec)+
                   F4(q_1_vec,q_2_vec,q_4_vec,q_3_vec)+F4(q_2_vec,q_3_vec,q_1_vec,q_4_vec)+F4(q_3_vec,q_4_vec,q_2_vec,q_1_vec)+F4(q_4_vec,q_1_vec,q_3_vec,q_2_vec)+
                   F4(q_1_vec,q_3_vec,q_2_vec,q_4_vec)+F4(q_2_vec,q_4_vec,q_3_vec,q_1_vec)+F4(q_3_vec,q_1_vec,q_4_vec,q_2_vec)+F4(q_4_vec,q_2_vec,q_1_vec,q_3_vec)+
                   F4(q_1_vec,q_4_vec,q_3_vec,q_2_vec)+F4(q_2_vec,q_1_vec,q_4_vec,q_3_vec)+F4(q_3_vec,q_2_vec,q_1_vec,q_4_vec)+F4(q_4_vec,q_3_vec,q_2_vec,q_1_vec));
}

double F4_sym(const double &q_1, const double &q_2, const double &q_3, const double &q_4,
              const double &mu_12, const double &mu_13, const double &mu_14, const double &mu_23, const double &mu_24, const double &mu_34)
{
    if(q_1 == 0 || q_2 == 0 || q_3 == 0 || q_4 == 0)
        return 0;

    std::vector<double> q_1_vec = {q_1, 0, 0, q_1};
    std::vector<double> q_2_vec = {q_2*mu_12, q_2*sqrt(1 - mu_12*mu_12), 0, q_2};
    double q_3_x = q_3*mu_13;
    double q_3_y = (q_2*q_3*mu_23 - q_3_x*q_2_vec.at(0)) / q_2_vec.at(1);
    double q_3_z = sqrt(q_3*q_3 - q_3_x*q_3_x - q_3_y*q_3_y);
    std::vector<double> q_3_vec = {q_3_x, q_3_y, q_3_z, q_3};

    double q_4_x = q_4*mu_14;
    double q_4_y = (q_2*q_4*mu_24 - q_4_x*q_2_vec.at(0)) / q_2_vec.at(1);
    double q_4_z = (q_3*q_4*mu_34 - q_4_x*q_3_vec.at(0) - q_4_y*q_3_vec.at(1)) / q_3_vec.at(2);

    double q_4_computed = sqrt(q_4_x*q_4_x + q_4_y*q_4_y + q_4_z*q_4_z);
    if (q_4_computed != q_4)
    {
        std::cout << q_4_computed << " is not equal to the length provided i.e. " << q_4 << std::endl;
        //assert (q_4_computed == q_4);
    }

    std::vector<double> q_4_vec = {q_4_x, q_4_y, q_4_z, q_4_computed};

    std::cout << "\n" << q_1_vec.at(0) << " " << q_1_vec.at(1) << " " << q_1_vec.at(2) << " " << q_1_vec.at(3) << std::endl;
    std::cout << q_2_vec.at(0) << " " << q_2_vec.at(1) << " " << q_2_vec.at(2) << " " << q_2_vec.at(3) << std::endl;
    std::cout << q_3_vec.at(0) << " " << q_3_vec.at(1) << " " << q_3_vec.at(2) << " " << q_3_vec.at(3) << std::endl;
    std::cout << q_4_vec.at(0) << " " << q_4_vec.at(1) << " " << q_4_vec.at(2) << " " << q_4_vec.at(3) << std::endl;

    return F4_sym(q_1_vec,q_2_vec,q_3_vec,q_4_vec);
}

double G4_sym(const std::vector<double> &q_1_vec, const std::vector<double> &q_2_vec, const std::vector<double> &q_3_vec, const std::vector<double> &q_4_vec)
{
    if(q_1_vec.at(3) == 0 || q_2_vec.at(3) == 0 || q_3_vec.at(3) == 0 || q_4_vec.at(3) == 0)
        return 0;

    return 1/24.0*(G4(q_1_vec,q_2_vec,q_3_vec,q_4_vec)+G4(q_2_vec,q_3_vec,q_4_vec,q_1_vec)+G4(q_3_vec,q_4_vec,q_1_vec,q_2_vec)+G4(q_4_vec,q_1_vec,q_2_vec,q_3_vec)+
                   G4(q_1_vec,q_3_vec,q_4_vec,q_2_vec)+G4(q_2_vec,q_4_vec,q_1_vec,q_3_vec)+G4(q_3_vec,q_1_vec,q_2_vec,q_4_vec)+G4(q_4_vec,q_2_vec,q_3_vec,q_1_vec)+
                   G4(q_1_vec,q_4_vec,q_2_vec,q_3_vec)+G4(q_2_vec,q_1_vec,q_3_vec,q_4_vec)+G4(q_3_vec,q_2_vec,q_4_vec,q_1_vec)+G4(q_4_vec,q_3_vec,q_1_vec,q_2_vec)+
                   G4(q_1_vec,q_2_vec,q_4_vec,q_3_vec)+G4(q_2_vec,q_3_vec,q_1_vec,q_4_vec)+G4(q_3_vec,q_4_vec,q_2_vec,q_1_vec)+G4(q_4_vec,q_1_vec,q_3_vec,q_2_vec)+
                   G4(q_1_vec,q_3_vec,q_2_vec,q_4_vec)+G4(q_2_vec,q_4_vec,q_3_vec,q_1_vec)+G4(q_3_vec,q_1_vec,q_4_vec,q_2_vec)+G4(q_4_vec,q_2_vec,q_1_vec,q_3_vec)+
                   G4(q_1_vec,q_4_vec,q_3_vec,q_2_vec)+G4(q_2_vec,q_1_vec,q_4_vec,q_3_vec)+G4(q_3_vec,q_2_vec,q_1_vec,q_4_vec)+G4(q_4_vec,q_3_vec,q_2_vec,q_1_vec));
}

