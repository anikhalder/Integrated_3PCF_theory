#include <halo_utils.hpp>
#include <math.h>
#include <nonlinear.h>
#include <integration_utils.h>
#include <assert.h>
#include <constants.h>

namespace
{
    const double delta_c = 1.686;
    const double Delta = 200;
}

double R_L(const double &M, ClassEngine *class_obj)
{
    // comoving Lagrangian radius [Mpc] of a halo of mass M [M_sun]
    // Note: as this is a comoving quantity it doesn't matter at which redshift z the halo is identified
    // physical radius of halo at z : R_halo(z) = (3M/(4 pi rho(z)))^(1/3)
    //                                          = (3M a(z)^3/(4 pi rho_0))^(1/3)
    //                                          = a(z) (3M/(4 pi rho_0))^(1/3)
    //                                         := a(z) R_L
    // ==> R_L = (3M / (4 pi rho_0) )^(1/3) in comoving [Mpc] e.g. see eqn (12.63) of Dodelson & Schmidt (2020)

    return pow(3.0*M*_M_SUN_/(4.0*M_PI*class_obj->get_rho_m_z(0)),1./3.0) / _Mpc_over_m_; // in comoving [Mpc]
}

double halo_b1_PBS(const double &M, const double &z, ClassEngine *class_obj)
{
    // see eqn (2) of Jing++ (1998) https://iopscience.iop.org/article/10.1086/311530/pdf
    double R = R_L(M, class_obj);
    double nu = delta_c / class_obj->get_sigma_R_z(R, z);
    //double nu = delta_c / class_obj->get_sigma_R_z_lin(R, z);
    double b1_L = (nu*nu-1.)/delta_c; // Lagrangian b1 bias

    return 1. + b1_L; // Eulerian b1 bias
}

double halo_b2_PBS(const double &M, const double &z, ClassEngine *class_obj)
{
    double R = R_L(M, class_obj);
    double nu = delta_c / class_obj->get_sigma_R_z(R, z);
    //double nu = delta_c / class_obj->get_sigma_R_z_lin(R, z);
    double b1_L = (nu*nu-1.)/delta_c; // Lagrangian b1 bias

    double nu2 = nu*nu;
    double b2_L = nu2*(nu2-3.)/(delta_c*delta_c); // Lagrangian b2 bias

    return 8./21.*b1_L + b2_L; // Eulerian b2 bias (the 8/21 factor appears from spherical symmetric collapse)
}

double halo_b2_Lazeyras(const double &M, const double &z, ClassEngine *class_obj)
{
    double b1 = halo_b1_Tinker2010(M, z, class_obj);
    //double b1 = halo_b1_PBS(M, z, class_obj);

    return 0.412 - 2.143*b1 + 0.929*b1*b1 + 0.008*b1*b1*b1;
}

double halo_bs2_coevolution(const double &M, const double &z, ClassEngine *class_obj)
{
    // e.g. see Table 1 Leicht++ 2020 and also Pandey++ 2020
    double b1 = halo_b1_Tinker2010(M, z, class_obj);
    return -4./7.*(b1-1); // Note that in Desjacques++ and many other literature one works with bK2 = 1/2*bs2 = -2/7(b1-1)
}

double halo_b1_nu_Tinker2010(const double &nu)
{
    // linear (Tinker) bias of a halo of peak height nu
    // eqn (6) and Table 2 of Tinker++ (2010) https://iopscience.iop.org/article/10.1088/0004-637X/724/2/878/pdf

    double y = log10(Delta);
    double xp = exp(-1.0*pow(4./y,4.));

    double A = 1. + 0.24*y*xp;
    double a = 0.44*y - 0.88;
    double B = 0.183;
    double b = 1.5;
    double C = 0.019 + 0.107*y + 0.19*xp;
    double c = 2.4;

    double nu_a = pow(nu,a);
    double nu_b = pow(nu,b);
    double nu_c = pow(nu,c);

    return 1. - A*nu_a/(nu_a+pow(delta_c,a)) + B*nu_b + C*nu_c;
}

double halo_b1_Tinker2010(const double &M, const double &z, ClassEngine *class_obj)
{
    // linear (Tinker) bias of a halo of mass M [M_sun]
    // eqn (6) and Table 2 of Tinker++ (2010) https://iopscience.iop.org/article/10.1088/0004-637X/724/2/878/pdf

    double R = R_L(M, class_obj);
    double nu = delta_c / class_obj->get_sigma_R_z(R, z);
    //double nu = delta_c / class_obj->get_sigma_R_z_lin(R, z);

    return halo_b1_nu_Tinker2010(nu);
}

double f_nu_PS(const double &nu)
{
    // eqn (12.73) of Dodelson & Schmidt (2020)
    return sqrt(2./M_PI)*nu*exp(-nu*nu*0.5);
}

double dn_dM_PS(const double &M, const double &z, ClassEngine *class_obj)
{
    // halo mass M is in units [M_sun]
    // eqn (12.73) of Dodelson & Schmidt (2020)

    // dn/dln M = (rho_0/M)*f_PS(nu)* |dln sigma(M,z) / dln M|    where nu = delta_c / sigma(M,z)
    // dn/dM = (rho_0/M)*f_PS(nu)* |dln sigma(M,z) / d M|
    //       = (rho_0/M)*f_PS(nu)* (dln sigma^-1 / d M)
    //       = (rho_0/M)*f_PS(nu)*(1/(3M))* (dln sigma^-1 / dln R)    where dM = 3M dln R
    //       = f_PS(nu)*(-1/3)*(rho_0/M^2)* (dln sigma / dln R)
    //       = f_PS(nu)*(-1/3)*(rho_0/M^2)* R/sigma (d sigma / d R)

    // dn/dln M = M dn/dM = f_PS(nu)*(-1/3)*(rho_0/M)* R/sigma (d sigma / d R)


    double R = R_L(M, class_obj);
    double sigma_R_z = class_obj->get_sigma_R_z(R, z);
    double sigma_prime_R_z = class_obj->get_sigma_prime_R_z(R, z); // in [1/Mpc]

    double nu = delta_c / sigma_R_z;

    double f_PS = f_nu_PS(nu);
    double rho_0 = class_obj->get_rho_m_z(0) / _M_SUN_ * pow(_Mpc_over_m_,3); // in [M_sun/Mpc^3]

    return f_PS*(-1./3.)*(rho_0/pow(M,2))*R/sigma_R_z*sigma_prime_R_z; // in [1/Mpc^3/M_sun]
}

double f_nu_Tinker2010(const double &nu, const double &z)
{
    // eqn (8) and Table 4 of Tinker++ (2010) https://iopscience.iop.org/article/10.1088/0004-637X/724/2/878/pdf

    assert(Delta == 200);

    // these parameters were calibrated at z=0
    double alpha_200 = 0.368;
    double beta_200 = 0.589;
    double gamma_200 = 0.864;
    double phi_200 = -0.729;
    double eta_200 = -0.243;

    double alpha = alpha_200;
    double beta = beta_200*pow(1.+z,0.20);
    double phi = phi_200*pow(1.+z,-0.08);
    double eta = eta_200*pow(1.+z,0.27);
    double gamma = gamma_200*pow(1.+z,-0.01);

    double f_nu = alpha*(1. + pow(beta*nu,-2*phi))*pow(nu,2*eta)*exp(-gamma*nu*nu*0.5);

    return f_nu;
}

double g_sigma_Tinker2010(const double &sigma_R_z, const double &z)
{
    // eqn (8) and Table 4 of Tinker++ (2010) https://iopscience.iop.org/article/10.1088/0004-637X/724/2/878/pdf

    double nu = delta_c / sigma_R_z;

    return nu*f_nu_Tinker2010(nu, z);
}

double f_sigma_Tinker2008(const double &sigma_R_z, const double &z)
{
    // eqn (3) and Table 2 of Tinker++ (2008) https://iopscience.iop.org/article/10.1086/591439/pdf

    assert(Delta == 200);

    double A0_200 = 0.186;
    double a0_200 = 1.47;
    double b0_200 = 2.57;
    double c0_200 = 1.19;

    double A = A0_200*pow(1.+z,-0.14);
    double a = a0_200*pow(1.+z,-0.06);
    double alpha = pow(10,-pow(0.75/(log10(Delta/75.0)),1.2));
    double b = b0_200*pow(1.+z,-alpha);
    double c = c0_200;

    double f_sigma = A*(pow(sigma_R_z/b,-a) + 1.)*exp(-c/pow(sigma_R_z,2));

    return f_sigma;
}

double dn_dM_Tinker2008(const double &M, const double &z, ClassEngine *class_obj)
{
    // halo mass M is in units [M_sun]
    // eqn (2) of Tinker++ (2008) https://iopscience.iop.org/article/10.1086/591439/pdf

    // dn/dM = f(sigma)*(rho_0/M)* (dln sigma^-1 / d M)
    //       = f(sigma)*(rho_0/M)*(1/(3M))* (dln sigma^-1 / dln R)    where dM = 3M dln R
    //       = f(sigma)*(-1/3)*(rho_0/M^2)* (dln sigma / dln R)
    //       = f(sigma)*(-1/3)*(rho_0/M^2)* R/sigma (d sigma / d R)

    // dn/dln M = M dn/dM = f(sigma)*(-1/3)*(rho_0/M)* R/sigma (d sigma / d R)

    double R = R_L(M, class_obj);
    double sigma_R_z = class_obj->get_sigma_R_z(R, z);
    double sigma_prime_R_z = class_obj->get_sigma_prime_R_z(R, z); // in [1/Mpc]

    //double f_sigma = f_sigma_Tinker2008(sigma_R_z, z);
    double f_sigma = g_sigma_Tinker2010(sigma_R_z, z);
    double rho_0 = class_obj->get_rho_m_z(0) / _M_SUN_ * pow(_Mpc_over_m_,3); // in [M_sun/Mpc^3]

    return f_sigma*(-1./3.)*(rho_0/pow(M,2))*R/sigma_R_z*sigma_prime_R_z; // in [1/Mpc^3/M_sun]
}

int N_h_z_M_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value)
{
    assert(ndim == 2);
    assert(fdim == 1);

    double z = k[0], M = k[1];

    params_N_h_z_M_integrand *p = static_cast<params_N_h_z_M_integrand *>(params);

    double chi_z = p->class_obj->get_chi_z(z);
    double H_inv = 1./p->class_obj->get_H_z(z);

    value[0] = 4.*M_PI*H_inv*pow(chi_z,2)*dn_dM_Tinker2008(M, z, p->class_obj);

    //M = exp(k[1]);
    //value[0] = 4.*M_PI*H_inv*pow(chi_z,2)*dn_dM_Tinker2008(M, z, p->class_obj)*M; // for dln M integrand

    return 0;
}

double N_h_hcubature(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max)
{
    // eqn (11) of https://iopscience.iop.org/article/10.3847/1538-4357/aa943d/pdf

    double result = 0;              // the result from the integration
    double error = 0;               // the estimated error from the integration

    params_N_h_z_M_integrand args = {class_obj};

    std::vector<double> lower_limits = { z_min, M_min };
    std::vector<double> upper_limits = { z_max, M_max };

    hcubature_integration(N_h_z_M_hcubature_integrand, static_cast<void *>(&args), lower_limits, upper_limits, 2, 0, result, error);

    return result;
}

// HOD functions

double N_gal_expected_Zacharegkas2020(const double &M, const double &z)
{
    // halo mass M is in units [M_sun]
    // this is using the Lens1-Source4 RedMagic best fit parameters from Table D1 of Zacharegkas 2020 (DES collaboration)

    double log_M_min = 12.13;
    double M_1 = pow(10,13.64);
    double sigma_logM = 0.50;
    double alpha = 2.06;
    double f_cen = 0.13;

    // eqn (1) and (2) of Zacharegkas 2020

    double N_cen_expected = f_cen/2.0*(1.+erf((log10(M) - log_M_min)/sigma_logM));
    double N_sat_expected = N_cen_expected*pow(M/M_1,alpha);

    //return N_cen_expected + N_sat_expected;
    //return N_cen_expected;
    return N_sat_expected;
}

double n_z_gal_M_qag_integrand(double M, void *params)
{
    params_HOD_M_integrand *p = static_cast<params_HOD_M_integrand *>(params);

    double N_galaxies_expected = N_gal_expected_Zacharegkas2020(M, p->z);

    return dn_dM_Tinker2008(M, p->z, p->class_obj)*N_galaxies_expected;
}

double n_z_gal_qag(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max)
{
    // eqn (8) of Zacharegkas 2020

    double result = 0, error = 0;

    params_HOD_M_integrand args = {z, class_obj};

    qag_1D_integration(&n_z_gal_M_qag_integrand, static_cast<void *>(&args), M_min, M_max,  calls_1e3, result, error);

    //std::cout << "Number density of galaxies [h^3/cMpc^3] for the given HOD parameters at redshift " << z << " and in the given mass bin --> " <<  result / pow(class_obj->get_h(),3) << " +/- " << error / pow(class_obj->get_h(),3) << std::endl;

    return result;
}

double N_gal_z_qag_integrand(double z, void *params)
{
    params_HOD_z_integrand *p = static_cast<params_HOD_z_integrand *>(params);

    double chi_z = p->class_obj->get_chi_z(z);
    double H_inv = 1./p->class_obj->get_H_z(z);

    return 4.*M_PI*H_inv*pow(chi_z,2)*n_z_gal_qag(z, p->class_obj, p->M_min, p->M_max);
}

double N_gal_qag(ClassEngine *class_obj, const double &z_min, const double &z_max, const double &M_min, const double &M_max)
{
    double result = 0, error = 0;

    params_HOD_z_integrand args = {class_obj, M_min, M_max};

    qag_1D_integration(&N_gal_z_qag_integrand, static_cast<void *>(&args), z_min, z_max,  calls_1e3, result, error);

    std::cout << "Number of galaxies for the given HOD parameters in the given redshift range and mass bin --> " <<  result << " +/- " << error << std::endl;

    return result;    
}

double b1_avg_z_gal_M_qag_integrand(double M, void *params)
{
    params_HOD_M_integrand *p = static_cast<params_HOD_M_integrand *>(params);

    double N_galaxies_expected = N_gal_expected_Zacharegkas2020(M, p->z);

    return dn_dM_Tinker2008(M, p->z, p->class_obj)*N_galaxies_expected*halo_b1_Tinker2010(M, p->z, p->class_obj); // mean halo bias which hosts the corresponding HOD galaxies
    //return dn_dM_Tinker2008(M, p->z, p->class_obj)*N_galaxies_expected*M; // mean halo mass which hosts the corresponding HOD galaxies
}

double b1_avg_z_gal_qag(const double &z, ClassEngine *class_obj, const double &M_min, const double &M_max)
{
    // eqn (7) of Zacharegkas 2020

    double result = 0, error = 0;

    params_HOD_M_integrand args = {z, class_obj};

    qag_1D_integration(&b1_avg_z_gal_M_qag_integrand, static_cast<void *>(&args), M_min, M_max, calls_1e3, result, error);

    double n_gal_z = n_z_gal_qag(z, class_obj, M_min, M_max);

    std::cout << "Average linear bias of galaxies for the given HOD parameters at redshift " << z << " and in the given mass bin --> " << result / n_gal_z << " +/- " << error / n_gal_z << std::endl;

    return result / n_gal_z;
}

int b1_avg_gal_z_M_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value)
{
    assert(ndim == 2);
    assert(fdim == 1);

    double z = k[0], M = k[1];

    params_b1_avg_gal_z_M_integrand *p = static_cast<params_b1_avg_gal_z_M_integrand *>(params);

    double N_galaxies_expected = N_gal_expected_Zacharegkas2020(M, z);
    double n_gal_z = n_z_gal_qag(z, p->class_obj, p->M_min, p->M_max);

    value[0] = p->n_of_z->interp(z)*dn_dM_Tinker2008(M, z, p->class_obj)*N_galaxies_expected*halo_b1_Tinker2010(M, z, p->class_obj) / n_gal_z;

    return 0;
}

double b1_avg_gal(ClassEngine *class_obj, Linear_interp_1D *n_of_z, const double &z_min, const double &z_max, const double &M_min, const double &M_max)
{
    double result = 0;              // the result from the integration
    double error = 0;               // the estimated error from the integration

    params_b1_avg_gal_z_M_integrand args = {class_obj, n_of_z, M_min, M_max};

    std::vector<double> lower_limits = { z_min, M_min };
    std::vector<double> upper_limits = { z_max, M_max };

    hcubature_integration(b1_avg_gal_z_M_integrand, static_cast<void *>(&args), lower_limits, upper_limits, 2, calls_1e3, result, error);

    std::cout << "Average galaxy bias for given n(z) --> " <<  result << " +/- " << error << std::endl;

    return result;
}

// bias and halo mass function integral constraints

double b1_integral_constraint_nu_qag_integrand(double nu, void *params)
{
    params_integral_constraints_integrand *p = static_cast<params_integral_constraints_integrand *>(params);

    return halo_b1_nu_Tinker2010(nu)*f_nu_Tinker2010(nu, p->z);
}

double b1_integral_constraint_qag(const double &z)
{
    double result = 0, error = 0;

    //parameters in integrand
    params_integral_constraints_integrand args = {z};

    qagiu_1D_integration(&b1_integral_constraint_nu_qag_integrand, static_cast<void *>(&args), 0, calls_1e5, result, error);

    std::cout << "Tinker 2010 Halo bias integral constraint at redshift " << z << " --> " <<  result << " +/- " << error << std::endl;

    return result;
}

double g_sigma_integral_constraint_sigma_qag_integrand(double sigma, void *params)
{
    params_integral_constraints_integrand *p = static_cast<params_integral_constraints_integrand *>(params);

    return g_sigma_Tinker2010(sigma, p->z) / sigma;
}

double g_sigma_integral_constraint_qag(const double &z)
{
    double result = 0, error = 0;

    //parameters in integrand
    params_integral_constraints_integrand args = {z};

    qagiu_1D_integration(&g_sigma_integral_constraint_sigma_qag_integrand, static_cast<void *>(&args), 0, calls_1e5, result, error);

    std::cout << "Tinker 2010 Halo mass function integral constraint ( \int g(sigma) / sigma dsigma ) at redshift " << z << " --> " <<  result << " +/- " << error << std::endl;

    return result;
}

double f_nu_integral_constraint_nu_qag_integrand(double nu, void *params)
{
    params_integral_constraints_integrand *p = static_cast<params_integral_constraints_integrand *>(params);

    return f_nu_Tinker2010(nu, p->z);
}

double f_nu_integral_constraint_qag(const double &z)
{
    double result = 0, error = 0;

    //parameters in integrand
    params_integral_constraints_integrand args = {z};

    qagiu_1D_integration(&f_nu_integral_constraint_nu_qag_integrand, static_cast<void *>(&args), 0, calls_1e5, result, error);

    std::cout << "Tinker 2010 Halo mass function integral constraint ( \int f(nu) dnu ) at redshift " << z << " --> " <<  result << " +/- " << error << std::endl;

    return result;
}