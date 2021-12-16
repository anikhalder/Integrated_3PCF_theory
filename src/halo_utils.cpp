#include <halo_utils.hpp>
#include <math.h>
#include <nonlinear.h>
#include <integration_utils.h>
#include <assert.h>

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
    double b1_L = (nu*nu-1.)/delta_c; // Lagrangian b1 bias

    return 1. + b1_L; // Eulerian b1 bias
}

double halo_b2_PBS(const double &M, const double &z, ClassEngine *class_obj)
{
    double R = R_L(M, class_obj);
    double nu = delta_c / class_obj->get_sigma_R_z(R, z);
    double b1_L = (nu*nu-1.)/delta_c; // Lagrangian b1 bias

    double nu2 = nu*nu;
    double b2_L = nu2*(nu2-3.)/(delta_c*delta_c); // Lagrangian b2 bias

    return 8./21.*b1_L + b2_L; // Eulerian b2 bias (the 8/21 factor appears from spherical symmetric collapse)
}

double halo_b2_Lazeyras(const double &M, const double &z, ClassEngine *class_obj)
{
    double b1 = halo_b1_Tinker2010(M, z, class_obj);
    //double b1 = halo_b1_PBS(M, z, class_obj);

    return 0.412 - 2.143*b1 + 0.929*pow(b1,2) + 0.008*pow(b1,3);
}

double halo_bs2_coevolution(const double &M, const double &z, ClassEngine *class_obj)
{
    // e.g. see Table 1 Leicht++ 2020 and also Pandey++ 2020
    double b1 = halo_b1_Tinker2010(M, z, class_obj);
    return -4./7.*(b1-1); // Note that in Desjacques++ and many other literature one works with bK2 = 1/2*bs2 = -2/7(b1-1)
}

double halo_b1_Tinker2010(const double &M, const double &z, ClassEngine *class_obj)
{
    // linear (Tinker) bias of a halo of mass M [M_sun]
    // eqn (6) and Table 2 of Tinker++ (2010) https://iopscience.iop.org/article/10.1088/0004-637X/724/2/878/pdf

    double R = R_L(M, class_obj);
    double nu = delta_c / class_obj->get_sigma_R_z(R, z);

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

    double R = R_L(M, class_obj);
    double sigma_R_z = class_obj->get_sigma_R_z(R, z);
    double sigma_prime_R_z = class_obj->get_sigma_prime_R_z(R, z); // in [1/Mpc]

    double nu = delta_c / sigma_R_z;

    double f_PS = f_nu_PS(nu);
    double rho_0 = class_obj->get_rho_m_z(0) / _M_SUN_ * pow(_Mpc_over_m_,3); // in [M_sun/Mpc^3]

    return f_PS*(-1./3.)*(rho_0/pow(M,2))*R/sigma_R_z*sigma_prime_R_z; // in [1/Mpc^3/M_sun]
}

double g_sigma_Tinker2010(const double &sigma_R_z, const double &z)
{
    // eqn (8) and Table 4 of Tinker++ (2010) https://iopscience.iop.org/article/10.1088/0004-637X/724/2/878/pdf

    double nu = delta_c / sigma_R_z;

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

    double f_sigma = alpha*(1. + pow(beta*nu,-2*phi))*pow(nu,2*eta)*exp(-gamma*nu*nu*0.5);

    return nu*f_sigma;
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

    double R = R_L(M, class_obj);
    double sigma_R_z = class_obj->get_sigma_R_z(R, z);
    double sigma_prime_R_z = class_obj->get_sigma_prime_R_z(R, z); // in [1/Mpc]

    double f_sigma = f_sigma_Tinker2008(sigma_R_z, z);
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
