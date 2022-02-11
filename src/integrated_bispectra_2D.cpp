#include <integrated_bispectra_2D.h>
#include <bispectrum.hpp>
#include <cosmology_utils.h>
#include <integration_utils.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <constants.h>
#include <assert.h>

namespace
{
    const double integration_pre_factor = 1./(16.*M_PI*M_PI*M_PI*M_PI); // 1/pow(2*M_PI,4)
}

// ######################################################################################

// Integrated 3-point function area pre-factors

// Gaussian-quadrature

double Adelta_theta_qag_integrand(double theta, void *params)
{
    params_Adelta_theta_integrand *p = static_cast<params_Adelta_theta_integrand *>(params);

    return theta*W2D_TH_RS_unnormalised(theta, p->theta_T);
}

double Adelta_qag(double theta_T)
{
    double result = 0;              // the result from the integration
    double error = 0;               // the estimated error from the integration

    //parameters in integrand
    params_Adelta_theta_integrand args = {theta_T};

    qag_1D_integration(&Adelta_theta_qag_integrand, static_cast<void *>(&args), 0, 4.0*theta_T, calls_1e5, result, error);

    return 2*M_PI*result;
}

double A2pt_phi_qag_integrand(double phi, void *params)
{
    params_A2pt_phi_integrand *p = static_cast<params_A2pt_phi_integrand *>(params);

    double length = l_ApB(p->theta, phi, p->alpha, p->phi_alpha);

    return W2D_TH_RS_unnormalised(length, p->theta_T);
}

double A2pt_theta_qag_integrand(double theta, void *params)
{
    params_A2pt_theta_phi_integrand *p = static_cast<params_A2pt_theta_phi_integrand *>(params);

    //parameters in integrand
    params_A2pt_phi_integrand args = {theta, p->phi_alpha, p->alpha, p->theta_T};

    double result = 0;              // the result from the integration
    double error = 0;               // the estimated error from the integration

    qag_1D_integration(&A2pt_phi_qag_integrand, static_cast<void *>(&args), 0, 2*M_PI, calls_1e5, result, error);

    return theta*W2D_TH_RS_unnormalised(theta, p->theta_T)*result;
}

double A2pt_qag(double alpha, double theta_T)
{
    double phi_alpha = 0;           // assuming angular-independency we can just fix the phi_alpha to any angle we want
    double result = 0;              // the result from the integration
    double error = 0;               // the estimated error from the integration

    //parameters in integrand
    params_A2pt_theta_phi_integrand args = {phi_alpha, alpha, theta_T};

    qag_1D_integration(&A2pt_theta_qag_integrand, static_cast<void *>(&args), 0, 4.0*theta_T, calls_1e5, result, error);

    return result;
}

// Monte-Carlo

double A2pt_theta_phi_mc_integrand(double *k, size_t dim, void *params)
{
    (void)(dim); // avoid unused parameter warnings

    struct params_A2pt_theta_phi_integrand *p = static_cast<params_A2pt_theta_phi_integrand *>(params);

    double theta = k[0];
    double phi = k[1];

    double length = l_ApB(theta, phi, p->alpha, p->phi_alpha);

    return theta*W2D_TH_RS_unnormalised(theta, p->theta_T)*W2D_TH_RS_unnormalised(length, p->theta_T);
}

double A2pt_mc(double alpha, double theta_T, const gsl_rng_type *T, const std::string &mc_integration_type)
{
    double phi_alpha = 0;           // assuming angular-independency we can just fix the phi_alpha to any angle we want
    double result = 0;              // the result from the integration
    double error = 0;               // the estimated error from the integration

    //parameters in integrand
    params_A2pt_theta_phi_integrand args = {phi_alpha, alpha, theta_T};

    gsl_monte_function G = { &A2pt_theta_phi_mc_integrand, 1, static_cast<void *>(&args)};

    std::vector<double> lower_limits = { 0, 0};
    std::vector<double> upper_limits = { 4.0*theta_T, 2*M_PI};

    if (mc_integration_type == "plain" )
        monte_carlo_plain_integration(&G, lower_limits, upper_limits, 2, calls_1e5, T, result, error);
    else if (mc_integration_type == "miser" )
        monte_carlo_miser_integration(&G, lower_limits, upper_limits, 2, calls_1e5, T, result, error);
    else if (mc_integration_type == "vegas" )
        monte_carlo_vegas_integration(&G, lower_limits, upper_limits, 2, calls_1e5, T, result, error);

    return result;
}

// h-cubature

int A2pt_theta_phi_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value)
{
    assert(ndim == 2);
    assert(fdim == 1);

    struct params_A2pt_theta_phi_integrand *p = static_cast<params_A2pt_theta_phi_integrand *>(params);

    double theta = k[0];
    double phi = k[1];

    double length = l_ApB(theta, phi, p->alpha, p->phi_alpha);

    value[0] = theta*W2D_TH_RS_unnormalised(theta, p->theta_T)*W2D_TH_RS_unnormalised(length, p->theta_T);

    return 0;
}

double A2pt_hcubature(double alpha, double theta_T)
{
    double phi_alpha = 0;           // assuming angular-independency we can just fix the phi_alpha to any angle we want
    double result = 0;              // the result from the integration
    double error = 0;               // the estimated error from the integration

    //parameters in integrand
    params_A2pt_theta_phi_integrand args = {phi_alpha, alpha, theta_T};

    std::vector<double> lower_limits = { 0, 0};
    std::vector<double> upper_limits = { 4.0*theta_T, 2*M_PI};

    hcubature_integration(A2pt_theta_phi_hcubature_integrand, static_cast<void *>(&args), lower_limits, upper_limits, 2, 0, result, error);

    return result;
}

// ##############################

// testing with phi_alpha angle average integration; this is redundant and time consuming ---> does not affect the results obtained only from A2pt_qag()

double A2pt_phi_alpha_qag_integrand(double phi_alpha, void *params)
{
    params_A2pt_phi_alpha_theta_phi_integrand *p = static_cast<params_A2pt_phi_alpha_theta_phi_integrand *>(params);

    //parameters in integrand
    params_A2pt_theta_phi_integrand args = {phi_alpha, p->alpha, p->theta_T};

    double result = 0;              // the result from the integration
    double error = 0;               // the estimated error from the integration

    //qag_1D_integration(&A2pt_theta_qag_integrand, static_cast<void *>(&args), 0, 4.0*p->theta_T, calls_1e5, result, error);
    qag_1D_integration_abs_rel(&A2pt_theta_qag_integrand, static_cast<void *>(&args), 0, 4.0*p->theta_T, calls_1e5, result, error);

    //return cos(4*phi_alpha)*result;
    return result;
}

double A2pt_angle_averaged_qag(double alpha, double theta_T)
{
    //parameters in integrand
    params_A2pt_phi_alpha_theta_phi_integrand args = {alpha, theta_T};

    double result = 0;              // the result from the integration
    double error = 0;               // the estimated error from the integration

    //qag_1D_integration(&A2pt_phi_alpha_qag_integrand, static_cast<void *>(&args), 0, 2*M_PI, calls_1e5, result, error);
    qag_1D_integration_abs_rel(&A2pt_phi_alpha_qag_integrand, static_cast<void *>(&args), 0, 2*M_PI, calls_1e5, result, error);

    return result/(2*M_PI);
}

// ######################################################################################
// ######################################################################################
// ######################################################################################

double W_products(const std::string &key, const double &l_1, const double &phi_1, const double &l_2, const double &phi_2, const double &l, const double &phi_l,
                  const struct_iB2D_W_FS &info_iB2D_W_FS)
{
    if (key == "X")
    {
        double W_1 = info_iB2D_W_FS.W_meandelta_FS(l_1,info_iB2D_W_FS.theta_T_meandelta);
        double W_2pl = info_iB2D_W_FS.W_2pt_FS(l_ApB(l_2,phi_2, l,phi_l), info_iB2D_W_FS.theta_T_2pt);
        double W_m1m2ml = info_iB2D_W_FS.W_2pt_FS(l_ApBpC(l_1,M_PI+phi_1, l_2,M_PI+phi_2, l,M_PI+phi_l), info_iB2D_W_FS.theta_T_2pt);

        return W_1*W_2pl*W_m1m2ml;
    }

    else if (key == "X_v2")
    {
        // this is for another version of the argument form of the window functions; TODO: Need to investigate this further
        double W_m1m2 = info_iB2D_W_FS.W_meandelta_FS(l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2), info_iB2D_W_FS.theta_T_meandelta);
        double W_1pl = info_iB2D_W_FS.W_2pt_FS(l_ApB(l_1,phi_1, l,phi_l), info_iB2D_W_FS.theta_T_2pt);
        double W_2ml = info_iB2D_W_FS.W_2pt_FS(l_ApB(l_2,phi_2, l,M_PI+phi_l), info_iB2D_W_FS.theta_T_2pt);
        
        return W_m1m2*W_1pl*W_2ml;
    }

    else if (key == "Y")
    {
        double W_1 = info_iB2D_W_FS.W_meandelta_FS(l_1, info_iB2D_W_FS.theta_T_meandelta);
        double W_m1m2pl = info_iB2D_W_FS.W_2pt_FS(l_ApBpC(l_1,M_PI+phi_1, l_2,M_PI+phi_2, l,phi_l), info_iB2D_W_FS.theta_T_2pt);
        double W_2ml = info_iB2D_W_FS.W_2pt_FS(l_ApB(l_2,phi_2, l,M_PI+phi_l), info_iB2D_W_FS.theta_T_2pt);

        return W_1*W_m1m2pl*W_2ml;
    }

    else if (key == "Z")
    {
        double W_m1m2 = info_iB2D_W_FS.W_meandelta_FS(l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2), info_iB2D_W_FS.theta_T_meandelta);
        double W_2pl = info_iB2D_W_FS.W_2pt_FS(l_ApB(l_2,phi_2, l,phi_l), info_iB2D_W_FS.theta_T_2pt);
        double W_1ml = info_iB2D_W_FS.W_2pt_FS(l_ApB(l_1,phi_1, l,M_PI+phi_l), info_iB2D_W_FS.theta_T_2pt);

        return W_m1m2*W_2pl*W_1ml;
    }

    return 0.0;
}

double evaluate_iB2D_phi_1_phi_2_integrand(const std::string &key, const double &l, const double &phi_l, const struct_iB2D_W_FS &info_iB2D_W_FS,
                                         const double &l_1, const double &l_2, const double &phi_1, const double &phi_2)
{
    if (key == "X_lF2")
    {
        double lF2_1_2 = kF2_EdS_angular(l_1,l_2,cos(phi_1-phi_2));

        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return lF2_1_2*X_W;
    }

    else if (key == "Y_lF2")
    {
        double lF2_1_2 = kF2_EdS_angular(l_1,l_2,cos(phi_1-phi_2));

        double Y_W = W_products("Y", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return lF2_1_2*Y_W;
    }

    else if (key == "Z_lF2")
    {
        double lF2_1_2 = kF2_EdS_angular(l_1,l_2,cos(phi_1-phi_2));

        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return lF2_1_2*Z_W;
    }

    else if (key == "lF2")
    {
        double lF2_1_2 = kF2_EdS_angular(l_1,l_2,cos(phi_1-phi_2));

        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        double Y_W = W_products("Y", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return lF2_1_2*(X_W + Y_W + Z_W);
    }

    else if (key == "X_W" || key == "P_1" || key == "P_2" || key == "A")
    {
        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return X_W;
    }

    else if (key == "Y_W")
    {
        double Y_W = W_products("Y", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Y_W;
    }

    else if (key == "Z_W" || key == "P_3")
    {
        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Z_W;
    }

    else if (key == "X_S2")
    {
        double S2_1_2 = S2_angular(cos(phi_1-phi_2));

        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return S2_1_2*X_W;
    }

    else if (key == "Y_S2")
    {
        double S2_1_2 = S2_angular(cos(phi_1-phi_2));

        double Y_W = W_products("Y", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return S2_1_2*Y_W;
    }

    else if (key == "Z_S2")
    {
        double S2_1_2 = S2_angular(cos(phi_1-phi_2));

        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return S2_1_2*Z_W;
    }

    return 0.0;
}

// h-cubature

int iB2D_phi_1_phi_2_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value)
{
    assert(ndim == 2);
    assert(fdim == 1);

    double phi_1 = k[0],  phi_2 = k[1];

    params_iB2D_phi_1_phi_2_integrand *p = static_cast<params_iB2D_phi_1_phi_2_integrand *>(params);

    value[0] = evaluate_iB2D_phi_1_phi_2_integrand(p->key, p->l, p->phi_l, p->info_iB2D_W_FS, p->l_1, p->l_2, phi_1, phi_2);

    return 0;
}

void iB2D_phi_1_phi_2_hcubature(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, const double &l_1, const double &l_2,
                              double &result, double &error, size_t max_evals)
{
    double phi_l = 0;           // assuming angular-independency we can just fix the phi_l to any angle we want

    params_iB2D_phi_1_phi_2_integrand args = {key, l, phi_l, info_iB2D_W_FS, l_1, l_2};

    std::vector<double> lower_limits = { 0, 0};
    std::vector<double> upper_limits = { 2*M_PI, 2*M_PI};

    hcubature_integration(iB2D_phi_1_phi_2_hcubature_integrand, static_cast<void *>(&args), lower_limits, upper_limits, 2, max_evals, result, error);
}

// ######################################################################################

double evaluate_iB2D_l_1_l_2_phi_1_phi_2_integrand(const std::string &key, const double &l, const double &phi_l, const double &z, const struct_iB2D_W_FS &info_iB2D_W_FS, 
                                                   ClassEngine *class_obj, bool use_pk_nl, const double &l_1, const double &l_2, const double &phi_1, const double &phi_2)
{
    double chi_inv = 1/class_obj->get_chi_z(z);
    double X_windows = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS); // default settings

    // for lll integared bispectra
    if (key == "B_nothing" || key == "A")
        return l_1*l_2*X_windows;

    else if (key == "B_P")
        return B_P(l_1*chi_inv, l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, z, class_obj, use_pk_nl)*l_1*l_2*X_windows;

    else if (key == "B_PP")
        return B_PP(l_1*chi_inv, l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, z, class_obj, use_pk_nl)*l_1*l_2*X_windows;

    else if (key == "B_S2PP")
        return B_S2PP(l_1*chi_inv, l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, z, class_obj, use_pk_nl)*l_1*l_2*X_windows;

    else if (key == "B_hhh_eps_eps_eps")
        return l_1*l_2*X_windows;

    else if (key == "B_hhh_delta_eps_eps")
        return B_P(l_1*chi_inv, l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, z, class_obj, use_pk_nl)*l_1*l_2*X_windows;

    else if (key == "B_P2P3")
        return B_PaPb(l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, z, class_obj, use_pk_nl)*l_1*l_2*X_windows;

    else if (key == "B_S2P2P3")
        return B_S2PaPb(l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, l_1*chi_inv, z, class_obj, use_pk_nl)*l_1*l_2*X_windows;

    // for sss integrated bispectra
    if (key == "B")
        X_windows *= 1;

    // default settings
    else if (key == "B_xip_cos")
        X_windows *= cos(2*(phi_2 - phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)));

    else if (key == "B_xip_sin")
        X_windows *= sin(2*(phi_2 - phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)));

    else if (key == "B_xim_cos")
        X_windows *= cos(2*(phi_2 + phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)));

    else if (key == "B_xim_sin")
        X_windows *= sin(2*(phi_2 + phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)));

    return B(l_1*chi_inv, l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, z, class_obj, use_pk_nl)*l_1*l_2*X_windows;
}

// Monte-Carlo

double iB2D_l_1_l_2_phi_1_phi_2_mc_integrand(double *k, size_t dim, void *params)
{
    (void)(dim); // avoid unused parameter warnings

    double l_1 = k[0], l_2 = k[1], phi_1 = k[2], phi_2 = k[3];

    params_iB2D_l_1_l_2_phi_1_phi_2_integrand *p = static_cast<params_iB2D_l_1_l_2_phi_1_phi_2_integrand *>(params);

    return evaluate_iB2D_l_1_l_2_phi_1_phi_2_integrand(p->key, p->l, p->phi_l, p->z, p->info_iB2D_W_FS, p->class_obj, p->use_pk_nl, l_1, l_2, phi_1, phi_2);

}

void iB2D_l_1_l_2_phi_1_phi_2_mc(const std::string &key, const double &l, const double &z, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                   std::vector<double> lower_limits, std::vector<double> upper_limits,
                                   const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls)
{
    double phi_l = 0.;           // assuming angular-independency we can just fix the phi_l to any angle we want
//    double phi_l = 0.625*M_PI;           // assuming angular-independency we can just fix the phi_l to any angle we want
//    double phi_l = 1.378*M_PI;           // assuming angular-independency we can just fix the phi_l to any angle we want

    params_iB2D_l_1_l_2_phi_1_phi_2_integrand args = {key, l, phi_l, z, info_iB2D_W_FS, class_obj, use_pk_nl};

    gsl_monte_function G = { &iB2D_l_1_l_2_phi_1_phi_2_mc_integrand, 4, static_cast<void *>(&args)};

    if (mc_integration_type == "plain" )
        monte_carlo_plain_integration(&G, lower_limits, upper_limits, 4, calls, T, result, error);
    else if (mc_integration_type == "miser" )
        monte_carlo_miser_integration(&G, lower_limits, upper_limits, 4, calls, T, result, error);
    else if (mc_integration_type == "vegas" )
        monte_carlo_vegas_integration(&G, lower_limits, upper_limits, 4, calls, T, result, error);
}

void iB2D_mc_4_dim(const std::string &key, const double &l, const double &z, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                    std::vector<double> lower_limits, std::vector<double> upper_limits,
                    const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls)
{
    result = 0.0;
    error = 0.0;

    if (key == "A")
        // passing l = 0 and z = 0explicitly as there is no l or z dependence for the integrals (depends only on l_1 and l_2)
        iB2D_l_1_l_2_phi_1_phi_2_mc("A", 0, 0, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_xip_cos")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B_xip_cos", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_xip_sin")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B_xip_sin", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_xim_cos")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B_xim_cos", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_xim_sin")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B_xim_sin", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    // components for tracer integrated bispectrum

    else if (key == "B_nothing")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B_nothing", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_P")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B_P", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_PP")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B_PP", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_S2PP")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B_S2PP", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_hhh_delta_eps_eps")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B_hhh_delta_eps_eps", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_hhh_eps_eps_eps")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B_hhh_eps_eps_eps", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_P2P3")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B_P2P3", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_S2P2P3")
        iB2D_l_1_l_2_phi_1_phi_2_mc("B_S2P2P3", l, z, info_iB2D_W_FS, class_obj, use_pk_nl, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    double patch_sqradians = spherical_cap_radius_2_sqradians(info_iB2D_W_FS.theta_T_2pt);
    result = integration_pre_factor * patch_sqradians * patch_sqradians * result;
    error = integration_pre_factor * patch_sqradians * patch_sqradians * error;
}

double iB2D_trapz_z(Linear_interp_1D *iB2D_z_interp, double &z_lower_limit, double &z_upper_limit, ClassEngine *class_obj, 
                    projection_kernel *q1, projection_kernel *q2, projection_kernel *q3)
{
    // uniform trapezoidal integration
    double delta_z = (z_upper_limit-z_lower_limit)/num_trapz_steps;

    double integral = 0.;

    for (int i = 0; i <= num_trapz_steps; i++)
    {
        double z = z_lower_limit + i*delta_z;

        double q_1 = q1->evaluate(z);
        double q_2 = q2->evaluate(z);
        double q_3 = q3->evaluate(z);

        if (q_1 == 0 || q_2 == 0 || q_3 == 0)
            continue;

        double chi_inv = 1/class_obj->get_chi_z(z);
        double Hz_inv = 1/class_obj->get_H_z(z);

        double val = Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*iB2D_z_interp->interp(z);

        if (i==0 || i==num_trapz_steps)
            integral += val;       
        else
            integral += 2.*val;
    }

     return integral*delta_z*0.5;
}

// ######################################################################################

double evaluate_iB2D_z_l_1_l_2_phi_1_phi_2_integrand(const std::string &key, const double &l, const double &phi_l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj,
                                                     bool use_pk_nl, const double &z, const double &l_1, const double &l_2, const double &phi_1, const double &phi_2,
                                                     projection_kernel *q1, projection_kernel *q2, projection_kernel *q3)
{
    double q_1 = q1->evaluate(z);
    double q_2 = q2->evaluate(z);
    double q_3 = q3->evaluate(z);

    if (q_1 == 0 || q_2 == 0 || q_3 == 0)
        return 0.0;

    double chi_inv = 1/class_obj->get_chi_z(z);
    double Hz_inv = 1/class_obj->get_H_z(z);
    double X_windows = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS); // default settings
    //double X_windows = W_products("X_v2", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

    // outdated code (To be DELETED)
    /*
    if (key == "X_lF2")
    {
        double lF2_1_2 = kF2_EdS_angular(l_1,l_2,cos(phi_1-phi_2));

        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                lF2_1_2*X_W;
    }

    else if (key == "Y_lF2")
    {
        double lF2_1_2 = kF2_EdS_angular(l_1,l_2,cos(phi_1-phi_2));

        double Y_W = W_products("Y", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                lF2_1_2*Y_W;
    }

    else if (key == "Z_lF2")
    {
        double lF2_1_2 = kF2_EdS_angular(l_1,l_2,cos(phi_1-phi_2));

        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                lF2_1_2*Z_W;
    }

    else if (key == "B_lF2")
    {
        double lF2_1_2 = kF2_EdS_angular(l_1,l_2,cos(phi_1-phi_2));

        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        double Y_W = W_products("Y", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return 2*Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                lF2_1_2*(X_W + Y_W + Z_W);
    }

    else if (key == "X_W")
    {
        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                l_1*class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                l_2*class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                X_W;
    }

    else if (key == "Y_W")
    {
        double Y_W = W_products("Y", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                l_1*class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                l_2*class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                Y_W;
    }

    else if (key == "Z_W")
    {
        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                l_1*class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                l_2*class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                Z_W;
    }

    else if (key == "X_S2")
    {
        double S2_1_2 = S2_angular(cos(phi_1-phi_2));

        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                l_1*class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                l_2*class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                S2_1_2*X_W;
    }

    else if (key == "Y_S2")
    {
        double S2_1_2 = S2_angular(cos(phi_1-phi_2));

        double Y_W = W_products("Y", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                l_1*class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                l_2*class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                S2_1_2*Y_W;
    }

    else if (key == "Z_S2")
    {
        double S2_1_2 = S2_angular(cos(phi_1-phi_2));

        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                l_1*class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                l_2*class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                S2_1_2*Z_W;
    }

    else if (key == "P_1")
    {
        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                l_1*class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                l_2*
                X_W;
    }

    else if (key == "P_2")
    {
        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                l_1*
                l_2*class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                X_W;
    }

    else if (key == "P_3")
    {
        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS);

        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                l_1*class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                l_2*
                Z_W;
    }

    else if (key == "B_lF2_xip_cos")
    {
        double lF2_1_2 = kF2_EdS_angular(l_1,l_2,cos(phi_1-phi_2));

        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS)*cos(2*(phi_2 - phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)));

        double Y_W = W_products("Y", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS)*cos(2*(phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2) - phi_2));

        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS)*cos(2*(phi_2 - phi_1));

        return 2*Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                lF2_1_2*(X_W + Y_W + Z_W);
    }

    else if (key == "B_lF2_xip_sin")
    {
        double lF2_1_2 = kF2_EdS_angular(l_1,l_2,cos(phi_1-phi_2));

        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS)*sin(2*(phi_2 - phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)));

        double Y_W = W_products("Y", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS)*sin(2*(phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2) - phi_2));

        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS)*sin(2*(phi_2 - phi_1));

        return 2*Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                lF2_1_2*(X_W + Y_W + Z_W);
    }

    else if (key == "B_lF2_xim_cos")
    {
        double lF2_1_2 = kF2_EdS_angular(l_1,l_2,cos(phi_1-phi_2));

        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS)*cos(2*(phi_2 + phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)));

        double Y_W = W_products("Y", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS)*cos(2*(phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2) + phi_2));

        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS)*cos(2*(phi_2 + phi_1));

        return 2*Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                lF2_1_2*(X_W + Y_W + Z_W);
    }

    else if (key == "B_lF2_xim_sin")
    {
        double lF2_1_2 = kF2_EdS_angular(l_1,l_2,cos(phi_1-phi_2));

        double X_W = W_products("X", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS)*sin(2*(phi_2 + phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)));

        double Y_W = W_products("Y", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS)*sin(2*(phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2) + phi_2));

        double Z_W = W_products("Z", l_1, phi_1, l_2, phi_2, l, phi_l, info_iB2D_W_FS)*sin(2*(phi_2 + phi_1));

        return 2*Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                class_obj->pk(l_1*chi_inv, z, use_pk_nl)*
                class_obj->pk(l_2*chi_inv, z, use_pk_nl)*
                lF2_1_2*(X_W + Y_W + Z_W);
    }
    */

    // for lll integared bispectra
    if (key == "B_nothing" || key == "A")
        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*l_1*l_2*X_windows;

    else if (key == "B_P")
        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                B_P(l_1*chi_inv, l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, z, class_obj, use_pk_nl)*
                l_1*l_2*X_windows;

    else if (key == "B_PP")
        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                B_PP(l_1*chi_inv, l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, z, class_obj, use_pk_nl)*
                l_1*l_2*X_windows;

    else if (key == "B_S2PP")
        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                B_S2PP(l_1*chi_inv, l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, z, class_obj, use_pk_nl)*
                l_1*l_2*X_windows;

    else if (key == "B_hhh_eps_eps_eps")
    {
        double n_h_z = dynamic_cast<projection_kernel_q_h*>(q3)->get_n_h_z(z);
        double _1_over_n_h_z_2 = 1./(n_h_z*n_h_z);
        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*_1_over_n_h_z_2*l_1*l_2*X_windows;
    }

    else if (key == "B_hhh_delta_eps_eps")
    {
        double _1_over_n_h_z = 1./dynamic_cast<projection_kernel_q_h*>(q3)->get_n_h_z(z);
        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                B_P(l_1*chi_inv, l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, z, class_obj, use_pk_nl)*_1_over_n_h_z*
                l_1*l_2*X_windows;
    }

    else if (key == "B_P2P3")
        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                B_PaPb(l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, z, class_obj, use_pk_nl)*
                l_1*l_2*X_windows;

    else if (key == "B_S2P2P3")
        return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                B_S2PaPb(l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, l_1*chi_inv, z, class_obj, use_pk_nl)*
                l_1*l_2*X_windows;

    // for sss integrated bispectra
    if (key == "B")
        X_windows *= 1;

    // default settings
    else if (key == "B_xip_cos")
        X_windows *= cos(2*(phi_2 - phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)));

    else if (key == "B_xip_sin")
        X_windows *= sin(2*(phi_2 - phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)));

    else if (key == "B_xim_cos")
        X_windows *= cos(2*(phi_2 + phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)));

    else if (key == "B_xim_sin")
        X_windows *= sin(2*(phi_2 + phi_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)));

    // this is for another version of the argument form of the window functions; TODO: Need to investigate this further
    // else if (key == "B_xip_cos")
    //     X_windows *= cos(2*(phi_1 - phi_2));

    // else if (key == "B_xip_sin")
    //     X_windows *= sin(2*(phi_1 - phi_2));

    // else if (key == "B_xim_cos")
    //     X_windows *= cos(2*(phi_1 + phi_2));

    // else if (key == "B_xim_sin")
    //     X_windows *= sin(2*(phi_1 + phi_2));

    return Hz_inv*chi_inv*chi_inv*chi_inv*chi_inv*q_1*q_2*q_3*
                B(l_1*chi_inv, l_2*chi_inv, l_ApB(l_1,M_PI+phi_1, l_2,M_PI+phi_2)*chi_inv, z, class_obj, use_pk_nl)*
                l_1*l_2*X_windows;
}

// Monte-Carlo

double iB2D_z_l_1_l_2_phi_1_phi_2_mc_integrand(double *k, size_t dim, void *params)
{
    (void)(dim); // avoid unused parameter warnings

    double z = k[0], l_1 = k[1], l_2 = k[2], phi_1 = k[3], phi_2 = k[4];

    params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand *p = static_cast<params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand *>(params);

    return evaluate_iB2D_z_l_1_l_2_phi_1_phi_2_integrand(p->key, p->l, p->phi_l, p->info_iB2D_W_FS, p->class_obj, p->use_pk_nl, z, l_1, l_2, phi_1, phi_2,
                                                         p->q1, p->q2, p->q3);

}

void iB2D_z_l_1_l_2_phi_1_phi_2_mc(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                   projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                   const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls)
{
    double phi_l = 0.;           // assuming angular-independency we can just fix the phi_l to any angle we want
//    double phi_l = 0.625*M_PI;           // assuming angular-independency we can just fix the phi_l to any angle we want
//    double phi_l = 1.378*M_PI;           // assuming angular-independency we can just fix the phi_l to any angle we want

    params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand args = {key, l, phi_l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3};

    gsl_monte_function G = { &iB2D_z_l_1_l_2_phi_1_phi_2_mc_integrand, 5, static_cast<void *>(&args)};

    if (mc_integration_type == "plain" )
        monte_carlo_plain_integration(&G, lower_limits, upper_limits, 5, calls, T, result, error);
    else if (mc_integration_type == "miser" )
        monte_carlo_miser_integration(&G, lower_limits, upper_limits, 5, calls, T, result, error);
    else if (mc_integration_type == "vegas" )
        monte_carlo_vegas_integration(&G, lower_limits, upper_limits, 5, calls, T, result, error);
}

void iB2D_mc(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
           projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
           const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls)
{
    result = 0.0;
    error = 0.0;

    if (key == "B_lF2")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_lF2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "X_W")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("X_W", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "Y_W")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("Y_W", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "Z_W")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("Z_W", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "X_S2")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("X_S2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "Y_S2")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("Y_S2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "Z_S2")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("Z_S2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "P_1")
        // passing l = 0 explicitly as there is no l dependence for the integrals (depends only on l_1 and l_2)
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("P_1", 0, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "P_2")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("P_2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "P_3")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("P_3", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "A")
        // passing l = 0 explicitly as there is no l dependence for the integrals (depends only on l_1 and l_2)
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("A", 0, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    if (key == "B")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_xip_cos")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_xip_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_xip_sin")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_xip_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_xim_cos")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_xim_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_xim_sin")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_xim_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    // components for tracer integrated bispectrum

    else if (key == "B_nothing")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_nothing", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_P")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_P", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_PP")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_PP", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_S2PP")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_S2PP", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_hhh_delta_eps_eps")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_hhh_delta_eps_eps", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_hhh_eps_eps_eps")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_hhh_eps_eps_eps", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_P2P3")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_P2P3", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_S2P2P3")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc("B_S2P2P3", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    double patch_sqradians = spherical_cap_radius_2_sqradians(info_iB2D_W_FS.theta_T_2pt);
    result = integration_pre_factor * patch_sqradians * patch_sqradians * result;
    error = integration_pre_factor * patch_sqradians * patch_sqradians * error;
}

// Monte-Carlo angle averaged

double iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc_integrand(double *k, size_t dim, void *params)
{
    (void)(dim); // avoid unused parameter warnings

    double phi_l = k[0], z = k[1], l_1 = k[2], l_2 = k[3], phi_1 = k[4], phi_2 = k[5];

    params_iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_integrand *p = static_cast<params_iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_integrand *>(params);

    return evaluate_iB2D_z_l_1_l_2_phi_1_phi_2_integrand(p->key, p->l, phi_l, p->info_iB2D_W_FS, p->class_obj, p->use_pk_nl, z, l_1, l_2, phi_1, phi_2,
                                                         p->q1, p->q2, p->q3);

}

void iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                         projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                         const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls)
{
    params_iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_integrand args = {key, l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3};

    gsl_monte_function G = { &iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc_integrand, 6, static_cast<void *>(&args)};

    if (mc_integration_type == "plain" )
        monte_carlo_plain_integration(&G, lower_limits, upper_limits, 6, calls, T, result, error);
    else if (mc_integration_type == "miser" )
        monte_carlo_miser_integration(&G, lower_limits, upper_limits, 6, calls, T, result, error);
    else if (mc_integration_type == "vegas" )
        monte_carlo_vegas_integration(&G, lower_limits, upper_limits, 6, calls, T, result, error);
}

void iB2D_mc_angle_averaged(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
           projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
           const gsl_rng_type *T, const std::string &mc_integration_type, double &result, double &error, size_t calls)
{
    result = 0.0;
    error = 0.0;

    if (key == "B_lF2")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_lF2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "X_W")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("X_W", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "Y_W")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("Y_W", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "Z_W")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("Z_W", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "X_S2")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("X_S2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "Y_S2")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("Y_S2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "Z_S2")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("Z_S2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "P_1")
        // passing l = 0 explicitly as there is no l dependence for the integrals (depends only on l_1 and l_2)
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("P_1", 0, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "P_2")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("P_2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "P_3")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("P_3", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "A")
        // passing l = 0 explicitly as there is no l dependence for the integrals (depends only on l_1 and l_2)
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("A", 0, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    if (key == "B")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_xip_cos")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_xip_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_xip_sin")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_xip_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_xim_cos")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_xim_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_xim_sin")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_xim_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    // components for tracer integrated bispectrum

    else if (key == "B_nothing")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_nothing", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_P")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_P", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_PP")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_PP", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_S2PP")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_S2PP", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_hhh_delta_eps_eps")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_hhh_delta_eps_eps", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_hhh_eps_eps_eps")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_hhh_eps_eps_eps", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_P2P3")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_P2P3", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    else if (key == "B_S2P2P3")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_mc("B_S2P2P3", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, T, mc_integration_type, result, error, calls);

    double patch_sqradians = spherical_cap_radius_2_sqradians(info_iB2D_W_FS.theta_T_2pt);
    result = integration_pre_factor * patch_sqradians * patch_sqradians * result;
    error = integration_pre_factor * patch_sqradians * patch_sqradians * error;
}

// Monte-Carlo CIGAR

double iB2D_z_l_1_l_2_phi_1_phi_2_mc_cigar_integrand(std::vector<double> k, void *params)
{
    params_iB2D_z_l_1_l_2_phi_1_phi_2_mc_cigar_integrand *p = static_cast<params_iB2D_z_l_1_l_2_phi_1_phi_2_mc_cigar_integrand *>(params);

    // need the limits of the integrand because the cigar implementation of vegas only performs interation within the [0,1] hypercube
    double z_ll = p->lower_limits[0];
    double l_1_ll = p->lower_limits[1];
    double l_2_ll = p->lower_limits[2];
    double phi_1_ll = p->lower_limits[3];
    double phi_2_ll = p->lower_limits[4];

    double z_ul_m_ll = p->upper_limits[0]-z_ll;
    double l_1_ul_m_ll = p->upper_limits[1]-l_1_ll;
    double l_2_ul_m_ll = p->upper_limits[2]-l_2_ll;
    double phi_1_ul_m_ll = p->upper_limits[3]-phi_1_ll;
    double phi_2_ul_m_ll = p->upper_limits[4]-phi_2_ll;

    double z = z_ll + z_ul_m_ll*k[0];
    double l_1 = l_1_ll + l_1_ul_m_ll*k[1];
    double l_2 = l_2_ll + l_2_ul_m_ll*k[2];
    double phi_1 = phi_1_ll + phi_1_ul_m_ll*k[3];
    double phi_2 = phi_2_ll + phi_2_ul_m_ll*k[4];

    return evaluate_iB2D_z_l_1_l_2_phi_1_phi_2_integrand(p->key, p->l, p->phi_l, p->info_iB2D_W_FS, p->class_obj, p->use_pk_nl, z, l_1, l_2, phi_1, phi_2,
                                                         p->q1, p->q2, p->q3)*z_ul_m_ll*l_1_ul_m_ll*l_2_ul_m_ll*phi_1_ul_m_ll*phi_2_ul_m_ll;
}

void iB2D_z_l_1_l_2_phi_1_phi_2_mc_cigar(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                         projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                         VEGAS_Integrator &cigar, double &result, double &error, int thread_count)
{
    double phi_l = 0.;           // assuming angular-independency we can just fix the phi_l to any angle we want
//    double phi_l = 0.625*M_PI;           // assuming angular-independency we can just fix the phi_l to any angle we want
//    double phi_l = 1.378*M_PI;           // assuming angular-independency we can just fix the phi_l to any angle we want

    params_iB2D_z_l_1_l_2_phi_1_phi_2_mc_cigar_integrand args = {key, l, phi_l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits};

    //VEGAS_Integrator cigar;
    //cigar.Set_Verbose(NONE);
    cigar.Set_Integrand(iB2D_z_l_1_l_2_phi_1_phi_2_mc_cigar_integrand, 5, static_cast<void *>(&args));
    cigar.Improve_Grid();
    cigar.Integration();

    result = cigar.Get_Result();
    error = cigar.Get_Error();
}

void iB2D_mc_cigar(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                 projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                 VEGAS_Integrator &cigar, double &result, double &error, int thread_count)
{
    result = 0.0;
    error = 0.0;

    //VEGAS_Integrator cigar;
    //cigar.Set_Verbose(NONE);

    if (key == "B")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc_cigar("B", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, cigar, result, error, thread_count);

    else if (key == "B_xip_cos")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc_cigar("B_xip_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, cigar, result, error, thread_count);

    else if (key == "B_xip_sin")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc_cigar("B_xip_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, cigar, result, error, thread_count);

    else if (key == "B_xim_cos")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc_cigar("B_xim_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, cigar, result, error, thread_count);

    else if (key == "B_xim_sin")
        iB2D_z_l_1_l_2_phi_1_phi_2_mc_cigar("B_xim_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, cigar, result, error, thread_count);

    double patch_sqradians = spherical_cap_radius_2_sqradians(info_iB2D_W_FS.theta_T_2pt);
    result = integration_pre_factor * patch_sqradians * patch_sqradians * result;
    error = integration_pre_factor * patch_sqradians * patch_sqradians * error;
}

// h-cubature

int iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value)
{
    assert(ndim == 5);
    assert(fdim == 1);

    double z = k[0], l_1 = k[1], l_2 = k[2], phi_1 = k[3], phi_2 = k[4];

    params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand *p = static_cast<params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand *>(params);

    value[0] = evaluate_iB2D_z_l_1_l_2_phi_1_phi_2_integrand(p->key, p->l, p->phi_l, p->info_iB2D_W_FS, p->class_obj, p->use_pk_nl, z, l_1, l_2, phi_1, phi_2,
                                                             p->q1, p->q2, p->q3);

    return 0;
}

void iB2D_z_l_1_l_2_phi_1_phi_2_hcubature(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                          projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                          double &result, double &error, size_t max_evals)
{
    double phi_l = 0.0;           // assuming angular-independency we can just fix the phi_l to any angle we want
//    double phi_l = 0.625*M_PI;           // assuming angular-independency we can just fix the phi_l to any angle we want
//    double phi_l = 1.378*M_PI;           // assuming angular-independency we can just fix the phi_l to any angle we want

    params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand args = {key, l, phi_l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3};

    hcubature_integration(iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_integrand, static_cast<void *>(&args), lower_limits, upper_limits, 5, max_evals, result, error);
}

void iB2D_hcubature(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                  projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                  double &result, double& error, size_t max_evals)
{
    result = 0.0;
    error = 0.0;

    if (key == "B_lF2")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("B_lF2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "X_W")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("X_W", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "Y_W")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("Y_W", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "Z_W")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("Z_W", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "X_S2")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("X_S2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "Y_S2")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("Y_S2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "Z_S2")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("Z_S2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "P_1")
        // passing l = 0 explicitly as there is no l dependence for the integrals (depends only on l_1 and l_2)
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("P_1", 0, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "P_2")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("P_2", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "P_3")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("P_3", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "A")
        // passing l = 0 explicitly as there is no l dependence for the integrals (depends only on l_1 and l_2)
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("A", 0, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_lF2_xip_cos")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("B_lF2_xip_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_lF2_xip_sin")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("B_lF2_xip_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_lF2_xim_cos")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("B_lF2_xim_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_lF2_xim_sin")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("B_lF2_xim_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    if (key == "B")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("B", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xip_cos")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("B_xip_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xip_sin")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("B_xip_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xim_cos")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("B_xim_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xim_sin")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature("B_xim_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    double patch_sqradians = spherical_cap_radius_2_sqradians(info_iB2D_W_FS.theta_T_2pt);
    result = integration_pre_factor * patch_sqradians * patch_sqradians * result;
    error = integration_pre_factor * patch_sqradians * patch_sqradians * error;
}

// h-cubature with 4 dimensions

int iB2D_z_l_1_l_2_phi_1_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value)
{
    assert(ndim == 4);
    assert(fdim == 1);

    double z = k[0], l_1 = k[1], l_2 = k[2], phi_2 = k[3];

    double phi_1 = 0;

    params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand *p = static_cast<params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand *>(params);

    value[0] = evaluate_iB2D_z_l_1_l_2_phi_1_phi_2_integrand(p->key, p->l, p->phi_l, p->info_iB2D_W_FS, p->class_obj, p->use_pk_nl, z, l_1, l_2, phi_1, phi_2,
                                                             p->q1, p->q2, p->q3);

    return 0;
}

void iB2D_z_l_1_l_2_phi_1_hcubature(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                    projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits,
                                    std::vector<double> upper_limits, double &result, double &error, size_t max_evals)
{
    double phi_l = 0;           // assuming angular-independency we can just fix the phi_l to any angle we want

    params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand args = {key, l, phi_l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3};

    hcubature_integration(iB2D_z_l_1_l_2_phi_1_hcubature_integrand, static_cast<void *>(&args), lower_limits, upper_limits, 4, max_evals, result, error);
}

void iB2D_hcubature_4dim(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                       projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                       double &result, double& error, size_t max_evals)
{
    result = 0.0;
    error = 0.0;

    if (key == "B")
        iB2D_z_l_1_l_2_phi_1_hcubature("B", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xip_cos")
        iB2D_z_l_1_l_2_phi_1_hcubature("B_xip_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xip_sin")
        iB2D_z_l_1_l_2_phi_1_hcubature("B_xip_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xim_cos")
        iB2D_z_l_1_l_2_phi_1_hcubature("B_xim_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xim_sin")
        iB2D_z_l_1_l_2_phi_1_hcubature("B_xim_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    double patch_sqradians = spherical_cap_radius_2_sqradians(info_iB2D_W_FS.theta_T_2pt);
    result = integration_pre_factor * 2*M_PI* patch_sqradians * patch_sqradians * result;
    error = integration_pre_factor * 2*M_PI* patch_sqradians * patch_sqradians * error;
}


// h-cubature vectorized

int iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_v_integrand(unsigned ndim, size_t npts, const double *k, void *params, unsigned fdim, double *value)
{
    assert(ndim == 5);
    assert(fdim == 1);

    double z, l_1, l_2, phi_1, phi_2;

    params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand *p = static_cast<params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand *>(params);

    for (unsigned j = 0; j < npts; ++j) { // evaluate the integrand for npts points

        z = k[j*ndim+0], l_1 = k[j*ndim+1], l_2 = k[j*ndim+2], phi_1 = k[j*ndim+3], phi_2 = k[j*ndim+4];

        value[j] = evaluate_iB2D_z_l_1_l_2_phi_1_phi_2_integrand(p->key, p->l, p->phi_l, p->info_iB2D_W_FS, p->class_obj, p->use_pk_nl, z, l_1, l_2, phi_1, phi_2,
                                                                 p->q1, p->q2, p->q3);
    }

    return 0;
}

void iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_v(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                            projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                            double &result, double &error, size_t max_evals)
{
    double phi_l = 0;           // assuming angular-independency we can just fix the phi_l to any angle we want

    params_iB2D_z_l_1_l_2_phi_1_phi_2_integrand args = {key, l, phi_l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3};

    hcubature_v_integration(iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_v_integrand, static_cast<void *>(&args), lower_limits, upper_limits, 5, max_evals, result, error);
}

double iB2D_hcubature_v(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                      projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                      size_t max_evals, double &result, double& error)
{
    result = 0.0;
    error = 0.0;

    if (key == "B")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_v("B", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xip_cos")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_v("B_xip_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xip_sin")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_v("B_xip_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xim_cos")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_v("B_xim_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xim_sin")
        iB2D_z_l_1_l_2_phi_1_phi_2_hcubature_v("B_xim_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    double patch_sqradians = spherical_cap_radius_2_sqradians(info_iB2D_W_FS.theta_T_2pt);
    result = integration_pre_factor * patch_sqradians * patch_sqradians * result;
    error = integration_pre_factor * patch_sqradians * patch_sqradians * error;

    return result;
}

// h-cubature angle averaged

int iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value)
{
    assert(ndim == 6);
    assert(fdim == 1);

    double phi_l = k[0], z = k[1], l_1 = k[2], l_2 = k[3], phi_1 = k[4], phi_2 = k[5];

    params_iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_integrand *p = static_cast<params_iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_integrand *>(params);

    value[0] = evaluate_iB2D_z_l_1_l_2_phi_1_phi_2_integrand(p->key, p->l, phi_l, p->info_iB2D_W_FS, p->class_obj, p->use_pk_nl, z, l_1, l_2, phi_1, phi_2,
                                                             p->q1, p->q2, p->q3);

    return 0;
}

void iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_hcubature(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                                projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                                double &result, double &error, size_t max_evals)
{
    params_iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_integrand args = {key, l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3};

    hcubature_integration(iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_hcubature_integrand, static_cast<void *>(&args), lower_limits, upper_limits, 6, max_evals, result, error);
}

void iB2D_hcubature_angle_averaged(const std::string &key, const double &l, const struct_iB2D_W_FS &info_iB2D_W_FS, ClassEngine *class_obj, const bool &use_pk_nl,
                                 projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> lower_limits, std::vector<double> upper_limits,
                                 double &result, double &error, size_t max_evals)
{
    result = 0.0;
    error = 0.0;

    if (key == "B")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_hcubature("B", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xip_cos")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_hcubature("B_xip_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xip_sin")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_hcubature("B_xip_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xim_cos")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_hcubature("B_xim_cos", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    else if (key == "B_xim_sin")
        iB2D_phi_l_z_l_1_l_2_phi_1_phi_2_hcubature("B_xim_sin", l, info_iB2D_W_FS, class_obj, use_pk_nl, q1, q2, q3, lower_limits, upper_limits, result, error, max_evals);

    double patch_sqradians = spherical_cap_radius_2_sqradians(info_iB2D_W_FS.theta_T_2pt);
    result = integration_pre_factor * patch_sqradians * patch_sqradians * result / (2*M_PI);
    error = integration_pre_factor * patch_sqradians * patch_sqradians * error / (2*M_PI);
}
