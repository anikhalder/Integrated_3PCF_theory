#include <power_spectra_2D.h>
#include <power_spectrum.hpp>
#include <cosmology_utils.h>
#include <integration_utils.h>
#include <constants.h>
#include <cubature.h>
#include <assert.h>

// ######################################################################################

// 2D power spectrum (integrated along the line-of-sight)

/* e.g. for Pmm (e.g. LSS)
 *
 * eqn (2.84) of https://arxiv.org/pdf/astro-ph/9912508.pdf
 *
 * P(l) = \int dx q_1(x)q_2(x)/(x^2) P_3D(l/x, eta0 - x)
 * => P(l) = \int dz/H(z) q_1(x(z))q_2(x(z))/(x(z)^2) P_3D(l/x(z), eta0 - x(z))
 * with q_1(x(z)) = q_2(x(z) = q_m(z, class_obj, n_m_of_z) */


/* e.g. for Pkk (e.g. CMB)
 *
 * eqn (2.84) of https://arxiv.org/pdf/astro-ph/9912508.pdf
 *
 * P(l) = \int dx q_1(x)q_2(x)/(x^2) P_3D(l/x, eta0 - x)
 * => P(l) = \int dz/H(z) q_1(x(z))q_2(x(z))/(x(z)^2) P_3D(l/x(z), eta0 - x(z))
 * with q_1(x(z)) = q_2(x(z) = q_k_zs_fixed(z, class_obj, z_cmb) i.e. eqn (6.22) of https://arxiv.org/pdf/astro-ph/9912508.pdf */


double evaluate_P_los_integrand(const std::string &key, const double &l, ClassEngine *class_obj, bool use_pk_nl, const double &z,
                                projection_kernel *q1, projection_kernel *q2)
{
    double q_1 = q1->evaluate(z);
    double q_2 = q2->evaluate(z);

    if (q_1 == 0 || q_2 == 0)
        return 0;

    double chi_inv = 1/class_obj->get_chi_z(z);
    double Hz_inv = 1/class_obj->get_H_z(z);

    if (key == "P")
        return Hz_inv*chi_inv*chi_inv*q_1*q_2*P(l*chi_inv,z,class_obj,use_pk_nl);

    else if (key == "P_nothing")
        return Hz_inv*chi_inv*chi_inv*q_1*q_2;

    else if (key == "P_hh_eps_eps")
    {
        double _1_over_n_h_z = 1./dynamic_cast<projection_kernel_q_h*>(q2)->get_n_h_z(z);
        return Hz_inv*chi_inv*chi_inv*q_1*q_2*_1_over_n_h_z;
    }

    else
        return Hz_inv*chi_inv*chi_inv*q_1*q_2*P(l*chi_inv,z,class_obj,use_pk_nl);
}

// Gaussian-quadrature

double P_los_qag_integrand(double z, void *params)
{
    params_P_los_integrand *p = static_cast<params_P_los_integrand *>(params);

    return evaluate_P_los_integrand(p->key, p->l, p->class_obj, p->use_pk_nl, z, p->q1, p->q2);
}

double P_los_qag(const std::string &key, const double &l, ClassEngine *class_obj, const bool &use_pk_nl,
                 projection_kernel *q1, projection_kernel *q2, const double &z_lower, const double &z_upper)
{
    double result = 0, error = 0;

    //parameters in integrand
    params_P_los_integrand args = {key, l, class_obj, use_pk_nl, q1, q2};

    qag_1D_integration(&P_los_qag_integrand, static_cast<void *>(&args), z_lower, z_upper, calls_1e5, result, error);

    return result;
}

// Monte-Carlo

double P_los_mc_integrand(double *k, size_t dim, void *params)
{
    (void)(dim); // avoid unused parameter warnings

    params_P_los_integrand *p = static_cast<params_P_los_integrand *>(params);

    double z = k[0];

    return evaluate_P_los_integrand(p->key, p->l, p->class_obj, p->use_pk_nl, z, p->q1, p->q2);
}

double P_los_mc(const std::string &key, const double &l, ClassEngine *class_obj, const bool &use_pk_nl,
                projection_kernel *q1, projection_kernel *q2, std::vector<double> &z_lower, std::vector<double> &z_upper,
                const gsl_rng_type *T, const std::string &mc_integration_type)
{
    double result = 0, error = 0;

    //parameters in integrand
    params_P_los_integrand args = {key, l, class_obj, use_pk_nl, q1, q2};

    gsl_monte_function G = { &P_los_mc_integrand, 1, static_cast<void *>(&args)};

    size_t calls = 2*calls_1e3;

    if (mc_integration_type == "plain" )
        monte_carlo_plain_integration(&G, z_lower, z_upper, 1, calls, T, result, error);
    else if (mc_integration_type == "miser" )
        monte_carlo_miser_integration(&G, z_lower, z_upper, 1, calls, T, result, error);
    else if (mc_integration_type == "vegas" )
        monte_carlo_vegas_integration(&G, z_lower, z_upper, 1, calls, T, result, error);

    return result;
}

// h-cubature

int P_los_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value)
{
    assert(ndim == 1);
    assert(fdim == 1);

    params_P_los_integrand *p = static_cast<params_P_los_integrand *>(params);

    double z = k[0];

    value[0] = evaluate_P_los_integrand(p->key, p->l, p->class_obj, p->use_pk_nl, z, p->q1, p->q2);

    return 0;
}

double P_los_hcubature(const std::string &key, const double &l, ClassEngine *class_obj, const bool &use_pk_nl,
                       projection_kernel *q1, projection_kernel *q2, std::vector<double> &z_lower, std::vector<double> &z_upper)
{
    double result = 0, error = 0;

    //parameters in integrand
    params_P_los_integrand args = {key, l, class_obj, use_pk_nl, q1, q2};

    hcubature_integration(P_los_hcubature_integrand, static_cast<void *>(&args), z_lower, z_upper, 1, 0, result, error);

    return result;
}
