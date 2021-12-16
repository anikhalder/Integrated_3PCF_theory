#include <bispectra_2D.h>
#include <bispectrum.hpp>
#include <cosmology_utils.h>
#include <integration_utils.h>
#include <constants.h>
#include <cubature.h>
#include <assert.h>

// ######################################################################################

// 2D bispectrum (integrated along the line-of-sight)

/*
 * B(l1,l2,l3) = \int dx q_1(x)q_2(x)q_3(x)/(x^4) B_3D(l1/x, l2/x, l3/x, eta0 - x)
 * => P(l) = \int dz/H(z) q_1(x(z))q_2(x(z))/(x(z)^3) B_3D(l1/x(z), l2/x(z), l3/x(z), eta0 - x(z)) */

double evaluate_B_los_integrand(const std::string &key, const double &l1, const double &l2, const double &l3, ClassEngine *class_obj, bool use_pk_nl, const double &z,
                                const double &q_1, const double &q_2, const double &q_3)
{
    if (q_1 == 0 || q_2 == 0 || q_3 == 0)
        return 0;

    double chi_inv = 1/class_obj->get_chi_z(z);
    double Hz_inv = 1/class_obj->get_H_z(z);

    if (key == "B")
        return Hz_inv*pow(chi_inv,4)*q_1*q_2*q_3*B(l1*chi_inv,l2*chi_inv,l3*chi_inv,z,class_obj,use_pk_nl);

    else if (key == "B_nothing")
        return Hz_inv*pow(chi_inv,4)*q_1*q_2*q_3;

    else
        return Hz_inv*pow(chi_inv,4)*q_1*q_2*q_3*B(l1*chi_inv,l2*chi_inv,l3*chi_inv,z,class_obj,use_pk_nl);

}

// Gaussian-quadrature

double B_los_qag_integrand(double z, void *params)
{
    params_B_los_integrand *p = static_cast<params_B_los_integrand *>(params);

    return evaluate_B_los_integrand(p->key, p->l1, p->l2, p->l3, p->class_obj, p->use_pk_nl, z, p->q1->evaluate(z), p->q2->evaluate(z), p->q3->evaluate(z));
}

double B_los_qag(const std::string &key, const double &l1, const double &l2, const double &l3, ClassEngine *class_obj, const bool &use_pk_nl,
                 projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, const double &z_lower, const double &z_upper)
{
    double result = 0, error = 0;

    //parameters in integrand
    params_B_los_integrand args = {key, l1, l2, l3, class_obj, use_pk_nl, q1, q2, q3};

    qag_1D_integration(&B_los_qag_integrand, static_cast<void *>(&args), z_lower, z_upper, calls_1e5, result, error);

    return result;
}

// Monte-Carlo

double B_los_mc_integrand(double *k, size_t dim, void *params)
{
    (void)(dim); // avoid unused parameter warnings

//    double z = k[0];

    params_B_los_integrand *p = static_cast<params_B_los_integrand *>(params);

    return evaluate_B_los_integrand(p->key, p->l1, p->l2, p->l3, p->class_obj, p->use_pk_nl, k[0], p->q1->evaluate(k[0]), p->q2->evaluate(k[0]), p->q3->evaluate(k[0]));
}

double B_los_mc(const std::string &key, const double &l1, const double &l2, const double &l3, ClassEngine *class_obj, const bool &use_pk_nl,
                projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> &z_lower, std::vector<double> &z_upper,
                const gsl_rng_type *T, const std::string &mc_integration_type)
{
    double result = 0, error = 0;

    //parameters in integrand
    params_B_los_integrand args = {key, l1, l2, l3, class_obj, use_pk_nl, q1, q2, q3};

    gsl_monte_function G = { &B_los_mc_integrand, 1, static_cast<void *>(&args)};

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

int B_los_hcubature_integrand(unsigned ndim, const double *k, void *params, unsigned fdim, double *value)
{
    assert(ndim == 1);
    assert(fdim == 1);

    params_B_los_integrand *p = static_cast<params_B_los_integrand *>(params);

    double z = k[0];

    value[0] = evaluate_B_los_integrand(p->key, p->l1, p->l2, p->l3, p->class_obj, p->use_pk_nl, z, p->q1->evaluate(z), p->q2->evaluate(z), p->q3->evaluate(z));

    return 0;
}

double B_los_hcubature(const std::string &key, const double &l1, const double &l2, const double &l3, ClassEngine *class_obj, const bool &use_pk_nl,
                       projection_kernel *q1, projection_kernel *q2, projection_kernel *q3, std::vector<double> &z_lower, std::vector<double> &z_upper)
{
    double result = 0, error = 0;

    //parameters in integrand
    params_B_los_integrand args = {key, l1, l2, l3, class_obj, use_pk_nl, q1, q2, q3};

    hcubature_integration(B_los_hcubature_integrand, static_cast<void *>(&args), z_lower, z_upper, 1, 0, result, error);

    return result;
}
