#include <power_spectrum.hpp>
#include <constants.h>
#include <cosmology_utils.h>

// ######################################################################################

// power spectrum prescription

double P(const double &k, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    if (z > class_obj->get_z_max_pk())
        return 0.0;

    double k_max_pk = class_obj->get_k_max_pk();
    if (k > k_max_pk)
        return 0.0;

    if (k == 0)
        return 0.0;

    // -----------------------
    // Tests for shell thickness integration

    //return class_obj->pk_nl(k,z); // without shell thickness

    //return class_obj->pk_nl(k,z)*pow(T17_shell_correction(k,class_obj),2); // shell thickness - T17 fitting function

    //double delta_r = delta_r_T17 / class_obj->get_h();
    //return Pk_shell_correction_qag(k, z, delta_r, class_obj); // shell thickness - T17 integration

    // -----------------------

    if (apply_T17_corrections && T17_box_correction(k,z,class_obj) == 0)
        return 0.0;

    double P;
    if (use_pk_nl)
        P = P_nl(k, z, class_obj);
    else
        P = P_tree(k, z, class_obj);

    if (apply_T17_corrections)
        P *= pow(T17_shell_correction(k,class_obj),2);

    return P;
}

// ######################################################################################

// tree-level SPT matter power spectrum

double P_tree(const double &k, const double &z, ClassEngine *class_obj)
{
    return class_obj->pk_lin(k,z);
}

// ######################################################################################

// non-linear matter power spectrum fitting formula

double P_nl(const double &k, const double &z, ClassEngine *class_obj)
{
    return class_obj->pk_nl(k,z);
}
