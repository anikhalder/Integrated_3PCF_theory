#ifndef POWER_SPECTRUM_HPP
#define POWER_SPECTRUM_HPP

#include <ClassEngine.hh>

// ######################################################################################

// power spectrum prescription

double P(const double &k, const double &z, ClassEngine *class_obj, bool use_pk_nl);

// ######################################################################################

// tree-level SPT matter power spectrum

double P_tree(const double &k, const double &z, ClassEngine *class_obj);

// ######################################################################################

// non-linear matter power spectrum fitting formula

double P_nl(const double &k, const double &z, ClassEngine *class_obj);

#endif // POWER_SPECTRUM_HPP
