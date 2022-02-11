#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <math.h>

const size_t calls_1e2  = 100;
const size_t calls_1e3  = 1000;
const size_t calls_1e4  = 10000;
const size_t calls_1e5  = 100000;
const size_t calls_1e6  = 1000000;
const size_t calls_1e7  = 10000000;
const size_t calls_1e8  = 100000000;
const size_t calls_1e9  = 1000000000;

const double error_m01 = 1.0e-1;
const double error_m02 = 1.0e-2;
const double error_m03 = 1.0e-3;
const double error_m04 = 1.0e-4;
const double error_m05 = 1.0e-5;
const double error_m06 = 1.0e-6;
const double error_m07 = 1.0e-7;
const double error_m08 = 1.0e-8;
const double error_m10 = 1.0e-10;

const double z_cmb = 1090;

const double arcmin = 1/60.0*M_PI/180.0; // expressed in radians

const double delta_z_step = 0.02; // for making grid
const int num_trapz_steps = 100; // for trapezoidal integration

// Takahashi et al. (2017) Appendix B eqn (28) fitting formula for finite shell thickness - parameter values
const double c1_T17 = 9.5171e-4;
const double c2_T17 = 5.1543e-3;
const double a1_T17 = 1.3063;
const double a2_T17 = 1.1475;
const double a3_T17 = 0.62793;
const double ell_res_T17 = 1.6*4096;
const double delta_r_T17 = 450.0; // [Mpc/h]

//const double delta_r_T17 = 126.0; // [Mpc/h] for MassiveNus

const bool apply_T17_corrections = false; // set this to false when calculating models for Fisher derivatives; when comparing to simulations set to true

const bool compute_bihalofit = false;

const bool verbose_print_outs = true;

#endif // CONSTANTS_H
