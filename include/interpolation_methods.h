#ifndef INTERPOLATION_METHODS_H
#define INTERPOLATION_METHODS_H

#include <iostream>
#include <vector>
#include <math.h>

/*

// ######################################################################################

// From Numerical Recipes interp_1d.h

class Base_interp_1D_NR
{
public:
    int n, mm, jsav, cor, dj;
    const double *xx, *yy;
    Base_interp_1D_NR(std::vector<double> &x, const double *y, int m);

    double interp(double x);

    int locate(const double x);
    int hunt(const double x);

    double virtual rawinterp(int jlo, double x) = 0;

    virtual ~Base_interp_1D_NR() = 0;
};


class Linear_interp_1D_NR : public Base_interp_1D_NR
{
public:
    Linear_interp_1D_NR(std::vector<double> &xv, std::vector<double> &vals);

    double rawinterp(int j, double x);

    ~Linear_interp_1D_NR();
};

*/

// ######################################################################################

// Using N-dim linear interpolation template from mlinterp --> see https://github.com/parsiad/mlinterp

// ######################################################################################

// Linear interpolation (over a rectilinear 1D grid)

class Linear_interp_1D
{
    std::vector<size_t> grid_dim; // number of points along the x axis of the rectilinear 1D grid
    std::vector<double> xv, vals; // the vertices along the x axis of the rectilinear 1D grid

public:
    Linear_interp_1D();
    Linear_interp_1D(const std::vector<double> &x, const std::vector<double> &v);

    std::vector<double> get_min_max_x_coordinates();

    double interp(const double &xi); // interpolation for a single point on the rectilinear 1D grid

    std::vector<double> interp_vec(const std::vector<double> &xi);  // interpolation for a vector of points on the rectilinear 1D grid

    ~Linear_interp_1D();
};

// ######################################################################################

// Bilinear interpolation (over a rectilinear 2D grid)

class Linear_interp_2D
{
    std::vector<size_t> grid_dim; // number of points along each of the x and y axes of the rectilinear 2D grid
    std::vector<double> xv, yv, vals; // the vertices along each of the x and y axes of the rectilinear 2D grid

public:
    Linear_interp_2D();
    Linear_interp_2D(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &v);

    std::vector<double> get_min_max_x_coordinates();
    std::vector<double> get_min_max_y_coordinates();

    double interp(const double &xi, const double &yi); // interpolation for a single point on the rectilinear 2D grid

    std::vector<double> interp_vec(const std::vector<double> &xi, const std::vector<double> &yi);  // interpolation for a vector of points on the rectilinear 2D grid

    ~Linear_interp_2D();
};

// ######################################################################################

// Trilinear interpolation (over a rectilinear 3D grid)

class Linear_interp_3D
{
    std::vector<size_t> grid_dim; // number of points along each of the x, y and z axes of the rectilinear 3D grid
    std::vector<double> xv, yv, zv, vals; // the vertices along each of the x, y and z axes of the rectilinear 3D grid

public:
    Linear_interp_3D();
    Linear_interp_3D(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, const std::vector<double> &v);

    std::vector<double> get_min_max_x_coordinates();
    std::vector<double> get_min_max_y_coordinates();
    std::vector<double> get_min_max_z_coordinates();

    double interp(const double &xi, const double &yi, const double &zi); // interpolation for a single point on the rectilinear 3D grid

    std::vector<double> interp_vec(const std::vector<double> &xi, const std::vector<double> &yi,
                                   const std::vector<double> &zi);  // interpolation for a vector of points on the rectilinear 3D grid

    ~Linear_interp_3D();
};

#endif // INTERPOLATION_METHODS_H
