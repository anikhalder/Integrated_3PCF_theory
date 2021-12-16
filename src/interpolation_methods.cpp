#include "interpolation_methods.h"
#include <mlinterp.hpp>
#include <assert.h>

/*
 *
// ######################################################################################

// From Numerical Recipes interp_1d.h

Base_interp_1D_NR::Base_interp_1D_NR(std::vector<double> &x, const double *y, int m) : n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y)
{
    dj = std::min(1,(int)pow((double)n,0.25));
}

double Base_interp_1D_NR::interp(double x)
{
    int jlo = cor ? hunt(x) : locate(x);
    return rawinterp(jlo,x);
}

int Base_interp_1D_NR::locate(const double x)
{
    int ju,jm,jl;
    if (n < 2 || mm < 2 || mm > n)
        throw("locate size error");

    bool ascnd=(xx[n-1] >= xx[0]);

    jl=0;
    ju=n-1;
    while (ju-jl > 1)
    {
        jm = (ju+jl) >> 1;
        if (x >= xx[jm] == ascnd)
            jl=jm;
        else
            ju=jm;
    }
    cor = abs(jl-jsav) > dj ? 0 : 1;
    jsav = jl;
    return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
}

int Base_interp_1D_NR::hunt(const double x)
{
    int jl=jsav, jm, ju, inc=1;
    if (n < 2 || mm < 2 || mm > n)
        throw("hunt size error");

    bool ascnd=(xx[n-1] >= xx[0]);

    if (jl < 0 || jl > n-1)
    {
        jl=0;
        ju=n-1;
    }
    else
    {
        if (x >= xx[jl] == ascnd)
        {
            for (;;)
            {
                ju = jl + inc;
                if (ju >= n-1)
                {
                    ju = n-1;
                    break;
                }
                else if (x < xx[ju] == ascnd)
                    break;
                else
                {
                    jl = ju;
                    inc += inc;
                }
            }
        }
        else
        {
            ju = jl;
            for (;;)
            {
                jl = jl - inc;
                if (jl <= 0)
                {
                    jl = 0;
                    break;
                }
                else if (x >= xx[jl] == ascnd)
                    break;
                else
                {
                    ju = jl;
                    inc += inc;
                }
            }
        }
    }

    while (ju-jl > 1)
    {
        jm = (ju+jl) >> 1;
        if (x >= xx[jm] == ascnd)
            jl=jm;
        else
            ju=jm;
    }

    cor = abs(jl-jsav) > dj ? 0 : 1;
    jsav = jl;
    return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
}

Base_interp_1D_NR::~Base_interp_1D_NR()
{

}

Linear_interp_1D_NR::Linear_interp_1D_NR(std::vector<double> &xv, std::vector<double> &vals)
    : Base_interp_1D_NR(xv,&vals[0],2)  {}

double Linear_interp_1D_NR::rawinterp(int j, double x)
{
    if (xx[j]==xx[j+1])
        return yy[j];
    else
        return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
}

Linear_interp_1D_NR::~Linear_interp_1D_NR()
{

}

*/

// ######################################################################################

// Using N-dim linear interpolation template from mlinterp --> see https://github.com/parsiad/mlinterp

// ######################################################################################

// Linear interpolation (over a 1D rectilinear 1D grid)

Linear_interp_1D::Linear_interp_1D()
{
    xv = {};
    vals = {};

    grid_dim.push_back(xv.size());
}

Linear_interp_1D::Linear_interp_1D(const std::vector<double> &x, const std::vector<double> &v) :
    xv(x), vals(v)
{
    grid_dim.push_back(xv.size());
}

std::vector<double> Linear_interp_1D::get_min_max_x_coordinates()
{
    std::vector<double> coordinates;
    coordinates.push_back(xv.front());
    coordinates.push_back(xv.back());

    return coordinates;
}

double Linear_interp_1D::interp(const double &xi)
{
    size_t ni = 1;

    std::vector<double> interpolated_value(ni,0);

    mlinterp::interp(grid_dim.data(), ni, vals.data(), interpolated_value.data(), xv.data(), &xi);

    return interpolated_value.at(0);
}

std::vector<double> Linear_interp_1D::interp_vec(const std::vector<double> &xi)
{

    size_t ni = xi.size();

    std::vector<double> interpolated_values(ni,0);

    mlinterp::interp(grid_dim.data(), ni, vals.data(), interpolated_values.data(), xv.data(), xi.data());

    return interpolated_values;
}

Linear_interp_1D::~Linear_interp_1D()
{

}

// ######################################################################################

// Bilinear interpolation (over a 2D rectilinear 2D grid)

Linear_interp_2D::Linear_interp_2D()
{
    xv = {};
    yv = {};
    vals = {};

    grid_dim.push_back(xv.size());
    grid_dim.push_back(yv.size());
}

Linear_interp_2D::Linear_interp_2D(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &v) :
    xv(x), yv(y), vals(v)
{
    grid_dim.push_back(xv.size());
    grid_dim.push_back(yv.size());
}

std::vector<double> Linear_interp_2D::get_min_max_x_coordinates()
{
    std::vector<double> coordinates;
    coordinates.push_back(xv.front());
    coordinates.push_back(xv.back());

    return coordinates;
}

std::vector<double> Linear_interp_2D::get_min_max_y_coordinates()
{
    std::vector<double> coordinates;
    coordinates.push_back(yv.front());
    coordinates.push_back(yv.back());

    return coordinates;
}

double Linear_interp_2D::interp(const double &xi, const double &yi)
{
    size_t ni = 1;

    std::vector<double> interpolated_value(ni,0);

    mlinterp::interp(grid_dim.data(), ni, vals.data(), interpolated_value.data(), xv.data(), &xi, yv.data(), &yi);

    return interpolated_value.at(0);
}

std::vector<double> Linear_interp_2D::interp_vec(const std::vector<double> &xi, const std::vector<double> &yi)
{
    assert(xi.size() == yi.size());

    size_t ni = xi.size();

    std::vector<double> interpolated_values(ni,0);

    mlinterp::interp(grid_dim.data(), ni, vals.data(), interpolated_values.data(), xv.data(), xi.data(), yv.data(), yi.data());

    return interpolated_values;
}

Linear_interp_2D::~Linear_interp_2D()
{

}

// ######################################################################################

// Trilinear interpolation (over a 3D rectilinear 3D grid)

Linear_interp_3D::Linear_interp_3D()
{
    xv = {};
    yv = {};
    zv = {};
    vals = {};

    grid_dim.push_back(xv.size());
    grid_dim.push_back(yv.size());
    grid_dim.push_back(zv.size());
}

Linear_interp_3D::Linear_interp_3D(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, const std::vector<double> &v) :
    xv(x), yv(y), zv(z), vals(v)
{
    grid_dim.push_back(xv.size());
    grid_dim.push_back(yv.size());
    grid_dim.push_back(zv.size());
}

std::vector<double> Linear_interp_3D::get_min_max_x_coordinates()
{
    std::vector<double> coordinates;
    coordinates.push_back(xv.front());
    coordinates.push_back(xv.back());

    return coordinates;
}

std::vector<double> Linear_interp_3D::get_min_max_y_coordinates()
{
    std::vector<double> coordinates;
    coordinates.push_back(yv.front());
    coordinates.push_back(yv.back());

    return coordinates;
}

std::vector<double> Linear_interp_3D::get_min_max_z_coordinates()
{
    std::vector<double> coordinates;
    coordinates.push_back(zv.front());
    coordinates.push_back(zv.back());

    return coordinates;
}

double Linear_interp_3D::interp(const double &xi, const double &yi, const double &zi)
{
    size_t ni = 1;

    std::vector<double> interpolated_value(ni,0);

    mlinterp::interp(grid_dim.data(), ni, vals.data(), interpolated_value.data(), xv.data(), &xi, yv.data(), &yi, zv.data(), &zi);

    return interpolated_value.at(0);
}

std::vector<double> Linear_interp_3D::interp_vec(const std::vector<double> &xi, const std::vector<double> &yi, const std::vector<double> &zi)
{
    assert(xi.size() == yi.size() && yi.size() == zi.size());

    size_t ni = xi.size();

    std::vector<double> interpolated_values(ni,0);

    mlinterp::interp(grid_dim.data(), ni, vals.data(), interpolated_values.data(), xv.data(), xi.data(), yv.data(), yi.data(), zv.data(), zi.data());

    return interpolated_values;
}

Linear_interp_3D::~Linear_interp_3D()
{

}
