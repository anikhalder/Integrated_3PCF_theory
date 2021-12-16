#include <real_space_2D.h>
#include <interpolation_methods.h>
#include <gsl/gsl_sf.h>
#include <constants.h>
#include <cosmology_utils.h>

void C_ells_spherical_sky(const std::vector<double> &l_array, const std::vector<double> &P, std::vector<int> &ells, std::vector<double> &C_ells)
{
    /* In this function, we use the following notation:
     * l_array and P array are computed for 2D flat sky at discrete 2D Fourier l values until some maximum l_max (which can be a floating number)
     * ells and C_ells array are computed for 2D spherical sky at every integer multipole ells = 0,1,...,l_max-1
     */

    // compute interpolated C_ells array for ells = 0,1,...,l_max-1

    Linear_interp_1D P_interpolated_obj(l_array, P); // l is assumed to start from 1

    int l_max = int(l_array.back());

    for (int l = 0; l <= l_max-1; l++)
    {
        ells.push_back(l);

        if (l == 0 || l == 1)
            C_ells.push_back(0.);
        else
        {
            //C_ells.push_back(P_interpolated_obj.interp(ell));
            C_ells.push_back((l+2.)*(l+1.)*l*(l-1.) / pow(l+0.5,4) * P_interpolated_obj.interp(l+0.5)); // Kitching++ correction
        }

        if(apply_T17_corrections)
            C_ells[l] *= T17_resolution_correction(l);
    }
}

double G_ell_2_x_p(const int &l, const double &x)
{   
    // eqn (A5) of Friedrich++ 2020 DES covariance modelling --- https://arxiv.org/pdf/2012.08568.pdf
    // here, G_ell_2_x_p = G_ell_2_+(x) + G_ell_2_-(x) i.e. LHS of Friedrich ++ for the 'p'lus middle sign

    double P_l_2 =  gsl_sf_legendre_Plm(l,2,x);
    double P_lm1_2 = gsl_sf_legendre_Plm(l-1,2,x) ;

    return P_l_2*( (4.-l+2.*x*(l-1.))/(1.-x*x) - l*(l-1.)/2. ) + P_lm1_2*( (l+2.)*(x-2.)/(1.-x*x));
}

double G_ell_2_x_m(const int &l, const double &x)
{
    // eqn (A5) of Friedrich++ 2020 DES covariance modelling --- https://arxiv.org/pdf/2012.08568.pdf
    // here, G_ell_2_x_m = G_ell_2_+(x) - G_ell_2_-(x) i.e. LHS of Friedrich ++ for the 'm'inus middle sign

    double P_l_2 =  gsl_sf_legendre_Plm(l,2,x);
    double P_lm1_2 = gsl_sf_legendre_Plm(l-1,2,x) ;

    return P_l_2*( (4.-l-2.*x*(l-1.))/(1.-x*x) - l*(l-1.)/2. ) + P_lm1_2*( (l+2.)*(x+2.)/(1.-x*x));
}

double xi_theta(const double &theta, const std::vector<int> &ells, const std::vector<double> &C_ells)
{
    // eqn (9) of Friedrich++ 2020 DES covariance modelling --- https://arxiv.org/pdf/2012.08568.pdf

    double xi = 0;
    double x = cos(theta); // theta should be in radians

    for (int l = 0; l <= ells.back() ; l++)
    {
        xi += (2.*l+1.) * gsl_sf_legendre_Pl(l, x) * C_ells.at(l);
    }

    return xi/(4.*M_PI);
}

double xip_theta(const double &theta, const std::vector<int> &ells, const std::vector<double> &C_ells)
{
    // eqn (9) of Friedrich++ 2020 DES covariance modelling --- https://arxiv.org/pdf/2012.08568.pdf

    double xip = 0;
    double x = cos(theta); // theta should be in radians

    for (int l = 0; l <= ells.back() ; l++)
    {
//        if (l > 150)
//            xip += (2*l+1) * G_ell_2_x_p(l,cos_theta)/(l*l*(l+1.)*(l+1.)) * C_ells.at(l);
//        else
//            continue;

        if (l > 2)
        {
            xip += (2.*l+1.) * G_ell_2_x_p(l,x)/(l*l*(l+1.)*(l+1.)) * C_ells.at(l);
        }
    }

    return xip*2./(4.*M_PI);
}

double xim_theta(const double &theta, const std::vector<int> &ells, const std::vector<double> &C_ells)
{
    // eqn (9) of Friedrich++ 2020 DES covariance modelling --- https://arxiv.org/pdf/2012.08568.pdf

    double xim = 0;
    double x = cos(theta); // theta should be in radians

    for (int l = 0; l <= ells.back() ; l++)
    {
        if (l > 2)
        {
            xim += (2.*l+1.) * G_ell_2_x_m(l,x)/(l*l*(l+1.)*(l+1.)) * C_ells.at(l);
        }
    }

    return xim*2./(4.*M_PI);
}

std::vector<double> xi_theta_array(std::vector<double> &alpha_angles_arcmins, std::vector<double> &l_array, std::vector<double> &P, const std::string &signal_type,
                                   int thread_count)
{
    std::vector<double> xi_array(alpha_angles_arcmins.size(), 0) ;

    std::vector<int> ells;
    std::vector<double> C_ells;
    C_ells_spherical_sky(l_array, P, ells, C_ells); // C_ells array for ells = 0,1,...,l_max-1

    #pragma omp parallel for num_threads(thread_count)
    for (size_t idx = 0; idx < alpha_angles_arcmins.size() ; idx++)
    {
        double theta = alpha_angles_arcmins.at(idx) * M_PI / 180.0 / 60.0;
        if (signal_type == "xi" || signal_type == "w")
            xi_array[idx] = xi_theta(theta, ells, C_ells);
        else if  (signal_type == "xip")
            xi_array[idx] = xip_theta(theta, ells, C_ells);
        else if  (signal_type == "xim")
            xi_array[idx] = xim_theta(theta, ells, C_ells);
    }

    return xi_array;
}

double P_ell_x_bin_averaged(const int &l, const double &x_min, const double &x_max, const double *P_l_x_min_array, const double *P_l_x_max_array)
{
    // eqn (65) of Friedrich++ 2020 DES covariance modelling --- https://arxiv.org/pdf/2012.08568.pdf

    return 1./(2.*l+1.)/(x_min - x_max)*(P_l_x_min_array[l+1]-P_l_x_min_array[l-1]-(P_l_x_max_array[l+1]-P_l_x_max_array[l-1]));
}

double G_ell_2_x_p_bin_averaged(const int &l, const double &x_min, const double &x_max,
                                const double *P_l_x_min_array, const double* ddx_P_l_x_min_array, const double *P_l_x_max_array, const double* ddx_P_l_x_max_array)
{
    // eqn (B5) of Friedrich++ 2020 DES covariance modelling --- https://arxiv.org/pdf/2012.08568.pdf

    return 1./(x_min - x_max)*(
            -l*(l-1.)/2.*(l+2./(2.*l+1.))*(P_l_x_min_array[l-1]-P_l_x_max_array[l-1])
            -l*(l-1.)*(2.-l)/2.*(x_min*P_l_x_min_array[l]-x_max*P_l_x_max_array[l])
            +l*(l-1.)/(2.*l+1.)*(P_l_x_min_array[l+1]-P_l_x_max_array[l+1])
            +(4.-l)*(ddx_P_l_x_min_array[l]-ddx_P_l_x_max_array[l])
            +(l+2.)*(x_min*ddx_P_l_x_min_array[l-1]-x_max*ddx_P_l_x_max_array[l-1]-(P_l_x_min_array[l-1]-P_l_x_max_array[l-1]))
            +2.*(l-1.)*(x_min*ddx_P_l_x_min_array[l]-x_max*ddx_P_l_x_max_array[l]-(P_l_x_min_array[l]-P_l_x_max_array[l]))
            -2.*(l+2.)*(ddx_P_l_x_min_array[l-1]-ddx_P_l_x_max_array[l-1])
             );

}

double G_ell_2_x_m_bin_averaged(const int &l, const double &x_min, const double &x_max,
                                const double *P_l_x_min_array, const double* ddx_P_l_x_min_array, const double *P_l_x_max_array, const double* ddx_P_l_x_max_array)
{
    // eqn (B5) of Friedrich++ 2020 DES covariance modelling --- https://arxiv.org/pdf/2012.08568.pdf

    return 1./(x_min - x_max)*(
            -l*(l-1.)/2.*(l+2./(2.*l+1.))*(P_l_x_min_array[l-1]-P_l_x_max_array[l-1])
            -l*(l-1.)*(2.-l)/2.*(x_min*P_l_x_min_array[l]-x_max*P_l_x_max_array[l])
            +l*(l-1.)/(2.*l+1.)*(P_l_x_min_array[l+1]-P_l_x_max_array[l+1])
            +(4.-l)*(ddx_P_l_x_min_array[l]-ddx_P_l_x_max_array[l])
            +(l+2.)*(x_min*ddx_P_l_x_min_array[l-1]-x_max*ddx_P_l_x_max_array[l-1]-(P_l_x_min_array[l-1]-P_l_x_max_array[l-1]))
            -2.*(l-1.)*(x_min*ddx_P_l_x_min_array[l]-x_max*ddx_P_l_x_max_array[l]-(P_l_x_min_array[l]-P_l_x_max_array[l]))
            +2.*(l+2.)*(ddx_P_l_x_min_array[l-1]-ddx_P_l_x_max_array[l-1])
             );

}

double xi_theta_bin_averaged(const double &theta_min, const double &theta_max, const std::vector<int> &ells, const std::vector<double> &C_ells)
{
    // eqn (9) of Friedrich++ 2020 DES covariance modelling replaced by bin averaged expression (B2) --- https://arxiv.org/pdf/2012.08568.pdf

    double xi = 0;

    double x_min = cos(theta_min); // theta_min should be in radians
    double x_max = cos(theta_max); // theta_max should be in radians

    int l_size = ells.back()+5;

    double P_l_x_min_array[l_size], ddx_P_l_x_min_array[l_size];
    gsl_sf_legendre_Pl_deriv_array(l_size, x_min, P_l_x_min_array, ddx_P_l_x_min_array);

    double P_l_x_max_array[l_size], ddx_P_l_x_max_array[l_size];
    gsl_sf_legendre_Pl_deriv_array(l_size, x_max, P_l_x_max_array, ddx_P_l_x_max_array);

    for (int l = 0; l <= ells.back() ; l++)
    {
        xi += (2.*l+1.) * P_ell_x_bin_averaged(l, x_min, x_max, P_l_x_min_array, P_l_x_max_array) * C_ells.at(l);
    }

    return xi/(4*M_PI);
}

double xip_theta_bin_averaged(const double &theta_min, const double &theta_max, const std::vector<int> &ells, const std::vector<double> &C_ells)
{
    // eqn (9) of Friedrich++ 2020 DES covariance modelling replaced by bin averaged expression (B5) --- https://arxiv.org/pdf/2012.08568.pdf

    double xip = 0.;

    double x_min = cos(theta_min); // theta_min should be in radians
    double x_max = cos(theta_max); // theta_max should be in radians

    int l_size = ells.back()+5;

    double P_l_x_min_array[l_size], ddx_P_l_x_min_array[l_size];
    gsl_sf_legendre_Pl_deriv_array(l_size, x_min, P_l_x_min_array, ddx_P_l_x_min_array);

    double P_l_x_max_array[l_size], ddx_P_l_x_max_array[l_size];
    gsl_sf_legendre_Pl_deriv_array(l_size, x_max, P_l_x_max_array, ddx_P_l_x_max_array);

    for (int l = 0; l <= ells.back() ; l++)
    {
        if (l > 2)
            xip += (2.*l+1.) * G_ell_2_x_p_bin_averaged(l, x_min, x_max, P_l_x_min_array, ddx_P_l_x_min_array,
                                                          P_l_x_max_array, ddx_P_l_x_max_array) / (l*l*(l+1.)*(l+1.)) * C_ells.at(l);
    }

    return xip*2./(4.*M_PI);
}

double xim_theta_bin_averaged(const double &theta_min, const double &theta_max, const std::vector<int> &ells, const std::vector<double> &C_ells)
{
    // eqn (9) of Friedrich++ 2020 DES covariance modelling replaced by bin averaged expression (B5) --- https://arxiv.org/pdf/2012.08568.pdf

    double xim = 0.;

    double x_min = cos(theta_min); // theta_min should be in radians
    double x_max = cos(theta_max); // theta_max should be in radians

    int l_max = ells.back()+1;

    double P_l_x_min_array[l_max], ddx_P_l_x_min_array[l_max];
    gsl_sf_legendre_Pl_deriv_array(l_max, x_min, P_l_x_min_array, ddx_P_l_x_min_array);

    double P_l_x_max_array[l_max], ddx_P_l_x_max_array[l_max];
    gsl_sf_legendre_Pl_deriv_array(l_max, x_max, P_l_x_max_array, ddx_P_l_x_max_array);

    for (int l = 0; l <= ells.back() ; l++)
    {
        if (l > 2)
            xim += (2.*l+1.) * G_ell_2_x_m_bin_averaged(l, x_min, x_max, P_l_x_min_array, ddx_P_l_x_min_array,
                                                          P_l_x_max_array, ddx_P_l_x_max_array) / (l*l*(l+1.)*(l+1.)) * C_ells.at(l);
    }

    return xim*2./(4.*M_PI);
}

std::vector<double> xi_theta_array_bin_averaged(const std::vector<double> &alpha_min_angles_arcmins, const std::vector<double> &alpha_max_angles_arcmins,
                                                const std::vector<double> &l_array, const std::vector<double> &P, const std::string &signal_type, int thread_count)
{
    std::vector<double> xi_array(alpha_min_angles_arcmins.size(), 0) ;

    std::vector<int> ells;
    std::vector<double> C_ells;
    C_ells_spherical_sky(l_array, P, ells, C_ells); // C_ells array for ells = 0,1,...,l_max-1

    #pragma omp parallel for num_threads(thread_count)
    for (size_t idx = 0; idx < alpha_min_angles_arcmins.size() ; idx++)
    {
        double theta_min = alpha_min_angles_arcmins.at(idx) * M_PI / 180.0 / 60.0;
        double theta_max = alpha_max_angles_arcmins.at(idx) * M_PI / 180.0 / 60.0;

        if (signal_type == "xi" || signal_type == "w")
            xi_array[idx] = xi_theta_bin_averaged(theta_min, theta_max, ells, C_ells);
        else if  (signal_type == "xip")
            xi_array[idx] = xip_theta_bin_averaged(theta_min, theta_max, ells, C_ells);
        else if  (signal_type == "xim")
            xi_array[idx] = xim_theta_bin_averaged(theta_min, theta_max, ells, C_ells);
    }

    return xi_array;
}
