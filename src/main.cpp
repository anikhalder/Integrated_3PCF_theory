#include <iostream>
#include <experimental/filesystem>
#include <omp.h>
#include <fstream>
#include <memory>
#include <constants.h>
#include <cosmology_utils.h>
#include <power_spectra_2D.h>
#include <bispectra_2D.h>
#include <power_spectrum.hpp>
#include <bispectrum.hpp>
#include <integrated_bispectra_2D.h>
#include <sys/time.h>
#include <assert.h>
#include <real_space_2D.h>
#include <sstream>
#include <halo_utils.hpp>
#include <VEGAS_Integrator.h>

/*
 * Things to check/change before compiling the code and performing an execution:
 *
 * constants.h:
 * - compute_bihalofit keyword
 * - apply_T17_corrections keyword
 *
 * bispectrum.cpp:
 * - bispectrum type (tree, GM, SC, bihalofit, GM+RF etc.)
 *
 * integrated_bispectra_2D.cpp:
 * - in the iB_z_l_1_l_2_phi_1_phi_2_mc() or iB_z_l_1_l_2_phi_1_phi_2_hcubature() check the phi_l value
 * - also check whether you are using X or X_v2 form of the window functions
 * 
 * main.cpp:
 * - thread_count
 * - index i
 * - k_max
 * - non linear halofit or hmcode
 * - use_pk_nl keyword
 * - spectra_folder and correlations_folder names
 * - l_array
 * - zs_bins
 * - filename_P and filename_iB names
 * - calls_iB_initial
 * - bin averaging or not for computing correlation functions
 *
 */

/*
double func_weight(std::vector<double> x, void* param)
{
    double dx = *((double *)param);
    double xmin = -1.0 + dx;
    double xmax = 1.0 - dx;
    double x_true = x[0]*2.0-1.0;
    if (x_true < xmin || x_true > xmax)
    {
        return 0;
    }
    return (1.0+x_true*x_true)/(1.0-x_true*x_true)*2;
}
double func_WW2hh(std::vector<double> x, void *param)
{
    double *par = (double *)param;
    double s = par[0]*par[0];
    double MH = 125.0;
    double MW = 80.385;
    double MH2 = MH*MH;
    double MW2 = MW*MW;
    double cth = x[0]*2.0-1.0;
    double cth2 = cth*cth;
    double num = 8*pow(MH2,3)*(cth2*s-2*MW2+s)+2*MH2*MH2*(16*(2*cth2+1)*MW2*MW2-40*cth2*MW2*s+(cth2-3)*s*s) + MH2*s*(16*(cth2-3)*MW2*MW2 + 4*(7*cth2+4)*MW2*s-(cth2-1)*s*s) - 8*(cth2-2)*MW2*MW2*s*s - 2*(cth2+3)*MW2*s*s*s;
    double den = 16*MW2*MW2*pow(s-MH2,2)*pow(-4*MH2*(cth2*(4*MW2-s)+s)+s*(cth2*(4*MW2-s)+s)+4*MH2*MH2,2);
    return pow(num,2)/den;
}
double func_4D(std::vector<double> x, void *param)
{
    double r1[4] = {0.33,0.5,0.5,0.5};
    double r2[4] = {0.67,0.5,0.5,0.5};
    double x1=0;
    double x2=0;
    for (int i = 0; i < 4; i++)
    {
        x1 += -100*pow(x[i]-r1[i],2);
        x2 += -100*pow(x[i]-r2[i],2);
    }
    return exp(x1) + exp(x2);
}
double func_2D(std::vector<double> x, void *param)
{
    double xl = -2;
    double xu = 3;
    double yl = -5;
    double yu = 15;
    
    return (pow(xl + (xu-xl)*x[0],2) + pow(yl + (yu-yl)*x[1],2))*(xu-xl)*(yu-yl);
}
*/

int main()
{  
    struct timeval start_file, end_file;
    gettimeofday(&start_file, nullptr);
    double time_taken_file;

    #pragma omp parallel
    {
        #pragma omp master
        {
            std::cout << "A maximum of " << omp_get_max_threads() << " OpenMP threads are present in this machine!" << std::endl;
        }
    }

    //const int thread_count = static_cast<int>(omp_get_max_threads())-1;
    const int thread_count = 120;
    
    if (verbose_print_outs)
        std::cout<<"Number of threads that will be used if parallelisation is requested = " << thread_count << std::endl;

    VEGAS_Integrator inter;
    inter.Set_Verbose(NONE);
    
    /*
    #pragma omp parallel for num_threads(thread_count)
    for (double dx = 0.02; dx < 0.31; dx += 0.02)
    {
        VEGAS_Integrator inter;
        inter.Set_Verbose(NONE);
        inter.Set_Integrand(func_weight, 1, &dx);
        inter.Improve_Grid();
        inter.Integration();
        std::cout<<"dx: "<<dx<<" res: "<<inter.Get_Result()<<" err: "<<inter.Get_Error()<<" chi2: "<<inter.Get_Chisq()<<std::endl;
    }
    */
    
    /*
    double energies[37] = {1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,5000,6000,7000,8000,9000,10000,12000,14000,16000,18000,20000,22000,24000,26000,28000,30000};
    //#pragma omp parallel for num_threads(thread_count)
    for (int i = 0; i < 37; i++)
    {
        VEGAS_Integrator inter;
        inter.Set_Verbose(NONE);
        inter.Set_Integrand(func_WW2hh, 1, &energies[i]);
        inter.Improve_Grid();
        inter.Integration();
        std::cout<<"Ecm: "<<energies[i]<<" res: "<<inter.Get_Result()/pow(energies[i],2)<<" err: "<<inter.Get_Error()/pow(energies[i],2)<<" chi2: "<<inter.Get_Chisq()<<std::endl;
    }
    
    inter.Set_Integrand(func_4D,4,NULL);
    inter.Improve_Grid();
    inter.Integration();
    std::cout<<"Result: "<<inter.Get_Result()<<" Error: "<<inter.Get_Error()<<" chi2: "<<inter.Get_Chisq()<<std::endl;

    inter.Set_Integrand(func_2D,2,NULL);
    inter.Improve_Grid();
    inter.Integration();
    std::cout<<"Result: "<<inter.Get_Result()<<" Error: "<<inter.Get_Error()<<" chi2: "<<inter.Get_Chisq()<<std::endl;
    */

    std::vector<std::string> filename_extension_array =
            { ".dat",
              "_minus_2step_Omega_cdm.dat", "_minus_1step_Omega_cdm.dat", "_plus_1step_Omega_cdm.dat", "_plus_2step_Omega_cdm.dat",
              "_minus_2step_sigma8.dat", "_minus_1step_sigma8.dat", "_plus_1step_sigma8.dat", "_plus_2step_sigma8.dat",
              "_minus_2step_n_s.dat", "_minus_1step_n_s.dat", "_plus_1step_n_s.dat", "_plus_2step_n_s.dat",
              "_minus_2step_w_0.dat", "_minus_1step_w_0.dat", "_plus_1step_w_0.dat", "_plus_2step_w_0.dat",
              "_minus_2step_w_a.dat", "_minus_1step_w_a.dat", "_plus_1step_w_a.dat", "_plus_2step_w_a.dat",
              "_minus_2step_eta_0.dat", "_minus_1step_eta_0.dat", "_plus_1step_eta_0.dat", "_plus_2step_eta_0.dat",
              "_minus_2step_c_min.dat", "_minus_1step_c_min.dat", "_plus_1step_c_min.dat", "_plus_2step_c_min.dat",
              "_minus_2step_h.dat", "_minus_1step_h.dat", "_plus_1step_h.dat", "_plus_2step_h.dat"};

    //for (size_t i=0; i<filename_extension_array.size(); i++)
    for (size_t i=0; i<1; i++)
    {
        // For normal runs only use the first iteration of the loop i.e. i=0 --> '.dat'
        // For Fisher forecast computations use the corresponding desired range of i as per the filename_extension_array above

        std::string filename_extension = filename_extension_array.at(i);

        // -------------------------------------------------------------------------------------

        struct timeval start, end;

        double time_taken;

        const gsl_rng_type *T = gsl_rng_default;

        size_t counter = 0;

        // ######################################################################################
        // ######################################################################################
        // ######################################################################################

        // Setting up cosmology and CLASS object

        gettimeofday(&start, nullptr);

        // -------------------------------------------------------------------------------------

        // Set cosmology (only comment out/enter the one that is needed)

        // Buzzard simulations cosmology

        //double omega_b = 0.02303; // Omega_b*h*h
        //double omega_cdm = 0.11711; // Omega_cdm*h*h
        //double h = 0.7;
        //double sigma8 = 0.82;
        //double n_s = 0.96;

        // MICE simulations cosmology

        double omega_b = 0.02156; // Omega_b*h*h
        double omega_cdm = 0.10094; // Omega_cdm*h*h
        double h = 0.7;
        double sigma8 = 0.8;
        double n_s = 0.95;

        // Planck 2015 (e.g. used in bihalofit)

        //double h=0.6727;
        //double omega_b = 0.0492*h*h;
        //double omega_cdm = 0.2664*h*h;
        //double sigma8 = 0.831;
        //double n_s = 0.9645;

        // Planck 2015 Colossus

        //double h=0.6774;
        //double omega_b = 0.0486*h*h;
        //double omega_cdm = (0.3089-0.0486)*h*h;
        //double sigma8 = 0.8159;
        //double n_s = 0.9667;

        // KiDS - Pierre Bürger cosmology

        //double h=0.6898;
        //double omega_b = 0.047*h*h;
        //double omega_cdm = 0.2905*h*h-omega_b;
        //double sigma8 = 0.813;
        //double n_s = 0.969;

        // MassiveNus simulations cosmology

        // double omega_b = 0.0223; // Omega_b*h*h
        // double omega_cdm = 0.1247; // Omega_cdm*h*h
        // double h = 0.7;
        // double A_s = 2.1e-9;
        // double n_s = 0.97;

        // TODO: add equations for massive neutrinos and set the CLASS massive neutrino parameters!!!!

        // Takahashi simulations cosmology

        ////double omega_b = 0.02254; // Omega_b*h*h i.e. the physical baryon density 
        ////double omega_cdm = 0.11417; // Omega_cdm*h*h i.e. the physical cdm density 
        // double Omega_b = 0.046;
        // double Omega_cdm = 0.233;
        // double h = 0.7;
        // double sigma8 = 0.82;
        // double n_s = 0.97;

        double w_0 = -1.0;
        double w_a = 0.0;
        double eta_0 = 0.603; // emu_dmonly (fiducial)
        double c_min = 3.13; // emu_dmonly (fiducial)

        //double eta_0 = 0.64; // owls_dmonly
        //double c_min = 3.43; // owls_dmonly

        //double eta_0 = 0.76; // owls_agn
        //double c_min = 2.32; // owls_agn

        double Omega_b = omega_b/h/h;
        double Omega_cdm = omega_cdm/h/h;

        // Useful only for Fisher forecasting i.e. i>0
//        if (i==1)
//            omega_cdm *= 0.92;
//        else if (i==2)
//            omega_cdm *= 0.96;
//        else if (i==3)
//            omega_cdm *= 1.04;
//        else if (i==4)
//            omega_cdm *= 1.08;
        if (i==1)
            Omega_cdm *= 0.92;
        else if (i==2)
            Omega_cdm *= 0.96;
        else if (i==3)
            Omega_cdm *= 1.04;
        else if (i==4)
            Omega_cdm *= 1.08;
        else if (i==5)
            sigma8 *= 0.96;
        else if (i==6)
            sigma8 *= 0.98;
        else if (i==7)
            sigma8 *= 1.02;
        else if (i==8)
            sigma8 *= 1.04;
        else if (i==9)
            n_s *= 0.8;
        else if (i==10)
            n_s *= 0.9;
        else if (i==11)
            n_s *= 1.1;
        else if (i==12)
            n_s *= 1.2;
        else if (i==13)
            w_0 += -0.16;
        else if (i==14)
            w_0 += -0.08;
        else if (i==15)
            w_0 += 0.08;
        else if (i==16)
            w_0 += 0.16;
        else if (i==17)
            w_a += -0.32;
        else if (i==18)
            w_a += -0.16;
        else if (i==19)
            w_a += 0.16;
        else if (i==20)
            w_a += 0.32;
        else if (i==21)
            eta_0 *= 0.9;
        else if (i==22)
            eta_0 *= 0.95;
        else if (i==23)
            eta_0 *= 1.05;
        else if (i==24)
            eta_0 *= 1.1;
        else if (i==25)
            c_min *= 0.9;
        else if (i==26)
            c_min *= 0.95;
        else if (i==27)
            c_min *= 1.05;
        else if (i==28)
            c_min *= 1.1;
        else if (i==29)
            h *= 0.96;
        else if (i==30)
            h *= 0.98;
        else if (i==31)
            h *= 1.02;
        else if (i==32)
            h *= 1.04;

        //double omega_b = Omega_b*h*h;
        //double omega_cdm = Omega_cdm*h*h;

        // -------------------------------------------------------------------------------------

        // CLASS config

        ClassParams pars;

        // -------------------------

        // set cosmological parameters within CLASS

        pars.add("h",h);
        pars.add("omega_b",omega_b); // Omega_b*h*h
        pars.add("omega_cdm",omega_cdm); // Omega_cdm*h*h
        pars.add("n_s",n_s);

        // use either sigma8 or A_s (comment out the one not being used)
        pars.add("sigma8",sigma8);
        //pars.add("A_s",A_s);

        // for dark energy
        // pars.add("Omega_Lambda",0.0);
        // pars.add("fluid_equation_of_state","CLP");
        // pars.add("w0_fld",w_0);
        // pars.add("wa_fld",w_a);

        // TODO: for massive neutrinos !!!!

        pars.add("N_ur",3.046);
        pars.add("N_ncdm",0.0);
        pars.add("Omega_k",0.0);
        pars.add("Omega_fld",0.0);
        pars.add("Omega_scf",0.0);
        pars.add("YHe",0.24);
        

        // -------------------------

        // for the computational output of the 3D matter power spectrum, total density and velocity transfer functions

        //pars.add("output","mPk, dTk, vTk"); // settings for papers
        pars.add("output","mPk");

        // -------------------------

        // for non-linear power spectrum use either halofit or hmcode (comment out the one not being used)

        pars.add("non linear","halofit");
        //pars.add("pk_eq","yes"); // can use this when using Halofit and Omega_fld != 0 & w_a != 0 (but this is a 'maybe' as suggested in class)

        //pars.add("non linear","hmcode");

        // -------------------------

        // baryonic parameters (only when using hmcode 2016; otherwise comment these lines out)

        //pars.add("eta_0",eta_0);
        //pars.add("c_min",c_min);

        // -------------------------

        // set k_max, z_max for 3D Pk computation by CLASS

        //pars.add("P_k_max_1/Mpc",10.0); // this is good for quick tests
        pars.add("P_k_max_1/Mpc",30.0); // this is good for quick tests
        //pars.add("P_k_max_1/Mpc",150.0); // this is good for quick tests

        //pars.add("P_k_max_1/Mpc",5000.0); // for l1 + l2 = 20000 this is good enough (for lowest z=0.001)
        //pars.add("P_k_max_1/Mpc",6000.0); // for l1 + l2 = 25000 this is good enough (for lowest z=0.001)
        //pars.add("P_k_max_h/Mpc",15853.0); // previous settings

        //pars.add("P_k_max_1/Mpc",3000.0); // for l1 + l2 = 50000 this is good enough (for lowest z=0.005) --> current settings for paper

        pars.add("z_max_pk",3.5);

        // -------------------------

        // set verbose parameters for showing output to the terminal

        pars.add("input_verbose",0);
        pars.add("background_verbose",0);
        pars.add("thermodynamics_verbose",0);
        pars.add("perturbations_verbose",0);
        pars.add("transfer_verbose",0);
        pars.add("primordial_verbose",0);
        pars.add("spectra_verbose",0);
        pars.add("nonlinear_verbose",0);
        pars.add("lensing_verbose",0);
        pars.add("output_verbose",0);

        // create CLASS object using parameters input in main.cpp
        std::unique_ptr<ClassEngine> class_obj (new ClassEngine(pars, false));

        // create CLASS object using parameters input in main.cpp along with .pre file
        //std::unique_ptr<ClassEngine> class_obj (new ClassEngine(pars, (char*)"pk_ref.pre", false));

        // create CLASS object using .ini file (can also provide .pre file)
        //char *argv[] = {(char*)"class", (char*)"takahashi_fiducial.ini", NULL};
        //int argc = 0; 
        //while(argv[++argc] != NULL);
        //std::unique_ptr<ClassEngine> class_obj (new ClassEngine(argc, argv));

        // -------------------------

        // TODO: set primordial non-Gaussianity (not really needed as of now and will investigate this in the future --> hence, setting them to zero)

        // class_obj->set_f_NL_local(0);
        // class_obj->set_f_NL_equilateral(0);
        // class_obj->set_f_NL_orthogonal(0);

        gettimeofday(&end, nullptr);
        time_taken = (end.tv_sec - start.tv_sec) * 1e6;
        time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
        std::cout << "\nTime taken for setting cosmology and CLASS object creation: " << time_taken << " sec" << std::endl;

        // -------------------------------------------------------------------------------------

        if (verbose_print_outs)
        {
            std::cout << "\n#########################################" << std::endl;
            std::cout << "#" << i << std::endl;
            std::cout << "Omega_cdm   sigma8   n_s   w_0   w_a   eta_0   c_min   h" << std::endl;
            std::cout << Omega_cdm << " " << sigma8 << " " << n_s << " " << w_0 << " " << w_a << " " << eta_0 << " " << c_min << " " << h << "\n" << std::endl;
        }

        // ######################################################################################
        // ######################################################################################
        // ######################################################################################

        // Computation settings

        bool compute_initial_checks = false;
        bool compute_sigma_quantities_tables = false;
        bool compute_Pk_shell_integration_table = false;
        bool compute_Pk_Bkkk_tables = false;
        bool compute_n_eff_table = false;
        bool compute_k_NL_table = false;
        bool compute_response_function_tables = false;
        bool compute_transfer_function_tables = false;
        bool compute_halo_quantities_tables = false;

        // bool compute_2D_integrated_3PCF_area_pre_factors = true;
        // bool compute_2D_power_spectra = true;
        // bool compute_2D_power_spectra_spherical_sky = false; // e.g. needed for FLASK (this can only be computed when compute_2D_power_spectra = true)
        // bool compute_2D_2PCF = true;
        // bool compute_2D_bispectra_equilateral = false;
        // bool compute_2D_integrated_bispectra = false; // OLD to be deleted
        // bool compute_2D_integrated_bispectra_v2 = true;
        // bool compute_2D_integrated_3PCF = true;

        bool compute_2D_integrated_3PCF_area_pre_factors = true;
        bool compute_2D_power_spectra = true;
        bool compute_2D_power_spectra_spherical_sky = false; // e.g. needed for FLASK (this can only be computed when compute_2D_power_spectra = true)
        bool compute_2D_2PCF = true;
        bool compute_2D_bispectra_equilateral = false;
        bool compute_2D_integrated_bispectra = false; // OLD to be deleted
        bool compute_2D_integrated_bispectra_v2 = false;
        bool compute_2D_integrated_bispectra_v3 = true;
        bool compute_2D_integrated_3PCF = true;

        // ######################################################################################
        // ######################################################################################
        // ######################################################################################

        // Folder settings

        //std::string spectra_folder = "./takahashi_B2D_GM_157_iBkxi_U70W75W75_cross_zs10_zs16_1e7_20000/";
        //std::string correlations_folder = "./takahashi_B2D_GM_157_iZkxi_U70W75W75_cross_zs10_zs16_1e7_20000/";

        //std::string spectra_folder = "./takahashi_B2D_GM_157_iBkxi_U70W75W75_cross_zs10_zs16_1e7_20000_with_hmcode_correct_change_params/";
        //std::string correlations_folder = "./takahashi_B2D_GM_157_iZkxi_U70W75W75_cross_zs10_zs16_1e7_20000_with_hmcode_correct_change_params/";

        //std::string spectra_folder = "./takahashi_B2D_GM_squeezed_RF_stitched_157_iBkxi_U70W75W75_cross_zs10_zs16_1e7_20000_squeezed_9_with_hmcode_correct_change_params/";
        //std::string correlations_folder = "./takahashi_B2D_GM_squeezed_RF_stitched_157_iZkxi_U70W75W75_cross_zs10_zs16_1e7_20000_squeezed_9_with_hmcode_correct_change_params/";

        //std::string spectra_folder = "./takahashi_B2D_GM_squeezed_RF_stitched_157_iBkxi_U70W75W75_zs10_mc_2e6_x2_150_20000_squeezed_9_with_baryons_change_params_nonsq_GM_sq_RF/";
        //std::string correlations_folder = "./takahashi_B2D_GM_squeezed_RF_stitched_157_iZkxi_U70W75W75_zs10_mc_2e6_x2_150_20000_squeezed_9_with_baryons_change_params_nonsq_GM_sq_RF/";

        // -------------------------
        // With 157 ells upto ell_max=20000; MC integration (2e6 and double of it for ell<150)

        // For validation against T17 with bsr corrections
        //std::string spectra_folder = "./takahashi_bsr_nonsq_GM_sq9_RF_ell157_iBkxi_U70W75W75_cross_zs10_zs16_mc_2e6_x2_150_20000/";
        //std::string correlations_folder = "./takahashi_bsr_nonsq_GM_sq9_RF_ell157_iZkxi_U70W75W75_cross_zs10_zs16_mc_2e6_x2_150_20000/";

        // For Fisher forecast
        //std::string spectra_folder = "./takahashi_nonsq_GM_sq_RF_ell157_iBkxi_U70W75W75_cross_zs10_zs16_mc_2e6_x2_150_20000_with_baryons_change_params/";
        //std::string correlations_folder = "./takahashi_nonsq_GM_sq_RF_ell157_iZkxi_U70W75W75_cross_zs10_zs16_mc_2e6_x2_150_20000_with_baryons_change_params/";

        // -------------------------
        // With 150 ells upto ell_max=12000; MC integration (2e6 and double of it for ell<150)

        // For validation against T17 with bsr corrections - GM / tree / bihalofit
        //std::string spectra_folder = "./takahashi_bsr_bihalofit_ell150_iBkxi_U70W75W75_cross_zs10_zs16_mc_2e6_x2_150_12000/";
        //std::string correlations_folder = "./takahashi_bsr_bihalofit_ell150_iZkxi_U70W75W75_cross_zs10_zs16_mc_2e6_x2_150_12000/";

        // For validation against T17 with bsr corrections - nonsq GM / tree + sq RF
        //std::string spectra_folder = "./takahashi_bsr_nonsq_GM_sq6_RF_ell150_iBkxi_U70W75W75_cross_zs10_zs16_mc_2e6_x2_150_12000/";
        //std::string correlations_folder = "./takahashi_bsr_nonsq_GM_sq6_RF_ell150_iZkxi_U70W75W75_cross_zs10_zs16_mc_2e6_x2_150_12000/";

        // For Fisher forecast - GM+RF
        //std::string spectra_folder = "./takahashi_nonsq_GM_sq9_RF_ell150_iBkxi_U70W75W75_cross_zs10_zs16_mc_2e6_x2_150_12000_with_baryons_change_params/";
        //std::string correlations_folder = "./takahashi_nonsq_GM_sq9_RF_ell150_iZkxi_U70W75W75_cross_zs10_zs16_mc_2e6_x2_150_12000_with_baryons_change_params/";

        // -------------------------
        // With 150 ells upto ell_max=20000; MC integration (1e6 and double of it for ell<=200)

        // For validation against T17 with bsr corrections - GM / tree / bihalofit
        //std::string spectra_folder = "./takahashi_bsr_bihalofit_ell150_iB_Mss_U70W75W75_cross_zs10_zs16_mc_1e6_x2_200_20000/";
        //std::string correlations_folder = "./takahashi_bsr_bihalofit_ell150_iZ_Mss_U70W75W75_cross_zs10_zs16_mc_1e6_x2_200_20000/";

        // For validation against T17 with bsr corrections - nonsq GM / tree + sq RF
        //std::string spectra_folder = "./takahashi_bsr_nonsq_GM_sq7_RF_ell150_iB_Mss_U70W75W75_cross_zs10_zs16_mc_1e6_x2_200_20000/";
        //std::string correlations_folder = "./takahashi_bsr_nonsq_GM_sq7_RF_ell150_iZ_Mss_U70W75W75_cross_zs10_zs16_mc_1e6_x2_200_20000/";

        // For Fisher forecast - GM+RF
        //*std::string spectra_folder = "./takahashi_nonsq_GM_sq7_RF_ell150_iB_Mss_U70W75W75_cross_zs10_zs16_mc_2.5e6_x2_200_20000_with_baryons_change_params/";
        //*std::string correlations_folder = "./takahashi_nonsq_GM_sq7_RF_ell150_iZ_Mss_U70W75W75_cross_zs10_zs16_mc_2.5e6_x2_200_20000_with_baryons_change_params/";

        // -------------------------
        // With 120 ells upto ell_max=20000; MC integration (2e6 and double of it for ell<=220)

//        std::string spectra_folder = "./takahashi_nonsq_GM_sq7_RF_ell120_iB_Mss_U70W75W75_cross_zs10_zs16_mc_2e6_x2_220_20000_owls_dmonly/";
//        std::string correlations_folder = "./takahashi_nonsq_GM_sq7_RF_ell120_iZ_Mss_U70W75W75_cross_zs10_zs16_mc_2e6_x2_220_20000_owls_dmonly/";
//        std::string correlations_folder = "./takahashi_nonsq_GM_sq7_RF_ell120_iZ_Mss_U70W75W75_cross_zs10_zs16_mc_2e6_x2_220_20000_owls_dmonly_bin_averaged/";

//        std::string spectra_folder = "./takahashi_nonsq_GM_sq7_RF_ell120_iB_Mss_U70W75W75_cross_zs10_zs16_mc_2e6_x2_220_20000_owls_agn/";
//        std::string correlations_folder = "./takahashi_nonsq_GM_sq7_RF_ell120_iZ_Mss_U70W75W75_cross_zs10_zs16_mc_2e6_x2_220_20000_owls_agn/";
//        std::string correlations_folder = "./takahashi_nonsq_GM_sq7_RF_ell120_iZ_Mss_U70W75W75_cross_zs10_zs16_mc_2e6_x2_220_20000_owls_agn_bin_averaged/";

        // For Fisher forecast - GM+RF
//        std::string spectra_folder = "./takahashi_nonsq_GM_sq7_RF_ell120_iB_Mss_U70W75W75_cross_zs10_zs16_mc_2e6_x2_220_20000_with_baryons_change_params_P_qag/";
//        //std::string correlations_folder = "./takahashi_nonsq_GM_sq7_RF_ell120_iZ_Mss_U70W75W75_cross_zs10_zs16_mc_2e6_x2_220_20000_with_baryons_change_params/";
//        std::string correlations_folder = "./takahashi_nonsq_GM_sq7_RF_ell120_iZ_Mss_U70W75W75_cross_zs10_zs16_mc_2e6_x2_220_20000_with_baryons_change_params_bin_averaged_P_qag/";

        // -------------------------
        // Halo bias
        // With 150 ells upto ell_max=20000; MC integration (1e6 and double of it for ell<=200)

        // For validation against T17 with bsr corrections - GM / tree / bihalofit
        //std::string spectra_folder = "./takahashi_tree_ell150_iB_hhh_W75W75W75_sample_2_mc_1e6_x2_200_20000_b1_b2_Tinker_bs2_coevolution_sn_HMF_Tinker/";
        //std::string correlations_folder = "./takahashi_tree_ell150_iZ_hhh_W75W75W75_sample_2_mc_1e6_x2_200_20000_b1_b2_Tinker_bs2_coevolution_sn_HMF_Tinker/";

        // With 120 ells upto ell_max=20000; MC integration (1e6 and double of it for ell<=220)

        // For validation against T17 with bsr corrections - GM / tree / bihalofit
        //std::string spectra_folder = "./takahashi_nl_ell120_iB_hhh_W75W75W75_redmagicLens1_mc_1e6_x2_220_20000_b1_b2_Tinker_bs2_coevolution_HMF_Tinker/";
        //std::string correlations_folder = "./takahashi_nl_ell120_iZ_hhh_W75W75W75_redmagicLens1_mc_1e6_x2_220_20000_b1_b2_Tinker_bs2_coevolution_HMF_Tinker/";

        // ########################
        // Test
        std::string spectra_folder = "./mice_tree_ell120_iB_kkk_W75W75W75_zs993_mc_1e6_20000_kmax_30_Mpc_Laurence_settings_v3/";
        std::string correlations_folder = "./mice_tree_ell120_iZ_kkk_W75W75W75_zs993_mc_1e6_20000_kmax_30_Mpc_Laurence_settings_v3/";

        //std::string spectra_folder = "./takahashi_bsr_nonsq_GM_sq7_RF_ell120_iB_Mss_U70W75W75_cross_zs10_zs16_mc_2e6_x2_220_20000_X_v2_bin_averaged_phi_l_zero/";
        //std::string correlations_folder = "./takahashi_bsr_nonsq_GM_sq7_RF_ell120_iZ_Mss_U70W75W75_cross_zs10_zs16_mc_2e6_x2_220_20000_X_v2_bin_averaged_phi_l_zero/";

        //std::string spectra_folder = "./takahashi_bsr_nonsq_tree_sq7_tree_ell120_iB_Mss_U70W75W75_cross_zs10_zs16_mc_cigar_20000_bin_averaged/";
        //std::string correlations_folder = "./takahashi_bsr_nonsq_tree_sq7_tree_ell120_iZ_Mss_U70W75W75_cross_zs10_zs16_mc_cigar_20000_bin_averaged/";

        // -------------------------------------------------------------------------------------

        // create folders

        std::experimental::filesystem::create_directory(spectra_folder);
        std::experimental::filesystem::create_directory(correlations_folder);

        // ######################################################################################
        // ######################################################################################
        // ######################################################################################

        // 2D window, projection kernel, spectra and correlation calculation settings

        double theta_U = 70*arcmin;

        //double patch_area_sq_degrees = 5.0;
        //double theta_T = spherical_cap_sqdeg_2_radius(patch_area_sq_degrees);

        double theta_T = 75*arcmin;
        double patch_area_sq_radians = spherical_cap_radius_2_sqradians(theta_T);

        if (verbose_print_outs)
            std::cout << "\nPatch area (square radians) = " << patch_area_sq_radians << std::endl;

        std::vector<double> alpha_table = read_1_column_table("../data/angular_bins/alpha_angles_arcmins_20_bins.tab");
        assert(!alpha_table.empty());

        std::vector<double> alpha_min_table = read_1_column_table("../data/angular_bins/alpha_min_angles_arcmins_20_bins.tab");
        assert(!alpha_min_table.empty());

        std::vector<double> alpha_max_table = read_1_column_table("../data/angular_bins/alpha_max_angles_arcmins_20_bins.tab");
        assert(!alpha_max_table.empty());

        std::vector<double> l_array;

        //for(int l=1; l<=10001; l++)
        //    l_array.push_back(l);

        //l_array = read_1_column_table("../data/ell_arrays/ell_array_157.tab"); // ell_max = 20000 where 157 points log-spaced (1-20000) -- previous settings
        //l_array = read_1_column_table("../data/ell_arrays/ell_array_150_disjoint.tab"); // ell_max = 20000 where 30 points log-spaced (1-200); 120 points log-spaced (230-20000)

        l_array = read_1_column_table("../data/ell_arrays/ell_array_120.tab"); // ell_max = 20000 where 120 points log-spaced (1-20000) -- current settings for papers

        //double a = log10(2);
        //double b = log10(20000);
        //double num_pts = 30;

        //for(int i=0; i<num_pts; i++)
        //    l_array.push_back(pow(10, a + i * (b - a) / (num_pts - 1)));

        if (verbose_print_outs)
            std::cout << "\nSize of l_array = " << l_array.size() << "\n" << std::endl;

        // --------------------------------------------------------
        // Source (shear/convergence) correlations settings
        // --------------------------------------------------------

        double zs1 = 0.993;  // MICE
        //double zs1 = z_cmb; // CMB

        double zs2 = 0.993;  // MICE

        //double zs1 =  0.5739; // Takahashi zs10
        //double zs2 =  1.0334; // Takahashi zs16

        //double zs1 =  1.0;
        //double zs2 =  1.5;

        //double zs1 = 2.0; // DESY3
        //double zs2 = 2.0; // DESY3

        //double zs1 = 0.5; // MassiveNus
        //double zs2 = 1.0; // MassiveNus

        double zs_lower = 0;
        double zs_upper = zs2;

        // redshift bins --> make sure that this array is in ascending order for multiple bins
        //std::vector<double> zs_bins{zs1, zs2}; // source (shear/convergence) redshift bin median --> for shear only

        std::vector<double> zs_bins{zs2}; // source (shear/convergence) redshift bin median --> for halos

        size_t num_2pt_ss_correlations = num_correlations(zs_bins.size(), 2);
        size_t num_3pt_sss_correlations = num_correlations(zs_bins.size(), 3);
        size_t num_i3pt_sss_correlations = num_correlations(zs_bins.size(), 3);

        bool do_cross_fields = true; // previous settings
        //bool do_cross_fields = false;

        std::shared_ptr<projection_kernel> qk1(new projection_kernel_q_k_zs_fixed(class_obj.get(), zs1));
        std::shared_ptr<projection_kernel> qk2(new projection_kernel_q_k_zs_fixed(class_obj.get(), zs2));

        //std::vector<std::vector<double>> nofz_s1_table = read_2_column_table("../data/nofz/DESY3_nofz/nofz_DESY3_source_BIN2.tab");
        //assert(!nofz_s1_table.empty());
        //normalise_nofz(nofz_s1_table);
        //Linear_interp_1D nofz_s1(nofz_s1_table.at(0), nofz_s1_table.at(1));

        //std::vector<std::vector<double>> nofz_s2_table = read_2_column_table("../data/nofz/DESY3_nofz/nofz_DESY3_source_BIN4.tab");
        //assert(!nofz_s2_table.empty());
        //normalise_nofz(nofz_s2_table);
        //Linear_interp_1D nofz_s2(nofz_s2_table.at(0), nofz_s2_table.at(1));

        //std::shared_ptr<projection_kernel> qk1(new projection_kernel_q_k_zs_distribution(class_obj.get(), &nofz_s1, zs1));
        //std::shared_ptr<projection_kernel> qk2(new projection_kernel_q_k_zs_distribution(class_obj.get(), &nofz_s2, zs2));

        std::vector<double> P_ss_lower_limit = { zs_lower };
        std::vector<double> P_ss_upper_limit = { zs_upper };

        std::vector<double> B2D_sss_lower_limit = { zs_lower };
        std::vector<double> B2D_sss_upper_limit = { zs_upper };

        std::vector<double> iB_sss_lower_limits = { zs_lower, 1, 1, 0, 0};
        std::vector<double> iB_sss_upper_limits = { zs_upper, 25000, 25000, 2*M_PI, 2*M_PI};
        //std::vector<double> iB_sss_upper_limits = { zs_upper, 30000, 30000, 2*M_PI, 2*M_PI}; // previous settings

        std::vector<double> iB_sss_angle_averaged_lower_limits = { 0, zs_lower, 1, 1, 0, 0};
        std::vector<double> iB_sss_angle_averaged_upper_limits = { 2*M_PI, zs_upper, 25000, 25000, 2*M_PI, 2*M_PI};

        std::vector<double> iB_sss_4_dim_lower_limits = { 1, 1, 0, 0};
        std::vector<double> iB_sss_4_dim_upper_limits = { 25000, 25000, 2*M_PI, 2*M_PI};

        std::vector<std::shared_ptr<projection_kernel>> qs_kernels;

        for (size_t i = 0; i < zs_bins.size(); i++)
            qs_kernels.emplace_back(new projection_kernel_q_k_zs_fixed(class_obj.get(), zs_bins.at(i)));

        // --------------------------------------------------------
        // Lens (halo/galaxy) correlations settings
        // --------------------------------------------------------

        // Takahashi sample 2
//        double zl_lower = 0.16;
//        double zl_upper = 0.36;

        // redmagicLens1
        double zl_lower = 0.15;
        double zl_upper = 0.35;

        double zl1 = (zl_upper - zl_lower) / 2.0;
        double M_halos_lower = pow(10,13.5)/h; // [M_sun]
        double M_halos_upper = pow(10,17.0)/h; // [M_sun]

        // redshift bins --> make sure that this array is in ascending order for multiple bins
        std::vector<double> zl_bins{zl1}; // lens redshift bin median

        size_t num_2pt_ll_correlations = zl_bins.size();
        size_t num_i3pt_lll_correlations =  zl_bins.size();

        std::vector<std::vector<double>> nofz_l1_table = read_2_column_table("../data/nofz/DESY3_nofz/nofz_DESY3_redmagic_lens_BIN1.tab");
        assert(!nofz_l1_table.empty());
        normalise_nofz(nofz_l1_table);

        //for (size_t idx = 0; idx < nofz_l1_table.at(0).size(); idx++)
        //     std::cout << nofz_l1_table.at(0).at(idx) << "     " << nofz_l1_table.at(1).at(idx) << std::endl;

        Linear_interp_1D nofz_l1(nofz_l1_table.at(0), nofz_l1_table.at(1));
        std::shared_ptr<projection_kernel> qg1(new projection_kernel_q_m(class_obj.get(), &nofz_l1));

        std::shared_ptr<projection_kernel> qh1_b1(new projection_kernel_q_h_b1(class_obj.get(), zl_lower, zl_upper, M_halos_lower, M_halos_upper)); // sample 2 in Takahashi et al.
        std::shared_ptr<projection_kernel> qh1_b2(new projection_kernel_q_h_b2(class_obj.get(), zl_lower, zl_upper, M_halos_lower, M_halos_upper));
        std::shared_ptr<projection_kernel> qh1_bs2(new projection_kernel_q_h_bs2(class_obj.get(), zl_lower, zl_upper, M_halos_lower, M_halos_upper));
        std::shared_ptr<projection_kernel> qh1(new projection_kernel_q_h(class_obj.get(), zl_lower, zl_upper, M_halos_lower, M_halos_upper));

        std::vector<double> P_ll_lower_limit = { zl_lower };
        std::vector<double> P_ll_upper_limit = { zl_upper };

        std::vector<double> iB_lll_lower_limits = { zl_lower, 1, 1, 0, 0};
        std::vector<double> iB_lll_upper_limits = { zl_upper, 25000, 25000, 2*M_PI, 2*M_PI};
        //std::vector<double> iB_lll_upper_limits = { zl_upper, 30000, 30000, 2*M_PI, 2*M_PI}; // previous settings

        std::vector<std::shared_ptr<projection_kernel>> ql_b1_kernels;
        std::vector<std::shared_ptr<projection_kernel>> ql_b2_kernels;
        std::vector<std::shared_ptr<projection_kernel>> ql_bs2_kernels;
        std::vector<std::shared_ptr<projection_kernel>> ql_kernels;

        for (size_t i = 0; i < zl_bins.size(); i++)
        {
            ql_b1_kernels.emplace_back(new projection_kernel_q_h_b1(class_obj.get(), zl_lower, zl_upper, M_halos_lower, M_halos_upper));
            ql_b2_kernels.emplace_back(new projection_kernel_q_h_b2(class_obj.get(), zl_lower, zl_upper, M_halos_lower, M_halos_upper));
            ql_bs2_kernels.emplace_back(new projection_kernel_q_h_bs2(class_obj.get(), zl_lower, zl_upper, M_halos_lower, M_halos_upper));
            ql_kernels.emplace_back(new projection_kernel_q_h(class_obj.get(), zl_lower, zl_upper, M_halos_lower, M_halos_upper));
        }
        double N = std::dynamic_pointer_cast<projection_kernel_q_h_b1>(ql_b1_kernels.at(0))->get_N_h();
        double n_bar = std::dynamic_pointer_cast<projection_kernel_q_h_b1>(ql_b1_kernels.at(0))->get_n_h(); // [Mpc^-3]

        // ----------------------------------------------------------------------------
        // Lens (halo/galaxy) x Source (shear/convergence) cross-correlations settings
        // ----------------------------------------------------------------------------

        std::vector<double> iB_lss_lower_limits = { zl_lower, 1, 1, 0, 0};
        std::vector<double> iB_lss_upper_limits = { zs_upper, 25000, 25000, 2*M_PI, 2*M_PI};

        // -------------------------------------------------------------------------------------

        bool use_pk_nl = true; // default settings
        //bool use_pk_nl = false; // for tree level computations (for halos and galaxies) 

        // -------------------------

        std::string filename_P;

        filename_P = "P_kk.dat"; // P_kk and P_ss (kappa/shear) power spectra are the same but not their 2PCF!
        //filename_P = "P_ss.dat"; // P_kk and P_ss (kappa/shear) power spectra are the same but not their 2PCF!
        //filename_P = "P_hh.dat"; // with bias

        std::string P_integration_algorithm = "qag"; // -- current settings for papers
        //std::string P_integration_algorithm = "mc";
        //std::string P_integration_algorithm = "hcubature"; // -- previous settings

        // -------------------------

        std::string filename_B2D;
        filename_B2D = "B2D_kkk.dat"; // B2D_kkk kappa equilateral bispectrum

        //std::string B2D_integration_algorithm = "qag";
        //std::string B2D_integration_algorithm = "mc";
        std::string B2D_integration_algorithm = "hcubature";

        // -------------------------

        std::string i3pt_area_pre_factors_integration_algorithm = "qag";
        //std::string i3pt_area_pre_factors_integration_algorithm = "mc";
        //std::string i3pt_area_pre_factors_integration_algorithm = "hcubature";

        // -------------------------

        std::string filename_iB;
        filename_iB = "iB_kkk.dat"; // kappa integrated bispectra (k stands for tophat kappa mass)
        //filename_iB = "iB_Mkk.dat"; // kappa integrated bispectra (M stands for aperture mass)
        //filename_iB = "iB_Mss.dat"; // shear integrated bispectra (M stands for aperture mass)
        //filename_iB = "iB_Mss_angle_averaged.dat"; // shear integrated bispectra (M stands for aperture mass) angle averaged
        //filename_iB = "iB_hkk.dat"; // halo integrated bispectra with bias
        //filename_iB = "iB_hhh.dat"; // halo integrated bispectra with bias

        // OLD STUFF --> TO BE DELETED
        // for cross correlations with halos
        //filename_iB = "iBhhh_mc.dat"; // with bias
        //filename_iB = "iBhhh_hcubature.dat"; // with bias

        // for cross correlations with galaxy - TODO: STILL NEEDS to be verified
        //filename_iB = "iBgkk_hcubature.dat"; // without bias
        //filename_iB = "iBggk_hcubature.dat"; // without bias
        //filename_iB = "iBggg_hcubature.dat"; // without bias

        //filename_iB = "iBkkk_hcubature.dat"; // for convergence integrated bispectra
        //filename_iB = "iBkxi_hcubature.dat"; // for shear integrated bispectra
        //filename_iB = "iBkxi_hcubature_v.dat"; // for shear integrated bispectra - vectorized integration TODO: still needs to be verified

        //filename_iB = "iBkxi_hcubature_angle_averaged.dat"; // for shear integrated bispectra - angle averaged TODO: still needs to be verified
        //filename_iB = "iBkxi_mc.dat"; // for shear integrated bispectra - using MC integration

        std::string iB_integration_algorithm = "mc";
        //std::string iB_integration_algorithm = "mc_cigar";
        //std::string iB_integration_algorithm = "hcubature";

        size_t calls_iB_initial;

        if (iB_integration_algorithm == "mc")
            //calls_iB_initial = 2*calls_1e6; // -- current settings for papers
            calls_iB_initial = 1*calls_1e6; // maybe useful when doing for halos or for faster checks
        else if (iB_integration_algorithm == "hcubature")
            calls_iB_initial = calls_1e7;

        // ######################################################################################
        // ######################################################################################
        // ######################################################################################

        // initial checks

        if (compute_initial_checks)
        {
            std::cout << "\nCosmological parameters that have been set in CLASS:" << std::endl;
            std::cout << "H0 [in units of c/Mpc]= " << class_obj->get_H0() << std::endl;
            std::cout << "H0 [in units of c/Mpc]= " << class_obj->get_H_z(0) << std::endl;
            std::cout << "H0 [in units of km/s/Mpc]= " << class_obj->get_H_z(0)*_c_/1.e3 << std::endl;
            std::cout << "h = " << class_obj->get_h() << std::endl;
            std::cout << "A_s = " << class_obj->get_A_s() << std::endl;
            std::cout << "n_s = " << class_obj->get_n_s() << std::endl;
            std::cout << "w0_fld = " << class_obj->get_w0_fld() << std::endl;
            std::cout << "wa_fld = " << class_obj->get_wa_fld() << std::endl;
            std::cout << "Omega0_b = " << class_obj->get_Omega0_b() << std::endl;
            std::cout << "Omega0_cdm = " << class_obj->get_Omega0_cdm() << std::endl;
            std::cout << "Omega0_m != Omega0_b + Omega0_cdm = " << class_obj->get_Omega0_m() << std::endl;
            std::cout << "Omega0_m != Omega0_b + Omega0_cdm = " << class_obj->get_Omega_m_z(0) << std::endl;
            std::cout << "Omega0_g = " << class_obj->get_Omega0_g() << std::endl;
            std::cout << "Omega0_ur = " << class_obj->get_Omega0_ur() << std::endl;
            std::cout << "Omega0_r != Omega0_g + Omega0_ur = " << class_obj->get_Omega0_r() << std::endl;
            std::cout << "Omega0_k = " << class_obj->get_Omega0_k() << std::endl;
            std::cout << "Omega0_Lambda = " << class_obj->get_Omega0_Lambda() << std::endl;
            std::cout << "Omega0_de = " << class_obj->get_Omega0_de() << std::endl;
            std::cout << "Omega0_fld = " << class_obj->get_Omega0_fld() << std::endl;

            std::cout << "\nk_max [1/Mpc]  " << class_obj->get_k_max_pk() << std::endl;

            std::cout << "\nchi [Mpc] = chi(0.4) = " << class_obj->get_chi_z(0.4) << std::endl;
            std::cout << "k_perp [1/Mpc] = l/chi = 100/chi(0.5) = " << 100/class_obj->get_chi_z(0.5) << std::endl;
            std::cout << "P_k_z_shell [Mpc^3]  " << Pk_shell_correction_qag(100/class_obj->get_chi_z(0.5), 0.5, 126./class_obj->get_h(), class_obj.get()) << std::endl;

            std::cout << "chi(z=0.001) [Mpc] = " << class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "chi(z=0.005) [Mpc] = " << class_obj->get_chi_z(0.005) << std::endl;
            std::cout << "chi(z=0.01) [Mpc] = " << class_obj->get_chi_z(0.01) << std::endl;
            std::cout << "chi(z=0.05) [Mpc] = " << class_obj->get_chi_z(0.05) << std::endl;

            std::cout << "k [1/Mpc] = l/chi = 20000/chi(0.001) = " << 20000/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 20000/chi(0.005) = " << 20000/class_obj->get_chi_z(0.005) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 20000/chi(0.01) = " << 20000/class_obj->get_chi_z(0.01) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 20000/chi(0.05) = " << 20000/class_obj->get_chi_z(0.05) << std::endl;

            std::cout << "k [1/Mpc] = l/chi = 0/chi(0.001) = " << 0/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 10000/chi(0.001) = " << 10000/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 20000/chi(0.001) = " << 20000/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 25000/chi(0.001) = " << 25000/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 26000/chi(0.001) = " << 26000/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 30000/chi(0.001) = " << 30000/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 33000/chi(0.001) = " << 33000/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 34000/chi(0.001) = " << 34000/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 35000/chi(0.001) = " << 35000/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 40000/chi(0.001) = " << 40000/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 41000/chi(0.001) = " << 41000/class_obj->get_chi_z(0.001) << std::endl;

            std::cout << "\nk [1/Mpc] = l/chi = 42000/chi(0.001) = " << 42000/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 42000/chi(0.005) = " << 42000/class_obj->get_chi_z(0.005) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 42000/chi(0.01) = " << 42000/class_obj->get_chi_z(0.01) << std::endl;

            std::cout << "\nk [1/Mpc] = l/chi = 50000/chi(0.001) = " << 50000/class_obj->get_chi_z(0.001) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 50000/chi(0.002) = " << 50000/class_obj->get_chi_z(0.002) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 50000/chi(0.003) = " << 50000/class_obj->get_chi_z(0.003) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 50000/chi(0.004) = " << 50000/class_obj->get_chi_z(0.004) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 50000/chi(0.005) = " << 50000/class_obj->get_chi_z(0.005) << std::endl;
            std::cout << "k [1/Mpc] = l/chi = 50000/chi(0.01) = " << 50000/class_obj->get_chi_z(0.01) << std::endl;

            std::cout << "\nqk1(z=0.001) = " << qk1->evaluate(0.001) << std::endl;
            std::cout << "qk1(z=0.005) = " << qk1->evaluate(0.005) << std::endl;
            std::cout << "qk1(z=0.01) = " << qk1->evaluate(0.01) << std::endl;
            std::cout << "qk1(z=0.05) = " << qk1->evaluate(0.05) << std::endl;
            std::cout << "qk1(z=0.23) = " << qk1->evaluate(0.23) << std::endl;
            std::cout << "qk1(z=0.57) = " << qk1->evaluate(0.57) << std::endl;

            std::cout << "\nqk2(z=0.001) = " << qk2->evaluate(0.001) << std::endl;
            std::cout << "qk2(z=0.005) = " << qk2->evaluate(0.005) << std::endl;
            std::cout << "qk2(z=0.01) = " << qk2->evaluate(0.01) << std::endl;
            std::cout << "qk2(z=0.05) = " << qk2->evaluate(0.05) << std::endl;
            std::cout << "qk2(z=0.51) = " << qk2->evaluate(0.51) << std::endl;
            std::cout << "qk2(z=1.02) = " << qk2->evaluate(1.02) << std::endl;

            std::cout << "\nqh1(z=0.01) = " << qh1_b1->evaluate(0.01) << std::endl;
            std::cout << "qh1(z=0.05) = " << qh1_b1->evaluate(0.05) << std::endl;
            std::cout << "qh1(z=0.16) = " << qh1_b1->evaluate(0.16) << std::endl;
            std::cout << "qh1(z=0.23) = " << qh1_b1->evaluate(0.23) << std::endl;
            std::cout << "qh1(z=0.36) = " << qh1_b1->evaluate(0.36) << std::endl;
            std::cout << "qh1(z=0.57) = " << qh1_b1->evaluate(0.57) << std::endl;

            std::cout << "\nN_h(z_min=0.16,z_max=0.36,M_min=1e13.5/h,M_max=1e17/h) = " << N << std::endl;
            std::cout << "\nn_h(z_min=0.16,z_max=0.36,M_min=1e13.5/h,M_max=1e17/h) [(h/Mpc)^3]= " << n_bar / pow(h,3) << std::endl;

            std::cout << "\ntest1" << std::endl;
            std::cout << "F4_sym(test) = " << F4_sym(1., 1., 1., 1., 1./sqrt(2.), 0., 1./sqrt(3.), 0, 2./sqrt(6.), 1./sqrt(3.)) << std::endl;
            std::cout << "\ntest2" << std::endl;
            std::cout << "F4_sym(test) = " << F4_sym(1., 1., 1., 1., 1./sqrt(2.), 0., 1./sqrt(3.), 0, 2./sqrt(6.), -1./sqrt(3.)) << std::endl;

            std::cout << "\nF3 and F4 mode coupling kernel tests" << std::endl;
            std::cout << "F3_sym(0.01, 1.6133604056037736, 1.6133604056037736, 0.02858865975884861, -0.02858865975884861, -0.999999) = " << F3_sym(0.01, 1.6133604056037736, 1.6133604056037736, 0.02858865975884861, -0.02858865975884861, -0.999999) << std::endl;
            std::cout << "F4_sym(0.33863256860613544, 0.33863256860613544, 0.021544346900318832, 0.021544346900318832, -0.999999, 0.04141981391515914, 0.8372570214291091, -0.04141981391515914, -0.8372570214291091, -0.5) = " << F4_sym(0.33863256860613544, 0.33863256860613544, 0.021544346900318832, 0.021544346900318832, -0.999999, 0.04141981391515914, 0.8372570214291091, -0.04141981391515914, -0.8372570214291091, -0.5) << std::endl;

            std::cout << "Comoving lagrangian radius [Mpc] for a halo of mass 1e12 Msun/h = " << R_L(1e12/h, class_obj.get()) << std::endl;

            //Compare with cosmology calculator (https://cosmocalc.icrar.org)

            std::cout << "\nCompare with cosmology calculator (https://cosmocalc.icrar.org) \n" << std::endl;

            double z = 0.5;

            std::cout << "The redshift z is = " << z << std::endl;
            std::cout << "The scale (expansion) factor a is = " << class_obj->get_a_z(z) << std::endl;

            std::cout << "The radial comoving radial distance along the LoS to z is [cMpc] = " << class_obj->get_chi_z(z) << std::endl;
            std::cout << "The luminosity distance to z is [Mpc] = " << class_obj->get_D_lum_z(z) << std::endl;
            std::cout << "The angular diameter distance to z is [pMpc] = " << class_obj->get_D_ang_z(z) << std::endl;
            std::cout << "The comoving volume to z is [cGpc^3] = " << 4.*M_PI/3.0*pow(class_obj->get_chi_z(z),3) / 1e9 << std::endl;

            std::cout << "Hubble's constant H at z is [km/s/pMpc] = " << class_obj->get_H_z(z)*_c_/1.e3 << std::endl;

            std::cout << "The Universe age at z is [Gyr] = " << class_obj->get_proper_time_z(z) / _Gyr_over_Mpc_ << std::endl;
            std::cout << "The look-back time at z is [Gyr] = " << (class_obj->get_proper_time_z(0) - class_obj->get_proper_time_z(z)) / _Gyr_over_Mpc_ << std::endl;

            std::cout << "The Universe age now is [Gyr] = " << class_obj->get_proper_time_z(0) / _Gyr_over_Mpc_ << std::endl;

            std::cout << "The comoving horizon (conformal time/age) of the Universe today is [Mpc] = " << class_obj->get_tau_z(0) << std::endl;
            std::cout << "The comoving horizon (conformal time/age) of the Universe today is [Gyr] = " << class_obj->get_tau_z(0) / _Gyr_over_Mpc_  << std::endl;
            std::cout << "The comoving horizon (conformal time/age) of the Universe at z is [Mpc] = " << class_obj->get_tau_z(z) << std::endl;
            std::cout << "The comoving horizon (conformal time/age) of the Universe at z is [Gyr] = " << class_obj->get_tau_z(z) / _Gyr_over_Mpc_  << std::endl;

            std::cout << "OmegaM at z is = " << class_obj->get_Omega_m_z(z) << std::endl;
            std::cout << "OmegaL at z is = " << class_obj->get_Omega_fld_z(z) << std::endl;
            std::cout << "OmegaR at z is = " << class_obj->get_Omega_r_z(z) << std::endl;
            std::cout << "OmegaK at z is = " << class_obj->get_Omega_k_z(z) << std::endl;

            std::cout << "The growth factor to z is = " << class_obj->get_D_plus_z(z) << std::endl;
            std::cout << "The growth rate at z is = " << class_obj->get_f_z(z) << std::endl;
            std::cout << "The sigma8 at z is = " << class_obj->get_sigma8_z(z) << std::endl;

            std::cout << "The critical energy density [10^10 M_sun/cMpc^3] at z is = " << class_obj->get_rho_crit_z(z) * _kg_m3_to_M_sun_Mpc3_units_ / pow(1.+z,3) /1e10 << std::endl;
            std::cout << "The mean energy (mass) density [10^10 M_sun/cMpc^3] at z is = " << class_obj->get_rho_m_z(z) * _kg_m3_to_M_sun_Mpc3_units_ / pow(1.+z,3) /1e10 << std::endl;
            std::cout << "The critical energy density (SI) at z is = " << class_obj->get_rho_crit_z(z) << std::endl;
            std::cout << "The mean energy (mass) density (SI) at z is = " << class_obj->get_rho_m_z(z) << std::endl;

            std::cout << "The critical energy density [10^10 M_sun/cMpc^3] at z is = " << class_obj->get_rho_crit_z(0) * _kg_m3_to_M_sun_Mpc3_units_ / pow(1.+z,3) /1e10 << std::endl;
            std::cout << "The mean energy (mass) density [10^10 M_sun/cMpc^3] at z is = " << class_obj->get_rho_m_z(0) * _kg_m3_to_M_sun_Mpc3_units_ / pow(1.+z,3) /1e10 << std::endl;
            std::cout << "The critical energy density (SI) today is = " << class_obj->get_rho_crit_z(0) << std::endl;
            std::cout << "The mean energy (mass) density (SI) today is = " << class_obj->get_rho_m_z(0) << std::endl;

            std::cout << "\nThe sigma8 at z=0.00 is = " << class_obj->get_sigma_R_z(8/h, 0.00) << std::endl;
            std::cout << "The sigma8 at z=0.25 is = " << class_obj->get_sigma_R_z(8/h, 0.25) << std::endl;
            std::cout << "The sigma8 at z=0.50 is = " << class_obj->get_sigma_R_z(8/h, 0.50) << std::endl;
            std::cout << "The sigma8 at z=0.75 is = " << class_obj->get_sigma_R_z(8/h, 0.75) << std::endl;
            std::cout << "The sigma8 at z=1.00 is = " << class_obj->get_sigma_R_z(8/h, 1.00) << std::endl;

            std::cout << "The sigma8 (manual integration) at z=0.00 is = " << sqrt(sigma_squared_R_z(8/h, 0.00, class_obj.get())) << std::endl;
            std::cout << "The sigma8 (manual integration) at z=0.25 is = " << sqrt(sigma_squared_R_z(8/h, 0.25, class_obj.get())) << std::endl;
            std::cout << "The sigma8 (manual integration) at z=0.50 is = " << sqrt(sigma_squared_R_z(8/h, 0.50, class_obj.get())) << std::endl;
            std::cout << "The sigma8 (manual integration) at z=0.75 is = " << sqrt(sigma_squared_R_z(8/h, 0.75, class_obj.get())) << std::endl;
            std::cout << "The sigma8 (manual integration) at z=1.00 is = " << sqrt(sigma_squared_R_z(8/h, 1.00, class_obj.get())) << std::endl;

            std::cout << "\nThe sigma squared prime at z=0.00 is = " << class_obj->get_sigma_squared_prime_R_z(8/h, 0.00) << std::endl;
            std::cout << "The sigma squared prime at z=0.25 is = " << class_obj->get_sigma_squared_prime_R_z(8/h, 0.25) << std::endl;
            std::cout << "The sigma squared prime at z=0.50 is = " << class_obj->get_sigma_squared_prime_R_z(8/h, 0.50) << std::endl;
            std::cout << "The sigma squared prime at z=0.75 is = " << class_obj->get_sigma_squared_prime_R_z(8/h, 0.75) << std::endl;
            std::cout << "The sigma squared prime at z=1.00 is = " << class_obj->get_sigma_squared_prime_R_z(8/h, 1.00) << std::endl;

            std::cout << "The sigma squared prime (manual integration) at z=0.00 is = " << sigma_squared_prime_R_z(8/h, 0.00, class_obj.get()) << std::endl;
            std::cout << "The sigma squared prime (manual integration) at z=0.25 is = " << sigma_squared_prime_R_z(8/h, 0.25, class_obj.get()) << std::endl;
            std::cout << "The sigma squared prime (manual integration) at z=0.50 is = " << sigma_squared_prime_R_z(8/h, 0.50, class_obj.get()) << std::endl;
            std::cout << "The sigma squared prime (manual integration) at z=0.75 is = " << sigma_squared_prime_R_z(8/h, 0.75, class_obj.get()) << std::endl;
            std::cout << "The sigma squared prime (manual integration) at z=1.00 is = " << sigma_squared_prime_R_z(8/h, 1.00, class_obj.get()) << std::endl;

        }

        // ######################################################################################

        // sigma(R,z) and sigma'(R,z) table in R=8 Mpc/h

        if (compute_sigma_quantities_tables)
        {
            std::ofstream sigma8_quantities_z_table;

            sigma8_quantities_z_table.open ("./output/sigma8_quantities_z_table.tab", std::ofstream::trunc);

            double a=0,b=1.5;
            int num_pts = 100;
            double z=0.0;
            double R = 8/h; // in [Mpc]
            for(int i=0; i<num_pts; i++)
            {
                z = i*(b-a)/num_pts;
                double sigma_R_z = class_obj->get_sigma_R_z(R,z); // dimensionless
                double dln_sigma8_dln_R_z = class_obj->get_sigma_prime_R_z(R,z)*R/sigma_R_z; // dimensionless
                sigma8_quantities_z_table << z << " " <<  sigma_R_z <<  " " <<  dln_sigma8_dln_R_z << "\n";
            }

            sigma8_quantities_z_table.close();
        }

        // ######################################################################################

        // 3D matter power spectrum and 3D equilateral bispectrum table

        if (compute_Pk_shell_integration_table)
        {
            gettimeofday(&start, nullptr);

            std::ofstream P_shell_integration_k_h_table;

            P_shell_integration_k_h_table.precision(20);

            P_shell_integration_k_h_table.open ("./output/P_shell_integration_k_h_table.tab", std::ofstream::trunc);

            double a=-5,b=2;
            int num_pts = 30;
            double z=0.4;
            double delta_r=delta_r_T17 / h; // in units of Mpc

            std::vector<std::vector<double>> P3D_k(num_pts, std::vector<double>(4, 0));

            #pragma omp parallel for num_threads(thread_count)
            for(int i=0; i<num_pts; i++)
            {
                double k_h = pow(10, a + i * (b - a) / (num_pts - 1)); // in units of h/Mpc

                double k = k_h*h;

                std::cout << k << std::endl;

                P3D_k[i][0] = k_h;
                P3D_k[i][1] = P_nl(k,z,class_obj.get())*pow(h,3);
                P3D_k[i][2] = Pk_shell_correction_qag(k, z, delta_r, class_obj.get())*pow(h,3);
                P3D_k[i][3] = P3D_k[i][2]/P3D_k[i][1];
            }

            for(int idx = 0; idx<num_pts; idx++)
            {
                P_shell_integration_k_h_table << P3D_k.at(idx).at(0) << " " <<  P3D_k.at(idx).at(1) <<  " " <<  P3D_k.at(idx).at(2) <<  " " <<  P3D_k.at(idx).at(3) << "\n";
            }

            P_shell_integration_k_h_table.close();

            gettimeofday(&end, nullptr);

            time_taken = (end.tv_sec - start.tv_sec) * 1e6;
            time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;

            std::cout << "Time taken : " << time_taken << " sec" << std::endl;
        }

        // ######################################################################################

        // 3D matter power spectrum and 3D equilateral bispectrum table

        if (compute_Pk_Bkkk_tables)
        {
            gettimeofday(&start, nullptr);

            std::ofstream P3D_k_h_table;
            std::ofstream B3D_kkk_h_table;

            P3D_k_h_table.precision(20);
            B3D_kkk_h_table.precision(20);

            P3D_k_h_table.open ("./output/P3D_k_h_table.tab", std::ofstream::trunc);
            B3D_kkk_h_table.open ("./output/B3D_kkk_h_table.tab", std::ofstream::trunc);
            //B3D_kkk_h_table.open ("./output/B3D_kh_kh_ks_h_table.tab", std::ofstream::trunc); // for hard and soft k - response function tests

            double a=-3,b=2;
            //double a=-3,b=4.20411998266;
            //double a=-3,b=3.65321251378;
            int num_pts = 100;
            double z=0.0;

            std::vector<std::vector<double>> P3D_k(num_pts, std::vector<double>(3, 0));
            std::vector<std::vector<double>> B3D_kkk(num_pts, std::vector<double>(3, 0));
            //std::vector<std::vector<double>> B3D_kkk(num_pts, std::vector<double>(5, 0));

            //#pragma omp parallel for num_threads(thread_count)
            for(int i=0; i<num_pts; i++)
            {
                double k_h = pow(10, a + i * (b - a) / (num_pts - 1)); // in units of h/Mpc

                double k = k_h*h;

                //std::cout << k << std::endl;

                P3D_k[i][0] = k_h;
                P3D_k[i][1] = P_tree(k,z,class_obj.get())*pow(h,3);
                P3D_k[i][2] = P_nl(k,z,class_obj.get())*pow(h,3);

                B3D_kkk[i][0] = k_h;
                B3D_kkk[i][1] = B_tree(k,k,k,z,class_obj.get(),false)*pow(h,6);
                B3D_kkk[i][2] = B_SC(k,k,k,z,class_obj.get(),true)*pow(h,6);
                //B3D_kkk[i][3] = B_bihalofit(k, k, k, z,class_obj.get(),true)*pow(h,6);
                //B3D_kkk[i][4] = B_1_loop_hcubature(k,k,k,z,class_obj.get(),false)*pow(h,6);

                 // for hard and soft k - response function tests
    //            double k_hard = k; // in units of 1/Mpc
    //            double k_soft = 0.1*h; // in units of 1/Mpc
    //            B3D_kkk[i][1] = B_tree(k_hard, k_hard, k_soft,z,class_obj.get(),false)*pow(h,6);
    //            B3D_kkk[i][2] = B_GM(k_hard, k_hard, k_soft,z,class_obj.get(),true)*pow(h,6);
    //            B3D_kkk[i][3] = B_squeezed_RF(k_hard, k_hard, k_soft, z,class_obj.get(),true)*pow(h,6);
            }

            for(int idx = 0; idx<num_pts; idx++)
            {
                P3D_k_h_table << P3D_k.at(idx).at(0) << " " <<  P3D_k.at(idx).at(1) <<  " " <<  P3D_k.at(idx).at(2) << "\n";
    //            B3D_kkk_h_table << B3D_kkk.at(idx).at(0) << " " <<  B3D_kkk.at(idx).at(1) <<  " " <<  B3D_kkk.at(idx).at(2) << " " <<  B3D_kkk.at(idx).at(3) <<
    //                             " " << B3D_kkk.at(idx).at(4) << "\n";
                B3D_kkk_h_table << B3D_kkk.at(idx).at(0) << " " <<  B3D_kkk.at(idx).at(1) <<  " " <<  B3D_kkk.at(idx).at(2) << "\n";
            }

            P3D_k_h_table.close();
            B3D_kkk_h_table.close();

            gettimeofday(&end, nullptr);

            time_taken = (end.tv_sec - start.tv_sec) * 1e6;
            time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;

            std::cout << "Time taken : " << time_taken << " sec" << std::endl;
        }

        // ######################################################################################

        // n_eff table

        if (compute_n_eff_table)
        {
            std::ofstream n_eff_k_h_table;

            n_eff_k_h_table.open ("./output/n_eff_k_h_table.tab", std::ofstream::trunc);

            double a=-3,b=2;
            int num_pts = 1000;
            double z=1.0;
            for(int i=0; i<=num_pts; i++)
            {
                double k_h = pow(10, a + i * (b - a) / (num_pts - 1)); // in units of h/Mpc
                double n_eff_lin_k_h = class_obj->get_n_eff_from_lin_Pk(k_h*h,z); // ideally, same as class_obj->get_n_eff_from_lin_Pk(k_h*h,0)
                double n_eff_nl_k_h = class_obj->get_n_eff_from_nl_Pk(k_h*h,z);
                n_eff_k_h_table << k_h << " " <<  n_eff_lin_k_h <<  " " <<  n_eff_nl_k_h << "\n";
            }

            n_eff_k_h_table.close();
        }

        // ######################################################################################

        // k_NL table

        if (compute_k_NL_table)
        {
            std::ofstream k_NL_h_z_table;

            k_NL_h_z_table.open ("./output/k_NL_h_z_table.tab", std::ofstream::trunc);

            double a=0,b=1.5;
            int num_pts = 100;
            double z=0.0;
            for(int i=0; i<num_pts; i++)
            {
                z = i*(b-a)/num_pts;
                double k_NL_h = class_obj->get_k_NL_from_lin_Pk_interp(z) / h; // [h/Mpc]
                k_NL_h_z_table << z << " " <<  k_NL_h << "\n";
            }

            k_NL_h_z_table.close();
        }

        // ######################################################################################

        // Response function tables

        if (compute_response_function_tables)
        {
            std::ofstream R_G_1_k_h_table;
            std::ofstream R_G_K_k_h_table;

            R_G_1_k_h_table.open ("./output/R_G_1_k_h_table.tab", std::ofstream::trunc);
            R_G_K_k_h_table.open ("./output/R_G_K_k_h_table.tab", std::ofstream::trunc);

            double a=-2,b=0;
            int num_pts = 1000;
            double z=0.0;
            for(int i=0; i<=num_pts; i++)
            {
                double k_h = pow(10, a + i * (b - a) / (num_pts - 1)); // in units of h/Mpc
                double n = class_obj->get_n_eff_from_nl_Pk(k_h*h,z);
                double R_1 = 1 - n/3.0 + class_obj->get_G_1_k_z_interp(k_h*h,z);
                double R_K = class_obj->get_G_K_k_z_interp(k_h*h,z) - n;

                double G_1 = class_obj->get_G_1_k_z_interp(k_h*h,z);
                double G_K = class_obj->get_G_K_k_z_interp(k_h*h,z);

                R_G_1_k_h_table << k_h << " " <<  R_1 <<  " " <<  G_1 <<"\n";
                R_G_K_k_h_table << k_h << " " <<  R_K <<  " " <<  G_K <<"\n";
            }

            R_G_1_k_h_table.close();
            R_G_K_k_h_table.close();
        }

        // ######################################################################################

        // Transfer functions tables from class

        if (compute_transfer_function_tables)
        {
            std::vector<double> k_Mpc_table; // output k table obtained internally from class in units of [1/Mpc]
            std::vector<double> d_cdm, d_b, d_ncdm, d_tot;
            std::vector<double> t_cdm, t_b, t_ncdm, t_tot;

            // get class transfer functions
            double z=0.0;
            class_obj->getTk(z, k_Mpc_table, d_cdm, d_b, d_ncdm, d_tot, t_cdm, t_b, t_ncdm, t_tot);

            std::ofstream dTk, vTk;

            dTk.open ("./output/dTk.tab", std::ofstream::trunc); // density transfer function
            dTk.precision(40);

            vTk.open ("./output/vTk.tab", std::ofstream::trunc); // velocity transfer function
            vTk.precision(40);

            for(size_t i=0; i<k_Mpc_table.size(); i++)
            {
                // Output table with k in units [h/Mpc]
                dTk << k_Mpc_table.at(i)/class_obj->get_h() << " " << d_cdm.at(i) << " " << d_b.at(i) << " " << d_ncdm.at(i) << " " << d_tot.at(i) << "\n";
                vTk << k_Mpc_table.at(i)/class_obj->get_h() << " " << t_cdm.at(i) << " " << t_b.at(i) << " " << t_ncdm.at(i) << " " << t_tot.at(i) << "\n";
            }

            dTk.close();
            vTk.close();

            // Normalise the transfer function d_tot by -1/k^2 with k in [Mpc] (camb format) and set it to 1 for largest scale (smallest k)
            std::vector<double> d_tot_normalised(k_Mpc_table.size(),0);

            for(size_t i=0; i<k_Mpc_table.size(); i++)
            {
                d_tot_normalised[i] = -d_tot[i]/pow(k_Mpc_table[i],2);
                d_tot_normalised[i] = d_tot_normalised[i] / (-d_tot[0]/pow(k_Mpc_table[0],2));

                std::cout << k_Mpc_table[i]/class_obj->get_h() << " " << d_tot_normalised[i] << std::endl;
            }

            std::ofstream dTk_tot_camb;

            dTk_tot_camb.open ("./output/dTk_tot_camb.tab", std::ofstream::trunc); // total matter density transfer function
            dTk_tot_camb.precision(40);

            for(size_t i=0; i<k_Mpc_table.size(); i++)
            {
                // Output table with k in units [h/Mpc]
                dTk_tot_camb << k_Mpc_table.at(i)/class_obj->get_h() << " " << d_tot_normalised.at(i) << "\n";
            }

            dTk_tot_camb.close();

            //Linear_interp_1D Tk_z_d_tot(k_Mpc_table, d_tot_normalised);
        }

        // ######################################################################################

        // Halo quantities (halo mass function, bias) tables

        if (compute_halo_quantities_tables)
        {
            std::ofstream dn_dlnM_Tinker_M_h_table;
            std::ofstream dn_dlnM_PS_M_h_table;

            std::ofstream halo_b1_Tinker_M_h_table;
            std::ofstream halo_b1_PBS_M_h_table;

            dn_dlnM_Tinker_M_h_table.open ("./output/dn_dlnM_Tinker_M_h_table.tab", std::ofstream::trunc);
            dn_dlnM_PS_M_h_table.open ("./output/dn_dlnM_PS_M_h_table.tab", std::ofstream::trunc);

            halo_b1_Tinker_M_h_table.open ("./output/halo_b1_Tinker_M_h_table.tab", std::ofstream::trunc);
            halo_b1_PBS_M_h_table.open ("./output/halo_b1_PBS_M_h_table.tab", std::ofstream::trunc);

            double a=11,b=15.5;
            int num_pts = 100;

            std::vector<std::vector<double>> dn_dlnM_Tinker_array(num_pts, std::vector<double>(5, 0));
            std::vector<std::vector<double>> dn_dlnM_PS_array(num_pts, std::vector<double>(5, 0));

            std::vector<std::vector<double>> halo_b1_Tinker_array(num_pts, std::vector<double>(5, 0));
            std::vector<std::vector<double>> halo_b1_PBS_array(num_pts, std::vector<double>(5, 0));

            for(int i=0; i<num_pts; i++)
            {
                double M_h = pow(10, a + i * (b - a) / (num_pts - 1)); // in units [M_sun/h]
                double M = M_h / h;

                dn_dlnM_Tinker_array[i][0] = M_h;
                dn_dlnM_Tinker_array[i][1] = M*dn_dM_Tinker2008(M,0.0,class_obj.get())/pow(h,3); // in units [h^3/Mpc^3]
                dn_dlnM_Tinker_array[i][2] = M*dn_dM_Tinker2008(M,1.0,class_obj.get())/pow(h,3); // in units [h^3/Mpc^3]
                dn_dlnM_Tinker_array[i][3] = M*dn_dM_Tinker2008(M,2.0,class_obj.get())/pow(h,3); // in units [h^3/Mpc^3]
                dn_dlnM_Tinker_array[i][4] = M*dn_dM_Tinker2008(M,4.0,class_obj.get())/pow(h,3); // in units [h^3/Mpc^3]

                dn_dlnM_PS_array[i][0] = M_h;
                dn_dlnM_PS_array[i][1] = M*dn_dM_PS(M,0.0,class_obj.get())/pow(h,3); // in units [h^3/Mpc^3]
                dn_dlnM_PS_array[i][2] = M*dn_dM_PS(M,1.0,class_obj.get())/pow(h,3); // in units [h^3/Mpc^3]
                dn_dlnM_PS_array[i][3] = M*dn_dM_PS(M,2.0,class_obj.get())/pow(h,3); // in units [h^3/Mpc^3]
                dn_dlnM_PS_array[i][4] = M*dn_dM_PS(M,4.0,class_obj.get())/pow(h,3); // in units [h^3/Mpc^3]


                halo_b1_Tinker_array[i][0] = M_h;
                halo_b1_Tinker_array[i][1] = halo_b1_Tinker2010(M,0.0,class_obj.get());
                halo_b1_Tinker_array[i][2] = halo_b1_Tinker2010(M,1.0,class_obj.get());
                halo_b1_Tinker_array[i][3] = halo_b1_Tinker2010(M,2.0,class_obj.get());
                halo_b1_Tinker_array[i][4] = halo_b1_Tinker2010(M,4.0,class_obj.get());

                halo_b1_PBS_array[i][0] = M_h;
                halo_b1_PBS_array[i][1] = halo_b1_PBS(M,0.0,class_obj.get());
                halo_b1_PBS_array[i][2] = halo_b1_PBS(M,1.0,class_obj.get());
                halo_b1_PBS_array[i][3] = halo_b1_PBS(M,2.0,class_obj.get());
                halo_b1_PBS_array[i][4] = halo_b1_PBS(M,4.0,class_obj.get());
            }

            for(int idx=0; idx<num_pts; idx++)
            {
                dn_dlnM_Tinker_M_h_table << dn_dlnM_Tinker_array.at(idx).at(0) << " " <<  dn_dlnM_Tinker_array.at(idx).at(1) << " " <<  dn_dlnM_Tinker_array.at(idx).at(2) <<  " " <<  dn_dlnM_Tinker_array.at(idx).at(3) <<  " " <<  dn_dlnM_Tinker_array.at(idx).at(4) << "\n";
                dn_dlnM_PS_M_h_table << dn_dlnM_PS_array.at(idx).at(0) << " " <<  dn_dlnM_PS_array.at(idx).at(1) << " " <<  dn_dlnM_PS_array.at(idx).at(2) <<  " " <<  dn_dlnM_PS_array.at(idx).at(3) <<  " " <<  dn_dlnM_PS_array.at(idx).at(4) << "\n";

                halo_b1_Tinker_M_h_table << halo_b1_Tinker_array.at(idx).at(0) << " " <<  halo_b1_Tinker_array.at(idx).at(1) << " " <<  halo_b1_Tinker_array.at(idx).at(2) <<  " " <<  halo_b1_Tinker_array.at(idx).at(3) <<  " " <<  halo_b1_Tinker_array.at(idx).at(4) << "\n";
                halo_b1_PBS_M_h_table << halo_b1_PBS_array.at(idx).at(0) << " " <<  halo_b1_PBS_array.at(idx).at(1) << " " <<  halo_b1_PBS_array.at(idx).at(2) <<  " " <<  halo_b1_PBS_array.at(idx).at(3) <<  " " <<  halo_b1_PBS_array.at(idx).at(4) << "\n";
            }

            dn_dlnM_Tinker_M_h_table.close();
            dn_dlnM_PS_M_h_table.close();

            halo_b1_Tinker_M_h_table.close();
            halo_b1_PBS_M_h_table.close();
        }

        // ######################################################################################

        /*

        // Test 3D interpolation

        l_array = read_1_column_table("../data/ell_array_test.tab");

        std::vector<double> iB_phi_1_phi_2_lF2_array = read_1_column_table("iB_phi_1_phi_2_lF2.dat");

        double exact_value, error;

        iB_phi_1_phi_2_hcubature("lF2", 344, theta_T, 344, 344, exact_value, error);

        Linear_interp_3D test(l_array, l_array, l_array, iB_phi_1_phi_2_lF2_array);

        std::cout << exact_value << " ; " << test.interp(344, 344, 344) << std::endl;

        iB_phi_1_phi_2_hcubature("X_W", 342, theta_T, 2, 502, exact_value, error);

        std::cout << exact_value << std::endl;

        iB_phi_1_phi_2_hcubature("X_W", 342, theta_T, 502, 2, exact_value, error);

        std::cout << exact_value << std::endl;

        iB_phi_1_phi_2_hcubature("A", 342, theta_T, 2, 502, exact_value, error);

        std::cout << exact_value << std::endl;

        iB_phi_1_phi_2_hcubature("A", 342, theta_T, 502, 2, exact_value, error);

        std::cout << exact_value << std::endl;

        */

        // ######################################################################################

        /*

        // Make interpolation table for iB_phi_phi_2 integrand

        gettimeofday(&start, nullptr);

        l_array = read_1_column_table("../data/ell_array.tab");

        size_t l_nd = l_array.size();
        size_t l_1_nd = l_array.size();
        size_t l_2_nd = l_array.size();
        size_t nd = l_nd*l_1_nd*l_2_nd;

        std::cout << "Size of l array: " << l_nd << std::endl;

        std::cout << "Total number of evaluations to perform: " << nd << std::endl;

        std::vector<double> iB_phi_1_phi_2_lF2_array(nd,0);

        filename = "iB_phi_1_phi_2_lF2.dat";

        counter = 0;

        #pragma omp parallel for num_threads(thread_count) shared(iB_phi_1_phi_2_lF2_array)
        for (size_t i = 0; i < l_nd; i++)
        {
            for (size_t j = 0; j < l_1_nd; j++)
            {
                for (size_t k = 0; k < l_2_nd; k++)
                {
                    size_t idx = k + j*l_2_nd + i*l_2_nd*l_1_nd;

                    double result = 0;
                    double error = 0;

    //                iB_phi_1_phi_2_mc("X_lF2", l, theta_T, l_1, l_2, T, "vegas", phi_1_phi_2_integrand_val, error, calls_1e5);
                    iB_phi_1_phi_2_hcubature("lF2", l_array.at(i), theta_T, l_array.at(j), l_array.at(k), result, error, calls_1e5);

                    iB_phi_1_phi_2_lF2_array[idx] = result;

    //                std::cout << l_array.at(i) << " " << l_array.at(j) << " " << l_array.at(k) << " " << iB_phi_1_phi_2_lF2_array[idx] << std::endl;

                    counter++;

                    if(counter % calls_1e6 == 0)
                        std::cout << counter << std::endl;
                }
            }
        }

        gettimeofday(&end, nullptr);

        time_taken = (end.tv_sec - start.tv_sec) * 1e6;
        time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;

        std::cout << "Time taken : " << time_taken << " sec" << std::endl;

        std::ofstream iB_phi_1_phi_2_lF2;

        iB_phi_1_phi_2_lF2.open (filename, std::ofstream::trunc);
        iB_phi_1_phi_2_lF2.precision(16);

        for (size_t idx = 0; idx < nd; idx++)
        {
            iB_phi_1_phi_2_lF2 << iB_phi_1_phi_2_lF2_array.at(idx) << "\n";
        }

        std::cout<<"\niB_phi_1_phi_2_lF2 output file created\n";

        iB_phi_1_phi_2_lF2.close();

        */

        // ######################################################################################

        // 2D power spectra

        if (compute_2D_power_spectra)
        {
            gettimeofday(&start, nullptr);

            if (verbose_print_outs)
                std::cout << "\nClock started (2D power spectra calculations started)" << std::endl;

            std::vector<std::vector<std::vector<double>>> P_ss_array;
            std::vector<std::vector<std::vector<double>>> P_ll_array;

            P_ss_array = std::vector<std::vector<std::vector<double>>>(1, std::vector<std::vector<double>>(num_2pt_ss_correlations, std::vector<double>(l_array.size(), 0)));
            P_ll_array = std::vector<std::vector<std::vector<double>>>(2, std::vector<std::vector<double>>(num_2pt_ll_correlations, std::vector<double>(l_array.size(), 0)));

            counter = 0;

            #pragma omp parallel for num_threads(thread_count) shared(l_array, class_obj, qs_kernels, ql_b1_kernels)
            for (size_t idx = 0; idx < l_array.size(); idx++)
            {
                if (filename_P == "P_kk.dat" || filename_P == "P_ss.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zs_bins.size() ; a++)
                    {
                        for (size_t b = a; b < zs_bins.size() ; b++)
                        {
                            assert(corr_idx != num_2pt_ss_correlations);

                            if (P_integration_algorithm == "qag")
                                P_ss_array[0][corr_idx][idx] = P2D_z_qag("P", l_array.at(idx), class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), P_ss_lower_limit.at(0), P_ss_upper_limit.at(0));

                            else if (P_integration_algorithm == "mc")
                                P_ss_array[0][corr_idx][idx] = P2D_z_mc("P", l_array.at(idx), class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), P_ss_lower_limit, P_ss_upper_limit, T, "vegas");

                            else if (P_integration_algorithm == "hcubature")
                                P_ss_array[0][corr_idx][idx] = P2D_z_hcubature("P", l_array.at(idx), class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), P_ss_lower_limit, P_ss_upper_limit);

                            corr_idx++;
                        }
                    }
                }

                else if (filename_P == "P_hh.dat")
                {
                    size_t corr_idx = 0;

                    for (size_t a = 0; a < zl_bins.size() ; a++) // only auto-P for lens bins
                    {
                        assert(corr_idx != num_2pt_ll_correlations);

                        if (P_integration_algorithm == "qag")
                        {
                            // b1 term
                            P_ll_array[0][corr_idx][idx] = P2D_z_qag("P", l_array.at(idx), class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), P_ll_lower_limit.at(0), P_ll_upper_limit.at(0));

                            // Shot noise term
                            P_ll_array[1][corr_idx][idx] = P2D_z_qag("P_hh_eps_eps", l_array.at(idx), class_obj.get(), use_pk_nl, ql_kernels.at(a).get(), ql_kernels.at(a).get(), P_ll_lower_limit.at(0), P_ll_upper_limit.at(0));
                        }

                        else if (P_integration_algorithm == "mc")
                        {
                            // b1 term
                            P_ll_array[0][corr_idx][idx] = P2D_z_mc("P", l_array.at(idx), class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), P_ll_lower_limit, P_ll_upper_limit, T, "vegas");

                            // Shot noise term
                            P_ll_array[1][corr_idx][idx] = P2D_z_mc("P_hh_eps_eps", l_array.at(idx), class_obj.get(), use_pk_nl, ql_kernels.at(a).get(), ql_kernels.at(a).get(), P_ll_lower_limit, P_ll_upper_limit, T, "vegas");
                        }

                        else if (P_integration_algorithm == "hcubature")
                        {
                            // b1 term
                            P_ll_array[0][corr_idx][idx] = P2D_z_hcubature("P", l_array.at(idx), class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), P_ll_lower_limit, P_ll_upper_limit);

                            // Shot noise term
                            P_ll_array[1][corr_idx][idx] = P2D_z_hcubature("P_hh_eps_eps", l_array.at(idx), class_obj.get(), use_pk_nl, ql_kernels.at(a).get(), ql_kernels.at(a).get(), P_ll_lower_limit, P_ll_upper_limit);
                        }

                        corr_idx++;
                    }
                }

                if (verbose_print_outs)
                    std::cout << ++counter << " " << l_array.at(idx) << std::endl;
            }

            gettimeofday(&end, nullptr);
            time_taken = (end.tv_sec - start.tv_sec) * 1e6;
            time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
            std::cout << "Time taken for 2D power spectra calculations: " << time_taken << " sec" << std::endl;

            if (filename_P == "P_kk.dat" || filename_P == "P_ss.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zs_bins.size() ; a++)
                {
                    for (size_t b = a; b < zs_bins.size() ; b++)
                    {
                        assert(corr_idx != num_2pt_ss_correlations);

                        std::stringstream s_P;
                        s_P << spectra_folder << "P_" << std::to_string(a+1) << std::to_string(b+1) << filename_extension;

                        std::ofstream file_P;
                        file_P.open (s_P.str(), std::ofstream::trunc);
                        file_P.precision(40);

                        assert(corr_idx != num_2pt_ss_correlations);

                        for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
                            file_P << l_array.at(l_idx) << " " << P_ss_array.at(0).at(corr_idx).at(l_idx) << "\n";

                        file_P.close();

                        corr_idx++;
                    }
                }
            }

           else if (filename_P == "P_hh.dat") // only auto-P for lens bins
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zl_bins.size() ; a++)
                {
                    assert(corr_idx != num_2pt_ll_correlations);

                    std::stringstream s_P;
                    s_P << spectra_folder << "P_" << std::to_string(a+1) << std::to_string(a+1) << filename_extension;

                    std::ofstream file_P;
                    file_P.open (s_P.str(), std::ofstream::trunc);
                    file_P.precision(40);

                    assert(corr_idx != num_2pt_ll_correlations);

                    for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
                        file_P << l_array.at(l_idx) << " " << P_ll_array.at(0).at(corr_idx).at(l_idx) << " " << P_ll_array.at(1).at(corr_idx).at(l_idx) << "\n";

                    file_P.close();

                    corr_idx++;
                }
            }

            if (verbose_print_outs)
                std::cout<<"2D power spectra output files created\n";

            // Spherical sky 2D power spectra (e.g. for FLASK config)

            if (compute_2D_power_spectra_spherical_sky)
            {
                if (filename_P == "P_kk.dat" || filename_P == "P_ss.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zs_bins.size() ; a++)
                    {
                        for (size_t b = a; b < zs_bins.size() ; b++)
                        {
                            assert(corr_idx != num_2pt_ss_correlations);

                            std::vector<int> ells;
                            std::vector<double> C_ells;
                            C_ells_spherical_sky(l_array, P_ss_array.at(0).at(corr_idx), ells, C_ells);

                            std::stringstream s_C_ells;
                            s_C_ells << spectra_folder << "Cl_" << std::to_string(a+1) << std::to_string(b+1) << filename_extension;

                            std::ofstream file_C_ells;
                            file_C_ells.open (s_C_ells.str(), std::ofstream::trunc);
                            file_C_ells.precision(40);

                            for (int l = 0; l <= ells.back(); l++)
                                file_C_ells << l << " " << C_ells.at(l) << "\n";

                            file_C_ells.close();

                            corr_idx++;
                        }
                    }
                }

                else if (filename_P == "P_hh.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zl_bins.size() ; a++)
                    {
                        assert(corr_idx != num_2pt_ll_correlations);

                        std::vector<int> ells;
                        std::vector<double> C_ells;
                        C_ells_spherical_sky(l_array, P_ll_array.at(0).at(corr_idx), ells, C_ells);

                        std::vector<int> ells_n_bar;
                        std::vector<double> C_ells_n_bar;
                        C_ells_spherical_sky(l_array, P_ll_array.at(1).at(corr_idx), ells_n_bar, C_ells_n_bar);

                        std::stringstream s_C_ells;
                        s_C_ells << spectra_folder << "Cl_" << std::to_string(a+1) << std::to_string(a+1) << filename_extension;

                        std::ofstream file_C_ells;
                        file_C_ells.open (s_C_ells.str(), std::ofstream::trunc);
                        file_C_ells.precision(40);

                        for (int l = 0; l <= ells.back(); l++)
                            file_C_ells << l << " " << C_ells.at(l) << " " << C_ells_n_bar.at(l) << "\n";

                        file_C_ells.close();

                        corr_idx++;
                    }
                }

                if (verbose_print_outs)
                    std::cout<<"2D power spectra in spherical sky output files created\n";
            }
        }

        // ######################################################################################

        // 2D 2PCF

        //for (size_t i=0; i < filename_extension_array.size(); i++)
        //{
        //filename_extension = filename_extension_array.at(i);

        if (compute_2D_2PCF)
        {
            gettimeofday(&start, nullptr);

            if (filename_P == "P_kk.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zs_bins.size() ; a++)
                {
                    for (size_t b = a; b < zs_bins.size() ; b++)
                    {
                        assert(corr_idx != num_2pt_ss_correlations);

                        std::vector<double> xi;

                        std::stringstream s_P;
                        s_P << spectra_folder << "P_" << std::to_string(a+1) << std::to_string(b+1) << filename_extension;
                        std::string P_flat_2D_spectrum =  s_P.str();

                        std::vector<std::vector<double>> P_matrix = read_2_column_table(P_flat_2D_spectrum);

                        assert (P_matrix.at(0).size() == l_array.size());

                        //xi = xi_theta_array(alpha_table, P_matrix.at(0), P_matrix.at(1), "xi", thread_count);

                        xi = xi_theta_array_bin_averaged(alpha_min_table, alpha_max_table, P_matrix.at(0), P_matrix.at(1), "xi", thread_count);

                        std::stringstream s_xi2D;
                        s_xi2D << correlations_folder << "xi_" << std::to_string(a+1) << std::to_string(b+1) << filename_extension;
                        std::ofstream file_xi( s_xi2D.str(), std::ofstream::trunc);

                        file_xi.precision(40);

                        for (size_t alpha_idx = 0; alpha_idx < alpha_table.size(); alpha_idx++)
                        {
                            file_xi << alpha_table.at(alpha_idx) << " " << xi.at(alpha_idx) << "\n";
                        }

                        file_xi.close();

                        corr_idx++;
                    }
                }
            }

            if (filename_P == "P_ss.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zs_bins.size() ; a++)
                {
                    for (size_t b = a; b < zs_bins.size() ; b++)
                    {
                        assert(corr_idx != num_2pt_ss_correlations);

                        std::vector<double> xip_Re, xim_Re;

                        std::stringstream s_P;
                        s_P << spectra_folder << "P_" << std::to_string(a+1) << std::to_string(b+1) << filename_extension;
                        std::string P_flat_2D_spectrum =  s_P.str();

                        std::vector<std::vector<double>> P_matrix = read_2_column_table(P_flat_2D_spectrum);

                        assert (P_matrix.at(0).size() == l_array.size());

                        //xip_Re = xi_theta_array(alpha_table, P_matrix.at(0), P_matrix.at(1), "xip", thread_count);
                        //xim_Re = xi_theta_array(alpha_table, P_matrix.at(0), P_matrix.at(1), "xim", thread_count);

                        xip_Re = xi_theta_array_bin_averaged(alpha_min_table, alpha_max_table, P_matrix.at(0), P_matrix.at(1), "xip", thread_count);
                        xim_Re = xi_theta_array_bin_averaged(alpha_min_table, alpha_max_table, P_matrix.at(0), P_matrix.at(1), "xim", thread_count);

                        std::stringstream s_xip_Re;
                        s_xip_Re << correlations_folder << "xip_Re_" << std::to_string(a+1) << std::to_string(b+1) << filename_extension;
                        std::ofstream file_xip_Re( s_xip_Re.str(), std::ofstream::trunc);

                        std::stringstream s_xim_Re;
                        s_xim_Re << correlations_folder << "xim_Re_" << std::to_string(a+1) << std::to_string(b+1) << filename_extension;
                        std::ofstream file_xim_Re( s_xim_Re.str(), std::ofstream::trunc);

                        file_xip_Re.precision(40);
                        file_xim_Re.precision(40);

                        for (size_t alpha_idx = 0; alpha_idx < alpha_table.size(); alpha_idx++)
                        {
                            file_xip_Re << alpha_table.at(alpha_idx) << " " << xip_Re.at(alpha_idx) << "\n";
                            file_xim_Re << alpha_table.at(alpha_idx) << " " << xim_Re.at(alpha_idx) << "\n";
                        }

                        file_xip_Re.close();
                        file_xim_Re.close();

                        corr_idx++;
                    }
                }
            }

            else if (filename_P == "P_hh.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zl_bins.size() ; a++)
                {
                    assert(corr_idx != num_2pt_ll_correlations);

                    std::vector<double> xi, xi_component;

                    std::stringstream s_P;
                    s_P << spectra_folder << "P_" << std::to_string(a+1) << std::to_string(a+1) << filename_extension;
                    std::string P_flat_2D_spectrum =  s_P.str();

                    std::vector<std::vector<double>> P_matrix = read_3_column_table(P_flat_2D_spectrum);

                    assert (P_matrix.at(0).size() == l_array.size());

                    for (size_t b_idx = 1; b_idx < P_matrix.size(); b_idx++)
                    {
                        // Each of the bias modelling components
                        std::stringstream s_xi_component;
                        s_xi_component << correlations_folder << "xi_b" << std::to_string(b_idx) << "_" << std::to_string(a+1) << std::to_string(a+1) << filename_extension;
                        std::ofstream file_xi_component( s_xi_component.str(), std::ofstream::trunc);
                        file_xi_component.precision(40);

                        //xi_component = xi_theta_array(alpha_table, P_matrix.at(0), P_matrix.at(b_idx), "xi", thread_count);

                        xi_component = xi_theta_array_bin_averaged(alpha_min_table, alpha_max_table, P_matrix.at(0), P_matrix.at(b_idx), "xi", thread_count);

                        for (size_t alpha_idx = 0; alpha_idx < alpha_table.size(); alpha_idx++)
                        {
                            file_xi_component << alpha_table.at(alpha_idx) << " " << xi_component.at(alpha_idx) << "\n";
                        }

                        file_xi_component.close();
                    }

                    std::vector<double> P(P_matrix.at(0).size(), 0);

                    // Total sum of the bias modelling components
                    for (size_t l_idx = 0; l_idx < P_matrix.at(0).size(); l_idx++)
                    {
                        for (size_t b_idx = 1; b_idx < P_matrix.size(); b_idx++)
                        {
                            if (b_idx == 2) // don't add Poisson shot noise term as it gives Dirac delta anyway
                                continue;

                            else
                                P[l_idx] += P_matrix[b_idx][l_idx];
                        }
                    }

                    //xi = xi_theta_array(alpha_table, P_matrix.at(0), P, "xi", thread_count);

                    xi = xi_theta_array_bin_averaged(alpha_min_table, alpha_max_table, P_matrix.at(0), P, "xi", thread_count);

                    std::stringstream s_xi2D;
                    s_xi2D << correlations_folder << "xi_" << std::to_string(a+1) << std::to_string(a+1) << filename_extension;
                    std::ofstream file_xi( s_xi2D.str(), std::ofstream::trunc);

                    file_xi.precision(40);

                    for (size_t alpha_idx = 0; alpha_idx < alpha_table.size(); alpha_idx++)
                    {
                        file_xi << alpha_table.at(alpha_idx) << " " << xi.at(alpha_idx) << "\n";
                    }

                    file_xi.close();

                    corr_idx++;
                }
            }

            gettimeofday(&end, nullptr);
            time_taken = (end.tv_sec - start.tv_sec) * 1e6;
            time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
            std::cout << "\nTime taken for 2D 2PCF calculations: " << time_taken << " sec" << std::endl;

            if (verbose_print_outs)
                std::cout<<"2D 2PCF output files created\n\n";

        }

        //}

        // ######################################################################################

        // 2D bispectra equilateral configurations

        if (compute_2D_bispectra_equilateral)
        {
            gettimeofday(&start, nullptr);

            std::cout << "\nClock started (2D bispectra calculations started)" << std::endl;

            std::vector<std::vector<std::vector<double>>> B2D_sss_array;

            B2D_sss_array = std::vector<std::vector<std::vector<double>>>(1, std::vector<std::vector<double>>(num_3pt_sss_correlations, std::vector<double>(l_array.size(), 0)));

            counter = 0;

            #pragma omp parallel for num_threads(thread_count) shared(l_array, class_obj, qs_kernels, ql_b1_kernels)
            for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
            {
                if (filename_B2D == "B2D_kkk.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zs_bins.size() ; a++)
                    {
                        for (size_t b = a; b < zs_bins.size() ; b++)
                        {
                            for (size_t c = b; c < zs_bins.size() ; c++)
                            {
                                assert(corr_idx != num_3pt_sss_correlations);

                                if (B2D_integration_algorithm == "qag")
                                    B2D_sss_array[0][corr_idx][l_idx] = B2D_z_qag("B", l_array.at(l_idx), l_array.at(l_idx), l_array.at(l_idx), class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), B2D_sss_lower_limit.at(0), B2D_sss_upper_limit.at(0));

                                else if (B2D_integration_algorithm == "mc")
                                    B2D_sss_array[0][corr_idx][l_idx] = B2D_z_mc("B", l_array.at(l_idx), l_array.at(l_idx), l_array.at(l_idx), class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), B2D_sss_lower_limit, B2D_sss_upper_limit, T, "vegas");

                                else if (B2D_integration_algorithm == "hcubature")
                                    B2D_sss_array[0][corr_idx][l_idx] = B2D_z_hcubature("B", l_array.at(l_idx), l_array.at(l_idx), l_array.at(l_idx), class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(),  qs_kernels.at(c).get(), B2D_sss_lower_limit, B2D_sss_upper_limit);

                                corr_idx++;
                            }
                        }
                    }
                }

                if (verbose_print_outs)
                    std::cout << ++counter << " " << l_array.at(l_idx) << std::endl;
            }

            gettimeofday(&end, nullptr);

            time_taken = (end.tv_sec - start.tv_sec) * 1e6;
            time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;

            std::cout << "\nTime taken for 2D bispectra calculations: " << time_taken << " sec" << std::endl;

            if (filename_B2D == "B2D_kkk.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zs_bins.size() ; a++)
                {
                    for (size_t b = a; b < zs_bins.size() ; b++)
                    {
                        for (size_t c = b; c < zs_bins.size() ; c++)
                        {
                            assert(corr_idx != num_3pt_sss_correlations);

                            std::stringstream s_B2D;
                            s_B2D << spectra_folder << "B2D_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;

                            std::ofstream file_B2D;
                            file_B2D.open (s_B2D.str(), std::ofstream::trunc);
                            file_B2D.precision(40);

                            for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
                                file_B2D << l_array.at(l_idx) << " " << B2D_sss_array.at(0).at(corr_idx).at(l_idx) << "\n";

                            file_B2D.close();

                            corr_idx++;
                        }
                    }
                }
            }

            if (verbose_print_outs)
                std::cout<<"2D bispectra output files created\n";
        }

        // ######################################################################################

        // 2D integrated 3PCF area pre-factors

        if (compute_2D_integrated_3PCF_area_pre_factors)
        {
            std::vector<double> A2pt(alpha_table.size(),0);

            double Adelta = Adelta_qag(theta_T);

            #pragma omp parallel for num_threads(thread_count)
            for (size_t alpha_idx = 0; alpha_idx < alpha_table.size(); alpha_idx++)
            {
                if (i3pt_area_pre_factors_integration_algorithm == "qag")
                {
                    A2pt[alpha_idx] = A2pt_qag(alpha_table.at(alpha_idx)*M_PI/180.0/60.0, theta_T);
                    //A2pt[alpha_idx] = A2pt_angle_averaged_qag(alpha_table.at(alpha_idx)*M_PI/180.0/60.0, theta_T); // angle-averaged phi_alpha (NOT NEEDED)
                }

                else if (i3pt_area_pre_factors_integration_algorithm == "mc")
                {
                    const gsl_rng_type *T = gsl_rng_default;
                    A2pt[alpha_idx] =  A2pt_mc(alpha_table.at(alpha_idx)*M_PI/180.0/60.0, theta_T, T, "vegas");
                }

                else if (i3pt_area_pre_factors_integration_algorithm == "hcubature")
                {
                    A2pt[alpha_idx] =  A2pt_hcubature(alpha_table.at(alpha_idx)*M_PI/180.0/60.0, theta_T);
                }

                std::cout << alpha_table.at(alpha_idx) << " " << Adelta << " " << A2pt.at(alpha_idx) << std::endl;
            }

            std::ofstream area_pre_factor;
            std::stringstream s_area_pre_factor;
            s_area_pre_factor << correlations_folder << "i3pt_area_pre_factors_" << i3pt_area_pre_factors_integration_algorithm << ".dat";
            area_pre_factor.open (s_area_pre_factor.str(), std::ofstream::trunc);
            area_pre_factor.precision(16);

            for (size_t alpha_idx = 0; alpha_idx < alpha_table.size(); alpha_idx++)
            {
                area_pre_factor << alpha_table.at(alpha_idx) << " " << Adelta << " " << A2pt.at(alpha_idx) << "\n";
            }

            if (verbose_print_outs)
                std::cout<<"2D integrated 3PCF area pre-factors output file created\n";

            area_pre_factor.close();
        }

        // ######################################################################################

        // 2D integrated bispectra

        /*
        if (compute_2D_integrated_bispectra)
        {
            std::vector<std::vector<double>> iB_111(l_array.size(), std::vector<double>(11, 0));
            std::vector<std::vector<double>> iB_112(l_array.size(), std::vector<double>(11, 0));
            std::vector<std::vector<double>> iB_122(l_array.size(), std::vector<double>(11, 0));
            std::vector<std::vector<double>> iB_222(l_array.size(), std::vector<double>(11, 0));

            gettimeofday(&start, nullptr);

            std::cout << "\nClock started (2D integrated bispectra calculations started)" << std::endl;

            double iB_P_1_val = 0, iB_A_val = 0;

            struct_iB2D_W_FS info_iB_W_FS = {&W2D_U_FS, theta_U, &W2D_TH_FS, theta_T};
            struct_iB2D_W_FS info_iB_W_FS_halos = {&W2D_TH_FS, theta_T, &W2D_TH_FS, theta_T};

            if (filename_iB == "iBggg_mc.dat")
            {
                iB_P_1_val = iB2D_mc("P_1", 0, info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", calls_iB_initial);
                std::cout << "\niB_P_1 = " << iB_P_1_val << std::endl;

                iB_A_val = iB2D_mc("A", 0, info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", calls_iB_initial);

                std::cout << "\niB_A = " << iB_A_val << std::endl;
            }

            else if (filename_iB == "iBggg_hcubature.dat")
            {
                iB_P_1_val = iB2D_hcubature("P_1", 0, info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_lll_lower_limits, iB_lll_upper_limits, calls_iB_initial);
                std::cout << "\niB_P_1 = " << iB_P_1_val << std::endl;

                iB_A_val = iB2D_hcubature("A", 0, info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_lll_lower_limits, iB_lll_upper_limits, calls_iB_initial);

                std::cout << "\niB_A = " << iB_A_val << std::endl;
            }

            counter = 0;

            #pragma omp parallel for num_threads(thread_count) shared(l_array, class_obj, qg1, qk1, qk2)
            for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
            {
                size_t calls_iB;

                calls_iB = calls_iB_initial;

                if (l_array.at(l_idx) <= 150)
                    calls_iB = 2*calls_iB_initial;
                else
                    calls_iB = calls_iB_initial;

                if (filename_iB == "iBhhh_mc.dat") // with bias
                {
                    iB_111[l_idx][0] = iB2D_mc("B", l_array.at(l_idx), info_iB_W_FS_halos, class_obj.get(), true, qh1_b1.get(), qh1_b1.get(), qh1_b1.get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", calls_iB);
                    iB_111[l_idx][1] = iB2D_mc("B_PP", l_array.at(l_idx), info_iB_W_FS_halos, class_obj.get(), true, qh1_b1.get(), qh1_b1.get(), qh1_b2.get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", calls_iB);
                }

                if (filename_iB == "iBhhh_hcubature.dat") // with bias
                {
                    iB_111[l_idx][0] = iB2D_hcubature("B", l_array.at(l_idx), info_iB_W_FS_halos, class_obj.get(), true, qh1_b1.get(), qh1_b1.get(), qh1_b1.get(), iB_lll_lower_limits, iB_lll_upper_limits, calls_iB);
                }

                if (filename_iB == "iBgkk_hcubature.dat")
                {
                     // without bias
                    iB_111[l_idx][0] = iB2D_hcubature("B_lF2", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qk1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][3] = iB2D_hcubature("Z_W", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qk1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][6] = iB2D_hcubature("Z_S2", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qk1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                }

                else if (filename_iB == "iBggk_hcubature.dat")
                {
                     // without bias
                    iB_111[l_idx][0] = iB2D_hcubature("B_lF2", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][2] = iB2D_hcubature("Y_W", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][3] = iB2D_hcubature("Z_W", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][5] = iB2D_hcubature("Y_S2", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][6] = iB2D_hcubature("Z_S2", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][9] = iB2D_hcubature("P_3", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                }

                else if (filename_iB == "iBggg_hcubature.dat")
                {
                     // without bias
                    iB_111[l_idx][0] = iB2D_hcubature("B_lF2", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][1] = iB2D_hcubature("X_W", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][2] = iB2D_hcubature("Y_W", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][3] = iB2D_hcubature("Z_W", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][4] = iB2D_hcubature("X_S2", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][5] = iB2D_hcubature("Y_S2", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][6] = iB2D_hcubature("Z_S2", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][7] = iB_P_1_val;
                    iB_111[l_idx][8] = iB2D_hcubature("P_2", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][9] = iB2D_hcubature("P_3", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][10] = iB_A_val;
                }

                if (filename_iB == "iBkkk_hcubature.dat")
                {
                    iB_111[l_idx][0] = iB2D_hcubature("B", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);

                    if (do_cross_fields)
                    {
                        iB_112[l_idx][0] = iB2D_hcubature("B", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                        iB_122[l_idx][0] = iB2D_hcubature("B", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                        iB_222[l_idx][0] = iB2D_hcubature("B", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                   }
                }

                else if (filename_iB == "iBkxi_hcubature.dat")
                {
                    iB_111[l_idx][0] = iB2D_hcubature("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][1] = iB2D_hcubature("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
        //            iB_111[idx][2] = iB2D_hcubature("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_lower_limits, iB_upper_limits, calls_iB);
        //            iB_111[idx][3] = iB2D_hcubature("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_lower_limits, iB_upper_limits, calls_iB);

                    if (do_cross_fields)
                    {
                        iB_112[l_idx][0] = iB2D_hcubature("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                        iB_112[l_idx][1] = iB2D_hcubature("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
        //                iB_112[idx][2] = iB2D_hcubature("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);
        //                iB_112[idx][3] = iB2D_hcubature("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);

                        iB_122[l_idx][0] = iB2D_hcubature("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                        iB_122[l_idx][1] = iB2D_hcubature("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
        //                iB_122[idx][2] = iB2D_hcubature("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);
        //                iB_122[idx][3] = iB2D_hcubature("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);

                        iB_222[l_idx][0] = iB2D_hcubature("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                        iB_222[l_idx][1] = iB2D_hcubature("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
        //                iB_222[idx][2] = iB2D_hcubature("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);
        //                iB_222[idx][3] = iB2D_hcubature("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);
                    }
                }

                else if (filename_iB == "iBkxi_hcubature_v.dat")
                {
                    iB_111[l_idx][0] = iB2D_hcubature_v("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][1] = iB2D_hcubature_v("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
        //            iB_111[idx][2] = iB2D_hcubature_v("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_lower_limits, iB_upper_limits, calls_iB);
        //            iB_111[idx][3] = iB2D_hcubature_v("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_lower_limits, iB_upper_limits, calls_iB);

                    if (do_cross_fields)
                    {
                        iB_112[l_idx][0] = iB2D_hcubature_v("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                        iB_112[l_idx][1] = iB2D_hcubature_v("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
        //                iB_112[idx][2] = iB2D_hcubature_v("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);
        //                iB_112[idx][3] = iB2D_hcubature_v("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);

                        iB_122[l_idx][0] = iB2D_hcubature_v("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                        iB_122[l_idx][1] = iB2D_hcubature_v("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
        //                iB_122[idx][2] = iB2D_hcubature_v("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);
        //                iB_122[idx][3] = iB2D_hcubature_v("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);

                        iB_222[l_idx][0] = iB2D_hcubature_v("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                        iB_222[l_idx][1] = iB2D_hcubature_v("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
        //                iB_222[idx][2] = iB2D_hcubature_v("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);
        //                iB_222[idx][3] = iB2D_hcubature_v("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);
                    }
                }

                else if (filename_iB == "iBkxi_hcubature_angle_averaged.dat")
                {
                    iB_111[l_idx][0] = iB2D_hcubature_angle_averaged("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                    iB_111[l_idx][1] = iB2D_hcubature_angle_averaged("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
        //            iB_111[idx][2] = iB2D_hcubature_angle_averaged("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_lower_limits, iB_upper_limits, calls_iB);
        //            iB_111[idx][3] = iB2D_hcubature_angle_averaged("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_lower_limits, iB_upper_limits, calls_iB);

                    if (do_cross_fields)
                    {
                        iB_112[l_idx][0] = iB2D_hcubature_angle_averaged("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                        iB_112[l_idx][1] = iB2D_hcubature_angle_averaged("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
        //                iB_112[idx][2] = iB2D_hcubature_angle_averaged("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);
        //                iB_112[idx][3] = iB2D_hcubature_angle_averaged("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);

                        iB_122[l_idx][0] = iB2D_hcubature_angle_averaged("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                        iB_122[l_idx][1] = iB2D_hcubature_angle_averaged("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
        //                iB_122[idx][2] = iB2D_hcubature_angle_averaged("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);
        //                iB_122[idx][3] = iB2D_hcubature_angle_averaged("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);

                        iB_222[l_idx][0] = iB2D_hcubature_angle_averaged("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
                        iB_222[l_idx][1] = iB2D_hcubature_angle_averaged("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, calls_iB);
        //                iB_222[idx][2] = iB2D_hcubature_angle_averaged("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);
        //                iB_222[idx][3] = iB2D_hcubature_angle_averaged("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, calls_iB);
                    }
                }

                else if (filename_iB == "iBkxi_mc.dat")
                {
                    iB_111[l_idx][0] = iB2D_mc("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", calls_iB);
                    iB_111[l_idx][1] = iB2D_mc("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", calls_iB);
        //            iB_111[idx][2] = iB2D_mc("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_lower_limits, iB_upper_limits, T, "vegas", calls_iB);
        //            iB_111[idx][3] = iB2D_mc("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk1.get(), iB_lower_limits, iB_upper_limits, T, "vegas", calls_iB);

                    if (do_cross_fields)
                    {
                        iB_112[l_idx][0] = iB2D_mc("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", calls_iB);
                        iB_112[l_idx][1] = iB2D_mc("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", calls_iB);
        //                iB_112[idx][2] = iB2D_mc("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_lower_limits, iB_upper_limits, T, "vegas", calls_iB);
        //                iB_112[idx][3] = iB2D_mc("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk1.get(), qk2.get(), iB_lower_limits, iB_upper_limits, T, "vegas", calls_iB);

                        iB_122[l_idx][0] = iB2D_mc("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", calls_iB);
                        iB_122[l_idx][1] = iB2D_mc("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", calls_iB);
        //                iB_122[idx][2] = iB2D_mc("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, T, "vegas", calls_iB);
        //                iB_122[idx][3] = iB2D_mc("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk1.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, T, "vegas", calls_iB);

                        iB_222[l_idx][0] = iB2D_mc("B_xip_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", calls_iB);
                        iB_222[l_idx][1] = iB2D_mc("B_xim_cos", l_array.at(l_idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", calls_iB);
        //                iB_222[idx][2] = iB2D_mc("B_xip_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, T, "vegas", calls_iB);
        //                iB_222[idx][3] = iB2D_mc("B_xim_sin", l_array.at(idx), info_iB_W_FS, class_obj.get(), true, qk2.get(), qk2.get(), qk2.get(), iB_lower_limits, iB_upper_limits, T, "vegas", calls_iB);

                    }
                }

                std::cout << ++counter << " " << l_array.at(l_idx) << " " << iB_111[l_idx][0] << std::endl;
            }

            gettimeofday(&end, nullptr);

            time_taken = (end.tv_sec - start.tv_sec) * 1e6;
            time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;

            std::cout << "\nTime taken for 2D integrated bispectra calculations: " << time_taken << " sec" << std::endl;

            std::ofstream file_iB_111, file_iB_222, file_iB_122, file_iB_211, file_iB_112, file_iB_212;

            std::stringstream s_iB_111;
            s_iB_111 << spectra_folder << "iB_111" << filename_extension;
            file_iB_111.open (s_iB_111.str(), std::ofstream::trunc);
            file_iB_111.precision(40);

            std::stringstream s_iB_112;
            s_iB_112 << spectra_folder << "iB_112" << filename_extension;
            file_iB_112.open (s_iB_112.str(), std::ofstream::trunc);
            file_iB_112.precision(40);

            std::stringstream s_iB_122;
            s_iB_122 << spectra_folder << "iB_122" << filename_extension;
            file_iB_122.open (s_iB_122.str(), std::ofstream::trunc);
            file_iB_122.precision(40);

            std::stringstream s_iB_222;
            s_iB_222 << spectra_folder << "iB_222" << filename_extension;
            file_iB_222.open (s_iB_222.str(), std::ofstream::trunc);
            file_iB_222.precision(40);

            for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
            {
                file_iB_111 << l_array.at(l_idx) << " " << iB_111.at(l_idx).at(0) << " " << iB_111.at(l_idx).at(1) << " " << iB_111.at(l_idx).at(2) <<
                              " " << iB_111.at(l_idx).at(3) << " " << iB_111.at(l_idx).at(4) << " " << iB_111.at(l_idx).at(5) <<
                              " " << iB_111.at(l_idx).at(6) << " " << iB_111.at(l_idx).at(7) << " " << iB_111.at(l_idx).at(8) <<
                              " " << iB_111.at(l_idx).at(9) << " " << iB_111.at(l_idx).at(10) << "\n";

                if (do_cross_fields)
                {
                    file_iB_112 << l_array.at(l_idx) << " " << iB_112.at(l_idx).at(0) << " " << iB_112.at(l_idx).at(1) << " " << iB_112.at(l_idx).at(2) <<
                                  " " << iB_112.at(l_idx).at(3) << " " << iB_112.at(l_idx).at(4) << " " << iB_112.at(l_idx).at(5) <<
                                  " " << iB_112.at(l_idx).at(6) << " " << iB_112.at(l_idx).at(7) << " " << iB_112.at(l_idx).at(8) <<
                                  " " << iB_112.at(l_idx).at(9) << " " << iB_112.at(l_idx).at(10) << "\n";

                    file_iB_122 << l_array.at(l_idx) << " " << iB_122.at(l_idx).at(0) << " " << iB_122.at(l_idx).at(1) << " " << iB_122.at(l_idx).at(2) <<
                                  " " << iB_122.at(l_idx).at(3) << " " << iB_122.at(l_idx).at(4) << " " << iB_122.at(l_idx).at(5) <<
                                  " " << iB_122.at(l_idx).at(6) << " " << iB_122.at(l_idx).at(7) << " " << iB_122.at(l_idx).at(8) <<
                                  " " << iB_122.at(l_idx).at(9) << " " << iB_122.at(l_idx).at(10) << "\n";

                    file_iB_222 << l_array.at(l_idx) << " " << iB_222.at(l_idx).at(0) << " " << iB_222.at(l_idx).at(1) << " " << iB_222.at(l_idx).at(2) <<
                                  " " << iB_222.at(l_idx).at(3) << " " << iB_222.at(l_idx).at(4) << " " << iB_222.at(l_idx).at(5) <<
                                  " " << iB_222.at(l_idx).at(6) << " " << iB_222.at(l_idx).at(7) << " " << iB_222.at(l_idx).at(8) <<
                                  " " << iB_222.at(l_idx).at(9) << " " << iB_222.at(l_idx).at(10) << "\n";
                }
            }

            file_iB_111.close();
            file_iB_112.close();
            file_iB_122.close();
            file_iB_222.close();

            std::cout<<"2D integrated bispectra output files created\n";
        }
        */

        if (compute_2D_integrated_bispectra_v2)
        {
            gettimeofday(&start, nullptr);

            std::cout << "\nClock started (2D integrated bispectra calculations started)" << std::endl;

    //        double iB_P_1_val = 0, iB_A_val = 0;

            struct_iB2D_W_FS info_iB_UWW_FS = {&W2D_U_FS, theta_U, &W2D_TH_FS, theta_T};
            struct_iB2D_W_FS info_iB_WWW_FS = {&W2D_TH_FS, theta_T, &W2D_TH_FS, theta_T};

    //        if (filename_iB == "iB_ggg_mc.dat")
    //        {
    //            iB_P_1_val = iB2D_mc("P_1", 0, info_iB_WWW_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_source_lower_limits, iB_source_upper_limits, T, "vegas", calls_iB_initial);
    //            std::cout << "\niB_P_1 = " << iB_P_1_val << std::endl;

    //            iB_A_val = iB2D_mc("A", 0, info_iB_WWW_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_source_lower_limits, iB_source_upper_limits, T, "vegas", calls_iB_initial);

    //            std::cout << "\niB_A = " << iB_A_val << std::endl;
    //        }

    //        else if (filename_iB == "iB_ggg_hcubature.dat")
    //        {
    //            iB_P_1_val = iB2D_hcubature("P_1", 0, info_iB_WWW_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_source_lower_limits, iB_source_upper_limits, calls_iB_initial);
    //            std::cout << "\niB_P_1 = " << iB_P_1_val << std::endl;

    //            iB_A_val = iB2D_hcubature("A", 0, info_iB_WWW_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_source_lower_limits, iB_source_upper_limits, calls_iB_initial);

    //            std::cout << "\niB_A = " << iB_A_val << std::endl;
    //        }

            std::vector<std::vector<std::vector<double>>> iB_sss_array(2, std::vector<std::vector<double>>(num_i3pt_sss_correlations, std::vector<double>(l_array.size(), 0))); // 2 columns to accommodate either (iB_Mkk) or (iB_Mxip and iB_Mxim) columns
            std::vector<std::vector<std::vector<double>>> iB_lll_array(5, std::vector<std::vector<double>>(num_i3pt_lll_correlations, std::vector<double>(l_array.size(), 0))); // 5 columns to accommodate (iB_hhh_b1, iB_hhh_b2, iB_hhh_bs2, iB_hhh_sn1, iB_hhh_sn2)
            //std::vector<std::vector<std::vector<double>>> iB_lls_array(5, std::vector<std::vector<double>>(num_i3pt_lls_correlations, std::vector<double>(l_array.size(), 0))); // 5 columns to accommodate (iB_hhh_b1, iB_hhh_b2, iB_hhh_bs2, iB_hhh_sn1, iB_hhh_sn2)
            std::vector<std::vector<std::vector<double>>> iB_lss_array(6, std::vector<std::vector<double>>(zl_bins.size()*num_2pt_ss_correlations, std::vector<double>(l_array.size(), 0))); // 6 columns to accommodate either (iB_hkk_b1, iB_hkk_b2, iB_hkk_bs2) or (iB_hxip_b1, iB_hxip_b2, iB_hxip_bs2, iB_hxim_b1, iB_hxim_b2, iB_hxim_bs2)

            std::vector<std::vector<std::vector<double>>> iB_sss_error_array(2, std::vector<std::vector<double>>(num_i3pt_sss_correlations, std::vector<double>(l_array.size(), 0)));
            std::vector<std::vector<std::vector<double>>> iB_lll_error_array(5, std::vector<std::vector<double>>(num_i3pt_lll_correlations, std::vector<double>(l_array.size(), 0)));
            std::vector<std::vector<std::vector<double>>> iB_lss_error_array(6, std::vector<std::vector<double>>(zl_bins.size()*num_2pt_ss_correlations, std::vector<double>(l_array.size(), 0)));

            counter = 0;

            VEGAS_Integrator cigar;
            cigar.Set_Verbose(NONE);

            #pragma omp parallel for num_threads(thread_count) shared(l_array, class_obj, qs_kernels, ql_b1_kernels, ql_b2_kernels, ql_bs2_kernels)
            for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
            //for (size_t l_idx = l_array.size()-1; l_idx >= 0; l_idx--)
            {
                size_t calls_iB;

                if (l_array.at(l_idx) <= 220)
                    //calls_iB = 2*calls_iB_initial;
                    calls_iB = calls_iB_initial;
                else
                    calls_iB = calls_iB_initial;

                // --------------------------------------------------------

                if (filename_iB == "iB_kkk.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zs_bins.size() ; a++)
                    {
                        for (size_t b = a; b < zs_bins.size() ; b++)
                        {
                            for (size_t c = b; c < zs_bins.size() ; c++)
                            {
                                assert(corr_idx != num_i3pt_sss_correlations);

                                double result = 0.0, error = 0.0;

                                if (iB_integration_algorithm == "mc")
                                {
                                    // iB term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;
                                }

                                if (iB_integration_algorithm == "mc_cigar")
                                {
                                    // iB term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc_cigar("B", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, cigar, result, error, thread_count);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;
                                }

                                else if (iB_integration_algorithm == "hcubature")
                                {
                                    // iB term
                                    result = 0.0, error = 0.0;
                                    iB2D_hcubature("B", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;
                                }

                                corr_idx++;
                            }
                        }
                    }
                }

                else if (filename_iB == "iB_Mkk.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zs_bins.size() ; a++)
                    {
                        for (size_t b = a; b < zs_bins.size() ; b++)
                        {
                            for (size_t c = b; c < zs_bins.size() ; c++)
                            {
                                assert(corr_idx != num_i3pt_sss_correlations);

                                double result = 0.0, error = 0.0;

                                if (iB_integration_algorithm == "mc")
                                {
                                    // iB term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;
                                }

                                if (iB_integration_algorithm == "mc_cigar")
                                {
                                    // iB term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc_cigar("B", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, cigar, result, error, thread_count);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;
                                }

                                else if (iB_integration_algorithm == "hcubature")
                                {
                                    // iB term
                                    result = 0.0, error = 0.0;
                                    iB2D_hcubature("B", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;
                                }

                                corr_idx++;
                            }
                        }
                    }
                }

                else if (filename_iB == "iB_Mss.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zs_bins.size() ; a++)
                    {
                        for (size_t b = a; b < zs_bins.size() ; b++)
                        {
                            for (size_t c = b; c < zs_bins.size() ; c++)
                            {
                                assert(corr_idx != num_i3pt_sss_correlations);

                                double result = 0.0, error = 0.0;

                                if (iB_integration_algorithm == "mc")
                                {
                                    // iBp term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B_xip_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;

                                    // iBm term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B_xim_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_sss_array[1][corr_idx][l_idx] = result;
                                    iB_sss_error_array[1][corr_idx][l_idx] = error;
                                }

                                if (iB_integration_algorithm == "mc_cigar")
                                {
                                    // iBp term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc_cigar("B_xip_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, cigar, result, error, thread_count);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;

                                    // iBm term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc_cigar("B_xim_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, cigar, result, error, thread_count);
                                    iB_sss_array[1][corr_idx][l_idx] = result;
                                    iB_sss_error_array[1][corr_idx][l_idx] = error;
                                }

                                else if (iB_integration_algorithm == "hcubature")
                                {
                                    // iBp term
                                    result = 0.0, error = 0.0;
                                    iB2D_hcubature("B_xip_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;

                                    // iBm term
                                    result = 0.0, error = 0.0;
                                    iB2D_hcubature("B_xim_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, result, error, calls_iB);
                                    iB_sss_array[1][corr_idx][l_idx] = result;
                                    iB_sss_error_array[1][corr_idx][l_idx] = error;
                                }

                                corr_idx++;
                            }
                        }
                    }
                }

                else if (filename_iB == "iB_Mss_angle_averaged.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zs_bins.size() ; a++)
                    {
                        for (size_t b = a; b < zs_bins.size() ; b++)
                        {
                            for (size_t c = b; c < zs_bins.size() ; c++)
                            {
                                assert(corr_idx != num_i3pt_sss_correlations);

                                double result = 0.0, error = 0.0;

                                if (iB_integration_algorithm == "mc")
                                {
                                    // iBp term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc_angle_averaged("B_xip_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_angle_averaged_lower_limits, iB_sss_angle_averaged_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;

                                    // iBm term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc_angle_averaged("B_xim_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_angle_averaged_lower_limits, iB_sss_angle_averaged_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_sss_array[1][corr_idx][l_idx] = result;
                                    iB_sss_error_array[1][corr_idx][l_idx] = error;
                                }

                                else if (iB_integration_algorithm == "hcubature")
                                {
                                    // iBp term
                                    result = 0.0, error = 0.0;
                                    iB2D_hcubature_angle_averaged("B_xip_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_angle_averaged_lower_limits, iB_sss_angle_averaged_upper_limits, result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;

                                    // iBm term
                                    result = 0.0, error = 0.0;
                                    iB2D_hcubature_angle_averaged("B_xim_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_angle_averaged_lower_limits, iB_sss_angle_averaged_upper_limits, result, error, calls_iB);
                                    iB_sss_array[1][corr_idx][l_idx] = result;
                                    iB_sss_error_array[1][corr_idx][l_idx] = error;
                                }

                                corr_idx++;
                            }
                        }
                    }
                }

                else if (filename_iB == "iB_hkk.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zl_bins.size() ; a++)
                    {
                        for (size_t b = 0; b < zs_bins.size() ; b++)
                        {
                            for (size_t c = b; c < zs_bins.size() ; c++)
                            {
                                assert(corr_idx != zl_bins.size()*num_2pt_ss_correlations);

                                double result = 0.0, error = 0.0;

                                if (iB_integration_algorithm == "mc")
                                {
                                    // b1 term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_lss_array[0][corr_idx][l_idx] = result;
                                    iB_lss_error_array[0][corr_idx][l_idx] = error;

                                    // b2 term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B_P2P3", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_b2_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_lss_array[1][corr_idx][l_idx] = result;
                                    iB_lss_error_array[1][corr_idx][l_idx] = error;

                                    // bs2 term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B_S2P2P3", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_bs2_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_lss_array[2][corr_idx][l_idx] = result;
                                    iB_lss_error_array[2][corr_idx][l_idx] = error;
                                }

                                corr_idx++;
                            }
                        }
                    }
                }

                else if (filename_iB == "iB_hhh.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zl_bins.size() ; a++)
                    {
                        assert(corr_idx != num_i3pt_lll_correlations);

                        double result = 0.0, error = 0.0;

                        if (iB_integration_algorithm == "mc")
                        {
                            // b1 term
                            result = 0.0, error = 0.0;
                            iB2D_mc("B", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                            iB_lll_array[0][corr_idx][l_idx] = result;
                            iB_lll_error_array[0][corr_idx][l_idx] = error;

                            // b2 term
                            result = 0.0, error = 0.0;
                            iB2D_mc("B_PP", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), ql_b2_kernels.at(a).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                            iB_lll_array[1][corr_idx][l_idx] = result;
                            iB_lll_error_array[1][corr_idx][l_idx] = error;

                            // bs2 term
                            result = 0.0, error = 0.0;
                            iB2D_mc("B_S2PP", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), ql_bs2_kernels.at(a).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                            iB_lll_array[2][corr_idx][l_idx] = result;
                            iB_lll_error_array[2][corr_idx][l_idx] = error;

                            // Shot noise 1 term
                            result = 0.0, error = 0.0;
                            iB2D_mc("B_hhh_delta_eps_eps", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), ql_kernels.at(a).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                            iB_lll_array[3][corr_idx][l_idx] = result;
                            iB_lll_error_array[3][corr_idx][l_idx] = error;

                            // Shot noise 2 term
                            result = 0.0, error = 0.0;
                            iB2D_mc("B_hhh_eps_eps_eps", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_kernels.at(a).get(), ql_kernels.at(a).get(), ql_kernels.at(a).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                            iB_lll_array[4][corr_idx][l_idx] = result;
                            iB_lll_error_array[4][corr_idx][l_idx] = error;
                        }

                        corr_idx++;
                    }
                }

                if (verbose_print_outs)
                    std::cout << ++counter << " " << l_array.at(l_idx) << std::endl;
            }

            gettimeofday(&end, nullptr);
            time_taken = (end.tv_sec - start.tv_sec) * 1e6;
            time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
            std::cout << "\nTime taken for 2D integrated bispectra calculations: " << time_taken << " sec" << std::endl;

            if (filename_iB == "iB_kkk.dat" || filename_iB == "iB_Mkk.dat" || filename_iB == "iB_Mss.dat" || filename_iB == "iB_Mss_angle_averaged.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zs_bins.size() ; a++)
                {
                    for (size_t b = a; b < zs_bins.size() ; b++)
                    {
                        for (size_t c = b; c < zs_bins.size() ; c++)
                        {
                            std::stringstream s_iB;
                            s_iB << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;

                            std::ofstream file_iB;
                            file_iB.open (s_iB.str(), std::ofstream::trunc);
                            file_iB.precision(40);

                            std::stringstream s_iB_error;
                            s_iB_error << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << "_error" << filename_extension;

                            std::ofstream file_iB_error;
                            file_iB_error.open (s_iB_error.str(), std::ofstream::trunc);
                            file_iB_error.precision(40);

                            assert(corr_idx != num_i3pt_sss_correlations);

                            for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
                            {
                                file_iB << l_array.at(l_idx) << " " << iB_sss_array[0][corr_idx][l_idx] <<  " " << iB_sss_array[1][corr_idx][l_idx] << " " << "\n";
                                file_iB_error << l_array.at(l_idx) << " " << iB_sss_error_array[0][corr_idx][l_idx] <<  " " << iB_sss_error_array[1][corr_idx][l_idx] << " " << "\n";
                            }

                            file_iB.close();
                            file_iB_error.close();

                            corr_idx++;
                        }
                    }
                }
            }

            else if (filename_iB == "iB_hkk.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zl_bins.size() ; a++)
                {
                    for (size_t b = 0; b < zs_bins.size() ; b++)
                    {
                        for (size_t c = b; c < zs_bins.size() ; c++)
                        {
                            std::stringstream s_iB;
                            s_iB << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;

                            std::ofstream file_iB;
                            file_iB.open (s_iB.str(), std::ofstream::trunc);
                            file_iB.precision(40);

                            std::stringstream s_iB_error;
                            s_iB_error << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << "_error" << filename_extension;

                            std::ofstream file_iB_error;
                            file_iB_error.open (s_iB_error.str(), std::ofstream::trunc);
                            file_iB_error.precision(40);

                            assert(corr_idx != num_i3pt_lll_correlations);

                            for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
                            {
                                file_iB << l_array.at(l_idx) << " " << iB_lss_array[0][corr_idx][l_idx] <<  " " << iB_lss_array[1][corr_idx][l_idx] << " " <<
                                           iB_lss_array[2][corr_idx][l_idx] << " " << iB_lss_array[3][corr_idx][l_idx] << " " << iB_lss_array[4][corr_idx][l_idx] << " " <<
                                           iB_lss_array[5][corr_idx][l_idx] << " " << "\n";

                                file_iB_error << l_array.at(l_idx) << " " << iB_lss_error_array[0][corr_idx][l_idx] <<  " " << iB_lss_error_array[1][corr_idx][l_idx] << " " <<
                                                 iB_lss_error_array[2][corr_idx][l_idx] << " " << iB_lss_error_array[3][corr_idx][l_idx] << " " << iB_lss_error_array[4][corr_idx][l_idx] << " " <<
                                                 iB_lss_error_array[5][corr_idx][l_idx] << " " << "\n";

                            }

                            file_iB.close();
                            file_iB_error.close();

                            corr_idx++;
                        }
                    }
                }
            }

            else if (filename_iB == "iB_hhh.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zl_bins.size() ; a++)
                {
                    std::stringstream s_iB;
                    s_iB << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(a+1) << std::to_string(a+1) << filename_extension;

                    std::ofstream file_iB;
                    file_iB.open (s_iB.str(), std::ofstream::trunc);
                    file_iB.precision(40);

                    std::stringstream s_iB_error;
                    s_iB_error << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(a+1) << std::to_string(a+1) << "_error" << filename_extension;

                    std::ofstream file_iB_error;
                    file_iB_error.open (s_iB_error.str(), std::ofstream::trunc);
                    file_iB_error.precision(40);

                    assert(corr_idx != num_i3pt_lll_correlations);

                    for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
                    {
                        file_iB << l_array.at(l_idx) << " " << iB_lll_array[0][corr_idx][l_idx] <<  " " << iB_lll_array[1][corr_idx][l_idx] << " " <<
                                   iB_lll_array[2][corr_idx][l_idx] << " " << iB_lll_array[3][corr_idx][l_idx] << " " << iB_lll_array[4][corr_idx][l_idx] << " " << "\n";

                        file_iB_error << l_array.at(l_idx) << " " << iB_lll_error_array[0][corr_idx][l_idx] <<  " " << iB_lll_error_array[1][corr_idx][l_idx] << " " <<
                                         iB_lll_error_array[2][corr_idx][l_idx] << " " << iB_lll_error_array[3][corr_idx][l_idx] << " " << iB_lll_error_array[4][corr_idx][l_idx] << " " << "\n";

                    }

                    file_iB.close();
                    file_iB_error.close();

                    corr_idx++;
                }
            }

            if (verbose_print_outs)
                std::cout<<"2D integrated bispectra output files created\n";
        }

        // ######################################################################################

        if (compute_2D_integrated_bispectra_v3)
        {
            gettimeofday(&start, nullptr);

            std::cout << "\nClock started (2D integrated bispectra calculations started)" << std::endl;

    //        double iB_P_1_val = 0, iB_A_val = 0;

            struct_iB2D_W_FS info_iB_UWW_FS = {&W2D_U_FS, theta_U, &W2D_TH_FS, theta_T};
            struct_iB2D_W_FS info_iB_WWW_FS = {&W2D_TH_FS, theta_T, &W2D_TH_FS, theta_T};

    //        if (filename_iB == "iB_ggg_mc.dat")
    //        {
    //            iB_P_1_val = iB2D_mc("P_1", 0, info_iB_WWW_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_source_lower_limits, iB_source_upper_limits, T, "vegas", calls_iB_initial);
    //            std::cout << "\niB_P_1 = " << iB_P_1_val << std::endl;

    //            iB_A_val = iB2D_mc("A", 0, info_iB_WWW_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_source_lower_limits, iB_source_upper_limits, T, "vegas", calls_iB_initial);

    //            std::cout << "\niB_A = " << iB_A_val << std::endl;
    //        }

    //        else if (filename_iB == "iB_ggg_hcubature.dat")
    //        {
    //            iB_P_1_val = iB2D_hcubature("P_1", 0, info_iB_WWW_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_source_lower_limits, iB_source_upper_limits, calls_iB_initial);
    //            std::cout << "\niB_P_1 = " << iB_P_1_val << std::endl;

    //            iB_A_val = iB2D_hcubature("A", 0, info_iB_WWW_FS, class_obj.get(), true, qg1.get(), qg1.get(), qg1.get(), iB_source_lower_limits, iB_source_upper_limits, calls_iB_initial);

    //            std::cout << "\niB_A = " << iB_A_val << std::endl;
    //        }

            std::vector<std::vector<std::vector<double>>> iB_sss_array(2, std::vector<std::vector<double>>(num_i3pt_sss_correlations, std::vector<double>(l_array.size(), 0))); // 2 columns to accommodate either (iB_Mkk) or (iB_Mxip and iB_Mxim) columns
            std::vector<std::vector<std::vector<double>>> iB_lll_array(5, std::vector<std::vector<double>>(num_i3pt_lll_correlations, std::vector<double>(l_array.size(), 0))); // 5 columns to accommodate (iB_hhh_b1, iB_hhh_b2, iB_hhh_bs2, iB_hhh_sn1, iB_hhh_sn2)
            //std::vector<std::vector<std::vector<double>>> iB_lls_array(5, std::vector<std::vector<double>>(num_i3pt_lls_correlations, std::vector<double>(l_array.size(), 0))); // 5 columns to accommodate (iB_hhh_b1, iB_hhh_b2, iB_hhh_bs2, iB_hhh_sn1, iB_hhh_sn2)
            std::vector<std::vector<std::vector<double>>> iB_lss_array(6, std::vector<std::vector<double>>(zl_bins.size()*num_2pt_ss_correlations, std::vector<double>(l_array.size(), 0))); // 6 columns to accommodate either (iB_hkk_b1, iB_hkk_b2, iB_hkk_bs2) or (iB_hxip_b1, iB_hxip_b2, iB_hxip_bs2, iB_hxim_b1, iB_hxim_b2, iB_hxim_bs2)

            std::vector<std::vector<std::vector<double>>> iB_sss_error_array(2, std::vector<std::vector<double>>(num_i3pt_sss_correlations, std::vector<double>(l_array.size(), 0)));
            std::vector<std::vector<std::vector<double>>> iB_lll_error_array(5, std::vector<std::vector<double>>(num_i3pt_lll_correlations, std::vector<double>(l_array.size(), 0)));
            std::vector<std::vector<std::vector<double>>> iB_lss_error_array(6, std::vector<std::vector<double>>(zl_bins.size()*num_2pt_ss_correlations, std::vector<double>(l_array.size(), 0)));

            counter = 0;

            VEGAS_Integrator cigar;
            cigar.Set_Verbose(NONE);

            #pragma omp parallel for num_threads(thread_count) shared(l_array, class_obj, qs_kernels, ql_b1_kernels, ql_b2_kernels, ql_bs2_kernels)
            for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
            //for (size_t l_idx = l_array.size()-1; l_idx >= 0; l_idx--)
            {
                size_t calls_iB;

                if (l_array.at(l_idx) <= 220)
                    //calls_iB = 2*calls_iB_initial;
                    calls_iB = calls_iB_initial;
                else
                    calls_iB = calls_iB_initial;

                // --------------------------------------------------------

                if (filename_iB == "iB_kkk.dat")
                {
                    std::vector<double> z_array;
                    for (double zi = 0.001; zi < zs_upper+0.001; zi+=0.1)
                        z_array.push_back(zi);

                    std::vector<double> iB2D_z_array(z_array.size(), 0);
                    for (size_t z_idx = 0; z_idx < z_array.size(); z_idx++)
                    {
                        double result = 0.0, error = 0.0;
                        iB2D_mc_4_dim("B", l_array.at(l_idx), z_array.at(z_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, iB_sss_4_dim_lower_limits, iB_sss_4_dim_upper_limits, T, "vegas", result, error, 1e5);
                        iB2D_z_array[z_idx] = result;
                    }

                    Linear_interp_1D iB2D_z_interp(z_array, iB2D_z_array);

                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zs_bins.size() ; a++)
                    {
                        for (size_t b = a; b < zs_bins.size() ; b++)
                        {
                            for (size_t c = b; c < zs_bins.size() ; c++)
                            {
                                assert(corr_idx != num_i3pt_sss_correlations);

                                if (iB_integration_algorithm == "mc")
                                {
                                    // iB term
                                    iB_sss_array[0][corr_idx][l_idx] = iB2D_trapz_z(&iB2D_z_interp, zs_lower, zs_upper, class_obj.get(), qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get());
                                }

                                corr_idx++;
                            }
                        }
                    }
                }

                else if (filename_iB == "iB_Mkk.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zs_bins.size() ; a++)
                    {
                        for (size_t b = a; b < zs_bins.size() ; b++)
                        {
                            for (size_t c = b; c < zs_bins.size() ; c++)
                            {
                                assert(corr_idx != num_i3pt_sss_correlations);

                                double result = 0.0, error = 0.0;

                                if (iB_integration_algorithm == "mc")
                                {
                                    // iB term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;
                                }

                                if (iB_integration_algorithm == "mc_cigar")
                                {
                                    // iB term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc_cigar("B", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, cigar, result, error, thread_count);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;
                                }

                                else if (iB_integration_algorithm == "hcubature")
                                {
                                    // iB term
                                    result = 0.0, error = 0.0;
                                    iB2D_hcubature("B", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;
                                }

                                corr_idx++;
                            }
                        }
                    }
                }

                /*

                else if (filename_iB == "iB_Mss.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zs_bins.size() ; a++)
                    {
                        for (size_t b = a; b < zs_bins.size() ; b++)
                        {
                            for (size_t c = b; c < zs_bins.size() ; c++)
                            {
                                assert(corr_idx != num_i3pt_sss_correlations);

                                double result = 0.0, error = 0.0;

                                if (iB_integration_algorithm == "mc")
                                {
                                    // iBp term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B_xip_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;

                                    // iBm term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B_xim_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_sss_array[1][corr_idx][l_idx] = result;
                                    iB_sss_error_array[1][corr_idx][l_idx] = error;
                                }

                                if (iB_integration_algorithm == "mc_cigar")
                                {
                                    // iBp term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc_cigar("B_xip_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, cigar, result, error, thread_count);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;

                                    // iBm term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc_cigar("B_xim_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, cigar, result, error, thread_count);
                                    iB_sss_array[1][corr_idx][l_idx] = result;
                                    iB_sss_error_array[1][corr_idx][l_idx] = error;
                                }

                                else if (iB_integration_algorithm == "hcubature")
                                {
                                    // iBp term
                                    result = 0.0, error = 0.0;
                                    iB2D_hcubature("B_xip_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;

                                    // iBm term
                                    result = 0.0, error = 0.0;
                                    iB2D_hcubature("B_xim_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_lower_limits, iB_sss_upper_limits, result, error, calls_iB);
                                    iB_sss_array[1][corr_idx][l_idx] = result;
                                    iB_sss_error_array[1][corr_idx][l_idx] = error;
                                }

                                corr_idx++;
                            }
                        }
                    }
                }

                else if (filename_iB == "iB_Mss_angle_averaged.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zs_bins.size() ; a++)
                    {
                        for (size_t b = a; b < zs_bins.size() ; b++)
                        {
                            for (size_t c = b; c < zs_bins.size() ; c++)
                            {
                                assert(corr_idx != num_i3pt_sss_correlations);

                                double result = 0.0, error = 0.0;

                                if (iB_integration_algorithm == "mc")
                                {
                                    // iBp term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc_angle_averaged("B_xip_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_angle_averaged_lower_limits, iB_sss_angle_averaged_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;

                                    // iBm term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc_angle_averaged("B_xim_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_angle_averaged_lower_limits, iB_sss_angle_averaged_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_sss_array[1][corr_idx][l_idx] = result;
                                    iB_sss_error_array[1][corr_idx][l_idx] = error;
                                }

                                else if (iB_integration_algorithm == "hcubature")
                                {
                                    // iBp term
                                    result = 0.0, error = 0.0;
                                    iB2D_hcubature_angle_averaged("B_xip_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_angle_averaged_lower_limits, iB_sss_angle_averaged_upper_limits, result, error, calls_iB);
                                    iB_sss_array[0][corr_idx][l_idx] = result;
                                    iB_sss_error_array[0][corr_idx][l_idx] = error;

                                    // iBm term
                                    result = 0.0, error = 0.0;
                                    iB2D_hcubature_angle_averaged("B_xim_cos", l_array.at(l_idx), info_iB_UWW_FS, class_obj.get(), use_pk_nl, qs_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_sss_angle_averaged_lower_limits, iB_sss_angle_averaged_upper_limits, result, error, calls_iB);
                                    iB_sss_array[1][corr_idx][l_idx] = result;
                                    iB_sss_error_array[1][corr_idx][l_idx] = error;
                                }

                                corr_idx++;
                            }
                        }
                    }
                }

                else if (filename_iB == "iB_hkk.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zl_bins.size() ; a++)
                    {
                        for (size_t b = 0; b < zs_bins.size() ; b++)
                        {
                            for (size_t c = b; c < zs_bins.size() ; c++)
                            {
                                assert(corr_idx != zl_bins.size()*num_2pt_ss_correlations);

                                double result = 0.0, error = 0.0;

                                if (iB_integration_algorithm == "mc")
                                {
                                    // b1 term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_lss_array[0][corr_idx][l_idx] = result;
                                    iB_lss_error_array[0][corr_idx][l_idx] = error;

                                    // b2 term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B_P2P3", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_b2_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_lss_array[1][corr_idx][l_idx] = result;
                                    iB_lss_error_array[1][corr_idx][l_idx] = error;

                                    // bs2 term
                                    result = 0.0, error = 0.0;
                                    iB2D_mc("B_S2P2P3", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_bs2_kernels.at(a).get(), qs_kernels.at(b).get(), qs_kernels.at(c).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                                    iB_lss_array[2][corr_idx][l_idx] = result;
                                    iB_lss_error_array[2][corr_idx][l_idx] = error;
                                }

                                corr_idx++;
                            }
                        }
                    }
                }

                else if (filename_iB == "iB_hhh.dat")
                {
                    size_t corr_idx = 0;
                    for (size_t a = 0; a < zl_bins.size() ; a++)
                    {
                        assert(corr_idx != num_i3pt_lll_correlations);

                        double result = 0.0, error = 0.0;

                        if (iB_integration_algorithm == "mc")
                        {
                            // b1 term
                            result = 0.0, error = 0.0;
                            iB2D_mc("B", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                            iB_lll_array[0][corr_idx][l_idx] = result;
                            iB_lll_error_array[0][corr_idx][l_idx] = error;

                            // b2 term
                            result = 0.0, error = 0.0;
                            iB2D_mc("B_PP", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), ql_b2_kernels.at(a).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                            iB_lll_array[1][corr_idx][l_idx] = result;
                            iB_lll_error_array[1][corr_idx][l_idx] = error;

                            // bs2 term
                            result = 0.0, error = 0.0;
                            iB2D_mc("B_S2PP", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), ql_bs2_kernels.at(a).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                            iB_lll_array[2][corr_idx][l_idx] = result;
                            iB_lll_error_array[2][corr_idx][l_idx] = error;

                            // Shot noise 1 term
                            result = 0.0, error = 0.0;
                            iB2D_mc("B_hhh_delta_eps_eps", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_b1_kernels.at(a).get(), ql_b1_kernels.at(a).get(), ql_kernels.at(a).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                            iB_lll_array[3][corr_idx][l_idx] = result;
                            iB_lll_error_array[3][corr_idx][l_idx] = error;

                            // Shot noise 2 term
                            result = 0.0, error = 0.0;
                            iB2D_mc("B_hhh_eps_eps_eps", l_array.at(l_idx), info_iB_WWW_FS, class_obj.get(), use_pk_nl, ql_kernels.at(a).get(), ql_kernels.at(a).get(), ql_kernels.at(a).get(), iB_lll_lower_limits, iB_lll_upper_limits, T, "vegas", result, error, calls_iB);
                            iB_lll_array[4][corr_idx][l_idx] = result;
                            iB_lll_error_array[4][corr_idx][l_idx] = error;
                        }

                        corr_idx++;
                    }
                }
                */

                if (verbose_print_outs)
                    std::cout << ++counter << " " << l_array.at(l_idx) << std::endl;
            }

            gettimeofday(&end, nullptr);
            time_taken = (end.tv_sec - start.tv_sec) * 1e6;
            time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
            std::cout << "\nTime taken for 2D integrated bispectra calculations: " << time_taken << " sec" << std::endl;

            if (filename_iB == "iB_kkk.dat" || filename_iB == "iB_Mkk.dat" || filename_iB == "iB_Mss.dat" || filename_iB == "iB_Mss_angle_averaged.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zs_bins.size() ; a++)
                {
                    for (size_t b = a; b < zs_bins.size() ; b++)
                    {
                        for (size_t c = b; c < zs_bins.size() ; c++)
                        {
                            std::stringstream s_iB;
                            s_iB << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;

                            std::ofstream file_iB;
                            file_iB.open (s_iB.str(), std::ofstream::trunc);
                            file_iB.precision(40);

                            std::stringstream s_iB_error;
                            s_iB_error << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << "_error" << filename_extension;

                            std::ofstream file_iB_error;
                            file_iB_error.open (s_iB_error.str(), std::ofstream::trunc);
                            file_iB_error.precision(40);

                            assert(corr_idx != num_i3pt_sss_correlations);

                            for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
                            {
                                file_iB << l_array.at(l_idx) << " " << iB_sss_array[0][corr_idx][l_idx] <<  " " << iB_sss_array[1][corr_idx][l_idx] << " " << "\n";
                                file_iB_error << l_array.at(l_idx) << " " << iB_sss_error_array[0][corr_idx][l_idx] <<  " " << iB_sss_error_array[1][corr_idx][l_idx] << " " << "\n";
                            }

                            file_iB.close();
                            file_iB_error.close();

                            corr_idx++;
                        }
                    }
                }
            }

            else if (filename_iB == "iB_hkk.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zl_bins.size() ; a++)
                {
                    for (size_t b = 0; b < zs_bins.size() ; b++)
                    {
                        for (size_t c = b; c < zs_bins.size() ; c++)
                        {
                            std::stringstream s_iB;
                            s_iB << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;

                            std::ofstream file_iB;
                            file_iB.open (s_iB.str(), std::ofstream::trunc);
                            file_iB.precision(40);

                            std::stringstream s_iB_error;
                            s_iB_error << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << "_error" << filename_extension;

                            std::ofstream file_iB_error;
                            file_iB_error.open (s_iB_error.str(), std::ofstream::trunc);
                            file_iB_error.precision(40);

                            assert(corr_idx != num_i3pt_lll_correlations);

                            for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
                            {
                                file_iB << l_array.at(l_idx) << " " << iB_lss_array[0][corr_idx][l_idx] <<  " " << iB_lss_array[1][corr_idx][l_idx] << " " <<
                                           iB_lss_array[2][corr_idx][l_idx] << " " << iB_lss_array[3][corr_idx][l_idx] << " " << iB_lss_array[4][corr_idx][l_idx] << " " <<
                                           iB_lss_array[5][corr_idx][l_idx] << " " << "\n";

                                file_iB_error << l_array.at(l_idx) << " " << iB_lss_error_array[0][corr_idx][l_idx] <<  " " << iB_lss_error_array[1][corr_idx][l_idx] << " " <<
                                                 iB_lss_error_array[2][corr_idx][l_idx] << " " << iB_lss_error_array[3][corr_idx][l_idx] << " " << iB_lss_error_array[4][corr_idx][l_idx] << " " <<
                                                 iB_lss_error_array[5][corr_idx][l_idx] << " " << "\n";

                            }

                            file_iB.close();
                            file_iB_error.close();

                            corr_idx++;
                        }
                    }
                }
            }

            else if (filename_iB == "iB_hhh.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zl_bins.size() ; a++)
                {
                    std::stringstream s_iB;
                    s_iB << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(a+1) << std::to_string(a+1) << filename_extension;

                    std::ofstream file_iB;
                    file_iB.open (s_iB.str(), std::ofstream::trunc);
                    file_iB.precision(40);

                    std::stringstream s_iB_error;
                    s_iB_error << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(a+1) << std::to_string(a+1) << "_error" << filename_extension;

                    std::ofstream file_iB_error;
                    file_iB_error.open (s_iB_error.str(), std::ofstream::trunc);
                    file_iB_error.precision(40);

                    assert(corr_idx != num_i3pt_lll_correlations);

                    for (size_t l_idx = 0; l_idx < l_array.size(); l_idx++)
                    {
                        file_iB << l_array.at(l_idx) << " " << iB_lll_array[0][corr_idx][l_idx] <<  " " << iB_lll_array[1][corr_idx][l_idx] << " " <<
                                   iB_lll_array[2][corr_idx][l_idx] << " " << iB_lll_array[3][corr_idx][l_idx] << " " << iB_lll_array[4][corr_idx][l_idx] << " " << "\n";

                        file_iB_error << l_array.at(l_idx) << " " << iB_lll_error_array[0][corr_idx][l_idx] <<  " " << iB_lll_error_array[1][corr_idx][l_idx] << " " <<
                                         iB_lll_error_array[2][corr_idx][l_idx] << " " << iB_lll_error_array[3][corr_idx][l_idx] << " " << iB_lll_error_array[4][corr_idx][l_idx] << " " << "\n";

                    }

                    file_iB.close();
                    file_iB_error.close();

                    corr_idx++;
                }
            }

            if (verbose_print_outs)
                std::cout<<"2D integrated bispectra output files created\n";
        }

        // ######################################################################################

        // 2D integrated 3PCF

        //for (size_t i=0; i < filename_extension_array.size(); i++)
        //{
        //filename_extension = filename_extension_array.at(i);

        if (compute_2D_integrated_3PCF)
        {
            gettimeofday(&start, nullptr);

            if (filename_iB == "iB_kkk.dat" || filename_iB == "iB_Mkk.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zs_bins.size() ; a++)
                {
                    for (size_t b = a; b < zs_bins.size() ; b++)
                    {
                        for (size_t c = b; c < zs_bins.size() ; c++)
                        {
                            assert(corr_idx != num_i3pt_sss_correlations);

                            std::vector<double> iZ;

                            std::stringstream s_iB_Re;
                            s_iB_Re << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;
                            std::string iB_flat_2D_spectrum =  s_iB_Re.str();

                            std::vector<std::vector<double>> iB_Re_matrix = read_3_column_table(iB_flat_2D_spectrum);
                            assert (iB_Re_matrix.at(0).size() == l_array.size());

                            //iZ = xi_theta_array(alpha_table, iB_Re_matrix.at(0), iB_Re_matrix.at(1), "xi", thread_count);

                            iZ = xi_theta_array_bin_averaged(alpha_min_table, alpha_max_table, iB_Re_matrix.at(0), iB_Re_matrix.at(1), "xi", thread_count);

                            std::stringstream s_iZ;
                            s_iZ << correlations_folder << "iZ_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;
                            std::ofstream file_iZ( s_iZ.str(), std::ofstream::trunc);
                            file_iZ.precision(40);


                            for (size_t alpha_idx = 0; alpha_idx < alpha_table.size(); alpha_idx++)
                            {
                                file_iZ << alpha_table.at(alpha_idx) << " " << iZ.at(alpha_idx) << "\n";
                            }

                            file_iZ.close();

                            corr_idx++;
                        }
                    }
                }
            }

            else if (filename_iB == "iB_Mss.dat" || filename_iB == "iB_Mss_angle_averaged.dat")
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zs_bins.size() ; a++)
                {
                    for (size_t b = a; b < zs_bins.size() ; b++)
                    {
                        for (size_t c = b; c < zs_bins.size() ; c++)
                        {
                            assert(corr_idx != num_i3pt_sss_correlations);

                            std::vector<double> iZp_Re, iZm_Re;

                            std::stringstream s_iB_Re;
                            s_iB_Re << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;
                            std::string iB_flat_2D_spectrum =  s_iB_Re.str();

                            std::vector<std::vector<double>> iB_Re_matrix = read_3_column_table(iB_flat_2D_spectrum);
                            assert (iB_Re_matrix.at(0).size() == l_array.size());

                            //iZp_Re = xi_theta_array(alpha_table, iB_Re_matrix.at(0), iB_Re_matrix.at(1), "xip", thread_count);
                            //iZm_Re = xi_theta_array(alpha_table, iB_Re_matrix.at(0), iB_Re_matrix.at(2), "xim", thread_count);

                            iZp_Re = xi_theta_array_bin_averaged(alpha_min_table, alpha_max_table, iB_Re_matrix.at(0), iB_Re_matrix.at(1), "xip", thread_count);
                            iZm_Re = xi_theta_array_bin_averaged(alpha_min_table, alpha_max_table, iB_Re_matrix.at(0), iB_Re_matrix.at(2), "xim", thread_count);

                            std::stringstream s_iZp_Re;
                            s_iZp_Re << correlations_folder << "iZp_Re_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;
                            std::ofstream file_iZp_Re( s_iZp_Re.str(), std::ofstream::trunc);
                            file_iZp_Re.precision(40);

                            std::stringstream s_iZm_Re;
                            s_iZm_Re << correlations_folder << "iZm_Re_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;
                            std::ofstream file_iZm_Re( s_iZm_Re.str(), std::ofstream::trunc);
                            file_iZm_Re.precision(40);

                            for (size_t alpha_idx = 0; alpha_idx < alpha_table.size(); alpha_idx++)
                            {
                                file_iZp_Re << alpha_table.at(alpha_idx) << " " << iZp_Re.at(alpha_idx) << "\n";
                                file_iZm_Re << alpha_table.at(alpha_idx) << " " << iZm_Re.at(alpha_idx) << "\n";
                            }

                            file_iZp_Re.close();
                            file_iZm_Re.close();

                            corr_idx++;
                        }
                    }
                }
            }

            if (filename_iB == "iB_hkk.dat") // with bias
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zl_bins.size() ; a++)
                {
                    for (size_t b = 0; b < zs_bins.size() ; b++)
                    {
                        for (size_t c = b; c < zs_bins.size() ; c++)
                        {
                            std::vector<double> iZ, iZ_component;

                            std::stringstream s_iB;
                            s_iB << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;
                            std::string iB_flat_2D_spectrum =  s_iB.str();

                            std::vector<std::vector<double>> iB_matrix = read_7_column_table(iB_flat_2D_spectrum);

                            assert (iB_matrix.at(0).size() == l_array.size());

                            std::vector<double> iB(iB_matrix.at(0).size(), 0);

                            for (size_t b_idx = 1; b_idx <= iB_matrix.size()/2; b_idx++)
                            {
                                // Each of the bias modelling components
                                std::stringstream s_iZ_component;
                                s_iZ_component << correlations_folder << "iZ_b" << std::to_string(b_idx) << "_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;
                                std::ofstream file_iZ_component( s_iZ_component.str(), std::ofstream::trunc);
                                file_iZ_component.precision(40);

                                //iZ_component = xi_theta_array(alpha_table, iB_matrix.at(0), iB_matrix.at(b_idx), "xi", thread_count);

                                iZ_component = xi_theta_array_bin_averaged(alpha_min_table, alpha_max_table, iB_matrix.at(0), iB_matrix.at(b_idx), "xi", thread_count);

                                for (size_t alpha_idx = 0; alpha_idx < alpha_table.size(); alpha_idx++)
                                {
                                    file_iZ_component << alpha_table.at(alpha_idx) << " " << iZ_component.at(alpha_idx) << "\n";
                                }

                                file_iZ_component.close();
                            }

                            // Total sum of the bias modelling components
                            for (size_t l_idx = 0; l_idx < iB_matrix.at(0).size(); l_idx++)
                            {
                                for (size_t b_idx = 1; b_idx <= iB_matrix.size()/2; b_idx++)
                                {
                                    iB[l_idx] += iB_matrix[b_idx][l_idx];
                                }
                            }

                            //iZ = xi_theta_array(alpha_table, iB_matrix.at(0), iB, "xi", thread_count);

                            iZ = xi_theta_array_bin_averaged(alpha_min_table, alpha_max_table, iB_matrix.at(0), iB, "xi", thread_count);

                            std::stringstream s_iZ;
                            s_iZ << correlations_folder << "iZ_" << std::to_string(a+1) << std::to_string(b+1) << std::to_string(c+1) << filename_extension;
                            std::ofstream file_iZ( s_iZ.str(), std::ofstream::trunc);
                            file_iZ.precision(40);

                            for (size_t alpha_idx = 0; alpha_idx < alpha_table.size(); alpha_idx++)
                            {
                                file_iZ << alpha_table.at(alpha_idx) << " " << iZ.at(alpha_idx) << "\n";
                            }

                            file_iZ.close();

                            corr_idx++;
                        }
                    }
                }
            }

            if (filename_iB == "iB_hhh.dat") // with bias
            {
                size_t corr_idx = 0;
                for (size_t a = 0; a < zl_bins.size() ; a++)
                {
                    std::vector<double> iZ, iZ_component;

                    std::stringstream s_iB;
                    s_iB << spectra_folder << "iB_" << std::to_string(a+1) << std::to_string(a+1) << std::to_string(a+1) << filename_extension;
                    std::string iB_flat_2D_spectrum =  s_iB.str();

                    std::vector<std::vector<double>> iB_matrix = read_6_column_table(iB_flat_2D_spectrum);

                    assert (iB_matrix.at(0).size() == l_array.size());

                    std::vector<double> iB(iB_matrix.at(0).size(), 0);

                    for (size_t b_idx = 1; b_idx < iB_matrix.size(); b_idx++)
                    {
                        // Each of the bias modelling components
                        std::stringstream s_iZ_component;
                        s_iZ_component << correlations_folder << "iZ_b" << std::to_string(b_idx) << "_" << std::to_string(a+1) << std::to_string(a+1) << std::to_string(a+1) << filename_extension;
                        std::ofstream file_iZ_component( s_iZ_component.str(), std::ofstream::trunc);
                        file_iZ_component.precision(40);

                        //iZ_component = xi_theta_array(alpha_table, iB_matrix.at(0), iB_matrix.at(b_idx), "xi", thread_count);

                        iZ_component = xi_theta_array_bin_averaged(alpha_min_table, alpha_max_table, iB_matrix.at(0), iB_matrix.at(b_idx), "xi", thread_count);

                        for (size_t alpha_idx = 0; alpha_idx < alpha_table.size(); alpha_idx++)
                        {
                            file_iZ_component << alpha_table.at(alpha_idx) << " " << iZ_component.at(alpha_idx) << "\n";
                        }

                        file_iZ_component.close();
                    }

                    // --------------------------------------------
                    // large-l (l > 5000) asymptote subtraction for iB_b4 term

                    std::vector<double> iB_b4_minus_large_l_asymptote(iB_matrix.at(4).size(), 0);

                    double iB_b4_large_l_asymptote = 0;
                    int num_l_modes_asymptote = 0;
                    for (size_t l_idx = 0; l_idx < iB_matrix.at(0).size(); l_idx++)
                    {
                        if (iB_matrix[0][l_idx] > 5000)
                        {
                            iB_b4_large_l_asymptote += iB_matrix[4][l_idx];
                            num_l_modes_asymptote++;
                        }
                        else
                            iB_b4_minus_large_l_asymptote[l_idx] =  iB_matrix[4][l_idx];
                    }
                    iB_b4_large_l_asymptote = iB_b4_large_l_asymptote / num_l_modes_asymptote;

                    for (size_t l_idx = 0; l_idx < iB_matrix.at(0).size(); l_idx++)
                    {
                        if (iB_matrix[0][l_idx] <= 5000)
                        {
                            iB_b4_minus_large_l_asymptote[l_idx] -= iB_b4_large_l_asymptote;
                        }
                        else
                            iB_b4_minus_large_l_asymptote[l_idx] = 0.0;
                    }

                    // --------------------------------------------

                    // Total sum of the bias modelling components
                    for (size_t l_idx = 0; l_idx < iB_matrix.at(0).size(); l_idx++)
                    {
                        for (size_t b_idx = 1; b_idx < iB_matrix.size(); b_idx++)
                        {
                            if (b_idx == 4)
                            {
                                iB[l_idx] += iB_b4_minus_large_l_asymptote[l_idx];
                            }

                            else if (b_idx == 5) // don't add Poisson shot noise term as it gives Dirac delta anyway
                                continue;

                            else
                                iB[l_idx] += iB_matrix[b_idx][l_idx];
                        }
                    }

                    //iZ = xi_theta_array(alpha_table, iB_matrix.at(0), iB, "xi", thread_count);

                    iZ = xi_theta_array_bin_averaged(alpha_min_table, alpha_max_table, iB_matrix.at(0), iB, "xi", thread_count);

                    std::stringstream s_iZ;
                    s_iZ << correlations_folder << "iZ_" << std::to_string(a+1) << std::to_string(a+1) << std::to_string(a+1) << filename_extension;
                    std::ofstream file_iZ( s_iZ.str(), std::ofstream::trunc);
                    file_iZ.precision(40);

                    for (size_t alpha_idx = 0; alpha_idx < alpha_table.size(); alpha_idx++)
                    {
                        file_iZ << alpha_table.at(alpha_idx) << " " << iZ.at(alpha_idx) << "\n";
                    }

                    file_iZ.close();

                    corr_idx++;
                }
            }

            gettimeofday(&end, nullptr);
            time_taken = (end.tv_sec - start.tv_sec) * 1e6;
            time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
            std::cout << "\nTime taken for 2D integrated 3PCF calculations: " << time_taken << " sec" << std::endl;

            if (verbose_print_outs)
                std::cout<<"2D integrated 3PCF output files created\n";
        }

        //}

    }

    gettimeofday(&end_file, nullptr);
    time_taken_file = (end_file.tv_sec - start_file.tv_sec) * 1e6;
    time_taken_file = (time_taken_file + (end_file.tv_usec - start_file.tv_usec)) * 1e-6;
    std::cout << "\nTotal time taken for the execution of main.cpp: " << time_taken_file << " sec" << std::endl;

    return 0;
}
