//--------------------------------------------------------------------------
//
// Description:
//  class ClassEngine :
//  encapsulation of class calls
//
//
// Author List:
//	Stephane Plaszczynski (plaszczy@lal.in2p3.fr)
//  Anik Halder
//
// History (add to end):
//	creation:   ven. nov. 4 11:02:20 CET 2011 
//  modified: 2020 onwards
//-----------------------------------------------------------------------

#ifndef ClassEngine_hh
#define ClassEngine_hh

//CLASS
#include"class.h"
#include"Engine.hh"
//STD
#include<string>
#include<vector>
#include<utility>
#include<ostream>

#include <interpolation_methods.h>

#define _rho_class_to_SI_units_ 3.*_c_*_c_/8./_PI_/_G_/_Mpc_over_m_/_Mpc_over_m_
#define _kg_m3_to_M_sun_Mpc3_units_ 1./_M_SUN_*_Mpc_over_m_*_Mpc_over_m_*_Mpc_over_m_

using std::string;

//general utility to convert safely numerical types to string
template<typename T> std::string str(const T &x);
//specialisations
template<> std::string str (const float &x);
template<> std::string str (const double &x);
template<> std::string str (const bool &x); //"yes" or "no"
template<> std::string str (const std::string &x);

std::string str(const char* x);
//////////////////////////////////////////////////////////////////////////
//class to encapsulate CLASS parameters from any type (numerical or string)
class ClassParams{
public:

  ClassParams(){}
  ClassParams( const ClassParams& o):pars(o.pars){}

  //use this to add a CLASS variable
  template<typename T> unsigned add(const string& key,const T& val){
  pars.push_back(make_pair(key,str(val)));
  return pars.size();
  }
  
  //accesors
  inline unsigned size() const {return pars.size();}
  inline string key(const unsigned& i) const {return pars[i].first;}
  inline string value(const unsigned& i) const {return pars[i].second;}


private:
  std::vector<std::pair<string,string> > pars;
};

///////////////////////////////////////////////////////////////////////////
class ClassEngine : public Engine
{

  friend class ClassParams;

public:
  //constructors
  ClassEngine(const ClassParams& pars, bool verbose=true);
  //with a class .pre file
  ClassEngine(const ClassParams& pars, const string & precision_file, bool verbose=true);
  //with a class .ini file and also .pre file
  ClassEngine(int argc, char **argv);

  // destructor
  ~ClassEngine();

  //modfiers: _FAILURE_ returned if CLASS pb:
  bool updateParValues(const std::vector<double>& par);

  //print content of file_content
  void printFC();

  double getCl(Engine::cltype t,const long &l); //get Cl value at l ( 2<l<lmax)

  void getCls(const std::vector<unsigned>& lVec, //input 
	      std::vector<double>& cltt, 
	      std::vector<double>& clte, 
	      std::vector<double>& clee, 
	      std::vector<double>& clbb);
  bool getLensing(const std::vector<unsigned>& lVec, //input 
	      std::vector<double>& clphiphi, 
	      std::vector<double>& cltphi, 
	      std::vector<double>& clephi);

  void call_perturb_sources_at_tau(
                           int index_md,
                           int index_ic,
                           int index_tp,
                           double tau,
                           double * psource
                           );

  void getTk( double z, 
        std::vector<double>& k,
        std::vector<double>& d_cdm,
        std::vector<double>& d_b,
        std::vector<double>& d_ncdm,
        std::vector<double>& d_tot,
        std::vector<double>& t_cdm,
        std::vector<double>& t_b,
        std::vector<double>& t_ncdm,
        std::vector<double>& t_tot );

  //for BAO
  inline double z_drag() const {return th.z_d;}
  inline double rs_drag() const {return th.rs_d;} 
  double getTauReio() const {return th.tau_reio;}
  inline int numCls() const {return sp.ct_size;}
  inline double Tcmb() const {return ba.T_cmb;}
  inline int l_max_scalars() const {return _lmax;}

  double get_a_z(const double &z); // scale (expansion) factor
  double get_tau_z(const double &z); // conformal time [Mpc]
  double get_chi_z(const double &z); // comoving distance [Mpc]
  double get_proper_time_z(const double &z); // proper time [Mpc]

  double get_A_z(const double &z);
  double get_Dv_z(const double &z);
  double get_D_ang_z(const double &z); // angular diameter distance [Mpc]
  double get_D_lum_z(const double &z); // luminosity distance [Mpc]
  double get_F_z(const double &z);
  double get_H_z(const double &z); // in units [c/Mpc]

  double get_D_plus_z(const double &z); // scale independent linear growth factor D+(z) for CDM perturbations
  double get_f_z(const double &z); // velocity growth factor
  double get_sigma8_z(const double &z);
  double get_sigma_R_z(const double &R, const double &z);
  double get_sigma_squared_prime_R_z(const double &R, const double &z);
  double get_sigma_prime_R_z(const double &R, const double &z);

  double get_H0(); // in units [c/Mpc]
  double get_h();
  double get_k_pivot();
  double get_A_s();
  double get_n_s();
  double get_Omega0_m();
  double get_Omega0_b();
  double get_Omega0_cdm();
  double get_Omega0_g();
  double get_Omega0_Lambda();
  double get_Omega0_k();
  double get_Omega0_r();
  double get_Omega0_ur();
  double get_Omega0_de();
  double get_Omega0_fld();
  double get_w0_fld();
  double get_wa_fld();
  double get_k_max_pk();
  double get_z_max_pk();

  double get_Omega_m_z(const double &z);
  double get_Omega_r_z(const double &z);
  double get_Omega_fld_z(const double &z);
  double get_Omega_k_z(const double &z);

  double get_rho_crit_z(const double &z); // in units [kg/m^3]
  double get_rho_m_z(const double &z); // in units [kg/m^3]

  double get_f_NL_local();
  void set_f_NL_local(const double &f_NL_local);

  double get_f_NL_equilateral();
  void set_f_NL_equilateral(const double &f_NL_equilateral);

  double get_f_NL_orthogonal();
  void set_f_NL_orthogonal(const double &f_NL_orthogonal);

  double get_k_NL_from_lin_Pk(const double &z);
  double get_k_NL_from_lin_Pk_interp(const double &z);
  double get_k_NL_from_nl_Pk_class(const double &z);

  double get_n_eff_from_lin_Pk(double k, double z); // logarithmic slope of linear P(k,z)
  double get_n_eff_from_nl_Pk(double k, double z);// logarithmic slope of non-linear P(k,z)

  double get_G_1_k_z_interp(double k, double z); // density growth-only response function G_1(k,z)
  double get_G_K_k_z_interp(double k, double z); // tidal growth-only response function G_K(k,z)

  double get_a_GM_k_z_interp(const double &k, const double &z); 
  double get_b_GM_k_z_interp(const double &k, const double &z);
  double get_c_GM_k_z_interp(const double &k, const double &z); 

  void bihalofit_compute_pk_norm();
  double bihalofit_pk_lin(const double &k_h, const double &z);
  double bihalofit_window(const double &x, const int &j);
  double bihalofit_sigma_j(const double &R, const int &j); // R [Mpc/h]
  double bihalofit_R_NL(const double &z);  // R_NL [Mpc/h]
  double bihalofit_n_eff(const double &z, const double &R_NL);

  double bihalofit_R_NL_interp(const double &z);
  double bihalofit_n_eff_interp(const double &z);

  double pk_primordial(const double &k); // P_R(k)
  double pk_gravitational_potential(const double &k); // P_phi(k)
  double M_poisson_factor(const double &k, const double &z); // Poisson factor relating delta_m to gravitational potential

  double pk_lin(const double &k, const double &z);
  double pk_nl(const double &k, const double &z);
  double pk(const double &k, const double &z, bool use_pk_nl);

private:
  //structures class en commun
  struct file_content fc;
  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */

  ErrorMsg _errmsg;            /* for error messages */
  double * cl;

  //helpers
  bool dofree;
  int freeStructs();

  bool input_init_already_done;

  //call once /model
  int computeCls();

  void compute_bispectrum_helpers();

  int class_main(
		 struct file_content *pfc,
		 struct precision * ppr,
		 struct background * pba,
		 struct thermo * pth,
		 struct perturbs * ppt,
		 struct transfers * ptr,
		 struct primordial * ppm,
		 struct spectra * psp,
		 struct nonlinear * pnl,
		 struct lensing * ple,
		 struct output * pop,
		 ErrorMsg errmsg);
  //parnames
  std::vector<std::string> parNames;

  Linear_interp_1D k_NL_z_from_lin_Pk_array;  /* k_NL of linear P(k,z) interpolation array */
  Linear_interp_2D G_1_k_z;  /* density growth-only response function G_1(k,z) */
  Linear_interp_2D G_K_k_z;  /* tidal growth-only response function G_K(k,z) */

  Linear_interp_2D a_GM_k_z;  
  Linear_interp_2D b_GM_k_z; 
  Linear_interp_2D c_GM_k_z; 

  double m_bihalofit_norm;

  Linear_interp_1D bihalofit_R_NL_z_array;  /* R_NL(z) bihalofit interpolation array */
  Linear_interp_1D bihalofit_n_eff_z_array;  /* n_eff(z) bihalofit interpolation array */

  // primordial non-Gaussianity

  double m_f_NL_local;
  double m_f_NL_equilateral;
  double m_f_NL_orthogonal;

protected:
 
};

#endif
