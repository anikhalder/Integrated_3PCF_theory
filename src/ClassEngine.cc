//--------------------------------------------------------------------------
//
// Description:
// 	class ClassEngine : see header file (ClassEngine.hh) for description.
//
//------------------------------------------------------------------------
//-----------------------
// This Class's Header --
//-----------------------
#include "ClassEngine.hh"
// ------------------------
// Collaborating classes --
//-------------------------
// C++
//--------------------
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <numeric>
#include <cassert>

#include <cosmology_utils.h>
#include <constants.h>
#include <bispectrum.hpp>

#include <sys/time.h>

//#define DBUG

using namespace std;

template<typename T> std::string str(const T &x)
{
  std::ostringstream os;
  os << x;
  return os.str();
}
//specilization
template<> std::string str (const float &x)
{
  std::ostringstream os;
  os << setprecision(8) << x;
  return os.str();
}
template<> std::string str (const double &x)
{
  std::ostringstream os;
  os << setprecision(16) << x;
  return os.str();
}
template<> std::string str (const bool &x)
{
  { return x ? "yes" : "no"; }
}

template<> std::string str (const std::string &x) {return x;}

std::string str (const char* s){return string(s);}

//instanciations
template string str(const int &x);
template string str(const signed char &x);
template string str(const unsigned char &x);
template string str(const short &x);
template string str(const unsigned short &x);
template string str(const unsigned int &x);
template string str(const long &x);
template string str(const unsigned long &x);
template string str(const long long &x);
template string str(const unsigned long long &x);

//---------------
// Constructors --
//----------------
ClassEngine::ClassEngine(const ClassParams& pars, bool verbose): cl(0),dofree(true)
{
  // struct timeval start, end;
  // double time_taken;
  // gettimeofday(&start, nullptr);
  
  //prepare fp structure
  size_t n=pars.size();
  //
  parser_init(&fc,n,"pipo",_errmsg);
  
  //config
  for (size_t i=0;i<pars.size();i++){
    strcpy(fc.name[i],pars.key(i).c_str());
    strcpy(fc.value[i],pars.value(i).c_str());
    //store
    parNames.push_back(pars.key(i));
    //identify lmax
    if(verbose) cout << pars.key(i) << "\t" << pars.value(i) <<endl;
    if (pars.key(i)=="l_max_scalars") {
      istringstream strstrm(pars.value(i));
      strstrm >> _lmax;
    }
  }
  if( verbose ) cout << __FILE__ << " : using lmax=" << _lmax <<endl;
  // assert(_lmax>0); // this collides with transfer function calculations

  //input
  if (input_init(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,_errmsg) == _FAILURE_) 
    throw invalid_argument(_errmsg);

  input_init_already_done = true;

  //proetction parametres mal defini
  for (size_t i=0;i<pars.size();i++){
    if (fc.read[i] !=_TRUE_) throw invalid_argument(string("invalid CLASS parameter: ")+fc.name[i]);
  }

  // gettimeofday(&end, nullptr);
  // time_taken = (end.tv_sec - start.tv_sec) * 1e6;
  // time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
  // std::cout << "Time taken for ClassEngine to do reading: " << time_taken << " sec" << std::endl;

  //calcul class
  computeCls(); // This seems to take the major amount of time
  
  //cout <<"creating " << sp.ct_size << " arrays" <<endl;
  if( pt.has_cl_cmb_temperature || pt.has_cl_cmb_polarization || pt.has_cl_lensing_potential ){
    cl=new double[sp.ct_size];
  }

  //printFC();

  compute_bispectrum_helpers();

  m_f_NL_local = 0;
  m_f_NL_equilateral = 0;
  m_f_NL_orthogonal = 0;

  // gettimeofday(&end, nullptr);
  // time_taken = (end.tv_sec - start.tv_sec) * 1e6;
  // time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
  // std::cout << "Time taken for ClassEngine to do reading and all computations: " << time_taken << " sec" << std::endl;
}

ClassEngine::ClassEngine(const ClassParams& pars,const string & precision_file, bool verbose): cl(0),dofree(true)
{
  struct file_content fc_precision;
  fc_precision.size = 0;
  //decode pre structure
  if (parser_read_file(const_cast<char*>(precision_file.c_str()),&fc_precision,_errmsg) == _FAILURE_){
    throw invalid_argument(_errmsg);
  }

  //pars
  struct file_content fc_input;
  fc_input.size = 0;
  fc_input.filename=new char[1];
 //prepare fc par structure
  size_t n=pars.size();
  parser_init(&fc_input,n,"pipo",_errmsg);
  //config
  for (size_t i=0;i<pars.size();i++){
    strcpy(fc_input.name[i],pars.key(i).c_str());
    strcpy(fc_input.value[i],pars.value(i).c_str());
    if (pars.key(i)=="l_max_scalars") {
      istringstream strstrm(pars.value(i));
      strstrm >> _lmax;
    }
  }
  if( verbose ) cout << __FILE__ << " : using lmax=" << _lmax <<endl;
  //assert(_lmax>0);
  
  //concatenate both
  if (parser_cat(&fc_input,&fc_precision,&fc,_errmsg) == _FAILURE_) throw invalid_argument(_errmsg);

  //parser_free(&fc_input);
  parser_free(&fc_precision);
  
  //input
  if (input_init(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,_errmsg) == _FAILURE_) 
    throw invalid_argument(_errmsg);

  input_init_already_done = true;

  //proetction parametres mal defini
  for (size_t i=0;i<pars.size();i++){
    if (fc.read[i] !=_TRUE_) throw invalid_argument(string("invalid CLASS parameter: ")+fc.name[i]);
  }

  //calcul class
  computeCls();
  
  //cout <<"creating " << sp.ct_size << " arrays" <<endl;
  if( pt.has_cl_cmb_temperature || pt.has_cl_cmb_polarization || pt.has_cl_lensing_potential ){
    cl=new double[sp.ct_size];
  }
  //printFC();

  compute_bispectrum_helpers();

  m_f_NL_local = 0;
  m_f_NL_equilateral = 0;
  m_f_NL_orthogonal = 0;
}

ClassEngine::ClassEngine(int argc, char **argv): cl(0),dofree(true)
{
  // struct timeval start, end;
  // double time_taken;
  // gettimeofday(&start, nullptr);

  //input
  if (input_init_from_arguments(argc,argv,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,_errmsg) == _FAILURE_) 
    throw invalid_argument(_errmsg);

  input_init_already_done = true;

  // gettimeofday(&end, nullptr);
  // time_taken = (end.tv_sec - start.tv_sec) * 1e6;
  // time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
  // std::cout << "Time taken for ClassEngine to do reading: " << time_taken << " sec" << std::endl;

  //calcul class
  computeCls();
  
  //cout <<"creating " << sp.ct_size << " arrays" <<endl;
  if( pt.has_cl_cmb_temperature || pt.has_cl_cmb_polarization || pt.has_cl_lensing_potential ){
    cl=new double[sp.ct_size];
  }
  //printFC();

  compute_bispectrum_helpers();

  m_f_NL_local = 0;
  m_f_NL_equilateral = 0;
  m_f_NL_orthogonal = 0;

  // gettimeofday(&end, nullptr);
  // time_taken = (end.tv_sec - start.tv_sec) * 1e6;
  // time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
  // std::cout << "Time taken for ClassEngine to do reading and all computations: " << time_taken << " sec" << std::endl;
}
  
//--------------
// Destructor --
//--------------
ClassEngine::~ClassEngine()
{
  //printFC();
  dofree && freeStructs();

  if( pt.has_cl_cmb_temperature || pt.has_cl_cmb_polarization || pt.has_cl_lensing_potential ){
    delete [] cl;
  }

  input_init_already_done = false;
}

//-----------------
// Member functions --
//-----------------
bool ClassEngine::updateParValues(const std::vector<double>& par)
{
  dofree && freeStructs();
  for (size_t i=0;i<par.size();i++) {
    double val=par[i];
    strcpy(fc.value[i],str(val).c_str());
    strcpy(fc.name[i],parNames[i].c_str());
  
  #ifdef DBUG
    cout << "update par values #" << i << "\t" <<  val << "\t" << str(val).c_str() << endl;
  #endif
  }

  input_init_already_done = false;
  
  int status=computeCls();
  #ifdef DBUG
    cout << "update par status=" << status << " succes=" << _SUCCESS_ << endl;
  #endif

  return (status==_SUCCESS_);
}

//print content of file_content
void ClassEngine::printFC()
{
  printf("FILE_CONTENT SIZE=%d\n",fc.size);
  for (int i=0;i<fc.size;i++) printf("%d : %s = %s\n",i,fc.name[i],fc.value[i]);
}

int ClassEngine::class_main(
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
                ErrorMsg errmsg)
{
  if (input_init_already_done == false)
  {
    if (input_init(pfc,ppr,pba,pth,ppt,ptr,ppm,psp,pnl,ple,pop,errmsg) == _FAILURE_) {
      printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
      dofree=false;
      return _FAILURE_;
    }
    input_init_already_done = true;
  }

  if (background_init(ppr,pba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",pba->error_message);
    dofree=false;
    return _FAILURE_;
  }

  if (thermodynamics_init(ppr,pba,pth) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",pth->error_message);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (perturb_init(ppr,pba,pth,ppt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",ppt->error_message);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (primordial_init(ppr,ppt,ppm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",ppm->error_message);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (nonlinear_init(ppr,pba,pth,ppt,ppm,pnl) == _FAILURE_)  {
    printf("\n\nError in nonlinear_init \n=>%s\n",pnl->error_message);
    primordial_free(&pm);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (transfer_init(ppr,pba,pth,ppt,pnl,ptr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",ptr->error_message);
    nonlinear_free(&nl);
    primordial_free(&pm);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (spectra_init(ppr,pba,ppt,ppm,pnl,ptr,psp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",psp->error_message);
    transfer_free(&tr);
    nonlinear_free(&nl);
    primordial_free(&pm);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (lensing_init(ppr,ppt,psp,pnl,ple) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",ple->error_message);
    spectra_free(&sp);
    transfer_free(&tr);
    nonlinear_free(&nl);
    primordial_free(&pm);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }


  dofree=true;
  return _SUCCESS_;
}


int ClassEngine::freeStructs()
{
  if (lensing_free(&le) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  if (spectra_free(&sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (perturb_free(&pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;
}

int ClassEngine::computeCls()
{

#ifdef DBUG
  cout <<"call computecls" << endl;
  //printFC();
#endif

  int status=this->class_main(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,_errmsg);
#ifdef DBUG
  cout <<"status=" << status << endl;
#endif
  return status;

}

void ClassEngine::compute_bispectrum_helpers()
{
  std::vector<double> m_z_array;
  std::vector<double> m_k_NL_z_array;
  std::vector<double> m_bihalofit_R_NL_z_array;
  std::vector<double> m_bihalofit_n_eff_z_array;

  if (compute_bihalofit)
    bihalofit_compute_pk_norm();

  double z=0.;
  while (z <= get_z_max_pk())
  {
    m_z_array.push_back(z);
    z += 0.05;
    //z += delta_z_step;
  }

  for (size_t z_idx=0; z_idx<m_z_array.size(); z_idx++)
  {
    z = m_z_array.at(z_idx);
    m_k_NL_z_array.push_back(get_k_NL_from_lin_Pk(z));

    double R_NL;

    if (compute_bihalofit)
    {
        R_NL = bihalofit_R_NL(z);
        m_bihalofit_R_NL_z_array.push_back(R_NL);
        m_bihalofit_n_eff_z_array.push_back(bihalofit_n_eff(z, R_NL));

        if (verbose_print_outs)
          std::cout << z << " " << get_Omega_m_z(z) << " " << R_NL << std::endl;
    }
    else
    {
        R_NL = 0;
        m_bihalofit_R_NL_z_array.push_back(0);
        m_bihalofit_n_eff_z_array.push_back(0);
    }
  }

  k_NL_z_from_lin_Pk_array = Linear_interp_1D(m_z_array, m_k_NL_z_array);
  bihalofit_R_NL_z_array = Linear_interp_1D(m_z_array, m_bihalofit_R_NL_z_array);
  bihalofit_n_eff_z_array = Linear_interp_1D(m_z_array, m_bihalofit_n_eff_z_array);

  std::vector<double> G_1_k_array = read_1_column_table("../data/response_functions/G_1_k_h.tab"); // in h/Mpc
  std::vector<double> G_1_z_array = read_1_column_table("../data/response_functions/G_1_z.tab");
  std::vector<double> G_1_vals_array = read_1_column_table("../data/response_functions/G_1_vals.tab");

  for (size_t i=0; i<G_1_k_array.size(); i++)
      G_1_k_array[i] *=  ba.h; // convert to 1/Mpc

  G_1_k_z = Linear_interp_2D(G_1_k_array, G_1_z_array, G_1_vals_array);

  std::vector<double> G_K_k_array = read_1_column_table("../data/response_functions/G_K_k_h.tab"); // in h/Mpc
  std::vector<double> G_K_z_array = read_1_column_table("../data/response_functions/G_K_z.tab");
  std::vector<double> G_K_vals_array = read_1_column_table("../data/response_functions/G_K_vals.tab");

  for (size_t i=0; i<G_K_k_array.size(); i++)
      G_K_k_array[i] *=  ba.h;// convert to 1/Mpc

  G_K_k_z = Linear_interp_2D(G_K_k_array, G_K_z_array, G_K_vals_array);

  // pre-computing GM bispectrum fitting functions
  std::vector<double> m_k_array;
  double a=-5,b=log10(get_k_max_pk());

  for (int i=0; i<num_k_pts; i++)
    m_k_array.push_back(pow(10, a + i * (b - a) / (num_k_pts - 1))); // in units of 1/Mpc

  std::vector<double> a_GM_k_z_vals(m_k_array.size()*m_z_array.size(), 0);
  std::vector<double> b_GM_k_z_vals(m_k_array.size()*m_z_array.size(), 0);
  std::vector<double> c_GM_k_z_vals(m_k_array.size()*m_z_array.size(), 0);

  for (size_t z_idx=0; z_idx<m_z_array.size(); z_idx++)
  {
    z = m_z_array.at(z_idx);
    double k_NL_z = get_k_NL_from_lin_Pk(z);
    double sigma8_z = get_sigma8_z(z); 

    for (size_t k_idx=0; k_idx<m_k_array.size(); k_idx++)
    {
      double k = m_k_array[k_idx];
      double n = get_n_eff_from_lin_Pk(k, 0);

      double q = k/k_NL_z;

      a_GM_k_z_vals[k_idx*m_z_array.size()+z_idx] = a_GM(n,q,sigma8_z);
      b_GM_k_z_vals[k_idx*m_z_array.size()+z_idx] = b_GM(n,q);
      c_GM_k_z_vals[k_idx*m_z_array.size()+z_idx] = c_GM(n,q);
    }
  }

  a_GM_k_z = Linear_interp_2D(m_k_array, m_z_array, a_GM_k_z_vals);
  b_GM_k_z = Linear_interp_2D(m_k_array, m_z_array, b_GM_k_z_vals);
  c_GM_k_z = Linear_interp_2D(m_k_array, m_z_array, c_GM_k_z_vals);
}

void ClassEngine::call_perturb_sources_at_tau(
                           int index_md,
                           int index_ic,
                           int index_tp,
                           double tau,
                           double * psource
                           )
{
  if( perturb_sources_at_tau( &pt, index_md, index_ic, index_tp, tau, psource ) == _FAILURE_){
    cerr << ">>>fail getting Tk type=" << (int)index_tp <<endl; 
    throw out_of_range(pt.error_message);
  }
}

void ClassEngine::getTk( double z, 
   std::vector<double>& k,
   std::vector<double>& d_cdm,
   std::vector<double>& d_b,
   std::vector<double>& d_ncdm,
   std::vector<double>& d_tot,
   std::vector<double>& t_cdm,
   std::vector<double>& t_b,
   std::vector<double>& t_ncdm,
   std::vector<double>& t_tot )
{
  
  if (!dofree) throw out_of_range("no Tk available because CLASS failed");

  double tau;
  int index;
  //transform redshift to conformal time
  background_tau_of_z(&ba,z,&tau);

  if(log(tau) < pt.ln_tau[0]){
    cerr << "Asking sources at a z bigger than z_max_pk, something probably went wrong\n";
    throw out_of_range(pt.error_message);
  }

  double *pvecback=new double[ba.bg_size];
  background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);
  double fHa = pvecback[ba.index_bg_f] * (pvecback[ba.index_bg_a]*pvecback[ba.index_bg_H]);
  delete[] pvecback;
  
  // copy transfer func data to temporary
  const size_t index_md = pt.index_md_scalars;
  d_cdm.assign( pt.k_size[index_md], 0.0 );
  d_b.assign( pt.k_size[index_md], 0.0 );
  d_ncdm.assign( pt.k_size[index_md], 0.0 );
  d_tot.assign( pt.k_size[index_md], 0.0 );
  t_cdm.assign( pt.k_size[index_md], 0.0 );
  t_b.assign( pt.k_size[index_md], 0.0 );
  t_ncdm.assign( pt.k_size[index_md], 0.0 );
  t_tot.assign( pt.k_size[index_md], 0.0 );

  if( pt.ic_size[index_md] > 1 ){
    cerr << ">>>have more than 1 ICs, will use first and ignore others" << endl;
  }

  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_delta_cdm, tau, &d_cdm[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_delta_b, tau, &d_b[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_delta_ncdm1, tau, &d_ncdm[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_delta_tot, tau, &d_tot[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_theta_cdm, tau, &t_cdm[0]); // THIS LINE WAS MISSING!!!
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_theta_b, tau, &t_b[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_theta_ncdm1, tau, &t_ncdm[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_theta_tot, tau, &t_tot[0]);

  //
  std::vector<double> h_prime(pt.k_size[index_md],0.0), eta_prime(pt.k_size[index_md],0.0);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_eta_prime, tau, &eta_prime[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_h_prime, tau, &h_prime[0]);

  // gauge trafo velocities, store k-vector
  for (int index_k=0; index_k<pt.k_size[index_md]; index_k++) 
  {
    auto ak = pt.k[index_md][index_k];

    // write data to vectors
    k.push_back( ak );

    // use the conformal Newtonian gauge for velocities
    // not correct, but N-body gauge currently not implemented
    double alphak2 = (h_prime[index_k]+6*eta_prime[index_k])/2; 

    t_cdm[index_k]  = (-alphak2) / fHa;
    t_b[index_k]    = (-alphak2 + t_b[index_k]) / fHa;
    t_ncdm[index_k] = (-alphak2 + t_ncdm[index_k]) / fHa;
    t_tot[index_k]  = (-alphak2 + t_tot[index_k]) / fHa;
  }
}

double ClassEngine::getCl(Engine::cltype t,const long &l)
{
  //get Cl value at l ( 2<l<lmax): in units = (micro-K)^2
  //don't call if FAILURE returned previously
  //throws std::execption if pb

  if (!dofree) throw out_of_range("no Cl available because CLASS failed");

  if (output_total_cl_at_l(&sp,&le,&op,static_cast<double>(l),cl) == _FAILURE_){
    cerr << ">>>fail getting Cl type=" << (int)t << " @l=" << l <<endl; 
    throw out_of_range(sp.error_message);
  }

  double zecl=-1;

  double tomuk=1e6*Tcmb();
  double tomuk2=tomuk*tomuk;

  switch(t)
    {
    case TT:
      (sp.has_tt==_TRUE_) ? zecl=tomuk2*cl[sp.index_ct_tt] : throw invalid_argument("no ClTT available");
      break;
    case TE:
      (sp.has_te==_TRUE_) ? zecl=tomuk2*cl[sp.index_ct_te] : throw invalid_argument("no ClTE available");
      break; 
    case EE:
      (sp.has_ee==_TRUE_) ? zecl=tomuk2*cl[sp.index_ct_ee] : throw invalid_argument("no ClEE available");
      break;
    case BB:
      (sp.has_bb==_TRUE_) ? zecl=tomuk2*cl[sp.index_ct_bb] : throw invalid_argument("no ClBB available");
      break;
    case PP:
      (sp.has_pp==_TRUE_) ? zecl=cl[sp.index_ct_pp] : throw invalid_argument("no ClPhi-Phi available");
      break;
    case TP:
      (sp.has_tp==_TRUE_) ? zecl=tomuk*cl[sp.index_ct_tp] : throw invalid_argument("no ClT-Phi available");
      break;
    case EP:
      (sp.has_ep==_TRUE_) ? zecl=tomuk*cl[sp.index_ct_ep] : throw invalid_argument("no ClE-Phi available");
      break;
    }
  
  return zecl;

}

void ClassEngine::getCls(const std::vector<unsigned>& lvec, //input
		      std::vector<double>& cltt, 
		      std::vector<double>& clte, 
		      std::vector<double>& clee, 
		      std::vector<double>& clbb)
{
  cltt.resize(lvec.size());
  clte.resize(lvec.size());
  clee.resize(lvec.size());
  clbb.resize(lvec.size());
  
  for (size_t i=0;i<lvec.size();i++){
    try{
      cltt[i]=getCl(ClassEngine::TT,lvec[i]);
      clte[i]=getCl(ClassEngine::TE,lvec[i]);
      clee[i]=getCl(ClassEngine::EE,lvec[i]);
      clbb[i]=getCl(ClassEngine::BB,lvec[i]);
    }
    catch(exception &e){
      throw e;
    }
  }

}
 
bool ClassEngine::getLensing(const std::vector<unsigned>& lvec, //input 
		std::vector<double>& clpp    , 
		std::vector<double>& cltp  , 
        std::vector<double>& clep  )
{
 
  clpp.resize(lvec.size());
  cltp.resize(lvec.size());
  clep.resize(lvec.size());
  
  for (size_t i=0;i<lvec.size();i++){
    try{
      clpp[i]=getCl(ClassEngine::PP,lvec[i]);
      cltp[i]=getCl(ClassEngine::TP,lvec[i]);
      clep[i]=getCl(ClassEngine::EP,lvec[i]);
    }
    catch(exception &e){
      cout << "plantage!" << endl;
      cout << __FILE__ << e.what() << endl;
      return false;
    }
  }
  return true;
}

double ClassEngine::get_a_z(const double &z)
{
    // scale factor of the Universe at z
    return 1./(1.+z);
}

double ClassEngine::get_tau_z(const double &z)
{
    // conformal time (a.k.a eta) or the comoving horizon of the Universe at z ==> comoving distance light has travelled unimpeded since t=0 to t(z)
    double tau;
    background_tau_of_z(&ba,z,&tau); //transform redshift to conformal time
    return tau; // in [Mpc]
}

double ClassEngine::get_chi_z(const double &z)
{
    // comoving distance of a source at z from us ==> comoving distance light has travelled unimpeded since t(z) to us at t(z=0)
    double tau = get_tau_z(z);
    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    int index;
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

    double chi_z = pvecback[ba.index_bg_conf_distance];
    free(pvecback);

    return chi_z; // in [Mpc]

    // alternative way of calculation using eqn (2.91) of Dodelson, Schmidt Modern Cosmology (2020)
    //return get_tau_z(0) - get_tau_z(z); // in [Mpc]
}

double ClassEngine::get_proper_time_z(const double &z)
{
    double tau = get_tau_z(z);
    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    int index;
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

    double proper_time = pvecback[ba.index_bg_time];
    free(pvecback);

    return proper_time; // in [Mpc] ==> age of the Universe at a given redshift in [Mpc]
}

// -------------------------
// ATTENTION FONCTION BIDON - GET omegam
double ClassEngine::get_A_z(const double &z)
{
  double Dv = get_Dv_z(z);
  // A(z)=100DV(z)sqrt(~mh2)/cz
  double omega_bidon = 0.12 ;
  double Az = 100.*Dv*sqrt(omega_bidon)/(3.e8*z); // is there speed of light somewhere ? 
  return Az;
}
// --------------------------

double ClassEngine::get_Dv_z(const double &z)
{
    double tau = get_tau_z(z);
    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    int index;
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

    double H_z = pvecback[ba.index_bg_H];
    double D_ang = pvecback[ba.index_bg_ang_distance];

    double D_v;
    D_v = pow(D_ang*(1+z),2)*z/H_z;
    D_v = pow(D_v,1./3.);

    free(pvecback);

    return D_v;
}

double ClassEngine::get_D_ang_z(const double &z)
{
    double tau = get_tau_z(z);
    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    int index;
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

    double D_ang = pvecback[ba.index_bg_ang_distance];
    free(pvecback);

    return D_ang; // angular diameter distance in [Mpc]. e.g. D_ang(z) = a(z)*chi(z) = chi(z)/(1.+z) in a spatially flat universe (K = 0) --> also the physical distance
}

double ClassEngine::get_D_lum_z(const double &z)
{
    double tau = get_tau_z(z);
    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    int index;
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

    double D_lum = pvecback[ba.index_bg_lum_distance];
    free(pvecback);

    return D_lum; // luminosity distance in [Mpc]
}

double ClassEngine::get_F_z(const double &z)
{
    double tau = get_tau_z(z);
    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    int index;
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

    double H_z = pvecback[ba.index_bg_H];
    double D_ang = pvecback[ba.index_bg_ang_distance];

    double F_z = (1.+z) * D_ang * H_z /(3.e8) ; // is there speed of light somewhere ?

    free(pvecback);

    return F_z;
}

double ClassEngine::get_H_z(const double &z)
{
    double tau = get_tau_z(z);
    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    int index;
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

    double H_z = pvecback[ba.index_bg_H];
    free(pvecback);

    return(H_z); // in units [c/Mpc] ==> get_H_z(z) * c /1.e3 will give the ACTUAL H(z) in [km/s/Mpc]
}

double ClassEngine::get_D_plus_z(const double &z)
{
    double tau = get_tau_z(z);
    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    int index;
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

    double D_plus_z = pvecback[ba.index_bg_D];
    free(pvecback);

    return D_plus_z; // scale independent growth factor D_+(a) for CDM perturbations
}

double ClassEngine::get_f_z(const double &z)
{
    double tau = get_tau_z(z);
    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    int index;
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

    double f_z = pvecback[ba.index_bg_f];
    free(pvecback);

    return f_z; // velocity growth factor dln D / dln a
}

double ClassEngine::get_sigma8_z(const double &z)
{
    double sigma8_z;
    nonlinear_sigmas_at_z(&pr, &ba, &nl, 8./ba.h, z, 0, out_sigma, &sigma8_z);

    return sigma8_z;

    // Alternative way to calculate sigma8(z) --> gives slightly different result than nonlinear_sigmas_at_z
    //return get_D_plus_z(z)*(*nl.sigma8);

//  double tau;
//  int index;
//  double *pvecback;
//  double sigma8 = 0.;
//  //transform redshift in conformal time
//  background_tau_of_z(&ba,z,&tau);

//  //pvecback must be allocated
//  pvecback=(double *)malloc(ba.bg_size*sizeof(double));

//  //call to fill pvecback
//  background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);
//  //background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback);

//  nonlinear_sigma_at_z(&ba,&nl,8./ba.h,z, 0, 10, &sigma8); // Outdated --> as suggested in class source code nonlinear.c rather use nonlinear_sigmas_at_z
//  //spectra_sigma(&ba,&pm,&sp,8./ba.h,z,&sigma8);

//#ifdef DBUG
//  cout << "sigma_8= "<< sigma8 <<endl;
//#endif

//  free(pvecback);

    //  return sigma8;
}

double ClassEngine::get_sigma_R_z(const double &R, const double &z)
{
    // R must be input in [Mpc]
    double sigma_R_z;
    nonlinear_sigmas_at_z(&pr, &ba, &nl, R, z, 0, out_sigma, &sigma_R_z);

    return sigma_R_z;
}

double ClassEngine::get_sigma_squared_prime_R_z(const double &R, const double &z)
{
    // R must be input in [Mpc]
    double sigma_squared_prime_R_z;
    nonlinear_sigmas_at_z(&pr, &ba, &nl, R, z, 0, out_sigma_prime, &sigma_squared_prime_R_z);

    // NOTE: out_sigma_prime corresponds to d sigma^2 / d R
    return sigma_squared_prime_R_z; // in [1/Mpc]
}

double ClassEngine::get_sigma_prime_R_z(const double &R, const double &z)
{
    // R must be input in [Mpc]
    // d sigma / d R = 1/(2 sigma)* d (sigma^2) / d R
    return get_sigma_squared_prime_R_z(R, z) / (2.*get_sigma_R_z(R,z)); // in [1/Mpc]

    // alternative way of calculation (crude numerical central difference derivative)
    //return (get_sigma_R_z(R*1.01, z) - get_sigma_R_z(R*0.99, z)) / (R*0.02); // in [1/Mpc]
}

double ClassEngine::get_H0()
{
    return ba.H0; // in units [c/Mpc] ==> get_H0 * c /.1e3 will give the ACTUAL H0 in [km/s/Mpc]
}

double ClassEngine::get_h()
{
    // reduced Hubble parameter
    return ba.h;
}

double ClassEngine::get_k_pivot()
{
    // k-pivot of primordial power spectrum [1/Mpc]
    return pm.k_pivot;
}

double ClassEngine::get_A_s()
{
     // amplitude of primordial power spectrum
    return pm.A_s;
}

double ClassEngine::get_n_s()
{
    // scalar spectral index of the primordial power spectrum (aka tilt)
    return pm.n_s;
}

double ClassEngine::get_Omega0_m()
{
    // total non relativistic matter = cdm + baryons + etc. (e.g. non-rel. neutrinos...)
    return ba.Omega0_m;
}

double ClassEngine::get_Omega0_b()
{
    // baryons
    return ba.Omega0_b;
}

double ClassEngine::get_Omega0_cdm()
{
    // cdm
    return ba.Omega0_cdm;
}

double ClassEngine::get_Omega0_g()
{
    // photon density fraction
    return ba.Omega0_g;
}

double ClassEngine::get_Omega0_Lambda()
{
    // cosmological constant i.e. Lambda defined for eqn of state w0=-1
    return ba.Omega0_lambda;
}

double ClassEngine::get_Omega0_k()
{
    // curvature
    return ba.Omega0_k;
}

double ClassEngine::get_Omega0_r()
{
    // total relativistic density fraction --> Omega0_g + Omega0_ur
    return ba.Omega0_r;
}

double ClassEngine::get_Omega0_ur()
{
    // ultra-relativistic species / massless neutrino density fraction
    return ba.Omega0_ur;
}

double ClassEngine::get_Omega0_de()
{
    // dark energy
    return ba.Omega0_de;
}

double ClassEngine::get_Omega0_fld()
{
    // fluid
    return ba.Omega0_fld;
}

double ClassEngine::get_k_max_pk()
{
    return pt.k_max_for_pk;
}

double ClassEngine::get_z_max_pk()
{
    return pt.z_max_pk;
}

double ClassEngine::get_w0_fld()
{
    // current fluid equation of state parameter w0
    return ba.w0_fld;
}

double ClassEngine::get_wa_fld()
{
    // current fluid equation of state parameter wa
    return ba.wa_fld;
}

double ClassEngine::get_Omega_m_z(const double &z)
{
    // total non relativistic matter (m) = cdm + baryons + etc. (e.g. non-rel. neutrinos...)
    double tau = get_tau_z(z);
    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    int index;
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

    double Omega_m_z = pvecback[ba.index_bg_Omega_m];
    free(pvecback);

    return(Omega_m_z);
}

double ClassEngine::get_Omega_r_z(const double &z)
{
    // total relativistic matter (r)
    double tau = get_tau_z(z);
    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    int index;
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

    double Omega_r_z = pvecback[ba.index_bg_Omega_r];
    free(pvecback);

    return(Omega_r_z);
}

double ClassEngine::get_Omega_fld_z(const double &z)
{
    // total relativistic matter (r)
    double tau = get_tau_z(z);
    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    int index;
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

    double Omega_fld_z = pvecback[ba.index_bg_rho_fld]/pvecback[ba.index_bg_rho_crit];
    free(pvecback);

    return(Omega_fld_z);
}

double ClassEngine::get_Omega_k_z(const double &z)
{
   return 1-get_Omega_m_z(z)-get_Omega_r_z(z)-get_Omega_fld_z(z);
}

double ClassEngine::get_rho_crit_z(const double &z)
{
//    double tau = get_tau_z(z);
//    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
//    int index;
//    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

//    double rho_crit_z = pvecback[ba.index_bg_rho_crit]; // rho_crit_z = get_H_z(z)^2 = H(z)^2/c^2 (this is the definition used in class)
//    free(pvecback);

//    return(rho_crit_z)*_rho_class_to_SI_units_; // = 3H(z)^2/(8*pi*G) in units [kg/m^3]

    // alternative way of calculation
    return get_H_z(z)*get_H_z(z)*_rho_class_to_SI_units_; // = 3H(z)^2/(8*pi*G) in units [kg/m^3]
}

double ClassEngine::get_rho_m_z(const double &z)
{
    // total non relativistic matter (m) = cdm + baryons + etc. (e.g. non-rel. neutrinos...)
    return get_Omega_m_z(z)*get_rho_crit_z(z); // in units [kg/m^3]

//    double tau = get_tau_z(z);
//    double *pvecback=(double *)malloc(ba.bg_size*sizeof(double));
//    int index;
//    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback); //call to fill pvecback

//    double rho_m_z = pvecback[ba.index_bg_rho_cdm] + pvecback[ba.index_bg_rho_b];
//    free(pvecback);

//    return(rho_m_z)*_rho_class_to_SI_units_; // in units [kg/m^3]
}

double ClassEngine::get_f_NL_local()
{
    return m_f_NL_local;
}

void ClassEngine::set_f_NL_local(const double &f_NL_local)
{
    m_f_NL_local = f_NL_local;
}

double ClassEngine::get_f_NL_equilateral()
{
    return m_f_NL_equilateral;
}

void ClassEngine::set_f_NL_equilateral(const double &f_NL_equilateral)
{
    m_f_NL_equilateral = f_NL_equilateral;
}

double ClassEngine::get_f_NL_orthogonal()
{
    return m_f_NL_orthogonal;
}

void ClassEngine::set_f_NL_orthogonal(const double &f_NL_orthogonal)
{
    m_f_NL_orthogonal = f_NL_orthogonal;
}

namespace
{
    const double inv_2PI2 = 1/(2.0*M_PI*M_PI);
}

double ClassEngine::get_k_NL_from_lin_Pk(const double &z)
{
    double k_NL, k1, k2, val_k_NL;

    k1 = k2 = 1;

    while ( k1*k1*k1*pk_lin(k1,z)*inv_2PI2 > 1 )
        k1 *= 0.5;

    while ( k2*k2*k2*pk_lin(k2,z)*inv_2PI2  < 1 )
        k2 *= 2.0;

    do
    {
        k_NL = 0.5*(k1+k2);

        val_k_NL = k_NL*k_NL*k_NL*pk_lin(k_NL,z)*inv_2PI2 ;

        if (val_k_NL < 1)
            k1 = k_NL;

        else if (val_k_NL > 1)
            k2 = k_NL;
    }
    while ( (val_k_NL != 1.) && (fabs(k2/k1-1.) > 1.0e-6) );

    return k_NL;
}

double ClassEngine::get_k_NL_from_lin_Pk_interp(const double &z)
{
    return k_NL_z_from_lin_Pk_array.interp(z);
}

double ClassEngine::get_k_NL_from_nl_Pk_class(const double &z)
{
    double k_NL, k_NL_cb;

    nonlinear_k_nl_at_z(&ba, &nl, z, &k_NL, &k_NL_cb);

    return k_NL;
}

double ClassEngine::get_n_eff_from_lin_Pk(double k, double z)
{
    // logarithmic slope of linear P(k,z)
    double n;

    nonlinear_pk_tilt_at_k_and_z(&ba, &pm, &nl, pk_linear, k, z, nl.index_pk_total, &n);

    return n;
}

double ClassEngine::get_n_eff_from_nl_Pk(double k, double z)
{
    // logarithmic slope of non-linear P(k,z)
    double n;

    nonlinear_pk_tilt_at_k_and_z(&ba, &pm, &nl, pk_nonlinear, k, z, nl.index_pk_total, &n);

    return n;
}

double ClassEngine::get_G_1_k_z_interp(double k, double z)
{
    // k in units of 1/Mpc
    std::vector<double> k_min_max = G_1_k_z.get_min_max_x_coordinates();
    double k_min = k_min_max.front();
    double k_max = k_min_max.back();

    if (k < k_min)
        return 26.0/21.0;

    else if (k > k_max)
    {
        double G_1_k_max = G_1_k_z.interp(k_max,z);

        double B_1 = -3.0/4.0;
        return B_1 + (G_1_k_max-B_1)*sqrt(k_max/k);
    }

    else
        return G_1_k_z.interp(k,z);
}

double ClassEngine::get_G_K_k_z_interp(double k, double z)
{
    // k in units of 1/Mpc
    std::vector<double> k_min_max = G_K_k_z.get_min_max_x_coordinates();
    double k_min = k_min_max.front();
    double k_max = k_min_max.back();

    if (k < k_min)
        return 8.0/7.0;

    else if (k > k_max)
    {
        double G_K_k_max = G_K_k_z.interp(k_max,z);

        double B_K = -9.0/4.0;
        return B_K + (G_K_k_max-B_K)*sqrt(k_max/k);
    }

    else
        return G_K_k_z.interp(k,z);
}

double ClassEngine::get_a_GM_k_z_interp(const double &k, const double &z)
{
  return a_GM_k_z.interp(k,z);
}

double ClassEngine::get_b_GM_k_z_interp(const double &k, const double &z)
{
  return b_GM_k_z.interp(k,z);
}

double ClassEngine::get_c_GM_k_z_interp(const double &k, const double &z)
{
  return c_GM_k_z.interp(k,z);
}


void ClassEngine::bihalofit_compute_pk_norm()
{
    m_bihalofit_norm = 1.;
    m_bihalofit_norm*= get_sigma8_z(0)/bihalofit_sigma_j(8.,0);   // P(k) amplitude normalized by sigma8
}

double ClassEngine::bihalofit_pk_lin(const double &k_h, const double &z)
{
    double h = get_h();
    return m_bihalofit_norm*m_bihalofit_norm*h*h*h*pk_lin(k_h*h,z);
}

double ClassEngine::bihalofit_window(const double &x, const int &j)
{
    if(j==0) return 3./(x*x*x)*(sin(x)-x*cos(x));  // top hat
    if(j==1) return exp(-0.5*x*x);   // gaussian
    if(j==2) return x*exp(-0.5*x*x);  // 1st derivative gaussian
    else printf("None of the possible options selected! \n");
}

double ClassEngine::bihalofit_sigma_j(const double &R, const int &j) // R [Mpc/h]
{
    double k1_h = 2.*M_PI/R; // [h/Mpc]
    double k2_h = 2.*M_PI/R; // [h/Mpc]

    double xxpp=-1.;
    double xx;

    for(;;)
    {
        k1_h=k1_h/10.;
        k2_h=k2_h*2.;

        double a=log(k1_h);
        double b=log(k2_h);

        double xxp=-1.;
        int n=2;

        for(;;)
        {
            n=n*2;
            double hh=(b-a)/(double)n;

            xx=0.;

            for(int i=1;i<n;i++)
            {
                double k_h=exp(a+hh*i); // [h/Mpc]
                xx+=k_h*k_h*k_h*bihalofit_pk_lin(k_h,0)*pow(bihalofit_window(k_h*R,j),2);
            }

            xx+=0.5*(k1_h*k1_h*k1_h*bihalofit_pk_lin(k1_h,0)*pow(bihalofit_window(k1_h*R,j),2)+k2_h*k2_h*k2_h*bihalofit_pk_lin(k2_h,0)*pow(bihalofit_window(k2_h*R,j),2));

            xx*=hh;

            if(fabs((xx-xxp)/xx) < 1.0e-6)
                break;

            xxp=xx;
        }

        if(fabs((xx-xxpp)/xx) < 1.0e-6)
            break;

        xxpp=xx;
    }

    return sqrt(xx/(2.0*M_PI*M_PI));
}

double ClassEngine::bihalofit_R_NL(const double &z)  // return r[Mpc/h] (=1/k_NL)
{
    // eqn (B1) of https://arxiv.org/pdf/1911.07886.pdf

    double D1 = get_D_plus_z(z);

    double k,k1,k2;

    k1=k2=1.; // [h/Mpc]

    for(;;)
    {
        if(D1*bihalofit_sigma_j(1./k1,1) < 1.)
            break;
        k1*=0.5;
    }

    for(;;)
    {
        if(D1*bihalofit_sigma_j(1./k2,1) > 1.)
            break;
        k2*=2.;
    }

    for(;;)
    {
        k=0.5*(k1+k2);

        if(D1*bihalofit_sigma_j(1./k,1) < 1.)
            k1=k;

        else if(D1*bihalofit_sigma_j(1./k,1) > 1.)
            k2=k;

        if(D1*bihalofit_sigma_j(1./k,1) == 1. || fabs(k2/k1-1.) < 1.0e-6)
            break;
    }

    return 1./k; // i.e. r_sigma = 1/k_NL [Mpc/h]
}

double ClassEngine::bihalofit_n_eff(const double &z, const double &R_NL)
{
    // eqn (B2) https://arxiv.org/pdf/1911.07886.pdf
    return -3. + 2.0*pow(get_D_plus_z(z)*bihalofit_sigma_j(R_NL,2),2);
}

double ClassEngine::bihalofit_R_NL_interp(const double &z)
{
    return bihalofit_R_NL_z_array.interp(z);
}

double ClassEngine::bihalofit_n_eff_interp(const double &z)
{
    return bihalofit_n_eff_z_array.interp(z);
}

//double ClassEngine::pk_primordial(double k)
//{
//    // P_R(k) = 2 pi^2/k^3 Delta_R^2(k) = 2 pi^2/k^3 A_s (k/k_pivot)^(n_s - 1)

//    // Units: k [1/Mpc] and P(k) [Mpc^3]
//    double pk;

//    primordial_spectrum_at_k(&pm, nl.index_md_scalars, linear, k, &pk);

//    return 2.*M_PI*M_PI/(k*k*k)*pk;
//}

double ClassEngine::pk_primordial(const double &k)
{
    // P_R(k) = 2 pi^2/k^3 Delta_R^2(k) = 2 pi^2/k^3 A_s (k/k_pivot)^(n_s - 1)

    // Units: k [1/Mpc] and P(k) [Mpc^3]

    return 2.*M_PI*M_PI/(k*k*k) * get_A_s() * pow(k/get_k_pivot(), get_n_s() - 1);
}

double ClassEngine::pk_gravitational_potential(const double &k)
{
    // P_phi(k) = (3/5)^2 P_R(k)

    // Units: k [1/Mpc] and P(k) [Mpc^3]

    return 0.36*pk_primordial(k);
}

double ClassEngine::M_poisson_factor(const double &k, const double &z)
{
    // Units:k [1/Mpc] and P(k) [Mpc^3]
    double M_k_z = sqrt(fabs(pk_lin(k,z)/pk_gravitational_potential(k))); // equation (18) Leicht, Baldauf et al. (2020)

    if (isnan(M_k_z))
        std::cout << "NaN encountered in the value of the Poisson factor" << std::endl;

    if (M_k_z == 0)
        std::cout << "The Poisson factor is zero and the value of k was!! = "<< k << std::endl;

    return M_k_z;

    //return 2*k*k*Tk_z0_d_tot.interp(k)*class_obj->D_plus_z(z) / (3*class_obj->get_Omega0_m()*pow(class_obj->get_Hz(0),2));
}

double ClassEngine::pk_lin(const double &k, const double &z)
{
    // Units: k [1/Mpc] and P(k) [Mpc^3]
    double pk;
    double pk_ic;

    nonlinear_pk_at_k_and_z(&ba, &pm, &nl, pk_linear, k, z, nl.index_pk_total, &pk, &pk_ic);

    return pk;
}

double ClassEngine::pk_nl(const double &k, const double &z)
{
    // Units: k [1/Mpc] and P(k) [Mpc^3]
    double pk;
    double pk_ic;

    nonlinear_pk_at_k_and_z(&ba, &pm, &nl, pk_nonlinear, k, z, nl.index_pk_total, &pk, &pk_ic);

    return pk;
}

double ClassEngine::pk(const double &k, const double &z, bool use_pk_nl)
{
    // Units: k [1/Mpc] and P(k) [Mpc^3]
    if (use_pk_nl == false)
        return pk_lin(k, z); // linear power spectrum
    else
        return pk_nl(k, z);
}
