//--------------------------------------------------------------------------
//
// Description:
// 	class Engine :
//base class for Boltzmann code
//
//
// Author List:
//	Stephane Plaszczynski (plaszczy@lal.in2p3.fr)
//
// History (add to end):
//	creation:   Tue Mar 13 15:28:50 CET 2012 
//
//------------------------------------------------------------------------

#ifndef Engine_hh
#define Engine_hh

#include<vector>
#include<ostream>

class Engine
{

public:

  enum cltype {TT=0,EE,TE,BB,PP,TP,EP}; //P stands for phi (lensing potential)

  //constructors
  Engine();

  //pure virtual:
  virtual bool updateParValues(const std::vector<double>& cosmopars)=0;

  // units = (micro-K)^2
  virtual void getCls(const std::vector<unsigned>& lVec, //input 
		      std::vector<double>& cltt, 
		      std::vector<double>& clte, 
		      std::vector<double>& clee, 
		      std::vector<double>& clbb)=0;
  

  virtual bool getLensing(const std::vector<unsigned>& lVec, //input 
			  std::vector<double>& clpp, 
			  std::vector<double>& cltp, 
			  std::vector<double>& clep)=0;


  virtual double z_drag() const=0;
  virtual double rs_drag() const =0;

  virtual double get_Dv_z(const double &z)=0;
  virtual double get_D_ang_z(const double &z)=0;
  virtual double get_sigma8_z(const double &z)=0;
  virtual double get_f_z(const double &z)=0;
  virtual double get_F_z(const double &z)=0;
  virtual double get_A_z(const double &z)=0;
  virtual double get_H_z(const double &z)=0;

  virtual double getTauReio() const=0;

  // destructor
  virtual ~Engine(){};

  //write Cl model+lensing in ostream
  virtual void writeCls(std::ostream &o);
  inline int lmax() {return _lmax;}

protected:
  int _lmax;

};

#endif

