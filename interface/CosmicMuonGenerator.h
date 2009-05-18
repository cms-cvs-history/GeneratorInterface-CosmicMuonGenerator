#ifndef CosmicMuonGenerator_h
#define CosmicMuonGenerator_h
//
// CosmicMuonGenerator by droll (04/DEC/2005)
// modified by P. Biallass 29.03.2006 to implement new cosmic generator (CMSCGEN.cc) 
//

// include files
#include "GeneratorInterface/CosmicMuonGenerator/interface/CMSCGENnorm.h"
#include "GeneratorInterface/CosmicMuonGenerator/interface/CMSCGEN.h"
#include "GeneratorInterface/CosmicMuonGenerator/interface/CosmicMuonParameters.h"
#include "GeneratorInterface/CosmicMuonGenerator/interface/SingleParticleEvent.h"
#include <iostream>
#include <string>
#include <vector>
#include "TRandom2.h"
#include "TFile.h"
#include "TTree.h"

#include "GeneratorInterface/CosmicMuonGenerator/interface/sim.h"


using namespace std;


// class definitions
class CosmicMuonGenerator{
public:
  // constructor
  CosmicMuonGenerator(){
    //initialize class which normalizes flux (added by P.Biallass 29.3.2006)
    Norm = new CMSCGENnorm();
    //initialize class which produces the cosmic muons  (modified by P.Biallass 29.3.2006)
    Cosmics = new CMSCGEN();
    // set default control parameters
    NumberOfEvents = 100;
    RanSeed = 135799468;
    MinP =     3.;
    MinP_CMS =     MinP;
    MaxP =   3000.;
    MinTheta =  0.*Deg2Rad;
    //MaxTheta = 84.26*Deg2Rad;
    MaxTheta = 89.0*Deg2Rad;
    MinPhi =    0.*Deg2Rad;
    MaxPhi =  360.*Deg2Rad;
    MinT0  = -12.5;
    MaxT0  =  12.5;
    ElossScaleFactor = 1.0;
    RadiusOfTarget = 8000.;
    ZDistOfTarget = 15000.;
    ZCentrOfTarget = 0.;
    TrackerOnly = false;
    MultiMuon = false;
    MultiMuonFileName = "dummy.root";
    TIFOnly_constant = false;
    TIFOnly_linear = false;
    MTCCHalf = false;
    EventRate = 0.;
    rateErr_stat = 0.;
    rateErr_syst = 0.;

    SumIntegrals = 0.;
    Ngen = 0.;
    Nsel = 0.;
    Ndiced = 0.;
    NotInitialized = true;
    Target3dRadius = 0.;
    SurfaceRadius = 0.;
    //set plug as default onto PX56 shaft
    PlugVx = PlugOnShaftVx;
    PlugVz = PlugOnShaftVz;


    std::cout << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << "***                                                   ***" << std::endl;
    std::cout << "***  C O S M I C  M U O N  G E N E R A T O R  (vC++)  ***" << std::endl;
    std::cout << "***                                                   ***" << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << std::endl;
  }
  // destructor
  ~CosmicMuonGenerator(){
    delete Norm; 
    delete Cosmics;
  }
  // event with one particle
  SingleParticleEvent OneMuoEvt;
  double EventWeight; //for multi muon events
  TFile* MultiIn; //file to be read in
  TTree* MultiTree; //tree of file with multi muon events
  sim* SimTree; //class to acces tree branches
  ULong64_t SimTreeEntries;
  ULong64_t SimTree_jentry;
  int NcloseMultiMuonEvents;

  int Id_at;
  double Px_at; double Py_at; double Pz_at; 
  double E_at; 
  //double M_at;
  double Vx_at; double Vy_at; double Vz_at; 
  double T0_at;
  
  vector<int> Id_sf;
  vector<double> Px_sf; vector<double> Py_sf; vector<double> Pz_sf; 
  vector<double> E_sf; 
  //vector<double> M_sf;
  vector<double> Vx_sf; vector<double> Vy_sf; vector<double> Vz_sf; 
  vector<double> T0_sf;
  
  vector<int> Id_ug;
  vector<double> Px_ug; vector<double> Py_ug; vector<double> Pz_ug;
  vector<double> E_ug; 
  //vector<double> M_ug;
  vector<double> Vx_ug; vector<double> Vy_ug; vector<double> Vz_ug;
  vector<double> T0_ug;


 

private:
  //initialize class which normalizes flux (added by P.Biallass 29.3.2006)
  CMSCGENnorm*  Norm ;
  //initialize class which produces the cosmic muons  (modified by P.Biallass 29.3.2006)
  CMSCGEN* Cosmics ; 
  // default control parameters
  unsigned int NumberOfEvents; // number of events to be generated
  int    RanSeed; // seed of random number generator
  double MinP;     // min. E     [GeV]
  double MinP_CMS; // min. E at CMS surface    [GeV]; default is MinE_CMS=MinE, thus no bias from access-shaft
  double MaxP;     // max. E     [GeV]
  double MinTheta; // min. theta [rad]
  double MaxTheta; // max. theta [rad]
  double MinPhi;   // min. phi   [rad]
  double MaxPhi;   // max. phi   [rad]
  double MinT0;    // min. t0   [ns]
  double MaxT0;    // max. t0   [ns]
  double ElossScaleFactor; // scale factor for energy loss
  double RadiusOfTarget; // Radius of target-cylinder which cosmics HAVE to hit [mm], default is CMS-dimensions
  double ZDistOfTarget; // z-length of target-cylinder which cosmics HAVE to hit [mm], default is CMS-dimensions
  double ZCentrOfTarget; // z-position of centre of target-cylinder which cosmics HAVE to hit [mm], default is Nominal Interaction Point (=0)
  bool   TrackerOnly; //if set to "true" detector with tracker-only setup is used, so no material or B-field outside is considerd
  bool   MultiMuon; //read in multi-muon events from file instead of generating single muon events
  std::string MultiMuonFileName; //file containing multi muon events, to be read in
  bool   TIFOnly_constant; //if set to "true" cosmics can also be generated below 2GeV with unphysical constant energy dependence
  bool   TIFOnly_linear; //if set to "true" cosmics can also be generated below 2GeV with unphysical linear energy dependence
  bool   MTCCHalf; //if set to "true" muons are sure to hit half of CMS important for MTCC, 
                   //still material and B-field of whole CMS is considered
  double EventRate; // number of muons per second [Hz]
  double rateErr_stat; // stat. error of number of muons per second [Hz]
  double rateErr_syst; // syst. error of number of muons per second [Hz] from error of known flux
  // other stuff needed
  double SumIntegrals; // sum of phase space integrals
  double Ngen; // number of generated events
  double Nsel; // number of selected events
  double Ndiced; // number of diced events
  double Target3dRadius; // radius of sphere around target (cylinder)
  double SurfaceRadius; // radius for area on surface that has to be considered (for event generation)
  double PlugVx; //Plug x position
  double PlugVz; //Plug z position

  // random number generator (periodicity > 10**14)
  TRandom2 RanGen; 
  // check user input
  bool NotInitialized;
  void checkIn();
  // check, if muon is pointing into target
  bool goodOrientation();
  // event display: initialize + display
  void initEvDis();
  void displayEv();

public:
  // set parameters
  void setNumberOfEvents(unsigned int N);
  void setRanSeed(int N);
  void setMinP(double P);
  void setMinP_CMS(double P);
  void setMaxP(double P);
  void setMinTheta(double Theta);
  void setMaxTheta(double Theta);
  void setMinPhi(double Phi);
  void setMaxPhi(double Phi);
  void setMinT0(double T0);
  void setMaxT0(double T0);
  void setElossScaleFactor(double ElossScaleFact);
  void setRadiusOfTarget(double R);
  void setZDistOfTarget(double Z);
  void setZCentrOfTarget(double Z);
  void setTrackerOnly(bool Tracker);
  void setMultiMuon(bool MultiMu);
  void setMultiMuonFileName(std::string MultiMuonFileName);
  void setTIFOnly_constant(bool TIF);
  void setTIFOnly_linear(bool TIF);
  void setMTCCHalf(bool MTCC);
  void setPlugVx(double PlugVtx);
  void setPlugVz(double PlugVtz);

  // initialize the generator
  void initialize();
   // prints rate + statistics
  void terminate();
  // initialize, generate and terminate the Cosmic Muon Generator
  void runCMG();
  // returns event rate
  double getRate();
  // generate next event/muon
  void nextEvent();
  // generate next multi muon event
  bool nextMultiEvent();
};
#endif
