#ifndef CosmicMuonParameters_h
#define CosmicMuonParameters_h
//
// Parameters for CosmicMuonGenerator by droll (05/DEC/2005)
//
//
// added plug and clay(moraine) specific constants, sonne (15/Jan/2009)
//
#include "TMath.h"

// flags
const bool Debug = false; // debugging printout
const bool EventDisplay = true; // display single events (if ROOT_INTERACTIVE is defined as true)

// algorithmic constants
const double MinStepSize = 10.; // minimal propagation step size [mm] must be small compared to target size
// mathematical constants
const double Pi = acos(-1.); // [rad]
const double TwoPi = 2.0*Pi; // [rad]
const double Deg2Rad = Pi/180.; // [deg] -> [rad]
const double Rad2Deg = 180./Pi; // [rad] -> [deg]
// physical constants
const double SpeedOfLight = 299.792458; // [mm/ns]
const double MuonMass = 0.105658357; // [GeV/c^2]
//const double ChargeFrac = 0.545454545; // n(mu+)/n(mu-) ~ 1.2 defined in CMSCGEN
// geometry
const double SurfaceOfEarth = 88874.; // Y-distance to surface of earth [mm]
const double Z_PX56 = -14000.; // [mm] Z position of PX56 centre [mm]
// densities of materials
const double RhoAir  = 0.00; // [g cm^-3]
const double RhoWall = 2.5; // [g cm^-3]
const double RhoRock = 2.20; // [g cm^-3]
const double RhoClay = 1.50; // [g cm^-3]
const double RhoPlug = 6.3; // [g cm^-3] 2-3 times concrete density
// width of clay layer between surface and rock
const double ClayWidth = 50000.; // [mm]
//plug constants
const double PlugWidth = 2250.; // [mm]
const double PlugXlength = 20600.; // [mm]
const double PlugZlength = 16000.; // [mm]
const double PlugNoseXlength = 6400.; // [mm]
const double PlugNoseZlength = 1800.; // [mm]
const double PlugOnShaftVx = 0.; // [mm]
const double PlugOnShaftVz = Z_PX56; // [mm]

// cylinder around CMS (with R, +-Z)
// WARNING: These values will be set to tracker-only setup if "TrackerOnly=true" in .cfg-file. 
// This means R=1200 and Z=2800, no material or B-field outside is considered
const double RadiusCMS =  8000.; // [mm]
const double Z_DistCMS = 15000.; // [mm]
const double RadiusTracker =  1200.; // [mm]
const double Z_DistTracker = 2800.; // [mm]
// cylinder actually used in the code
//const double RadiusTarget = RadiusCMS; // [mm]  //now controlled by cfg-file!!!
//const double Z_DistTarget = Z_DistCMS; // [mm]  //now controlled by cfg-file!!!

const double MuProdAlt = 7.e6; // 7km in [mm]

//define different materials
enum {Unknown=0, Plug, Wall, Air, Clay, Rock};

//Multi Muon relevant parameters
const double NorthCMSzDeltaPhi = 3./8.*Pi; //rad (Pi/2 if CMS -x = North)
const double maxMultiMuDist = 30000.; //30000.; //30m [mm]
const int max_trials = 1000000;


#endif
