#include "GeneratorInterface/CosmicMuonGenerator/interface/CosMuoGenSource.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"


edm::CosMuoGenSource::CosMuoGenSource( const ParameterSet & pset, InputSourceDescription const& desc ) :
  GeneratedInputSource(pset, desc ) ,  
  //RanS(pset.getParameter<int>("RanSeed", 123456)), //get seed now from Framework
  MinP(pset.getParameter<double>("MinP")),
  MinP_CMS(pset.getParameter<double>("MinP_CMS")),
  MaxP(pset.getParameter<double>("MaxP")),
  MinT(pset.getParameter<double>("MinTheta")),
  MaxT(pset.getParameter<double>("MaxTheta")),
  MinPh(pset.getParameter<double>("MinPhi")),
  MaxPh(pset.getParameter<double>("MaxPhi")),
  MinS(pset.getParameter<double>("MinT0")),
  MaxS(pset.getParameter<double>("MaxT0")),
  ELSF(pset.getParameter<double>("ElossScaleFactor")),
  RTarget(pset.getParameter<double>("RadiusOfTarget")),
  ZTarget(pset.getParameter<double>("ZDistOfTarget")),
  ZCTarget(pset.getParameter<double>("ZCentrOfTarget")),
  TrackerOnly(pset.getParameter<bool>("TrackerOnly")),
  MultiMuon(pset.getParameter<bool>("MultiMuon")),
  MultiMuonFileName(pset.getParameter<std::string>("MultiMuonFileName")),
  MultiMuonFileFirstEvent(pset.getParameter<int>("MultiMuonFileFirstEvent")),
  TIFOnly_constant(pset.getParameter<bool>("TIFOnly_constant")),
  TIFOnly_linear(pset.getParameter<bool>("TIFOnly_linear")),
  MTCCHalf(pset.getParameter<bool>("MTCCHalf")),
  PlugVtx(pset.getParameter<double>("PlugVx")),
  PlugVtz(pset.getParameter<double>("PlugVz")),
  cmVerbosity_(pset.getParameter<bool>("Verbosity"))
  {

    //if not specified (i.e. negative) then use MinP also for MinP_CMS
    if(MinP_CMS < 0) MinP_CMS = MinP;

    //get seed now from Framework
    edm::Service<edm::RandomNumberGenerator> rng;
    RanS = rng->mySeed();
    // set up the generator
    CosMuoGen = new CosmicMuonGenerator();
    CosMuoGen->setNumberOfEvents(numberEventsInRun());
    CosMuoGen->setRanSeed(RanS);
    CosMuoGen->setMinP(MinP);
    CosMuoGen->setMinP_CMS(MinP_CMS);
    CosMuoGen->setMaxP(MaxP);
    CosMuoGen->setMinTheta(MinT);
    CosMuoGen->setMaxTheta(MaxT);
    CosMuoGen->setMinPhi(MinPh);
    CosMuoGen->setMaxPhi(MaxPh);
    CosMuoGen->setMinT0(MinS);
    CosMuoGen->setMaxT0(MaxS);
    CosMuoGen->setElossScaleFactor(ELSF);
    CosMuoGen->setRadiusOfTarget(RTarget);
    CosMuoGen->setZDistOfTarget(ZTarget);
    CosMuoGen->setZCentrOfTarget(ZCTarget);
    CosMuoGen->setTrackerOnly(TrackerOnly);
    CosMuoGen->setMultiMuon(MultiMuon);
    CosMuoGen->setMultiMuonFileName(MultiMuonFileName);
    CosMuoGen->setMultiMuonFileFirstEvent(MultiMuonFileFirstEvent);
    CosMuoGen->setTIFOnly_constant(TIFOnly_constant);
    CosMuoGen->setTIFOnly_linear(TIFOnly_linear);
    CosMuoGen->setMTCCHalf(MTCCHalf);
    CosMuoGen->setPlugVx(PlugVtx);
    CosMuoGen->setPlugVz(PlugVtz);    
    CosMuoGen->initialize();
    produces<HepMCProduct>();
    //  fEvt = new HepMC::GenEvent();
  }

edm::CosMuoGenSource::~CosMuoGenSource(){
  CosMuoGen->terminate();
  delete CosMuoGen;
  //  delete fEvt;
  clear();
}

void edm::CosMuoGenSource::clear(){}

bool edm::CosMuoGenSource::produce(Event &e)
{
  // generate event
  if (!MultiMuon) {
    CosMuoGen->nextEvent();
  }
  else {
    bool success = CosMuoGen->nextMultiEvent();
    if (!success) return false;
  }
 
  std::cout << "CosMuoGenSource.cc: CosMuoGen->EventWeight=" << CosMuoGen->EventWeight 
	    << "  CosMuoGen: Nmuons=" << CosMuoGen->Id_sf.size() << std::endl; 
  std::cout << "CosMuoGen->Id_at=" << CosMuoGen->Id_at
	    << "  CosMuoGen->Vx_at=" << CosMuoGen->Vx_at 
	    << "  CosMuoGen->Vy_at=" << CosMuoGen->Vy_at
	    << "  CosMuoGen->Vz_at=" << CosMuoGen->Vz_at << std::endl;
  std::cout << "  Px=" << CosMuoGen->Px_at
	    << "  Py=" << CosMuoGen->Py_at
	    << "  Pz=" << CosMuoGen->Pz_at << std::endl;
  for (unsigned int i=0; i<CosMuoGen->Id_sf.size(); ++i) {
    std::cout << "Id_sf[" << i << "]=" << CosMuoGen->Id_sf[i]
	      << "  Vx_sf[" << i << "]=" << CosMuoGen->Vx_sf[i]
	      << "  Vy_sf=" << CosMuoGen->Vy_sf[i]
	      << "  Vz_sf=" << CosMuoGen->Vz_sf[i]
	      << "  Px_sf=" << CosMuoGen->Px_sf[i]
	      << "  Py_sf=" << CosMuoGen->Py_sf[i]
	      << "  Pz_sf=" << CosMuoGen->Pz_sf[i] << std::endl;
    std::cout << "Id_ug[" << i << "]=" << CosMuoGen->Id_ug[i] 
	      << "  Vx_ug[" << i << "]=" << CosMuoGen->Vx_ug[i] 
	      << "  Vy_ug=" << CosMuoGen->Vy_ug[i]
	      << "  Vz_ug=" << CosMuoGen->Vz_ug[i]
	      << "  Px_ug=" << CosMuoGen->Px_ug[i]
	      << "  Py_ug=" << CosMuoGen->Py_ug[i]
	      << "  Pz_ug=" << CosMuoGen->Pz_ug[i] << std::endl;
  }


  fEvt = new HepMC::GenEvent();

  HepMC::GenVertex* Vtx_at = new  HepMC::GenVertex(HepMC::FourVector(CosMuoGen->Vx_at/10., //[mm->cm]
  							     CosMuoGen->Vy_at/10., //[mm->cm]
  							     CosMuoGen->Vz_at/10., //[mm->cm]
  							     CosMuoGen->T0_at/10.)); //[mm->cm]
  HepMC::FourVector p_at(CosMuoGen->Px_at,CosMuoGen->Py_at,CosMuoGen->Pz_at,CosMuoGen->E_at);
  HepMC::GenParticle* Part_at =
    new HepMC::GenParticle(p_at,CosMuoGen->Id_at, 3);//Comment mother particle in
  Vtx_at->add_particle_in(Part_at);


  //loop here in case of multi muon events (else just one iteration)
  for (unsigned int i=0; i<CosMuoGen->Id_sf.size(); ++i) {

    HepMC::FourVector p_sf(CosMuoGen->Px_sf[i],CosMuoGen->Py_sf[i],CosMuoGen->Pz_sf[i],CosMuoGen->E_sf[i]);
    HepMC::GenParticle* Part_sf_in =
      new HepMC::GenParticle(p_sf,CosMuoGen->Id_sf[i], 3); //Comment daughter particle
    Vtx_at->add_particle_out(Part_sf_in);
    
    HepMC::GenVertex* Vtx_sf = new HepMC::GenVertex(HepMC::FourVector(CosMuoGen->Vx_sf[i]/10.,                             CosMuoGen->Vy_sf[i]/10., CosMuoGen->Vz_sf[i]/10., CosMuoGen->T0_sf[i]/10.));
    HepMC::GenParticle* Part_sf_out =
      new HepMC::GenParticle(p_sf,CosMuoGen->Id_sf[i], 3); //Comment daughter particle
    
    Vtx_sf->add_particle_in(Part_sf_in);
    Vtx_sf->add_particle_out(Part_sf_out);
    
    fEvt->add_vertex(Vtx_sf); //one per muon

    HepMC::GenVertex* Vtx_ug = new HepMC::GenVertex(HepMC::FourVector(CosMuoGen->Vx_ug[i]/10.,                             CosMuoGen->Vy_ug[i]/10., CosMuoGen->Vz_ug[i]/10., CosMuoGen->T0_ug[i]/10.));
    
    HepMC::FourVector p_ug(CosMuoGen->Px_ug[i],CosMuoGen->Py_ug[i],CosMuoGen->Pz_ug[i],CosMuoGen->E_ug[i]);
    HepMC::GenParticle* Part_ug =
      new HepMC::GenParticle(p_ug,CosMuoGen->Id_ug[i], 1);//Final state daughter particle

    Vtx_ug->add_particle_in(Part_sf_out);
    Vtx_ug->add_particle_out(Part_ug);

    fEvt->add_vertex(Vtx_ug); //one per muon

  }

  fEvt->add_vertex(Vtx_at);
  fEvt->set_signal_process_vertex(Vtx_at);

  fEvt->set_event_number(event());
  fEvt->set_signal_process_id(13);

  fEvt->weights().push_back( CosMuoGen->EventWeight ); // just one event weight 
  fEvt->weights().push_back( CosMuoGen->Trials ); // int Trials number (unweighted) 


  if (cmVerbosity_) fEvt->print();

  std::auto_ptr<HepMCProduct> CMProduct(new HepMCProduct());
  CMProduct->addHepMCData( fEvt );
  e.put(CMProduct);

  return true;
}
