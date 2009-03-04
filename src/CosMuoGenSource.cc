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
  CosMuoGen->nextEvent();

  // delete and re-create fEvt (memory)
  // delete fEvt;
  fEvt = new HepMC::GenEvent();
  HepMC::GenVertex* Vtx = new  HepMC::GenVertex(HepMC::FourVector(CosMuoGen->OneMuoEvt.vx(),
								  CosMuoGen->OneMuoEvt.vy(),
								  CosMuoGen->OneMuoEvt.vz(),
								  CosMuoGen->OneMuoEvt.t0()));
  HepMC::FourVector p(CosMuoGen->OneMuoEvt.px(),CosMuoGen->OneMuoEvt.py(),CosMuoGen->OneMuoEvt.pz(),CosMuoGen->OneMuoEvt.e());
  HepMC::GenParticle* Part = 
    new HepMC::GenParticle(p,CosMuoGen->OneMuoEvt.id(),1);
  Vtx->add_particle_out(Part); 
  fEvt->add_vertex(Vtx);
  fEvt->set_event_number(event());
  fEvt->set_signal_process_id(13);

  if (cmVerbosity_) fEvt->print();

  std::auto_ptr<HepMCProduct> CMProduct(new HepMCProduct());
  CMProduct->addHepMCData( fEvt );
  e.put(CMProduct);
     
  return true;
}

