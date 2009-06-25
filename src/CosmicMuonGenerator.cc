///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// modified by P. Biallass 29.03.2006 to implement new cosmic generator (CMSCGEN.cc) and new normalization of flux (CMSCGENnorm.cc)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 04.12.2008 sonne: replaced Min/MaxE by Min/MaxP to get cos_sf/ug scripts working again
// 20.04.2009 sonne: Implemented mechanism to read in multi muon events and propagate each muon
#define sim_cxx

#include "GeneratorInterface/CosmicMuonGenerator/interface/CosmicMuonGenerator.h"

void CosmicMuonGenerator::runCMG(){
  initialize();
  for (unsigned int iGen=0; iGen<NumberOfEvents; ++iGen){ nextEvent(); }
  terminate();
}

void CosmicMuonGenerator::initialize(){
  checkIn();
  if (NumberOfEvents > 0){
    RanGen.SetSeed(RanSeed); //set seed for Random Generator (seed can be controled by config-file)
    // set up "surface geometry" dimensions
    double RadiusTargetEff = RadiusOfTarget; //get this from cfg-file
    double Z_DistTargetEff = ZDistOfTarget;  //get this from cfg-file
    //double Z_CentrTargetEff = ZCentrOfTarget;  //get this from cfg-file
    if(TrackerOnly==true){
    RadiusTargetEff = RadiusTracker;
    Z_DistTargetEff = Z_DistTracker;
    }
    Target3dRadius = sqrt(RadiusTargetEff*RadiusTargetEff + Z_DistTargetEff*Z_DistTargetEff) + MinStepSize;
    if (Debug) std::cout << "  radius of sphere  around  target = " << Target3dRadius << " mm" << std::endl;
    SurfaceRadius = (SurfaceOfEarth+PlugWidth+RadiusTargetEff)*tan(MaxTheta) + Target3dRadius;  
    if (Debug) std::cout << "  starting point radius at Surface + PlugWidth = " << SurfaceRadius << " mm" << std::endl;

    OneMuoEvt.PlugVx = PlugVx;
    OneMuoEvt.PlugVz = PlugVz;
    //set energy and angle limits for CMSCGEN, give same seed as above 
    Cosmics->initialize(MinP, MaxP, MinTheta, MaxTheta, RanSeed, TIFOnly_constant, TIFOnly_linear);
   
#if ROOT_INTERACTIVE
  // book histos
  TH1D* ene = new TH1D("ene","generated energy",210,0.,1050.);
  TH1D* the = new TH1D("the","generated theta",90,0.,90.);
  TH1D* phi = new TH1D("phi","generated phi",120,0.,360.);
  TH3F* ver = new TH3F("ver","Z-X-Y coordinates",50,-25.,25.,20,-10.,10.,20,-10.,10.);
#endif
    if (EventDisplay) initEvDis();
    std::cout << std::endl;
    std::cout << "  generating " << NumberOfEvents << " events with random seed " << RanSeed << std::endl;

    if (MultiMuon) {
      MultiIn = 0;
      MultiIn = new TFile( MultiMuonFileName.c_str() );
      if (!MultiIn) std::cout << "MultiMuon=True: MultiMuonFileName='"
                              << MultiMuonFileName.c_str() << "' does not exist" << std::endl;
      else std::cout << "MultiMuonFile: " << MultiMuonFileName.c_str() << " opened!" << std::endl;
      //MultiTree = (TTree*) gDirectory->Get("sim");
      MultiTree = (TTree*) MultiIn->Get("sim");
      SimTree = new sim(MultiTree);
      SimTree->Init(MultiTree);
      SimTreeEntries = SimTree->fChain->GetEntriesFast();
      std::cout << "SimTreeEntries=" << SimTreeEntries << std::endl;

      if (MultiMuonFileFirstEvent <= 0)
	SimTree_jentry = 0;
      else 
	SimTree_jentry = MultiMuonFileFirstEvent - 1; //1=1st evt (SimTree_jentry=0)

      NcloseMultiMuonEvents = 0;
      NskippedMultiMuonEvents = 0;
    }

    if (!MultiMuon || (MultiMuon && MultiIn)) NotInitialized = false;


  }
}

void CosmicMuonGenerator::nextEvent(){

  double E = 0.; double Theta = 0.; double Phi = 0.; double RxzV = 0.; double PhiV = 0.;
  if (int(Nsel)%100 == 0) std::cout << "    generated " << int(Nsel) << " events" << std::endl;
  // generate cosmic (E,theta,phi)
  bool   notSelected = true;
  while (notSelected){
	bool   badMomentumGenerated = true;
	while (badMomentumGenerated){
	  Cosmics->generate(); //dice one event now
	  E = sqrt(Cosmics->momentum_times_charge()*Cosmics->momentum_times_charge() + MuonMass*MuonMass);
	  Theta = TMath::ACos( Cosmics->cos_theta() ) ; //angle has to be in RAD here
	  Ngen+=1.;   //count number of initial cosmic events (in surface area), vertices will be added later
	  badMomentumGenerated = false;
	  Phi = RanGen.Rndm()*(MaxPhi-MinPhi) + MinPhi;
	}
	Norm->events_n100cos(E, Theta); //test if this muon is in normalization range
	Ndiced += 1; //one more cosmic is diced
  
    // generate vertex
    double Nver = 0.;
    bool   badVertexGenerated = true;
    while (badVertexGenerated){
      RxzV = sqrt(RanGen.Rndm())*SurfaceRadius;
      PhiV = RanGen.Rndm()*TwoPi;
      // check phi range (for a sphere with Target3dRadius around the target)
      double dPhi = Pi; if (RxzV > Target3dRadius) dPhi = asin(Target3dRadius/RxzV);
      double rotPhi = PhiV + Pi; if (rotPhi > Pi) rotPhi -= TwoPi;
      double disPhi = fabs(rotPhi - Phi); if (disPhi > Pi) disPhi = TwoPi - disPhi;
      if (disPhi < dPhi) badVertexGenerated = false;
      Nver+=1.;
    }
    Ngen += (Nver-1.); //add number of generated vertices to initial cosmic events
    
    // complete event at surface
    int                             id =  13; // mu-
    if (Cosmics->momentum_times_charge() >0.) id = -13; // mu+
    double absMom = sqrt(E*E - MuonMass*MuonMass);
    double verMom = absMom*cos(Theta);
    double horMom = absMom*sin(Theta);
    double Px = horMom*sin(Phi); // [GeV/c]
    double Py = -verMom;         // [GeV/c]
    double Pz = horMom*cos(Phi); // [GeV/c]
    double Vx = RxzV*sin(PhiV);  // [mm]
    double Vy = SurfaceOfEarth + PlugWidth;  // [mm]
    double Vz = RxzV*cos(PhiV);  // [mm]
    double T0 = (RanGen.Rndm()*(MaxT0-MinT0) + MinT0)*SpeedOfLight; // [mm/c];

    Id_at = id;
    Px_at = Px; Py_at = Py; Pz_at = Pz; E_at = E; //M_at = MuonMass;
    Vx_at = Vx; Vy_at = Vy; Vz_at = Vz; T0_at = T0;

    //std::cout << "px=" << Px_at << " py=" << Py_at << " pz=" << Pz_at 
    //      << " theta=" << Theta << " phi=" << Phi << std::endl;
    //std::cout << "vx=" << Vx << " vy=" << Vy << " vz=" << Vz << std::endl;

    OneMuoEvt.create(id, Px, Py, Pz, E, MuonMass, Vx, Vy, Vz, T0); 
    // if angles are ok, propagate to target
    if (goodOrientation()) OneMuoEvt.propagate(ElossScaleFactor, RadiusOfTarget, ZDistOfTarget, ZCentrOfTarget, TrackerOnly, MTCCHalf);
    // if cosmic hits target test also if P>MinP_CMS; the default is MinP_surface=MinP_CMS, thus no bias from access shaft
    if (OneMuoEvt.hitTarget() && sqrt(OneMuoEvt.e()*OneMuoEvt.e() - MuonMass*MuonMass) > MinP_CMS){
      Nsel+=1.; //count number of generated and accepted events  
      notSelected = false;
    }
  }

  EventWeight = 1.;

  //just one outgoing particle at SurFace
  Id_sf.resize(1); 
  Px_sf.resize(1); 
  Py_sf.resize(1); 
  Pz_sf.resize(1);
  E_sf.resize(1); 
  //M_sf.resize(1);
  Vx_sf.resize(1);
  Vy_sf.resize(1);
  Vz_sf.resize(1);
  T0_sf.resize(1);

  Id_sf[0] = Id_at;
  Px_sf[0] = Px_at; Py_sf[0] = Py_at; Pz_sf[0] = Pz_at; E_sf[0] = E_at; //M_fs[0] = MuonMass;
  Vx_sf[0] = Vx_at; Vy_sf[0] = Vy_at; Vz_sf[0] = Vz_at; T0_sf[0] = T0_at;
  

  //just one particle at UnderGround  
  Id_ug.resize(1); 
  Px_ug.resize(1); 
  Py_ug.resize(1); 
  Pz_ug.resize(1);
  E_ug.resize(1); 
  //M_ug.resize(1);
  Vx_ug.resize(1);
  Vy_ug.resize(1);
  Vz_ug.resize(1);
  T0_ug.resize(1);

  Id_ug[0] = OneMuoEvt.id();
  Px_ug[0] = OneMuoEvt.px(); 
  Py_ug[0] = OneMuoEvt.py(); 
  Pz_ug[0] = OneMuoEvt.pz();
  E_ug[0] = OneMuoEvt.e(); 
  //M_ug[0] = OneMuoEvt.m();
  Vx_ug[0] = OneMuoEvt.vx(); 
  Vy_ug[0] = OneMuoEvt.vy(); 
  Vz_ug[0] = OneMuoEvt.vz();
  T0_ug[0] = OneMuoEvt.t0();

  // plot variables of selected events
#if ROOT_INTERACTIVE
  ene->Fill(OneMuoEvt.e());
  the->Fill((OneMuoEvt.theta()*Rad2Deg));
  phi->Fill((OneMuoEvt.phi()*Rad2Deg));
  ver->Fill((OneMuoEvt.vz()/1000.),(OneMuoEvt.vx()/1000.),(OneMuoEvt.vy()/1000.));
#endif
  if (Debug){
    std::cout << "new event" << std::endl;
    std::cout << "  Px,Py,Pz,E,m = " << OneMuoEvt.px() << ", " << OneMuoEvt.py() << ", "
         << OneMuoEvt.pz() << ", " << OneMuoEvt.e() << ", " << OneMuoEvt.m() << " GeV" << std::endl;
    std::cout << "  Vx,Vy,Vz,t0  = " << OneMuoEvt.vx() << ", " << OneMuoEvt.vy() << ", " 
         << OneMuoEvt.vz() << ", " << OneMuoEvt.t0() << " mm" << std::endl;
  }
  if (EventDisplay) displayEv();
  
  //std::cout << "CosmicMuonGenerator.cc: Vy_sf=" << Vy_sf[0] << "   Vy_ug[0]=" << Vy_ug[0] << std::endl;

}




bool CosmicMuonGenerator::nextMultiEvent() {

  std::cout << "\nEntered CosmicMuonGenerator::nextMultiEvent()" << std::endl;
  bool EvtRejected = true;
  bool MuInMaxDist = false;
  double MinDist; //[mm] 

  while (EvtRejected) {

    bool newEvt = true;

    //read in event from SimTree
    //ULong64_t ientry = SimTree->LoadTree(SimTree_jentry);
    Long64_t ientry = SimTree->GetEntry(SimTree_jentry);
    std::cout << "CosmicMuonGenerator::nextMultiEvent(): SimTree_jentry=" << SimTree_jentry
      //<< " ientry=" << ientry 
	      << " SimTreeEntries=" << SimTreeEntries << std::endl;
    if (ientry < 0) return false; //stop run
    if (SimTree_jentry < SimTreeEntries) {
      SimTree_jentry++;
    }
    else {
      std::cout << "CosmicMuonGenerator.cc::nextMultiEvent: No more events in file!" << std::endl;
      return false; //stop run
    }
    


    int nmuons = SimTree->shower_nParticlesWritten;
    if (nmuons < 2) {
      std::cout << "CosmicMuonGenerator.cc: Warning!  Less than two muons in event: Nmuons=" 
		<< nmuons << std::endl;
      std::cout << "trying next event from file" << std::endl;
      NskippedMultiMuonEvents++;
      //break; //EvtRejected while loop
      continue; //EvtRejected while loop: get next event from file
    }



    Px_mu.resize(nmuons); Py_mu.resize(nmuons); Pz_mu.resize(nmuons);
    P_mu.resize(nmuons);

    MinDist = 99999.e9; //[mm] 
    double MuMuDist;
    MuInMaxDist = false;
    //check if at least one muon pair closer than 30m at surface
    int NmuPmin = 0;
    for (int imu=0; imu<nmuons; ++imu) {

      Px_mu[imu] = -SimTree->particle__Px[imu]*sin(NorthCMSzDeltaPhi)
	+ SimTree->particle__Py[imu]*cos(NorthCMSzDeltaPhi);
      Pz_mu[imu] = SimTree->particle__Px[imu]*cos(NorthCMSzDeltaPhi)
	+ SimTree->particle__Py[imu]*sin(NorthCMSzDeltaPhi);
      Py_mu[imu] = -SimTree->particle__Pz[imu]; //Corsika down going particles defined in -z direction!
      P_mu[imu] = sqrt(Px_mu[imu]*Px_mu[imu] + Py_mu[imu]*Py_mu[imu] + Pz_mu[imu]*Pz_mu[imu]);
      
      if (P_mu[imu] < MinP_CMS) continue;
      else if (SimTree->particle__ParticleID[imu] != 5 &&
	       SimTree->particle__ParticleID[imu] != 6) continue;
      else NmuPmin++;
      
      for (int jmu=0; jmu<imu; ++jmu) {
	if (P_mu[jmu] < MinP_CMS) continue;
	if (SimTree->particle__ParticleID[imu] != 5 &&
	    SimTree->particle__ParticleID[imu] != 6) continue;
	MuMuDist = 0.1*sqrt( (SimTree->particle__x[imu]-SimTree->particle__x[jmu])*
			     (SimTree->particle__x[imu]-SimTree->particle__x[jmu]) 
			     +(SimTree->particle__y[imu]-SimTree->particle__y[jmu])*
			     (SimTree->particle__y[imu]-SimTree->particle__y[jmu])
			     ); //[cm] -> [mm]
	if (MuMuDist < MinDist) MinDist = MuMuDist;
	if (MuMuDist < 2.*Target3dRadius) MuInMaxDist = true;
      }
    }
    if (MuInMaxDist) {
      NcloseMultiMuonEvents++;
    }
    else {
      std::cout << "CosmicMuonGenerator.cc: Warning! No muon pair closer than " 
		<< Target3dRadius/1000. << "m   MinDist=" << MinDist/1000. << "m at surface" << std::endl;
      std::cout << "Fraction of too wide opening angle multi muon events: "
		<< 1 - double(NcloseMultiMuonEvents)/SimTree_jentry << std::endl;
      std::cout << "NcloseMultiMuonEvents=" << NcloseMultiMuonEvents << std::endl;
      std::cout << "trying next event from file" << std::endl;
      NskippedMultiMuonEvents++;
      continue; //EvtRejected while loop: get next event from file
    }
  
    if (NmuPmin < 2) continue; //EvtRejected while loop: get next event from file

    //if (Debug) 
    std::cout << "start trial do loop: MuMuDist=" << MinDist/1000. << "[m]   Nmuons=" 
	      << nmuons << "  NcloseMultiMuonEvents=" << NcloseMultiMuonEvents 
	      << "  NskippedMultiMuonEvents=" << NskippedMultiMuonEvents << std::endl;  

    
    //int primary_id = SimTree->run_ParticleID;
    Id_at = SimTree->shower_EventID;
    
    double M_at = 0.;
    //if (Id_at == 13) {
    Id_at = 2212; //convert from Corsika to HepPDT
    M_at = 938.272e-3; //[GeV] mass
    //}
    
    E_at = SimTree->shower_Energy;
    Theta_at = SimTree->shower_Theta;
    double phi_at = SimTree->shower_Phi - NorthCMSzDeltaPhi; //rotate by almost 90 degrees
    if (phi_at < -Pi) phi_at +=TwoPi; //bring into interval (-Pi,Pi]
    else if (phi_at > Pi) phi_at -= TwoPi; 
    double P_at = sqrt(E_at*E_at - M_at*M_at);
    //need to rotate about 90degrees around x->N axis => y<=>z, 
    //then rotate new x-z-plane from x->North to x->LHC centre
    Px_at = P_at*sin(Theta_at)*sin(phi_at);
    Py_at = -P_at*cos(Theta_at);
    Pz_at = P_at*sin(Theta_at)*cos(phi_at);

    //compute maximal theta of secondary muons
    double theta_mu_max = Theta_at;
    double theta_mu_min = Theta_at;

    double phi_rel_min = 0.; //phi_mu_min - phi_at
    double phi_rel_max = 0.; //phi_mu_max - phi_at

    //cout << "SimTree->shower_Energy=" << SimTree->shower_Energy << endl;

    /*
    if (newEvt) cout << "Theta_at=" << Theta_at << " phi_at=" << phi_at << " Px_at=" << Px_at
		     << " Pz_at=" << Pz_at << " Py_at=" << Py_at 
		     << " Vx_at=" << Vx_at << " Vy_at=" << Vy_at << " Vz_at=" << Vz_at 
		     << endl;
    */

    Theta_mu.resize(nmuons);
    for (int imu=0; imu<nmuons; ++imu) {
      Theta_mu[imu] = acos(-Py_mu[imu]/P_mu[imu]);
      if (Theta_mu[imu]>theta_mu_max) theta_mu_max = Theta_mu[imu];
      if (Theta_mu[imu]<theta_mu_min) theta_mu_min = Theta_mu[imu];

      double phi_mu = atan2(Px_mu[imu],Pz_mu[imu]); // in (-Pi,Pi]
      //if (phi_mu>phi_mu_max) phi_mu_max = phi_mu;
      //if (phi_mu<phi_mu_min) phi_mu_min = phi_mu;
      double phi_rel = phi_mu - phi_at;
      //cout << "phi_rel=" << phi_rel;
      if (phi_rel < -Pi) phi_rel += TwoPi; //bring into interval (-Pi,Pi] 
      else if (phi_rel > Pi) phi_rel -= TwoPi;
      //cout << " after centering: phi_rel=" << phi_rel;
      if (phi_rel < phi_rel_min) phi_rel_min = phi_rel;
      else if (phi_rel > phi_rel_max) phi_rel_max =phi_rel;

      //if (newEvt) cout << "Theta_mu[" << imu << "]=" << Theta_mu[imu] 
      //	       << " phi=" << atan2(Px_mu[imu],Pz_mu[imu]) 
      //	       << " Px=" << Px_mu[imu] << " Pz=" << Pz_mu[imu] 
      //	       << " Py=" << Py_mu[imu] << endl
;
      //cout << "SimTree->particle__Pz[" << imu << "]=" << SimTree->particle__Pz[imu] << endl;

    }

    //if (newEvt) cout << "theta_mu_min=" << theta_mu_min << " theta_mu_max=" << theta_mu_max << endl;
    //if (newEvt) cout << "phi_at=" << phi_at << " phi_rel_min=" << phi_rel_min << " phi_rel_max=" << phi_rel_max << endl;
    //newEvt = false;
    
    
    double h_sf = SurfaceOfEarth + PlugWidth; //[mm]

    double R_at = h_sf*tan(Theta_at);
    
    double JdRxzV_dR_trans = 1.; 
    double JdPhiV_dPhi_trans = 1.; 
    double JdR_trans_sqrt = 1.;
    
    //chose random vertex Phi and Rxz weighted to speed up and smoothen
    double R_mu_max = (h_sf+Target3dRadius)*tan(theta_mu_max);
    double R_max = min(SurfaceRadius, R_mu_max);
    double R_mu_min = (h_sf-Target3dRadius)*tan(theta_mu_min);
    double R_min = max(0., R_mu_min);
    
    if (R_at>SurfaceRadius) {
      cout << "CosmicMuonGenerator.cc: Warning! R_at=" << R_at 
	   << " > SurfaceRadius=" << SurfaceRadius << endl;
    }
    
    //do phase space transformation for horizontal radius R
    
    //determine first phase space limits
    
    double psR1min = R_min + 0.25*(R_max-R_min); 
    double psR1max = min(SurfaceRadius,R_max-0.25*(R_max-R_min)); //no R's beyond R_max
    double psR1 = psR1max - psR1min;

    double psR2min = R_min;
    double psR2max = R_max;
    double psR2 = psR2max - psR2min;

    double psR3min = 0.;
    double psR3max = SurfaceRadius;
    double psR3 = psR3max - psR3min;

    //double psall = psR1+psR2+psR3;
    double psRall = psR3;

    double fR1=psR1/psRall, fR2=psR2/psRall, fR3=psR3/psRall; //f1+f2+f3=130%
    double pR1=0.25, pR2=0.7, pR3=0.05;


    //do phase space transformation for azimuthal angle phi
    double psPh1 = 0.5*(phi_rel_max - phi_rel_min);
    double psPh2 = phi_rel_max - phi_rel_min;
    double psPh3 = TwoPi;
    //double psPhall = psPh1+psPh2+psPh3;
    double psPhall = psPh3;
    
    double fPh1=psPh1/psPhall, fPh2=psPh2/psPhall, fPh3=psPh3/psPhall; //(f1+f2+f3=TwoPi+f1+f2)/(TwoPi+f1+f2) 

    double pPh1=0.25, pPh2=0.7, pPh3=0.05;

    int temp = 0;
    bool btemp = false;
  
    double trials = 0;
    Vx_mu.resize(nmuons); Vy_mu.resize(nmuons); Vz_mu.resize(nmuons);
    //do { //while (Id_sf.size() < 2 && trials < max_trials)
    int NmuHitTarget = 0;
    while (NmuHitTarget < 2) { 

      NmuHitTarget = 0; //re-initialize every loop iteration
      double Nver = 0.;

      
      //chose phase space class
      double RxzV;
      double which_R_class = RanGen.Rndm();
      if (which_R_class < pR1) { //pR1% in psR1
	RxzV = psR1min + psR1 * RanGen.Rndm();
	//JdRxzV_dR_trans = fR1/pR1;
	JdRxzV_dR_trans = fR1/pR1 * SurfaceRadius/psR1;
      }
      else if (which_R_class < pR1+pR2) { //further pR2% in psR2 
	RxzV = psR2min + psR2 * RanGen.Rndm();	  
	//JdRxzV_dR_trans = fR2/pR2;
	JdRxzV_dR_trans = fR2/pR2 * SurfaceRadius/psR2;
      }
      else { //remaining pR3% in psR3=[0., R_max]
	RxzV = psR3min + psR3 * RanGen.Rndm();
	//JdRxzV_dR_trans = fR3/pR3;
	JdRxzV_dR_trans = fR3/pR3 * SurfaceRadius/psR3;
      }
      
      //double R_trans_sqrt = RanGen.Rndm()*R_max;
      //JdR_trans_sqrt = 2.*R_trans_sqrt/R_max;

      //JdR_trans_sqrt = 2.*RxzV/R_max; //flat in sqrt(r) space

      //JdR_trans_sqrt = 2.*(RxzV-R_min)/(R_max-R_min); //flat in sqrt(r) space
      
      JdR_trans_sqrt = 2.*RxzV/SurfaceRadius; //flat in sqrt(r) space

      
      // cout << "JddRxzV_trans_dR=" <<  JdRxzV_dR_trans << " f1/p1=" << f1/p1 << " f2/p2=" << f2/p2 << " f3/p3=" << f3/p3 << endl; 


      //chose phase space class
      double PhiV;
      double which_phi_class = RanGen.Rndm();
      if (which_phi_class < pPh1) { //pPh1% in psPh1 
	PhiV = phi_at + phi_rel_min + psPh1 * RanGen.Rndm();
	//JdPhiV_dPhi_trans = fPh1/pPh1;
	JdPhiV_dPhi_trans = fPh1/pPh1 * TwoPi/psPh1;
      }
      else if (which_phi_class < pPh1+pPh2) { //further pPh2% in psPh2
	PhiV = phi_at + phi_rel_min + psPh2 * RanGen.Rndm();	  
	//JdPhiV_dPhi_trans = fPh2/pPh2;
	JdPhiV_dPhi_trans = fPh2/pPh2 * TwoPi/psPh2;
      }
      else { //remaining pPh3% in psPh3=[-Pi,Pi]
	PhiV = phi_at + phi_rel_min + psPh3 * RanGen.Rndm();
	//JdPhiV_dPhi_trans = fPh3/pPh3;
	JdPhiV_dPhi_trans = fPh3/pPh3 * TwoPi/psPh3;
      }

      //shuffle PhiV into [-Pi,+Pi] interval
      if (PhiV < -Pi) PhiV+=TwoPi;
      else if (PhiV > Pi) PhiV-=TwoPi;
      
      
      Nver++;
      trials += JdR_trans_sqrt * JdRxzV_dR_trans * JdPhiV_dPhi_trans;
      if (trials > max_trials) break; //while (Id_sf.size() < 2) loop 
      Ngen += (Nver-1.); //add number of generated vertices to initial cosmic events
         
   
      Vx_at = RxzV*sin(PhiV); // [mm]
      Vy_at = h_sf + 10.e3; // [mm] (SurfaceOfEarth + PlugWidth + atmosphere height eventually)
      Vz_at = RxzV*cos(PhiV); // [mm]
      
      int NmuHitTargetSphere = 0;
      for (int imu=0; imu<nmuons; ++imu) {
	
	Vx_mu[imu] = Vx_at + (-SimTree->particle__x[imu]*sin(NorthCMSzDeltaPhi)
			      +SimTree->particle__y[imu]*cos(NorthCMSzDeltaPhi) )*10; //[mm] (Corsika cm to CMSCGEN mm)
	Vy_mu[imu] = h_sf; //[mm] fixed at surface + PlugWidth
	Vz_mu[imu] = Vz_at + ( SimTree->particle__x[imu]*cos(NorthCMSzDeltaPhi)
			       +SimTree->particle__y[imu]*sin(NorthCMSzDeltaPhi) )*10; //[mm] (Corsika cm to CMSCGEN mm)


	//only muons
	if (SimTree->particle__ParticleID[imu] != 5 &&
	    SimTree->particle__ParticleID[imu] != 6) continue;

	if (P_mu[imu] < MinP_CMS) continue; 

	//check here if at least 2 muons make it to the target sphere
	double Vxz_mu = sqrt(Vx_mu[imu]*Vx_mu[imu] + Vz_mu[imu]*Vz_mu[imu]);
	theta_mu_max = atan((Vxz_mu+Target3dRadius)/(h_sf-Target3dRadius));
	theta_mu_min = atan((Vxz_mu-Target3dRadius)/(h_sf+Target3dRadius));
	if (Theta_mu[imu] > theta_mu_min && Theta_mu[imu] < theta_mu_max) {
	  
	  // check phi range (for a sphere with Target3dRadius around the target)
	  double dPhi = Pi; if (Vxz_mu > Target3dRadius) dPhi = asin(Target3dRadius/Vxz_mu);
	  double PhiPmu = atan2(Px_mu[imu],Pz_mu[imu]); //muon phi
	  double PhiVmu = atan2(Vx_mu[imu],Vz_mu[imu]); //muon phi
	  double rotPhi = PhiVmu + Pi; if (rotPhi > Pi) rotPhi -= TwoPi;
	  double disPhi = fabs(rotPhi - PhiPmu); if (disPhi > Pi) disPhi = TwoPi - disPhi;
	  if (disPhi < dPhi) {
	    NmuHitTargetSphere++;
	  }

	  /*
	    else if (btemp == false) {
	      temp++;
	      cout << "PhiPmu[" << imu << "]=" << PhiPmu 
		   << " PxMu=" << Px_mu[imu] << " PzMu=" << Pz_mu[imu]
		   << " PhiVmu=" << PhiVmu
		   << " VxMu=" << Vx_mu[imu] << " VzMu=" << Vz_mu[imu]
		   << " rotPhiVxz=" << rotPhi 
		   << " dPhi=" << dPhi << " disPhi=" << disPhi << endl;  
	    }
	  */
	}

      } //end imu for loop
      
      if (temp > 0) btemp = true;


      if (NmuHitTargetSphere < 2) continue; //while (Id_sf.size() < 2) loop
    
      
      T0_at = (RanGen.Rndm()*(MaxT0-MinT0) + MinT0)*SpeedOfLight; // [mm/c];
      
      //std::cout << "primary: ID=" << Id_at << " E_at=" << E_at << " theta_at="
      //	<< theta_at << " phi_at=" << phi_at << " Vx_at=" << Vx_at
      //	<< " Vy_at=" << Vy_at << " Vz_at=" << Vz_at << std::endl;
      //std::cout << "vx(vy=0) = " << Vx_at - Vy_at/Py_at*Px_at 
      //	<< "   vz(vy=0) = " << Vz_at - Vy_at/Py_at*Pz_at << std::endl;
      
      //nmuons outgoing particle at SurFace
      Id_sf.clear(); 
      Px_sf.clear(); 
      Py_sf.clear(); 
      Pz_sf.clear();
      E_sf.clear(); 
      //M_sf_out.clear();
      Vx_sf.clear(); 
      Vy_sf.clear(); 
      Vz_sf.clear();
      T0_sf.clear();
      
      //nmuons particles at UnderGround  
      Id_ug.clear(); 
      Px_ug.clear(); 
      Py_ug.clear(); 
      Pz_ug.clear();
      E_ug.clear(); 
      //M_ug.clear();
      Vx_ug.clear();
      Vy_ug.clear();
      Vz_ug.clear();
      T0_ug.clear();
      
      int Id_sf_this =0;
      double Px_sf_this =0., Py_sf_this=0., Pz_sf_this=0.;
      double E_sf_this=0.;
      //double M_sf_this=0.;
      double Vx_sf_this=0., Vy_sf_this=0., Vz_sf_this=0.;
      double T0_sf_this=0.;
      

      for (int imu=0; imu<nmuons; ++imu) {

	if (P_mu[imu] < MinP_CMS) continue;
	//for the time being only muons
	if (SimTree->particle__ParticleID[imu] != 5 &&
	    SimTree->particle__ParticleID[imu] != 6) continue;
	
	Id_sf_this = SimTree->particle__ParticleID[imu];
	if (Id_sf_this == 5) Id_sf_this = -13; //mu+
	else if (Id_sf_this == 6) Id_sf_this = 13; //mu-

	else if (Id_sf_this == 1) Id_sf_this = 22; //gamma
	else if (Id_sf_this == 2) Id_sf_this = -11; //e+
	else if (Id_sf_this == 3) Id_sf_this = 11; //e-
	else if (Id_sf_this == 7) Id_sf_this = 111; //pi0
	else if (Id_sf_this == 8) Id_sf_this = 211; //pi+
	else if (Id_sf_this == 9) Id_sf_this = -211; //pi-
	else if (Id_sf_this == 10) Id_sf_this = 130; //KL0
	else if (Id_sf_this == 11) Id_sf_this = 321; //K+
	else if (Id_sf_this == 12) Id_sf_this = -321; //K-
	else if (Id_sf_this == 13) Id_sf_this = 2112; //n
	else if (Id_sf_this == 14) Id_sf_this = 2212; //p
	else if (Id_sf_this == 15) Id_sf_this = -2212; //pbar
	else if (Id_sf_this == 16) Id_sf_this = 310; //Ks0
	else if (Id_sf_this == 17) Id_sf_this = 221; //eta
	else if (Id_sf_this == 18) Id_sf_this = 3122; //Lambda
	
	else {
	  cout << "CosmicMuonGenerator.cc: Warning! Muon Id=" << Id_sf_this 
	       << " from file read in" << endl;
	  Id_sf_this = 99999; //trouble
	}
	
	Px_sf_this = Px_mu[imu];
	Py_sf_this = Py_mu[imu]; //Corsika down going particles defined in -z direction!
	Pz_sf_this = Pz_mu[imu];
	//double P_sf_this = P_mu[imu];
	E_sf_this = sqrt(P_mu[imu]*P_mu[imu] + MuonMass*MuonMass);
	double theta_sf_this = Theta_mu[imu];
	double phi_sf_this = atan2(Px_sf_this,Pz_sf_this);
	Vx_sf_this = Vx_mu[imu];
	Vy_sf_this = Vy_mu[imu]; //[mm] fixed at surface + PlugWidth
	Vz_sf_this = Vz_mu[imu];
	T0_sf_this = T0_at + SimTree->particle__Time[imu]*SpeedOfLight; // Corsika [ns] (*mm/ns) to [mm]
	OneMuoEvt.create(Id_sf_this, Px_sf_this, Py_sf_this, Pz_sf_this, E_sf_this, MuonMass, Vx_sf_this, Vy_sf_this, Vz_sf_this, T0_sf_this); 
	// if angles are ok, propagate to target
	if (goodOrientation()) {
	  OneMuoEvt.propagate(ElossScaleFactor, RadiusOfTarget, ZDistOfTarget, ZCentrOfTarget, TrackerOnly, MTCCHalf);
	}
	
	/*
	else {
	  std::cout << "CosmicMuonGenerator.cc: Orientation failed!" << std::endl; 
	}
	std::cout << "px=" << OneMuoEvt.px() << " py=" << OneMuoEvt.py() << " pz=" << OneMuoEvt.pz() 
		  << " theta=" << theta_sf_this << " phi=" << phi_sf_this << std::endl;
	std::cout << "vx=" << OneMuoEvt.vx() << " vy=" << OneMuoEvt.vy() << " Vz=" << OneMuoEvt.vz() << std::endl;
	std::cout << "vx(vy=0) = " << OneMuoEvt.vx() - OneMuoEvt.vy()/OneMuoEvt.py()*OneMuoEvt.px() 
		  << "   vz(vy=0) = " << OneMuoEvt.vz() - OneMuoEvt.vy()/OneMuoEvt.py()*OneMuoEvt.pz() << std::endl;
	*/
	//	     << std::endl;
	// if cosmic hits target test also if P>MinP_CMS; the default is MinP_surface=MinP_CMS, thus no bias from access shaft
	/*
	if (OneMuoEvt.hitTarget() ) 
	  std::cout << "Target hit!" << std::endl;
	std::cout << " MinP_CMS=" << MinP_CMS
		  << "   oneMuoEvt.p()=" << sqrt(OneMuoEvt.e()*OneMuoEvt.e() - MuonMass*MuonMass) 
		  << std::endl;
	*/
	if (OneMuoEvt.hitTarget() && sqrt(OneMuoEvt.e()*OneMuoEvt.e() - MuonMass*MuonMass) > MinP_CMS){
	  
	  Id_sf.push_back(Id_sf_this);
	  Px_sf.push_back(Px_sf_this);
	  Py_sf.push_back(Py_sf_this);
	  Pz_sf.push_back(Pz_sf_this);
	  E_sf.push_back(E_sf_this);
	  //M_sf.push_back(M_sf_this);
	  Vx_sf.push_back(Vx_sf_this);
	  Vy_sf.push_back(Vy_sf_this);
	  Vz_sf.push_back(Vz_sf_this);
	  T0_sf.push_back(T0_sf_this);
	  
	  Id_ug.push_back(OneMuoEvt.id());
	  Px_ug.push_back(OneMuoEvt.px());
	  Py_ug.push_back(OneMuoEvt.py());
	  Pz_ug.push_back(OneMuoEvt.pz());
	  E_ug.push_back(OneMuoEvt.e());
	  //M_sf.push_back(OneMuoEvt.m());
	  Vx_ug.push_back(OneMuoEvt.vx());
	  Vy_ug.push_back(OneMuoEvt.vy());
	  Vz_ug.push_back(OneMuoEvt.vz());
	  T0_ug.push_back(OneMuoEvt.t0());

	  NmuHitTarget++;
	}
      }
      
      //tmp_theta = acos(-tmp_py/tmp_p);
      //tmp_phi = atan2(tmp_px,tmp_pz);
      //std::cout << "P_at=" << P_at << " theta=" << theta_at << " phi=" << phi_at << endl;
      //std::cout << "Nmu's p=" << tmp_p << " theta=" << tmp_theta << " phi=" << tmp_phi << endl;
      
    } // while (Id_sf.size() < 2); //end of do loop
    

    if (trials > max_trials) {
      std::cout << "CosmicMuonGenerator.cc: Warning! trials reach max_trials=" << max_trials
		<< " without accepting event!" << std::endl;
      std::cout << "trying next event (" << SimTree_jentry << ") from file" << std::endl;
      NskippedMultiMuonEvents++;
      continue; //EvtRejected while loop: get next event from file
    }
    else {
      if (NmuHitTarget < 2) {
	std::cout << "CosmicMuonGenerator.cc: Warning! less than two muons hit target: N(mu=)" 
		  << NmuHitTarget << std::endl;
	std::cout << "trying next event (" << SimTree_jentry << ") from file" << std::endl;
	NskippedMultiMuonEvents++;
	continue; //EvtRejected while loop: get next event from file
      }
      else { //if (MuInMaxDist) {
	Nsel += 1;
	//EventWeight = 1./trials;
	EventWeight = JdR_trans_sqrt * JdRxzV_dR_trans * JdPhiV_dPhi_trans / trials;
	EvtRejected = false;
	std::cout << "CosmicMuonGenerator.cc: Theta_at=" << Theta_at << " phi_at=" << phi_at << " Px_at=" << Px_at << " Py_at=" << Py_at << " Pz_at=" << Pz_at << " Vx_at=" << Vx_at << " Vy_at=" << Vy_at << " Vz_at=" << Vz_at 
		  << " EventWeight=" << EventWeight << " Nmuons=" << Id_sf.size() << std::endl;
      }
    }
    
    
  } //while loop EvtRejected
  

  return true; //write event to HepMC;

}




void CosmicMuonGenerator::terminate(){
  if (NumberOfEvents > 0){
    std::cout << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << "***                                                   ***" << std::endl;
    std::cout << "***    C O S M I C   M U O N   S T A T I S T I C S    ***" << std::endl;
    std::cout << "***                                                   ***" << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << std::endl;  
    std::cout << "       number of initial cosmic events:  " << int(Ngen) << std::endl;
    std::cout << "       number of actually diced events:  " << int(Ndiced) << std::endl;
    std::cout << "       number of generated and accepted events:  " << int(Nsel) << std::endl;
    double selEff = Nsel/Ngen; // selection efficiency
    std::cout << "       event selection efficiency:  " << selEff*100. << "%" << std::endl;
    int n100cos =  Norm->events_n100cos(0., 0.); //get final amount of cosmics in defined range for normalisation of flux
    std::cout << "       events with ~100 GeV and 1 - cos(theta) < 1/2pi: " << n100cos << std::endl;
    std::cout << std::endl;
    std::cout << "       momentum range: " << MinP             << " ... " << MaxP << " GeV" << std::endl;
    std::cout << "       theta  range:   " << MinTheta*Rad2Deg << " ... " << MaxTheta*Rad2Deg << " deg" << std::endl; 
    std::cout << "       phi    range:   " << MinPhi*Rad2Deg   << " ... " << MaxPhi*Rad2Deg << " deg" << std::endl;
    std::cout << "       time   range:   " << MinT0            << " ... " << MaxT0 << " ns" << std::endl;
    std::cout << "       energy  loss:   " << ElossScaleFactor*100. << "%" << std::endl;
    std::cout << std::endl;
    double area = 1.e-6*Pi*SurfaceRadius*SurfaceRadius; // area on surface [m^2] 
    std::cout << "       area of initial cosmics on Surface + PlugWidth:   " << area << " m^2" << std::endl;
    std::cout << "       depth of CMS detector (from Surface, without PlugWidth)):   " << SurfaceOfEarth/1000 << " m" << std::endl;
    //OneMuoEvt.theta()*Rad2Deg
    if(n100cos>0 && MaxTheta<84.26*Deg2Rad){
      // rate: corrected for area and selection-Eff. and normalized to known flux, integration over solid angle (dOmega) is implicit
      // flux is normalised with respect to known flux of vertical 100GeV muons in area at suface level 
      // rate seen by detector is lower than rate at surface area, so has to be corrected for selection-Eff.
      // normalisation factor has unit [1/s/m^2] 
      // rate = N/time --> normalization factor gives 1/runtime/area 
      // normalization with respect to number of actually diced events (Ndiced)
      EventRate= (Ndiced * Norm->norm(n100cos)) * area * selEff;
      rateErr_stat = EventRate/sqrt( (double) n100cos);  // stat. rate error 
      rateErr_syst = EventRate/2.63e-3 * 0.06e-3;  // syst. rate error, from error of known flux 

      // normalisation in region 1.-cos(theta) < 1./(2.*Pi), if MaxTheta even lower correct for this
      if(MaxTheta<0.572){
	double spacean = 2.*Pi*(1.-cos(MaxTheta));
	EventRate= (Ndiced * Norm->norm(n100cos)) * area * selEff * spacean;
	rateErr_stat = EventRate/sqrt( (double) n100cos);  // rate error 
	rateErr_syst = EventRate/2.63e-3 * 0.06e-3;  // syst. rate error, from error of known flux 
      }

    }else{
      EventRate=Nsel; //no info as no muons at 100 GeV
      rateErr_stat =Nsel;
      rateErr_syst =Nsel;
      std::cout << std::endl;
      if (MinP > 100.)
	std::cout << " !!! MinP > 100 GeV. Cannot apply normalisation!" << std::endl;
      else if (MaxTheta > 84.26*Deg2Rad)
	std::cout << " !!! Note: generated cosmics exceed parameterisation. No flux calculated!" << std::endl;
      else 
	std::cout << " !!! Not enough statistics to apply normalisation (rate=1 +- 1) !!!" << std::endl;
    } 
    
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "       rate is " << EventRate << " +-" << rateErr_stat <<" (stat) " << "+-" << 
      rateErr_syst << " (syst) " <<" muons per second" << std::endl;
    if(EventRate!=0) std::cout << "       number of events corresponds to " << Nsel/EventRate << " s" << std::endl;  //runtime at CMS = Nsel/rate
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << "*********************************************************" << std::endl;
  }
}

void CosmicMuonGenerator::checkIn(){
  if (MinP < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: min.energy is out of range (0 GeV ... inf]" << std::endl << std::endl; }
  if (MaxP < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.energy is out of range (0 GeV ... inf]" << std::endl << std::endl; }
  if (MaxP <= MinP){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.energy is not greater than min.energy" << std::endl << std::endl; }
  if (MinTheta < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: min.theta is out of range [0 deg ... 90 deg)" << std::endl << std::endl; }
  if (MaxTheta < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.theta is out of range [0 deg ... 90 deg)" << std::endl << std::endl; }
  if (MaxTheta <= MinTheta){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.theta is not greater than min.theta" << std::endl << std::endl; }
  if (MinPhi < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: min.phi is out of range [0 deg ... 360 deg]" << std::endl << std::endl; }
  if (MaxPhi < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.phi is out of range [0 deg ... 360 deg]" << std::endl << std::endl; }
  if (MaxPhi <= MinPhi){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.phi is not greater than min.phi" << std::endl << std::endl; }
  if (MaxT0 <= MinT0){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.t0 is not greater than min.t0" << std::endl << std::endl; }
  if (ElossScaleFactor < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: E-loss scale factor is out of range [0 ... inf)" << std::endl << std::endl; }
}

bool CosmicMuonGenerator::goodOrientation(){
  // check angular range (for a sphere with Target3dRadius around the target)
  bool goodAngles = false;
  bool phiaccepted = false;
  bool thetaaccepted = false;
  double RxzV = sqrt(OneMuoEvt.vx()*OneMuoEvt.vx() + OneMuoEvt.vz()*OneMuoEvt.vz());
  double rVY = sqrt(RxzV*RxzV + (SurfaceOfEarth+PlugWidth)*(SurfaceOfEarth+PlugWidth));
  double Phi = OneMuoEvt.phi();
  double PhiV = atan2(OneMuoEvt.vx(),OneMuoEvt.vz()) + Pi; if (PhiV > TwoPi) PhiV -= TwoPi;
  double disPhi = fabs(PhiV - Phi); if (disPhi > Pi) disPhi = TwoPi - disPhi;
  double dPhi = Pi; if (RxzV > Target3dRadius) dPhi = asin(Target3dRadius/RxzV);
  if (disPhi < dPhi) phiaccepted = true;
  double Theta = OneMuoEvt.theta();
  double ThetaV = asin(RxzV/rVY);
  double dTheta = Pi; if (rVY > Target3dRadius) dTheta = asin(Target3dRadius/rVY);
  //std::cout << "    dPhi = " <<   dPhi << "  (" <<   Phi << " <p|V> " <<   PhiV << ")" << std::endl;
  //std::cout << "  dTheta = " << dTheta << "  (" << Theta << " <p|V> " << ThetaV << ")" << std::endl;

  /*
  if (!phiaccepted && RxzV < Target3dRadius)
  //if (RxzV < Target3dRadius)
    std::cout << "Rejected phi=" << Phi << "  PhiV="  << PhiV 
	 << "  dPhi=" << dPhi << "  disPhi=" << disPhi
	 << "  RxzV=" << RxzV << "  Target3dRadius=" << Target3dRadius 
	 << "  Theta=" << Theta << std::endl;
  */

  if (fabs(Theta-ThetaV) < dTheta) thetaaccepted = true;
  //if (thetaaccepted) std::cout << "Bingo! thetaaccepted=true!" << std::endl;
  if (phiaccepted && thetaaccepted) goodAngles = true;
  //if (phiaccepted) std::cout << "Bingo! phiaccepted=true!" << std::endl;

  return goodAngles;
}

void CosmicMuonGenerator::initEvDis(){
#if ROOT_INTERACTIVE
  float rCMS = RadiusCMS/1000.;
  float zCMS = Z_DistCMS/1000.;
  if(TrackerOnly==true){
    rCMS = RadiusTracker/1000.;
    zCMS = Z_DistTracker/1000.;
}
  TH2F* disXY = new TH2F("disXY","X-Y view",160,-rCMS,rCMS,160,-rCMS,rCMS);
  TH2F* disZY = new TH2F("disZY","Z-Y view",150,-zCMS,zCMS,160,-rCMS,rCMS);
  gStyle->SetPalette(1,0);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1.5);
  TCanvas *disC = new TCanvas("disC","Cosmic Muon Event Display",0,0,800,410);
  disC->Divide(2,1);
  disC->cd(1);
  gPad->SetTicks(1,1);
  disXY->SetMinimum(log10(MinP));
  disXY->SetMaximum(log10(MaxP));
  disXY->GetXaxis()->SetLabelSize(0.05);
  disXY->GetXaxis()->SetTitleSize(0.05);
  disXY->GetXaxis()->SetTitleOffset(1.0);
  disXY->GetXaxis()->SetTitle("X [m]");
  disXY->GetYaxis()->SetLabelSize(0.05);
  disXY->GetYaxis()->SetTitleSize(0.05);
  disXY->GetYaxis()->SetTitleOffset(0.8);
  disXY->GetYaxis()->SetTitle("Y [m]");
  disC->cd(2);
  gPad->SetGrid(1,1);
  gPad->SetTicks(1,1);
  disZY->SetMinimum(log10(MinP));
  disZY->SetMaximum(log10(MaxP));
  disZY->GetXaxis()->SetLabelSize(0.05);
  disZY->GetXaxis()->SetTitleSize(0.05);
  disZY->GetXaxis()->SetTitleOffset(1.0);
  disZY->GetXaxis()->SetTitle("Z [m]");
  disZY->GetYaxis()->SetLabelSize(0.05);
  disZY->GetYaxis()->SetTitleSize(0.05);
  disZY->GetYaxis()->SetTitleOffset(0.8);
  disZY->GetYaxis()->SetTitle("Y [m]");
#endif
}

void CosmicMuonGenerator::displayEv(){
#if ROOT_INTERACTIVE
  double RadiusDet=RadiusCMS;
  double Z_DistDet=Z_DistCMS;
  if(TrackerOnly==true){
    RadiusDet = RadiusTracker;
    Z_DistDet = Z_DistTracker;
  }
  disXY->Reset();
  disZY->Reset();
  TMarker* InteractionPoint = new TMarker(0.,0.,2);
  TArc* r8m = new TArc(0.,0.,(RadiusDet/1000.));
  TLatex* logEaxis = new TLatex(); logEaxis->SetTextSize(0.05);
  float energy = float(OneMuoEvt.e());
  float verX = float(OneMuoEvt.vx()/1000.); // [m]
  float verY = float(OneMuoEvt.vy()/1000.); // [m]
  float verZ = float(OneMuoEvt.vz()/1000.); // [m]
  float dirX = float(OneMuoEvt.px())/fabs(OneMuoEvt.py());
  float dirY = float(OneMuoEvt.py())/fabs(OneMuoEvt.py());
  float dirZ = float(OneMuoEvt.pz())/fabs(OneMuoEvt.py());
  float yStep = disXY->GetYaxis()->GetBinWidth(1);
  int   NbinY = disXY->GetYaxis()->GetNbins();
  for (int iy=0; iy<NbinY; ++iy){
    verX += dirX*yStep;
    verY += dirY*yStep;
    verZ += dirZ*yStep;
    float rXY = sqrt(verX*verX + verY*verY)*1000.; // [mm]
    float absZ = fabs(verZ)*1000.;                 // [mm]
    if (rXY < RadiusDet && absZ < Z_DistDet){
      disXY->Fill(verX,verY,log10(energy));
      disZY->Fill(verZ,verY,log10(energy));
      disC->cd(1); disXY->Draw("COLZ"); InteractionPoint->Draw("SAME"); r8m->Draw("SAME");
      logEaxis->DrawLatex((0.65*RadiusDet/1000.),(1.08*RadiusDet/1000.),"log_{10}E(#mu^{#pm})");
      disC->cd(2); disZY->Draw("COL"); InteractionPoint->Draw("SAME");
      gPad->Update();
    }
  }
#endif
}

void CosmicMuonGenerator::setNumberOfEvents(unsigned int N){ if (NotInitialized) NumberOfEvents = N; }

void CosmicMuonGenerator::setRanSeed(int N){ if (NotInitialized) RanSeed = N; }

void CosmicMuonGenerator::setMinP(double P){ if (NotInitialized) MinP = P; }

void CosmicMuonGenerator::setMinP_CMS(double P){ if (NotInitialized) MinP_CMS = P; }

void CosmicMuonGenerator::setMaxP(double P){ if (NotInitialized) MaxP = P; }

void CosmicMuonGenerator::setMinTheta(double Theta){ if (NotInitialized) MinTheta = Theta*Deg2Rad; }

void CosmicMuonGenerator::setMaxTheta(double Theta){ if (NotInitialized) MaxTheta = Theta*Deg2Rad; } 

void CosmicMuonGenerator::setMinPhi(double Phi){ if (NotInitialized) MinPhi = Phi*Deg2Rad; }

void CosmicMuonGenerator::setMaxPhi(double Phi){ if (NotInitialized) MaxPhi = Phi*Deg2Rad; }

void CosmicMuonGenerator::setMinT0(double T0){ if (NotInitialized) MinT0 = T0; }

void CosmicMuonGenerator::setMaxT0(double T0){ if (NotInitialized) MaxT0 = T0; }

void CosmicMuonGenerator::setElossScaleFactor(double ElossScaleFact){ if (NotInitialized) ElossScaleFactor = ElossScaleFact; }

void CosmicMuonGenerator::setRadiusOfTarget(double R){ if (NotInitialized) RadiusOfTarget = R; }

void CosmicMuonGenerator::setZDistOfTarget(double Z){ if (NotInitialized) ZDistOfTarget = Z; }

void CosmicMuonGenerator::setZCentrOfTarget(double Z){ if (NotInitialized) ZCentrOfTarget = Z; }

void CosmicMuonGenerator::setTrackerOnly(bool Tracker){ if (NotInitialized) TrackerOnly = Tracker; }

void CosmicMuonGenerator::setMultiMuon(bool MultiMu){ if (NotInitialized) MultiMuon = MultiMu; }

void CosmicMuonGenerator::setMultiMuonFileName(std::string MultiMuFile){ if (NotInitialized) MultiMuonFileName = MultiMuFile; }
void CosmicMuonGenerator::setMultiMuonFileFirstEvent(int MultiMuFile1stEvt){ if (NotInitialized) MultiMuonFileFirstEvent = MultiMuFile1stEvt; }

void CosmicMuonGenerator::setTIFOnly_constant(bool TIF){ if (NotInitialized) TIFOnly_constant = TIF; }

void CosmicMuonGenerator::setTIFOnly_linear(bool TIF){ if (NotInitialized) TIFOnly_linear = TIF; }

void CosmicMuonGenerator::setMTCCHalf(bool MTCC){ if (NotInitialized) MTCCHalf = MTCC; }

void CosmicMuonGenerator::setPlugVx(double PlugVtx){ if (NotInitialized) PlugVx = PlugVtx; }
void CosmicMuonGenerator::setPlugVz(double PlugVtz){ if (NotInitialized) PlugVz = PlugVtz; }


double CosmicMuonGenerator::getRate(){ return EventRate; }
