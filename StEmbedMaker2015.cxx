/*
// class:    StEmbedMaker2015
// author:   Hannah Harrison 
// RCF user: harrison
*/

#include "math.h"
#include "StEmbedMaker2015.h"
#include "TDataSetIter.h"
#include "TChain.h"
#include "TTree.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include "TVector.h"
#include "TLorentzVector.h"

// jet
#include "StSpinPool/StJetEvent/StJetEvent.h"
#include "StSpinPool/StJetEvent/StJetElement.h"
#include "StSpinPool/StJetEvent/StJetCandidate.h"
#include "StSpinPool/StJetEvent/StJetParticle.h"
#include "StSpinPool/StJetEvent/StJetTrack.h"
#include "StSpinPool/StJetEvent/StJetTower.h"
#include "StSpinPool/StJetEvent/StJetVertex.h"
#include "StSpinPool/StJetEvent/StJetEventTypes.h"

#include "StSpinPool/StJetSkimEvent/StJetSkimEvent.h"
#include "StSpinPool/StJetSkimEvent/StPythiaEvent.h"


// underlying event
#include "StSpinPool/StUeEvent/StUeJet.h"
#include "StSpinPool/StUeEvent/StUeVertex.h"
#include "StSpinPool/StUeEvent/StUeOffAxisCones.h"
#include "StSpinPool/StUeEvent/StUeOffAxisConesEvent.h"
#include "StSpinPool/StUeEvent/StUeOffAxisConesJet.h"

// MC Event & Pythia
#include "StMcEvent/StMcEvent.hh"
#include "StMcEvent/StMcVertex.hh"
#include "StMcEvent/StMcTrack.hh"


// startless TOF
#include "StBTofHeader.h"
#include "StBTofCollection.h"
#include "StBTofHit.h"
#include "StBTofPidTraits.h"
#include "StEventTypes.h"


ClassImp(StEmbedMaker2015)

//_________________________________________
// StEmbedMaker2015 constructor  

StEmbedMaker2015::StEmbedMaker2015(const char *name, const char *outputfile, TChain *jetChain, TChain *skimChain, TChain *UeChain):StMaker(name), mJetChain(jetChain), mSkimChain(skimChain), mUeChain(UeChain){
//save output ROOT file with name outputfile
  FFRootFile = TString(outputfile);
}

//______________________________________ /// StEmbedMaker2015 destructor                                                                                       
StEmbedMaker2015::~StEmbedMaker2015(){
}   

//____________________________________ /// Init- StChain initializes makers                                                     
Int_t StEmbedMaker2015::Init(){
  
  FFfile= new TFile(FFRootFile, "RECREATE");  
  //create tree with (name, title)
  FFTree=new TTree("Embed","Embed");

  //***** Histograms *****
  float pi= TMath::Pi();

  NumEvts=new TH1F("NumEvts","NumEvts",100,0,100);
  DeltRJet=new TH1F("DeltRJetMatch","DeltRJetMatch",10000,0,2*pi);
  DeltRJet_best=new TH1F("DeltRJet_bestmatch","DeltRJet_bestmatch",10000,0,2*pi);

  DetJetPt_all=new TH1F("DetJetPt_all","DetJetPt_all",100,0,100);
  PartJetPt_all=new TH1F("PartJetPt_all","PartJetPt_all",100,0,100);
  DetJetPt_nomatch=new TH1F("DetJetPt_nomatch","DetJetPt_nomatch",100,0,100);
  PartJetPt_nomatch=new TH1F("PartJetPt_nomatch","PartJetPt_nomatch",100,0,100);

  PartJetEta_all=new TH1F("PartJetEta_all","PartJetEta_all",10000,-1*pi,pi);
  PartJetEta_cut=new TH1F("PartJetEta_cut","PartJetEta_cut",10000,-1*pi,pi);


  //diagnostic hists to see why particle jT is weird?
  DeltRParticle=new TH1F("DeltRParticle","DeltRParticle",1000,-2*pi,2*pi);
  DeltEtaParticle=new TH1F("DeltEtaParticle","DeltEtaParticle",10000,-10,10);
  DeltPhiParticle=new TH1F("DeltPhiParticle","DeltPhiParticle",1000,-2*pi,2*pi);
  DeltEtaPhi_Ratio=new TH1F("DeltEta/Phi","DeltEta/Phi",10000,-100,100);


  //initialize jet events, skim events
  mJetEvent =0;
  mJetEventParticle=0;
  mSkimEvent =0;
  mUeEvent =0;
  nEventsAnalyzed=0;

  
  //********* add branches to TChain ************
  //names must match input trees!

  //Jet File
  mJetChain->SetBranchAddress("AntiKtR060NHits12", &mJetEvent);
  mJetChain->SetBranchAddress("AntiKtR060Particle", &mJetEventParticle);

  //Skim File
  mSkimChain->SetBranchAddress("skimEventBranch", &mSkimEvent);

  //Skim/Particle File
  mUeChain->SetBranchAddress("AntiKtR060NHits12OffAxisConesR060", &mUeEvent);
  mUeChain->SetBranchAddress("AntiKtR060ParticleOffAxisConesR060", &mUeEvent);


  //********* add branches to TTree ************
  // FFTree->Branch("RunNumber",&runNum, "runNum/I");
  FFTree->Branch("PartonicPt", "std::vector<Float_t>",&PartonicpT);
  
  FFTree->Branch("PionpT", "std::vector<Float_t>",&PionpT);
  // FFTree->Branch("PionMomentum","std::vector<Float_t>",&PionMomentum);
  //FFTree->Branch("PionEta", "std::vector<Float_t>",&PionEta);
  //FFTree->Branch("PionPhi", "std::vector<Float_t>",&PionPhi);
  //FFTree->Branch("ChargeSign","std::vector<Float_t>",&ChargeSign);

  FFTree->Branch("JetPt", "std::vector<Float_t>",&JetpT);
  //FFTree->Branch("JetEta", "std::vector<Float_t>",&JetEta);
  //FFTree->Branch("JetPhi", "std::vector<Float_t>",&JetPhi);
  //FFTree->Branch("JetRt", "std::vector<Float_t>", &JetRt);

  FFTree->Branch("DeltaR", "std::vector<Float_t>",&DeltaR);
  FFTree->Branch("DetZh", "std::vector<Float_t>",&DetZh);
  FFTree->Branch("DetjT", "std::vector<Float_t>",&DetjT);
  //FFTree->Branch("JetMomentum","std::vector<Float_t>",&JetMomentum);
  //FFTree->Branch("JetDetEta", "std::vector<Float_t>",&JetDetEta);
  //FFTree->Branch("Eta", "std::vector<Float_t>",&Eta);
  
  //FFTree->Branch("VertexZ","std::vector<Float_t>",&VertexZ);
  //FFTree->Branch("SumTrackpT","std::vector<Float_t>",&SumTrackpT);
  //FFTree->Branch("NSigTOFPi","std::vector<Float_t>",&NSigTOFPi);
  //FFTree->Branch("NSigTPCPi", "std::vector<Float_t>",&NSigTPCPi);
  //FFTree->Branch("NHitsFit", "std::vector<Float_t>",&NHitsFit);

  // UE INFO -- DETECTOR
  FFTree->Branch("UeDensity","std::vector<Float_t>",&UeRho);
  FFTree->Branch("UeArea","std::vector<Float_t>",&UeA);

  // PARTICLE-LEVEL INFO
  //FFTree->Branch("DeltaRjet", "std::vector<Float_t>",&DeltaRjet);
  FFTree->Branch("ParticleDeltaR", "std::vector<Float_t>",&ParticleDeltaR);
  FFTree->Branch("ParticleJetpT", "std::vector<Float_t>",&ParticleJetpT);
  FFTree->Branch("ParticleZh", "std::vector<Float_t>",&ParticleZh);
  FFTree->Branch("ParticlejT", "std::vector<Float_t>",&ParticlejT);


  // UE INFO -- PARTICLE
  FFTree->Branch("ParticleUeDensity","std::vector<Float_t>",&ParticleUeRho);
  FFTree->Branch("ParticleUeArea","std::vector<Float_t>",&ParticleUeA);

  //initialize tree values
  InitTreeVals();
  
  return StMaker::Init();
} //END INIT


//___________________________________/// Make - this method is called in loop for each event
Int_t StEmbedMaker2015::Make(){


  nEventsAnalyzed++;
  NumEvts->Fill(1);

  //initialize non-struct variables
  VerZ = jetPt = jetPhi = jetEta = jetRt=jetE =jetEt = 0;
  NTracks = NTow = TrackPt = TrackEta = TrackPhi = 0.;
  partonicPt=0;
  TowPt = TowEta = TowPhi = TowEt= nSigPi= nSigTOFPi= 0.;
  NVert = NJets = deltaR = HitsFit= 0;
  Zh = jT= jetDetEta= SumTrackPt=UeTrackPt= 0;

  deltaRjet=particlejetPt=0;
  UeDensity=UeJetArea=0;
  UeDensity_particle= UeJetArea_particle =0;
 
  runNum=0;
  float pi= TMath::Pi();

  assert(mJetEvent && mSkimEvent && mUeEvent &&mJetEventParticle);
  assert(mJetEvent->runId()==mSkimEvent->runId());

  runNum= mJetEvent->runId();
 
   //OLD DET LEVEL VERT
  //if (ivertex>1) continue;
  //mJetVertex =mJetEvent->vertex(ivertex);
  //    if(mJetVertex->ranking()> pow(10,6)){

  //********************************************
  // DETECTOR LEVEL
  mJetVertex =mJetEvent->vertex(0); //best primary vert is at pos 0
  NJets =mJetVertex->numberOfJets();  

  // PARTICLE LEVEL
  mJetVertexsim= mJetEventParticle->vertex(0);

  NJetssim =mJetVertexsim->numberOfJets();
  
  //************ JET CUTS **********************

  // DETECTOR JET LOOP
  for(Int_t ijet=0; ijet<NJets; ++ijet){
    InitTreeVals();
    
    mJet =mJetVertex->jet(ijet);
    VerZ =mJetVertex->position().z();
    
    //require Abs(VerZ)<30 & JetRt<0.95
    if(TMath::Abs(VerZ)>30) continue;
    jetRt= mJet->neutralFraction();
    
    //if(jetRt >0.95) continue;
    jetEta =mJet->eta();
    if(TMath::Abs(jetEta)>1) continue;
    
    jetPt =mJet->pt(); 
    //pThist-> part all jets
    DetJetPt_all->Fill(jetPt);

    jetPhi =mJet->phi();
    jetE =mJet->E();
    jetEt =(jetE)/cosh(jetEta);
    
    UeJetArea=mJet->area();
    UeDensity=mJet->ueDensity()["OffAxisConesR060"];
    
    jetDetEta= mJet->detEta();
    
    //require abs(JetDetEta) <0.8
    if(TMath::Abs(jetDetEta)>0.8) continue;
    
    //require sum of all jet pTtracks >0.5
    SumTrackPt= mJet->sumTrackPt();
    if(SumTrackPt<0.5) continue;
    


    //************ PARTICLE JET LOOP ****************
    Float_t deltRbest=10000;
    Float_t njet_bestmatch=1000;
    Float_t particlejetPt_best=10000;
    TVector3 particlejetP_best;

    Float_t particleZh_best=10000;
    Float_t particlejT_best=10000;

    Float_t particleUeDensity_best=10000;
    Float_t particleUeA_best  =10000;


    for(Int_t partjet=0; partjet<NJetssim; ++partjet){

      mJetSim =mJetVertexsim->jet(partjet);
      VerZPart =mJetVertexsim->position().z();
      jetEtaPart=mJetSim->eta();

      
      // TEST EFFECT OF PART JET ETA CUT
      PartJetEta_all->Fill(jetEtaPart);
      //cout<<"particle jet eta: "<< jetEtaPart<<endl;

      //if(TMath::Abs(jetEtaPart)>1){
      //PartJetEta_cut->Fill(jetEtaPart);
      //}
      particlejetPt= mJetSim->pt();
      
      //pThist-> part all jets
      PartJetPt_all->Fill(particlejetPt);

      TVector3 particlejetP= mJetSim->momentum();

      jetPhiPart=mJetSim->phi();

      UeJetArea_particle=mJetSim->area();
      UeDensity_particle=mJetSim->ueDensity()["OffAxisConesR060"];


      // *** FIND BEST MATCH ***
      //get deltaRparticle here
      //this is the dist btwn particle jet/detector jet
      //sqrt( deltaR^2 +deltaeta^2)
      deltaEtaPart =(jetEta-jetEtaPart);
      deltaPhiPart =(jetPhi-jetPhiPart);
      
      deltaRjet= sqrt((deltaPhiPart*deltaPhiPart)+(deltaEtaPart*deltaEtaPart));
      
      if(deltaRjet<= 0.2){
	DeltRJet->Fill(deltaRjet);	

      //***** Match particle jet --> detector jet *****      
	//save best deltR match
	// if deltaR=bestmatch, save part jet info
	if(deltaRjet < deltRbest){
	  deltRbest=deltaRjet;
	  njet_bestmatch=partjet;
	  
	  particlejetPt_best= particlejetPt;
	  particlejetP_best= particlejetP;

	  particleUeA_best= UeJetArea_particle;
	  particleUeDensity_best= UeDensity_particle;
	}//deltaRbest loop
      }//deltaR<=0.2

    }//part jet loop

    
    //**** PARTONIC PT ****
    mPythiaEvent= mSkimEvent->mcEvent(); 
    partonicPt= mPythiaEvent->pt();

    NTracks =mJet->numberOfTracks();    


    //if no best match meets criteria, fill hist
    if(njet_bestmatch==1000){
      DetJetPt_nomatch->Fill(jetPt);
      PartJetPt_nomatch->Fill(particlejetPt);
    }

   
    //**** IF PART JET HAS BEST MATCH... ****
    else{
      // ******** PARTICLE-LEVEL EVENT: PART LOOP *******************
      //now only get particle info for best match jet

      mJetSim =mJetVertexsim->jet(njet_bestmatch);//get best match jet
      NParticlesSim=mJetSim->numberOfParticles();//get part in best match jet
      
      //loop over particles in best-match PARTICLE jet
      for(Int_t part=0; part<NParticlesSim; ++part){
	mJetParticleSim =mJetEventParticle->particle(part);
	
	TLorentzVector particleFourMom =mJetParticleSim->fourMomentum();
	TVector3 particleP =particleFourMom.Vect();
	float_t particlePt= particleFourMom.Pt();

	//cout<<"test part eta: "<< particleFourMom.Eta()<<endl;
	float_t particleEta=particleFourMom.Eta();
	float_t particlePhi=particleFourMom.Phi();
	

	


	deltaR_particlebest=1000;
	//loop over particles in best-match DETECTOR jet
	for(Int_t track=0; track<NTracks; track++){
    
	  //det track info
	  mJetTrack =mJet->track(track);
	  TrackPt =mJetTrack->pt();
	  TrackEta =mJetTrack->eta();
	  TrackPhi =mJetTrack->phi();
	  chargeSign = mJetTrack->charge();
	  nSigPi =mJetTrack->nSigmaPion();
	  nSigTOFPi= mJetTrack->nSigmaTofPion();
	  HitsFit= mJetTrack->nHitsFit();
	  

	  //dist btwn DETECTOR jet axis and track in jet
	  deltaEta=(jetEta-TrackEta);
	  deltaPhi=(jetPhi-TrackPhi);
	  deltaR=sqrt((deltaPhi*deltaPhi)+(deltaEta*deltaEta));


	  if(deltaPhi > pi){
	    deltaR= deltaR-2*pi;
	  }
	  
	  //detector jet jT & z  
	  TVector3 JetP =mJet->momentum();
	  TVector3 TrackP=mJetTrack->momentum();
	  jT= TrackP.Perp(JetP);
	  Zh=(TrackPt)/(jetPt);
	  
	  //find best match particle within best match jet
	  deltaEta_particlematch =(TrackEta-particleEta);
	  deltaPhi_particlematch =(TrackPhi-particlePhi);
	  
	
	  //cout<<"partEta: "<<particleEta<<" trackEta: "<<TrackEta<<endl;
	  //cout<<"partPhi: "<<particlePhi<<" trackPhi: "<<TrackPhi<<endl;
  
	  if(deltaPhi_particlematch > pi){
	    //cout<<"deltaPhi: "<<deltaPhi_particlematch<<endl;
	    deltaPhi_particlematch=deltaPhi_particlematch-pi;
	    //cout<<"deltaPhicorr: "<<deltaPhi_particlematch<<endl;
	    
	  }
	  if(deltaPhi_particlematch < -1*pi){
	    //cout<<"deltaPhi: "<<deltaPhi_particlematch<<endl;
	    deltaPhi_particlematch=deltaPhi_particlematch+pi;
	    //cout<<"deltaPhicorr: "<<deltaPhi_particlematch<<endl;
	    
	  }
	  
	  deltaR_particlematch=sqrt((deltaPhi_particlematch*deltaPhi_particlematch)+(deltaEta_particlematch*deltaEta_particlematch));
	  
	  //diagnostic hists to see why particle jT is weird?
	  DeltEtaParticle->Fill(deltaEta_particlematch);
	  DeltPhiParticle->Fill(deltaPhi_particlematch);
	  DeltRParticle->Fill(deltaR_particlematch);


	  float_t ratio= deltaEta_particlematch/deltaPhi_particlematch;
	  DeltEtaPhi_Ratio->Fill( ratio );	  
	  if(ratio>=10){
	    //cout<<"****"<<endl;
	    //cout<<"phi: "<<deltaPhi_particlematch<<endl;
	    //cout<<"eta: "<<deltaEta_particlematch<<endl;
	    
	    //cout<<"ratio: "<<ratio<<endl;
	  }

	  particlejT= particleP.Perp(particlejetP_best);
	  particleZh=(particlePt)/(particlejetPt_best);
	  
	  
	  // if deltaR_particlcematch=bestmatch, save info
	  if(deltaR_particlematch < deltaR_particlebest){
	    deltaR_particlebest= deltaR_particlematch;
	    
	    particleZh_best= particleZh;
	    particlejT_best= particlejT;
	  }
	  
	  
	}//det track loop
	
	// SAVE TO TREES
	PartonicpT.push_back(partonicPt);
	
	//*** DET LEVEL ***
	DeltaR.push_back(deltaR);
	JetpT.push_back(jetPt);
	PionpT.push_back(TrackPt);
	UeA.push_back(UeJetArea);
	UeRho.push_back(UeDensity);
	
	DetjT.push_back(jT);
	DetZh.push_back(Zh);


	//*** PART LEVEL ***
	
	ParticleJetpT.push_back(particlejetPt_best);
	ParticleUeA.push_back(particleUeA_best);
	ParticleUeRho.push_back(particleUeDensity_best);

	ParticleZh.push_back(particleZh_best);
	ParticlejT.push_back(particlejT_best);

    
      }//part-level loop
    }//loop over jets w a match
    
   
    //CURRENTLY SAVING INFO FOR ALL TRACKS- WANT TO ONLY SAVE FOR BEST MATCH
    /*for(Int_t track=0; track<NTracks; track++){
    
      mJetTrack =mJet->track(track);
      TrackPt =mJetTrack->pt();
      TrackEta =mJetTrack->eta();
      TrackPhi =mJetTrack->phi();
      chargeSign = mJetTrack->charge();
      nSigPi =mJetTrack->nSigmaPion();
      nSigTOFPi= mJetTrack->nSigmaTofPion();
      HitsFit= mJetTrack->nHitsFit();

      //vector projections to get jT                                       
      TVector3 JetP =mJet->momentum();
      TVector3 TrackP=mJetTrack->momentum();
    
      //detector jet jT & z  
      jT= TrackP.Perp(JetP);
      Zh=(TrackPt)/(jetPt);

      deltaEta=(jetEta-TrackEta);
      deltaPhi=(jetPhi-TrackPhi);
    
      //*** only take closest matches (smallest NPartSim # of deltaR's)
      //dist btwn jet axis and track in jet
      deltaR= sqrt((deltaPhi*deltaPhi)+(deltaEta*deltaEta));
      
      if(deltaPhi > pi){
	deltaR= deltaR-2*pi;
      }
    }//loop over tracks
    /*


    /*     
      //******** INFO TO SAVE TO TREES ********
      JetpT.push_back(jetPt);
      PartonicpT.push_back(partonicPt);
      VertexZ.push_back(VerZ);
      SumTrackpT.push_back(SumTrackPt);
      NSigTPCPi.push_back(nSigPi);
      NHitsFit.push_back(HitsFit);
      PionpT.push_back(TrackPt);
      PionEta.push_back(TrackEta);
      PionPhi.push_back(TrackPhi);
      ChargeSign.push_back(chargeSign);
      //JetMomentum.push_back(JetP.Mag());
      //PionMomentum.push_back(TrackP.Mag());
      JetEta.push_back(jetEta);
      JetPhi.push_back(jetPhi);
      JetRt.push_back(jetRt);
      DeltaR.push_back(deltaR);
      DetjT.push_back(jT);
      JetDetEta.push_back(jetDetEta);
      Eta.push_back(eta);
      NSigTOFPi.push_back(nSigTOFPi);
      UeA.push_back(UeJetArea);
      UeRho.push_back(UeDensity);
      DeltaRjet.push_back(deltaRjet);
      ParticleJetpT.push_back(particlejetPt_best);
      //ParticleZh.push_back(particleZh_best);
      //ParticlejT.push_back(particlejT_best);
      //    ParticleUeA.push_back(particleUeA_best);
      ParticleUeRho.push_back(particleUeDensity_best);
  }//det track loop
    */    

  
    //fill tree once per jet
    FFTree->Fill();
  }//NJETS LOOP
  
  return kStOK;
} //END MAKE


//_____________________________________________________________________________ /// Finish
Int_t StEmbedMaker2015::Finish(){

  FFfile->Write();
  FFfile->Close();

  return kStOK;
}//END FINISH


//_____________________________________________________________________________ /// initialize tree values

void StEmbedMaker2015::InitTreeVals(){

  JetpT.clear();
  PartonicpT.clear();
 
  PionpT.clear();
  PionEta.clear();
  PionPhi.clear();
  ChargeSign.clear();
  JetpT.clear();
  JetEta.clear();
  JetPhi.clear();
  JetRt.clear();
  DeltaR.clear();

  DetZh.clear();
  DetjT.clear();
  JetDetEta.clear();
  Eta.clear();
  VertexZ.clear();
  SumTrackpT.clear();
  NSigTOFPi.clear();
  NSigTPCPi.clear();
  NHitsFit.clear();
  
  UeRho.clear();
  UeA.clear();
  
  JetMomentum.clear();
  PionMomentum.clear();
  
  //DeltaRjet.clear();
  ParticleDeltaR.clear();
  ParticleJetpT.clear();
  ParticleZh.clear();
  ParticlejT.clear();

  
  ParticleUeRho.clear();
  ParticleUeA.clear();
}


