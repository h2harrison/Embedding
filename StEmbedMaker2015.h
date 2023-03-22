
//****************************

// class:    StEmbedMaker2015
// author:   Hannah Harrison 
// RCF user: harrison

//****************************

#ifndef STAR_StEmbedMaker2015
#define STAR_StEmbedMaker2015

#ifndef StMaker_H
#include "StThreeVectorF.hh"
#include "StMaker.h"
#endif
#include <vector>
#include "StRoot/StSpinPool/StJetSkimEvent/StPythiaEvent.h"

class TTree;
class TFile;
class TChain;
class TH1F;
class TH2F;
class StEvent;
class StJetEvent;
class StJetVertex;
class StJetCandidate;
class StJetElement;
class StJetTrack;
class StJetTower;
class StJetParticle;
class StJetTower;

class StJetSkimEvent;
class StJetSkimTrig;
class StMcEvent;

class StUeOffAxisConesEvent;
class StUeVertex;
class StUeEvent;

class StEmbedMaker2015 : public StMaker {
 private:
  //
 protected:
  //

 public:

  StEmbedMaker2015(const char *name, const char *outputfile, TChain *jetChain, TChain *skimChain, TChain *UeChain);
  virtual       ~StEmbedMaker2015();
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();


  //Root TTrees and output TFile
  TString outfile;
  TString FFRootFile;
  TFile *FFfile;
  TTree *FFTree;

  //Jet chain info
  TChain *mJetChain, *mSkimChain, *mUeChain;
  StJetEvent *mJetEvent, *mJetEventParticle;
  StJetElement *mJetElement;
  StJetParticle *mJetParticle;
  StUeEvent *mUeEvent;

  const StPythiaEvent *mPythiaEvent;
  //StUeOffAxisConesEvent *mUeEvent;
  StJetSkimEvent *mSkimEvent;
  StJetVertex *mJetVertex, *mJetVertexsim;
  StUeVertex *mUeVertex;
  StJetCandidate *mJet, *mJetSim;
  StJetTrack *mJetTrack, *mJetTrackSim;
  StJetTower *mJetTower;
  StJetParticle *mJetParticleSim;

  StJetSkimTrig *pAu_ssdmb;
  StJetSkimTrig *pp_ssdmb1;
  StJetSkimTrig *pp_ssdmb2;
  StJetSkimTrig *pp_ssdmb3;
  StJetSkimTrig *pp_ssdmb4;
  Bool_t is_pplong1,is_pplong2,is_pplong3,is_pptrans,is_pAu;
  Bool_t is_minbias;

  Int_t runNum, nEventsAnalyzed;
  Int_t NVert, NJets;
  Int_t NVertsim, NJetssim, NParticlesSim;

  TH1F* NumEvts;
  TH1F* DeltRJet, *DeltRJet_best;

  TH1F* DeltRParticle, *DeltRParticle_best;
  TH1F* DeltEtaParticle, *DeltPhiParticle;
  TH1F* DeltEtaPhi_Ratio;

  TH1F* DetJetPt_all, *PartJetPt_all;
  TH1F* DetJetPt_nomatch, *PartJetPt_nomatch;
  TH1F* DetJetPt_match, *PartJetPt_match;

  TH1F* PartJetEta_all, *PartJetEta_cut;
  Float_t partonicPt;


  // DETECTOR JET INFO
  Float_t VerZ;
  Float_t jetPt, jetEt;
  Float_t jetEta, jetPhi, jetE;
  Float_t jetRt; // (jet calorimeter E)/(total jet E)
  Float_t Zh, jT;
  Float_t jetDetEta, eta, SumTrackPt;
  Float_t Charge, chargeSign;
  Float_t NTracks, NTow; //Jet Track in TPC
  Float_t TrackPt, TrackEta, TrackPhi, nSigPi; //TPC
  Float_t nSigTOFPi; //TOF
  Float_t TowEta, TowPhi, TowPt, TowEt, TowE; //Jet Tower
  Float_t deltaPhi, deltaEta;
  Float_t deltaR;//dist btwn jet axis and track in jet                         
  Float_t HitsFit; //nhits used for reconstruction of jet track                 

  // PARTICLE JET INFO
  Float_t particleNTracks, particleTrackPt;
  Float_t deltaRjet;//dist btwn detector jet and particle jet                   
  Float_t VerZPart, jetRtPart, jetEtaPart;
  Float_t jetPhiPart, particlejetPt;
  Float_t deltaPhiPart, deltaEtaPart;

  Float_t particleZh, particlejT;
  Float_t deltaEta_particlematch, deltaPhi_particlematch;
  Float_t deltaR_particlematch;//dist btwn det part and particle part
  Float_t deltaR_particlebest;

  //UE info -- PARTICLE
  Float_t UeJetArea_particle;
  float UeDensity_particle;

  //UE info -- DETECTOR
  Float_t UeJetArea;
  float UeDensity;
  Float_t UeTrackPt;

  //SAVE TO OUTPUT TREE
  std:: vector <Float_t> PartonicpT;
  std:: vector <Float_t> PionpT;
  std:: vector <Float_t> PionMomentum;
  std:: vector <Float_t> PionEta;
  std:: vector <Float_t> PionPhi;
  std:: vector <Float_t> ChargeSign;

  std:: vector <Float_t> JetpT;
  std:: vector <Float_t> JetEta;
  std:: vector <Float_t> JetPhi;
  std:: vector <Float_t> JetRt;
  std:: vector <Float_t> DeltaR;

  std:: vector <Float_t> DetZh;
  std:: vector <Float_t> DetjT;
  std:: vector <Float_t> JetMomentum;
  std:: vector <Float_t> JetDetEta;
  std:: vector <Float_t> Eta;//global eta                                                               
  std:: vector <Float_t> VertexZ;
  std:: vector <Float_t> SumTrackpT;
  std:: vector <Float_t> NSigTOFPi;
  std:: vector <Float_t> NSigTPCPi; //pion cut                                 
  std:: vector <Float_t> NHitsFit;  //pion cut                                                          
  std:: vector <Float_t> UeRho;
  std:: vector <Float_t> UeA;

  std:: vector <Float_t> DeltaRjet;
  std:: vector <Float_t> ParticleDeltaR;
  std:: vector <Float_t> ParticleJetpT;
  std:: vector <Float_t> ParticleZh;
  std:: vector <Float_t> ParticlejT;

  std:: vector <Float_t> ParticleUeRho;
  std:: vector <Float_t> ParticleUeA;


  struct Embed{
    //***** used but not saved in tree *************************
    Float_t totalNumTracks;
    Float_t totalNumTow;
    Float_t totalNumJets;
    Float_t totalNumVert;
    Float_t totalNumEvts; //tot #events in jet
  };
  

  //***** declare members of struct ************
  Embed EmbedData;

  //Additional Functions                                                       
  void InitTreeVals();
  //inline void CountNumEvts(FF &mode){if(is_minbias==true){mode.totalNumEvts++;}}
  
  //******************************************************                        /// Displayed on session exit, leave it as-is please ...                                                                                
  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $ $Id: StEmbedMaker2015.h,  Exp $ built " __DATE__ " " __TIME__ ;
    return cvs;
  }
  ClassDef(StEmbedMaker2015,0)   //StAF chain virtual base class for Makers                                                               
    };

#endif
