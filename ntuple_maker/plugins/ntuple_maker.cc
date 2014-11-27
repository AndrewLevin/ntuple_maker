// -*- C++ -*-
//
// Package:    ntuple_maker
// Class:      ntuple_maker
// 
/**\class ntuple_maker ntuple_maker.cc ntuple_maker/ntuple_maker/plugins/ntuple_maker.cc

 Description: makes smaller ntuples from miniaod

 Implementation:
 writes a flat ttree
*/
//
// Original Author:  Andrew Michael Levin
//         Created:  Fri, 25 Jul 2014 16:58:43 GMT
// $Id$
//
//



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/LorentzVector.h"

//
// class declaration
//



class ntuple_maker : public edm::EDAnalyzer {
public:
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
  
      explicit ntuple_maker(const edm::ParameterSet&);
      ~ntuple_maker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  enum Cuts {
    Lep1FullSelectionV1  = 1UL<<1, 
    Lep1FullSelectionV2  = 1UL<<2, 
    Lep1FullSelectionV3  = 1UL<<3,
    Lep1FullSelectionV4  = 1UL<<4, 
    Lep1FullSelectionV5  = 1UL<<5, 
    Lep2FullSelectionV1  = 1UL<<6, 
    Lep2FullSelectionV2  = 1UL<<7, 
    Lep2FullSelectionV3  = 1UL<<8, 
    Lep2FullSelectionV4  = 1UL<<9, 
    Lep2FullSelectionV5  = 1UL<<10 
  };


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<pat::TauCollection> tauToken_;
  edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;


  TH1F * n_events_run_over;
  UInt_t cuts;
  UInt_t event;
  UInt_t run;
  UInt_t lumi;
  UInt_t nvtx;
  TTree * tree;
  Float_t jetpt;
  Float_t jet1pujetid;
  Float_t jet1btag;
  Float_t jet1btagincl;
  Float_t jet2pujetid;
  Float_t jet2btag;
  Float_t jet2btagincl;
  Float_t metpt;
  Float_t metphi;
  Float_t metsumet;
  Float_t metgenmetpt;
  Float_t metptshiftup;
  Float_t metptshiftdown;
  LorentzVector jet1;
  LorentzVector jet2;
  LorentzVector lep1;
  LorentzVector lep2;
  Int_t lep1id;
  Int_t lep2id;
  Int_t lep1q;
  Int_t lep2q;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ntuple_maker::ntuple_maker(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
{
  //now do what ever initialization is needed

}


ntuple_maker::~ntuple_maker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ntuple_maker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  n_events_run_over->Fill(0.5);

  std::vector<UInt_t> loose_muon_indices;
  std::vector<UInt_t> loose_electron_indices;

  cuts = 0;

   using namespace edm;

   run=iEvent.eventAuxiliary().run(); 
   lumi=iEvent.eventAuxiliary().luminosityBlock();
   event=iEvent.eventAuxiliary().event();

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if (vertices->empty()) return; // skip the event if no PV found
   const reco::Vertex &PV = vertices->front();

   nvtx=vertices->size();

   edm::Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);

   //   for (const pat::Electron &el : *electrons) {
   for(UInt_t i = 0; i < electrons->size(); i++){

     if( (*electrons)[i].pt() < 10 || (*electrons)[i].chargedHadronIso()/(*electrons)[i].pt() > 0.3) 
       continue;

     loose_electron_indices.push_back(i);

     //if (el.chargedHadronIso()/el.pt() > 0.3)
     //  continue;

   }

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);

   //Bool_t found_muon = kFALSE;

   

   for(UInt_t i = 0; i < muons->size(); i++){
     if ((*muons)[i].pt() < 10 || !(*muons)[i].isLooseMuon() || (*muons)[i].chargedHadronIso()/(*muons)[i].pt() > 0.3)
       continue;

     //std::cout << "(*muons)[i].pt() = " << (*muons)[i].pt() << std::endl;

     loose_muon_indices.push_back(i);

    }

   if(loose_muon_indices.size() >= 2){

     UInt_t i1 = loose_muon_indices[0];
     UInt_t i2 = loose_muon_indices[1];

     Float_t relative_isolation_1 = ((*muons)[i1].pfIsolationR04().sumChargedHadronPt+ std::max(0.0,(*muons)[i1].pfIsolationR04().sumNeutralHadronEt + (*muons)[i1].pfIsolationR04().sumPhotonEt - 0.5 * (*muons)[i1].pfIsolationR04().sumChargedHadronPt))/(*muons)[i1].pt();
     Float_t relative_isolation_2 = ((*muons)[i2].pfIsolationR04().sumChargedHadronPt+ std::max(0.0,(*muons)[i2].pfIsolationR04().sumNeutralHadronEt + (*muons)[i2].pfIsolationR04().sumPhotonEt - 0.5 * (*muons)[i2].pfIsolationR04().sumChargedHadronPt))/(*muons)[i2].pt();

     if ((*muons)[i1].isTightMuon(PV) && relative_isolation_1 < 0.12) 
       cuts = cuts | Lep1FullSelectionV1;

     if ((*muons)[i2].isTightMuon(PV) && relative_isolation_2 < 0.12) 
       cuts = cuts | Lep2FullSelectionV1;

     if ((*muons)[i1].isTightMuon(PV) && relative_isolation_1 < 0.2) 
       cuts = cuts | Lep1FullSelectionV2;

     if ((*muons)[i2].isTightMuon(PV) && relative_isolation_2 < 0.2) 
       cuts = cuts | Lep2FullSelectionV2;

     if ((*muons)[i1].isTightMuon(PV) && relative_isolation_1 < 0.3) 
       cuts = cuts | Lep1FullSelectionV3;

     if ((*muons)[i2].isTightMuon(PV) && relative_isolation_2 < 0.3) 
       cuts = cuts | Lep2FullSelectionV3;

     if ((*muons)[i1].isTightMuon(PV) && relative_isolation_1 < 0.4) 
       cuts = cuts | Lep1FullSelectionV4;

     if ((*muons)[i2].isTightMuon(PV) && relative_isolation_2 < 0.4) 
       cuts = cuts | Lep2FullSelectionV4;

     if (relative_isolation_1 < 0.2) 
       cuts = cuts | Lep1FullSelectionV5;

     if (relative_isolation_2 < 0.2) 
       cuts = cuts | Lep2FullSelectionV5;
     
     lep1 = (*muons)[i1].p4();
     lep1q = (*muons)[i1].charge();
     lep1id = (*muons)[i1].pdgId();

     lep2 = (*muons)[i2].p4();
     lep2q = (*muons)[i2].charge();
     lep2id = (*muons)[i2].pdgId();
     
   }
   else if (loose_muon_indices.size() >=1 && loose_electron_indices.size() >= 1){

     UInt_t im = loose_muon_indices[0];
     UInt_t ie = loose_electron_indices[0];


     Float_t relative_isolation_1 = ((*muons)[im].pfIsolationR04().sumChargedHadronPt+ std::max(0.0,(*muons)[im].pfIsolationR04().sumNeutralHadronEt + (*muons)[im].pfIsolationR04().sumPhotonEt - 0.5 * (*muons)[im].pfIsolationR04().sumChargedHadronPt))/(*muons)[im].pt();

     if ((*muons)[im].isTightMuon(PV) && relative_isolation_1 < 0.12) 
       cuts = cuts | Lep1FullSelectionV1;

     if ((*muons)[im].isTightMuon(PV) && relative_isolation_1 < 0.2) 
       cuts = cuts | Lep1FullSelectionV2;

     if ((*muons)[im].isTightMuon(PV) && relative_isolation_1 < 0.3) 
       cuts = cuts | Lep1FullSelectionV3;

     if ((*muons)[im].isTightMuon(PV) && relative_isolation_1 < 0.4) 
       cuts = cuts | Lep1FullSelectionV4;

     if (relative_isolation_1 < 0.2) 
       cuts = cuts | Lep1FullSelectionV5;

     lep1 = (*muons)[im].p4();
     lep1q = (*muons)[im].charge();
     lep1id = (*muons)[im].pdgId();

     reco::GsfElectron::PflowIsolationVariables pfIso = (*electrons)[ie].pfIsolationVariables();

     Float_t absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
     Float_t relIsoWithDBeta = absiso/(*electrons)[ie].pt();
     Float_t ooEmooP = 0;

     if( (*electrons)[ie].ecalEnergy() == 0 ){
       std::cout << "Electron energy is zero!" << std::endl;
       ooEmooP = 1e30;
     }
     else if (!std::isfinite((*electrons)[ie].ecalEnergy())){
       std::cout << "Electron energy is not finite!" << std::endl;
       ooEmooP = 1e30;
     }
     else
       ooEmooP = fabs(1.0/(*electrons)[ie].ecalEnergy() - (*electrons)[ie].eSuperClusterOverP()/(*electrons)[ie].ecalEnergy() );

     //loose working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if((*electrons)[ie].superCluster()->eta() < 2.5 && (*electrons)[ie].superCluster()->eta() > 1.479 ){
       if(
	  (fabs((*electrons)[ie].deltaEtaSuperClusterTrackAtVtx()) < 0.0124)
	  &&
	  ( fabs((*electrons)[ie].deltaPhiSuperClusterTrackAtVtx()) < 0.0642)
	  &&
	  ((*electrons)[ie].full5x5_sigmaIetaIeta() < 0.035)
	  &&
	  ( (*electrons)[ie].hcalOverEcal() < 0.1115)
	  &&
	  ( fabs((-1) * (*electrons)[ie].gsfTrack()->dxy(PV.position())) < 0.098)
	  &&
	  (  fabs((*electrons)[ie].gsfTrack()->dz( PV.position() )) < 0.9187)
	  &&
	  (fabs(ooEmooP) < 0.1443)
	  &&
	  (relIsoWithDBeta < 0.3529)
	  &&
	  ((*electrons)[ie].passConversionVeto())
	  &&
	  ((*electrons)[ie].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep2FullSelectionV3;
     } else if ((*electrons)[ie].superCluster()->eta() < 1.479) {
       if(
	  (fabs((*electrons)[ie].deltaEtaSuperClusterTrackAtVtx()) < 0.0181)
	  &&
	  ( fabs((*electrons)[ie].deltaPhiSuperClusterTrackAtVtx()) < 0.0936)
	  &&
	  ((*electrons)[ie].full5x5_sigmaIetaIeta() < 0.0123)
	  &&
	  ( (*electrons)[ie].hcalOverEcal() < 0.141)
	  &&
	  ( fabs((-1) * (*electrons)[ie].gsfTrack()->dxy(PV.position())) < 0.0166)
	  &&
	  (  fabs((*electrons)[ie].gsfTrack()->dz( PV.position() )) < 0.54342)
	  &&
	  (fabs(ooEmooP) < 0.1353)
	  &&
	  (relIsoWithDBeta < 0.24)
	  &&
	  ((*electrons)[ie].passConversionVeto())
	  &&
	  ((*electrons)[ie].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep2FullSelectionV3;
     } 
       
     //medium working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if((*electrons)[ie].superCluster()->eta() < 2.5 && (*electrons)[ie].superCluster()->eta() > 1.479 ){
       if(
	  (fabs((*electrons)[ie].deltaEtaSuperClusterTrackAtVtx()) < 0.0108)
	  &&
	  ( fabs((*electrons)[ie].deltaPhiSuperClusterTrackAtVtx()) < 0.0455)
	  &&
	  ((*electrons)[ie].full5x5_sigmaIetaIeta() < 0.0318)
	  &&
	  ( (*electrons)[ie].hcalOverEcal() < 0.097)
	  &&
	  ( fabs((-1) * (*electrons)[ie].gsfTrack()->dxy(PV.position())) < 0.0845)
	  &&
	  (  fabs((*electrons)[ie].gsfTrack()->dz( PV.position() )) < 0.7523)
	  &&
	  (fabs(ooEmooP) < 0.1201)
	  &&
	  (relIsoWithDBeta < 0.254)
	  &&
	  ((*electrons)[ie].passConversionVeto())
	  &&
	  ((*electrons)[ie].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep2FullSelectionV1;
     } else if ((*electrons)[ie].superCluster()->eta() < 1.479) {
       if(
	  (fabs((*electrons)[ie].deltaEtaSuperClusterTrackAtVtx()) < 0.0106)
	  &&
	  ( fabs((*electrons)[ie].deltaPhiSuperClusterTrackAtVtx()) < 0.0323)
	  &&
	  ((*electrons)[ie].full5x5_sigmaIetaIeta() < 0.0107)
	  &&
	  ( (*electrons)[ie].hcalOverEcal() < 0.067)
	  &&
	  ( fabs((-1) * (*electrons)[ie].gsfTrack()->dxy(PV.position())) < 0.0131)
	  &&
	  (  fabs((*electrons)[ie].gsfTrack()->dz( PV.position() )) < 0.22310)
	  &&
	  (fabs(ooEmooP) < 0.1043)
	  &&
	  (relIsoWithDBeta < 0.2179)
	  &&
	  ((*electrons)[ie].passConversionVeto())
	  &&
	  ((*electrons)[ie].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep2FullSelectionV1;
     } 

     //tight working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if((*electrons)[ie].superCluster()->eta() < 2.5 && (*electrons)[ie].superCluster()->eta() > 1.479 ){
       if(
	  (fabs((*electrons)[ie].deltaEtaSuperClusterTrackAtVtx()) < 0.0106)
	  &&
	  ( fabs((*electrons)[ie].deltaPhiSuperClusterTrackAtVtx()) < 0.0359)
	  &&
	  ((*electrons)[ie].full5x5_sigmaIetaIeta() < 0.0305)
	  &&
	  ( (*electrons)[ie].hcalOverEcal() < 0.0835)
	  &&
	  ( fabs((-1) * (*electrons)[ie].gsfTrack()->dxy(PV.position())) < 0.0163)
	  &&
	  (  fabs((*electrons)[ie].gsfTrack()->dz( PV.position() )) < 0.5999)
	  &&
	  (fabs(ooEmooP) < 0.1126)
	  &&
	  (relIsoWithDBeta < 0.2075)
	  &&
	  ((*electrons)[ie].passConversionVeto())
	  &&
	  ((*electrons)[ie].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep2FullSelectionV2;
     } else if ((*electrons)[ie].superCluster()->eta() < 1.479) {
       if(
	  (fabs((*electrons)[ie].deltaEtaSuperClusterTrackAtVtx()) < 0.0091)
	  &&
	  ( fabs((*electrons)[ie].deltaPhiSuperClusterTrackAtVtx()) < 0.031)
	  &&
	  ((*electrons)[ie].full5x5_sigmaIetaIeta() < 0.0106)
	  &&
	  ( (*electrons)[ie].hcalOverEcal() < 0.0532)
	  &&
	  ( fabs((-1) * (*electrons)[ie].gsfTrack()->dxy(PV.position())) < 0.0126)
	  &&
	  (  fabs((*electrons)[ie].gsfTrack()->dz( PV.position() )) < 0.0116)
	  &&
	  (fabs(ooEmooP) < 0.0609)
	  &&
	  (relIsoWithDBeta < 0.1649)
	  &&
	  ((*electrons)[ie].passConversionVeto())
	  &&
	  ((*electrons)[ie].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep2FullSelectionV2;
     } 

     lep2= (*electrons)[ie].p4();
     lep2q = (*electrons)[ie].charge();
     lep2id = (*electrons)[ie].pdgId();
     
     
   } else if (loose_electron_indices.size() >= 2){

     UInt_t i1 = loose_electron_indices[0];
     UInt_t i2 = loose_electron_indices[1];

     reco::GsfElectron::PflowIsolationVariables pfIso_1 = (*electrons)[i1].pfIsolationVariables();

     Float_t absiso_1 = pfIso_1.sumChargedHadronPt + std::max(0.0 , pfIso_1.sumNeutralHadronEt + pfIso_1.sumPhotonEt - 0.5 * pfIso_1.sumPUPt );
     Float_t relIsoWithDBeta_1 = absiso_1/(*electrons)[i1].pt();
     Float_t ooEmooP_1 = 0;

     if( (*electrons)[i1].ecalEnergy() == 0 ){
       std::cout << "Electron energy is zero!" << std::endl;
       ooEmooP_1 = 1e30;
     }
     else if (!std::isfinite((*electrons)[i1].ecalEnergy())){
       std::cout << "Electron energy is not finite!" << std::endl;
       ooEmooP_1 = 1e30;
     }
     else
       ooEmooP_1 = fabs(1.0/(*electrons)[i1].ecalEnergy() - (*electrons)[i1].eSuperClusterOverP()/(*electrons)[i1].ecalEnergy() );

     //loose working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if((*electrons)[i1].superCluster()->eta() < 2.5 && (*electrons)[i1].superCluster()->eta() > 1.479 ){
       if(
	  (fabs((*electrons)[i1].deltaEtaSuperClusterTrackAtVtx()) < 0.0124)
	  &&
	  ( fabs((*electrons)[i1].deltaPhiSuperClusterTrackAtVtx()) < 0.0642)
	  &&
	  ((*electrons)[i1].full5x5_sigmaIetaIeta() < 0.035)
	  &&
	  ( (*electrons)[i1].hcalOverEcal() < 0.1115)
	  &&
	  ( fabs((-1) * (*electrons)[i1].gsfTrack()->dxy(PV.position())) < 0.098)
	  &&
	  (  fabs((*electrons)[i1].gsfTrack()->dz( PV.position() )) < 0.9187)
	  &&
	  (fabs(ooEmooP_1) < 0.1443)
	  &&
	  (relIsoWithDBeta_1 < 0.3529)
	  &&
	  ((*electrons)[i1].passConversionVeto())
	  &&
	  ((*electrons)[i1].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep1FullSelectionV3;
     } else if ((*electrons)[i1].superCluster()->eta() < 1.479) {
       if(
	  (fabs((*electrons)[i1].deltaEtaSuperClusterTrackAtVtx()) < 0.0181)
	  &&
	  ( fabs((*electrons)[i1].deltaPhiSuperClusterTrackAtVtx()) < 0.0936)
	  &&
	  ((*electrons)[i1].full5x5_sigmaIetaIeta() < 0.0123)
	  &&
	  ( (*electrons)[i1].hcalOverEcal() < 0.141)
	  &&
	  ( fabs((-1) * (*electrons)[i1].gsfTrack()->dxy(PV.position())) < 0.0166)
	  &&
	  (  fabs((*electrons)[i1].gsfTrack()->dz( PV.position() )) < 0.54342)
	  &&
	  (fabs(ooEmooP_1) < 0.1353)
	  &&
	  (relIsoWithDBeta_1 < 0.24)
	  &&
	  ((*electrons)[i1].passConversionVeto())
	  &&
	  ((*electrons)[i1].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep1FullSelectionV3;
     } 

     //medium working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if((*electrons)[i1].superCluster()->eta() < 2.5 && (*electrons)[i1].superCluster()->eta() > 1.479 ){
       if(
	  (fabs((*electrons)[i1].deltaEtaSuperClusterTrackAtVtx()) < 0.0108)
	  &&
	  ( fabs((*electrons)[i1].deltaPhiSuperClusterTrackAtVtx()) < 0.0455)
	  &&
	  ((*electrons)[i1].full5x5_sigmaIetaIeta() < 0.0318)
	  &&
	  ( (*electrons)[i1].hcalOverEcal() < 0.097)
	  &&
	  ( fabs((-1) * (*electrons)[i1].gsfTrack()->dxy(PV.position())) < 0.0845)
	  &&
	  (  fabs((*electrons)[i1].gsfTrack()->dz( PV.position() )) < 0.7523)
	  &&
	  (fabs(ooEmooP_1) < 0.1201)
	  &&
	  (relIsoWithDBeta_1 < 0.254)
	  &&
	  ((*electrons)[i1].passConversionVeto())
	  &&
	  ((*electrons)[i1].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep1FullSelectionV1;
     } else if ((*electrons)[i1].superCluster()->eta() < 1.479) {
       if(
	  (fabs((*electrons)[i1].deltaEtaSuperClusterTrackAtVtx()) < 0.0106)
	  &&
	  ( fabs((*electrons)[i1].deltaPhiSuperClusterTrackAtVtx()) < 0.0323)
	  &&
	  ((*electrons)[i1].full5x5_sigmaIetaIeta() < 0.0107)
	  &&
	  ( (*electrons)[i1].hcalOverEcal() < 0.067)
	  &&
	  ( fabs((-1) * (*electrons)[i1].gsfTrack()->dxy(PV.position())) < 0.0131)
	  &&
	  (  fabs((*electrons)[i1].gsfTrack()->dz( PV.position() )) < 0.22310)
	  &&
	  (fabs(ooEmooP_1) < 0.1043)
	  &&
	  (relIsoWithDBeta_1 < 0.2179)
	  &&
	  ((*electrons)[i1].passConversionVeto())
	  &&
	  ((*electrons)[i1].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep1FullSelectionV1;
     } 

     //tight working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if((*electrons)[i1].superCluster()->eta() < 2.5 && (*electrons)[i1].superCluster()->eta() > 1.479 ){
       if(
	  (fabs((*electrons)[i1].deltaEtaSuperClusterTrackAtVtx()) < 0.0106)
	  &&
	  ( fabs((*electrons)[i1].deltaPhiSuperClusterTrackAtVtx()) < 0.0359)
	  &&
	  ((*electrons)[i1].full5x5_sigmaIetaIeta() < 0.0305)
	  &&
	  ( (*electrons)[i1].hcalOverEcal() < 0.0835)
	  &&
	  ( fabs((-1) * (*electrons)[i1].gsfTrack()->dxy(PV.position())) < 0.0163)
	  &&
	  (  fabs((*electrons)[i1].gsfTrack()->dz( PV.position() )) < 0.5999)
	  &&
	  (fabs(ooEmooP_1) < 0.1126)
	  &&
	  (relIsoWithDBeta_1 < 0.2075)
	  &&
	  ((*electrons)[i1].passConversionVeto())
	  &&
	  ((*electrons)[i1].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep1FullSelectionV2;
     } else if ((*electrons)[i1].superCluster()->eta() < 1.479) {
       if(
	  (fabs((*electrons)[i1].deltaEtaSuperClusterTrackAtVtx()) < 0.0091)
	  &&
	  ( fabs((*electrons)[i1].deltaPhiSuperClusterTrackAtVtx()) < 0.031)
	  &&
	  ((*electrons)[i1].full5x5_sigmaIetaIeta() < 0.0106)
	  &&
	  ( (*electrons)[i1].hcalOverEcal() < 0.0532)
	  &&
	  ( fabs((-1) * (*electrons)[i1].gsfTrack()->dxy(PV.position())) < 0.0126)
	  &&
	  (  fabs((*electrons)[i1].gsfTrack()->dz( PV.position() )) < 0.0116)
	  &&
	  (fabs(ooEmooP_1) < 0.0609)
	  &&
	  (relIsoWithDBeta_1 < 0.1649)
	  &&
	  ((*electrons)[i1].passConversionVeto())
	  &&
	  ((*electrons)[i1].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep1FullSelectionV2;
     } 
     
     lep1= (*electrons)[i1].p4();
     lep1q = (*electrons)[i1].charge();
     lep1id = (*electrons)[i1].pdgId();

     reco::GsfElectron::PflowIsolationVariables pfIso_2 = (*electrons)[i2].pfIsolationVariables();

     Float_t absiso_2 = pfIso_2.sumChargedHadronPt + std::max(0.0 , pfIso_2.sumNeutralHadronEt + pfIso_2.sumPhotonEt - 0.5 * pfIso_2.sumPUPt );
     Float_t relIsoWithDBeta_2 = absiso_2/(*electrons)[i2].pt();
     Float_t ooEmooP_2 = 0;

     if( (*electrons)[i2].ecalEnergy() == 0 ){
       std::cout << "Electron energy is zero!" << std::endl;
       ooEmooP_2 = 1e30;
     }
     else if (!std::isfinite((*electrons)[i2].ecalEnergy())){
       std::cout << "Electron energy is not finite!" << std::endl;
       ooEmooP_2 = 1e30;
     }
     else
       ooEmooP_2 = fabs(1.0/(*electrons)[i2].ecalEnergy() - (*electrons)[i2].eSuperClusterOverP()/(*electrons)[i2].ecalEnergy() );


     //loose working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if((*electrons)[i2].superCluster()->eta() < 2.5 && (*electrons)[i2].superCluster()->eta() > 1.479 ){
       if(
	  (fabs((*electrons)[i2].deltaEtaSuperClusterTrackAtVtx()) < 0.0124)
	  &&
	  ( fabs((*electrons)[i2].deltaPhiSuperClusterTrackAtVtx()) < 0.0642)
	  &&
	  ((*electrons)[i2].full5x5_sigmaIetaIeta() < 0.035)
	  &&
	  ( (*electrons)[i2].hcalOverEcal() < 0.1115)
	  &&
	  ( fabs((-1) * (*electrons)[i2].gsfTrack()->dxy(PV.position())) < 0.098)
	  &&
	  (  fabs((*electrons)[i2].gsfTrack()->dz( PV.position() )) < 0.9187)
	  &&
	  (fabs(ooEmooP_2) < 0.1443)
	  &&
	  (relIsoWithDBeta_2 < 0.3529)
	  &&
	  ((*electrons)[i2].passConversionVeto())
	  &&
	  ((*electrons)[i2].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep2FullSelectionV3;
     } else if ((*electrons)[i2].superCluster()->eta() < 1.479) {
       if(
	  (fabs((*electrons)[i2].deltaEtaSuperClusterTrackAtVtx()) < 0.0181)
	  &&
	  ( fabs((*electrons)[i2].deltaPhiSuperClusterTrackAtVtx()) < 0.0936)
	  &&
	  ((*electrons)[i2].full5x5_sigmaIetaIeta() < 0.0123)
	  &&
	  ( (*electrons)[i2].hcalOverEcal() < 0.141)
	  &&
	  ( fabs((-1) * (*electrons)[i2].gsfTrack()->dxy(PV.position())) < 0.0166)
	  &&
	  (  fabs((*electrons)[i2].gsfTrack()->dz( PV.position() )) < 0.54342)
	  &&
	  (fabs(ooEmooP_2) < 0.1353)
	  &&
	  (relIsoWithDBeta_2 < 0.24)
	  &&
	  ((*electrons)[i2].passConversionVeto())
	  &&
	  ((*electrons)[i2].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep2FullSelectionV3;
     } 

     //medium working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if((*electrons)[i2].superCluster()->eta() < 2.5 && (*electrons)[i2].superCluster()->eta() > 1.479 ){
       if(
	  (fabs((*electrons)[i2].deltaEtaSuperClusterTrackAtVtx()) < 0.0108)
	  &&
	  ( fabs((*electrons)[i2].deltaPhiSuperClusterTrackAtVtx()) < 0.0455)
	  &&
	  ((*electrons)[i2].full5x5_sigmaIetaIeta() < 0.0318)
	  &&
	  ( (*electrons)[i2].hcalOverEcal() < 0.097)
	  &&
	  ( fabs((-1) * (*electrons)[i2].gsfTrack()->dxy(PV.position())) < 0.0845)
	  &&
	  (  fabs((*electrons)[i2].gsfTrack()->dz( PV.position() )) < 0.7523)
	  &&
	  (fabs(ooEmooP_2) < 0.1201)
	  &&
	  (relIsoWithDBeta_2 < 0.254)
	  &&
	  ((*electrons)[i2].passConversionVeto())
	  &&
	  ((*electrons)[i2].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep2FullSelectionV1;
     } else if ((*electrons)[i2].superCluster()->eta() < 1.479) {
       if(
	  (fabs((*electrons)[i2].deltaEtaSuperClusterTrackAtVtx()) < 0.0106)
	  &&
	  ( fabs((*electrons)[i2].deltaPhiSuperClusterTrackAtVtx()) < 0.0323)
	  &&
	  ((*electrons)[i2].full5x5_sigmaIetaIeta() < 0.0107)
	  &&
	  ( (*electrons)[i2].hcalOverEcal() < 0.067)
	  &&
	  ( fabs((-1) * (*electrons)[i2].gsfTrack()->dxy(PV.position())) < 0.0131)
	  &&
	  (  fabs((*electrons)[i2].gsfTrack()->dz( PV.position() )) < 0.22310)
	  &&
	  (fabs(ooEmooP_2) < 0.1043)
	  &&
	  (relIsoWithDBeta_2 < 0.2179)
	  &&
	  ((*electrons)[i2].passConversionVeto())
	  &&
	  ((*electrons)[i2].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep2FullSelectionV1;
     } 

     //tight working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if((*electrons)[i2].superCluster()->eta() < 2.5 && (*electrons)[i2].superCluster()->eta() > 1.479 ){
       if(
	  (fabs((*electrons)[i2].deltaEtaSuperClusterTrackAtVtx()) < 0.0106)
	  &&
	  ( fabs((*electrons)[i2].deltaPhiSuperClusterTrackAtVtx()) < 0.0359)
	  &&
	  ((*electrons)[i2].full5x5_sigmaIetaIeta() < 0.0305)
	  &&
	  ( (*electrons)[i2].hcalOverEcal() < 0.0835)
	  &&
	  ( fabs((-1) * (*electrons)[i2].gsfTrack()->dxy(PV.position())) < 0.0163)
	  &&
	  (  fabs((*electrons)[i2].gsfTrack()->dz( PV.position() )) < 0.5999)
	  &&
	  (fabs(ooEmooP_2) < 0.1126)
	  &&
	  (relIsoWithDBeta_2 < 0.2075)
	  &&
	  ((*electrons)[i2].passConversionVeto())
	  &&
	  ((*electrons)[i2].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep2FullSelectionV2;
     } else if ((*electrons)[i2].superCluster()->eta() < 1.479) {
       if(
	  (fabs((*electrons)[i2].deltaEtaSuperClusterTrackAtVtx()) < 0.0091)
	  &&
	  ( fabs((*electrons)[i2].deltaPhiSuperClusterTrackAtVtx()) < 0.031)
	  &&
	  ((*electrons)[i2].full5x5_sigmaIetaIeta() < 0.0106)
	  &&
	  ( (*electrons)[i2].hcalOverEcal() < 0.0532)
	  &&
	  ( fabs((-1) * (*electrons)[i2].gsfTrack()->dxy(PV.position())) < 0.0126)
	  &&
	  (  fabs((*electrons)[i2].gsfTrack()->dz( PV.position() )) < 0.0116)
	  &&
	  (fabs(ooEmooP_2) < 0.0609)
	  &&
	  (relIsoWithDBeta_2 < 0.1649)
	  &&
	  ((*electrons)[i2].passConversionVeto())
	  &&
	  ((*electrons)[i2].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1)
	  )
	 cuts = cuts | Lep2FullSelectionV2;
     } 

     lep2= (*electrons)[i2].p4();
     lep2q = (*electrons)[i2].charge();
     lep2id = (*electrons)[i2].pdgId();

   } else 
     return;
   
   

   /*


     found_muon=kTRUE;
     break;

     //printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
     //	    mu.pt(), mu.muonBestTrack()->dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV));
   }

   

   if(!found_muon)
     return;

   Bool_t found_electron = kFALSE;

   */



   //   std::cout << "loose_electron_indices.size() = " << loose_electron_indices.size() << std::endl;
   //   std::cout << "loose_muon_indices.size() = " << loose_muon_indices.size() << std::endl;

   //std::cout << "sum = " << loose_muon_indices.size()+ loose_electron_indices.size() << std::endl;

   /*



     //printf("elec with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), lost hits %d, pass conv veto %d\n",
     //	    el.pt(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(), el.passConversionVeto());
   }

  if(!found_electron)
     return;

   edm::Handle<pat::PhotonCollection> photons;
   iEvent.getByToken(photonToken_, photons);
   for (const pat::Photon &pho : *photons) {
     if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
     printf("phot with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes)\n",
	    pho.pt(), pho.superCluster()->eta(), pho.sigmaIetaIeta(), pho.full5x5_sigmaIetaIeta());
   }


   edm::Handle<pat::TauCollection> taus;
   iEvent.getByToken(tauToken_, taus);
   for (const pat::Tau &tau : *taus) {
     if (tau.pt() < 20) continue;
     printf("tau  with pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",
	    tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());
   }

   */


   edm::Handle<pat::JetCollection> jets;
   iEvent.getByToken(jetToken_, jets);

   if (jets->size() < 2) 
     return;
   else if ((*jets)[1].pt() < 20)
     return;
   jet1=(*jets)[0].p4();
   jet2=(*jets)[1].p4();   
   jet1btag = std::max(0.f,(*jets)[0].bDiscriminator("combinedSecondaryVertexBJetTags"));
   jet1btagincl = std::max(0.f,(*jets)[0].bDiscriminator("combinedInclusiveSecondaryVertexBJetTags"));
   jet2btag = std::max(0.f,(*jets)[1].bDiscriminator("combinedSecondaryVertexBJetTags"));
   jet2btagincl = std::max(0.f,(*jets)[1].bDiscriminator("combinedInclusiveSecondaryVertexBJetTags"));
   jet1pujetid = (*jets)[0].userFloat("pileupJetId:fullDiscriminant");
   jet2pujetid = (*jets)[1].userFloat("pileupJetId:fullDiscriminant");


   for (const pat::Jet &j : *jets) {
     if (j.pt() < 20) continue;

     if( std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")) > 0.5)
       return;

   }

   /*

   edm::Handle<pat::JetCollection> fatjets;
   iEvent.getByToken(fatjetToken_, fatjets);
   for (const pat::Jet &j : *fatjets) {
     printf("AK8j with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed, %5.1f pruned, %5.1f trimmed, %5.1f filtered. CMS TopTagger %.1f\n",
            j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), j.mass(), j.userFloat("ak8PFJetsCHSPrunedLinks"), j.userFloat("ak8PFJetsCHSTrimmedLinks"), j.userFloat("ak8PFJetsCHSFilteredLinks"), j.userFloat("cmsTopTagPFJetsCHSLinksAK8"));
   }

   */
 
   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();
   //printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
   //	  met.pt(), met.phi(), met.sumEt(),
   //	  met.genMET()->pt(),
   //	  met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));

   //   printf("\n");

   metpt = met.pt();
   metphi = met.phi();
   metsumet = met.sumEt();
   metgenmetpt = met.genMET()->pt();
   metptshiftup = met.shiftedPt(pat::MET::JetEnUp);
   metptshiftdown = met.shiftedPt(pat::MET::JetEnDown);

   tree->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
ntuple_maker::beginJob()
{

  edm::Service<TFileService> fs;

  n_events_run_over= fs->make<TH1F>("n_events_run_over","n_events_run_over",1,0,1);

  tree = fs->make<TTree>( "events"  , "events");

  tree->Branch("cuts",&cuts);

  tree->Branch("event",&event);
  tree->Branch("lumi",&lumi);
  tree->Branch("run",&run);

  tree->Branch("metpt",&metpt);
  tree->Branch("metphi",&metphi);
  tree->Branch("metsumet",&metsumet);
  tree->Branch("metgenmetpt",&metgenmetpt);
  tree->Branch("metptshiftup",&metptshiftup);
  tree->Branch("metptshiftdown",&metptshiftdown);
  tree->Branch("jet1",&jet1);
  tree->Branch("jet2",&jet2);
  tree->Branch("jet1pujetid",&jet1pujetid);
  tree->Branch("jet2pujetid",&jet2pujetid);
  tree->Branch("jet1btag",&jet1btag);
  tree->Branch("jet2btag",&jet2btag);
  tree->Branch("jet1btagincl",&jet1btagincl);
  tree->Branch("jet2btagincl",&jet2btagincl);
  tree->Branch("lep1",&lep1);
  tree->Branch("lep2",&lep2);
  tree->Branch("nvtx",&nvtx);
  tree->Branch("lep1q",&lep1q);
  tree->Branch("lep1id",&lep1id);
  tree->Branch("lep2q",&lep2q);
  tree->Branch("lep2id",&lep2id);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ntuple_maker::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

/*
void 
ntuple_maker::beginRun(edm::Run const&, edm::EventSetup const&)
{



}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ntuple_maker::endRun(edm::Run const&, edm::EventSetup const&)
{



}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ntuple_maker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ntuple_maker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ntuple_maker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ntuple_maker);
