// -*- C++ -*-
//
// Package:    ll_ntuple_maker
// Class:      ll_ntuple_maker
// 
/**\class ll_ntuple_maker ll_ntuple_maker.cc ntuple_maker/ntuple_maker/plugins/ll_ntuple_maker.cc

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
#include <limits>

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

#include "ntuple_maker/ntuple_maker/interface/enum_definition.h"
#include "ntuple_maker/ntuple_maker/interface/lepton_ids.h"
#include "ntuple_maker/ntuple_maker/interface/triggers.h"

#include "ntuple_maker/ntuple_maker/interface/lhe_and_gen.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//
// class declaration
//

const Float_t z_mass = 91.18800;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class ll_ntuple_maker : public edm::EDAnalyzer {
public:

  
      explicit ll_ntuple_maker(const edm::ParameterSet&);
      ~ll_ntuple_maker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
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
  edm::EDGetTokenT< edm::TriggerResults > triggerResultsToken_;
  edm::EDGetTokenT< pat::TriggerObjectStandAloneCollection > triggerObjectToken_;
  edm::EDGetTokenT< pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<double> rhoHLTElectronSelectionToken_;

  TH1F * n_events_run_over;
  UInt_t flags;
  ULong64_t event;
  UInt_t run;
  UInt_t lumi;
  UInt_t nvtx;
  TTree * tree;

  Float_t maxbtagevent;

  Float_t metpt;
  Float_t metphi;
  Float_t metsumet;
  Float_t metgenmetpt;
  //Float_t metptshiftup;
  //Float_t metptshiftdown;
  LorentzVector lep1;
  LorentzVector lep2;
  Int_t lep1id;
  Int_t lep2id;
  Int_t lep1q;
  Int_t lep2q;
  Int_t n_pu_interactions;
  lhe_and_gen lhe_and_gen_object; //separate the part that runs over the generator and lhe information

  bool syscalcinfo_;
  bool lheinfo_;

  bool isMC_;

  std::string which_triggers_;

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
ll_ntuple_maker::ll_ntuple_maker(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  triggerResultsToken_(consumes< edm::TriggerResults >(edm::InputTag("TriggerResults","","HLT"))),
  triggerObjectToken_( consumes< pat::TriggerObjectStandAloneCollection >(edm::InputTag("selectedPatTrigger"))),
  //  lheRunInfoLabel_(iConfig.getParameter<edm::InputTag>("lheruninfo")),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  pileupSummaryToken_(consumes <std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileup_summary"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  rhoHLTElectronSelectionToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoHLTElectronSelection"))),
  syscalcinfo_(iConfig.getUntrackedParameter<bool>("syscalcinfo")),
  lheinfo_(iConfig.getUntrackedParameter<bool>("lheinfo")),
  isMC_(iConfig.getUntrackedParameter<bool>("isMC")),
  which_triggers_(iConfig.getUntrackedParameter<std::string>("which_triggers"))
{
  //now do what ever initialization is needed

  //if lhe_info_ is false, syscalcinfo_ should also be false
  assert(!syscalcinfo_ || lheinfo_);

  //lhe_and_gen_object.isMC_ = isMC_;
  //lhe_and_gen_object.prunedGenToken_ = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedgenparticles"));
  //lhe_and_gen_object.packedGenToken_ = consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedgenparticles"));
  //lhe_and_gen_object.lheEvtToken_ = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheevent"));
  //lhe_and_gen_object.syscalcinfo_ = syscalcinfo_;
  //lhe_and_gen_object.lheinfo_ = lheinfo_;
  //lhe_and_gen_object.lheRunInfoLabel_ = iConfig.getParameter<edm::InputTag>("lheruninfo");
  
}


ll_ntuple_maker::~ll_ntuple_maker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ll_ntuple_maker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //if (iEvent.eventAuxiliary().luminosityBlock() != 11 || iEvent.eventAuxiliary().event() != 2092)
  //  return;

  flags = 0;

  edm::Handle<double> rhoHandle;

  iEvent.getByToken(rhoToken_,rhoHandle);

  float rho    =  *rhoHandle;

  edm::Handle<double> rhoHLTElectronSelectionHandle;

  iEvent.getByToken(rhoHLTElectronSelectionToken_,rhoHLTElectronSelectionHandle);

  float rhoHLTElectronSelection  =  *rhoHLTElectronSelectionHandle;

  if (isMC_)
    n_events_run_over->Fill(0.5);

  if (isMC_){

    edm::Handle< std::vector<PileupSummaryInfo> > pileupSummaryHandle;

    iEvent.getByToken(pileupSummaryToken_,pileupSummaryHandle);

    for(const auto & pu : *pileupSummaryHandle)
      {

	if (pu.getBunchCrossing() == 0)
	  n_pu_interactions = pu.getTrueNumInteractions();
      }

  }


     edm::Handle<pat::TauCollection> taus;
     iEvent.getByToken(tauToken_, taus);

     for (const pat::Tau &tau : *taus) {
       if (tau.pt() < 20) continue;

       if (fabs(tau.eta()) > 2.4) continue;

       if(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"))
	 flags = flags | WLLJJVetoV3;
     }

     edm::Handle<pat::PackedCandidateCollection> pfs;
     iEvent.getByToken(pfToken_, pfs);

     for (unsigned int i = 0; i < pfs->size(); ++i) {
       for (unsigned int j = i+1; j < pfs->size(); ++j) {
	 const pat::PackedCandidate &pf1 = (*pfs)[i];
	 const pat::PackedCandidate &pf2 = (*pfs)[j];

	 if ( ! (abs(pf1.pdgId()) == 11 || abs(pf1.pdgId()) == 13 || abs(pf1.pdgId()) == 15 || abs(pf2.pdgId()) == 11 || abs(pf2.pdgId()) == 13 || abs(pf2.pdgId()) == 15) )
	   continue;

	 if ( abs(z_mass - (pf1.p4() + pf2.p4()).mass()) > 15 )
	   continue;

	 if (fabs(pf1.dz()) > 0.1 && fabs(pf2.dz()) > 0.1)
	   continue;

	 double charged = 0, neutral = 0, pileup  = 0;

	 //calculate an isolation variable for pf1
	 for (unsigned int k = 0; k < pfs->size(); ++k) {

	   const pat::PackedCandidate &pf3 = (*pfs)[k];

	   if(i == k)
	     continue;

	   if (deltaR(pf3,pf1) > 0.3) continue;

	   if (pf3.charge() == 0) {

	     if (pf3.pt() > 0.5) neutral += pf3.pt();
	     
	   } else if (pf3.fromPV() >= 2) {
	     
	     charged += pf3.pt();
	     
	   } else {
	     
	     if (pf3.pt() > 0.5) pileup += pf3.pt();

	   }

	 }

	 double iso_pf1 = charged + std::max(0.0, neutral-0.5*pileup);

	 charged = 0; neutral = 0; pileup  = 0;

	 //calculate an isolation variable for pf1
	 for (unsigned int k = 0; k < pfs->size(); ++k) {

	   const pat::PackedCandidate &pf3 = (*pfs)[k];

	   if(j == k)
	     continue;

	   if (deltaR(pf3,pf2) > 0.3) continue;

	   if (pf3.charge() == 0) {

	     if (pf3.pt() > 0.5) neutral += pf3.pt();
	     
	   } else if (pf3.fromPV() >= 2) {
	     
	     charged += pf3.pt();
	     
	   } else {
	     
	     if (pf3.pt() > 0.5) pileup += pf3.pt();

	   }

	 }

	 double iso_pf2 = charged + std::max(0.0, neutral-0.5*pileup);


     if (iso_pf2/pf2.pt() > 0.1 || iso_pf1/pf1.pt() > 0.1)
       continue;

	 flags = flags | WLLJJVetoV4;

       }
     }     



  edm::Handle< edm::TriggerResults> triggerResultsHandle;

  iEvent.getByToken(triggerResultsToken_,triggerResultsHandle);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResultsHandle);

  if (! trigger_fired(names, triggerResultsHandle,"doublelepton"))
    return;

  /*

  iEvent.getByToken(triggerObjectToken_,triggerObjectHandle);

  for (pat::TriggerObjectStandAlone obj : *triggerObjectHandle) { 

    obj.unpackPathNames(names);

    for (unsigned h = 0; h < obj.filterIds().size(); ++h){

      std::cout << "obj.filterIds()[h] = " << obj.filterIds()[h] << std::endl;

    }

  }

  */

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);

  maxbtagevent = std::numeric_limits<Float_t>::min();
   
  for (const pat::Jet &j : *jets) {
    if (j.pt() < 20) continue;

    //std::cout << j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;

    if (std::max(0.f,j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")) > maxbtagevent)
      maxbtagevent = std::max(0.f,j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

    //see here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X
    //if( std::max(0.f,j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")) > 0.814)
    //  return;

    
    
  }

  std::vector<UInt_t> veryloose_muon_indices;
  std::vector<UInt_t> veryloose_electron_indices;

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

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);

   for (UInt_t i = 0; i < muons->size(); i++){
     for (UInt_t j = i+1; j < muons->size(); j++){

       if ((*muons)[i].pt() < 10 || (*muons)[j].pt() < 10)
	 continue;

       if ((*muons)[j].charge() != (*muons)[i].charge() && abs(((*muons)[i].p4() + (*muons)[j].p4()).mass() - z_mass)  < 15){

	 flags = flags | WLLJJVetoV1;
	 
       }

     }
   }

   for(UInt_t i = 0; i < electrons->size(); i++){
     for(UInt_t j = i+1; j < electrons->size(); j++){

       if ((*electrons)[i].pt() < 10 || (*electrons)[j].pt() < 10)
	 continue;

       if ((*electrons)[j].charge() != (*electrons)[i].charge() && abs(((*electrons)[i].p4() + (*electrons)[j].p4()).mass() - z_mass)  < 15){

	 flags = flags | WLLJJVetoV1;
	 
       }
     }
   }

   for(UInt_t i = 0; i < electrons->size(); i++){

     if( (*electrons)[i].pt() < 10) 
       continue;

     if (!passVeryLooseElectronSelection((*electrons)[i],PV,rho,rhoHLTElectronSelection))
       continue;

     veryloose_electron_indices.push_back(i);

   }

   for(UInt_t i = 0; i < muons->size(); i++){
     for(UInt_t j = i+1; j < muons->size(); j++){

       if ((*muons)[i].pt() < 10 || (*muons)[j].pt() < 10)
	 continue;

       if ((*muons)[j].charge() != (*muons)[i].charge() && abs(((*muons)[i].p4() + (*muons)[j].p4()).mass() - z_mass)  < 15){

	 flags = flags | WLLJJVetoV6;
	 
       }
     }
   }

   for(UInt_t i = 0; i < muons->size(); i++){

     if ((*muons)[i].pt() < 10)
       continue;

     if (! passVeryLooseMuonSelection( (*muons)[i],PV ) )
       continue;

     veryloose_muon_indices.push_back(i);

    }

   if(veryloose_muon_indices.size() >= 2){

     UInt_t i1 = veryloose_muon_indices[0];
     UInt_t i2 = veryloose_muon_indices[1];

     if (passLooseMuonSelectionV1((*muons)[i1],PV) )
       flags = flags | Lep1LooseSelectionV1;

     if (passLooseMuonSelectionV1((*muons)[i2],PV) )
       flags = flags | Lep2LooseSelectionV1;

     if (passTightMuonSelectionV1((*muons)[i1],PV)) 
       flags = flags | Lep1TightSelectionV1;

     if (passTightMuonSelectionV1((*muons)[i2],PV)) 
       flags = flags | Lep2TightSelectionV1;

     if (passTightMuonSelectionV2((*muons)[i1],PV)) 
       flags = flags | Lep1TightSelectionV2;

     if (passTightMuonSelectionV2((*muons)[i2],PV)) 
       flags = flags | Lep2TightSelectionV2;

     if (passTightMuonSelectionV3((*muons)[i1],PV)) 
       flags = flags | Lep1TightSelectionV3;

     if (passTightMuonSelectionV3((*muons)[i2],PV)) 
       flags = flags | Lep2TightSelectionV3;
     
     lep1 = (*muons)[i1].p4();
     lep1q = (*muons)[i1].charge();
     lep1id = (*muons)[i1].pdgId();

     lep2 = (*muons)[i2].p4();
     lep2q = (*muons)[i2].charge();
     lep2id = (*muons)[i2].pdgId();
     
     for(UInt_t i = 0; i < muons->size(); i++){

       if (i == i1 || i == i2)
	 continue;

       if ( passWLLJJVetoMuonId( (*muons)[i],PV ) ) 
	 flags = flags | WLLJJVetoV2;

       if (passSoftMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 3)
	 flags = flags | WLLJJVetoV5;

     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if ( passWLLJJVetoElectronId( (*electrons)[i], PV , rho, rhoHLTElectronSelection) ) 
	 flags = flags | WLLJJVetoV2;

     }


   }
   else if (veryloose_muon_indices.size() >=1 && veryloose_electron_indices.size() >= 1){

     UInt_t im = veryloose_muon_indices[0];
     UInt_t ie = veryloose_electron_indices[0];

     if (passLooseMuonSelectionV1((*muons)[im],PV))
       flags = flags | Lep1LooseSelectionV1;

     if (passTightMuonSelectionV1((*muons)[im],PV)) 
       flags = flags | Lep1TightSelectionV1;

     if (passTightMuonSelectionV2((*muons)[im],PV)) 
       flags = flags | Lep1TightSelectionV2;

     if (passTightMuonSelectionV3((*muons)[im],PV) ) 
       flags = flags | Lep1TightSelectionV3;

     lep1 = (*muons)[im].p4();
     lep1q = (*muons)[im].charge();
     lep1id = (*muons)[im].pdgId();

     if (passTightElectronSelectionV2((*electrons)[ie], PV,rho))
       flags = flags | Lep2TightSelectionV1;
     if (passLooseElectronSelectionV1((*electrons)[ie],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep2LooseSelectionV1;
     if (passLooseElectronSelectionV2((*electrons)[ie],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep2LooseSelectionV2;
       
     lep2= (*electrons)[ie].p4();
     lep2q = (*electrons)[ie].charge();
     lep2id = (*electrons)[ie].pdgId();
     
     for(UInt_t i = 0; i < muons->size(); i++){

       if (i == im)
	 continue;

       if ( passWLLJJVetoMuonId( (*muons)[i],PV ) ) 
	 flags = flags | WLLJJVetoV2;

       if (passSoftMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 3)
	 flags = flags | WLLJJVetoV5;

     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if (i == ie)
	 continue;

       if ( passWLLJJVetoElectronId( (*electrons)[i], PV , rho, rhoHLTElectronSelection) ) 
	 flags = flags | WLLJJVetoV2;


     }

     
   } else if (veryloose_electron_indices.size() >= 2){

     UInt_t i1 = veryloose_electron_indices[0];
     UInt_t i2 = veryloose_electron_indices[1];

     lep1= (*electrons)[i1].p4();
     lep1q = (*electrons)[i1].charge();
     lep1id = (*electrons)[i1].pdgId();

     lep2= (*electrons)[i2].p4();
     lep2q = (*electrons)[i2].charge();
     lep2id = (*electrons)[i2].pdgId();

     if (passTightElectronSelectionV2((*electrons)[i1],PV,rho))
	 flags = flags | Lep1TightSelectionV1;
     if (passTightElectronSelectionV2((*electrons)[i2],PV,rho))
	 flags = flags | Lep2TightSelectionV1;
     if (passLooseElectronSelectionV1((*electrons)[i1],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep1LooseSelectionV1;
     if (passLooseElectronSelectionV1((*electrons)[i2],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep2LooseSelectionV1;
     if (passLooseElectronSelectionV2((*electrons)[i1],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep1LooseSelectionV2;
     if (passLooseElectronSelectionV2((*electrons)[i2],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep2LooseSelectionV2;

     for(UInt_t i = 0; i < muons->size(); i++){

       if ( passWLLJJVetoMuonId( (*muons)[i],PV ) ) 
	 flags = flags | WLLJJVetoV2;

       if (passSoftMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 3)
	 flags = flags | WLLJJVetoV5;

     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if (i == i1 || i == i2)
	 continue;

       if ( passWLLJJVetoElectronId( (*electrons)[i],PV, rho , rhoHLTElectronSelection) ) 
	 flags = flags | WLLJJVetoV2;

     }


   } else 
     return;

   //   lhe_and_gen_object.analyze(iEvent,lep1,lep2);

   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();

   metpt = met.pt();
   metphi = met.phi();
   metsumet = met.sumEt();

   if (isMC_)
     metgenmetpt = met.genMET()->pt();
   else
     metgenmetpt = 0;

   //metptshiftup = met.shiftedPt(pat::MET::JetEnUp);

   //metptshiftdown = met.shiftedPt(pat::MET::JetEnDown);

   tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
ll_ntuple_maker::beginJob()
{

  edm::Service<TFileService> fs;

  if (isMC_)
    n_events_run_over= fs->make<TH1F>("n_events_run_over","n_events_run_over",1,0,1);

  //  lhe_and_gen_object.initrwgt_header_tree_ = fs->make<TTree >("initrwgt_header","initrwgt_header");

  //  lhe_and_gen_object.slha_header_tree_ = fs->make<TTree >("slha_header","slha_header");

  tree = fs->make<TTree>( "events"  , "events");

  tree->Branch("flags",&flags);

  tree->Branch("event",&event);
  tree->Branch("lumi",&lumi);
  tree->Branch("run",&run);

  tree->Branch("metpt",&metpt);
  tree->Branch("metphi",&metphi);
  tree->Branch("metsumet",&metsumet);
  tree->Branch("metgenmetpt",&metgenmetpt);
  //tree->Branch("metptshiftup",&metptshiftup);
  //tree->Branch("metptshiftdown",&metptshiftdown);
  tree->Branch("lep1",&lep1);
  tree->Branch("lep2",&lep2);
  tree->Branch("nvtx",&nvtx);
  tree->Branch("lep1q",&lep1q);
  tree->Branch("lep1id",&lep1id);
  tree->Branch("lep2q",&lep2q);
  tree->Branch("lep2id",&lep2id);

  tree->Branch("n_pu_interactions",&n_pu_interactions);

  tree->Branch("maxbtagevent",&maxbtagevent);

  //  lhe_and_gen_object.defineBranches(tree);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ll_ntuple_maker::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------


void 
ll_ntuple_maker::beginRun(edm::Run const& iRun, edm::EventSetup const&)
{

  //  lhe_and_gen_object.beginRun(iRun);

}


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
ll_ntuple_maker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ll_ntuple_maker);
