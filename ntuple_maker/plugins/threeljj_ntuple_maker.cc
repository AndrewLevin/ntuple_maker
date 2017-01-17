// -*- C++ -*-
//
// Package:    threeljj_ntuple_maker
// Class:      threeljj_ntuple_maker
// 
/**\class threeljj_ntuple_maker threeljj_ntuple_maker.cc ntuple_maker/ntuple_maker/plugins/threeljj_ntuple_maker.cc

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

#include "ntuple_maker/ntuple_maker/interface/threeljj_enum_definition.h"
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

//
// class declaration
//

const Float_t z_mass = 91.18800;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class threeljj_ntuple_maker : public edm::EDAnalyzer {
public:

  
      explicit threeljj_ntuple_maker(const edm::ParameterSet&);
      ~threeljj_ntuple_maker();

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
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<LHEEventProduct> lheEvtToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEvtToken_;
  edm::EDGetTokenT<double> rhoHLTElectronSelectionToken_;


  TH1F * n_events_run_over;
  TH1F * n_weighted_events_run_over;
  UInt_t flags;
  UInt_t event;
  UInt_t run;
  UInt_t lumi;
  UInt_t nvtx;
  TTree * tree;

  Float_t maxbtagevent;

  Float_t jetpt;
  Float_t jet1pujetid;
  Float_t jet1btag;
  Float_t jet2pujetid;
  Float_t jet2btag;
  Float_t metpt;
  Float_t metphi;
  Float_t metsumet;
  Float_t metgenmetpt;
  //Float_t metptshiftup;
  //Float_t metptshiftdown;
  LorentzVector jet1;
  LorentzVector jet2;
  LorentzVector lep1;
  LorentzVector lep2;
  LorentzVector lep3;

  Float_t lep1iso;
  Float_t lep2iso;
  Float_t lep3iso;

  Int_t lep1id;
  Int_t lep2id;
  Int_t lep3id;
  Int_t lep1q;
  Int_t lep2q;
  Int_t lep3q;
  //  lhe_and_gen lhe_and_gen_object; //separate the part that runs over the generator and lhe information

  bool syscalcinfo_;
  bool mgreweightinfo_;

  bool apply_trigger_;
  bool isMC_;

  Float_t gen_weight;

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
threeljj_ntuple_maker::threeljj_ntuple_maker(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
    triggerResultsToken_(consumes< edm::TriggerResults >(edm::InputTag("TriggerResults","","HLT"))),
  //  triggerResultsToken_(consumes< edm::TriggerResults >(edm::InputTag("TriggerResults","","HLT2"))),
  triggerObjectToken_( consumes< pat::TriggerObjectStandAloneCollection >(edm::InputTag("selectedPatTrigger"))),
  //  lheRunInfoLabel_(iConfig.getParameter<edm::InputTag>("lheruninfo")),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  lheEvtToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheevent"))),
  genEvtToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genevent"))),
  rhoHLTElectronSelectionToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoHLTElectronSelection"))),

  syscalcinfo_(iConfig.getUntrackedParameter<bool>("syscalcinfo")),
  mgreweightinfo_(iConfig.getUntrackedParameter<bool>("mgreweightinfo")),
  apply_trigger_(iConfig.getUntrackedParameter<bool>("apply_trigger")),
  isMC_(iConfig.getUntrackedParameter<bool>("isMC"))

{
  //now do what ever initialization is needed

  consumes< LHERunInfoProduct, edm::InRun > (iConfig.getParameter<edm::InputTag>("lheruninfo"));

  
}


threeljj_ntuple_maker::~threeljj_ntuple_maker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
threeljj_ntuple_maker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //if (iEvent.eventAuxiliary().luminosityBlock() != 11 || iEvent.eventAuxiliary().event() != 2092)
  //  return;

  if (isMC_){
    n_events_run_over->Fill(0.5);

    edm::Handle<GenEventInfoProduct> hGenEvt;
    iEvent.getByToken(genEvtToken_,hGenEvt);

    if (hGenEvt->weight() > 0)
      n_weighted_events_run_over->Fill(0.5,1);
    else
      n_weighted_events_run_over->Fill(0.5,-1);

    gen_weight = hGenEvt->weight();
    

  }


  edm::Handle<double> rhoHandle;

  iEvent.getByToken(rhoToken_,rhoHandle);

  float rho    =  *rhoHandle;

  edm::Handle<double> rhoHLTElectronSelectionHandle;
  
  iEvent.getByToken(rhoHLTElectronSelectionToken_,rhoHLTElectronSelectionHandle);

  float rhoHLTElectronSelection  =  *rhoHLTElectronSelectionHandle;

  flags = 0;


     edm::Handle<pat::TauCollection> taus;
     iEvent.getByToken(tauToken_, taus);


  edm::Handle< edm::TriggerResults> triggerResultsHandle;

  iEvent.getByToken(triggerResultsToken_,triggerResultsHandle);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResultsHandle);

  //  if (apply_trigger_ && ! trigger_fired(names,triggerResultsHandle,"soup"))
      if (apply_trigger_ && ! trigger_fired(names,triggerResultsHandle,"doublelepton"))
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
  std::vector<UInt_t> tight_muon_indices;
  std::vector<UInt_t> tight_electron_indices;

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


   for(UInt_t i = 0; i < electrons->size(); i++){

     if( (*electrons)[i].pt() < 20) 
       continue;

     bool should_be_cleaned = false;

     for (UInt_t j = 0; j < muons->size(); j++){
       if (reco::deltaR((*electrons)[i],(*muons)[j]) < 0.05 && passTightMuonIdV2((*muons)[j],PV))
         should_be_cleaned = true;
     }

     if (should_be_cleaned) continue;

     if (passTightElectronSelectionV2((*electrons)[i],PV,rho))
       tight_electron_indices.push_back(i);
     else  if (passVeryLooseElectronSelection((*electrons)[i],PV,rho, rhoHLTElectronSelection))
       veryloose_electron_indices.push_back(i);

   }

   for(UInt_t i = 0; i < muons->size(); i++){

     if ((*muons)[i].pt() < 20)
       continue;

     if (passTightMuonSelectionV1( (*muons)[i],PV ) )
       tight_muon_indices.push_back(i);
     else if (passVeryLooseMuonSelection( (*muons)[i],PV ) )
       veryloose_muon_indices.push_back(i);

    }

   /*

      std::cout << "veryloose_muon_indices.size() = " << veryloose_muon_indices.size() << std::endl;
      std::cout << "veryloose_electron_indices.size() = " << veryloose_electron_indices.size() << std::endl;
     std::cout << "tight_muon_indicies.size() = " << tight_muon_indices.size() << std::endl;
     std::cout << "tight_electron_indicies.size() = " << tight_electron_indices.size() << std::endl;

   */

   if(tight_muon_indices.size() >= 3){

     UInt_t i1 = tight_muon_indices[0];
     UInt_t i2 = tight_muon_indices[1];
     UInt_t i3 = tight_muon_indices[2];

     if (passLooseMuonSelectionV5((*muons)[i1],PV) )
       flags = flags | Lep1LooseSelectionV5;

     if (passLooseMuonSelectionV5((*muons)[i2],PV) )
       flags = flags | Lep2LooseSelectionV5;

     if (passLooseMuonSelectionV5((*muons)[i3],PV) )
       flags = flags | Lep3LooseSelectionV5;

     if (passTightMuonSelectionV1((*muons)[i1],PV)) 
       flags = flags | Lep1TightSelectionV1;

     if (passTightMuonSelectionV1((*muons)[i2],PV)) 
       flags = flags | Lep2TightSelectionV1;

     if (passTightMuonSelectionV1((*muons)[i3],PV)) 
       flags = flags | Lep3TightSelectionV1;

     if (passTightMuonSelectionV2((*muons)[i1],PV)) 
       flags = flags | Lep1TightSelectionV2;

     if (passTightMuonSelectionV2((*muons)[i2],PV)) 
       flags = flags | Lep2TightSelectionV2;

     if (passTightMuonSelectionV2((*muons)[i3],PV)) 
       flags = flags | Lep3TightSelectionV2;

     if (passTightMuonSelectionV3((*muons)[i1],PV)) 
       flags = flags | Lep1TightSelectionV3;

     if (passTightMuonSelectionV3((*muons)[i2],PV)) 
       flags = flags | Lep2TightSelectionV3;

     if (passTightMuonSelectionV3((*muons)[i3],PV)) 
       flags = flags | Lep3TightSelectionV3;
     
     lep1 = (*muons)[i1].p4();
     lep1q = (*muons)[i1].charge();
     lep1id = (*muons)[i1].pdgId();
     lep1iso = muon_isolation((*muons)[i1],PV);

     lep2 = (*muons)[i2].p4();
     lep2q = (*muons)[i2].charge();
     lep2id = (*muons)[i2].pdgId();
     lep2iso = muon_isolation((*muons)[i2],PV);

     lep3 = (*muons)[i3].p4();
     lep3q = (*muons)[i3].charge();
     lep3id = (*muons)[i3].pdgId();
     lep3iso = muon_isolation((*muons)[i3],PV);
     
   }
   else if (tight_muon_indices.size() >=2 && tight_electron_indices.size() >= 1){

     UInt_t im1 = tight_muon_indices[0];
     UInt_t im2 = tight_muon_indices[1];
     UInt_t ie = tight_electron_indices[0];

     if (passLooseMuonSelectionV5((*muons)[im1],PV))
       flags = flags | Lep1LooseSelectionV5;

     if (passLooseMuonSelectionV5((*muons)[im2],PV))
       flags = flags | Lep2LooseSelectionV5;

     if (passTightMuonSelectionV1((*muons)[im1],PV)) 
       flags = flags | Lep1TightSelectionV1;

     if (passTightMuonSelectionV1((*muons)[im2],PV)) 
       flags = flags | Lep2TightSelectionV1;

     if (passTightMuonSelectionV2((*muons)[im1],PV)) 
       flags = flags | Lep1TightSelectionV2;

     if (passTightMuonSelectionV2((*muons)[im2],PV)) 
       flags = flags | Lep2TightSelectionV2;

     if (passTightMuonSelectionV3((*muons)[im1],PV) ) 
       flags = flags | Lep1TightSelectionV3;

     if (passTightMuonSelectionV3((*muons)[im2],PV) ) 
       flags = flags | Lep2TightSelectionV3;

     lep1 = (*muons)[im1].p4();
     lep1q = (*muons)[im1].charge();
     lep1id = (*muons)[im1].pdgId();
     lep1iso = muon_isolation((*muons)[im1],PV);

     lep2 = (*muons)[im2].p4();
     lep2q = (*muons)[im2].charge();
     lep2id = (*muons)[im2].pdgId();
     lep2iso = muon_isolation((*muons)[im2],PV);

     if (passTightElectronSelectionV1((*electrons)[ie], PV,rho))
       flags = flags | Lep3TightSelectionV1;
     if (passTightElectronSelectionV2((*electrons)[ie], PV,rho))
       flags = flags | Lep3TightSelectionV2;
     if (passLooseElectronSelectionV5((*electrons)[ie],PV,rho,rhoHLTElectronSelection))
	 flags = flags | Lep3LooseSelectionV5;
       
     lep3= (*electrons)[ie].p4();
     lep3q = (*electrons)[ie].charge();
     lep3id = (*electrons)[ie].pdgId();
     lep3iso = electron_isolation((*electrons)[ie],PV,rho);
     
   } else if (tight_muon_indices.size() >=1 && tight_electron_indices.size() >= 2) {

     UInt_t im = tight_muon_indices[0];
     UInt_t ie1 = tight_electron_indices[0];
     UInt_t ie2 = tight_electron_indices[1];

     if (passLooseMuonSelectionV5((*muons)[im],PV))
       flags = flags | Lep1LooseSelectionV5;

     if (passTightMuonSelectionV1((*muons)[im],PV)) 
       flags = flags | Lep1TightSelectionV1;

     if (passTightMuonSelectionV2((*muons)[im],PV)) 
       flags = flags | Lep1TightSelectionV2;

     if (passTightMuonSelectionV3((*muons)[im],PV) ) 
       flags = flags | Lep1TightSelectionV3;

     lep1 = (*muons)[im].p4();
     lep1q = (*muons)[im].charge();
     lep1id = (*muons)[im].pdgId();
     lep1iso = muon_isolation((*muons)[im],PV);

     if (passTightElectronSelectionV1((*electrons)[ie1], PV,rho))
       flags = flags | Lep2TightSelectionV1;
     if (passTightElectronSelectionV2((*electrons)[ie1], PV,rho))
       flags = flags | Lep2TightSelectionV2;
     if (passLooseElectronSelectionV5((*electrons)[ie1],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep2LooseSelectionV5;

     if (passTightElectronSelectionV1((*electrons)[ie2], PV,rho))
       flags = flags | Lep3TightSelectionV1;
     if (passTightElectronSelectionV2((*electrons)[ie2], PV,rho))
       flags = flags | Lep3TightSelectionV2;
     if (passLooseElectronSelectionV5((*electrons)[ie2],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep3LooseSelectionV5;

     lep2= (*electrons)[ie1].p4();
     lep2q = (*electrons)[ie1].charge();
     lep2id = (*electrons)[ie1].pdgId();
     lep2iso = electron_isolation((*electrons)[ie1],PV,rho);

     lep3= (*electrons)[ie2].p4();
     lep3q = (*electrons)[ie2].charge();
     lep3id = (*electrons)[ie2].pdgId();
     lep2iso = electron_isolation((*electrons)[ie2],PV,rho);

   } else if (tight_electron_indices.size() >= 3){

     UInt_t i1 = tight_electron_indices[0];
     UInt_t i2 = tight_electron_indices[1];
     UInt_t i3 = tight_electron_indices[2];


     lep1= (*electrons)[i1].p4();
     lep1q = (*electrons)[i1].charge();
     lep1id = (*electrons)[i1].pdgId();
     lep1iso = electron_isolation((*electrons)[i1],PV,rho);

     lep2= (*electrons)[i2].p4();
     lep2q = (*electrons)[i2].charge();
     lep2id = (*electrons)[i2].pdgId();
     lep2iso = electron_isolation((*electrons)[i2],PV,rho);

     lep3= (*electrons)[i3].p4();
     lep3q = (*electrons)[i3].charge();
     lep3id = (*electrons)[i3].pdgId();
     lep3iso = electron_isolation((*electrons)[i3],PV,rho);

     if (passTightElectronSelectionV1((*electrons)[i1],PV,rho))
	 flags = flags | Lep1TightSelectionV1;
     if (passTightElectronSelectionV1((*electrons)[i2],PV,rho))
	 flags = flags | Lep2TightSelectionV1;
     if (passTightElectronSelectionV1((*electrons)[i3],PV,rho))
	 flags = flags | Lep3TightSelectionV1;

     if (passTightElectronSelectionV2((*electrons)[i1],PV,rho))
	 flags = flags | Lep1TightSelectionV2;
     if (passTightElectronSelectionV2((*electrons)[i2],PV,rho))
	 flags = flags | Lep2TightSelectionV2;
     if (passTightElectronSelectionV2((*electrons)[i3],PV,rho))
	 flags = flags | Lep3TightSelectionV2;


     if (passLooseElectronSelectionV5((*electrons)[i1],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep1LooseSelectionV5;
     if (passLooseElectronSelectionV5((*electrons)[i2],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep2LooseSelectionV5;
     if (passLooseElectronSelectionV5((*electrons)[i3],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep3LooseSelectionV5;

   } else if(tight_muon_indices.size() >= 2 && veryloose_muon_indices.size() >= 1){

     UInt_t i1 = tight_muon_indices[0];
     UInt_t i2 = tight_muon_indices[1];
     UInt_t i3 = veryloose_muon_indices[0];

     if (passLooseMuonSelectionV5((*muons)[i1],PV) )
       flags = flags | Lep1LooseSelectionV5;

     if (passLooseMuonSelectionV5((*muons)[i2],PV) )
       flags = flags | Lep2LooseSelectionV5;

     if (passLooseMuonSelectionV5((*muons)[i3],PV) )
       flags = flags | Lep3LooseSelectionV5;

     if (passTightMuonSelectionV1((*muons)[i1],PV)) 
       flags = flags | Lep1TightSelectionV1;

     if (passTightMuonSelectionV1((*muons)[i2],PV)) 
       flags = flags | Lep2TightSelectionV1;

     if (passTightMuonSelectionV1((*muons)[i3],PV)) 
       flags = flags | Lep3TightSelectionV1;

     if (passTightMuonSelectionV2((*muons)[i1],PV)) 
       flags = flags | Lep1TightSelectionV2;

     if (passTightMuonSelectionV2((*muons)[i2],PV)) 
       flags = flags | Lep2TightSelectionV2;

     if (passTightMuonSelectionV2((*muons)[i3],PV)) 
       flags = flags | Lep3TightSelectionV2;

     if (passTightMuonSelectionV3((*muons)[i1],PV)) 
       flags = flags | Lep1TightSelectionV3;

     if (passTightMuonSelectionV3((*muons)[i2],PV)) 
       flags = flags | Lep2TightSelectionV3;

     if (passTightMuonSelectionV3((*muons)[i3],PV)) 
       flags = flags | Lep3TightSelectionV3;

     lep1 = (*muons)[i1].p4();
     lep1q = (*muons)[i1].charge();
     lep1id = (*muons)[i1].pdgId();
     lep1iso = muon_isolation((*muons)[i1],PV);

     lep2 = (*muons)[i2].p4();
     lep2q = (*muons)[i2].charge();
     lep2id = (*muons)[i2].pdgId();
     lep2iso = muon_isolation((*muons)[i2],PV);

     lep3 = (*muons)[i3].p4();
     lep3q = (*muons)[i3].charge();
     lep3id = (*muons)[i3].pdgId();
     lep3iso = muon_isolation((*muons)[i3],PV);
     
   } else if (tight_muon_indices.size() >=2 && veryloose_electron_indices.size() >= 1){

     UInt_t im1 = tight_muon_indices[0];
     UInt_t im2 = tight_muon_indices[1];
     UInt_t ie = veryloose_electron_indices[0];

     if (passLooseMuonSelectionV5((*muons)[im1],PV))
       flags = flags | Lep1LooseSelectionV5;

     if (passLooseMuonSelectionV5((*muons)[im2],PV))
       flags = flags | Lep2LooseSelectionV5;

     if (passTightMuonSelectionV1((*muons)[im1],PV)) 
       flags = flags | Lep1TightSelectionV1;

     if (passTightMuonSelectionV1((*muons)[im2],PV)) 
       flags = flags | Lep2TightSelectionV1;

     if (passTightMuonSelectionV2((*muons)[im1],PV)) 
       flags = flags | Lep1TightSelectionV2;

     if (passTightMuonSelectionV2((*muons)[im2],PV)) 
       flags = flags | Lep2TightSelectionV2;

     if (passTightMuonSelectionV3((*muons)[im1],PV) ) 
       flags = flags | Lep1TightSelectionV3;

     if (passTightMuonSelectionV3((*muons)[im2],PV) ) 
       flags = flags | Lep2TightSelectionV3;

     lep1 = (*muons)[im1].p4();
     lep1q = (*muons)[im1].charge();
     lep1id = (*muons)[im1].pdgId();
     lep1iso = muon_isolation((*muons)[im1],PV);

     lep2 = (*muons)[im2].p4();
     lep2q = (*muons)[im2].charge();
     lep2id = (*muons)[im2].pdgId();
     lep2iso = muon_isolation((*muons)[im2],PV);

     if (passTightElectronSelectionV1((*electrons)[ie], PV,rho))
       flags = flags | Lep3TightSelectionV1;
     if (passTightElectronSelectionV2((*electrons)[ie], PV,rho))
       flags = flags | Lep3TightSelectionV2;
     if (passLooseElectronSelectionV5((*electrons)[ie],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep3LooseSelectionV5;
       
     lep3= (*electrons)[ie].p4();
     lep3q = (*electrons)[ie].charge();
     lep3id = (*electrons)[ie].pdgId();
     lep3iso = electron_isolation((*electrons)[ie],PV,rho);


   } else if (veryloose_muon_indices.size() >=1 && tight_electron_indices.size() >= 2) {

     UInt_t im = veryloose_muon_indices[0];
     UInt_t ie1 = tight_electron_indices[0];
     UInt_t ie2 = tight_electron_indices[1];

     if (passLooseMuonSelectionV5((*muons)[im],PV))
       flags = flags | Lep1LooseSelectionV5;

     if (passTightMuonSelectionV1((*muons)[im],PV)) 
       flags = flags | Lep1TightSelectionV1;

     if (passTightMuonSelectionV2((*muons)[im],PV)) 
       flags = flags | Lep1TightSelectionV2;

     if (passTightMuonSelectionV3((*muons)[im],PV) ) 
       flags = flags | Lep1TightSelectionV3;

     lep1 = (*muons)[im].p4();
     lep1q = (*muons)[im].charge();
     lep1id = (*muons)[im].pdgId();
     lep1iso = muon_isolation((*muons)[im],PV);

     if (passTightElectronSelectionV1((*electrons)[ie1], PV,rho))
       flags = flags | Lep2TightSelectionV1;
     if (passTightElectronSelectionV2((*electrons)[ie1], PV,rho))
       flags = flags | Lep2TightSelectionV2;
     if (passLooseElectronSelectionV5((*electrons)[ie1],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep2LooseSelectionV5;


     if (passTightElectronSelectionV1((*electrons)[ie2], PV,rho))
       flags = flags | Lep3TightSelectionV1;
     if (passTightElectronSelectionV2((*electrons)[ie2], PV,rho))
       flags = flags | Lep3TightSelectionV2;
     if (passLooseElectronSelectionV5((*electrons)[ie2],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep3LooseSelectionV5;

       
     lep2= (*electrons)[ie1].p4();
     lep2q = (*electrons)[ie1].charge();
     lep2id = (*electrons)[ie1].pdgId();
     lep2iso = electron_isolation((*electrons)[ie1],PV,rho);

     lep3= (*electrons)[ie2].p4();
     lep3q = (*electrons)[ie2].charge();
     lep3id = (*electrons)[ie2].pdgId();
     lep3iso = electron_isolation((*electrons)[ie2],PV,rho);

   } else if (tight_electron_indices.size() >= 2 && veryloose_electron_indices.size() >= 1){

     UInt_t i1 = tight_electron_indices[0];
     UInt_t i2 = tight_electron_indices[1];
     UInt_t i3 = veryloose_electron_indices[0];

     lep1= (*electrons)[i1].p4();
     lep1q = (*electrons)[i1].charge();
     lep1id = (*electrons)[i1].pdgId();
     lep1iso = electron_isolation((*electrons)[i1],PV,rho);

     lep2= (*electrons)[i2].p4();
     lep2q = (*electrons)[i2].charge();
     lep2id = (*electrons)[i2].pdgId();
     lep2iso = electron_isolation((*electrons)[i2],PV,rho);

     lep3= (*electrons)[i3].p4();
     lep3q = (*electrons)[i3].charge();
     lep3id = (*electrons)[i3].pdgId();
     lep3iso = electron_isolation((*electrons)[i3],PV,rho);

     if (passTightElectronSelectionV1((*electrons)[i1],PV,rho))
	 flags = flags | Lep1TightSelectionV1;
     if (passTightElectronSelectionV1((*electrons)[i2],PV,rho))
	 flags = flags | Lep2TightSelectionV1;
     if (passTightElectronSelectionV1((*electrons)[i3],PV,rho))
	 flags = flags | Lep3TightSelectionV1;

     if (passTightElectronSelectionV2((*electrons)[i1],PV,rho))
	 flags = flags | Lep1TightSelectionV2;
     if (passTightElectronSelectionV2((*electrons)[i2],PV,rho))
	 flags = flags | Lep2TightSelectionV2;
     if (passTightElectronSelectionV2((*electrons)[i3],PV,rho))
	 flags = flags | Lep3TightSelectionV2;


     if (passLooseElectronSelectionV5((*electrons)[i1],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep1LooseSelectionV5;
     if (passLooseElectronSelectionV5((*electrons)[i2],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep2LooseSelectionV5;
     if (passLooseElectronSelectionV5((*electrons)[i3],PV,rho, rhoHLTElectronSelection))
	 flags = flags | Lep3LooseSelectionV5;

 } else 
     return;

   //   lhe_and_gen_object.analyze(iEvent,lep1,lep2);

   std::vector<const pat::Jet *> cleaned_jets;

   for (const pat::Jet &j : *jets) {

     if ( reco::deltaR(j.p4(),lep1) < 0.5 || reco::deltaR(j.p4(),lep2) < 0.5 || reco::deltaR(j.p4(),lep3) < 0.5 ){
       continue;

     }
     cleaned_jets.push_back(&j);
   }

   if (cleaned_jets.size() < 2) 
     return;
   else if (cleaned_jets[1]->pt() < 20)
     return;

   jet1=cleaned_jets[0]->p4();
   jet2=cleaned_jets[1]->p4();   
   jet1btag = std::max(0.f,cleaned_jets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
   jet2btag = std::max(0.f,cleaned_jets[1]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
   jet1pujetid = cleaned_jets[0]->userFloat("pileupJetId:fullDiscriminant");
   jet2pujetid = cleaned_jets[1]->userFloat("pileupJetId:fullDiscriminant");

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
threeljj_ntuple_maker::beginJob()
{

  edm::Service<TFileService> fs;

  if (isMC_){
    n_events_run_over= fs->make<TH1F>("n_events_run_over","n_events_run_over",1,0,1);
    n_weighted_events_run_over= fs->make<TH1F>("n_weighted_events_run_over","n_weighted_events_run_over",1,0,1);
  }

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
  tree->Branch("jet1",&jet1);
  tree->Branch("jet2",&jet2);
  tree->Branch("jet1pujetid",&jet1pujetid);
  tree->Branch("jet2pujetid",&jet2pujetid);
  tree->Branch("jet1btag",&jet1btag);
  tree->Branch("jet2btag",&jet2btag);
  tree->Branch("lep1",&lep1);
  tree->Branch("lep2",&lep2);
  tree->Branch("lep3",&lep3);
  tree->Branch("nvtx",&nvtx);
  tree->Branch("lep1q",&lep1q);
  tree->Branch("lep1id",&lep1id);
  tree->Branch("lep2q",&lep2q);
  tree->Branch("lep2id",&lep2id);
  tree->Branch("lep3q",&lep3q);
  tree->Branch("lep3id",&lep3id);
  tree->Branch("lep1iso",&lep1iso);
  tree->Branch("lep2iso",&lep2iso);
  tree->Branch("lep3iso",&lep3iso);

  tree->Branch("maxbtagevent",&maxbtagevent);

  if (isMC_)
    tree->Branch("gen_weight",&gen_weight);

  //  lhe_and_gen_object.defineBranches(tree);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
threeljj_ntuple_maker::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------


void 
threeljj_ntuple_maker::beginRun(edm::Run const& iRun, edm::EventSetup const&)
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
threeljj_ntuple_maker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(threeljj_ntuple_maker);
