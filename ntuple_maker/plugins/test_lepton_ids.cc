// -*- C++ -*-
//
// Package:    ntuple_maker
// Class:      test_lepton_ids
// 
/**\class test_lepton_ids test_lepton_ids.cc ntuple_maker/ntuple_maker/plugins/test_lepton_ids.cc

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

#include "ntuple_maker/ntuple_maker/interface/enum_definition.h"
#include "ntuple_maker/ntuple_maker/interface/lepton_ids.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


//
// class declaration
//

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class test_lepton_ids : public edm::EDAnalyzer {
public:

  
      explicit test_lepton_ids(const edm::ParameterSet&);
      ~test_lepton_ids();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

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
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;

  TH1F * n_events_run_over;
  UInt_t flags;
  UInt_t event;
  UInt_t run;
  UInt_t lumi;
  UInt_t nvtx;
  TTree * tree;
  Float_t jetpt;
  Float_t jet1pujetid;
  Float_t jet1btag;
  Float_t jet2pujetid;
  Float_t jet2btag;
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
  LorentzVector lep1_nearestparton_4mom;
  LorentzVector lep2_nearestparton_4mom;
  Int_t lep1_nearestparton_pdgid;
  Int_t lep2_nearestparton_pdgid;
  Int_t lep1_matching_real_gen_lepton_pdgid;
  Int_t lep2_matching_real_gen_lepton_pdgid;
  Int_t lep1_matching_real_gen_lepton_q;
  Int_t lep2_matching_real_gen_lepton_q;
  Float_t lep1_nearest_gen_electron_dr;
  Float_t lep2_nearest_gen_electron_dr;
  Float_t lep1_nearest_gen_muon_dr;
  Float_t lep2_nearest_gen_muon_dr;

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
test_lepton_ids::test_lepton_ids(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedgenparticles"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedgenparticles")))
{
  //now do what ever initialization is needed

}


test_lepton_ids::~test_lepton_ids()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
test_lepton_ids::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  n_events_run_over->Fill(0.5);

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
   
  flags = 0;

   using namespace edm;

   run=iEvent.eventAuxiliary().run(); 
   lumi=iEvent.eventAuxiliary().luminosityBlock();
   event=iEvent.eventAuxiliary().event();

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if (vertices->empty()) return; // skip the event if no PV found
   const reco::Vertex &PV = vertices->front();

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);

   for(UInt_t i = 0; i < muons->size(); i++){

     std::cout << passTightMuonId((*muons)[i], PV) << " " << (*muons)[i].isTightMuon(PV) << std::endl;

     assert(passTightMuonId((*muons)[i], PV) == (*muons)[i].isTightMuon(PV));
     
   }


}


// ------------ method called once each job just before starting event loop  ------------
void 
test_lepton_ids::beginJob()
{

  edm::Service<TFileService> fs;

  n_events_run_over= fs->make<TH1F>("n_events_run_over","n_events_run_over",1,0,1);

  tree = fs->make<TTree>( "events"  , "events");

  tree->Branch("flags",&flags);

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
  tree->Branch("lep1",&lep1);
  tree->Branch("lep2",&lep2);
  tree->Branch("nvtx",&nvtx);
  tree->Branch("lep1q",&lep1q);
  tree->Branch("lep1id",&lep1id);
  tree->Branch("lep2q",&lep2q);
  tree->Branch("lep2id",&lep2id);
  tree->Branch("lep1_nearestparton_4mom",&lep1_nearestparton_4mom);
  tree->Branch("lep1_nearestparton_pdgid",&lep1_nearestparton_pdgid);
  tree->Branch("lep2_nearestparton_4mom",&lep2_nearestparton_4mom);
  tree->Branch("lep2_nearestparton_pdgid",&lep2_nearestparton_pdgid);
  tree->Branch("lep1_matching_real_gen_lepton_pdgid",&lep1_matching_real_gen_lepton_pdgid);
  tree->Branch("lep2_matching_real_gen_lepton_pdgid",&lep2_matching_real_gen_lepton_pdgid);
  tree->Branch("lep1_matching_real_gen_lepton_q",&lep1_matching_real_gen_lepton_q);
  tree->Branch("lep2_matching_real_gen_lepton_q",&lep2_matching_real_gen_lepton_q);
  tree->Branch("lep1_nearest_gen_electron_dr",&lep1_nearest_gen_electron_dr);
  tree->Branch("lep2_nearest_gen_electron_dr",&lep2_nearest_gen_electron_dr);
  tree->Branch("lep1_nearest_gen_muon_dr",&lep1_nearest_gen_muon_dr);
  tree->Branch("lep2_nearest_gen_muon_dr",&lep2_nearest_gen_muon_dr);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
test_lepton_ids::endJob() 
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
test_lepton_ids::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(test_lepton_ids);
