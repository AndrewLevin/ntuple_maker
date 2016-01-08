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
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  edm::EDGetTokenT<double> rhoToken_;

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
  //electronToken_(mayConsume<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedgenparticles"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedgenparticles"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho")))

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

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
   
   using namespace edm;

   edm::Handle<double> rhoHandle;

   iEvent.getByToken(rhoToken_,rhoHandle);

   float rho    =  *rhoHandle;


   //run=iEvent.eventAuxiliary().run(); 
   //lumi=iEvent.eventAuxiliary().luminosityBlock();
   //event=iEvent.eventAuxiliary().event();

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if (vertices->empty()) return; // skip the event if no PV found
   const reco::Vertex &PV = vertices->front();

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);

   //for(UInt_t i = 0; i < muons->size(); i++){
   //
   //     std::cout << passTightMuonIdV1((*muons)[i], PV) << " " << (*muons)[i].isTightMuon(PV) << std::endl;
   //
   //     assert(passTightMuonIdV1((*muons)[i], PV) == (*muons)[i].isTightMuon(PV));
   //     
   //   }

   //edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
   //edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 

   //iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
   //iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);
   
   edm::Handle<pat::ElectronCollection> electrons;


   //edm::Handle<edm::View<reco::GsfElectron> > electrons;
   iEvent.getByToken(electronToken_, electrons);

   for(UInt_t i = 0; i < electrons->size(); i++){

     assert( (*electrons)[i].isElectronIDAvailable("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") );
     assert( (*electrons)[i].isElectronIDAvailable("cutBasedElectronID-Spring15-25ns-V1-standalone-tight") );
     //std::cout << (*electrons)[i].electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight") << std::endl;
     //std::cout << (*electrons)[i].electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight") << std::endl;

     //std::cout << "andrew debug 0" << std::endl;

     std::cout << "(*electrons)[i].pt() = " << (*electrons)[i].pt() << std::endl;

     std::cout << passTightElectronSelectionV1((*electrons)[i], PV,rho) << " " << ((*electrons)[i].chargeInfo().isGsfCtfScPixConsistent && (*electrons)[i].electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium")) << std::endl;

     assert(passTightElectronSelectionV1((*electrons)[i], PV,rho) == ((*electrons)[i].chargeInfo().isGsfCtfScPixConsistent && (*electrons)[i].electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium")) );

     //std::cout << "andrew debug 1" << std::endl;

     //assert(passTightElectronSelectionV1((*electrons)[i], PV,rho) == (*electrons)[i].electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium"));

     //bool isPassMedium = (*medium_id_decisions)[el];
     //bool isPassTight  = (*tight_id_decisions)[el];

     
     //std::cout << "isPassMedium = " << isPassMedium << std::endl;
     //std::cout << "isPassTight = " << isPassTight << std::endl;

     
     //std::cout << (*electrons)[i].userFloat("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose") << std::endl;

   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
test_lepton_ids::beginJob()
{

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
