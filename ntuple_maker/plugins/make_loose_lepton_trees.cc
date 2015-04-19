// -*- C++ -*-
//
// Package:    make_loose_lepton_trees
// Class:      make_loose_lepton_trees
// 
/**\class make_loose_lepton_trees make_loose_lepton_trees.cc ntuple_maker/ntuple_maker/plugins/make_loose_lepton_trees.cc

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

#include "DataFormats/Math/interface/deltaR.h"

#include "ntuple_maker/ntuple_maker/interface/fr_enum_definition.h"
#include "ntuple_maker/ntuple_maker/interface/lepton_ids.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


//
// class declaration
//



class make_loose_lepton_trees : public edm::EDAnalyzer {
public:
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
  
      explicit make_loose_lepton_trees(const edm::ParameterSet&);
      ~make_loose_lepton_trees();

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
  UInt_t njets;
  UInt_t nmuons;
  UInt_t n_loose_muons;
  UInt_t n_loose_electrons;
  UInt_t nelectrons;
  TTree * muon_tree;
  TTree * electron_tree;
  Float_t jetpt;
  Float_t drnearestgenmuon;
  Float_t drnearestgenelectron;

  //the pt of the highest pt jet that is at least delta R = 1 away from the lepton
  Float_t ptjetaway;

  Float_t maxjetbtag;
  Float_t metpt;
  Float_t metphi;
  Float_t metsumet;
  Float_t metgenmetpt;
  Float_t metptshiftup;
  Float_t metptshiftdown;
  LorentzVector jet1;
  LorentzVector jet2;
  LorentzVector nearestparton_4mom;
  Int_t nearestparton_pdgid;
  LorentzVector electron_4mom;
  LorentzVector muon_4mom;
  LorentzVector lep2;
  Int_t lep1id;
  Int_t lep2id;
  Int_t lep1q;
  Int_t lep2q;
  Bool_t pass_tight_muon_id;
  Bool_t pass_tight_electron_id;


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
make_loose_lepton_trees::make_loose_lepton_trees(const edm::ParameterSet& iConfig):
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


make_loose_lepton_trees::~make_loose_lepton_trees()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
make_loose_lepton_trees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  n_events_run_over->Fill(0.5);

  std::vector<UInt_t> loose_muon_indices;
  std::vector<UInt_t> loose_electron_indices;

  flags = 0;

   using namespace edm;

   run=iEvent.eventAuxiliary().run(); 
   lumi=iEvent.eventAuxiliary().luminosityBlock();
   event=iEvent.eventAuxiliary().event();

   //if (event != 6612489)
   //  return;

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if (vertices->empty()) return; // skip the event if no PV found
   const reco::Vertex &PV = vertices->front();

   nvtx=vertices->size();

   edm::Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);

   nelectrons = electrons->size();

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);

   nmuons = muons->size();

   n_loose_muons = 0;
   n_loose_electrons = 0;

   pass_tight_muon_id = kFALSE;
   pass_tight_electron_id = kFALSE;

   for(UInt_t i = 0; i < muons->size(); i++){

     if( (*muons)[i].pt() < 10)
       continue;

     Float_t relative_isolation = ((*muons)[i].pfIsolationR04().sumChargedHadronPt+ std::max(0.0,(*muons)[i].pfIsolationR04().sumNeutralHadronEt + (*muons)[i].pfIsolationR04().sumPhotonEt - 0.5 * (*muons)[i].pfIsolationR04().sumPUPt))/(*muons)[i].pt();

     if ( passLooseMuonId((*muons)[i],PV) && relative_isolation < 1.0)
       flags = flags | LepLooseSelectionV1;

     if ( !(flags & LepLooseSelectionV1)  )
       continue;

     n_loose_muons++;

     muon_4mom = (*muons)[i].p4();

     if (passTightMuonId((*muons)[i],PV) && relative_isolation < 0.4) {
       flags = flags | LepTightSelectionV1;
     }
     if (passTightMuonId((*muons)[i],PV) && relative_isolation < 0.2) {
       flags = flags | LepTightSelectionV2;
     }
     if (passTightMuonId((*muons)[i],PV) && relative_isolation < 0.12) {
       flags = flags | LepTightSelectionV3;
     }

   }

   for(UInt_t i = 0; i < electrons->size(); i++){

     if( (*electrons)[i].pt() < 10)
       continue;

     if (! passLooseElectronId((*electrons)[i], PV))
       continue;

     flags = flags | LepLooseSelectionV1;

     n_loose_electrons++;

     electron_4mom = (*electrons)[i].p4();
     
     if (passTightElectronId((*electrons)[i], PV))
       flags = flags | LepTightSelectionV1;
     
   }
   
   edm::Handle<pat::JetCollection> jets;
   iEvent.getByToken(jetToken_, jets);

   njets = jets->size();

   Float_t maxbtag = 0;

   for (const pat::Jet &j : *jets) {

       if (maxbtag < j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags")) 
	 maxbtag = j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
     }

   maxjetbtag = maxbtag;

   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();

   metpt = met.pt();
   metphi = met.phi();
   metsumet = met.sumEt();
   metgenmetpt = met.genMET()->pt();
   metptshiftup = met.shiftedPt(pat::MET::JetEnUp);
   metptshiftdown = met.shiftedPt(pat::MET::JetEnDown);

   Handle<edm::View<reco::GenParticle> > pruned;
   iEvent.getByToken(prunedGenToken_,pruned);
   // Packed particles are all the status 1, so usable to remake jets
   // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
   Handle<edm::View<pat::PackedGenParticle> > packed;
   iEvent.getByToken(packedGenToken_,packed);
   
   nearestparton_pdgid=0;

   if (n_loose_electrons + n_loose_muons > 1)
     std::cout  << "warning: n_loose_leptons > 1, not saving the event" << std::endl;
   else if (n_loose_electrons == 1){

     drnearestgenmuon = std::numeric_limits<Float_t>::max();
     drnearestgenelectron = std::numeric_limits<Float_t>::max();

     for(size_t j=0; j<packed->size();j++){
       if (abs((*packed)[j].pdgId()) == 11 && (*packed)[j].status() == 1){
	 if (reco::deltaR((*packed)[j].p4(),electron_4mom) < drnearestgenelectron)
	   drnearestgenelectron = reco::deltaR((*packed)[j].p4(),electron_4mom);
	 flags = flags | GenElectronInTheEvent;
	 }
       if (abs((*packed)[j].pdgId()) == 13 && (*packed)[j].status() == 1){
	 if (reco::deltaR((*packed)[j].p4(),electron_4mom) < drnearestgenmuon)
	   drnearestgenmuon = reco::deltaR((*packed)[j].p4(),electron_4mom);
	 flags = flags | GenMuonInTheEvent;
	 }
     }

     Float_t minpartondR = std::numeric_limits<Float_t>::max();
     
     for(size_t j=0; j<pruned->size();j++){
       if ( abs((*pruned)[j].pdgId()) == 1 || abs((*pruned)[j].pdgId()) == 2 || abs((*pruned)[j].pdgId()) == 3 || abs((*pruned)[j].pdgId()) == 4 || abs((*pruned)[j].pdgId()) == 5 || abs((*pruned)[j].pdgId()) == 21 ){
	 
	 //status 23 means outgoing: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
	 if (   (     (*pruned)[j].status() == 23 || (*pruned)[j].status() == 71 )  && reco::deltaR((*pruned)[j].p4(),electron_4mom) < minpartondR){
	   minpartondR = reco::deltaR((*pruned)[j].p4(),electron_4mom);
	   nearestparton_pdgid=(*pruned)[j].pdgId();
	   nearestparton_4mom=(*pruned)[j].p4();
	   
	 }
	 
       //std::cout << (*pruned)[j].pdgId()  << " " << (*pruned)[j].status() << std::endl;
       }
     }

     Float_t maxptjetaway = -1;

     for (const pat::Jet &j : *jets) {

       if (j.pt() > maxptjetaway && reco::deltaR(j,electron_4mom) > 1)
	 maxptjetaway = j.pt();
     }

     ptjetaway=maxptjetaway;
     electron_tree->Fill();
   }
   else if (n_loose_muons == 1){

     drnearestgenelectron = std::numeric_limits<Float_t>::max();
     drnearestgenmuon = std::numeric_limits<Float_t>::max();

     for(size_t j=0; j<packed->size();j++){
       if (abs((*packed)[j].pdgId()) == 11 && (*packed)[j].status() == 1){
	 if (reco::deltaR((*packed)[j].p4(),muon_4mom) < drnearestgenelectron)
	   drnearestgenelectron = reco::deltaR((*packed)[j].p4(),muon_4mom);
	 flags = flags | GenElectronInTheEvent;
	 }
       if (abs((*packed)[j].pdgId()) == 13 && (*packed)[j].status() == 1){
	 if (reco::deltaR((*packed)[j].p4(),muon_4mom) < drnearestgenmuon)
	   drnearestgenmuon = reco::deltaR((*packed)[j].p4(),muon_4mom);
	 flags = flags | GenMuonInTheEvent;
	 }
     }

     Float_t minpartondR = std::numeric_limits<Float_t>::max();
     
     for(size_t j=0; j<pruned->size();j++){
       
       if ( abs((*pruned)[j].pdgId()) == 1 || abs((*pruned)[j].pdgId()) == 2 || abs((*pruned)[j].pdgId()) == 3 || abs((*pruned)[j].pdgId()) == 4 || abs((*pruned)[j].pdgId()) == 5 || abs((*pruned)[j].pdgId()) == 21 ){

	 //status 23 means outgoing: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
	 if ( (  (*pruned)[j].status() == 23 || (*pruned)[j].status() == 71 )  && reco::deltaR((*pruned)[j].p4(),muon_4mom) < minpartondR){
	   
	   minpartondR = reco::deltaR((*pruned)[j].p4(),muon_4mom);
	   nearestparton_pdgid=(*pruned)[j].pdgId();
	   nearestparton_4mom=(*pruned)[j].p4();
	 }
	 
       //std::cout << (*pruned)[j].pdgId()  << " " << (*pruned)[j].status() << std::endl;
       }
     }

     Float_t maxptjetaway = -1;

     for (const pat::Jet &j : *jets) {

       if (j.pt() > maxptjetaway && reco::deltaR(j,muon_4mom) > 1)
	 maxptjetaway = j.pt();
     }

     ptjetaway=maxptjetaway;

     muon_tree->Fill();

   }
       


}


// ------------ method called once each job just before starting event loop  ------------
void 
make_loose_lepton_trees::beginJob()
{

  edm::Service<TFileService> fs;

  n_events_run_over= fs->make<TH1F>("n_events_run_over","n_events_run_over",1,0,1);

  muon_tree = fs->make<TTree>( "loose_muons"  , "loose_muons");
  electron_tree = fs->make<TTree>( "loose_electrons"  , "loose_electrons");

  muon_tree->Branch("event",&event);
  muon_tree->Branch("lumi",&lumi);
  muon_tree->Branch("run",&run);
  muon_tree->Branch("maxjetbtag",&maxjetbtag);
  muon_tree->Branch("nearestparton_4mom",&nearestparton_4mom);
  muon_tree->Branch("nearestparton_pdgid",&nearestparton_pdgid);
  muon_tree->Branch("muon_4mom",&muon_4mom);
  muon_tree->Branch("pass_tight_muon_id",&pass_tight_muon_id);
  muon_tree->Branch("ptjetaway",&ptjetaway);
  muon_tree->Branch("metpt",&metpt);
  muon_tree->Branch("flags",&flags);
  muon_tree->Branch("drnearestgenmuon",&drnearestgenmuon);
  muon_tree->Branch("drnearestgenelectron",&drnearestgenelectron);

  electron_tree->Branch("event",&event);
  electron_tree->Branch("lumi",&lumi);
  electron_tree->Branch("run",&run);
  electron_tree->Branch("maxjetbtag",&maxjetbtag);
  electron_tree->Branch("nearestparton_4mom",&nearestparton_4mom);
  electron_tree->Branch("nearestparton_pdgid",&nearestparton_pdgid);
  electron_tree->Branch("electron_4mom",&electron_4mom);
  electron_tree->Branch("pass_tight_electron_id",&pass_tight_electron_id);
  electron_tree->Branch("ptjetaway",&ptjetaway);
  electron_tree->Branch("metpt",&metpt);
  electron_tree->Branch("flags",&flags);
  electron_tree->Branch("drnearestgenelectron",&drnearestgenelectron);
  electron_tree->Branch("drnearestgenmuon",&drnearestgenmuon);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
make_loose_lepton_trees::endJob() 
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
make_loose_lepton_trees::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(make_loose_lepton_trees);
