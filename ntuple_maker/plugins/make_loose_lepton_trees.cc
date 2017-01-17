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

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "ntuple_maker/ntuple_maker/interface/triggers.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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
  edm::EDGetTokenT< edm::TriggerResults > triggerResultsToken_;
  edm::EDGetTokenT< pat::TriggerObjectStandAloneCollection > triggerObjectToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEvtToken_;
  edm::EDGetTokenT<double> rhoToken_;
  //  edm::EDGetTokenT<LHEEventProduct> lheEvtToken_;
  edm::EDGetTokenT<double> rhoHLTElectronSelectionToken_;

  TH1F * n_events_run_over;
  TH1F * n_weighted_events_run_over;

  UInt_t flags;
  UInt_t event;
  UInt_t run;
  Float_t iso;

  UInt_t lumi;
  UInt_t nvtx;
  UInt_t njets;
  UInt_t nmuons;
  UInt_t n_veryloose_muons;
  UInt_t n_veryloose_electrons;
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

  Float_t gen_weight;
  //  Float_t lhe_weight_orig;

  Bool_t isMC_;
  std::string lepton_flavor_;
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
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedgenparticles"))),
  //  triggerResultsToken_(consumes< edm::TriggerResults >(edm::InputTag("TriggerResults","","HLT"))),
  triggerResultsToken_(consumes< edm::TriggerResults >(edm::InputTag("TriggerResults","","HLT2"))),
  triggerObjectToken_( consumes< pat::TriggerObjectStandAloneCollection >(edm::InputTag("selectedPatTrigger"))),
  genEvtToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genevent"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  rhoHLTElectronSelectionToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoHLTElectronSelection"))),
  //  lheEvtToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheevent"))),
  isMC_(iConfig.getUntrackedParameter<bool>("isMC")),
  lepton_flavor_(iConfig.getUntrackedParameter<std::string>("lepton_flavor")),
  which_triggers_(iConfig.getUntrackedParameter<std::string>("which_triggers"))

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


  if(isMC_){

    n_events_run_over->Fill(0.5);

    edm::Handle<GenEventInfoProduct> hGenEvt;
    iEvent.getByToken(genEvtToken_,hGenEvt);

    //edm::Handle<LHEEventProduct> hLheEvt;                                                                                                                                 //iEvent.getByToken(lheEvtToken_,hLheEvt);   
  
    if (hGenEvt->weight() > 0)
      n_weighted_events_run_over->Fill(0.5,1);
    else
      n_weighted_events_run_over->Fill(0.5,-1);


    gen_weight = hGenEvt->weight();

    //lhe_weight_orig = hLheEvt->originalXWGTUP();
    
  }

  /*


  for(unsigned int i = 0; i < hLheEvt->hepeup().IDUP.size(); i++){
    std::cout << hLheEvt->hepeup().IDUP[i] << std::endl;

    LorentzVector v;

    v.SetPxPyPzE(hLheEvt->hepeup().PUP.at(i)[0],hLheEvt->hepeup().PUP.at(i)[1],hLheEvt->hepeup().PUP.at(i)[2],hLheEvt->hepeup().PUP.at(i)[3]);

      std::cout << "v.eta() = " << v.eta() << std::endl;
      std::cout << "v.pt() = " << v.pt() << std::endl;

  }

  */

  //std::cout << "hLheEvt->hepeup().IDUP.size() = " << hLheEvt->hepeup().IDUP.size() << std::endl;

  edm::Handle<double> rhoHandle;

  iEvent.getByToken(rhoToken_,rhoHandle);

  float rho    =  *rhoHandle;

  edm::Handle<double> rhoHLTElectronSelectionHandle;

  iEvent.getByToken(rhoHLTElectronSelectionToken_,rhoHLTElectronSelectionHandle);

  float rhoHLTElectronSelection  =  *rhoHLTElectronSelectionHandle;

  edm::Handle< edm::TriggerResults> triggerResultsHandle;

  iEvent.getByToken(triggerResultsToken_,triggerResultsHandle);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResultsHandle);

  if (! trigger_fired(names,triggerResultsHandle,which_triggers_))
    return;


  edm::Handle<pat::TriggerObjectStandAloneCollection > triggerObjectHandle;

  iEvent.getByToken(triggerObjectToken_,triggerObjectHandle);

  for (pat::TriggerObjectStandAlone obj : *triggerObjectHandle) {

    obj.unpackPathNames(names);

    //std::cout << "obj.pathNames(false).size() = " << obj.pathNames(false).size() << std::endl;

    for (unsigned h = 0, n =  obj.pathNames(false).size(); h < n; ++h) {
      bool isBoth = obj.hasPathName( obj.pathNames(false)[h], true, true ); 
      bool isL3   = obj.hasPathName( obj.pathNames(false)[h], false, true ); 
      bool isLF   = obj.hasPathName( obj.pathNames(false)[h], true, false ); 
      bool isNone = obj.hasPathName( obj.pathNames(false)[h], false, false ); 

      if (isNone && !isBoth && !isL3 && !isLF) 
      	continue;

      if (obj.filterIds()[h] == trigger::TriggerMuon){

	if( obj.pathNames(false)[h].find("HLT_Mu8_v1") != std::string::npos ){

	  std::cout << "obj.pt() = " << obj.pt() << std::endl;

	 }

      }

    }

  }



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

   n_veryloose_muons = 0;
   n_veryloose_electrons = 0;

   UInt_t ie = 0;
   UInt_t im = 0;

   for(UInt_t i = 0; i < muons->size(); i++){

     if( (*muons)[i].pt() < 10)
       continue;

     if ( ! passVeryLooseMuonSelection((*muons)[i],PV)  )
       continue;

     n_veryloose_muons++;

     if ( passLooseMuonSelectionV1((*muons)[i],PV))
       flags = flags | LepLooseSelectionV1;
     if ( passLooseMuonSelectionV2((*muons)[i],PV))
       flags = flags | LepLooseSelectionV2;
     if ( passLooseMuonSelectionV3((*muons)[i],PV))
       flags = flags | LepLooseSelectionV3;
     if ( passLooseMuonSelectionV4((*muons)[i],PV))
       flags = flags | LepLooseSelectionV4;
     if ( passLooseMuonSelectionV5((*muons)[i],PV))
       flags = flags | LepLooseSelectionV5;


     muon_4mom = (*muons)[i].p4();
     im = i;


     if (passTightMuonSelectionV1((*muons)[i],PV)) {
       flags = flags | LepTightSelectionV1;
     }
     if (passTightMuonSelectionV2((*muons)[i],PV)) {
       flags = flags | LepTightSelectionV2;
     }
     if (passTightMuonSelectionV3((*muons)[i],PV)) {
       flags = flags | LepTightSelectionV3;
     }

   }

   for(UInt_t i = 0; i < electrons->size(); i++){

     if( (*electrons)[i].pt() < 10)
       continue;

     if (! passVeryLooseElectronSelection((*electrons)[i], PV,rho,rhoHLTElectronSelection))
       continue;

     n_veryloose_electrons++;

     if (passLooseElectronSelectionV1((*electrons)[i], PV,rho,rhoHLTElectronSelection))
       flags = flags | LepLooseSelectionV1;
     if (passLooseElectronSelectionV2((*electrons)[i], PV,rho,rhoHLTElectronSelection))
       flags = flags | LepLooseSelectionV2;
     if (passLooseElectronSelectionV3((*electrons)[i], PV,rho,rhoHLTElectronSelection))
       flags = flags | LepLooseSelectionV3;
     if (passLooseElectronSelectionV4((*electrons)[i], PV,rho,rhoHLTElectronSelection))
       flags = flags | LepLooseSelectionV4;
     if (passLooseElectronSelectionV5((*electrons)[i], PV,rho,rhoHLTElectronSelection))
       flags = flags | LepLooseSelectionV5;

     electron_4mom = (*electrons)[i].p4();
     ie = i;

     if (passTightElectronSelectionV1((*electrons)[i], PV,rho))
       flags = flags | LepTightSelectionV1;
     if (passTightElectronSelectionV2((*electrons)[i], PV,rho))
       flags = flags | LepTightSelectionV2;
     if (passTightElectronSelectionV3((*electrons)[i], PV,rho))
       flags = flags | LepTightSelectionV3;
     if (passTightElectronSelectionV4((*electrons)[i], PV,rho))
       flags = flags | LepTightSelectionV4;
     if (passTightElectronSelectionV5((*electrons)[i], PV,rho))
       flags = flags | LepTightSelectionV5;
     
   }


   //std::cout << "n_veryloose_electrons = " << n_veryloose_electrons << std::endl;
   //std::cout << "n_veryloose_muon = " << n_veryloose_muons << std::endl;
   
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

   Handle<edm::View<reco::GenParticle> > pruned;
   Handle<edm::View<pat::PackedGenParticle> > packed;

   if (isMC_){

     iEvent.getByToken(prunedGenToken_,pruned);
     // Packed particles are all the status 1, so usable to remake jets
     // The navigation from status 1 to pruned is possible (the other direction should be made by hand)

     iEvent.getByToken(packedGenToken_,packed);
     
     nearestparton_pdgid=0;
   }

   if (lepton_flavor_ == "electron" && n_veryloose_electrons == 1){

     if (isMC_){

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
       
     }

     Float_t maxptjetaway = -1;

     for (const pat::Jet &j : *jets) {

       if (j.pt() > maxptjetaway && reco::deltaR(j,electron_4mom) > 1)
	 maxptjetaway = j.pt();
     }

     ptjetaway=maxptjetaway;
     
     iso = electron_isolation((*electrons)[ie], PV,rho);

     electron_tree->Fill();
   }
   else if (lepton_flavor_ == "muon" && n_veryloose_muons == 1){

     if(isMC_){

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
       
     }

     Float_t maxptjetaway = -1;

     for (const pat::Jet &j : *jets) {

       if (j.pt() > maxptjetaway && reco::deltaR(j,muon_4mom) > 1)
	 maxptjetaway = j.pt();
     }

     ptjetaway=maxptjetaway;

     iso = muon_isolation((*muons)[im], PV);

     muon_tree->Fill();

   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
make_loose_lepton_trees::beginJob()
{

  edm::Service<TFileService> fs;

  if (isMC_){

    n_events_run_over= fs->make<TH1F>("n_events_run_over","n_events_run_over",1,0,1);
    n_weighted_events_run_over= fs->make<TH1F>("n_weighted_events_run_over","n_weighted_events_run_over",1,0,1);

  }

  assert(lepton_flavor_ == "muon" || lepton_flavor_ == "electron");

  if (lepton_flavor_ == "muon"){

    muon_tree = fs->make<TTree>( "loose_muons"  , "loose_muons");

    muon_tree->Branch("event",&event);
    muon_tree->Branch("lumi",&lumi);
    muon_tree->Branch("run",&run);
    muon_tree->Branch("maxjetbtag",&maxjetbtag);
    muon_tree->Branch("nearestparton_4mom",&nearestparton_4mom);
    muon_tree->Branch("nearestparton_pdgid",&nearestparton_pdgid);
    muon_tree->Branch("muon_4mom",&muon_4mom);
    muon_tree->Branch("ptjetaway",&ptjetaway);
    muon_tree->Branch("metpt",&metpt);
    muon_tree->Branch("iso",&iso);
    muon_tree->Branch("metphi",&metphi);
    muon_tree->Branch("flags",&flags);
    muon_tree->Branch("drnearestgenmuon",&drnearestgenmuon);
    muon_tree->Branch("drnearestgenelectron",&drnearestgenelectron);

    if (isMC_) {

      muon_tree->Branch("gen_weight",&gen_weight);
      //      muon_tree->Branch("lhe_weight_orig",&lhe_weight_orig);

    }

  } else {

    electron_tree = fs->make<TTree>( "loose_electrons"  , "loose_electrons");
    
    electron_tree->Branch("event",&event);
    electron_tree->Branch("lumi",&lumi);
    electron_tree->Branch("run",&run);
    electron_tree->Branch("maxjetbtag",&maxjetbtag);
    electron_tree->Branch("nearestparton_4mom",&nearestparton_4mom);
    electron_tree->Branch("nearestparton_pdgid",&nearestparton_pdgid);
    electron_tree->Branch("electron_4mom",&electron_4mom);
    electron_tree->Branch("ptjetaway",&ptjetaway);
    electron_tree->Branch("metpt",&metpt);
    electron_tree->Branch("metphi",&metphi);
    electron_tree->Branch("iso",&iso);
    electron_tree->Branch("flags",&flags);
    electron_tree->Branch("drnearestgenelectron",&drnearestgenelectron);
    electron_tree->Branch("drnearestgenmuon",&drnearestgenmuon);

    if (isMC_) {

      electron_tree->Branch("gen_weight",&gen_weight);
      //      electron_tree->Branch("lhe_weight_orig",&lhe_weight_orig);
    
    }

  }

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
