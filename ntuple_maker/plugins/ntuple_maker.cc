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
#include <limits>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"

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

#include "ntuple_maker/ntuple_maker/interface/lhe_and_gen.h"
#include "ntuple_maker/ntuple_maker/interface/triggers.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "TRandom3.h"

//
// class declaration
//

const Float_t z_mass = 91.18800;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;


inline void fill_electron_selection_flags(const pat::Electron & el, const reco::Vertex &PV, const double &rho, const double &rhoHLTElectronSelection, UInt_t &flags, const int & flag_index) {

  if (flag_index == 1) {

    if (passTightElectronSelectionV1(el, PV,rho))
      flags = flags | Lep1TightSelectionV1;
    if (passTightElectronSelectionV2(el, PV,rho))
      flags = flags | Lep1TightSelectionV2;
    if (passTightElectronSelectionV3(el, PV,rho))
      flags = flags | Lep1TightSelectionV3;
    if (passTightElectronSelectionV4(el, PV,rho))
      flags = flags | Lep1TightSelectionV4;
    if (passTightElectronSelectionV5(el, PV,rho))
      flags = flags | Lep1TightSelectionV5;
    if (passLooseElectronSelectionV1(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep1LooseSelectionV1;
    if (passLooseElectronSelectionV2(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep1LooseSelectionV2;
    if (passLooseElectronSelectionV3(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep1LooseSelectionV3;
    if (passLooseElectronSelectionV4(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep1LooseSelectionV4;
    if (passLooseElectronSelectionV5(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep1LooseSelectionV5;

  }
  else if (flag_index == 2) {
    
    if (passTightElectronSelectionV1(el, PV,rho))
      flags = flags | Lep2TightSelectionV1;
    if (passTightElectronSelectionV2(el, PV,rho))
      flags = flags | Lep2TightSelectionV2;
    if (passTightElectronSelectionV3(el, PV,rho))
      flags = flags | Lep2TightSelectionV3;
    if (passTightElectronSelectionV4(el, PV,rho))
      flags = flags | Lep2TightSelectionV4;
    if (passTightElectronSelectionV5(el, PV,rho))
      flags = flags | Lep2TightSelectionV5;
    if (passLooseElectronSelectionV1(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep2LooseSelectionV1;
    if (passLooseElectronSelectionV2(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep2LooseSelectionV2;
    if (passLooseElectronSelectionV3(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep2LooseSelectionV3;
    if (passLooseElectronSelectionV4(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep2LooseSelectionV4;
    if (passLooseElectronSelectionV5(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep2LooseSelectionV5;    

  }
  else {
    assert(flag_index == 3);
    
    if (passTightElectronSelectionV1(el, PV,rho))
      flags = flags | Lep3TightSelectionV1;
    if (passTightElectronSelectionV2(el, PV,rho))
      flags = flags | Lep3TightSelectionV2;
    if (passTightElectronSelectionV3(el, PV,rho))
      flags = flags | Lep3TightSelectionV3;
    if (passTightElectronSelectionV4(el, PV,rho))
      flags = flags | Lep3TightSelectionV4;
    if (passTightElectronSelectionV5(el, PV,rho))
      flags = flags | Lep3TightSelectionV5;
    if (passLooseElectronSelectionV1(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep3LooseSelectionV1;
    if (passLooseElectronSelectionV2(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep3LooseSelectionV2;
    if (passLooseElectronSelectionV3(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep3LooseSelectionV3;
    if (passLooseElectronSelectionV4(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep3LooseSelectionV4;
    if (passLooseElectronSelectionV5(el,PV,rho,rhoHLTElectronSelection))
      flags = flags | Lep3LooseSelectionV5;    

 }

}

inline void fill_muon_selection_flags(const pat::Muon & muon, const reco::Vertex &PV,UInt_t &flags,const int & flag_index) {

  if (flag_index == 1){

    if (passLooseMuonSelectionV1(muon,PV) )
      flags = flags | Lep1LooseSelectionV1;
    
    if (passLooseMuonSelectionV2(muon,PV) )
      flags = flags | Lep1LooseSelectionV2;
    
    if (passLooseMuonSelectionV3(muon,PV) )
      flags = flags | Lep1LooseSelectionV3;
    
    if (passLooseMuonSelectionV4(muon,PV) )
      flags = flags | Lep1LooseSelectionV4;
    
    if (passLooseMuonSelectionV5(muon,PV) )
      flags = flags | Lep1LooseSelectionV5;
    
    if (passTightMuonSelectionV1(muon,PV))
      flags = flags | Lep1TightSelectionV1;
    
    if (passTightMuonSelectionV2(muon,PV)) 
      flags = flags | Lep1TightSelectionV2;
    
    if (passTightMuonSelectionV3(muon,PV)) 
      flags = flags | Lep1TightSelectionV3;
  }
  else if (flag_index == 2) {
    
    if (passLooseMuonSelectionV1(muon,PV) )
      flags = flags | Lep2LooseSelectionV1;
    
    if (passLooseMuonSelectionV2(muon,PV))
      flags = flags | Lep2LooseSelectionV2;
    
    if (passLooseMuonSelectionV3(muon,PV) )
      flags = flags | Lep2LooseSelectionV3;
    
    if (passLooseMuonSelectionV4(muon,PV) )
      flags = flags | Lep2LooseSelectionV4;
    
    if (passLooseMuonSelectionV5(muon,PV) )
      flags = flags | Lep2LooseSelectionV5;
    
    if (passTightMuonSelectionV1(muon,PV))
      flags = flags | Lep2TightSelectionV1;
    
    if (passTightMuonSelectionV2(muon,PV)) 
      flags = flags | Lep2TightSelectionV2;
    
    if (passTightMuonSelectionV3(muon,PV)) 
      flags = flags | Lep2TightSelectionV3;
    
  }   else  {
    
    assert(flag_index == 3);
    
    if (passLooseMuonSelectionV1(muon,PV) )
      flags = flags | Lep3LooseSelectionV1;
    
    if (passLooseMuonSelectionV2(muon,PV) )
      flags = flags | Lep3LooseSelectionV2;
	
    if (passLooseMuonSelectionV3(muon,PV) )
      flags = flags | Lep3LooseSelectionV3;
    
    if (passLooseMuonSelectionV4(muon,PV) )
      flags = flags | Lep3LooseSelectionV4;
    
    if (passLooseMuonSelectionV5(muon,PV) )
      flags = flags | Lep3LooseSelectionV5;
    
    if (passTightMuonSelectionV1(muon,PV))
      flags = flags | Lep3TightSelectionV1;
    
    if (passTightMuonSelectionV2(muon,PV)) 
      flags = flags | Lep3TightSelectionV2;
  
    if (passTightMuonSelectionV3(muon,PV)) 
      flags = flags | Lep3TightSelectionV3;
    
  }

}

std::pair<int,int> get_two_highest_pt_jet_indices(const std::vector<float> &corrected_jet_pts) {

  assert(corrected_jet_pts.size() >= 2);

  int index1 = -1;
  int index2 = -1;
  

  for (unsigned int i = 0; i < corrected_jet_pts.size(); i++){

    if (index1 == -1){
      index1 = i;
      continue;
    }

    if (index2 == -1){
      index2 = i;
      continue;
    }
    
    if (corrected_jet_pts[index1] < corrected_jet_pts[i]){
      index1 = i;
      continue;
    }

    if (corrected_jet_pts[index2] < corrected_jet_pts[i]){
      index2 = i;
      continue;
    }
    

  }

  assert(index1 != -1 && index2 != -1);

  return std::pair<int,int>(index1,index2);

}

class ntuple_maker : public edm::EDAnalyzer {
public:

  
      explicit ntuple_maker(const edm::ParameterSet&);
      ~ntuple_maker();

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
  edm::InputTag lheRunInfoLabel_;
  edm::EDGetTokenT< pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::EDGetTokenT<LHEEventProduct> lheEvtToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEvtToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<double> rhoHLTElectronSelectionToken_;


  TH1F * n_events_run_over;
  TH1F * n_weighted_events_run_over;
  UInt_t lepton_selection_flags;
  UInt_t flags;
  ULong64_t event;
  UInt_t run;
  UInt_t lumi;
  UInt_t nvtx;
  TTree * tree;

  Float_t maxbtagevent;

  Float_t jetpt;
  Float_t jet1pujetid;
  bool jet1loosejetid;
  Float_t jet1btag;
  Float_t jet2pujetid;
  Float_t jet2btag;
  bool jet2loosejetid;
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
  Int_t lep1id;
  Int_t lep2id;
  Int_t lep3id;
  Float_t lep1iso;
  Float_t lep2iso;
  Float_t lep3iso;
  Int_t lep1q;
  Int_t lep2q;
  Int_t lep3q;
  lhe_and_gen lhe_and_gen_object; //separate the part that runs over the generator and lhe information

  Int_t n_pu_interactions;

  bool syscalcinfo_;
  bool mgreweightinfo_;

  bool apply_trigger_;


  bool isMC_;

  JetCorrectionUncertainty *jecUnc_;

  bool isJecUncSet_ = false;

  std::string jes_;
  std::string jer_;

  bool lhe_lepton_info_;

private:

  TRandom3*rnd_{0};

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
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  triggerResultsToken_(consumes< edm::TriggerResults >(edm::InputTag("TriggerResults","",iConfig.getUntrackedParameter<std::string>("trigger_results_process")))),
  triggerObjectToken_( consumes< pat::TriggerObjectStandAloneCollection >(edm::InputTag("selectedPatTrigger"))),
  lheRunInfoLabel_(iConfig.getParameter<edm::InputTag>("lheruninfo")),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  pileupSummaryToken_(consumes <std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileup_summary"))),
  lheEvtToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheevent"))),
  genEvtToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genevent"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  rhoHLTElectronSelectionToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoHLTElectronSelection"))),
  syscalcinfo_(iConfig.getUntrackedParameter<bool>("syscalcinfo")),
  mgreweightinfo_(iConfig.getUntrackedParameter<bool>("mgreweightinfo")),
  apply_trigger_(iConfig.getUntrackedParameter<bool>("apply_trigger")),
  isMC_(iConfig.getUntrackedParameter<bool>("isMC")),
  jes_(iConfig.getUntrackedParameter<std::string>("jes")),
  jer_(iConfig.getUntrackedParameter<std::string>("jer")),
  lhe_lepton_info_(iConfig.getUntrackedParameter<bool>("lheleptoninfo"))
{
  //now do what ever initialization is needed

  rnd_=new TRandom3();
  rnd_->SetSeed((unsigned)time(NULL));

  lhe_and_gen_object.isMC_ = isMC_;
  lhe_and_gen_object.prunedGenToken_ = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedgenparticles"));
  lhe_and_gen_object.packedGenToken_ = consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedgenparticles"));
  lhe_and_gen_object.lheEvtToken_ = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheevent"));
  lhe_and_gen_object.genEvtToken_ = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genevent"));
  lhe_and_gen_object.syscalcinfo_ = syscalcinfo_;
  lhe_and_gen_object.mgreweightinfo_ = mgreweightinfo_;
  lhe_and_gen_object.lhe_lepton_info_ = lhe_lepton_info_;
  lhe_and_gen_object.lheRunInfoLabel_ = iConfig.getParameter<edm::InputTag>("lheruninfo");

  consumes< LHERunInfoProduct, edm::InRun > (iConfig.getParameter<edm::InputTag>("lheruninfo"));
  
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

  //if (iEvent.eventAuxiliary().luminosityBlock() != 11 || iEvent.eventAuxiliary().event() != 2092)
  //  return;

  lepton_selection_flags = 0;
  flags = 0;
  lep3id = 0;
  lep3q = 0;

  if (!isJecUncSet_ && isMC_){

    std::string payload = "AK4PFchs";
    
    edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    
    iSetup.get<JetCorrectionsRecord>().get(payload ,JetCorParColl); 
    
    JetCorrectorParameters const& par = (*JetCorParColl)["Uncertainty"];
    
    jecUnc_ = new JetCorrectionUncertainty(par);

    isJecUncSet_ = true;

  }

  edm::Handle<double> rhoHandle;

  iEvent.getByToken(rhoToken_,rhoHandle);

  float rho    =  *rhoHandle;

  edm::Handle<double> rhoHLTElectronSelectionHandle;

  iEvent.getByToken(rhoHLTElectronSelectionToken_,rhoHLTElectronSelectionHandle);

  float rhoHLTElectronSelection  =  *rhoHLTElectronSelectionHandle;

  
  if (isMC_){
    n_events_run_over->Fill(0.5);

    //edm::Handle<LHEEventProduct> hLheEvt;
    //iEvent.getByToken(lheEvtToken_,hLheEvt);

    edm::Handle<GenEventInfoProduct> hGenEvt;
    iEvent.getByToken(genEvtToken_,hGenEvt);
    
    //std::cout << "hGenEvt->weight() = " << hGenEvt->weight() << std::endl;
    //std::cout << "hLheEvt->originalXWGTUP() = " << hLheEvt->originalXWGTUP() << std::endl;

    //assert(hGenEvt->weight() == 1 || hGenEvt->weight() == -1);

    //assert(hGenEvt->weight() * hLheEvt->originalXWGTUP() > 0);

    if (hGenEvt->weight() > 0)
      n_weighted_events_run_over->Fill(0.5,1);
    else
      n_weighted_events_run_over->Fill(0.5,-1);
  }


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

  edm::Handle< edm::TriggerResults> triggerResultsHandle;

  iEvent.getByToken(triggerResultsToken_,triggerResultsHandle);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResultsHandle);

  if (trigger_fired(names,triggerResultsHandle,"doublelepton"))
    flags = flags | PassTriggerV1;

  if (trigger_fired(names,triggerResultsHandle,"soup"))
    flags = flags | PassTriggerV2;
    
  if (trigger_fired(names,triggerResultsHandle,"doublemu"))
    flags = flags | PassTriggerV3;

  if (trigger_fired(names,triggerResultsHandle,"doubleeg"))
    flags = flags | PassTriggerV4;

  if (trigger_fired(names,triggerResultsHandle,"muoneg"))
    flags = flags | PassTriggerV5;

  if (trigger_fired(names,triggerResultsHandle,"doublemudz"))
    flags = flags | PassTriggerV6;

  //    if (apply_trigger_ && ! trigger_fired(names,triggerResultsHandle,"soup"))
  //        return;
  
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

   //   std::cout << "muons->size() = " << muons->size() << std::endl;
   //   std::cout << "electrons->size() = " << electrons->size() << std::endl;


   for(UInt_t i = 0; i < electrons->size(); i++){
     for(UInt_t j = 0; j < electrons->size(); j++){

       if (i == j)
	 continue;

       if ((*electrons)[i].pt() < 10 || (*electrons)[j].pt() < 10 || abs((*electrons)[i].eta()) > 2.5 || abs((*electrons)[j].eta()) > 2.5)
	 continue;

       if ((*electrons)[j].charge() != (*electrons)[i].charge() && abs(((*electrons)[i].p4() + (*electrons)[j].p4()).mass() - z_mass)  < 15){

	 flags = flags | VetoV8;

	 if( passTightElectronSelectionV4((*electrons)[i],PV,rho) && passVeryLooseElectronSelection((*electrons)[j],PV,rho,rhoHLTElectronSelection)){

	   flags = flags | VetoV9;
	 
	 }

       }
     }
   }

   for(UInt_t i = 0; i < electrons->size(); i++){

     if( (*electrons)[i].pt() < 10 || abs((*electrons)[i].eta()) > 2.5) 
       continue;

     /*

     bool should_be_cleaned = false;

     for (UInt_t j = 0; j < muons->size(); j++){
       if (reco::deltaR((*electrons)[i],(*muons)[j]) < 0.05 && passTightMuonIdV2((*muons)[j],PV) && (*muons)[j].pt() > 20 && abs((*muons)[j].eta()) < 2.4 ){

	 should_be_cleaned = true;
       }
     }

     if (should_be_cleaned) continue;

     */

     if (passTightElectronSelectionV4((*electrons)[i],PV,rho))
       tight_electron_indices.push_back(i);
     else if (passVeryLooseElectronSelection((*electrons)[i],PV,rho,rhoHLTElectronSelection))
       veryloose_electron_indices.push_back(i);



   }

   for(UInt_t i = 0; i < muons->size(); i++){
     for(UInt_t j = i+1; j < muons->size(); j++){

       if ((*muons)[i].pt() < 10 || (*muons)[j].pt() < 10)
	 continue;

       if ((*muons)[j].charge() != (*muons)[i].charge() && abs(((*muons)[i].p4() + (*muons)[j].p4()).mass() - z_mass)  < 15){

	 flags = flags | VetoV6;
	 
       }
     }
   }

   for(UInt_t i = 0; i < muons->size(); i++){

     if ((*muons)[i].pt() < 10 || abs((*muons)[i].eta()) > 2.4)
       continue;

     //if (!passTightMuonSelectionV1((*muons)[i],PV) )
     //  continue;

     if (passTightMuonSelectionV1( (*muons)[i],PV ) )
       tight_muon_indices.push_back(i);
     else if (passVeryLooseMuonSelection( (*muons)[i],PV ) )
       veryloose_muon_indices.push_back(i);
     

    }

   //         std::cout << "tight_muon_indices.size() = " << tight_muon_indices.size() << std::endl;
   //         std::cout << "tight_electron_indices.size() = " << tight_electron_indices.size() << std::endl;
   //         std::cout << "veryloose_muon_indices.size() = " << veryloose_muon_indices.size() << std::endl;
   //         std::cout << "veryloose_electron_indices.size() = " << veryloose_electron_indices.size() << std::endl;

   if(tight_muon_indices.size() >= 2){

     UInt_t i1 = tight_muon_indices[0];
     UInt_t i2 = tight_muon_indices[1];

     fill_muon_selection_flags((*muons)[i1],PV,lepton_selection_flags,1);
     fill_muon_selection_flags((*muons)[i2],PV,lepton_selection_flags,2);

     lep1 = (*muons)[i1].p4();
     lep1q = (*muons)[i1].charge();
     lep1id = (*muons)[i1].pdgId();
     lep1iso = muon_isolation((*muons)[i1],PV);

     lep2 = (*muons)[i2].p4();
     lep2q = (*muons)[i2].charge();
     lep2id = (*muons)[i2].pdgId();
     lep2iso = muon_isolation((*muons)[i2],PV);

     if (tight_muon_indices.size() >= 3){

       UInt_t i3 = tight_muon_indices[2];

       fill_muon_selection_flags((*muons)[i3],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[i3].p4();
       lep3q = (*muons)[i3].charge();
       lep3id = (*muons)[i3].pdgId();
       lep3iso = muon_isolation((*muons)[i3],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 4 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == i1 || i == i2 || i == i3)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (tight_electron_indices.size() >= 1) {

       UInt_t ie = tight_electron_indices[0];

       fill_electron_selection_flags((*electrons)[ie],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie].p4();
       lep3q = (*electrons)[ie].charge();
       lep3id = (*electrons)[ie].pdgId();
       lep3iso = electron_isolation((*electrons)[ie],PV,rho);
       
       //fourth lepton veto
       if (tight_muon_indices.size() >= 3 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == i1 || i == i2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_muon_indices.size() >= 1){

       UInt_t i3 = veryloose_muon_indices[0];

       fill_muon_selection_flags((*muons)[i3],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[i3].p4();
       lep3q = (*muons)[i3].charge();
       lep3id = (*muons)[i3].pdgId();
       lep3iso = muon_isolation((*muons)[i3],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 3 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == i1 || i == i2 || i == i3)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_electron_indices.size() >= 1){

       UInt_t ie = veryloose_electron_indices[0];

       fill_electron_selection_flags((*electrons)[ie],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie].p4();
       lep3q = (*electrons)[ie].charge();
       lep3id = (*electrons)[ie].pdgId();
       lep3iso = electron_isolation((*electrons)[ie],PV,rho);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 3 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == i1 || i == i2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     }

     for(UInt_t i = 0; i < muons->size(); i++){

       if (i == i1 || i == i2)
	 continue;

       if ((passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4) && !(flags & VetoV5)  )
	 flags = flags | VetoV5;

       else if (passLooseMuonSelectionV1((*muons)[i],PV)  && (*muons)[i].pt() > 10 && abs((*muons)[i].eta()) < 2.4 )
	 flags = flags | VetoV4;

       
     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if ( passWLLJJVetoElectronId( (*electrons)[i], PV, rho, rhoHLTElectronSelection ) && (*electrons)[i].pt() > 10 && abs((*electrons)[i].eta()) < 2.5) {

	 flags = flags | VetoV2;
       }
     }


   }
   else if (tight_muon_indices.size() >=1 && tight_electron_indices.size() >= 1){

     UInt_t im = tight_muon_indices[0];
     UInt_t ie = tight_electron_indices[0];

     fill_muon_selection_flags((*muons)[im],PV,lepton_selection_flags,1);
     fill_electron_selection_flags((*electrons)[ie],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,2);

     lep1 = (*muons)[im].p4();
     lep1q = (*muons)[im].charge();
     lep1id = (*muons)[im].pdgId();
     lep1iso = muon_isolation((*muons)[im],PV);

     lep2= (*electrons)[ie].p4();
     lep2q = (*electrons)[ie].charge();
     lep2id = (*electrons)[ie].pdgId();
     lep2iso = electron_isolation((*electrons)[ie],PV,rho);

     //     std::cout << "lep2.pt() = " << lep2.pt() << std::endl;

     if (tight_muon_indices.size() >= 2){

       UInt_t im2  = tight_muon_indices[1];

       fill_muon_selection_flags((*muons)[im2],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im2].p4();
       lep3q = (*muons)[im2].charge();
       lep3id = (*muons)[im2].pdgId();
       lep3iso = muon_isolation((*muons)[im2],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 3 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im || i == im2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (tight_electron_indices.size() >= 2) {

       UInt_t ie2 = tight_electron_indices[1];

       fill_electron_selection_flags((*electrons)[ie2],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie2].p4();
       lep3q = (*electrons)[ie2].charge();
       lep3id = (*electrons)[ie2].pdgId();
       lep3iso = electron_isolation((*electrons)[ie2],PV,rho);
       
       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 3)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_muon_indices.size() >= 1){

       UInt_t im2 = veryloose_muon_indices[0];

       fill_muon_selection_flags((*muons)[im2],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im2].p4();
       lep3q = (*muons)[im2].charge();
       lep3id = (*muons)[im2].pdgId();
       lep3iso = muon_isolation((*muons)[im2],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im || i == im2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
 
     } else if (veryloose_electron_indices.size() >= 1){

       UInt_t ie2 = veryloose_electron_indices[0];

       fill_electron_selection_flags((*electrons)[ie2],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie2].p4();
       lep3q = (*electrons)[ie2].charge();
       lep3id = (*electrons)[ie2].pdgId();
       lep3iso = electron_isolation((*electrons)[ie2],PV,rho);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     }


     for(UInt_t i = 0; i < muons->size(); i++){

       if (i == im)
	 continue;

       if ((passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4) && !(flags & VetoV5))
	 flags = flags | VetoV5;

       else if (passLooseMuonSelectionV1((*muons)[i],PV)  && (*muons)[i].pt() > 10 && abs((*muons)[i].eta()) < 2.4 )
	 flags = flags | VetoV4;

     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if (i == ie)
	 continue;

       if ( passWLLJJVetoElectronId( (*electrons)[i], PV, rho, rhoHLTElectronSelection)  && (*electrons)[i].pt() > 10 && abs((*electrons)[i].eta()) < 2.5) {

	 flags = flags | VetoV2;
       }

     }

     
   } else if (tight_electron_indices.size() >= 2){

     UInt_t i1 = tight_electron_indices[0];
     UInt_t i2 = tight_electron_indices[1];

     lep1= (*electrons)[i1].p4();
     lep1q = (*electrons)[i1].charge();
     lep1id = (*electrons)[i1].pdgId();
     lep1iso = electron_isolation((*electrons)[i1],PV,rho);

     lep2= (*electrons)[i2].p4();
     lep2q = (*electrons)[i2].charge();
     lep2id = (*electrons)[i2].pdgId();
     lep2iso = electron_isolation((*electrons)[i2],PV,rho);

     fill_electron_selection_flags((*electrons)[i1],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,1);
     fill_electron_selection_flags((*electrons)[i2],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,2);

     if (tight_muon_indices.size() >= 1){

       UInt_t im = tight_muon_indices[0];

       fill_muon_selection_flags((*muons)[im],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im].p4();
       lep3q = (*muons)[im].charge();
       lep3id = (*muons)[im].pdgId();
       lep3iso = muon_isolation((*muons)[im],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 3)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (tight_electron_indices.size() >= 3) {

       UInt_t i3 = tight_electron_indices[2];

       fill_electron_selection_flags((*electrons)[i3],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[i3].p4();
       lep3q = (*electrons)[i3].charge();
       lep3id = (*electrons)[i3].pdgId();
       lep3iso = electron_isolation((*electrons)[i3],PV,rho);
       
       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 4)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_muon_indices.size() >= 1){

       UInt_t im = veryloose_muon_indices[0];

       fill_muon_selection_flags((*muons)[im],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im].p4();
       lep3q = (*muons)[im].charge();
       lep3id = (*muons)[im].pdgId();
       lep3iso = muon_isolation((*muons)[im],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 3)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
 
     } else if (veryloose_electron_indices.size() >= 1){

       UInt_t i3 = veryloose_electron_indices[0];

       fill_electron_selection_flags((*electrons)[i3],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[i3].p4();
       lep3q = (*electrons)[i3].charge();
       lep3id = (*electrons)[i3].pdgId();
       lep3iso = electron_isolation((*electrons)[i3],PV,rho);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 3)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     }
    
     //     std::cout << "lepton_selection_flags & Lep2LooseSelectionV5 = " << bool(lepton_selection_flags & Lep2LooseSelectionV5) << std::endl;
     //     std::cout << "lepton_selection_flags & Lep2LooseSelectionV3 = " << bool(lepton_selection_flags & Lep2LooseSelectionV3) << std::endl;

     for(UInt_t i = 0; i < muons->size(); i++){

       if ((passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4) && !(flags & VetoV5))
	 flags = flags | VetoV5;

       else if (passLooseMuonSelectionV1((*muons)[i],PV)  && (*muons)[i].pt() > 10 && abs((*muons)[i].eta()) < 2.4 )
	 flags = flags | VetoV4;

     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if (i == i1 || i == i2)
	 continue;

       if ( passWLLJJVetoElectronId( (*electrons)[i],PV, rho, rhoHLTElectronSelection )  && (*electrons)[i].pt() > 10 && abs((*electrons)[i].eta()) < 2.5) 
	 flags = flags | VetoV2;

     }


   } else if (tight_muon_indices.size() >= 1 && veryloose_muon_indices.size() >= 1){

     UInt_t i1 = tight_muon_indices[0];
     UInt_t i2 = veryloose_muon_indices[0];

     fill_muon_selection_flags((*muons)[i1],PV,lepton_selection_flags,1);
     fill_muon_selection_flags((*muons)[i2],PV,lepton_selection_flags,2);

     lep1 = (*muons)[i1].p4();
     lep1q = (*muons)[i1].charge();
     lep1id = (*muons)[i1].pdgId();
     lep1iso = muon_isolation((*muons)[i1],PV);

     lep2 = (*muons)[i2].p4();
     lep2q = (*muons)[i2].charge();
     lep2id = (*muons)[i2].pdgId();
     lep2iso = muon_isolation((*muons)[i2],PV);

     if (tight_muon_indices.size() >= 2){

       UInt_t i3 = tight_muon_indices[1];

       fill_muon_selection_flags((*muons)[i3],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[i3].p4();
       lep3q = (*muons)[i3].charge();
       lep3id = (*muons)[i3].pdgId();
       lep3iso = muon_isolation((*muons)[i3],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 3 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == i1 || i == i2 || i == i3)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (tight_electron_indices.size() >= 1) {

       UInt_t ie = tight_electron_indices[0];

       fill_electron_selection_flags((*electrons)[ie],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie].p4();
       lep3q = (*electrons)[ie].charge();
       lep3id = (*electrons)[ie].pdgId();
       lep3iso = electron_isolation((*electrons)[ie],PV,rho);
       
       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == i1 || i == i2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_muon_indices.size() >= 2){

       UInt_t i3 = veryloose_muon_indices[1];

       fill_muon_selection_flags((*muons)[i3],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[i3].p4();
       lep3q = (*muons)[i3].charge();
       lep3id = (*muons)[i3].pdgId();
       lep3iso = muon_isolation((*muons)[i3],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 3 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == i1 || i == i2 || i == i3)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_electron_indices.size() >= 1){

       UInt_t ie = veryloose_electron_indices[0];

       fill_electron_selection_flags((*electrons)[ie],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie].p4();
       lep3q = (*electrons)[ie].charge();
       lep3id = (*electrons)[ie].pdgId();
       lep3iso = electron_isolation((*electrons)[ie],PV,rho);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == i1 || i == i2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     }

     
     for(UInt_t i = 0; i < muons->size(); i++){

       if (i == i1 || i == i2)
	 continue;

       if ((passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4) && !(flags & VetoV5))
	 flags = flags | VetoV5;

       else if (passLooseMuonSelectionV1((*muons)[i],PV)  && (*muons)[i].pt() > 10 && abs((*muons)[i].eta()) < 2.4 )
	 flags = flags | VetoV4;

     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if ( passWLLJJVetoElectronId( (*electrons)[i], PV,  rho, rhoHLTElectronSelection )  && (*electrons)[i].pt() > 10 && abs((*electrons)[i].eta()) < 2.5) 
	 flags = flags | VetoV2;

     }

   } else if (tight_muon_indices.size() >= 1 && veryloose_electron_indices.size() >= 1){

     UInt_t im = tight_muon_indices[0];
     UInt_t ie = veryloose_electron_indices[0];

     fill_muon_selection_flags((*muons)[im],PV,lepton_selection_flags,1);
     fill_electron_selection_flags((*electrons)[ie],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,2);

     lep1 = (*muons)[im].p4();
     lep1q = (*muons)[im].charge();
     lep1id = (*muons)[im].pdgId();
     lep1iso = muon_isolation((*muons)[im],PV);

     lep2= (*electrons)[ie].p4();
     lep2q = (*electrons)[ie].charge();
     lep2id = (*electrons)[ie].pdgId();
     lep2iso = electron_isolation((*electrons)[ie],PV,rho);

     if (tight_muon_indices.size() >= 2){

       UInt_t im2 = tight_muon_indices[1];

       fill_muon_selection_flags((*muons)[im2],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im2].p4();
       lep3q = (*muons)[im2].charge();
       lep3id = (*muons)[im2].pdgId();
       lep3iso = muon_isolation((*muons)[im2],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 3 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im || i == im2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (tight_electron_indices.size() >= 1) {

       UInt_t ie2 = tight_electron_indices[0];

       fill_electron_selection_flags((*electrons)[ie2],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie2].p4();
       lep3q = (*electrons)[ie2].charge();
       lep3id = (*electrons)[ie2].pdgId();
       lep3iso = electron_isolation((*electrons)[ie2],PV,rho);
       
       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_muon_indices.size() >= 1){

       UInt_t im2 = veryloose_muon_indices[0];

       fill_muon_selection_flags((*muons)[im2],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im2].p4();
       lep3q = (*muons)[im2].charge();
       lep3id = (*muons)[im2].pdgId();
       lep3iso = muon_isolation((*muons)[im2],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im || i == im2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_electron_indices.size() >= 2){

       UInt_t ie2 = veryloose_electron_indices[1];

       fill_electron_selection_flags((*electrons)[ie2],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie2].p4();
       lep3q = (*electrons)[ie2].charge();
       lep3id = (*electrons)[ie2].pdgId();
       lep3iso = electron_isolation((*electrons)[ie2],PV,rho);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 3 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     }


     for(UInt_t i = 0; i < muons->size(); i++){

       if (i == im)
	 continue;

       if ((passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4) && !(flags & VetoV5))
	 flags = flags | VetoV5;

       else if (passLooseMuonSelectionV1((*muons)[i],PV)  && (*muons)[i].pt() > 10 && abs((*muons)[i].eta()) < 2.4 )
	 flags = flags | VetoV4;

     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if (i == ie)
	 continue;

       if ( passWLLJJVetoElectronId( (*electrons)[i], PV ,rho,  rhoHLTElectronSelection )  && (*electrons)[i].pt() > 10 && abs((*electrons)[i].eta()) < 2.5) 
	 flags = flags | VetoV2;


     }

     

   } else if (tight_electron_indices.size() >= 1 && veryloose_muon_indices.size() >= 1){


     UInt_t im = veryloose_muon_indices[0];
     UInt_t ie = tight_electron_indices[0];

     fill_muon_selection_flags((*muons)[im],PV,lepton_selection_flags,1);
     fill_electron_selection_flags((*electrons)[ie],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,2);

     lep1 = (*muons)[im].p4();
     lep1q = (*muons)[im].charge();
     lep1id = (*muons)[im].pdgId();
     lep1iso = muon_isolation((*muons)[im],PV);
       
     lep2= (*electrons)[ie].p4();
     lep2q = (*electrons)[ie].charge();
     lep2id = (*electrons)[ie].pdgId();
     lep2iso = electron_isolation((*electrons)[ie],PV,rho);

     if (tight_muon_indices.size() >= 1){

       UInt_t im2 = tight_muon_indices[0];

       fill_muon_selection_flags((*muons)[im2],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im2].p4();
       lep3q = (*muons)[im2].charge();
       lep3id = (*muons)[im2].pdgId();
       lep3iso = muon_isolation((*muons)[im2],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im || i == im2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (tight_electron_indices.size() >= 2) {

       UInt_t ie2 = tight_electron_indices[1];

       fill_electron_selection_flags((*electrons)[ie2],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie2].p4();
       lep3q = (*electrons)[ie2].charge();
       lep3id = (*electrons)[ie2].pdgId();
       lep3iso = electron_isolation((*electrons)[ie2],PV,rho);
       
       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 3)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_muon_indices.size() >= 2){

       UInt_t im2 = veryloose_muon_indices[1];

       fill_muon_selection_flags((*muons)[im2],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im2].p4();
       lep3q = (*muons)[im2].charge();
       lep3id = (*muons)[im2].pdgId();
       lep3iso = muon_isolation((*muons)[im2],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 3 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im || i == im2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_electron_indices.size() >= 1){

       UInt_t ie2 = veryloose_electron_indices[0];

       fill_electron_selection_flags((*electrons)[ie2],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie2].p4();
       lep3q = (*electrons)[ie2].charge();
       lep3id = (*electrons)[ie2].pdgId();
       lep3iso = electron_isolation((*electrons)[ie2],PV,rho);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     }

     for(UInt_t i = 0; i < muons->size(); i++){

       if (i == im)
	 continue;

       if ((passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4) && !(flags & VetoV5))
	 flags = flags | VetoV5;

       else if (passLooseMuonSelectionV1((*muons)[i],PV)  && (*muons)[i].pt() > 10 && abs((*muons)[i].eta()) < 2.4 )
	 flags = flags | VetoV4;

     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if (i == ie)
	 continue;

       if ( passWLLJJVetoElectronId( (*electrons)[i], PV , rho, rhoHLTElectronSelection) && (*electrons)[i].pt() > 10 && abs((*electrons)[i].eta()) < 2.5  ) 
	 flags = flags | VetoV2;


     }

   } else if (tight_electron_indices.size() >= 1 && veryloose_electron_indices.size() >= 1){

     UInt_t i1 = tight_electron_indices[0];
     UInt_t i2 = veryloose_electron_indices[0];

     lep1= (*electrons)[i1].p4();
     lep1q = (*electrons)[i1].charge();
     lep1id = (*electrons)[i1].pdgId();
     lep1iso = electron_isolation((*electrons)[i1],PV,rho);


     lep2= (*electrons)[i2].p4();
     lep2q = (*electrons)[i2].charge();
     lep2id = (*electrons)[i2].pdgId();
     lep2iso = electron_isolation((*electrons)[i2],PV,rho);

     fill_electron_selection_flags((*electrons)[i1],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,1);
     fill_electron_selection_flags((*electrons)[i2],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,2);


     if (tight_muon_indices.size() >= 1){

       UInt_t im = tight_muon_indices[0];

       fill_muon_selection_flags((*muons)[im],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im].p4();
       lep3q = (*muons)[im].charge();
       lep3id = (*muons)[im].pdgId();
       lep3iso = muon_isolation((*muons)[im],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (tight_electron_indices.size() >= 2) {

       UInt_t i3 = tight_electron_indices[1];

       fill_electron_selection_flags((*electrons)[i3],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[i3].p4();
       lep3q = (*electrons)[i3].charge();
       lep3id = (*electrons)[i3].pdgId();
       lep3iso = electron_isolation((*electrons)[i3],PV,rho);
       
       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 3)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_muon_indices.size() >= 1){

       UInt_t im = veryloose_muon_indices[0];

       fill_muon_selection_flags((*muons)[im],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im].p4();
       lep3q = (*muons)[im].charge();
       lep3id = (*muons)[im].pdgId();
       lep3iso = muon_isolation((*muons)[im],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_electron_indices.size() >= 2){

       UInt_t i3 = veryloose_electron_indices[1];

       fill_electron_selection_flags((*electrons)[i3],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[i3].p4();
       lep3q = (*electrons)[i3].charge();
       lep3id = (*electrons)[i3].pdgId();
       lep3iso = electron_isolation((*electrons)[i3],PV,rho);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 3 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     }

     //     std::cout << "lepton_selection_flags & Lep2LooseSelectionV5 = " << bool(lepton_selection_flags & Lep2LooseSelectionV5) << std::endl;
     //     std::cout << "lepton_selection_flags & Lep2LooseSelectionV3 = " << bool(lepton_selection_flags & Lep2LooseSelectionV3) << std::endl;

     for(UInt_t i = 0; i < muons->size(); i++){

       if ((passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4) && !(flags & VetoV5))
	 flags = flags | VetoV5;

       else if (passLooseMuonSelectionV1((*muons)[i],PV)  && (*muons)[i].pt() > 10 && abs((*muons)[i].eta()) < 2.4 )
	 flags = flags | VetoV4;

     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if (i == i1 || i == i2)
	 continue;

       if ( passWLLJJVetoElectronId( (*electrons)[i],PV , rho,  rhoHLTElectronSelection)  && (*electrons)[i].pt() > 10 && abs((*electrons)[i].eta()) < 2.5) 
	 flags = flags | VetoV2;

     }


   } else if (veryloose_muon_indices.size() >= 2){

     UInt_t i1 = veryloose_muon_indices[0];
     UInt_t i2 = veryloose_muon_indices[1];

     fill_muon_selection_flags((*muons)[i1],PV,lepton_selection_flags,1);
     fill_muon_selection_flags((*muons)[i2],PV,lepton_selection_flags,2);
     
     lep1 = (*muons)[i1].p4();
     lep1q = (*muons)[i1].charge();
     lep1id = (*muons)[i1].pdgId();
     lep1iso = muon_isolation((*muons)[i1],PV);

     lep2 = (*muons)[i2].p4();
     lep2q = (*muons)[i2].charge();
     lep2id = (*muons)[i2].pdgId();
     lep2iso = muon_isolation((*muons)[i2],PV);

     if (tight_muon_indices.size() >= 1){

       UInt_t i3 = tight_muon_indices[0];

       fill_muon_selection_flags((*muons)[i3],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[i3].p4();
       lep3q = (*muons)[i3].charge();
       lep3id = (*muons)[i3].pdgId();
       lep3iso = muon_isolation((*muons)[i3],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 3 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == i1 || i == i2 || i == i3)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (tight_electron_indices.size() >= 1) {

       UInt_t ie = tight_electron_indices[0];

       fill_electron_selection_flags((*electrons)[ie],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie].p4();
       lep3q = (*electrons)[ie].charge();
       lep3id = (*electrons)[ie].pdgId();
       lep3iso = electron_isolation((*electrons)[ie],PV,rho);
       
       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 3 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == i1 || i == i2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_muon_indices.size() >= 3){

       UInt_t i3 = veryloose_muon_indices[2];

       fill_muon_selection_flags((*muons)[i3],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[i3].p4();
       lep3q = (*muons)[i3].charge();
       lep3id = (*muons)[i3].pdgId();
       lep3iso = muon_isolation((*muons)[i3],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 4 || veryloose_electron_indices.size() >= 1 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == i1 || i == i2 || i == i3)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_electron_indices.size() >= 1){

       UInt_t ie = veryloose_electron_indices[0];

       fill_electron_selection_flags((*electrons)[ie],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie].p4();
       lep3q = (*electrons)[ie].charge();
       lep3id = (*electrons)[ie].pdgId();
       lep3iso = electron_isolation((*electrons)[ie],PV,rho);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 3 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == i1 || i == i2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     }
     
     for(UInt_t i = 0; i < muons->size(); i++){

       if (i == i1 || i == i2)
	 continue;

       if ((passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4) && !(flags & VetoV5))
	 flags = flags | VetoV5;

       else if (passLooseMuonSelectionV1((*muons)[i],PV)  && (*muons)[i].pt() > 10 && abs((*muons)[i].eta()) < 2.4 )
	 flags = flags | VetoV4;

     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if ( passWLLJJVetoElectronId( (*electrons)[i], PV,  rho, rhoHLTElectronSelection )  && (*electrons)[i].pt() > 10 && abs((*electrons)[i].eta()) < 2.5) 
	 flags = flags | VetoV2;

     }


   } else if (veryloose_muon_indices.size() >= 1 && veryloose_electron_indices.size() >= 1){

     UInt_t im = veryloose_muon_indices[0];
     UInt_t ie = veryloose_electron_indices[0];

     fill_muon_selection_flags((*muons)[im],PV,lepton_selection_flags,1);
     fill_electron_selection_flags((*electrons)[ie],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,2);

     lep1 = (*muons)[im].p4();
     lep1q = (*muons)[im].charge();
     lep1id = (*muons)[im].pdgId();
     lep1iso = muon_isolation((*muons)[im],PV);

     lep2= (*electrons)[ie].p4();
     lep2q = (*electrons)[ie].charge();
     lep2id = (*electrons)[ie].pdgId();
     lep2iso = electron_isolation((*electrons)[ie],PV,rho);

     if (tight_muon_indices.size() >= 1){

       UInt_t im2 = tight_muon_indices[0];

       fill_muon_selection_flags((*muons)[im2],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im2].p4();
       lep3q = (*muons)[im2].charge();
       lep3id = (*muons)[im2].pdgId();
       lep3iso = muon_isolation((*muons)[im2],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im || i == im2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (tight_electron_indices.size() >= 1) {

       UInt_t ie2 = tight_electron_indices[0];

       fill_electron_selection_flags((*electrons)[ie2],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie2].p4();
       lep3q = (*electrons)[ie2].charge();
       lep3id = (*electrons)[ie2].pdgId();
       lep3iso = electron_isolation((*electrons)[ie2],PV,rho);
       
       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_muon_indices.size() >= 2){

       UInt_t im2 = veryloose_muon_indices[1];

       fill_muon_selection_flags((*muons)[im2],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im2].p4();
       lep3q = (*muons)[im2].charge();
       lep3id = (*muons)[im2].pdgId();
       lep3iso = muon_isolation((*muons)[im2],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 3 || veryloose_electron_indices.size() >= 2 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im || i == im2)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_electron_indices.size() >= 2){

       UInt_t ie2 = veryloose_electron_indices[1];

       fill_electron_selection_flags((*electrons)[ie2],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie2].p4();
       lep3q = (*electrons)[ie2].charge();
       lep3id = (*electrons)[ie2].pdgId();
       lep3iso = electron_isolation((*electrons)[ie2],PV,rho);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 3 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     }

     for(UInt_t i = 0; i < muons->size(); i++){

       if (i == im)
	 continue;

       if ((passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4) && !(flags & VetoV5))
	 flags = flags | VetoV5;

       else if (passLooseMuonSelectionV1((*muons)[i],PV)  && (*muons)[i].pt() > 10 && abs((*muons)[i].eta()) < 2.4 )
	 flags = flags | VetoV4;
     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if (i == ie)
	 continue;

       if ( passWLLJJVetoElectronId( (*electrons)[i], PV ,rho,  rhoHLTElectronSelection )  && (*electrons)[i].pt() > 10 && abs((*electrons)[i].eta()) < 2.5) 
	 flags = flags | VetoV2;


     }

   } else if (veryloose_electron_indices.size() >= 2){

     UInt_t i1 = veryloose_electron_indices[0];
     UInt_t i2 = veryloose_electron_indices[1];

     lep1= (*electrons)[i1].p4();
     lep1q = (*electrons)[i1].charge();
     lep1id = (*electrons)[i1].pdgId();
     lep1iso = electron_isolation((*electrons)[i1],PV,rho);


     lep2= (*electrons)[i2].p4();
     lep2q = (*electrons)[i2].charge();
     lep2id = (*electrons)[i2].pdgId();
     lep2iso = electron_isolation((*electrons)[i2],PV,rho);

     fill_electron_selection_flags((*electrons)[i1],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,1);
     fill_electron_selection_flags((*electrons)[i2],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,2);

     if (tight_muon_indices.size() >= 1){

       UInt_t im = tight_muon_indices[0];

       fill_muon_selection_flags((*muons)[im],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im].p4();
       lep3q = (*muons)[im].charge();
       lep3id = (*muons)[im].pdgId();
       lep3iso = muon_isolation((*muons)[im],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 2 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 3 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (tight_electron_indices.size() >= 1) {

       UInt_t ie = tight_electron_indices[0];

       fill_electron_selection_flags((*electrons)[ie],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[ie].p4();
       lep3q = (*electrons)[ie].charge();
       lep3id = (*electrons)[ie].pdgId();
       lep3iso = electron_isolation((*electrons)[ie],PV,rho);
       
       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 3 || tight_electron_indices.size() >= 2)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_muon_indices.size() >= 1){

       UInt_t im = veryloose_muon_indices[0];

       fill_muon_selection_flags((*muons)[im],PV,lepton_selection_flags,3);
       
       lep3 = (*muons)[im].p4();
       lep3q = (*muons)[im].charge();
       lep3id = (*muons)[im].pdgId();
       lep3iso = muon_isolation((*muons)[im],PV);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 2 || veryloose_electron_indices.size() >= 3 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (i == im)
	   continue;
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     } else if (veryloose_electron_indices.size() >= 3){

       UInt_t i3 = veryloose_electron_indices[2];

       fill_electron_selection_flags((*electrons)[i3],PV,rho,rhoHLTElectronSelection,lepton_selection_flags,3);
       
       lep3= (*electrons)[i3].p4();
       lep3q = (*electrons)[i3].charge();
       lep3id = (*electrons)[i3].pdgId();
       lep3iso = electron_isolation((*electrons)[i3],PV,rho);

       //fourth lepton veto
       if (tight_muon_indices.size() >= 1 || veryloose_muon_indices.size() >= 1 || veryloose_electron_indices.size() >= 4 || tight_electron_indices.size() >= 1)
	 flags = flags | VetoV12;

       for(UInt_t i = 0; i < muons->size(); i++){
	 
	 if (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4)
	   flags = flags | VetoV13;
       }
       
     }

     //     std::cout << "lepton_selection_flags & Lep2LooseSelectionV5 = " << bool(lepton_selection_flags & Lep2LooseSelectionV5) << std::endl;
     //     std::cout << "lepton_selection_flags & Lep2LooseSelectionV3 = " << bool(lepton_selection_flags & Lep2LooseSelectionV3) << std::endl;

     for(UInt_t i = 0; i < muons->size(); i++){

       if ( (passWLLJJVetoMuonId( (*muons)[i],PV) && (*muons)[i].pt() > 5 && abs((*muons)[i].eta()) < 2.4) && !(flags & VetoV5) )
	 flags = flags | VetoV5;

       else if (passLooseMuonSelectionV1((*muons)[i],PV)  && (*muons)[i].pt() > 10 && abs((*muons)[i].eta()) < 2.4 )
	 flags = flags | VetoV4;

     }

     for(UInt_t i = 0; i < electrons->size(); i++){

       if (i == i1 || i == i2)
	 continue;

       if ( passWLLJJVetoElectronId( (*electrons)[i],PV , rho,  rhoHLTElectronSelection)  && (*electrons)[i].pt() > 10 && abs((*electrons)[i].eta()) < 2.5) 
	 flags = flags | VetoV2;

     }

     

   } else
     return;

   lhe_and_gen_object.analyze(iEvent,lep1,lep2);

   std::vector<const pat::Jet *> cleaned_jets;

     for (const pat::Tau &tau : *taus) {
       if (tau.pt() < 18) continue;

       if (fabs(tau.eta()) > 2.3) continue;

       if ( reco::deltaR(tau.p4(),lep1) < 0.4 || reco::deltaR(tau.p4(),lep2) < 0.4 ) continue;

       //std::cout << tau.tauID("byMediumIsolationMVArun2v1DBnewDMwLT") << std::endl;
       //std::cout << tau.charge() << std::endl;

       if(tau.tauID("byMediumIsolationMVArun2v1DBnewDMwLT")){
	 flags = flags | VetoV3;
       }

       float iso_sum = 0;

       for(const reco::CandidatePtr& gamma_cand: tau.isolationGammaCands() ){
	 iso_sum = iso_sum + gamma_cand->pt();
       }

       for(const reco::CandidatePtr& charged_hadron_cand: tau.isolationChargedHadrCands()){
	   iso_sum = iso_sum + charged_hadron_cand->pt();
	 }

       for(const reco::CandidatePtr& neutral_hadron_cand: tau.isolationNeutrHadrCands()){
	 iso_sum = iso_sum + neutral_hadron_cand->pt();
       }

       if(tau.tauID("decayModeFinding") && tau.tauID("decayModeFindingNewDMs") && iso_sum < 5.0){
	 flags = flags | VetoV7;
       }
     }


    
   std::vector<float> corrected_jet_pts;

   for (const pat::Jet &j : *jets) {

     if (abs(j.eta()) > 4.7)
       continue;

     if ( reco::deltaR(j.p4(),lep1) < 0.4 || reco::deltaR(j.p4(),lep2) < 0.4 ){
       continue;

     }

     if (lep3id != 0){       
       if ( reco::deltaR(j.p4(),lep3) < 0.4 ){
	 continue;       
       }       
     }
       
     if (isMC_){

       jecUnc_->setJetEta(j.eta());
       jecUnc_->setJetPt(j.pt());
       
       float jecunc = jecUnc_->getUncertainty(true);
       
       float oldpt = j.pt();
       
       if (jes_ == "up")
	 oldpt = oldpt*(1+jecunc);
       else if (jes_ == "down")
	 oldpt = oldpt*(1-jecunc);
       else
	 assert(jes_ == "nominal");
       
       JME::JetResolution resolution = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
       JME::JetResolutionScaleFactor resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");
       JME::JetParameters jpar;
       jpar.setJetPt( j.pt()).setJetEta( j.eta() ).setRho( rho)  ;
       float res = resolution.getResolution(jpar);
       float sf = resolution_sf.getScaleFactor(jpar);
       float sf_up = resolution_sf.getScaleFactor(jpar, Variation::UP);
       float sf_down = resolution_sf.getScaleFactor(jpar, Variation::DOWN);
       
       // https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
       const reco::GenJet* genjet= j.genJet() ;
       float newpt = j.pt();
       
       float newptUp = j.pt();
       float newptDown = j.pt();
       
       //
       // -- if (genjet==NULL) cout <<"DEBUG: NO GENJET"<<endl;
       // -- else cout <<"DEBUG: YES GENJET"<<"DR="<<reco::deltaR(*genjet,j) <<" DPT="<<fabs(j.pt()-genjet->pt())<< " 3Sig="<<3*res*genjet->pt()<<endl;
       
       // additional matching dR(reco jet, gen jet)<Rcone/2 and dpt=abs(pT-pTgen)<3*sigma_MC
       if (genjet!=NULL and reco::deltaR(*genjet,j) <0.2 and fabs(j.pt()-genjet->pt())<3*res*genjet->pt()) { //scaling
	 float genpt=genjet->pt();
	 newpt     = std::max(float(0.),float(genpt+sf*(j.pt()-genpt)      ));
	 newptUp   = std::max(float(0.),float(genpt+sf_up*(j.pt()-genpt)   ));
	 newptDown = std::max(float(0.),float(genpt+sf_down*(j.pt()-genpt) ));
       }
       else{ // smearing
	 float sigma=TMath::Sqrt( sf*sf -1) * res;
	 float sigmaUp=TMath::Sqrt( sf_up*sf_up -1) * res;
	 float sigmaDown=TMath::Sqrt( sf_down*sf_down -1) * res;
	 float s = rnd_->Gaus(0,1);
	 newpt = s * sigma + oldpt;
	 newptUp = s * sigmaUp + oldpt;
	 newptDown = s * sigmaDown +oldpt;
       }
       //*(TLorentzVector*)(*p4)[p4->GetEntriesFast()-1] *= newpt /oldpt ; // Update the pt stored in p4 at the last position (N-1)
       //ptResUncUp->push_back(newptUp);
       //ptResUncDown->push_back(newptDown) ;
       
       if (jer_ == "nominal")
	 corrected_jet_pts.push_back(newpt);
       else if (jer_ == "up")
	 corrected_jet_pts.push_back(newptUp);
       else if (jer_ == "down")
	 corrected_jet_pts.push_back(newptDown);
       else if (jer_ == "no_correction")
	 corrected_jet_pts.push_back(j.pt());
       else
	 assert(0);
       
     }
     else
       corrected_jet_pts.push_back(j.pt());     

     cleaned_jets.push_back(&j);
   }

  maxbtagevent = std::numeric_limits<Float_t>::min();
   
  for (UInt_t i = 0; i < cleaned_jets.size(); i++) {
    if (cleaned_jets[i]->pt() < 20) continue;

    //std::cout << cleaned_jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;

    if (std::max(0.f,cleaned_jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")) > maxbtagevent)
      maxbtagevent = std::max(0.f,cleaned_jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

    //see here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X
    //if( std::max(0.f,j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")) > 0.814)
    //  return;
    
  }


   if (cleaned_jets.size() < 2) 
     return;


   std::pair<int,int> highest_pt_jet_indices = get_two_highest_pt_jet_indices(corrected_jet_pts);

   if (cleaned_jets[highest_pt_jet_indices.first]->pt() < 30 || cleaned_jets[highest_pt_jet_indices.second]->pt() < 30)
     return;

   float NHF0    = cleaned_jets[highest_pt_jet_indices.first]->neutralHadronEnergyFraction();
   float NEMF0   = cleaned_jets[highest_pt_jet_indices.first]->neutralEmEnergyFraction();
   float CHF0    = cleaned_jets[highest_pt_jet_indices.first]->chargedHadronEnergyFraction();
   //float MUF0    = cleaned_jets[0]->muonEnergyFraction();
   float CEMF0   = cleaned_jets[highest_pt_jet_indices.first]->chargedEmEnergyFraction();
   int NumConst0 = cleaned_jets[highest_pt_jet_indices.first]->chargedMultiplicity()+cleaned_jets[highest_pt_jet_indices.first]->neutralMultiplicity();
   int CHM0      = cleaned_jets[highest_pt_jet_indices.first]->chargedMultiplicity();
   
   float NHF1    = cleaned_jets[1]->neutralHadronEnergyFraction();
   float NEMF1   = cleaned_jets[1]->neutralEmEnergyFraction();
   float CHF1    = cleaned_jets[1]->chargedHadronEnergyFraction();
   //float MUF1    = cleaned_jets[1]->muonEnergyFraction();
   float CEMF1   = cleaned_jets[1]->chargedEmEnergyFraction();
   int NumConst1 = cleaned_jets[1]->chargedMultiplicity()+cleaned_jets[1]->neutralMultiplicity();
   int CHM1      = cleaned_jets[1]->chargedMultiplicity();

   if (abs(cleaned_jets[highest_pt_jet_indices.first]->eta()) <= 3.0)
     jet1loosejetid = (NHF0<0.99 && NEMF0<0.99 && NumConst0>1) && ((abs(cleaned_jets[highest_pt_jet_indices.first]->eta())<=2.4 && CHF0>0 && CHM0>0 && CEMF0<0.99) || abs(cleaned_jets[highest_pt_jet_indices.first]->eta())>2.4);
   else
     jet1loosejetid = NEMF0<0.90 && cleaned_jets[highest_pt_jet_indices.first]->neutralMultiplicity()>10;

   if (abs(cleaned_jets[highest_pt_jet_indices.second]->eta()) <= 3.0)
     jet2loosejetid = (NHF1<0.99 && NEMF1<1.99 && NumConst1>1) && ((abs(cleaned_jets[highest_pt_jet_indices.second]->eta())<=2.4 && CHF1>0 && CHM1>0 && CEMF1<0.99) || abs(cleaned_jets[highest_pt_jet_indices.second]->eta())>2.4);
   else
     jet2loosejetid = NEMF1<0.90 && cleaned_jets[highest_pt_jet_indices.second]->neutralMultiplicity()>10;
   
   //   std::cout << "andrew debug 1" << std::endl;
   //   std::cout << corrected_jet_pts[highest_pt_jet_indices.first]/cleaned_jets[highest_pt_jet_indices.first]->pt() << std::endl;
   //   std::cout << corrected_jet_pts[highest_pt_jet_indices.second]/cleaned_jets[highest_pt_jet_indices.second]->pt() << std::endl;
   //   std::cout << "andrew debug 2" << std::endl;

   jet1=cleaned_jets[highest_pt_jet_indices.first]->p4();
   jet2=cleaned_jets[highest_pt_jet_indices.second]->p4();

   jet1 *= (corrected_jet_pts[highest_pt_jet_indices.first]/cleaned_jets[highest_pt_jet_indices.first]->pt());
   jet2 *= (corrected_jet_pts[highest_pt_jet_indices.second]/cleaned_jets[highest_pt_jet_indices.second]->pt());

   jet1btag = std::max(0.f,cleaned_jets[highest_pt_jet_indices.first]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
   jet2btag = std::max(0.f,cleaned_jets[highest_pt_jet_indices.second]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
   jet1pujetid = cleaned_jets[highest_pt_jet_indices.first]->userFloat("pileupJetId:fullDiscriminant");
   jet2pujetid = cleaned_jets[highest_pt_jet_indices.second]->userFloat("pileupJetId:fullDiscriminant");

   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);

   assert(mets->size() == 1);

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
ntuple_maker::beginJob()
{

  edm::Service<TFileService> fs;

  if (isMC_){
    n_events_run_over= fs->make<TH1F>("n_events_run_over","n_events_run_over",1,0,1);
    n_weighted_events_run_over= fs->make<TH1F>("n_weighted_events_run_over","n_weighted_events_run_over",1,0,1);
  }

  lhe_and_gen_object.initrwgt_header_tree_ = fs->make<TTree >("initrwgt_header","initrwgt_header");

  lhe_and_gen_object.slha_header_tree_ = fs->make<TTree >("slha_header","slha_header");

  tree = fs->make<TTree>( "events"  , "events");

  tree->Branch("lepton_selection_flags",&lepton_selection_flags);
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
  tree->Branch("lep1iso",&lep1iso);
  tree->Branch("lep2iso",&lep2iso);
  tree->Branch("lep3iso",&lep3iso);
  tree->Branch("nvtx",&nvtx);
  tree->Branch("lep1q",&lep1q);
  tree->Branch("lep1id",&lep1id);
  tree->Branch("lep2q",&lep2q);
  tree->Branch("lep2id",&lep2id);
  tree->Branch("lep3q",&lep3q);
  tree->Branch("lep3id",&lep3id);

  tree->Branch("maxbtagevent",&maxbtagevent);

  if (isMC_)
    tree->Branch("n_pu_interactions",&n_pu_interactions);

  lhe_and_gen_object.defineBranches(tree);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ntuple_maker::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------


void 
ntuple_maker::beginRun(edm::Run const& iRun, edm::EventSetup const&)
{

  lhe_and_gen_object.beginRun(iRun);

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
ntuple_maker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ntuple_maker);
