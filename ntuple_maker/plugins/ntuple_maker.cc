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

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

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

  //using a token instead of label does not work for the lheruninfo, see here: https://hypernews.cern.ch/HyperNews/CMS/get/edmFramework/3319/2.html
  edm::InputTag lheRunInfoLabel_;

  TH1F * n_events_run_over;
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
  lhe_and_gen lhe_and_gen_object; //separate the part that runs over the generator and lhe information

  std::vector<int> pdf_weight_indices;
  int qcd_weight_up_index;
  int qcd_weight_down_index;

  bool syscalcinfo_;

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
  triggerResultsToken_(consumes< edm::TriggerResults >(edm::InputTag("TriggerResults","","HLT"))),
  triggerObjectToken_( consumes< pat::TriggerObjectStandAloneCollection >(edm::InputTag("selectedPatTrigger"))),
  lheRunInfoLabel_(iConfig.getParameter<edm::InputTag>("lheruninfo")),
  syscalcinfo_(iConfig.getUntrackedParameter<bool>("syscalcinfo"))
{
  //now do what ever initialization is needed

  lhe_and_gen_object.prunedGenToken_ = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedgenparticles"));
  lhe_and_gen_object.packedGenToken_ = consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedgenparticles"));
  lhe_and_gen_object.lheEvtToken_ = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheevent"));
  lhe_and_gen_object.syscalcinfo_ = syscalcinfo_;

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

  n_events_run_over->Fill(0.5);

  edm::Handle< edm::TriggerResults> triggerResultsHandle;

  iEvent.getByToken(triggerResultsToken_,triggerResultsHandle);

  std::vector<std::string> triggerNames;

  triggerNames.push_back("HLT_Mu17_Mu8");
  triggerNames.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL");
  triggerNames.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL");
  triggerNames.push_back("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL");
  triggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL");

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResultsHandle);

  edm::Handle<pat::TriggerObjectStandAloneCollection > triggerObjectHandle;

  Bool_t trigger_fired = kFALSE;

  for (unsigned int i = 0; i < names.size(); i++) {

    for(unsigned int j=0;j< triggerNames.size() ;++j){

      std::string name = names.triggerName(i);

      if (name.find( (triggerNames)[j]) != std::string::npos){

	trigger_fired = kTRUE;
	
      }
    }
  }

  if (! trigger_fired)
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

  std::vector<UInt_t> loose_muon_indices;
  std::vector<UInt_t> loose_electron_indices;

  flags = 0;

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

   for(UInt_t i = 0; i < electrons->size(); i++){

     if( (*electrons)[i].pt() < 10) 
       continue;

     if (!passLooseElectronId((*electrons)[i],PV))
       continue;

     loose_electron_indices.push_back(i);

   }

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);

   for(UInt_t i = 0; i < muons->size(); i++){
     if ((*muons)[i].pt() < 10)
       continue;

     Float_t relative_isolation = ((*muons)[i].pfIsolationR04().sumChargedHadronPt+ std::max(0.0,(*muons)[i].pfIsolationR04().sumNeutralHadronEt + (*muons)[i].pfIsolationR04().sumPhotonEt - 0.5 * (*muons)[i].pfIsolationR04().sumPUPt))/(*muons)[i].pt();

     if (! (passLooseMuonId((*muons)[i],PV) && relative_isolation < 1.0) )
       continue;

     loose_muon_indices.push_back(i);

    }

   if(loose_muon_indices.size() >= 2){

     UInt_t i1 = loose_muon_indices[0];
     UInt_t i2 = loose_muon_indices[1];

     Float_t relative_isolation_1 = ((*muons)[i1].pfIsolationR04().sumChargedHadronPt+ std::max(0.0,(*muons)[i1].pfIsolationR04().sumNeutralHadronEt + (*muons)[i1].pfIsolationR04().sumPhotonEt - 0.5 * (*muons)[i1].pfIsolationR04().sumPUPt))/(*muons)[i1].pt();
     Float_t relative_isolation_2 = ((*muons)[i2].pfIsolationR04().sumChargedHadronPt+ std::max(0.0,(*muons)[i2].pfIsolationR04().sumNeutralHadronEt + (*muons)[i2].pfIsolationR04().sumPhotonEt - 0.5 * (*muons)[i2].pfIsolationR04().sumPUPt))/(*muons)[i2].pt();

     if (passTightMuonId((*muons)[i1],PV) && relative_isolation_1 < 0.12) 
       flags = flags | Lep1TightSelectionV1;

     if (passTightMuonId((*muons)[i2],PV) && relative_isolation_2 < 0.12) 
       flags = flags | Lep2TightSelectionV1;

     
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

     if (passTightMuonId((*muons)[im],PV) && relative_isolation_1 < 0.12) 
       flags = flags | Lep1TightSelectionV1;

     lep1 = (*muons)[im].p4();
     lep1q = (*muons)[im].charge();
     lep1id = (*muons)[im].pdgId();

     if (passTightElectronId((*electrons)[ie], PV))
       flags = flags | Lep2TightSelectionV1;
       
     lep2= (*electrons)[ie].p4();
     lep2q = (*electrons)[ie].charge();
     lep2id = (*electrons)[ie].pdgId();
     
     
   } else if (loose_electron_indices.size() >= 2){

     UInt_t i1 = loose_electron_indices[0];
     UInt_t i2 = loose_electron_indices[1];

     lep1= (*electrons)[i1].p4();
     lep1q = (*electrons)[i1].charge();
     lep1id = (*electrons)[i1].pdgId();

     lep2= (*electrons)[i2].p4();
     lep2q = (*electrons)[i2].charge();
     lep2id = (*electrons)[i2].pdgId();

     if (passTightElectronId((*electrons)[i1],PV))
	 flags = flags | Lep1TightSelectionV1;
     if (passTightElectronId((*electrons)[i2],PV))
	 flags = flags | Lep2TightSelectionV1;

   } else 
     return;


   lhe_and_gen_object.analyze(iEvent,lep1,lep2,pdf_weight_indices,qcd_weight_up_index,qcd_weight_down_index);

   std::vector<const pat::Jet *> cleaned_jets;

   for (const pat::Jet &j : *jets) {

     if ( reco::deltaR(j.p4(),lep1) < 0.5 || reco::deltaR(j.p4(),lep2) < 0.5 ){
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
   metgenmetpt = met.genMET()->pt();
   metptshiftup = met.shiftedPt(pat::MET::JetEnUp);
   metptshiftdown = met.shiftedPt(pat::MET::JetEnDown);

   tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
ntuple_maker::beginJob()
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

  tree->Branch("maxbtagevent",&maxbtagevent);

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

  edm::Handle<LHERunInfoProduct> hLheRun;
  iRun.getByLabel(lheRunInfoLabel_,hLheRun);


  for ( LHERunInfoProduct::headers_const_iterator lheruniter = hLheRun.product()->headers_begin(); lheruniter != hLheRun.product()->headers_end(); lheruniter++ ) {
    
    std::cout << "lheruniter->tag() = " << lheruniter->tag() << std::endl;

    if (lheruniter->tag() != "initrwgt")
      continue;

    bool in_NNPDF23_lo_as_0130_qed = false;

    for ( LHERunInfoProduct::Header::const_iterator iter = lheruniter->begin(); iter != lheruniter->end(); iter++ ) {

      if ( (*iter).find("mur=2 muf=2") != std::string::npos){

	std::stringstream ss;
	ss << (*iter).substr((*iter).find("<weight id=")+12,(*iter).find("> mu") - (*iter).find("<weight id=") - 13);
	ss >> qcd_weight_up_index;
	qcd_weight_up_index=qcd_weight_up_index-1;
	continue;
      }

      if ( (*iter).find("mur=0.5 muf=0.5") != std::string::npos){

	std::stringstream ss;
	ss << (*iter).substr((*iter).find("<weight id=")+12,(*iter).find("> mu") - (*iter).find("<weight id=") - 13);
	ss >> qcd_weight_down_index;
	qcd_weight_down_index=qcd_weight_down_index-1;
	continue;
      }

      if ( (*iter).find("NNPDF23_lo_as_0130_qed.LHgrid") != std::string::npos){
	in_NNPDF23_lo_as_0130_qed = true;
	continue;
      }

      if (in_NNPDF23_lo_as_0130_qed){
	if ( (*iter).find("/weightgroup") != std::string::npos){
	  in_NNPDF23_lo_as_0130_qed = false;
	  continue;
	}

	assert((*iter).find("Member") != std::string::npos);

	//std::cout << (*iter).substr((*iter).find("<weight id=")+12,(*iter).find(">Member") - (*iter).find("<weight id=") - 13) << std::endl;

	int weight_index;

	std::stringstream ss;
	ss << (*iter).substr((*iter).find("<weight id=")+12,(*iter).find(">Member") - (*iter).find("<weight id=") - 13);
	ss >> weight_index;

	//need to subtract one because the weight numbers start from 1
	pdf_weight_indices.push_back(weight_index-1);

      }

      //std::cout << (*iter) << std::endl;
    }

  }

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
