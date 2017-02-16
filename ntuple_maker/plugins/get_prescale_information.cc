// -*- C++ -*-
//
// Package:    get_prescale_information
// Class:      get_prescale_information
// 
/**\class get_prescale_information get_prescale_information.cc get_prescale_information/get_prescale_information/plugins/get_prescale_information.cc

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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/LorentzVector.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//
// class declaration
//

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class get_prescale_information : public edm::EDAnalyzer {
public:

  
      explicit get_prescale_information(const edm::ParameterSet&);
      ~get_prescale_information();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

  // ----------member data ---------------------------

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesL1max_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesL1min_;

  TTree * tree;

  Int_t l1min;
  Int_t l1max;
  Int_t prescale;

  ULong64_t event;
  UInt_t run;
  UInt_t lumi;

  UInt_t current_run;
  UInt_t current_lumi;

  Bool_t first_event;

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
get_prescale_information::get_prescale_information(const edm::ParameterSet& iConfig):
  triggerBits_(consumes< edm::TriggerResults >(edm::InputTag("TriggerResults","","HLT"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  triggerPrescalesL1max_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("l1max"))),
  triggerPrescalesL1min_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("l1min"))),
  first_event(true)
{
  //now do what ever initialization is needed

}


get_prescale_information::~get_prescale_information()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
get_prescale_information::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if ( (iEvent.eventAuxiliary().run() == run && iEvent.eventAuxiliary().luminosityBlock() == lumi) && !first_event)
    return;

  first_event = false;

  run=iEvent.eventAuxiliary().run();
  lumi=iEvent.eventAuxiliary().luminosityBlock();
  event=iEvent.eventAuxiliary().event();

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1max;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1min;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  iEvent.getByToken(triggerPrescalesL1max_, triggerPrescalesL1max);
  iEvent.getByToken(triggerPrescalesL1min_, triggerPrescalesL1min);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

  int n_matches = 0;

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {

    //std::cout << "names.triggerName(i) = " << names.triggerName(i) << std::endl;

    //if (names.triggerName(i).find("HLT_Mu17_TrkIsoVVL_v") == std::string::npos)
    //  continue;

    if (names.triggerName(i).find("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") == std::string::npos)
      continue;

    n_matches++;

    //if (names.triggerName(i).find("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v") == std::string::npos)
    // continue;

    prescale = triggerPrescales->getPrescaleForIndex(i);
    l1max = triggerPrescalesL1max->getPrescaleForIndex(i);
    l1min = triggerPrescalesL1min->getPrescaleForIndex(i);

    /*

    std::cout << "Trigger " << names.triggerName(i) << 
      ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
      ", l1max " << triggerPrescalesL1max->getPrescaleForIndex(i) <<
      ", l1min " << triggerPrescalesL1min->getPrescaleForIndex(i) <<
      ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") 
	      << std::endl;

    */

  }

  if (n_matches != 1){

    std::cout << "n_matches = " << n_matches << std::endl;

    assert(n_matches == 1);

  }

  tree->Fill();

  //triggerPrescales->

  //prescale = triggerPrescales->getPrescaleForName("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v",true);
  //l1min = triggerPrescalesL1min->getPrescaleForName("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v",true);
  //l1max = triggerPrescalesL1max->getPrescaleForName("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v",true);



}


// ------------ method called once each job just before starting event loop  ------------
void 
get_prescale_information::beginJob()
{

  edm::Service<TFileService> fs;

  tree = fs->make<TTree>( "prescales"  , "prescales");

  tree->Branch("prescale",&prescale);
  tree->Branch("l1min",&l1min);
  tree->Branch("l1max",&l1max);  
  //tree->Branch("event",&event);
  tree->Branch("lumi",&lumi);
  tree->Branch("run",&run);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
get_prescale_information::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------


void 
get_prescale_information::beginRun(edm::Run const& iRun, edm::EventSetup const&)
{

}


// ------------ method called when ending the processing of a run  ------------
 /*
void 
ntuple_maker::endRun(edm::Run const&, edm::EventSetup const&)
{



}
 */
// ------------ method called when starting to processes a luminosity block  ------------

void 
get_prescale_information::beginLuminosityBlock(edm::LuminosityBlock const&iLumi, edm::EventSetup const&iSetup)
{

  /*

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1max;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1min;

  //iLumi.getByToken(triggerBits_, triggerBits);
  iLumi.getByToken(triggerPrescales_, triggerPrescales);
  iLumi.getByToken(triggerPrescalesL1max_, triggerPrescalesL1max);
  iLumi.getByToken(triggerPrescalesL1min_, triggerPrescalesL1min);

  //const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  std::cout << "\n === TRIGGER PATHS === " << std::endl;
    //    if (names.triggerName(i).find("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v") == std::string::npos)
    //      continue;

    std::cout << // "Trigger " << names.triggerName(i) << 
      ", prescale " << triggerPrescales->getPrescaleForName("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v",true) <<
      ", l1max " << triggerPrescalesL1max->getPrescaleForName("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v",true)  <<
      ", l1min " << triggerPrescalesL1min->getPrescaleForName("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v",true) 
	      << std::endl;

  */

}

// ------------ method called when ending the processing of a luminosity block  ------------

void 
get_prescale_information::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
get_prescale_information::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(get_prescale_information);
