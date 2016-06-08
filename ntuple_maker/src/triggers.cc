#ifndef triggers_H
#define triggers_H


#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "ntuple_maker/ntuple_maker/interface/triggers.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

bool trigger_fired(const edm::TriggerNames &names, const edm::Handle< edm::TriggerResults> &triggerResultsHandle, std::string which_triggers){

  assert(which_triggers == "doubleeg" || which_triggers == "doublemu" || which_triggers == "muoneg" || which_triggers == "doublelepton" || which_triggers == "electron_fake_rate" || which_triggers == "muon_fake_rate");
  
  std::vector<std::string> triggerNames;

  if (which_triggers == "doubleeg" || which_triggers == "doublelepton"){

    triggerNames.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

  } 


if (which_triggers == "doublemu" || which_triggers == "doublelepton") {

    triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
    triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");

  } 

 if (which_triggers == "muoneg" || which_triggers == "doublelepton") {

    triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    triggerNames.push_back("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v");

  }

 if (which_triggers == "muon_fake_rate") {

   //triggerNames.push_back("HLT_Mu8_v");
   
   //triggerNames.push_back("HLT_Mu17_v");
   triggerNames.push_back("HLT_Mu17_TrkIsoVVL_v");
   
   //triggerNames.push_back("HLT_Mu24_v"); //not present for some of the 2015 data
   //triggerNames.push_back("HLT_Mu34_v"); //not present for some of the 2015 data

  }

 if (which_triggers == "electron_fake_rate") {

   triggerNames.push_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v");

   //triggerNames.push_back("HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v"); //not present for some of the 2015 data

   //triggerNames.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");


   //triggerNames.push_back("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v");

  }

  edm::Handle<pat::TriggerObjectStandAloneCollection > triggerObjectHandle;
  
  Bool_t trigger_fired = kFALSE;
  
  for (unsigned int i = 0; i < names.size(); i++) {
    
    //if (triggerResultsHandle->accept(i))
    //	std::cout << "names.triggerName(i) = " << names.triggerName(i) << std::endl;                                                                                     
      

    
    for(unsigned int j=0;j< triggerNames.size() ;++j){
      
      std::string name = names.triggerName(i);
      

      if (name.find( (triggerNames)[j]) != std::string::npos && triggerResultsHandle->accept(i)){
	//      if (name.find( (triggerNames)[j]) != std::string::npos){                                                                                               	
	trigger_fired = kTRUE;
	
      }
    }
  }

  return trigger_fired;

}

#endif
