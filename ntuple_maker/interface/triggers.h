#ifndef triggers_H
#define triggers_H


#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"


bool trigger_fired(const edm::TriggerNames &names, const edm::Handle< edm::TriggerResults> &triggerResultsHandle, std::string which_triggers);

#endif
