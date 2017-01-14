import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

process.maxEvents = cms.untracked.PSet( 
#input = cms.untracked.int32(10000) ,
input = cms.untracked.int32(-1) ,

)

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",

#skipEvents= cms.untracked.uint32(1),

    # replace 'myfile.root' with the source file you want to use
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),

    fileNames = cms.untracked.vstring(
'/store/data/Run2016B/ZeroBias/MINIAOD/23Sep2016-v3/00000/02213046-F097-E611-9F99-0242AC130004.root'

    ),
#eventsToProcess = cms.untracked.VEventRange('1:2437:158389-1:2437:158389'),

)

process.demo = cms.EDAnalyzer('get_prescale_information',

prescales = cms.InputTag("patTrigger"),
l1min = cms.InputTag("patTrigger","l1min"),
l1max = cms.InputTag("patTrigger","l1max")

)

process.p = cms.Path(process.demo)



