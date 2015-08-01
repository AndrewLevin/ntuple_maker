import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

process.maxEvents = cms.untracked.PSet( 
input = cms.untracked.int32(100) ,
#input = cms.untracked.int32(-1) ,

)

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",

#skipEvents= cms.untracked.uint32(1),

    # replace 'myfile.root' with the source file you want to use
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(

'/store/mc/RunIISpring15DR74/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/003F747A-D835-E511-BFC8-002590D9D976.root'

#'/store/mc/RunIISpring15DR74/WW_DoubleScattering_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/30000/04DCEB17-A12F-E511-8BFB-782BCB161FC2.root'


    ),
#eventsToProcess = cms.untracked.VEventRange('1:11:2092-1:11:2092'),

)

process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
   src = cms.InputTag("slimmedMuons"),
   preselection = cms.string("track.isNonnull"),
   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
   fractionOfSharedSegments = cms.double(0.499))

process.demo = cms.EDAnalyzer('ntuple_maker',

  syscalcinfo = cms.untracked.bool (True), #fill the information from syscalc
  lheinfo = cms.untracked.bool (True),
  vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  muons = cms.InputTag("cleanedMu"),
  lheevent = cms.InputTag("externalLHEProducer"),
  lheruninfo = cms.InputTag("externalLHEProducer"),
#  lheruninfo = cms.InputTag("source"),
#  lheevent = cms.InputTag("source"),
  #muons = cms.InputTag("slimmedMuons"),
  electrons = cms.InputTag("slimmedElectrons"),
  taus = cms.InputTag("slimmedTaus"),
  photons = cms.InputTag("slimmedPhotons"),
  jets = cms.InputTag("slimmedJets"),
  fatjets = cms.InputTag("slimmedJetsAK8"),
  mets = cms.InputTag("slimmedMETs"),  
  prunedgenparticles = cms.InputTag("prunedGenParticles"),  
  packedgenparticles = cms.InputTag("packedGenParticles"),  
                              
)

process.p = cms.Path(process.cleanedMu*process.demo)


