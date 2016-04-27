
import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

process.maxEvents = cms.untracked.PSet( 
input = cms.untracked.int32(-1) ,
#input = cms.untracked.int32(-1) ,

)

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",

#skipEvents= cms.untracked.uint32(1),

    # replace 'myfile.root' with the source file you want to use
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),

    fileNames = cms.untracked.vstring(

'/store/mc/RunIIFall15MiniAODv2/WmWmJJ_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/80000/84FEA6AF-ABED-E511-AD3F-001E67581494.root'

#'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/00DF0A73-17C2-E511-B086-E41D2D08DE30.root'

#'file:/afs/cern.ch/work/a/anlevin/VBS/13_tev/Merged.root'

#'/store/mc/RunIISpring15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/9E0EDDCD-BE70-E511-9A85-6C3BE5B5F210.root'

    ),
#eventsToProcess = cms.untracked.VEventRange('1:11:2092-1:11:2092'),

)

process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
   src = cms.InputTag("slimmedMuons"),
   preselection = cms.string("track.isNonnull"),
   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
   fractionOfSharedSegments = cms.double(0.499))

process.demo = cms.EDAnalyzer('ntuple_maker',

  syscalcinfo = cms.untracked.bool (False), #fill the information from syscalc
  mgreweightinfo = cms.untracked.bool (False),
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
  isMC = cms.untracked.bool(True),  
  prunedgenparticles = cms.InputTag("prunedGenParticles"),  
  packedgenparticles = cms.InputTag("packedGenParticles"),  
  pfCands = cms.InputTag("packedPFCandidates"),
  pileup_summary = cms.InputTag("slimmedAddPileupInfo"),
  genevent = cms.InputTag("generator"),
  rho = cms.InputTag("fixedGridRhoFastjetAll")                              
)

process.p = cms.Path(process.cleanedMu*process.demo)

