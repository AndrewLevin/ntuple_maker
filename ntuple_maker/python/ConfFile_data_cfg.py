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

#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v4/000/258/703/00000/A4DA0EA2-2172-E511-986D-02163E011AB8.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v4/000/258/703/00000/B21074A3-2172-E511-8603-02163E014135.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v4/000/258/703/00000/EA5BA99F-2172-E511-A486-02163E011A96.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v4/000/258/703/00000/F04EF5EF-AA71-E511-ABE7-02163E0126EE.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v4/000/258/703/00000/F61432C0-2172-E511-AE8B-02163E014237.root'

'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v4/000/258/159/00000/027612B0-306C-E511-BD47-02163E014496.root'

#'/store/data/Run2015D/DoubleMuon/MINIAOD/PromptReco-v4/000/258/159/00000/0C6D4AB0-6F6C-E511-8A64-02163E0133CD.root'

#'/store/mc/RunIISpring15DR74/WLLJJToLNu_M-4to60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/0EC87D48-E861-E511-B9FB-002590A3C97E.root'


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
  lheinfo = cms.untracked.bool (False),
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
  isMC = cms.untracked.bool(False),  
  prunedgenparticles = cms.InputTag("prunedGenParticles"),  
  packedgenparticles = cms.InputTag("packedGenParticles"),  
  pfCands = cms.InputTag("packedPFCandidates"),                              
pileup_summary = cms.InputTag("slimmedAddPileupInfo"),
genevent = cms.InputTag("generator"),
rho = cms.InputTag("fixedGridRhoFastjetAll")

)

process.p = cms.Path(process.cleanedMu*process.demo)



