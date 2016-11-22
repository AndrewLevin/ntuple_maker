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

'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/275/310/00000/F8C0DD37-AD37-E611-A8B8-02163E01198C.root',
'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/274/316/00000/A0C12CED-262A-E611-AD94-02163E011DF4.root',
'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/275/068/00000/4C361213-D234-E611-946E-02163E011DC3.root',
'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/275/375/00000/B8C3D0C5-8039-E611-BC04-02163E0143CF.root',
'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/425/00000/0E9B07A7-CC1B-E611-AF52-02163E012247.root',
'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/274/319/00000/26F48DC7-462A-E611-8C0E-02163E0139A0.root'


#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v4/000/258/703/00000/A4DA0EA2-2172-E511-986D-02163E011AB8.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v4/000/258/703/00000/B21074A3-2172-E511-8603-02163E014135.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v4/000/258/703/00000/EA5BA99F-2172-E511-A486-02163E011A96.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v4/000/258/703/00000/F04EF5EF-AA71-E511-ABE7-02163E0126EE.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v4/000/258/703/00000/F61432C0-2172-E511-AE8B-02163E014237.root'

#'/store/data/Run2015D/MuonEG/MINIAOD/16Dec2015-v1/00000/1ED2F05D-E8AF-E511-BDF3-0002C94D575E.root'

#'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v1/000/272/774/00000/C4E46C54-DE15-E611-A30E-02163E011E58.root'

#'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v1/000/272/728/00000/00375458-4215-E611-9618-02163E0138DD.root'

#'/store/data/Run2015D/DoubleMuon/MINIAOD/PromptReco-v4/000/258/159/00000/0C6D4AB0-6F6C-E511-8A64-02163E0133CD.root'

#'/store/mc/RunIISpring15DR74/WLLJJToLNu_M-4to60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/0EC87D48-E861-E511-B9FB-002590A3C97E.root'

'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/2227B110-EE19-E611-9B2C-02163E0145A3.root'

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
  isMC = cms.untracked.bool(False),  
  prunedgenparticles = cms.InputTag("prunedGenParticles"),  
  packedgenparticles = cms.InputTag("packedGenParticles"),  
  pfCands = cms.InputTag("packedPFCandidates"),                              
apply_trigger = cms.untracked.bool(True),
pileup_summary = cms.InputTag("slimmedAddPileupInfo"),
genevent = cms.InputTag("generator"),
rho = cms.InputTag("fixedGridRhoFastjetAll"),
jes = cms.untracked.string("nominal"),
jer = cms.untracked.string("nominal"),
lheleptoninfo = cms.untracked.bool(False)
)

process.p = cms.Path(process.cleanedMu*process.demo)



