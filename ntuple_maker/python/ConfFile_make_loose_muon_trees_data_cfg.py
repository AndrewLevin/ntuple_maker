import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/00000/00269AA1-479B-E611-A359-0025905A6082.root'

    )
)

process.demo = cms.EDAnalyzer('make_loose_lepton_trees',
  vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  muons = cms.InputTag("slimmedMuons"),
  electrons = cms.InputTag("slimmedElectrons"),
  taus = cms.InputTag("slimmedTaus"),
  photons = cms.InputTag("slimmedPhotons"),
  jets = cms.InputTag("slimmedJets"),
  fatjets = cms.InputTag("slimmedJetsAK8"),
  mets = cms.InputTag("slimmedMETs"),  
  prunedgenparticles = cms.InputTag("prunedGenParticles"),  
  packedgenparticles = cms.InputTag("packedGenParticles"),  
  isMC = cms.untracked.bool(False),
  lepton_flavor = cms.untracked.string("muon"),
  rho = cms.InputTag("fixedGridRhoFastjetAll"),
rhoHLTElectronSelection = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
  genevent = cms.InputTag("generator"),
  lheevent = cms.InputTag("externalLHEProducer"),
  which_triggers = cms.untracked.string('muon_fake_rate')

)


process.p = cms.Path(process.demo)
