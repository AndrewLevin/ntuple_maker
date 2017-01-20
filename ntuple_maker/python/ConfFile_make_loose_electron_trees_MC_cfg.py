import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

#'/store/mc/RunIISpring15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/008865AA-596D-E511-92F1-0025905A6110.root'

#'/store/mc/RunIISpring16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/127DD5B1-DC1A-E611-93AE-0025905AC824.root'

#'/store/mc/RunIISpring16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/20000/02652A75-C654-E611-9734-0CC47A4D7686.root'

#'/store/mc/RunIISpring16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/60000/DACF1211-B756-E611-AD9C-0CC47A4D76AC.root'

'/store/mc/RunIISpring16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/60000/E8769581-5655-E611-8AF9-14DDA9D4F20C.root'

    ),

#eventsToProcess = cms.untracked.VEventRange('1:4605:3206531-1:4605:3206531'),



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
  isMC = cms.untracked.bool(True),
  lepton_flavor = cms.untracked.string("electron"),
  rho = cms.InputTag("fixedGridRhoFastjetAll"),
rhoHLTElectronSelection = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
  lheevent = cms.InputTag("externalLHEProducer"),
  which_triggers = cms.untracked.string("electron_fake_rate"),
  genevent = cms.InputTag("generator"),

)


process.p = cms.Path(process.demo)
