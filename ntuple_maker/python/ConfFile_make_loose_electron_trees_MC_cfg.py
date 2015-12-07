import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

'/store/mc/RunIISpring15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/008865AA-596D-E511-92F1-0025905A6110.root'

#'/store/data/Run2015D/DoubleMuon/MINIAOD/PromptReco-v4/000/260/627/00000/6294465C-2285-E511-A6F9-02163E014431.root',
#'/store/data/Run2015D/DoubleMuon/MINIAOD/PromptReco-v4/000/260/627/00000/648DE552-1E85-E511-944F-02163E0143CF.root',
#'/store/data/Run2015D/DoubleMuon/MINIAOD/PromptReco-v4/000/260/627/00000/66962EAE-1E85-E511-8F42-02163E01472E.root',
#'/store/data/Run2015D/DoubleMuon/MINIAOD/PromptReco-v4/000/260/627/00000/6CA5132D-1E85-E511-ABA0-02163E0118CC.root',
#'/store/data/Run2015D/DoubleMuon/MINIAOD/PromptReco-v4/000/260/627/00000/703A0D2C-1E85-E511-A827-02163E0118CC.root',

#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v3/000/257/400/00000/587D2C79-4C65-E511-A018-02163E0142E9.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v3/000/257/400/00000/64CFDD24-5065-E511-9BAB-02163E011B22.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v3/000/257/400/00000/7A458316-5065-E511-8FFC-02163E01421A.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v3/000/257/400/00000/84057F14-5065-E511-BE28-02163E01412F.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v3/000/257/400/00000/94B0A799-5065-E511-969F-02163E0133BA.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v3/000/257/400/00000/98063C1D-4D65-E511-95E1-02163E011A31.root'


#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v3/000/256/584/00000/D4DF072A-855D-E511-B91B-02163E01468C.root',
#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v3/000/256/587/00000/560F124D-925D-E511-BC70-02163E011D99.root',


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
  isMC = cms.untracked.bool(True),
  lepton_flavor = cms.untracked.string("electron"),
  rho = cms.InputTag("fixedGridRhoFastjetAll"),
  lheevent = cms.InputTag("externalLHEProducer"),
  which_triggers = cms.untracked.string("electron_fake_rate"),
  genevent = cms.InputTag("generator"),

)


process.p = cms.Path(process.demo)
