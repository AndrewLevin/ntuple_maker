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

'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/457/00000/4CF346F1-9646-E511-94B6-02163E0143CC.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/458/00000/5074B052-9C46-E511-AECC-02163E0128D5.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/459/00000/3606102A-5B47-E511-8A32-02163E012AAF.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/459/00000/4A56412B-5B47-E511-8BE3-02163E012B2D.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/459/00000/5EFCB82B-5B47-E511-8739-02163E015603.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/459/00000/B265E62C-5B47-E511-ABAC-02163E01350F.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/459/00000/B8FAD227-5B47-E511-AB7D-02163E0140F0.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/459/00000/F8864D27-5B47-E511-9E5A-02163E014222.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/500/00000/A293A84A-D646-E511-B8FD-02163E01339E.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/511/00000/6AFB99EB-D246-E511-94F3-02163E0146E8.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/512/00000/4A22A5C7-2547-E511-B34C-02163E014648.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/513/00000/4A85DFEC-D146-E511-B236-02163E0144AE.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/530/00000/4C97D391-0B47-E511-A99E-02163E014434.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/532/00000/002755AE-6F47-E511-9874-02163E015529.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/532/00000/4204CEC0-6F47-E511-89E5-02163E01448B.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/602/00000/32BFEAE1-A447-E511-8787-02163E0143DD.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/607/00000/B6695F05-BB47-E511-AA80-02163E0139A4.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/608/00000/5C2C4A19-0249-E511-A008-02163E014414.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/608/00000/AEF0521D-0249-E511-8A61-02163E01559C.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/608/00000/B484B020-0249-E511-878A-02163E01343B.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/608/00000/B4FCAC1A-0249-E511-9A4E-02163E011DB9.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/608/00000/FE58F71A-0249-E511-BBDC-02163E011E27.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/712/00000/3CBD8F3C-CE48-E511-B616-02163E011B1D.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/768/00000/3EA24229-1049-E511-AFB8-02163E0141C4.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/769/00000/40FFD4DE-1849-E511-AFD0-02163E011AB8.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/780/00000/CE2BECA1-2949-E511-AC4D-02163E014738.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/785/00000/C82351CD-3D49-E511-834E-02163E014414.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/790/00000/22890CFA-034A-E511-A23D-02163E011EC1.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/790/00000/52A6F3F2-034A-E511-AE67-02163E014204.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/790/00000/B6F44AEE-034A-E511-95DF-02163E014176.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/790/00000/BA9FD8E8-034A-E511-BB81-02163E0141C7.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/790/00000/C29C6EFD-034A-E511-8A4A-02163E015603.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/790/00000/DEDA90EF-034A-E511-9F59-02163E014589.root',
'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/824/00000/A4E2DEAF-364A-E511-91A9-02163E011A47.root',

#/store/mc/RunIISpring15DR74/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/003F747A-D835-E511-BFC8-002590D9D976.root'

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

  syscalcinfo = cms.untracked.bool (False), #fill the information from syscalc
  lheinfo = cms.untracked.bool (False),
  vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  muons = cms.InputTag("cleanedMu"),
#  lheevent = cms.InputTag("externalLHEProducer"),
#  lheruninfo = cms.InputTag("externalLHEProducer"),
  lheruninfo = cms.InputTag("source"),
  lheevent = cms.InputTag("source"),
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
                              
)

process.p = cms.Path(process.cleanedMu*process.demo)



