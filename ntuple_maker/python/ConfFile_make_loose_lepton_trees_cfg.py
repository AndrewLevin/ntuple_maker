import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/790/00000/F0EC5CBE-3F4A-E511-A3F1-02163E011D72.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/790/00000/F82645C6-3F4A-E511-9E60-02163E011CF7.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/824/00000/D2E7EDBC-504A-E511-97C4-02163E014288.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/0EFC4F55-FB4A-E511-BFDE-02163E0141FD.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/1C007744-FB4A-E511-A22C-02163E0136A2.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/22CC8856-FB4A-E511-AB3D-02163E011F63.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/30B68949-FB4A-E511-BE6C-02163E0118A6.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/3C506B4E-FB4A-E511-8996-02163E0118A1.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/4016C44E-FB4A-E511-BC40-02163E011E62.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/58B60B54-FB4A-E511-BC56-02163E011FE6.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/5A911356-FB4A-E511-AA28-02163E0119FB.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/70B2C251-FB4A-E511-BCC8-02163E012687.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/720E4A53-FB4A-E511-8348-02163E014145.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/727A7B0C-044B-E511-948B-02163E012AAA.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/74E2404A-FB4A-E511-87CE-02163E01412A.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/76138555-FB4A-E511-9934-02163E013463.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/80018D50-FB4A-E511-9262-02163E011DC6.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/82768752-FB4A-E511-94CD-02163E013680.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/86BDC843-FB4A-E511-81D1-02163E0141C7.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/888EFC55-FB4A-E511-AB0E-02163E014566.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/90181353-FB4A-E511-92D7-02163E011CB3.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/9063FA52-FB4A-E511-B115-02163E01454C.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/9AF7284F-FB4A-E511-83F0-02163E012144.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/AA50BB53-FB4A-E511-B028-02163E0137AA.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/AE9CAE47-FB4A-E511-B53B-02163E0141D7.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/BA97E25C-FB4A-E511-AF0A-02163E011DE0.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/BEBC064C-FB4A-E511-9ED0-02163E014237.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/C45CCE47-FB4A-E511-9AEB-02163E014459.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/D648B648-FB4A-E511-8B88-02163E01438E.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/DC61DCBA-FD4A-E511-91EE-02163E0138C2.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/DCE3F751-FB4A-E511-B01B-02163E011CDF.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/E2B2DD64-FB4A-E511-A7A5-02163E013776.root',
'/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/833/00000/F05FA947-FB4A-E511-8A3A-02163E013408.root',

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
)


process.p = cms.Path(process.demo)
