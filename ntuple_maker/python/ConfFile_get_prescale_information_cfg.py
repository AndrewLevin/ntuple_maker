import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

process.maxEvents = cms.untracked.PSet( 
input = cms.untracked.int32(10000) ,
#input = cms.untracked.int32(-1) ,

)

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",

#skipEvents= cms.untracked.uint32(1),

    # replace 'myfile.root' with the source file you want to use
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),

    fileNames = cms.untracked.vstring(
'/store/data/Run2015D/ZeroBias/MINIAOD/PromptReco-v3/000/258/136/00000/526CC2AD-FC6A-E511-8E65-02163E0145AC.root'

#'/store/data/Run2015D/ZeroBias/MINIAOD/PromptReco-v4/000/258/159/00000/040923F6-D16B-E511-A319-02163E014780.root'

#'/store/data/Run2015B/ZeroBias/MINIAOD/05Oct2015-v1/30000/00FE726D-8A6E-E511-822E-0025905A6084.root'

#' /store/data/Run2015D/ZeroBias/MINIAOD/PromptReco-v3/000/256/584/00000/4AB13C16-865D-E511-A9DA-02163E011A40.root'

#'/store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v4/000/258/159/00000/027612B0-306C-E511-BD47-02163E014496.root'

#'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/20FF6359-AD70-E511-9CDF-003048FFCBB0.root',
#'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/22F3621E-B770-E511-8E67-002590596498.root',
#'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/22FBD2F7-DD70-E511-982B-0025905A48F2.root',
#'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/24723F7A-DF70-E511-A2F6-002618943954.root',
#'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/248920C4-AD6F-E511-A133-6CC2173BB810.root',
#'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/248C8B7D-AD70-E511-BE6C-0025905A60DA.root'

#'/store/mc/RunIISpring15DR74/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/001E8D25-8112-E511-A426-001517E74088.root',
#'/store/mc/RunIISpring15DR74/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/18E763F0-8012-E511-BA75-0002C90A3456.root',
#'/store/mc/RunIISpring15DR74/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/72306A00-8112-E511-85EA-485B3919F121.root',
#'/store/mc/RunIISpring15DR74/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/B001A6FB-8012-E511-B8AA-001E675A58D9.root'

#'/store/mc/RunIISpring15DR74/WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/60A6DBF4-8566-E511-9DA1-002618943975.root'

#'/store/mc/RunIISpring15DR74/WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/102AAAF1-0966-E511-B871-24BE05C6D5B1.root'

#'/store/mc/RunIISpring15DR74/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/001E8D25-8112-E511-A426-001517E74088.root'

    ),
#eventsToProcess = cms.untracked.VEventRange('1:2437:158389-1:2437:158389'),

)

process.demo = cms.EDAnalyzer('get_prescale_information',

prescales = cms.InputTag("patTrigger"),
l1min = cms.InputTag("patTrigger","l1min"),
l1max = cms.InputTag("patTrigger","l1max")

)

process.p = cms.Path(process.demo)



