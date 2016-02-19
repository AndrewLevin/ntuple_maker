import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

process.maxEvents = cms.untracked.PSet( 
#input = cms.untracked.int32(1) ,
input = cms.untracked.int32(-1) ,

)

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",

#skipEvents= cms.untracked.uint32(1),

    # replace 'myfile.root' with the source file you want to use
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),

    fileNames = cms.untracked.vstring(

# '/store/mc/RunIIFall15MiniAODv2/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/10000/044F8A3A-43B8-E511-8F98-0025904CF75A.root'

'/store/mc/RunIIFall15MiniAODv2/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/E89EB2E2-A7B8-E511-B30A-001E67A402C1.root'

#'file:eos/cms/store/user/anlevin/data/MINIAOD/wgjj_ewk_qcd_v1/step4_output_1.root',
#'file:eos/cms/store/user/anlevin/data/MINIAOD/wgjj_ewk_qcd_v1/step4_output_101.root',
#'file:eos/cms/store/user/anlevin/data/MINIAOD/wgjj_ewk_qcd_v1/step4_output_201.root',
#'file:eos/cms/store/user/anlevin/data/MINIAOD/wgjj_ewk_qcd_v1/step4_output_301.root',
#'file:eos/cms/store/user/anlevin/data/MINIAOD/wgjj_ewk_qcd_v1/step4_output_401.root',
#'file:eos/cms/store/user/anlevin/data/MINIAOD/wgjj_ewk_qcd_v1/step4_output_501.root',
#'file:eos/cms/store/user/anlevin/data/MINIAOD/wgjj_ewk_qcd_v1/step4_output_601.root',
#'file:eos/cms/store/user/anlevin/data/MINIAOD/wgjj_ewk_qcd_v1/step4_output_701.root',
#'file:eos/cms/store/user/anlevin/data/MINIAOD/wgjj_ewk_qcd_v1/step4_output_801.root',
#'file:eos/cms/store/user/anlevin/data/MINIAOD/wgjj_ewk_qcd_v1/step4_output_901.root',


#'/store/mc/RunIIFall15MiniAODv1/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/006EF3A3-34B1-E511-B8FD-B083FED42FC4.root'

#'/store/mc/RunIISpring15MiniAODv2/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/80000/88743C5E-BE7C-E511-8568-0025905A48B2.root'

#'file:/afs/cern.ch/work/a/anlevin/VBS/13_tev/Merged.root'

#'file:/afs/cern.ch/work/a/anlevin/tmp/306D0BE6-F57A-E511-BF8B-003048FF9AA6.root',
#'file:/afs/cern.ch/work/a/anlevin/tmp/C40BEFD4-F57A-E511-8427-0025905B85B2.root'

#'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/00759690-D16E-E511-B29E-00261894382D.root'

#'/store/mc/RunIISpring15DR74/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/60000/42CA2568-3D78-E511-9982-A01D48EE83AE.root'
#'/store/mc/RunIISpring15DR74/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/60000/AA453EED-1F79-E511-9AF3-001E67A4069F.root'
#'/store/mc/RunIISpring15DR74/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/80000/306D0BE6-F57A-E511-BF8B-003048FF9AA6.root'
#'/store/mc/RunIISpring15DR74/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/80000/66238FE4-F57A-E511-B1CA-0025905964B2.root'
#'/store/mc/RunIISpring15DR74/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/80000/C40BEFD4-F57A-E511-8427-0025905B85B2.root'

#'/store/mc/RunIISpring15DR74/WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/26DBC5F1-8566-E511-8C82-20CF3027A607.root',
#'/store/mc/RunIISpring15DR74/WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/3A5140F2-8566-E511-8242-B083FED42488.root',
#'/store/mc/RunIISpring15DR74/WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/4818A9C3-8566-E511-9BFB-008CFA11139C.root',
#'/store/mc/RunIISpring15DR74/WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/48249135-8966-E511-A11D-0002C92A1030.root'


#'/store/mc/RunIISpring15DR74/WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/60A6DBF4-8566-E511-9DA1-002618943975.root'

#'/store/mc/RunIISpring15DR74/WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/102AAAF1-0966-E511-B871-24BE05C6D5B1.root'

#'/store/mc/RunIISpring15DR74/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/001E8D25-8112-E511-A426-001517E74088.root'

    ),
#eventsToProcess = cms.untracked.VEventRange('1:2437:158389-1:2437:158389'),

)

process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
   src = cms.InputTag("slimmedMuons"),
   preselection = cms.string("track.isNonnull"),
   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
   fractionOfSharedSegments = cms.double(0.499))

process.demo = cms.EDAnalyzer('ntuple_maker',

  syscalcinfo = cms.untracked.bool (True), #fill the information from syscalc
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
#  pileup_summary = cms.InputTag("addPileupInfo"),                              
  pileup_summary = cms.InputTag("slimmedAddPileupInfo"),
  genevent = cms.InputTag("generator"),
  rho = cms.InputTag("fixedGridRhoFastjetAll")
)

process.p = cms.Path(process.cleanedMu*process.demo)
