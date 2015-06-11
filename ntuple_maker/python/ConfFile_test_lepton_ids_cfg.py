import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
#'file:/afs/cern.ch/work/a/anlevin/tmp/MySkim_1.electrons.root'
#'file:/afs/cern.ch/work/a/anlevin/delete_this_5/CMSSW_7_2_0/src/DiJet.root'
        '/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root'
#        '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU40bx25_POSTLS170_V5-v2/00000/00CAA728-D626-E411-9112-00215AD4D6E2.root'
#        '/store/user/anlevin/data/MINIAOD/wpwp_13_tev_qed_4_qcd_0_v3/step5_output_9601.root'
#        'file:/afs/cern.ch/work/a/anlevin/VBS/13_tev/Merged.root'
#'file:/afs/cern.ch/work/a/anlevin/tmp/step5_output_4.root'


    )
)

process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
   src = cms.InputTag("slimmedMuons"),
   preselection = cms.string("track.isNonnull"),
   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
   fractionOfSharedSegments = cms.double(0.499))

process.demo = cms.EDAnalyzer('test_lepton_ids',
  vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  muons = cms.InputTag("cleanedMu"),
  #muons = cms.InputTag("slimmedMuons"),
  electrons = cms.InputTag("slimmedElectrons"),
  taus = cms.InputTag("slimmedTaus"),
  photons = cms.InputTag("slimmedPhotons"),
  jets = cms.InputTag("slimmedJets"),
  fatjets = cms.InputTag("slimmedJetsAK8"),
  mets = cms.InputTag("slimmedMETs"),  
  prunedgenparticles = cms.InputTag("prunedGenParticles"),  
  packedgenparticles = cms.InputTag("packedGenParticles"),  
                              
)

process.p = cms.Path(process.cleanedMu*process.demo)


