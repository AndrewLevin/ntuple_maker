
import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

#process.options = cms.untracked.PSet(
#    numberOfStreams = cms.untracked.uint32( 0 ),
#   numberOfThreads = cms.untracked.uint32( 8 ),
#    wantSummary     = cms.untracked.bool( True )
#)

process.maxEvents = cms.untracked.PSet( 
input = cms.untracked.int32(-1) ,
#input = cms.untracked.int32(10) ,

)

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",

#skipEvents= cms.untracked.uint32(3146),

    # replace 'myfile.root' with the source file you want to use
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),

    fileNames = cms.untracked.vstring(

'/store/mc/RunIISummer16MiniAODv2/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/044AEFD6-92DE-E611-B842-0CC47A4D7616.root',

'/store/mc/RunIISummer16MiniAODv2/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/CA3039C2-9CDE-E611-9E74-0CC47A4D7662.root',

'/store/mc/RunIISummer16MiniAODv2/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/E40448C3-94DE-E611-9C31-0CC47A78A340.root'

#'/store/mc/RunIISummer16MiniAODv2/WLLJJ_WToLNu_MLL-4To60_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/0E8F00D6-DFEB-E611-99F6-001E67504475.root'

#'/store/mc/RunIISummer16MiniAODv2/WLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1A95882D-64C6-E611-9437-0CC47A4C8F30.root'

#'/store/mc/RunIISummer16MiniAODv2/WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/1ECE8AAB-E7DA-E611-A1D1-001E67E71B96.root'

#'/store/mc/RunIISummer16MiniAODv2/WLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/C27E1FD3-C8C5-E611-8E3C-0CC47A4D7692.root'

#'/store/mc/RunIISummer16MiniAODv2/WLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/B849F6A3-5BC5-E611-94A6-6C3BE5B59150.root'

#'/store/mc/RunIISummer16MiniAODv2/WLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/0A6FE368-F5C5-E611-9F7C-782BCB20E48C.root'

    ),

#eventsToProcess = cms.untracked.VEventRange('1:1747:239313-1:1747:239313'),

)

process.load("Configuration.Geometry.GeometryIdeal_cff") 
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")

#mc https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Global_Tags_for_Run2_MC_Producti
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'

from CondCore.DBCommon.CondDBSetup_cfi import *

###################################################################

#jecString = cms.string('sqlite:src/ntuple_maker/ntuple_maker/data/jec/Spring16_25nsV1_MC.db')

jecString = cms.string('sqlite:src/ntuple_maker/ntuple_maker/data/jec/Summer16_23Sep2016V4_MC.db')

process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK4PFchs'),
#            tag    = cms.string('JetCorrectorParametersCollection_Spring16_25nsV1_MC_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs')
            ),
      ## here you add as many jet types as you need
      ## note that the tag name is specific for the particular sqlite file 
      ), 
      connect = jecString
     # uncomment above tag lines and this comment to use MC JEC
)

process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring('L1FastJet',  'L2Relative', 'L3Absolute'), 'None')
)

process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

jerString = cms.string('sqlite:src/ntuple_maker/ntuple_maker/data/jer/Spring16_25nsV6_MC.db')

process.jer = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        toGet = cms.VPSet(
            # Resolution
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Spring16_25nsV6_MC_PtResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_pt')
                ),

            cms.PSet(
                 record = cms.string('JetResolutionRcd'),
                 tag    = cms.string('JR_Spring16_25nsV6_MC_PhiResolution_AK4PFchs'),
                 label  = cms.untracked.string('AK4PFchs_phi')
                 ),

            # Scale factors
            cms.PSet(
                record = cms.string('JetResolutionScaleFactorRcd'),
                tag    = cms.string('JR_Spring16_25nsV6_MC_SF_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs')
                ),
            ),
        connect = jerString
        )

process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')




process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
   src = cms.InputTag("slimmedMuons"),
   preselection = cms.string("track.isNonnull"),
   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
   fractionOfSharedSegments = cms.double(0.499))

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
           isData=False,

           )

process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')
from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)

process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')

process.demo = cms.EDAnalyzer('ntuple_maker',

  syscalcinfo = cms.untracked.bool (True), #fill the information from syscalc
  mgreweightinfo = cms.untracked.bool (False),
  vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#  muons = cms.InputTag("cleanedMu"),
  lheevent = cms.InputTag("externalLHEProducer"),
  lheruninfo = cms.InputTag("externalLHEProducer"),
#  lheruninfo = cms.InputTag("source"),
#  lheevent = cms.InputTag("source"),
  muons = cms.InputTag("slimmedMuons"),
  electrons = cms.InputTag("calibratedPatElectrons"),
#  electrons = cms.InputTag("slimmedElectrons"),
  taus = cms.InputTag("slimmedTaus"),
  photons = cms.InputTag("slimmedPhotons"),
jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
#  jets = cms.InputTag("slimmedJets"),
  fatjets = cms.InputTag("slimmedJetsAK8"),
#  mets = cms.InputTag("slimmedMETs"),
mets = cms.InputTag('slimmedMETs','','Demo'),

  isMC = cms.untracked.bool(True),
  prunedgenparticles = cms.InputTag("prunedGenParticles"),  
  packedgenparticles = cms.InputTag("packedGenParticles"),
  pfCands = cms.InputTag("packedPFCandidates"),
#  pileup_summary = cms.InputTag("addPileupInfo"),       
  pileup_summary = cms.InputTag("slimmedAddPileupInfo"),
  genevent = cms.InputTag("generator"),
  rho = cms.InputTag("fixedGridRhoFastjetAll"),
  rhoHLTElectronSelection = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
  #jes = cms.untracked.string("up"),
  #jes = cms.untracked.string("down"),
  jes = cms.untracked.string("nominal"),
#  jer = cms.untracked.string("up"),
  #jer = cms.untracked.string("down"),
#  jer = cms.untracked.string("nominal"),
  jer = cms.untracked.string("no_correction"),
  lheleptoninfo = cms.untracked.bool(False),
  apply_trigger = cms.untracked.bool(False),
  trigger_results_process = cms.untracked.string("HLT")                              
)

#process.p = cms.Path( process.regressionApplication * process.calibratedPatElectrons * process.fullPatMetSequence  * process.jecSequence*  process.cleanedMu*process.demo)
process.p = cms.Path( process.regressionApplication * process.calibratedPatElectrons * process.fullPatMetSequence  * process.jecSequence* process.demo)

