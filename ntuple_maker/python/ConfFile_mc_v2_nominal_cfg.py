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

#skipEvents= cms.untracked.uint32(1),

    # replace 'myfile.root' with the source file you want to use
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),

    fileNames = cms.untracked.vstring(

#'/store/mc/RunIIFall15MiniAODv2/WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/C69294D5-4CB8-E511-9FF5-001E67E95C40.root'

#'/store/mc/RunIISummer16MiniAODv2/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/044AEFD6-92DE-E611-B842-0CC47A4D7616.root'

#'/store/mc/RunIISpring16MiniAODv2/WmWmJJ_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/00000/1C2655B1-5A38-E611-A3C0-C4346BC8C638.root'

#'/store/mc/RunIISpring16MiniAODv2/WLLJJToLNu_M-60_EWK_13TeV-madgraph-pythia8/MINIAODSIM/premix_withHLT_80X_mcRun2_asymptotic_v14-v1/80000/064D833C-537F-E611-987E-20CF3027A5BC.root'

#'/store/mc/RunIISpring16MiniAODv2/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/08277AB0-F225-E611-8046-D4AE526A33F5.root'

#'/store/mc/RunIISpring16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/20000/02652A75-C654-E611-9734-0CC47A4D7686.root'

'/store/mc/RunIISpring16MiniAODv2/WLLJJToLNu_M-60_EWK_13TeV-madgraph-pythia8/MINIAODSIM/premix_withHLT_80X_mcRun2_asymptotic_v14-v1/80000/064D833C-537F-E611-987E-20CF3027A5BC.root'

    ),

#eventsToProcess = cms.untracked.VEventRange('1:985:121697-1:985:121697'),

)

process.load("Configuration.Geometry.GeometryIdeal_cff") 
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")

#mc https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Global_Tags_for_Run2_MC_Producti
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'

from CondCore.DBCommon.CondDBSetup_cfi import *

###################################################################

#jecString = cms.string('sqlite:src/ntuple_maker/ntuple_maker/data/jec/Spring16_25nsV1_MC.db')

jecString = cms.string('sqlite:src/ntuple_maker/ntuple_maker/data/jec/Spring16_23Sep2016V1_MC.db')

process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Spring16_23Sep2016V1_MC_AK4PFchs'),
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
  electrons = cms.InputTag("slimmedElectrons"),
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

process.p = cms.Path(process.fullPatMetSequence*process.jecSequence*  process.cleanedMu*process.demo)

