import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

process.maxEvents = cms.untracked.PSet( 
#input = cms.untracked.int32(100) ,
input = cms.untracked.int32(-1) ,

)

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",

#skipEvents= cms.untracked.uint32(1),

    # replace 'myfile.root' with the source file you want to use
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),

    fileNames = cms.untracked.vstring(

'/store/data/Run2016H/MuonEG/MINIAOD/03Feb2017_ver3-v1/110000/88C33F8D-96EA-E611-9767-0025905B85FC.root'

#'/store/data/Run2016H/MuonEG/MINIAOD/PromptReco-v2/000/283/876/00000/AA7CCD56-C19C-E611-8F3B-02163E011872.root'

#'/store/data/Run2016G/DoubleMuon/MINIAOD/03Feb2017-v1/100000/00182C13-EEEA-E611-8897-001E675A6C2A.root'

#'/store/data/Run2016H/DoubleMuon/MINIAOD/03Feb2017_ver2-v1/110000/E6E5089B-C4EC-E611-BC44-549F3525C0BC.root'




#'/store/data/Run2016G/DoubleEG/MINIAOD/03Feb2017-v1/100000/7EC857A9-C1EA-E611-982E-0025905A60F4.root',

#'/store/data/Run2016G/DoubleEG/MINIAOD/03Feb2017-v1/100000/0AC03D12-D1EA-E611-A8A1-008CFA197448.root',

#'/store/data/Run2016G/DoubleEG/MINIAOD/03Feb2017-v1/100000/9E7A1BFF-D0EA-E611-BCF5-0CC47AD99176.root',

#'/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50000/4E70C3B0-21EB-E611-9F02-90B11C050371.root',

#'/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50000/B02F2A53-15EB-E611-8D90-002590E39D8A.root',

#'/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50000/EEE22B96-1BEB-E611-A204-90B11C2ACF20.root',

#'/store/data/Run2016E/SingleElectron/MINIAOD/03Feb2017-v1/110000/90776164-7DEA-E611-8C40-0CC47AD9914A.root',

#'/store/data/Run2016B/SingleElectron/MINIAOD/03Feb2017_ver2-v2/110001/26EEF3AF-FAEB-E611-BA16-0CC47ADAF60A.root'

),

#eventsToProcess = cms.untracked.VEventRange('284043:89:139875961-284043:89:139875961'),

)

process.load("Configuration.Geometry.GeometryIdeal_cff") 
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")

#mc https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Global_Tags_for_Run2_MC_Producti
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'

from CondCore.DBCommon.CondDBSetup_cfi import *

#
jecString = cms.string('sqlite:src/ntuple_maker/ntuple_maker/data/jec/Summer16_23Sep2016AllV4_DATA.db')

#jecString = cms.string('sqlite:src/ntuple_maker/ntuple_maker/data/Spring16_23Sep2016AllV1_DATA.db')

process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
#            tag    = cms.string('JetCorrectorParametersCollection_Spring16_23Sep2016AllV1_DATA_AK4PFchs'),
            tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016AllV4_DATA_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs')
            ),
      ## here you add as many jet types as you need
      ## note that the tag name is specific for the particular sqlite file 
      ), 
      connect = jecString
     # uncomment above tag lines and this comment to use MC JEC
)

process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

#see here https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution and https://github.com/cms-jet/JRDatabase/tree/master/SQLiteFiles and https://twiki.cern.ch/twiki/bin/view/CMS/JERCReference

jerString = cms.string('sqlite:src/ntuple_maker/ntuple_maker/data/jer/Spring16_25nsV6_DATA.db')

#this is needed for met recalculation
process.jer = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        toGet = cms.VPSet(
            # Resolution
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Spring16_25nsV6_DATA_PtResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_pt')
                ),

            cms.PSet(
                 record = cms.string('JetResolutionRcd'),
                 tag    = cms.string('JR_Spring16_25nsV6_DATA_PhiResolution_AK4PFchs'),
                 label  = cms.untracked.string('AK4PFchs_phi')
                 ),     

            # Scale factors
            cms.PSet(
                record = cms.string('JetResolutionScaleFactorRcd'),
                tag    = cms.string('JR_Spring16_25nsV6_DATA_SF_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs')
                ),
            ),
        connect = jerString
        )

process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring('L1FastJet',  'L2Relative', 'L3Absolute', 'L2L3Residual'), 'None') 
)

process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)


process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
   src = cms.InputTag("slimmedMuons"),
   preselection = cms.string("track.isNonnull"),
   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
   fractionOfSharedSegments = cms.double(0.499))

#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#runMetCorAndUncFromMiniAOD(process,
#           isData=True,
                           
#           )

process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')
from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)

process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process, isData=True)

process.demo = cms.EDAnalyzer('ntuple_maker',

  syscalcinfo = cms.untracked.bool (False), #fill the information from syscalc
  mgreweightinfo = cms.untracked.bool (False),
  vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#  muons = cms.InputTag("cleanedMu"),
  muons = cms.InputTag("slimmedMuons"),
  lheevent = cms.InputTag("externalLHEProducer"),
  lheruninfo = cms.InputTag("externalLHEProducer"),
#  lheruninfo = cms.InputTag("source"),
#  lheevent = cms.InputTag("source"),
  #muons = cms.InputTag("slimmedMuons"),
#  electrons = cms.InputTag("slimmedElectrons"),
  electrons = cms.InputTag("calibratedPatElectrons"),
  taus = cms.InputTag("slimmedTaus"),
  photons = cms.InputTag("slimmedPhotons"),
#  jets = cms.InputTag("slimmedJets"),
jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
  fatjets = cms.InputTag("slimmedJetsAK8"),
#  mets = cms.InputTag("slimmedMETs"),  
 mets = cms.InputTag('slimmedMETs','','Demo'),

  isMC = cms.untracked.bool(False),  
  prunedgenparticles = cms.InputTag("prunedGenParticles"),  
  packedgenparticles = cms.InputTag("packedGenParticles"),  
  pfCands = cms.InputTag("packedPFCandidates"),                              
apply_trigger = cms.untracked.bool(True),
pileup_summary = cms.InputTag("slimmedAddPileupInfo"),
genevent = cms.InputTag("generator"),
rho = cms.InputTag("fixedGridRhoFastjetAll"),
rhoHLTElectronSelection = cms.InputTag("fixedGridRhoFastjetCentralCalo"),

jes = cms.untracked.string("nominal"),
jer = cms.untracked.string("nominal"),
lheleptoninfo = cms.untracked.bool(False),
trigger_results_process = cms.untracked.string("HLT")
)


#process.p = cms.Path( process.regressionApplication * process.calibratedPatElectrons * process.fullPatMetSequence  * process.jecSequence*  process.cleanedMu*process.demo)
process.p = cms.Path( process.regressionApplication * process.calibratedPatElectrons * process.fullPatMetSequence  * process.jecSequence*  process.demo)


