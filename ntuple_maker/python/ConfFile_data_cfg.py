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

#'/store/data/Run2016D/DoubleMuon/MINIAOD/23Sep2016-v1/90000/1EEEE4BB-0B8C-E611-B87F-3417EBE480D1.root'

#'/store/data/Run2016E/DoubleEG/MINIAOD/23Sep2016-v1/90000/1AF5373A-3C91-E611-AFC4-0CC47A4C8E8A.root'

#'/store/data/Run2016E/DoubleEG/MINIAOD/23Sep2016-v1/80000/5E6F8078-BE88-E611-9D50-0CC47A4D764C.root'

#'/store/data/Run2016B/DoubleEG/MINIAOD/23Sep2016-v3/110000/408DD8D7-3398-E611-BDC5-02163E014AD3.root'

#'/store/data/Run2016B/DoubleEG/MINIAOD/23Sep2016-v3/120000/18E0EC8C-3598-E611-BF53-0CC47A78A3B4.root'

#'/store/data/Run2016B/DoubleEG/MINIAOD/23Sep2016-v3/00000/84FFB190-E697-E611-AB0A-003048F5E840.root'

#'/store/data/Run2016B/DoubleEG/MINIAOD/23Sep2016-v3/120000/18E0EC8C-3598-E611-BF53-0CC47A78A3B4.root'

#'/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v3/60000/C81A53E1-8598-E611-9D91-F45214939730.root'

#'/store/data/Run2016B/SingleElectron/MINIAOD/23Sep2016-v3/60000/F0FEEBFC-959A-E611-8E9B-0025905AA9CC.root'

#'/store/data/Run2016B/DoubleEG/MINIAOD/23Sep2016-v3/110000/3E3AEF51-4898-E611-85CA-0025905A6084.root'

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/00000/CC154516-589B-E611-A533-008CFA197E0C.root'

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/00000/00269AA1-479B-E611-A359-0025905A6082.root'

'/store/data/Run2016H/MuonEG/MINIAOD/PromptReco-v2/000/283/876/00000/AA7CCD56-C19C-E611-8F3B-02163E011872.root'

),

#eventsToProcess = cms.untracked.VEventRange('276525:394:592892349-276525:394:592892349'),

)

#process.load("Configuration.Geometry.GeometryIdeal_cff") 
#process.load('Configuration.StandardSequences.Services_cff')
#process.load("Configuration.StandardSequences.MagneticField_cff")

#mc https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Global_Tags_for_Run2_MC_Producti
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

#process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v3'

from CondCore.DBCommon.CondDBSetup_cfi import *

#
jecString = cms.string('sqlite:src/ntuple_maker/ntuple_maker/data/jec/Spring16_25nsV10All_DATA.db')

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
            tag    = cms.string('JetCorrectorParametersCollection_Spring16_25nsV10All_DATA_AK4PFchs'),
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

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
           isData=True,
                           
           )

process.demo = cms.EDAnalyzer('ntuple_maker',

  syscalcinfo = cms.untracked.bool (False), #fill the information from syscalc
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
lheleptoninfo = cms.untracked.bool(False)
)

process.p = cms.Path(process.fullPatMetSequence*process.jecSequence*  process.cleanedMu*process.demo)



