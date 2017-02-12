import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")


process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("output_tree.root") )

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/00000/00269AA1-479B-E611-A359-0025905A6082.root'

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/60000/BC46B8EE-599B-E611-81CE-0025905A48B2.root'

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/60000/7E7E5DCD-4F9B-E611-A66B-0CC47A4D7602.root'

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/60000/3EBD40B4-599B-E611-B916-003048F5B2A0.root'

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/60000/B879E731-9A9B-E611-B10E-0242AC130004.root'

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/120000/76A06DDA-A99B-E611-9BBB-0025905A60AA.root'

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/60000/8A706036-909B-E611-B52E-0025905B8598.root'

#        '/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/60000/A463B263-6D9B-E611-B2C7-0CC47A4C8E56.root',

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/60000/5A95AF67-619B-E611-B8D7-0025905B8592.root',

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/60000/D66EC9FD-4C9B-E611-A060-0CC47A78A4B0.root',

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/60000/EC00F87E-459B-E611-AC94-0025905B856C.root'

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/00000/7A7B15D8-2F9B-E611-B470-0025905B85B2.root'

#'/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v3/120000/B8B1D53E-9C98-E611-804D-FA163EBA41B0.root'

#'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/00000/D64161BC-3C9B-E611-83D8-B083FED4239C.root'

#'/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v3/00000/C82F4E71-D297-E611-96F8-0025907B4EE2.root'

#'/store/data/Run2016B/DoubleEG/MINIAOD/23Sep2016-v3/110000/549ADDCD-5D98-E611-82C0-FA163EC3811B.root'

'/store/data/Run2016B/MuonEG/MINIAOD/23Sep2016-v3/120000/B66DCFF5-759A-E611-9379-0025907B4EC8.root'

    ),

eventsToProcess = cms.untracked.VEventRange('275073:286:516519502-275073:286:516519502'),

)

from CondCore.DBCommon.CondDBSetup_cfi import *

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

            # Phi
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


from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
           isData=True,
           )

process.demo = cms.EDAnalyzer('make_loose_lepton_trees',
  vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  muons = cms.InputTag("slimmedMuons"),
  electrons = cms.InputTag("slimmedElectrons"),
  taus = cms.InputTag("slimmedTaus"),
  photons = cms.InputTag("slimmedPhotons"),
  jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
#  jets = cms.InputTag("slimmedJets"),
  fatjets = cms.InputTag("slimmedJetsAK8"),
  mets = cms.InputTag("slimmedMETs",'',"Demo"),  
#mets = cms.InputTag('slimmedMETs'),

  prunedgenparticles = cms.InputTag("prunedGenParticles"),  
  packedgenparticles = cms.InputTag("packedGenParticles"),  
  isMC = cms.untracked.bool(False),
  lepton_flavor = cms.untracked.string("muon"),
  rho = cms.InputTag("fixedGridRhoFastjetAll"),
  rhoHLTElectronSelection = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
  genevent = cms.InputTag("generator"),
  lheevent = cms.InputTag("externalLHEProducer"),
  which_triggers = cms.untracked.string('muon_fake_rate')

)

process.p = cms.Path(process.fullPatMetSequence *process.jecSequence* process.demo)
#process.p = cms.Path(process.fullPatMetSequence* process.demo)
#process.p = cms.Path(process.demo)
