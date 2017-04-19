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
#input = cms.untracked.int32(10000) ,

)

#process.load('EgammaAnalysis/ElectronTools/egmGsfElectronIDs_cff')

process.source = cms.Source("PoolSource",

#skipEvents= cms.untracked.uint32(1),

    # replace 'myfile.root' with the source file you want to use
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),

    fileNames = cms.untracked.vstring(

#'/store/mc/RunIISummer16MiniAODv2/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/044AEFD6-92DE-E611-B842-0CC47A4D7616.root',

#'/store/mc/RunIISummer16MiniAODv2/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/CA3039C2-9CDE-E611-9E74-0CC47A4D7662.root',

#'/store/mc/RunIISummer16MiniAODv2/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/E40448C3-94DE-E611-9C31-0CC47A78A340.root'

#'/store/mc/RunIISummer16MiniAODv2/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0008345B-C0C9-E611-91E6-A4BADB1E763D.root'

#'/store/mc/RunIISummer16MiniAODv2/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0023AF2C-D7CD-E611-9247-002590E7D7CE.root'

#'/store/mc/RunIISummer16MiniAODv2/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0030B9D6-72C1-E611-AE49-02163E00E602.root'

'/store/mc/RunIISummer16MiniAODv2/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0030B9D6-72C1-E611-AE49-02163E00E602.root'

#'/store/mc/RunIISpring16MiniAODv2/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/2E6533A7-F225-E611-9C92-0025905B85BA.root'

#'/store/mc/RunIISpring16MiniAODv2/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/CEB823ED-E224-E611-85CC-3417EBE64444.root'

#'/store/mc/RunIIFall15MiniAODv2/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/E89EB2E2-A7B8-E511-B30A-001E67A402C1.root'

#'/store/mc/RunIIFall15MiniAODv2/WZJJ_QCD_13TeV-madgraph-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/48811764-39B8-E511-BCC5-F46D043B3CE5.root'

# '/store/mc/RunIIFall15MiniAODv2/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/10000/044F8A3A-43B8-E511-8F98-0025904CF75A.root'

#'/store/mc/RunIIFall15MiniAODv2/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/20000/420FA4F7-1CB9-E511-948E-02163E013C5C.root'

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
#''

#'/store/mc/RunIIFall15MiniAODv2/WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/C69294D5-4CB8-E511-9FF5-001E67E95C40.root'

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

  syscalcinfo = cms.untracked.bool (False), #fill the information from syscalc
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

#process.p = cms.Path(process.fullPatMetSequence*process.jecSequence*  process.cleanedMu*process.demo)
process.p = cms.Path(process.fullPatMetSequence*process.jecSequence* process.demo)

