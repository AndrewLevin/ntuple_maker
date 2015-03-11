# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: skim --step SKIM:DiJet --conditions auto --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('SKIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Skims_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root'),
    eventsToProcess = cms.untracked.VEventRange(
"1:10634:1063353-1:10634:1063353",
"1:203643:20364295-1:203643:20364295",
"1:48267:4826612-1:48267:4826612",
"1:398901:39890067-1:398901:39890067",
"1:269044:26904362-1:269044:26904362",
"1:245340:24533918-1:245340:24533918",
"1:559615:55961437-1:559615:55961437",
"1:169175:16917437-1:169175:16917437",
"1:218871:21887025-1:218871:21887025",
"1:320159:32015882-1:320159:32015882",
"1:335709:33570833-1:335709:33570833",
"1:439345:43934439-1:439345:43934439",
"1:532468:53246760-1:532468:53246760",
"1:793863:79386223-1:793863:79386223",
"1:631302:63130184-1:631302:63130184",
"1:1094230:109422986-1:1094230:109422986",
"1:771785:77178473-1:771785:77178473",
"1:983713:98371293-1:983713:98371293",
"1:858054:85805306-1:858054:85805306",
"1:1237603:123760274-1:1237603:123760274",
)


)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('skim nevts:1'),
    name = cms.untracked.string('Applications')
)

# Output definition

# Additional output definition
process.SKIMStreamMySkim = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('MySkim.root'),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880)
)

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PHYS14_25_V1', '')

# Path and EndPath definitions
process.SKIMStreamMySkimOutPath = cms.EndPath(process.SKIMStreamMySkim)

# Schedule definition
process.schedule = cms.Schedule(process.SKIMStreamMySkimOutPath)
