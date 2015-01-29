from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'flattish_nuples_v9'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuple_maker/ntuple_maker/python/ConfFile_cfg.py'

config.section_("Data")
#config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU40bx25_POSTLS170_V7-v2/MINIAODSIM'
#config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU40bx25_POSTLS170_V5-v2/MINIAODSIM'
#config.Data.inputDataset= '/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset= '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDataset= '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1000
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'flattish_ntuples_v9'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
