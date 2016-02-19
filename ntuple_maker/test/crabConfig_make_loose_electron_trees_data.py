from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'fr_electron_ntuples_data_v401'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuple_maker/ntuple_maker/python/ConfFile_make_loose_electron_trees_data_cfg.py'

config.section_("Data")

config.Data.inputDataset='/DoubleEG/Run2015C_25ns-16Dec2015-v1/MINIAOD'
#config.Data.inputDataset='/DoubleEG/Run2015D-16Dec2015-v2/MINIAOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'fr_electron_ntuples_data_v400'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt'
