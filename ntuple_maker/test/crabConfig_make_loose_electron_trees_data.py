from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'fr_electron_ntuples_data_v10006'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'src/ntuple_maker/ntuple_maker/python/ConfFile_make_loose_electron_trees_data_cfg.py'

config.section_("Data")

#config.Data.inputDataset = '/DoubleEG/Run2016B-23Sep2016-v3/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016C-23Sep2016-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016D-23Sep2016-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016E-23Sep2016-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016F-23Sep2016-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016G-23Sep2016-v1/MINIAOD'
config.Data.inputDataset = '/DoubleEG/Run2016H-PromptReco-v2/MINIAOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'fr_electron_ntuples_data_v10000'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'

#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'

config.Site.blacklist = ["T2_US_Vanderbilt","T2_US_Florida","T2_US_Wisconsin"]
