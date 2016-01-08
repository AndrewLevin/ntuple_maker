from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'll_flattish_nuples_v503'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuple_maker/ntuple_maker/python/ConfFile_ll_data_cfg.py'

config.section_("Data")

#config.Data.inputDataset = '/DoubleEG/Run2015D-PromptReco-v3/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2015D-PromptReco-v4/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2015B-PromptReco-v1/MINIAOD'
config.Data.inputDataset = '/DoubleEG/Run2015C-PromptReco-v1/MINIAOD'

#config.Data.inputDataset = '/DoubleMuon/Run2015D-PromptReco-v3/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2015D-PromptReco-v4/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2015B-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2015C-PromptReco-v1/MINIAOD'


config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'll_flattish_ntuples_v500'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
