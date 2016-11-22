from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'threeljj_flattish_nuples_80X_v2011'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'src/ntuple_maker/ntuple_maker/python/ConfFile_threeljj_data_cfg.py'

config.section_("Data")

#config.Data.inputDataset = '/DoubleEG/Run2016B-PromptReco-v2/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016B-PromptReco-v2/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2016B-PromptReco-v2/MINIAOD'

#config.Data.inputDataset = '/DoubleEG/Run2016C-PromptReco-v2/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016C-PromptReco-v2/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2016C-PromptReco-v2/MINIAOD'

config.Data.inputDataset = '/DoubleEG/Run2016D-PromptReco-v2/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016D-PromptReco-v2/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2016D-PromptReco-v2/MINIAOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'threeljj_flattish_ntuples_80X_v2000'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt'
