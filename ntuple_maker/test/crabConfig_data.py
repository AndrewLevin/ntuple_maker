from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'flattish_nuples_v808'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuple_maker/ntuple_maker/python/ConfFile_data_cfg.py'

config.section_("Data")

#config.Data.inputDataset = '/DoubleMuon/Run2015B-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2015B-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2015B-PromptReco-v1/MINIAOD'


#config.Data.inputDataset = '/DoubleMuon/Run2015C_25ns-05Oct2015-v1/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2015C_25ns-05Oct2015-v1/MINIAOD'
config.Data.inputDataset = '/DoubleEG/Run2015C_25ns-05Oct2015-v1/MINIAOD'

#config.Data.inputDataset = '/DoubleMuon/Run2015D-PromptReco-v4/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2015D-PromptReco-v4/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2015D-PromptReco-v4/MINIAOD'

#config.Data.inputDataset = '/DoubleMuon/Run2015D-PromptReco-v3/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2015D-PromptReco-v3/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2015D-PromptReco-v3/MINIAOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'flattish_ntuples_v800'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
