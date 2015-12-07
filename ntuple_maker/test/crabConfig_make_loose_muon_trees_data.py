from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'fr_muon_ntuples_v11002'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuple_maker/ntuple_maker/python/ConfFile_make_loose_muon_trees_data_cfg.py'

config.section_("Data")

#config.Data.inputDataset='/DoubleMuon/Run2015B-PromptReco-v1/MINIAOD'
#config.Data.inputDataset='/DoubleMuon/Run2015C-PromptReco-v1/MINIAOD'
#config.Data.inputDataset='/DoubleMuon/Run2015D-PromptReco-v3/MINIAOD'
config.Data.inputDataset='/DoubleMuon/Run2015D-PromptReco-v4/MINIAOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'fr_muon_ntuples_v11000'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'

config.Site.blacklist = ["T1_US_FNAL","T2_DE_DESY","T2_US_UCSD","T2_BE_UCL"]
