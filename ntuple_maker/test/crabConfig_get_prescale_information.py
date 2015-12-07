from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'prescale_ntuples_v1802'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuple_maker/ntuple_maker/python/ConfFile_get_prescale_information_cfg.py'

config.section_("Data")


#config.Data.inputDataset='/ZeroBias/Run2015B-05Oct2015-v1/MINIAOD'

config.Data.inputDataset='/ZeroBias/Run2015C_25ns-05Oct2015-v1/MINIAOD'

#config.Data.inputDataset='/ZeroBias/Run2015D-PromptReco-v3/MINIAOD'
#config.Data.inputDataset='/ZeroBias/Run2015D-PromptReco-v4/MINIAOD'




config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'prescale_ntuples_v1800'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.blacklist = ["T1_US_FNAL","T2_DE_DESY","T2_US_UCSD","T2_UK_SGrid_Bristol"]

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
