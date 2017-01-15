from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'prescale_ntuples_v800'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'src/ntuple_maker/ntuple_maker/python/ConfFile_get_prescale_information_cfg.py'

config.section_("Data")

config.Data.inputDataset='/ZeroBias/Run2016B-23Sep2016-v3/MINIAOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'prescale_ntuples_v800'
config.Data.ignoreLocality = True


config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.blacklist = ["T1_US_FNAL","T2_DE_DESY","T2_US_UCSD","T2_UK_SGrid_Bristol"]

#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
