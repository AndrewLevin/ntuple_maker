from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'flattish_nuples_80X_v6031'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'src/ntuple_maker/ntuple_maker/python/ConfFile_data_cfg.py'

config.section_("Data")


#config.Data.inputDataset = '/DoubleMuon/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016C-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016D-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016E-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016F-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016G-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016H-03Feb2017_ver3-v1/MINIAOD'

#config.Data.inputDataset = '/DoubleEG/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016C-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016D-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016E-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016F-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016G-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016H-03Feb2017_ver2-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016H-03Feb2017_ver3-v1/MINIAOD'

#config.Data.inputDataset = '/MuonEG/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2016C-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2016D-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2016E-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2016F-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2016G-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2016H-03Feb2017_ver2-v1/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2016H-03Feb2017_ver3-v1/MINIAOD'

#config.Data.inputDataset = '/SingleElectron/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#config.Data.inputDataset = '/SingleElectron/Run2016C-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/SingleElectron/Run2016D-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/SingleElectron/Run2016E-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/SingleElectron/Run2016F-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/SingleElectron/Run2016G-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/SingleElectron/Run2016H-03Feb2017_ver2-v1/MINIAOD'
config.Data.inputDataset = '/SingleElectron/Run2016H-03Feb2017_ver3-v1/MINIAOD'



config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'flattish_ntuples_80X_v6000'
config.Data.ignoreLocality = True
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'

config.Site.blacklist = ['T3_*','T2_US_Vanderbilt']





