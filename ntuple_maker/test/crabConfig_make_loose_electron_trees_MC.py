from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'fr_electron_ntuples_mc_v200'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuple_maker/ntuple_maker/python/ConfFile_make_loose_electron_trees_MC_cfg.py'

config.section_("Data")

#config.Data.inputDataset='/DoubleEG/Run2015B-PromptReco-v1/MINIAOD'
#config.Data.inputDataset='/DoubleEG/Run2015C-PromptReco-v1/MINIAOD'
#config.Data.inputDataset='/DoubleEG/Run2015D-PromptReco-v3/MINIAOD'
#config.Data.inputDataset='/DoubleEG/Run2015D-PromptReco-v4/MINIAOD'

config.Data.inputDataset='/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'fr_electron_ntuples_mc_v200'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'

#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
