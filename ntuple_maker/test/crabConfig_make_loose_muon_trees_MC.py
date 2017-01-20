from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'fr_muon_ntuples_mc_v400'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'src/ntuple_maker/ntuple_maker/python/ConfFile_make_loose_muon_trees_MC_cfg.py'

config.section_("Data")

config.Data.inputDataset='/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'fr_muon_ntuples_mc_v400'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'

