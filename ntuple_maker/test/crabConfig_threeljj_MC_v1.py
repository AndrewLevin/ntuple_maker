from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'threeljj_flattish_nuples_mc_v1_v300'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuple_maker/ntuple_maker/python/ConfFile_threeljj_mc_v1_cfg.py'

config.section_("Data")

#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/MINIAODSIM'
#config.Data.inputDataset = '/WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#config.Data.inputDataset='/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset='/WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1000
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'threeljj_flattish_ntuples_mc_v1_v300'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Data.allowNonValidInputDataset = True
