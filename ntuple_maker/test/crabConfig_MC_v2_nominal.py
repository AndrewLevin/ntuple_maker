from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'flattish_nuples_mc_v2_nominal_v3006'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'src/ntuple_maker/ntuple_maker/python/ConfFile_mc_v2_nominal_cfg.py'

config.section_("Data")

config.Data.inputDataset = '/WWJJ_SS_WToLNu_EWK_aQGC-FT-FS-FM_TuneCUETP8M1_13TeV_madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'

#config.Data.inputDataset = '/WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-premix_withHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'

#config.Data.inputDataset = '/WGJJToLNu_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2-premix_withHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'

#config.Data.inputDataset = '/WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'

#config.Data.inputDataset = '/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'

#config.Data.inputDataset = '/WLLJJToLNu_M-60_EWK_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2-premix_withHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1000
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'flattish_ntuples_mc_v2_nominal_v3000'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Data.allowNonValidInputDataset = True

config.Data.ignoreLocality = True
