from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'flattish_nuples_mc_v1_nominal_v203'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'src/ntuple_maker/ntuple_maker/python/ConfFile_mc_v1_nominal_cfg.py'

config.section_("Data")

config.Data.inputDataset = '/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'

#config.Data.inputDataset ='/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'

#config.Data.inputDataset = '/WW_DoubleScattering_13TeV-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'

#config.Data.inputDataset= '/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset= '/WW_DoubleScattering_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'

#config.Data.inputDataset = '/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'

#config.Data.inputDataset = '/WpWpJJ_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1000
#config.Data.lumiMask = 'lumi_mask_JSON.txt'
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'flattish_ntuples_mc_v1_nominal_v200'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Data.allowNonValidInputDataset = True
