from CRABClient.UserUtilities import config
config = config()


config.General.requestName = 'therequest'
config.General.workArea = 'crabworkarea'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'

config.JobType.psetName = 'JMEanalysis.py'
config.JobType.outputFiles = ['output.root']
config.JobType.maxJobRuntimeMin = 500
config.JobType.maxMemoryMB = 2500
config.JobType.inputFiles = [
'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt',
'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt',
'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt',
'Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt',
'RochesterCorrections',
'UnprefireableEventList'
]
config.Data.inputDataset = '/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/RunIISummer19UL17MiniAOD-FlatPU0to70_106X_mc2017_realistic_v6-v3/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.allowNonValidInputDataset = True
config.Data.publication = True
config.Data.outputDatasetTag = 'RunIISummer19UL17MiniAOD-FlatPU0to70_106X_mc2017_realistic_v6-v3_PuppiV15_Pt2'

config.Site.storageSite = 'T2_BE_IIHE'
config.Site.blacklist = ['T2_US_Vanderbilt','T1_IT_CNAF']

config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
config.User.voGroup = 'becms'
