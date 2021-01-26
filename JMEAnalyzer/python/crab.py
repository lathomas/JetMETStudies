from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'qcdflatsummer20_2016apv'
config.General.workArea = 'crabworkarea'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'

config.JobType.psetName = 'JMEanalysis.py'
config.JobType.outputFiles = ['output.root']
config.JobType.maxJobRuntimeMin = 300
config.JobType.maxMemoryMB = 2500
config.JobType.inputFiles = [
'Summer16_07Aug2017_V11_MC.db','Summer16_07Aug2017All_V11_DATA.db','Summer16_25nsV1b_MC.db','Summer16_25nsV1b_DATA.db',
'Fall17_17Nov2017_V32_102X_MC.db','Fall17_17Nov2017_V32_102X_DATA.db','Fall17_V3b_MC.db','Fall17_V3b_DATA.db',
'Autumn18_V19_MC.db','Autumn18_RunABCD_V19_DATA.db','Autumn18_V7b_MC.db','Autumn18_V7b_DATA.db',
'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt',
'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt',
'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt',
'Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt',
'RochesterCorrections',
'Summer19UL18_JRV1_MC.db',
'Summer19UL18_RunA_V4_DATA.db',
'Summer19UL18_RunB_V4_DATA.db',
'Summer19UL18_RunC_V4_DATA.db',
'Summer19UL18_RunD_V4_DATA.db',
'Summer19UL18_V2_MC.db',
'UnprefireableEventList'
]
#Make sure all the needed files are included. For example, running on UL, you should add files such as Summer19UL17_JRV2_DATA.db Summer19UL17_RunF_V5_DATA.db
#config.Data.inputDataset = '/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/RunIISummer20UL16MiniAOD-FlatPU0to70_106X_mcRun2_asymptotic_v13-v1/MINIAODSIM'
config.Data.inputDataset = '/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/RunIISummer20UL16MiniAODAPV-FlatPU0to75_106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.allowNonValidInputDataset = True
config.Data.publication = True
config.Data.outputDatasetTag = 'RunIISummer20UL16MiniAODAPV-FlatPU0to75_106X_mcRun2_asymptotic_preVFP_v8-v2'


config.Site.storageSite = 'T2_BE_IIHE'
config.Site.blacklist = ['T2_US_Vanderbilt','T1_IT_CNAF']

config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
config.User.voGroup = 'becms'
