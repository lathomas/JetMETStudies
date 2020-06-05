from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'dy2016_correct'
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
'Autumn18_V19_MC.db','Autumn18_RunABCD_V19_DATA.db','Autumn18_V7b_MC.db','Autumn18_V7b_DATA.db'
]
#config.JobType.inputFiles = ['Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt']

config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1


config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'RunIISummer16MiniAODv3_PUMoriond17_94X_mcRun2_asymptotic_v3_ext2_v1'
#config.Data.outputDatasetTag = 'RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1_v1'
#config.Data.outputDatasetTag = 'RunIIAutumn18MiniAOD_102X_upgrade2018_realistic_v15_ext2_v1'

config.Site.storageSite = 'T2_BE_IIHE'
config.Site.blacklist = ['T2_US_Vanderbilt','T1_IT_CNAF']

config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
config.User.voGroup = 'becms'
