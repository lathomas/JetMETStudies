from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'THEREQUESTNAME'
config.General.workArea = 'crabworkarea_10nov2022_new'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'

config.JobType.psetName = 'theconfig.py'
config.JobType.outputFiles = ['output.root']
#config.JobType.maxJobRuntimeMin = 1000
#config.JobType.maxMemoryMB = 2500
config.JobType.inputFiles = [
#'Summer16_07Aug2017_V11_MC.db','Summer16_07Aug2017All_V11_DATA.db','Summer16_25nsV1b_MC.db','Summer16_25nsV1b_DATA.db',
#'Fall17_17Nov2017_V32_102X_MC.db','Fall17_17Nov2017_V32_102X_DATA.db','Fall17_V3b_MC.db','Fall17_V3b_DATA.db',
#'Autumn18_V19_MC.db','Autumn18_RunABCD_V19_DATA.db','Autumn18_V7b_MC.db','Autumn18_V7b_DATA.db',
'*UL*.db',
'*.json',
'Winter22Run3_RunA_V1_DATA.db',
'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt',
'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt',
'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt',
'Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt',
'RochesterCorrections',
'UnprefireableEventList',
'Cert_Fill7733.txt'
]
#You may need to add .db file to the list above if you fetch JECs from there instead of global tag 
#Make sure all the needed files are included. For example, running on UL, you should add files such as Summer19UL17_JRV2_DATA.db Summer19UL17_RunF_V5_DATA.db
config.Data.inputDataset = 'THEDATASET'

config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 2
config.Data.allowNonValidInputDataset = True 


#config.Data.publication = True
config.Data.publication = False
config.Data.outputDatasetTag = 'CAMPAIGN_RUNERA_THESKIM'


config.Data.lumiMask = 'Cert_Collisions2016to2022_273158_357900_Golden.json'
#config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/L1Trigger/lathomas/Run3Commissioning'
#config.Site.storageSite = 'T2_CH_CERN'

config.Site.storageSite = 'T2_BE_IIHE'
#config.Site.blacklist = ['T2_US_Vanderbilt','T1_IT_CNAF']
#config.Site.whitelist = ['T1_US_FNAL']
#config.Site.ignoreGlobalBlacklist = True
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
config.User.voGroup = 'becms'
