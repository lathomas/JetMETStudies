
#TheSkim = "MCJECs"
TheSkim = "ZJetsResiduals"
TheSkim = "MCJECs" 
TheSkim = "HFJet"
TheSkim = "L1Unprefirable"
TheSkim = "ZJetsResiduals"
TheSkim = "MCJECs"
TheSkim = "L1Study"
#TheSkim = "L1Study_ZToMuMu"
TheSkim = "L1Study_ZToEE"
#TheSkim = "L1Study_SingleMuforJME"
ReclusterCHSJets = False
ReclusterGenJets = False
#TheSkim = ""

#runEra="DataUL2017F"
#runEra="MCUL2017"
#runEra="MCRun3"
runEra="DataRun3"
#runEra="DataUL2018D"
ISMC=True
if "MC" in runEra:
    ISMC=True
if not ISMC:
    ReclusterGenJets = False

UseSQLiteFiles=True


if TheSkim == "MCJECs":
   ReclusterCHSJets = True
   ReclusterGenJets = True


import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )  
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(23374) )  

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)


process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(

#2018
#'/store/mc/RunIISummer19UL18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/100000/00556E5A-E6CD-6C4B-83BA-0B5C00C0A5E4.root'
#2017
#'/store/mc/RunIISummer19UL17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/100000/002C691B-A0CE-A24F-8805-03B4C52C9004.root'
#2016
#'/store/mc/RunIISummer20UL16MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v13-v2/260000/003F1A76-9CDA-7644-A34E-923C4B1C0E5E.root'

#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_HighStat-v2/2580000/1026a87e-952c-4996-8ee7-57f5f0d63bb9.root',
#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_HighStat-v2/2580000/1fa066fd-adf9-440d-8f43-db3f68c28d2e.root',
#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_HighStat-v2/2580000/4fe44364-63b4-49a2-9577-38a6278cdf70.root',
#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_HighStat-v2/2580000/9f84ce7c-602d-465f-a4eb-6ce431a00de0.root'
#'/store/mc/Run3Winter22MiniAOD/SingleNeutrino_E-10-gun/MINIAODSIM/L1TPU0to99FEVT_SNB_122X_mcRun3_2021_realistic_v9-v2/2520000/04c888dc-4034-4fed-8616-5aa0663f84f6.root'
#'file:/user/lathomas/L1Studies/SampleGeneration/SingleNeutrinoPU1/CMSSW_12_2_1/src/SingleNuPU1_MINIAOD_500.root'
#'file:/user/lathomas/L1Studies/SampleGeneration/SingleNeutrinoPU1/CMSSW_12_2_1/src/jobsubmission/output_MINIAOD.root'
#'file:/user/lathomas/Run3DQM/CMSSW_12_3_4_patch2/src/JetMETStudies/JMEAnalyzer/python/files2022/0634e291-f4b1-49fc-a542-dde43480831a.root'
#'/store/data/Run2022C/EGamma/MINIAOD/PromptReco-v1/000/356/071/00000/861766b9-f568-4c74-96be-9a3d037b087d.root'
'/store/mc/Run3Winter22MiniAOD/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/MINIAODSIM/122X_mcRun3_2021_realistic_v9_ext1-v1/2830000/006d687f-44da-49f0-b0c0-f46b8600731f.root'
#'file:MINIAODRun3.root'
#'/store/data/Run2022C/SingleMuon/MINIAOD/PromptReco-v1/000/356/071/00000/40fbb456-859b-4f95-880a-15491af226fe.root'
#'/store/data/Run2022C/SingleMuon/MINIAOD/PromptReco-v1/000/356/071/00000/41ea401b-9ada-4481-b948-9c831fcafeb3.root',
#'/store/data/Run2022C/SingleMuon/MINIAOD/PromptReco-v1/000/356/071/00000/42073b59-45a3-4af0-ba64-4ad328e2f91e.root',
#'/store/data/Run2022C/SingleMuon/MINIAOD/PromptReco-v1/000/356/071/00000/5206522d-1e9b-4a67-be27-29f816054ca4.root',
#'/store/data/Run2022C/SingleMuon/MINIAOD/PromptReco-v1/000/356/071/00000/6706e10d-3638-4527-baab-b3598f667512.root',
#'/store/data/Run2022C/SingleMuon/MINIAOD/PromptReco-v1/000/356/071/00000/6f4eb734-caff-4e5d-b522-3145e7d6ff38.root',
#'/store/data/Run2022C/SingleMuon/MINIAOD/PromptReco-v1/000/356/071/00000/72f1d890-6801-4f60-b965-d20ac4a34f65.root'
#'/store/data/Run2022B/SingleMuon/MINIAOD/PromptReco-v1/000/355/102/00000/3a5f2a35-f621-4e05-9f9a-e259db2c4c8c.root'

##'file:/user/lathomas/Run3DQM/CMSSW_12_3_4_patch2/src/JetMETStudies/JMEAnalyzer/python/files2022/JetHT/6c336f02-ea66-4df5-bdf9-4088f7999493.root',
##'file:/user/lathomas/Run3DQM/CMSSW_12_3_4_patch2/src/JetMETStudies/JMEAnalyzer/python/files2022/JetHT/99ab6f95-1100-4615-9f5c-45d0117fd084.root',
##'file:/user/lathomas/Run3DQM/CMSSW_12_3_4_patch2/src/JetMETStudies/JMEAnalyzer/python/files2022/JetHT/b281a525-3f21-49c0-bb51-e84a540121ae.root',
#'file:/user/lathomas/Run3DQM/CMSSW_12_3_4_patch2/src/JetMETStudies/JMEAnalyzer/python/files2022/JetHT/e5d5dd34-1d4a-4478-857e-d62e3e60ae63.root',
#'file:/user/lathomas/Run3DQM/CMSSW_12_3_4_patch2/src/JetMETStudies/JMEAnalyzer/python/files2022/JetHT/3a7aad45-2acd-4b0e-900b-b53fd84eab99.root'

#'/store/relval/CMSSW_12_3_1/RelValZpToEE_m6000_14TeV/MINIAODSIM/PU_123X_mcRun3_2021_realistic_v13-v1/2580000/2ef5a51d-34cf-4341-9236-031697fc1711.root'
#'/store/data/Run2022A/ZeroBias10/MINIAOD/PromptReco-v1/000/353/689/00000/718b27d7-376a-4cab-a569-39671c815605.root'
#'file:files/05ddba8a-e911-476e-a25e-36b63740d3ab.root',
#'file:files/118b95e1-347e-475b-8fdd-c393bc14f18e.root',
#'file:files/171d91a8-8eb5-4738-9a33-56df8968385f.root',
#'file:files/37f35043-a4cc-46ac-845b-9e359947f6b0.root',
#'file:files/3ef6099d-116b-468d-b41b-8a6cb96c7ff2.root',
#'file:files/489dd159-007f-4591-becb-ec32fbb0d383.root',
#'file:files/5bb4f831-bdaa-417d-9126-b59644ee6517.root',
#'file:files/6b76e233-a175-4932-9456-897537a29e90.root',
#'file:files/6f816e2e-6cb8-4c43-8b55-0737c6a28ca3.root',
#'file:files/74df2149-16bc-4975-87c2-d07a2c5cf6a3.root',
#'file:files/9d38101d-d57a-40ba-9452-dcb76759bd50.root',
#'file:files/c12aa1ed-8a35-4aeb-96a5-572c07e5bb9b.root',
#'file:files/c55a1631-7234-4874-a997-c096f3cacf98.root',
#'file:files/f80d97eb-3615-4de8-bcba-0e2144546a7f.root'
#'file:miniaod.root'
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/10e5f0ef-3852-498e-bf3a-a4a073319b56.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/20da8959-7aec-44d0-b6b7-76f5d2fb4e17.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/3d9ad4fa-3781-475c-9257-4612bcd9d894.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/49b7b016-0286-477c-8b99-245774356fb8.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/54f285f6-82e6-4886-a514-3322a1d9d299.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/6a32adbf-fdc9-450c-af6c-ee9c3105878f.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/7548aa25-fda9-49e3-a68b-f637b4efca33.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/7f89bc06-949c-4572-9b8d-a500c6f9006a.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/8457731f-aba8-4350-b04d-e007788d0f68.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/92c22aa9-4b26-4891-9f5e-abae200c3a6c.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/9fa281ca-3c8b-4315-96fb-27d2d89efd38.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/b5071877-90b9-46c4-914f-6e7d92a20624.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/cdb7b0f0-f9d7-4266-a080-036814104f67.root',
#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/f046e8f8-1d0b-4394-89c7-528475756bef.root',

#'/store/data/Commissioning2021/ZeroBias1/MINIAOD/PromptReco-v1/000/346/512/00000/f1e3f57d-9fba-40e7-aa9d-f4444c10ade0.root'
#'/store/data/Commissioning2022/ZeroBias2/MINIAOD/PromptReco-v1/000/352/135/00000/f28f464f-50f1-4789-9741-b81fa125b5e0.root'
#'/store/data/Commissioning2022/ZeroBias2/MINIAOD/PromptReco-v1/000/352/173/00000/6d5c8157-a4c5-456c-8742-054f9a348363.root'
#'/store/data/Run2022A/ZeroBias2/MINIAOD/PromptReco-v1/000/352/417/00000/5ed54d9c-04ed-4eb6-b9a0-6b76a8703ae0.root'
#'/store/data/Run2022A/ZeroBias5/MINIAOD/PromptReco-v1/000/352/425/00000/2c45e4e9-3aed-4a5a-a68c-700461ffa202.root'
#'/store/data/Run2022A/ZeroBias2/MINIAOD/PromptReco-v1/000/352/417/00000/9b44dc16-758e-46ff-90dc-f547cd70d505.root'
#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_TkmkFitHighStat-v2/2580000/130135b7-81a8-4bb8-8348-20960f027375.root',
#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_TkmkFitHighStat-v2/2580000/99eec549-5ce3-4bdc-a93a-b7c4bfd6745e.root',
#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_TkmkFitHighStat-v2/2580000/d2c99988-cf49-4ba2-8b72-e9d6d4761ad2.root',
#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_TkmkFitHighStat-v2/2580000/d6d8c549-727f-4073-88f1-1042167bff67.root'



#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_HighStat-v2/2580000/1026a87e-952c-4996-8ee7-57f5f0d63bb9.root'
#'/store/mc/RunIISummer20UL17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/250000/EBB47EA5-4B28-7B42-BFC8-C0F75800C7AC.root'
#'/store/mc/RunIISummer20UL18MiniAODv2/QCD_bEnriched_HT100to200_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/120000/09FB91A8-7902-AF42-AFEC-AFB6A294955E.root'
#'/store/mc/RunIISummer19UL18MiniAODv2/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/FlatPU0to70_106X_upgrade2018_realistic_v16_L1v1-v1/280000/0024EE42-B032-5E49-A78D-152B70540F9C.root'


#'/store/mc/RunIISummer20UL18MiniAODv2/QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/110000/C46792AE-B560-3541-B67A-DD73050C49A8.root'

#'/store/data/Run2018D/DoubleMuon/MINIAOD/UL2018_MiniAODv2-v1/270002/152AABA5-91A3-9F4B-82BE-654A96DEAE98.root'
#Data
#'/store/data/Run2017F/DoubleMuon/MINIAOD/09Aug2019_UL2017-v1/270000/527C5A3A-7C09-4F42-B9DA-A84871504EBF.root'
#'/store/data/Run2018A/DoubleMuon/MINIAOD/12Nov2019_UL2018-v2/270000/F0CC36DA-D61E-0F4B-BAF4-1FF8913DBB78.root'


#'file:/user/lathomas/SUSYSSDilepton/new/New/CMSSW_10_6_23/src/pickevents_2018a.root',
#'file:/user/lathomas/SUSYSSDilepton/new/New/CMSSW_10_6_23/src/pickevents_2018b.root',
#'file:/user/lathomas/SUSYSSDilepton/new/New/CMSSW_10_6_23/src/pickevents_2018c.root',
#'file:/user/lathomas/SUSYSSDilepton/new/New/CMSSW_10_6_23/src/pickevents_2018d.root'


#'/store/data/Run2016E/JetHT/MINIAOD/21Feb2020_UL2016_HIPM-v1/10000/30434744-4792-9C4F-A1D3-87ACE1D2E6E0.root'
#Here's a MINIAOD file on lxplus in case you want to use a local sample
#'file:/afs/cern.ch/work/l/lathomas/public/qcdht1000to1500_1.root'
        )
                            )

#process.TFileService = cms.Service("TFileService", fileName = cms.string("outputQCDHT1000to1500_puppiv16_200kevts.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("outputQCDHT1000to1500_puppiv16_23374evts.root") )
process.TFileService = cms.Service("TFileService", fileName = cms.string("output.root") )

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag="102X_dataRun2_v8"


process.GlobalTag.globaltag="123X_mcRun3_2021_realistic_v13"

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")


process.load('Configuration.StandardSequences.Reconstruction_cff')

if ReclusterCHSJets: 
    CHSJetCollectionName ="selectedPatJetsCHS"
else: 
    CHSJetCollectionName ="updatedPatJetsUpdatedJEC"

CHSJetCollectionName ="slimmedJets"

if ReclusterGenJets:
    GenJetCollectionName="ak4GenJetsNoNuNEW"
else:
    GenJetCollectionName="slimmedGenJets"


RochesterCorrectionFile="RochesterCorrections/"




#EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
#EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
#EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
#PhotonTightWP='mvaPhoID-RunIIFall17-v2-wp80'
#RochesterCorrectionFile+="RoccoR2018UL.txt" #Muon POG hasn't released Rochester corrections for UL18 yet


#Starting with pre UL data
if "Data2016" in runEra:
    process.GlobalTag.globaltag="102X_dataRun2_v13" #2016
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2016.txt"

if "Data2017" in runEra:
    process.GlobalTag.globaltag="102X_dataRun2_v13" #2017      
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2017.txt"

if "Data2018" in runEra:
    if "2018D" in runEra: 
        process.GlobalTag.globaltag="102X_dataRun2_Prompt_v16" #2018D
    else:
        process.GlobalTag.globaltag="102X_dataRun2_v13" #2018ABC
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2018.txt"

#Now pre UL MC
if "MC2016" in runEra:
    process.GlobalTag.globaltag="102X_mcRun2_asymptotic_v8" #2016
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2016.txt"

if "MC2017" in runEra:
    process.GlobalTag.globaltag="102X_mc2017_realistic_v8" #2017
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2017.txt"

if "MC2018" in runEra:
    process.GlobalTag.globaltag="102X_upgrade2018_realistic_v21" #2018     
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2018.txt"


#Now UL data
#2016 is divided into two parts
if "DataUL2016B" in runEra or "DataUL2016C" in runEra or "DataUL2016D" in runEra or "DataUL2016E" in runEra or "DataUL2016F" in runEra:
    process.GlobalTag.globaltag="106X_dataRun2_v32" #UL2016
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v2-wp80'
    RochesterCorrectionFile+="RoccoR2016aUL.txt" 

if "DataUL2016Flate" in runEra or "DataUL2016G" in runEra or "DataUL2016H" in runEra:
    process.GlobalTag.globaltag="106X_dataRun2_v32" #UL2016
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v2-wp80'
    RochesterCorrectionFile+="RoccoR2016bUL.txt" 

if "DataUL2017" in runEra:
    process.GlobalTag.globaltag="106X_dataRun2_v32" #UL2017      
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v2-wp80'
    RochesterCorrectionFile+="RoccoR2017UL.txt"

if "DataUL2018" in runEra:
    process.GlobalTag.globaltag="106X_dataRun2_v32" #UL2018      
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v2-wp80'
    RochesterCorrectionFile+="RoccoR2018UL.txt"

if "DataRun3" in runEra:
    process.GlobalTag.globaltag="123X_dataRun3_Prompt_v12" #Run 3
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v2-wp80'
    RochesterCorrectionFile+="RoccoR2018UL.txt"


#Now UL MC

if "MCUL2016APV" in runEra:
    process.GlobalTag.globaltag="106X_mcRun2_asymptotic_preVFP_v11" #2016
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v2-wp80'
    RochesterCorrectionFile+="RoccoR2016aUL.txt"

if "MCUL2016nonAPV" in runEra:
    process.GlobalTag.globaltag="106X_mcRun2_asymptotic_v17" #2016
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v2-wp80'
    RochesterCorrectionFile+="RoccoR2016bUL.txt"

if "MCUL2017" in runEra:
    process.GlobalTag.globaltag="106X_mc2017_realistic_v8" #UL2017
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v2-wp80'
    RochesterCorrectionFile+="RoccoR2017UL.txt"

if "MCUL2018" in runEra:
    process.GlobalTag.globaltag="106X_upgrade2018_realistic_v15" #UL2018
    EleVetoWP='cutBasedElectronID-Fall17-94X-V2-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    EleLooseWP='mvaEleID-Fall17-iso-V2-wpHZZ'
    PhotonTightWP='mvaPhoID-RunIIFall17-v2-wp80'
    RochesterCorrectionFile+="RoccoR2018UL.txt" #Muon POG hasn't released Rochester corrections for UL18 yet   

print("Roch corr file: ")
print(RochesterCorrectionFile)

process.jmeanalyzer = cms.EDAnalyzer('JMEAnalyzer',
                                     METFiltersPAT = cms.InputTag("TriggerResults::PAT"),
                                     METFiltersRECO = cms.InputTag("TriggerResults::RECO"),
                                     ECALBadCalibFilterUpdate=cms.InputTag("ecalBadCalibReducedMINIAOD2019Filter"),
                                     ECALLaserCorrFilterUpdate=cms.InputTag("ecalLaserCorrFilter"),
                                     ECALDeadCellBoundaryEnergyFilterUpdate=cms.InputTag("ecalDeadCellBoundaryEnergyFilterUpdate"),
                                     BadChargedCandidateFilterUpdate=cms.InputTag("BadChargedCandidateFilterUpdate"),
                                     Vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                     Jets=cms.InputTag(CHSJetCollectionName),
                                     JetsAK8=cms.InputTag("selectedPatJetsAK8CHS"),
#                                     JetsPuppi=cms.InputTag("updatedPatJetsUpdatedJECPuppi"),
                                     JetsPuppi=cms.InputTag("slimmedJetsPuppi"),
                                     JetsPuppiAK8=cms.InputTag("slimmedJetsAK8"),
                                     JetsCalo=cms.InputTag("slimmedCaloJets"),
                                     JetsPFnoCHS=cms.InputTag("selectedPatJetsPlain"),
                                     GenJets=cms.InputTag(GenJetCollectionName),
                                     GenAK8Jets=cms.InputTag("slimmedGenJetsAK8"),
                                     pileupJetIdDiscriminantUpdate = cms.InputTag('pileupJetIdUpdate:fullDiscriminant'),
                                     pileupJetIdDiscriminantUpdate2017 = cms.InputTag('pileupJetIdUpdate2017:fullDiscriminant'),
                                     pileupJetIdDiscriminantUpdate2018 = cms.InputTag('pileupJetIdUpdate2018:fullDiscriminant'),
                                     pileupJetIdVariablesUpdate = cms.InputTag('pileupJetIdUpdate'),
                                     QuarkGluonLikelihood = cms.InputTag('QGTagger:qgLikelihood'),
                                     PFCandidates=cms.InputTag("packedPFCandidates"),
                                     PuppiWeights=cms.InputTag("puppi"),
                                     PULabel = cms.InputTag("slimmedAddPileupInfo"),
                                     Triggers = cms.InputTag("TriggerResults::HLT"),
                                     l1GtSrc = cms.InputTag("gtStage2Digis"),
                                     GenParticles=cms.InputTag("prunedGenParticles"),
                                     GenInfo=cms.InputTag("generator"),
                                     LHELabel = cms.InputTag("externalLHEProducer"),
                                     LHELabelALT = cms.InputTag("source"),
                                     GenJetMatchCHS= cms.InputTag("patJetGenJetMatchUpdate"),
                                     GenJetWithNuMatchCHS= cms.InputTag("patJetGenWithNuJetMatchUpdate"),
                                     GenJetMatchPuppi= cms.InputTag("patJetGenJetMatchUpdatePuppi"),
                                     GenJetWithNuMatchPuppi= cms.InputTag("patJetGenWithNuJetMatchUpdatePuppi"),
                                     GenJetMatchCalo= cms.InputTag("patJetGenJetMatchUpdateCalo"),
                                     PFMet=cms.InputTag("slimmedMETs"),
                                     PuppiMet=cms.InputTag("slimmedMETsPuppi"),
                                     Electrons=cms.InputTag("slimmedElectrons"),
                                     Muons=cms.InputTag("slimmedMuons"),
                                     Photons=cms.InputTag("slimmedPhotons"),
                                     JetPtCut=cms.double(2),
                                     AK8JetPtCut=cms.double(10),
                                     ElectronPtCut=cms.double(2),
                                     ElectronVetoWorkingPoint=cms.string(EleVetoWP),
                                     ElectronLooseWorkingPoint=cms.string(EleLooseWP),
                                     ElectronTightWorkingPoint=cms.string(EleTightWP),
                                     MuonPtCut=cms.double(1),
                                     RochCorrFile=cms.string(RochesterCorrectionFile),
                                     PhotonPtCut=cms.double(5),
                                     PhotonTightWorkingPoint=cms.string(PhotonTightWP),
                                     PFCandPtCut=cms.double(25000),
                                     SaveTree=cms.bool(True),
                                     IsMC=cms.bool(ISMC),
                                     SavePUIDVariables=cms.bool(False),
                                     SaveAK8Jets=cms.bool(False),
                                     SaveCaloJets=cms.bool(True),
                                     SavenoCHSJets=cms.bool(True),
                                     DropUnmatchedJets=cms.bool(True),
                                     DropBadJets=cms.bool(False),
                                     SavePFinJets=cms.bool(False),
                                     ApplyPhotonID=cms.bool(False),
                                     Skim=cms.string(TheSkim),
#                                     Skim=cms.string("FourLeptons"),
#                                     Skim=cms.string("L1Unprefirable"),
#                                     Skim=cms.string("HighHT"),
                                     Debug=cms.bool(False)
                              )


if TheSkim == "MCJECs":
    process.jmeanalyzer.JetPtCut=cms.double(-1)
    process.jmeanalyzer.AK8JetPtCut=cms.double(10)
    process.jmeanalyzer.SaveAK8Jets=cms.bool(True)
    process.jmeanalyzer.SaveCaloJets=cms.bool(True)
    process.jmeanalyzer.SavenoCHSJets=cms.bool(True)
    process.jmeanalyzer.DropUnmatchedJets=cms.bool(True)



if TheSkim == "ZJetsResiduals" or TheSkim == "GammaJetsResiduals":
    process.jmeanalyzer.JetPtCut=cms.double(-1)
    process.jmeanalyzer.AK8JetPtCut=cms.double(1000)
    process.jmeanalyzer.PhotonPtCut=cms.double(20)
    process.jmeanalyzer.ElectronPtCut=cms.double(10)
    process.jmeanalyzer.MuonPtCut=cms.double(10)
    process.jmeanalyzer.ApplyPhotonID=cms.bool(True)
    process.jmeanalyzer.SaveAK8Jets=cms.bool(False)
    process.jmeanalyzer.SaveCaloJets=cms.bool(False)
    process.jmeanalyzer.SavenoCHSJets=cms.bool(False)
    process.jmeanalyzer.DropUnmatchedJets=cms.bool(False)
    process.jmeanalyzer.DropBadJets=cms.bool(True)

if TheSkim == "FourLeptons":
    process.jmeanalyzer.ElectronPtCut=cms.double(10)
    process.jmeanalyzer.MuonPtCut=cms.double(10)
    process.jmeanalyzer.JetPtCut=cms.double(25)
    process.jmeanalyzer.AK8JetPtCut=cms.double(1000)
    process.jmeanalyzer.PhotonPtCut=cms.double(2000)
    process.jmeanalyzer.SaveAK8Jets=cms.bool(False)
    process.jmeanalyzer.SaveCaloJets=cms.bool(False)
    process.jmeanalyzer.SavenoCHSJets=cms.bool(False)
    process.jmeanalyzer.DropUnmatchedJets=cms.bool(False)
    process.jmeanalyzer.DropBadJets=cms.bool(True)

if TheSkim == "PhotonHFJet" or TheSkim == "HFJet" or TheSkim == "ZHFJet":
    process.jmeanalyzer.JetPtCut=cms.double(30)
    process.jmeanalyzer.AK8JetPtCut=cms.double(1000)
    process.jmeanalyzer.PhotonPtCut=cms.double(20)
    process.jmeanalyzer.ElectronPtCut=cms.double(10)
    process.jmeanalyzer.MuonPtCut=cms.double(10)
    process.jmeanalyzer.ApplyPhotonID=cms.bool(True)
    process.jmeanalyzer.SaveAK8Jets=cms.bool(False)
    process.jmeanalyzer.SaveCaloJets=cms.bool(False)
    process.jmeanalyzer.SavenoCHSJets=cms.bool(False)
    process.jmeanalyzer.DropUnmatchedJets=cms.bool(False)
    process.jmeanalyzer.DropBadJets=cms.bool(False)

if TheSkim == "L1Unprefirable":
    process.jmeanalyzer.JetPtCut=cms.double(20)
    process.jmeanalyzer.AK8JetPtCut=cms.double(1000)
    process.jmeanalyzer.PhotonPtCut=cms.double(15)
    process.jmeanalyzer.ElectronPtCut=cms.double(10)
    process.jmeanalyzer.MuonPtCut=cms.double(10)
    process.jmeanalyzer.ApplyPhotonID=cms.bool(False)
    process.jmeanalyzer.SaveAK8Jets=cms.bool(False)
    process.jmeanalyzer.SaveCaloJets=cms.bool(False)
    process.jmeanalyzer.SavenoCHSJets=cms.bool(False)
    process.jmeanalyzer.DropUnmatchedJets=cms.bool(False)
    process.jmeanalyzer.DropBadJets=cms.bool(False)

if TheSkim == "L1Study_ZToMuMu" or TheSkim == "L1Study_ZToEE":
    process.jmeanalyzer.JetPtCut=cms.double(20000)
    process.jmeanalyzer.AK8JetPtCut=cms.double(20000)
    process.jmeanalyzer.PhotonPtCut=cms.double(20000)
    process.jmeanalyzer.ElectronPtCut=cms.double(10)
    process.jmeanalyzer.MuonPtCut=cms.double(5)
    process.jmeanalyzer.ApplyPhotonID=cms.bool(False)
    process.jmeanalyzer.SaveAK8Jets=cms.bool(False)
    process.jmeanalyzer.SaveCaloJets=cms.bool(False)
    process.jmeanalyzer.SavenoCHSJets=cms.bool(False)
    process.jmeanalyzer.DropUnmatchedJets=cms.bool(False)
    process.jmeanalyzer.DropBadJets=cms.bool(False)
    
if TheSkim == "L1Study_SingleMuforJME" or TheSkim == "L1Study_SinglePhotonforJME":
    process.jmeanalyzer.JetPtCut=cms.double(20)
    process.jmeanalyzer.AK8JetPtCut=cms.double(20000)
    process.jmeanalyzer.PhotonPtCut=cms.double(20)
    process.jmeanalyzer.ElectronPtCut=cms.double(10)
    process.jmeanalyzer.MuonPtCut=cms.double(10)
    process.jmeanalyzer.ApplyPhotonID=cms.bool(False)
    process.jmeanalyzer.SaveAK8Jets=cms.bool(False)
    process.jmeanalyzer.SaveCaloJets=cms.bool(False)
    process.jmeanalyzer.SavenoCHSJets=cms.bool(False)
    process.jmeanalyzer.DropUnmatchedJets=cms.bool(False)
    process.jmeanalyzer.DropBadJets=cms.bool(False)
    
#Rerunning the ecalbadcalibration filter
from RecoMET.METFilters.ecalBadCalibFilter_cfi import ecalBadCalibFilter

baddetEcallistnew2019 = cms.vuint32(
    [872439604,872422825,872420274,872423218,872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,872437052,872420649,872421950,872437185,
     872422564,872421566,872421695,872421955,872421567,872437184,872421951,872421694,
     872437056,872437057,872437313,872438182,872438951,872439990,872439864,872439609,
     872437181,872437182,872437053,872436794,872436667,872436536,872421541,872421413,
     872421414,872421031,872423083,872421439,872423224,872421438,872420397,872421566,
     872422589,872423096,872422717,872423214,872421415,872422311,872421926,872439469,
     872438567,872436659,872439731,872438311,872438078,872438438,872439601,872437951,
     872437950,872439729,872436792,872438183,872439468,872436663,872439728,872439727,
     872437694,872437823,872438845,872438973,872439354,872438566,872439733,872436530,
     872436655,872439600,872439730]
    )

process.ecalBadCalibReducedMINIAOD2019Filter = ecalBadCalibFilter.clone(
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallistnew2019,
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )

#Rerunning the laser correction filter
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
process.ecalLaserCorrFilter = cms.EDFilter(
    "EcalLaserCorrFilter",
    EBRecHitSource = cms.InputTag("reducedEgamma:reducedEBRecHits"),
    EERecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    EBLaserMIN     = cms.double(0.3),
    EELaserMIN     = cms.double(0.3),
    EBLaserMAX     = cms.double(5.0), #this was updated wrt default
    EELaserMAX     = cms.double(100.0), #this was updated wrt default
    EBEnegyMIN     = cms.double(10.0),
    EEEnegyMIN     = cms.double(10.0),
    taggingMode    = cms.bool(True), #updated wrt default
    Debug          = cms.bool(False)
    )

#Rerunning EcalDeadCellBoundaryEnergyFilter
from RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi import EcalDeadCellBoundaryEnergyFilter
process.ecalDeadCellBoundaryEnergyFilterUpdate=EcalDeadCellBoundaryEnergyFilter.clone(
    recHitsEB = cms.InputTag("reducedEgamma:reducedEBRecHits"),
    recHitsEE = cms.InputTag("reducedEgamma:reducedEERecHits"),
    cutBoundEnergyDeadCellsEE=cms.untracked.double(10),
    taggingMode    = cms.bool(True)
    )

#Rerunning BadChargedCandidateFilter
from RecoMET.METFilters.BadChargedCandidateFilter_cfi import BadChargedCandidateFilter 
process.BadChargedCandidateFilterUpdate=BadChargedCandidateFilter.clone(
    muons = cms.InputTag("slimmedMuons"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    vtx = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taggingMode    = cms.bool(True)
)



import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())


JSONfile =''

if "DataUL2017" in runEra:
    JSONfile = 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
if "Data2018" in runEra or "DataUL2018" in runEra:
    JSONfile = 'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
if "Data2017" in runEra:
    JSONfile = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'
if "Data2016" in runEra or "DataUL2016" in runEra:
    JSONfile = 'Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt'
if "DataRun3" in runEra:
    JSONfile = 'Cert_Collisions2022_355100_357900_Golden.json' 

myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
if not ISMC:
    process.source.lumisToProcess.extend(myLumis)
print( "json" )
print( JSONfile )


#Updating JECs
JECsVersion=""
if "MC2018" in runEra:
    JECsVersion='Autumn18_V19_MC'
if "MC2017" in runEra:
    JECsVersion='Fall17_17Nov2017_V32_102X_MC'
if "MC2016" in runEra:
    JECsVersion='Summer16_07Aug2017_V11_MC'

if "MCUL2016APV" in runEra:
    JECsVersion='Summer19UL16APV_V7_MC'
if "MCUL2016nonAPV" in runEra:
    JECsVersion='Summer19UL16_V7_MC'
if "MCUL2017" in runEra:
    JECsVersion='Summer19UL17_V6_MC'
if "MCUL2018" in runEra:
    JECsVersion='Summer19UL18_V5_MC'

if "MCRun3" in runEra:
    JECsVersion='Summer19UL18_V5_MC'
if "DataRun3" in runEra:
    JECsVersion='Winter22Run3_RunA_V1_DATA'

if "Data2018" in runEra:
    JECsVersion='Autumn18_RunABCD_V19_DATA'
if "Data2017" in runEra:
    JECsVersion='Fall17_17Nov2017_V32_102X_DATA'
if "Data2016" in runEra or "DataUL2016" in runEra:
    JECsVersion='Summer16_07Aug2017All_V11_DATA'

if "DataUL2017B" in runEra:
    JECsVersion='Summer19UL17_RunB_V6_DATA'
if "DataUL2017C" in runEra:
    JECsVersion='Summer19UL17_RunC_V6_DATA'
if "DataUL2017D" in runEra:
    JECsVersion='Summer19UL17_RunD_V6_DATA'
if "DataUL2017E" in runEra:
    JECsVersion='Summer19UL17_RunE_V6_DATA'
if "DataUL2017F" in runEra:
    JECsVersion='Summer19UL17_RunF_V6_DATA'

if "DataUL2018A" in runEra:
    JECsVersion='Summer19UL18_RunA_V5_DATA'
if "DataUL2018B" in runEra:
    JECsVersion='Summer19UL18_RunB_V5_DATA'
if "DataUL2018C" in runEra:
    JECsVersion='Summer19UL18_RunC_V5_DATA'
if "DataUL2018D" in runEra:
    JECsVersion='Summer19UL18_RunD_V5_DATA'

if "DataUL2016B" in runEra:
    JECsVersion='Summer19UL16APV_RunBCD_V7_DATA'
if "DataUL2016C" in runEra:
    JECsVersion='Summer19UL16APV_RunBCD_V7_DATA'
if "DataUL2016D" in runEra:
    JECsVersion='Summer19UL16APV_RunBCD_V7_DATA'
if "DataUL2016E" in runEra:
    JECsVersion='Summer19UL16APV_RunEF_V7_DATA'
if "DataUL2016F" in runEra:
    JECsVersion='Summer19UL16APV_RunEF_V7_DATA'
if "DataUL2016Flate" in runEra:
    JECsVersion='Summer19UL16_RunFGH_V7_DATA'
if "DataUL2016G" in runEra:
    JECsVersion='Summer19UL16_RunFGH_V7_DATA'
if "DataUL2016H" in runEra:
    JECsVersion='Summer19UL16_RunFGH_V7_DATA'
 


SQLiteFile='sqlite:'+JECsVersion+'.db'

TagForAK4CHSJet='JetCorrectorParametersCollection_'+JECsVersion+'_AK4PFchs'
TagForAK4PuppiJet='JetCorrectorParametersCollection_'+JECsVersion+'_AK4PFPuppi'

from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup


process.jec = cms.ESSource('PoolDBESSource',
                           CondDBSetup,
                           connect = cms.string(SQLiteFile),
                           toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string(TagForAK4CHSJet),
            label  = cms.untracked.string('AK4PFchs')
            ),
        cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string(TagForAK4PuppiJet),
            label  = cms.untracked.string('AK4PFPuppi')
            ) 
        )
                           )

# Add an ESPrefer to override JEC that might be available from the global tag
if UseSQLiteFiles:
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')
else:
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'GlobalTag')


JERVersion=''

if "MC2018" in runEra:
    JERVersion='Autumn18_V7b_MC'
if "MC2017" in runEra:
    JERVersion='Fall17_V3b_MC'
if "MC2016" in runEra:
    JERVersion='Summer16_25nsV1b_MC'
if "MCUL2016" in runEra:
    JERVersion='Summer16_25nsV1b_MC'

if "Data2018" in runEra:
    JERVersion='Autumn18_V7b_DATA'
if "Data2017" in runEra:
    JERVersion='Fall17_V3b_DATA'
if "Data2016" in runEra or "DataUL2016" in runEra:
    JERVersion='Summer16_25nsV1b_DATA'



if "MCUL2017" in runEra:
    JERVersion='Summer19UL17_JRV3_MC'
if "DataUL2017" in runEra:
    JERVersion='Summer19UL17_JRV3_DATA'

if "MCUL2018" in runEra:
    JERVersion='Summer19UL18_JRV2_MC'
if "DataUL2018" in runEra:
    JERVersion='Summer19UL18_JRV2_DATA'

if "MCRun3" in runEra:
    JERVersion='Summer19UL18_JRV2_MC'
if "DataRun3" in runEra:
    JERVersion='Summer19UL18_JRV2_DATA'


if "MCUL2016APV" in runEra:
    JERVersion='Summer20UL16APV_JRV3_MC'
if "MCUL2016nonAPV" in runEra:
    JERVersion='Summer20UL16_JRV3_MC'

if "DataUL2016B" in runEra:
    JERVersion='Summer20UL16APV_JRV3_DATA'
if "DataUL2016C" in runEra:
    JERVersion='Summer20UL16APV_JRV3_DATA'
if "DataUL2016D" in runEra:
    JERVersion='Summer20UL16APV_JRV3_DATA'
if "DataUL2016E" in runEra:
    JERVersion='Summer20UL16APV_JRV3_DATA'
if "DataUL2016F" in runEra:
    JERVersion='Summer20UL16APV_JRV3_DATA'
if "DataUL2016Flate" in runEra:
    JERVersion='Summer20UL16_JRV3_DATA'
if "DataUL2016G" in runEra:
    JERVersion='Summer20UL16_JRV3_DATA'
if "DataUL2016H" in runEra:
    JERVersion='Summer20UL16_JRV3_DATA'


SQLiteFileJER='sqlite:'+JERVersion+'.db'


print(SQLiteFileJER)
process.jer = cms.ESSource("PoolDBESSource",
                           CondDBSetup,
                           toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('JetResolutionRcd'),
            tag    = cms.string('JR_'+JERVersion+'_PtResolution_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs_pt')
            ),
        cms.PSet(
            record = cms.string('JetResolutionScaleFactorRcd'),
            tag    = cms.string('JR_'+JERVersion+'_SF_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs')
            ),
#        cms.PSet(                                                                                                                                                                                            
#            record = cms.string('JetResolutionRcd'),                                                                                                                                                         
#            tag    = cms.string('JR_'+JERVersion+'_PtResolution_AK4PFPuppi'),
#            label  = cms.untracked.string('AK4PFPuppi_pt')                                                                                                                                                   
#            ),                                                                                                                                                                                               
#        cms.PSet(                                                                                                                                                                                            
#            record = cms.string('JetResolutionScaleFactorRcd'),                                                                                                                                              
#            tag    = cms.string('JR_'+JERVersion+'_SF_AK4PFPuppi'),                                                                                                                                    
#            label  = cms.untracked.string('AK4PFPuppi')                                                                                                                                                      
#            ),                                                                                                                                                                                               
        ),
                           connect = cms.string(SQLiteFileJER)
                           )

if UseSQLiteFiles: 
    process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')
else:
    process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'GlobalTag')


from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    labelName = 'UpdatedJEC',
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  
)

updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJetsPuppi'),
    labelName = 'UpdatedJEC',
    postfix = 'Puppi',
    jetCorrections = ('AK4PFPuppi', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None') 
)
process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.patJetCorrFactorsUpdatedJECPuppi * process.updatedPatJetsUpdatedJECPuppi)


#Recluster gen jets
## Filter out neutrinos from packed GenParticles 
#The filter on pdgid 2101, 2103, 2203 and 1103 should be harmless for standard samples. I added it because some private samples mistakenly added those unstable states as stable products. 
process.packedGenParticlesForJetsNoNuNEW = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16 && abs(pdgId) != 2101 &&abs(pdgId) != 2103 && abs(pdgId) != 2203  && abs(pdgId) != 1103 "))
process.packedGenParticlesForJetsWithNuNEW = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 2101 &&abs(pdgId) != 2103 && abs(pdgId) != 2203  && abs(pdgId) != 1103 "))
## Define GenJets 
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsNoNuNEW = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNuNEW')
process.ak4GenJetsWithNuNEW = ak4GenJets.clone(src = 'packedGenParticlesForJetsWithNuNEW')
#I didn't manage to create a new jet collection on top of MINIAOD with the matching to this updated gen jet collection   
#Work around: do the matching by hand                                                                                                                                                                                       
#Now redo the matching. The patJetGenJetMatch produces a matching between the gen jets and the reco jets. 
from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import patJetGenJetMatch

process.patJetGenJetMatchUpdate = patJetGenJetMatch.clone(
src         = cms.InputTag(CHSJetCollectionName),
matched     = cms.InputTag("ak4GenJetsNoNuNEW")
)
process.patJetGenJetMatchUpdatePuppi = patJetGenJetMatch.clone(
src         = cms.InputTag("slimmedJetsPuppi"),
matched     = cms.InputTag("ak4GenJetsNoNuNEW")
)
process.patJetGenWithNuJetMatchUpdate = patJetGenJetMatch.clone(
src         = cms.InputTag(CHSJetCollectionName),
matched     = cms.InputTag("ak4GenJetsWithNuNEW")
)
process.patJetGenWithNuJetMatchUpdatePuppi = patJetGenJetMatch.clone(
src         = cms.InputTag("slimmedJetsPuppi"),
matched     = cms.InputTag("ak4GenJetsWithNuNEW")
)
process.patJetGenJetMatchUpdateCalo = patJetGenJetMatch.clone(
src         = cms.InputTag("slimmedCaloJets"),
matched     = cms.InputTag("ak4GenJetsNoNuNEW")
)



#Update MET
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

#Recompute PFMET (with updated JECs)

if TheSkim == "MCJECs":

    runMetCorAndUncFromMiniAOD (
        process,
        isData = not ISMC
    )



'''
#Rerunning PUPPI (standard approach, default tuning)
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD( process, True );
#Set to false if you want to recompute PUPPI weights
process.puppiNoLep.useExistingWeights = True
process.puppi.useExistingWeights = True

#Recompute PUPPI MET (with updated JECs and possibly updated PUPPI weights)
runMetCorAndUncFromMiniAOD(process,
                           jetCollUnskimmed="updatedPatJetsUpdatedJECPuppi",
                           isData= not ISMC,
                           metType="Puppi",
                           postfix="Puppi",
                           jetFlavor="AK4PFPuppi",
                           reapplyJEC = False
                           )
'''



from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection

#Rerunning PUPPi with v15 tune
'''
process.ak4PuppiJets  = ak4PFJets.clone (src = 'puppi', doAreaFastjet = True, jetPtMin = 2.)
addJetCollection(process,labelName = 'Puppi', jetSource = cms.InputTag('ak4PuppiJets'), algo = 'AK', rParam=0.4, genJetCollection=cms.InputTag(GenJetCollectionName), jetCorrections = ('AK4PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'], 'None'),pfCandidates = cms.InputTag('packedPFCandidates'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 muSource =cms.InputTag( 'slimmedMuons'),
                 elSource = cms.InputTag('slimmedElectrons'),
                 genParticles= cms.InputTag('prunedGenParticles'),
                 getJetMCFlavour=ISMC
)

process.patJetsPuppi.addGenPartonMatch = cms.bool(ISMC)
process.patJetsPuppi.addGenJetMatch = cms.bool(ISMC)

from CommonTools.PileupAlgos.customizePuppiTune_cff import UpdatePuppiTuneV15

patAlgosToolsTask.add(process.ak4PuppiJets)
UpdatePuppiTuneV15(process,ISMC)
'''

#Now doing PF jets (not CHS). Those are not in MINIAOD
process.ak4PFJetsBis  = ak4PFJets.clone (src = 'packedPFCandidates', doAreaFastjet = True, jetPtMin = 2.)
patAlgosToolsTask.add(process.ak4PFJetsBis)
addJetCollection(process,labelName = 'Plain', jetSource = cms.InputTag('ak4PFJetsBis'), algo = 'AK', rParam=0.4, genJetCollection=cms.InputTag(GenJetCollectionName), jetCorrections = ('AK4PF', ['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'], 'None'),pfCandidates = cms.InputTag('packedPFCandidates'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 muSource =cms.InputTag( 'slimmedMuons'),
                 elSource = cms.InputTag('slimmedElectrons'),
                 genParticles= cms.InputTag('prunedGenParticles'),
                 getJetMCFlavour=ISMC
)
process.patJetsPlain.addGenPartonMatch = cms.bool(ISMC)
process.patJetsPlain.addGenJetMatch = cms.bool(ISMC)


#Recluster CHS jets to go lower in pt
if ReclusterCHSJets:
    process.ak4PFJetsCHSBis  = ak4PFJets.clone (src = 'pfCHS', doAreaFastjet = True, jetPtMin = 2.)
    patAlgosToolsTask.add(process.ak4PFJetsCHSBis)
    addJetCollection(process,labelName = 'CHS', jetSource = cms.InputTag('ak4PFJetsCHSBis'), algo = 'AK', rParam=0.4, genJetCollection=cms.InputTag(GenJetCollectionName), jetCorrections = ('AK4PFchs', ['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'], 'None'),pfCandidates = cms.InputTag('pfCHS'),
                     pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                     svSource = cms.InputTag('slimmedSecondaryVertices'),
                     muSource =cms.InputTag( 'slimmedMuons'),
                     elSource = cms.InputTag('slimmedElectrons'),
                     genParticles= cms.InputTag('prunedGenParticles'),
                     getJetMCFlavour=ISMC
                     )
    process.patJetsCHS.addGenPartonMatch = cms.bool(ISMC)
    process.patJetsCHS.addGenJetMatch = cms.bool(ISMC)

    #Now AK8 CHS
    from RecoJets.JetProducers.ak8PFJets_cfi import ak8PFJets
    process.ak8PFJetsCHSBis  = ak8PFJets.clone (src = 'pfCHS', doAreaFastjet = True, jetPtMin = 10.)
    patAlgosToolsTask.add(process.ak8PFJetsCHSBis)
    addJetCollection(process,labelName = 'AK8CHS', jetSource = cms.InputTag('ak8PFJetsCHSBis'), algo = 'AK', rParam=0.8, genJetCollection=cms.InputTag('slimmedGenJetsAK8'), jetCorrections = ('AK8PFchs', ['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'], 'None'),pfCandidates = cms.InputTag('pfCHS'),
                     pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                     svSource = cms.InputTag('slimmedSecondaryVertices'),
                     muSource =cms.InputTag( 'slimmedMuons'),
                     elSource = cms.InputTag('slimmedElectrons'),
                     genParticles= cms.InputTag('prunedGenParticles'),
                     getJetMCFlavour=ISMC
                     )
    process.patJetsAK8CHS.addGenPartonMatch = cms.bool(ISMC)
    process.patJetsAK8CHS.addGenJetMatch = cms.bool(ISMC)
#The produced collection is called selectedPatJetsAK8CHS 



#Recompute pile up ID
from RecoJets.JetProducers.PileupJetID_cfi import  _chsalgos_81x, _chsalgos_94x, _chsalgos_102x
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdate = process.pileupJetId.clone()
process.pileupJetIdUpdate.jets = cms.InputTag(CHSJetCollectionName)
process.pileupJetIdUpdate.inputIsCorrected = True
process.pileupJetIdUpdate.applyJec = False
process.pileupJetIdUpdate.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pileupJetIdUpdate.algos = cms.VPSet(_chsalgos_81x) 


process.pileupJetIdUpdate2017 = process.pileupJetId.clone()
process.pileupJetIdUpdate2017.jets = cms.InputTag(CHSJetCollectionName)
process.pileupJetIdUpdate2017.inputIsCorrected = True
process.pileupJetIdUpdate2017.applyJec = False
process.pileupJetIdUpdate2017.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pileupJetIdUpdate2017.algos = cms.VPSet(_chsalgos_94x)

process.pileupJetIdUpdate2018 = process.pileupJetId.clone()
process.pileupJetIdUpdate2018.jets = cms.InputTag(CHSJetCollectionName)
process.pileupJetIdUpdate2018.inputIsCorrected = True
process.pileupJetIdUpdate2018.applyJec = False
process.pileupJetIdUpdate2018.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pileupJetIdUpdate2018.algos = cms.VPSet(_chsalgos_102x)


#Compute QGL 
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets          = cms.InputTag(CHSJetCollectionName)
process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs')   


#

eraforEGMSmearing=''

if "UL2017" in runEra or "UL2018" in runEra:
    if "UL2017" in runEra:
        eraforEGMSmearing='2017-UL'
    if "UL2018" in runEra:
        eraforEGMSmearing='2018-UL'
    if "UL2016B" in runEra or "UL2016C" in runEra or "UL2016D" in runEra or "UL2016E" in runEra or "UL2016F" in runEra  or "UL2016APV" in runEra :
        eraforEGMSmearing='2016-UL-preVFP'
    if "UL2016Flate" in runEra or "UL2016G" in runEra or "UL2016H" in runEra  or "UL2016nonAPV" in runEra :
        eraforEGMSmearing='2016-UL-postVFP'
    from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
    setupEgammaPostRecoSeq(process,
                           runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                           era=eraforEGMSmearing)    
    process.ApplyEGMScaleSmearing=cms.Path(process.egammaPostRecoSeq)


from RecoMET.METFilters.primaryVertexFilter_cfi import primaryVertexFilter


#process.applyjecs =  cms.Path( process.jecSequence )
if ISMC and ReclusterGenJets: 
    process.reclustergenjets = cms.Path(process.packedGenParticlesForJetsNoNuNEW * process.packedGenParticlesForJetsWithNuNEW *process.ak4GenJetsNoNuNEW * process.ak4GenJetsWithNuNEW * process.patJetGenJetMatchUpdate *process.patJetGenJetMatchUpdatePuppi  * process.patJetGenWithNuJetMatchUpdate  * process.patJetGenWithNuJetMatchUpdatePuppi *process.patJetGenJetMatchUpdateCalo)


#You may want to comment out some of the following lines to speed things up
process.ApplyPatAlgos  = cms.Path(process.patAlgosToolsTask)

#process.rerunmetfilters = cms.Path( process.ecalBadCalibReducedMINIAOD2019Filter * process.ecalLaserCorrFilter * process.ecalDeadCellBoundaryEnergyFilterUpdate * process.BadChargedCandidateFilterUpdate ) 
#process.computepuid = cms.Path(process.pileupJetIdUpdate  * process.pileupJetIdUpdate2017 * process.pileupJetIdUpdate2018)
#process.computeqgl = cms.Path(process.QGTagger)

#This one obviously shouldn't be commented out
process.endpath = cms.EndPath( process.jmeanalyzer  )



