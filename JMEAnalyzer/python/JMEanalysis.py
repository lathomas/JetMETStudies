#TheSkim = "MCJECs"
TheSkim = "ZJetsResiduals"
TheSkim = "MCJECs" 
TheSkim = "HFJet"
TheSkim = "L1Unprefirable"
TheSkim = "ZJetsResiduals"
TheSkim = ""
#TheSkim = "L1Unprefirable"

ReclusterCHSJets = False
ReclusterGenJets = False


#runEra="DataUL2017F"
#runEra="MCUL2017"
runEra="MCUL2018"
#runEra="DataUL2018A"
ISMC=False
if "MC" in runEra:
    ISMC=True
if not ISMC:
    ReclusterGenJets = False

UseSQLiteFiles=False


if TheSkim == "MCJECs":
   ReclusterCHSJets = True
   ReclusterGenJets = True


import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )  
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(23374) )  

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)


process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(

#2018
'/store/mc/RunIISummer19UL18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/100000/00556E5A-E6CD-6C4B-83BA-0B5C00C0A5E4.root'
#2017
#'/store/mc/RunIISummer19UL17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/100000/002C691B-A0CE-A24F-8805-03B4C52C9004.root'
#2016
#'/store/mc/RunIISummer20UL16MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v13-v2/260000/003F1A76-9CDA-7644-A34E-923C4B1C0E5E.root'

#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/1266EB1A-9309-FD47-905B-DF9F4A3C4BFD.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/2AB745ED-157C-264C-988A-E6112679D39D.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/3EADBED7-D184-8446-B790-6287BB3C28F4.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/4EB946BE-9305-E741-8084-E805369734FF.root'

#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/00CD9760-02A6-8E41-9533-005920ED7062.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/0486D96E-5CC4-1545-A5CA-1FE0A79066D5.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/083B4F90-23D0-A44E-87BB-1E48570E3EA9.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/0AEE46BA-2522-8A41-89D4-FD141414A82E.root'

#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/4cores5k_106X_upgrade2018_realistic_v16_L1v1-v1/00000/1C829C47-045D-C249-9C17-FDBB6FC8C579.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/4cores5k_106X_upgrade2018_realistic_v16_L1v1-v1/00000/3EB73A23-BB61-FE4E-9103-505ADE78C08B.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/4cores5k_106X_upgrade2018_realistic_v16_L1v1-v1/00000/3FAD9883-07BF-394D-901E-0AF59D3B5557.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/4cores5k_106X_upgrade2018_realistic_v16_L1v1-v1/00000/45A3A67B-C8B4-AD4A-ADFE-7D3EFB8E0AE5.root'


#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/0AE62F14-AD3D-AD41-8F36-FCB332F3D861.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/0B03D729-090C-2E45-906A-226E6889E56D.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/0DAE62E9-D8D2-2041-943D-804CE3E02057.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/0EC3E8BD-4338-0646-9F99-1B22F328F2A9.root',
#'/store/mc/RunIISummer19UL18MiniAODv2/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/1225A64C-2381-D244-9BA9-3C4AE5D2D1B5.root'

#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0173DB9C-3E46-4741-B80B-9671E9E917CA.root'

#'/store/mc/RunIISummer20UL18MiniAODv2/VBF_HToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/110000/760A04EE-56EC-B14B-A9B4-BFCD87A61B7B.root'

#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0173DB9C-3E46-4741-B80B-9671E9E917CA.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0CB669FF-60C0-F84F-8165-6D75030DE776.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/2A5DB637-6AE8-874D-B4FA-B2CBBB079936.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/58324CFD-66FC-0746-B982-9E46F16D17E7.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/814DA621-8793-D343-AA2E-10A3B766CBA4.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/8F1158C5-5215-B14C-B6B6-BBBAD461A698.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/A02D8765-0B73-4A40-AEBE-2A966BCA0833.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/E8398719-E3DC-8743-9BC6-42FA554EE131.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/FA720C61-4D95-C046-9A21-3B6C6BC80AEE.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/270000/23266930-C608-104B-BB3B-64C66A052D1E.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/270000/6B41FEB4-0F68-8141-98D6-2FE8199B518B.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/270000/8D3F290E-CC30-8347-84E3-CE2103C9F516.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/270000/9E953C27-DBE5-EC4E-914E-A4D4470F829C.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/270000/C08A9B63-4FB6-354A-8C95-D1E647DD4239.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/006B8112-CA62-B245-B1E2-6B0B65146A38.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/0B7251EA-15C0-804D-A4DC-D7BE764218C2.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/39789092-8B9A-6646-9031-A52564DCC2FF.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/3EF1957A-6959-7F4D-A12F-E6780C4151E5.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/430EA7A9-7D28-C64F-949F-8D6B07BB20AC.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/44BA33AE-00CB-B444-BA38-77F18DD99ECB.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/46502098-8938-1647-B807-1B9ADEA80F44.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/4DB78BA2-1F7D-9842-B243-A127EF6B0944.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/4F0CC6B0-9745-024F-AF2F-125EC4D3893A.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/50CA582D-256F-494B-93A4-01302F66EF99.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/5AE89B6A-1B2C-AF47-A967-2E36626317B3.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/77B1AC88-F649-0441-93DC-C6682AEB6046.root',
#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/280000/7BC50BB6-5C9A-674A-B847-2D24EF59ACC0.root'

#'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/E8398719-E3DC-8743-9BC6-42FA554EE131.root'
#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_HighStat-v2/2580000/1026a87e-952c-4996-8ee7-57f5f0d63bb9.root',
#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_HighStat-v2/2580000/1fa066fd-adf9-440d-8f43-db3f68c28d2e.root',
#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_HighStat-v2/2580000/4fe44364-63b4-49a2-9577-38a6278cdf70.root',
#'/store/relval/CMSSW_12_1_0_pre4/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/PU_121X_mcRun3_2021_realistic_v10_HighStat-v2/2580000/9f84ce7c-602d-465f-a4eb-6ce431a00de0.root'

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
#process.TFileService = cms.Service("TFileService", fileName = cms.string("output_GJetsHT200to400.root") )
process.TFileService = cms.Service("TFileService", fileName = cms.string("output.root") )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag="102X_dataRun2_v8"

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")


process.load('Configuration.StandardSequences.Reconstruction_cff')

if ReclusterCHSJets: 
    CHSJetCollectionName ="selectedPatJetsCHS"
else: 
    CHSJetCollectionName ="updatedPatJetsUpdatedJEC"

if ReclusterGenJets:
    GenJetCollectionName="ak4GenJetsNoNuNEW"
else:
    GenJetCollectionName="slimmedGenJets"


EleVetoWP=''
EleTightWP=''
PhotonTightWP=''
#Rochester corrections folder: 
RochesterCorrectionFile="RochesterCorrections/"


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

print "Roch corr file: " 
print RochesterCorrectionFile

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
                                     PackedGenParticles=cms.InputTag("packedGenParticles"),
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
                                     JetPtCut=cms.double(20),
                                     AK8JetPtCut=cms.double(1000),
                                     ElectronPtCut=cms.double(1000),
                                     ElectronVetoWorkingPoint=cms.string(EleVetoWP),
                                     ElectronLooseWorkingPoint=cms.string(EleLooseWP),
                                     ElectronTightWorkingPoint=cms.string(EleTightWP),
                                     MuonPtCut=cms.double(10),
                                     RochCorrFile=cms.string(RochesterCorrectionFile),
                                     PhotonPtCut=cms.double(15),
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

myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
if not ISMC:
    process.source.lumisToProcess.extend(myLumis)
print "json" 
print JSONfile 

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


print SQLiteFileJER
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
#if TheSkim == "MCJECs":
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



process.applyjecs =  cms.Path( process.jecSequence )
if ISMC and ReclusterGenJets: 
    process.reclustergenjets = cms.Path(process.packedGenParticlesForJetsNoNuNEW * process.packedGenParticlesForJetsWithNuNEW *process.ak4GenJetsNoNuNEW * process.ak4GenJetsWithNuNEW * process.patJetGenJetMatchUpdate *process.patJetGenJetMatchUpdatePuppi  * process.patJetGenWithNuJetMatchUpdate  * process.patJetGenWithNuJetMatchUpdatePuppi *process.patJetGenJetMatchUpdateCalo)


#You may want to comment out some of the following lines to speed things up
process.ApplyPatAlgos  = cms.Path(process.patAlgosToolsTask)

process.rerunmetfilters = cms.Path( process.ecalBadCalibReducedMINIAOD2019Filter * process.ecalLaserCorrFilter * process.ecalDeadCellBoundaryEnergyFilterUpdate * process.BadChargedCandidateFilterUpdate ) 

process.computepuid = cms.Path(process.pileupJetIdUpdate  * process.pileupJetIdUpdate2017 * process.pileupJetIdUpdate2018)
process.computeqgl = cms.Path(process.QGTagger)

#This one obviously shouldn't be commented out
process.endpath = cms.EndPath( process.jmeanalyzer)



