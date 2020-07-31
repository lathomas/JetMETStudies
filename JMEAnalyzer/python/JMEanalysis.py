import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200000) )  
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(23374) )  

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)


process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
#                                    '/store/mc/PhaseIITDRSpring19MiniAOD/TTbar_14TeV_TuneCP5_Pythia8/MINIAODSIM/PU200_106X_upgrade2023_realistic_v3_ext1-v3/60000/E73765BC-F41D-6349-9394-2858D7CF81E1.root'
#'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext2-v1/50000/EE271696-E60F-1743-8DB2-AC50FD899528.root'
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/EA38FF5C-D444-E811-BFA6-001E6779262E.root'
#'/store/mc/RunIIAutumn18MiniAOD/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/20000/0E772C3C-0A8F-8341-886E-21A550D8E283.root'
#'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/90000/FD805C2C-92E0-4044-8061-770C77462EBC.root'
#'/store/mc/RunIIAutumn18MiniAOD/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/20000/435256CF-C2B1-4646-B28A-EE1A23DFD4DC.root'
#'/store/mc/RunIIAutumn18MiniAOD/TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/10000/014F2D77-CD31-8348-B6CC-B09134F7E8E4.root'
#'/store/mc/RunIISummer19UL17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/FCBBCC59-954C-3143-9634-DF8825F76369.root'
#'/store/data/Run2017F/SingleElectron/MINIAOD/09Aug2019_UL2017_rsb-v2/280000/114CE0B4-5CA0-FA47-BCE6-C13EE15A7191.root'

#'/store/data/Run2017F/ZeroBias/MINIAOD/09Aug2019_UL2017-v1/50000/CD375D55-3419-BD44-891D-0B3396DC9183.root',
#'/store/data/Run2017F/ZeroBias/MINIAOD/09Aug2019_UL2017-v1/50000/CAA1C836-1A53-1140-8DAE-C3C428EDB304.root',
#'/store/data/Run2017F/ZeroBias/MINIAOD/09Aug2019_UL2017-v1/50000/C93984AD-C2B3-084D-8735-FE533D3FC735.root',
#'/store/data/Run2017F/ZeroBias/MINIAOD/09Aug2019_UL2017-v1/50000/C7F52CEB-C7C1-6B45-B016-4731B8FF000E.root',
#'/store/data/Run2017F/ZeroBias/MINIAOD/09Aug2019_UL2017-v1/50000/C764B195-A6AD-984A-9A44-CA92289342DB.root',
#'/store/data/Run2017F/ZeroBias/MINIAOD/09Aug2019_UL2017-v1/50000/C722A73C-E67C-E445-BC9D-F287C2C962D1.root'


#'/store/data/Run2017F/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/FAAFF59A-2340-4641-88E0-8C468D6F56AD.root',
#'/store/data/Run2017F/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/F47C4B7D-B18A-3D4C-A6B9-601003D5AFF6.root',
#'/store/data/Run2017F/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/F4750D03-36C4-AD4A-8CAB-1501A5EE8EF3.root',
#'/store/data/Run2017F/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/F441F52A-9F68-BE42-B91B-DB4A8D1ABA99.root',
#'/store/data/Run2017F/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/F36E09AC-E7E6-EB43-B61A-A59227282B33.root',
#'/store/data/Run2017F/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/F18D0162-FF02-7E4A-BBB7-8D2A13CE7E87.root',
#'/store/data/Run2017F/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/F0E52585-2679-E542-BA08-EB22E63443A0.root',
#'/store/data/Run2017F/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/F0BDD302-F7CF-F94A-902B-8DFEF92CF3E8.root',
#'/store/data/Run2017F/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/F0496411-475D-9345-9A07-41BD1DD19F6D.root',
#'/store/data/Run2017F/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/F002BB55-D805-9B49-8A68-EFEB41E223EB.root',
#'/store/data/Run2017F/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/EF575753-1880-D047-92BC-DC8FFB908595.root',
#'/store/data/Run2017F/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/EE195096-4700-6149-AE6F-56BC9150D4E5.root'

'/store/data/Run2017F/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/3CD3E778-D5BF-394E-9D46-98D77A3CB58D.root'
#'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FC03A7B5-7FC0-E811-99BA-B496910A9A2C.root',
#'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FADA5029-2CC7-E811-BECD-0025907D1D6C.root',
#'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/F855E2C9-30C7-E811-A621-0025901AC0F8.root',
#'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/F6A1640F-65C7-E811-BA45-0025901F8740.root',
#'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/F6989592-7DC7-E811-96B7-0CC47A0AD3BC.root',
#'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/F6670B03-D5C4-E811-AAE4-008CFA1CBB34.root',
#'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/F640A4A9-55C7-E811-A17E-0025901AA5AE.root',
#'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/F6072CD4-40C7-E811-9E08-0CC47AD24CF8.root'
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/EC1F3053-63C3-6245-814E-74759808D0C7.root'


#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/230000/323B80E2-C415-694A-9A24-3E0F68A170D7.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/230000/3BAD0ABF-DB86-9B4C-BD26-D5DED96E54B2.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/230000/3BE12E39-1ADB-F444-ABE2-08E9F3765A11.root'
#'file:/user/lathomas/GamGamToLLStudy/CMSSW_10_6_13/src/output_met200_partial.root'

#'file:qcdht1000to1500_1.root',
#'file:qcdht1000to1500_2.root'
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/270000/000C0A31-7750-7A4D-A883-4676E3204108.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/270000/EDBBDDF4-C086-6E4C-AEC1-F95894C26521.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/270000/ED4BB7A9-8CCA-2848-AED9-071C95A97161.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/270000/945B5A22-C9E6-3D4E-A998-322431FD8759.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/270000/94007307-8399-1F4E-B59A-B32C70FC6375.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/270000/92BE7189-3636-3149-BFCD-91ED2498573C.root'


#'file:pickevents_ZeynepVBF.root'

#Here's a MINIAOD file on lxplus in case you want to use a local sample
#'file:/afs/cern.ch/work/l/lathomas/public/qcdht1000to1500_1.root'
        )
                            )

#process.TFileService = cms.Service("TFileService", fileName = cms.string("outputQCDHT1000to1500_puppiv16_200kevts.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("outputQCDHT1000to1500_puppiv16_23374evts.root") )
process.TFileService = cms.Service("TFileService", fileName = cms.string("output.root") )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag="102X_dataRun2_v8"

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")


process.load('Configuration.StandardSequences.Reconstruction_cff')


ISMC=False

runEra="DataUL2017F"
#runEra="MCUL2017"
UseSQLiteFiles=True


if "MC" in runEra:
    ISMC=True


EleVetoWP=''
EleTightWP=''
PhotonTightWP=''
#Rochester corrections folder: 
RochesterCorrectionFile="RochesterCorrections/"


if "Data2018" in runEra:
    if "2018D" in runEra: 
        process.GlobalTag.globaltag="102X_dataRun2_Prompt_v16" #2018D
    else:
        process.GlobalTag.globaltag="102X_dataRun2_v13" #2018ABC
    EleVetoWP='cutBasedElectronID-Fall17-94X-V1-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2018.txt"

if "Data2017" in runEra:
    process.GlobalTag.globaltag="102X_dataRun2_v13" #2017      
    EleVetoWP='cutBasedElectronID-Fall17-94X-V1-veto'
    EleTightWP='mvaEleID-Fall17-iso-V1-wp90'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2017.txt"

if "DataUL2017" in runEra:
    process.GlobalTag.globaltag="102X_dataRun2_v13" #2017      
    EleVetoWP='cutBasedElectronID-Fall17-94X-V1-veto'
    EleTightWP='mvaEleID-Fall17-iso-V1-wp90'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2017UL.txt"

if "Data2016" in runEra:
    process.GlobalTag.globaltag="102X_dataRun2_v13" #2016
    EleVetoWP='cutBasedElectronID-Fall17-94X-V1-veto'
    EleTightWP='mvaEleID-Fall17-iso-V1-wp90'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2016.txt"

if "MC2018" in runEra:
    process.GlobalTag.globaltag="102X_upgrade2018_realistic_v21" #2018     
    EleVetoWP='cutBasedElectronID-Fall17-94X-V1-veto'
    EleTightWP='mvaEleID-Fall17-iso-V2-wp90'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2018.txt"

if "MC2017" in runEra:
    process.GlobalTag.globaltag="102X_mc2017_realistic_v8" #2017
    EleVetoWP='cutBasedElectronID-Fall17-94X-V1-veto'
    EleTightWP='mvaEleID-Fall17-iso-V1-wp90'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2017.txt"

if "MC2016" in runEra:
    process.GlobalTag.globaltag="102X_mcRun2_asymptotic_v8" #2016
    EleVetoWP='cutBasedElectronID-Fall17-94X-V1-veto'
    EleTightWP='mvaEleID-Fall17-iso-V1-wp90'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2016.txt"

if "MCUL2017" in runEra:
    process.GlobalTag.globaltag="102X_mc2017_realistic_v8" #2017
    EleVetoWP='cutBasedElectronID-Fall17-94X-V1-veto'
    EleTightWP='mvaEleID-Fall17-iso-V1-wp90'
    PhotonTightWP='mvaPhoID-RunIIFall17-v1p1-wp80'
    RochesterCorrectionFile+="RoccoR2017UL.txt"

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
                                     Jets=cms.InputTag("updatedPatJetsUpdatedJEC"),
#                                     JetsPuppi=cms.InputTag("updatedPatJetsUpdatedJECPuppi"),
                                     JetsPuppi=cms.InputTag("slimmedJetsPuppi"),
                                     JetsPuppiAK8=cms.InputTag("slimmedJetsAK8"),
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
                                     PFMet=cms.InputTag("slimmedMETs"),
                                     PuppiMet=cms.InputTag("slimmedMETsPuppi"),
                                     Electrons=cms.InputTag("slimmedElectrons"),
                                     Muons=cms.InputTag("slimmedMuons"),
                                     Photons=cms.InputTag("slimmedPhotons"),
                                     JetPtCut=cms.double(20),
                                     AK8JetPtCut=cms.double(200),
                                     ElectronPtCut=cms.double(20),
                                     ElectronVetoWorkingPoint=cms.string(EleVetoWP),
                                     ElectronTightWorkingPoint=cms.string(EleTightWP),
                                     MuonPtCut=cms.double(20),
                                     RochCorrFile=cms.string(RochesterCorrectionFile),
                                     PhotonPtCut=cms.double(20),
                                     PhotonTightWorkingPoint=cms.string(PhotonTightWP),
                                     PFCandPtCut=cms.double(1000),
                                     SaveTree=cms.bool(True),
                                     IsMC=cms.bool(ISMC),
                                     SavePUIDVariables=cms.bool(True),
                                     SaveAK8Jets=cms.bool(True),
                                     DropUnmatchedJets=cms.bool(False),
                                     DropBadJets=cms.bool(False),
                                     ApplyPhotonID=cms.bool(False),
#                                     Skim=cms.string("ZToEEorMuMu"),
                                     Skim=cms.string(""),
#                                     Skim=cms.string("HighHT"),
                                     Debug=cms.bool(False)
                              )





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
    taggingMode    = cms.bool(True)
)



import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())


JSONfile =''

if "DataUL2017" in runEra:
    JSONfile = 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
if "Data2018" in runEra:
    JSONfile = 'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
if "Data2017" in runEra:
    JSONfile = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'
if "Data2016" in runEra:
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
if "MCUL2017" in runEra:
    JECsVersion='Summer19UL17_V5_MC'

if "Data2018" in runEra:
    JECsVersion='Autumn18_RunABCD_V19_DATA'
if "Data2017" in runEra:
    JECsVersion='Fall17_17Nov2017_V32_102X_DATA'
if "Data2016" in runEra:
    JECsVersion='Summer16_07Aug2017All_V11_DATA'

if "DataUL2017B" in runEra:
    JECsVersion='Summer19UL17_RunB_V5_DATA'
if "DataUL2017C" in runEra:
    JECsVersion='Summer19UL17_RunC_V5_DATA'
if "DataUL2017D" in runEra:
    JECsVersion='Summer19UL17_RunD_V5_DATA'
if "DataUL2017E" in runEra:
    JECsVersion='Summer19UL17_RunE_V5_DATA'
if "DataUL2017F" in runEra:
    JECsVersion='Summer19UL17_RunF_V5_DATA'



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

if "Data2018" in runEra:
    JERVersion='Autumn18_V7b_DATA'
if "Data2017" in runEra:
    JERVersion='Fall17_V3b_DATA'
if "Data2016" in runEra:
    JERVersion='Summer16_25nsV1b_DATA'

if "MCUL2017" in runEra:
    JERVersion='Summer19UL17_JRV2_MC'
if "DataUL2017" in runEra:
    JERVersion='Summer19UL17_JRV2_DATA'


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
src         = cms.InputTag("updatedPatJetsUpdatedJEC"),
matched     = cms.InputTag("ak4GenJetsNoNuNEW")
)
process.patJetGenJetMatchUpdatePuppi = patJetGenJetMatch.clone(
src         = cms.InputTag("slimmedJetsPuppi"),
matched     = cms.InputTag("ak4GenJetsNoNuNEW")
)
process.patJetGenWithNuJetMatchUpdate = patJetGenJetMatch.clone(
src         = cms.InputTag("updatedPatJetsUpdatedJEC"),
matched     = cms.InputTag("ak4GenJetsWithNuNEW")
)
process.patJetGenWithNuJetMatchUpdatePuppi = patJetGenJetMatch.clone(
src         = cms.InputTag("slimmedJetsPuppi"),
matched     = cms.InputTag("ak4GenJetsWithNuNEW")
)






#Update MET
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

#Recompute PFMET (with updated JECs)
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

#Rerunning PUPPi with v14-Chihuahua tune
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak4PuppiJets  = ak4PFJets.clone (src = 'puppi', doAreaFastjet = True, jetPtMin = 2.)

from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
addJetCollection(process,labelName = 'Puppi', jetSource = cms.InputTag('ak4PuppiJets'), algo = 'AK', rParam=0.4, genJetCollection=cms.InputTag('slimmedGenJets'), jetCorrections = ('AK4PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'], 'None'),pfCandidates = cms.InputTag('packedPFCandidates'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 muSource =cms.InputTag( 'slimmedMuons'),
                 elSource = cms.InputTag('slimmedElectrons'),
                 genParticles= cms.InputTag('prunedGenParticles'),
                 getJetMCFlavour=ISMC
)

process.patJetsPuppi.addGenPartonMatch = cms.bool(ISMC)
process.patJetsPuppi.addGenJetMatch = cms.bool(ISMC)

from CommonTools.PileupAlgos.customizePuppiTune_cff import UpdatePuppiTuneV14

patAlgosToolsTask.add(process.ak4PuppiJets)
UpdatePuppiTuneV14(process,ISMC)



#Recompute pile up ID
from RecoJets.JetProducers.PileupJetID_cfi import  _chsalgos_81x, _chsalgos_94x, _chsalgos_102x
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdate = process.pileupJetId.clone()
process.pileupJetIdUpdate.jets = cms.InputTag("updatedPatJetsUpdatedJEC")
process.pileupJetIdUpdate.inputIsCorrected = True
process.pileupJetIdUpdate.applyJec = False
process.pileupJetIdUpdate.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pileupJetIdUpdate.algos = cms.VPSet(_chsalgos_81x) 


process.pileupJetIdUpdate2017 = process.pileupJetId.clone()
process.pileupJetIdUpdate2017.jets = cms.InputTag("updatedPatJetsUpdatedJEC")
process.pileupJetIdUpdate2017.inputIsCorrected = True
process.pileupJetIdUpdate2017.applyJec = False
process.pileupJetIdUpdate2017.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pileupJetIdUpdate2017.algos = cms.VPSet(_chsalgos_94x)

process.pileupJetIdUpdate2018 = process.pileupJetId.clone()
process.pileupJetIdUpdate2018.jets = cms.InputTag("updatedPatJetsUpdatedJEC")
process.pileupJetIdUpdate2018.inputIsCorrected = True
process.pileupJetIdUpdate2018.applyJec = False
process.pileupJetIdUpdate2018.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pileupJetIdUpdate2018.algos = cms.VPSet(_chsalgos_102x)


#Compute QGL 
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets          = cms.InputTag('updatedPatJetsUpdatedJEC')
process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs')   


#

eraforEGMSmearing=''

if "UL2017" in runEra or "UL2018" in runEra:
    if "UL2017" in runEra:
        eraforEGMSmearing='2017-UL'
    if "UL2018" in runEra:
        eraforEGMSmearing='2018-UL'
    from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
    setupEgammaPostRecoSeq(process,
                           runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                           era=eraforEGMSmearing)    
    process.ApplyEGMScaleSmearing=cms.Path(process.egammaPostRecoSeq)



process.applyjecs =  cms.Path( process.jecSequence )
if ISMC: 
    process.reclustergenjets = cms.Path(process.packedGenParticlesForJetsNoNuNEW * process.packedGenParticlesForJetsWithNuNEW *process.ak4GenJetsNoNuNEW * process.ak4GenJetsWithNuNEW * process.patJetGenJetMatchUpdate *process.patJetGenJetMatchUpdatePuppi  * process.patJetGenWithNuJetMatchUpdate  * process.patJetGenWithNuJetMatchUpdatePuppi)


#You may want to comment out some of the following lines to speed things up
process.ApplyPatAlgos  = cms.Path(process.patAlgosToolsTask)

process.rerunmetfilters = cms.Path( process.ecalBadCalibReducedMINIAOD2019Filter * process.ecalLaserCorrFilter * process.ecalDeadCellBoundaryEnergyFilterUpdate * process.BadChargedCandidateFilterUpdate ) 
process.computepuid = cms.Path(process.pileupJetIdUpdate  * process.pileupJetIdUpdate2017 * process.pileupJetIdUpdate2018)
process.computeqgl = cms.Path(process.QGTagger)

#This one obviously shouldn't be commented out
process.endpath = cms.EndPath( process.jmeanalyzer)



