# JetMETStudies

```
cmsrel CMSSW_10_6_25
cd CMSSW_10_6_25/src
cmsenv
git cms-addpkg RecoMET/METFilters
git clone https://github.com/lathomas/JetMETStudies.git 
scram b -j4
```

To get the EGM scaling/smearing correction (enabled by default), the following is also needed, as instructed in: <br>
https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaUL2016To2018
```
git cms-addpkg RecoEgamma/EgammaTools  ### essentially just checkout the package from CMSSW
git clone https://github.com/cms-egamma/EgammaPostRecoTools.git
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.
git clone -b ULSSfiles_correctScaleSysMC https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/
git cms-addpkg EgammaAnalysis/ElectronTools
scram b -j 4
```


Make sure you adapt the **runEra** string and the **UseSQLiteFiles** boolean in JMEanalysis.py. 
In general, when available, global tag should be preferred to SQLite files. 

You can make test with a local input file currently on lxplus: 
'file:/afs/cern.ch/work/l/lathomas/public/qcdht1000to1500_1.root'

For that file, please set
```
runEra=MCUL2017
```


One can send jobs through crabs by using the SubmitToCrab.sh. For example: 
```
sh SubmitToCrab.sh /DoubleMuon/Run2018A-UL2018_MiniAODv2-v1/MINIAOD doublemuon2018a 1 ZJetsResiduals DataUL2018A
sh SubmitToCrab.sh /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM dyamcatnlo2018  1 ZJetsResiduals MCUL2018
```

Arguments are (in that order):
- dataset name
- local folder for crab 
- files per job 
- skim name
- Label specifying era (for reapplying JECs etc)