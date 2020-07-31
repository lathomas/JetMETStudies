# JetMETStudies

```
cmsrel CMSSW_10_6_13
cd CMSSW_10_6_13/src
git cms-addpkg RecoMET/METFilters
git clone https://github.com/lathomas/JetMETStudies.git 
scram b -j4
cd JetMETStudies/JMEAnalyzer/python/
cmsRun JMEanalysis.py
```

To get the EGM scaling/smearing correction (enabled by default), the following is also needed, as instructed in: <br>
https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaUL2016To2018
```
git cms-merge-topic jainshilpi:ULV1_backport106X_forUsers
git clone https://github.com/jainshilpi/EgammaPostRecoTools.git -b ULV0  
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.
git clone https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git -b UL2017SSV2 EgammaAnalysis/ElectronTools/data/
scram b -j 8
```


Make sure you adapt the **runEra** string and the **UseSQLiteFiles** boolean in JMEanalysis.py. 
In general, when available, global tag should be preferred to SQLite files. For UL17, however, no global tag is currently available so one needs to rely on .db files. 

You can make test with a local input file currently on lxplus: 

For that file, please set
```
runEra=MCUL2017
```
