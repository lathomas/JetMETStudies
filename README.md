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
