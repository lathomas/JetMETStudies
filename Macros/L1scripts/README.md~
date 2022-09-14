# JetMETStudies

```
cmsrel CMSSW_12_4_8
cd CMSSW_12_4_8/src
cmsenv
git cms-addpkg RecoMET/METFilters
git clone https://github.com/lathomas/JetMETStudies.git 
scram b -j4
```



Make sure you adapt the **runEra** string and the **UseSQLiteFiles** boolean in JMEanalysis.py. 
In general, when available, global tag should be preferred to SQLite files. 


Currently there are two separate config files, one to run locally, one to run on crab. 
Make sure they keep in synch ! 


To run locally:

```
cmsRun JMEAnalyzer/python/JMEanalysis.py
```


One can send jobs through crabs by using the SubmitToCrab.sh script. For example: 
```
sh SubmitToCrab.sh /Muon/Run2022D-PromptReco-v2/MINIAOD muon_zmumu_2022dv2  1 L1Study_ZToMuMu DataRun3 
```

Arguments are (in that order):
- dataset name
- local folder for crab 
- files per job 
- skim name
- Label specifying era (for reapplying JECs etc)