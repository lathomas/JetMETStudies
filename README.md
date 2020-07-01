# JetMETStudies

cmsrel CMSSW_10_6_13 <br>
cd CMSSW_10_6_13/src <br>
git cms-addpkg RecoMET/METFilters <br>
git clone https://github.com/lathomas/JetMETStudies.git -b feb2020<br>
scram b -j4 <br>
cd JetMETStudies/JMEAnalyzer/python/ <br>
cmsRun JMEanalysis.py <br>
<br>

Make sure you adapt the run era and the UseSQLiteFiles boolean in JMEanalysis.py. 
In general, when available, global tag should be preferred to SQLite files. For UL17, however, no global is currently available. 