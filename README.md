# JetMETStudies

cmsrel CMSSW_10_6_3
cd CMSSW_10_6_3/src 
git cms-addpkg RecoMET/METFilters 
git clone https://github.com/lathomas/JetMETStudies.git
scram b -j4
cd JetMETStudies/JMEAnalyzer/python/
cmsRun JMEanalysis.py
