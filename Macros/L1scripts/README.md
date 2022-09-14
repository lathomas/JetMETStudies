# Script for L1 studies


- Script to analyze ntuples and produce histos: 

```
python3 analyzel1.py  --maxEvents 100000 -c ZToMuMu
```

(relies on helper.py and Helper.h)

- Script to submit jobs on condor:

```
condor_submit   scriptcondorZToMuMu.sub
```
(the script itself is scriptcondor.sh)


- Script to draw histos: 
```
python3 drawplots.py -t efficiency --saveplot True -i input.root --den h_Qual12_plots_eta1p24to2p4 --num h_Qual12_plots_eta1p24to2p4_l1thrgeq10 h_Qual12_plots_eta1p24to2p4_l1thrgeq15 h_Qual12_plots_eta1p24to2p4_l1thrgeq20 h_Qual12_plots_eta1p24to2p4_l1thrgeq22 h_Qual12_plots_eta1p24to2p4_l1thrgeq25 --xtitle 'p_{T}^{#mu}(reco) (GeV)' --ytitle 'L1 Efficiency (Qual>=12)' --legend 'p_{T}^{L1}>10 GeV' 'p_{T}^{L1}>15 GeV' 'p_{T}^{L1}>20 GeV' 'p_{T}^{L1}>22 GeV' 'p_{T}^{L1}>25 GeV' --extralabel '#splitline{Z#rightarrow#mu#mu}{1.24<|#eta^{#mu}(reco)|<2.4}' --setlogx True --plotname L1Mu_TurnOnQual12_EMTF

```
