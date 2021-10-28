#define templateforJERC_cxx
#include "templateforJERC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TrueNvtxReweighting.h"
#include <iostream>


/*  !!!!!


Important caveat. Currently the filename argument in the loop method needs to include the following strings: 
For MC: 

ZtoMuMu/ZtoEE  => to pick the mumu/ee channel
MC => to tell the code this is MC
UL2017/UL2018  => to tell the code that this is UL reprocessing of a given year


For data: 
DoubleMuon/DoubleEG =>  to pick the mumu/ee channel 
UL2018A/UL2018B/... => to tell the code that this is UL reprocessing of a given run era

*/

//Jet pt bins
const int Nptbins = 20; 
double ptbins[Nptbins+1] =   {20,22,24,26,28,30,35,40,50,60,85,105,130,175,230,300,400,500,700,1000,1500};
const int Netabins = 32;

//Jet eta bins
double etabins[Netabins+1] = {-5.191, -4.013, -3.489, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.043, -1.930, -1.740, -1.305, -1.131, -0.783, -0.522, 0., 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 4.013, 5.191 };

//jet Pt threshold (corrected jet pt)
double ThePtMin  =12;

//Fraction of events to use. Use <1 for tests
double fractouse =1;

//Minimum Pt of the Z boson
double ptZmin = 0.;

//Label for plots splitted in bins of PU
TString pubin[3] = {"pult20","pu20to40","pugt40"};



//Option to make a second fit iteration (closure test)
bool seconditeration = false;

//Histogram to import existing customized JECs
TH2F * hjecs; 

//Apply or not some PU jet ID by default. Do not do that right now.
bool applypuid=false;


//For PU ID scale factor measurements (not related to JER)
const int nbinseta_forPUID =12 ;
double bineta_forPUID[nbinseta_forPUID+1]={-5,-3.0,-2.75,-2.5,-2.0,-1.479,0.,1.479,2.0,2.5,2.75,3.0,5.0};

const int nbinspt_forPUID =5 ;
double binpt_forPUID[nbinspt_forPUID+1]={15,20,25,30,40,50};

 
/*

To do: better PU template ! 

For each zpt bin, generate 10000 values for the relevant zpt range and then 10000k values based on  h_PtRecoJetPU

generate random jetpt from PU with phi(jet) = 3.14159+Zphi

Loop (possibly several times) over events, add that jet to the jet list and run the standard code. ( in particular, check that njetspt20eta2p4 <=1 and njetspt30<=1 (including that jet) ) 
If that jet is not the selected jet (rare but can happen if it has pt < 20/30) => continue

Maybe it's good enough to make the following approximation:
 check that njetspt20eta2p4 ==0 and njetspt30 ==0


 Maybe one could also require no additional jet passing PU ID with pt>15 with |eta|<2.5? 


One more possible thing: 
compute a MC correction factor by dividing h_PtRecoJetoverPtRecoZ_NoGenMatchedJet and h_PtRecoJetoverPtRecoZ_PU (using a wide eta range: barrel, endcap1, endcap2, hf)

Split low pt bins and go down to 15 GeV?
 */




void templateforJERC::Loop(TString samplename, bool useresiduals,  TString unctytype)
{
  
  //Only relevant if you apply custom JECs
  TFile * filewithjecs = new TFile("resultsjecszjets.root","open");
  hjecs = (TH2F*) filewithjecs->Get("mean_diff");
  
  
  
  //PUID stuff (not relevant for JER)
  TH1F * h_PUMVAID_PUjets[nbinseta_forPUID][nbinspt_forPUID]; 
  TH1F * h_PUMVAID_Realjets[nbinseta_forPUID][nbinspt_forPUID]; 
  TH1F * h_PUMVAID_RealjetsPUTemp[nbinseta_forPUID][nbinspt_forPUID]; 
  TH1F * h_PUMVAID_PUgenjets[nbinseta_forPUID][nbinspt_forPUID]; 
  TH1F * h_PUMVAID_Realgenjets[nbinseta_forPUID][nbinspt_forPUID]; 
  TH1F * h_dphiZjets_L[nbinseta_forPUID][nbinspt_forPUID][3]; 
  TH1F * h_dphiZjets_M[nbinseta_forPUID][nbinspt_forPUID][3];
  TH1F * h_dphiZjets_T[nbinseta_forPUID][nbinspt_forPUID][3];
  TH1F * h_dphiZjets_ALL[nbinseta_forPUID][nbinspt_forPUID][3];
  TH1F * h_dphiZjets_genmatched_L[nbinseta_forPUID][nbinspt_forPUID][3]; 
  TH1F * h_dphiZjets_genmatched_M[nbinseta_forPUID][nbinspt_forPUID][3];
  TH1F * h_dphiZjets_genmatched_T[nbinseta_forPUID][nbinspt_forPUID][3];
  TH1F * h_dphiZjets_genmatched_ALL[nbinseta_forPUID][nbinspt_forPUID][3];
  TH1F * h_dphiZjets_genunmatched_L[nbinseta_forPUID][nbinspt_forPUID][3]; 
  TH1F * h_dphiZjets_genunmatched_M[nbinseta_forPUID][nbinspt_forPUID][3];
  TH1F * h_dphiZjets_genunmatched_T[nbinseta_forPUID][nbinspt_forPUID][3];
  TH1F * h_dphiZjets_genunmatched_ALL[nbinseta_forPUID][nbinspt_forPUID][3];
  
  for(int i = 0; i < nbinseta_forPUID ; i++){
    for(int j = 0; j < nbinspt_forPUID ; j++){
      TString suffix = "eta_" +(TString)Form("%4.3f",bineta_forPUID[i]) + "to"+(TString)Form("%4.3f",bineta_forPUID[i+1])  +  "pt_" +(TString)Form("%2.0f",binpt_forPUID[j]) + "to"+(TString)Form("%2.0f",binpt_forPUID[j+1]) ;
      suffix.ReplaceAll(".","p");
      suffix.ReplaceAll("-","min");
      h_PUMVAID_PUjets[i][j]= new TH1F("h_PUMVAID_PUjets"+suffix,"",100,-1,1);
      h_PUMVAID_PUjets[i][j]->Sumw2();
      h_PUMVAID_Realjets[i][j]= new TH1F("h_PUMVAID_Realjets"+suffix,"",100,-1,1);
      h_PUMVAID_Realjets[i][j]->Sumw2();
      h_PUMVAID_RealjetsPUTemp[i][j]= new TH1F("h_PUMVAID_RealjetsPUTemp"+suffix,"",100,-1,1);
      h_PUMVAID_RealjetsPUTemp[i][j]->Sumw2();
      h_PUMVAID_PUgenjets[i][j]= new TH1F("h_PUMVAID_PUgenjets"+suffix,"",100,-1,1);
      h_PUMVAID_PUgenjets[i][j]->Sumw2();
      h_PUMVAID_Realgenjets[i][j]= new TH1F("h_PUMVAID_Realgenjets"+suffix,"",100,-1,1);
      h_PUMVAID_Realgenjets[i][j]->Sumw2();
      
      for(int k=0; k<3;k++){
	TString suffixbal = suffix;
	if (k ==0) suffixbal += "_goodbalance" ;
	if (k ==1) suffixbal += "_baddbalance" ;
	
	h_dphiZjets_L[i][j][k]= new TH1F("h_dphiZjets_L"+suffixbal,"",100,0,2);; h_dphiZjets_L[i][j][k]->Sumw2();
	h_dphiZjets_M[i][j][k]= new TH1F("h_dphiZjets_M"+suffixbal,"",100,0,2);;h_dphiZjets_M[i][j][k]->Sumw2();
	h_dphiZjets_T[i][j][k]= new TH1F("h_dphiZjets_T"+suffixbal,"",100,0,2);;h_dphiZjets_T[i][j][k]->Sumw2();
	h_dphiZjets_ALL[i][j][k]= new TH1F("h_dphiZjets_ALL"+suffixbal,"",100,0,2);;h_dphiZjets_ALL[i][j][k]->Sumw2();
	h_dphiZjets_genmatched_L[i][j][k]= new TH1F("h_dphiZjets_genmatched_L"+suffixbal,"",100,0,2);; h_dphiZjets_genmatched_L[i][j][k]->Sumw2();
	h_dphiZjets_genmatched_M[i][j][k]= new TH1F("h_dphiZjets_genmatched_M"+suffixbal,"",100,0,2);;h_dphiZjets_genmatched_M[i][j][k]->Sumw2();
	h_dphiZjets_genmatched_T[i][j][k]= new TH1F("h_dphiZjets_genmatched_T"+suffixbal,"",100,0,2);;h_dphiZjets_genmatched_T[i][j][k]->Sumw2();
	h_dphiZjets_genmatched_ALL[i][j][k]= new TH1F("h_dphiZjets_genmatched_ALL"+suffixbal,"",100,0,2);;h_dphiZjets_genmatched_ALL[i][j][k]->Sumw2();
	h_dphiZjets_genunmatched_L[i][j][k]= new TH1F("h_dphiZjets_genunmatched_L"+suffixbal,"",100,0,2);; h_dphiZjets_genunmatched_L[i][j][k]->Sumw2();
	h_dphiZjets_genunmatched_M[i][j][k]= new TH1F("h_dphiZjets_genunmatched_M"+suffixbal,"",100,0,2);;h_dphiZjets_genunmatched_M[i][j][k]->Sumw2();
	h_dphiZjets_genunmatched_T[i][j][k]= new TH1F("h_dphiZjets_genunmatched_T"+suffixbal,"",100,0,2);;h_dphiZjets_genunmatched_T[i][j][k]->Sumw2();
	h_dphiZjets_genunmatched_ALL[i][j][k]= new TH1F("h_dphiZjets_genunmatched_ALL"+suffixbal,"",100,0,2);;h_dphiZjets_genunmatched_ALL[i][j][k]->Sumw2();
      }
    }
  }
  
  TH2F *  h_PUID_T_PUjets = new TH2F("h_PUID_T_PUjets","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_M_PUjets = new TH2F("h_PUID_M_PUjets","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_L_PUjets = new TH2F("h_PUID_L_PUjets","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_PUjetsden = new TH2F("h_PUID_PUjetsden","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_T_Realjets = new TH2F("h_PUID_T_Realjets","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_M_Realjets = new TH2F("h_PUID_M_Realjets","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_L_Realjets = new TH2F("h_PUID_L_Realjets","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_Realjetsden = new TH2F("h_PUID_Realjetsden","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  
  TH2F *  h_PUID_T_RealjetsPUTemp = new TH2F("h_PUID_T_RealjetsPUTemp","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_M_RealjetsPUTemp = new TH2F("h_PUID_M_RealjetsPUTemp","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_L_RealjetsPUTemp = new TH2F("h_PUID_L_RealjetsPUTemp","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_RealjetsPUTempden = new TH2F("h_PUID_RealjetsPUTempden","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  
  TH2F *  h_PUID_T_PUgenjets = new TH2F("h_PUID_T_PUgenjets","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_M_PUgenjets = new TH2F("h_PUID_M_PUgenjets","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_L_PUgenjets = new TH2F("h_PUID_L_PUgenjets","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_PUgenjetsden = new TH2F("h_PUID_PUgenjetsden","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_T_Realgenjets = new TH2F("h_PUID_T_Realgenjets","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_M_Realgenjets = new TH2F("h_PUID_M_Realgenjets","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_L_Realgenjets = new TH2F("h_PUID_L_Realgenjets","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
  TH2F *  h_PUID_Realgenjetsden = new TH2F("h_PUID_Realgenjetsden","",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);

  h_PUID_M_PUjets->Sumw2();
  h_PUID_L_PUjets->Sumw2();
  h_PUID_PUjetsden->Sumw2();
  h_PUID_T_Realjets->Sumw2();
  h_PUID_M_Realjets->Sumw2();
  h_PUID_L_Realjets->Sumw2();
  h_PUID_Realjetsden->Sumw2();

  h_PUID_T_RealjetsPUTemp->Sumw2();
  h_PUID_M_RealjetsPUTemp->Sumw2();
  h_PUID_L_RealjetsPUTemp->Sumw2();
  h_PUID_RealjetsPUTempden->Sumw2();

  h_PUID_T_PUgenjets->Sumw2();
  h_PUID_M_PUgenjets->Sumw2();
  h_PUID_L_PUgenjets->Sumw2();
  h_PUID_PUgenjetsden->Sumw2();
  h_PUID_T_Realgenjets->Sumw2();
  h_PUID_M_Realgenjets->Sumw2();
  h_PUID_L_Realgenjets->Sumw2();
  h_PUID_Realgenjetsden->Sumw2();
  

  const  int nbinsfornvtx = 12;
  int nvtxbins[nbinsfornvtx+1] = {0,5,10,15,20,25,30,35,40,45,50,55,60};

  TH1F *  h_nvtx=new TH1F("h_nvtx","",100,0,100);
  TH1F *  h_nvtx_nvtxbinned[nbinsfornvtx];
  TH2F *  h_PUID_T_PUjets_nvtxbinned[nbinsfornvtx];
  TH2F *  h_PUID_M_PUjets_nvtxbinned[nbinsfornvtx];
  TH2F *  h_PUID_L_PUjets_nvtxbinned[nbinsfornvtx];
  TH2F *  h_PUID_PUjetsden_nvtxbinned[nbinsfornvtx];

  for(int i = 0; i<nbinsfornvtx;i++){ 
    TString suffix = "_nvtx"+(TString)Form("%d",nvtxbins[i])+"to"+(TString)Form("%d",nvtxbins[i+1]);
    h_PUID_T_PUjets_nvtxbinned[i]= new TH2F("h_PUID_T_PUjets"+suffix,"",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
    h_PUID_M_PUjets_nvtxbinned[i]= new TH2F("h_PUID_M_PUjets"+suffix,"",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
    h_PUID_L_PUjets_nvtxbinned[i]= new TH2F("h_PUID_L_PUjets"+suffix,"",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
    h_PUID_PUjetsden_nvtxbinned[i]= new TH2F("h_PUID_PUjetsden"+suffix,"",nbinspt_forPUID,binpt_forPUID,nbinseta_forPUID,bineta_forPUID);
    h_nvtx_nvtxbinned[i]=new TH1F("h_nvtx"+suffix,"",100,0,100);

    h_PUID_T_PUjets_nvtxbinned[i]->Sumw2();
    h_PUID_M_PUjets_nvtxbinned[i]->Sumw2();
    h_PUID_L_PUjets_nvtxbinned[i]->Sumw2();
    h_PUID_PUjetsden_nvtxbinned[i]->Sumw2();
    h_nvtx_nvtxbinned[i]->Sumw2();
  }
  
  //All the above is related to PU ID SF (not JER SF)


  //Input samples
  TString fname ="";
  
  if(samplename.Index("MCUL2018MadGraph")>=0)fname = "/pnfs/iihe/cms/store/user/lathomas/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAOD_106X_upgrade2018_realistic_v11_L1v1_v1_MCUL2018_SkimDileptonJERC_Sept2020/200930_150455/TOTAL.root";
  //else if(samplename.Index("MCUL2018")>=0)fname = "/pnfs/iihe/cms/store/user/lathomas/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer19UL18MiniAOD_106X_upgrade2018_realistic_v11_L1v1_v2_MCUL2018_SkimDileptonJERC_Sept2020/200922_235006/TOTAL.root";
  else if(samplename.Index("MCUL2018")>=0)fname = "/pnfs/iihe/cms/store/user/lathomas/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2_106X_upgrade2018_realistic_v16_L1v1_v2_MCUL2018_ZJetsResiduals_May2021/211015_002042/TOTAL.root";
  else if(samplename.Index("MCUL2017")>=0)fname = "/pnfs/iihe/cms/store/user/lathomas/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer19UL17MiniAOD_106X_mc2017_realistic_v6_v2_MCUL2017_SkimDileptonJERC_June2020/200709_104206/TOTAL.root";
  else if(samplename.Index("MC2016")>=0)fname = "/pnfs/iihe/cms/store/user/lathomas/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3_PUMoriond17_94X_mcRun2_asymptotic_v3_ext2_v1_MC2016_SkimDileptonJERC_Nov2019_Autumnv18_18b/191101_193724/TOTAL.root";
  else if(samplename.Index("MC2017")>=0)fname = "/pnfs/iihe/cms/store/user/lathomas/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1_v1_MC2017_SkimDileptonJERC_Nov2019_Autumnv18_18b/191101_193602/TOTAL.root";
  else if(samplename.Index("MC2018")>=0)fname = "/pnfs/iihe/cms/store/user/lathomas/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD_102X_upgrade2018_realistic_v15_ext2_v1_MC2018_SkimDileptonJERC_Nov2019_Autumnv18_18b/191101_193444/TOTAL.root";
  
  
  if(samplename=="DoubleEG2016B") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2016B_17Jul2018_ver2_v1_Run2016B_SkimDileptonJERC_July2019_LowJetPt/190802_100107/TOTAL.root"; 
  if(samplename=="DoubleEG2016C") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2016C_17Jul2018_v1_Run2016C_SkimDileptonJERC_July2019_LowJetPt/190802_101120/TOTAL.root";
  if(samplename=="DoubleEG2016D") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2016D_17Jul2018_v1_Run2016D_SkimDileptonJERC_July2019_LowJetPt/190802_101246/TOTAL.root";
  if(samplename=="DoubleEG2016E") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2016E_17Jul2018_v1_Run2016E_SkimDileptonJERC_July2019_LowJetPt/190802_101407/TOTAL.root";
  if(samplename=="DoubleEG2016F") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2016F_17Jul2018_v1_Run2016F_SkimDileptonJERC_July2019_LowJetPt/190802_101534/TOTAL.root";
  if(samplename=="DoubleEG2016G") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2016G_17Jul2018_v1_Run2016G_SkimDileptonJERC_July2019_LowJetPt/190802_101658/TOTAL.root";
  if(samplename=="DoubleEG2016H") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2016H_17Jul2018_v1_Run2016H_SkimDileptonJERC_July2019_LowJetPt/190802_101823/TOTAL.root";
  
  if(samplename=="DoubleEG2017B") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2017B_31Mar2018_v1_Run2017B_SkimDileptonJERC_July2019_LowJetPt/190802_083309/TOTAL.root";
  if(samplename=="DoubleEG2017C") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2017C_31Mar2018_v1_Run2017C_SkimDileptonJERC_July2019_LowJetPt/190802_085818/TOTAL.root";
  if(samplename=="DoubleEG2017D") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2017D_31Mar2018_v1_Run2017D_SkimDileptonJERC_July2019_LowJetPt/190802_084205/TOTAL.root";
  if(samplename=="DoubleEG2017E") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2017E_31Mar2018_v1_Run2017E_SkimDileptonJERC_July2019_LowJetPt/190802_084328/TOTAL.root";
  if(samplename=="DoubleEG2017F") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2017F_31Mar2018_v1_Run2017F_SkimDileptonJERC_July2019_LowJetPt/190802_084625/TOTAL.root";
  
  if(samplename=="DoubleMuon2016B") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2016B_17Jul2018_ver2_v1_Run2016B_SkimDileptonJERC_July2019_LowJetPt/190802_101959/TOTAL.root";
  if(samplename=="DoubleMuon2016C") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2016C_17Jul2018_v1_Run2016C_SkimDileptonJERC_July2019_LowJetPt/190802_102119/TOTAL.root";
  if(samplename=="DoubleMuon2016D") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2016D_17Jul2018_v1_Run2016D_SkimDileptonJERC_July2019_LowJetPt/190802_102241/TOTAL.root";
  if(samplename=="DoubleMuon2016E") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2016E_17Jul2018_v1_Run2016E_SkimDileptonJERC_July2019_LowJetPt/190802_102403/TOTAL.root";
  if(samplename=="DoubleMuon2016F") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2016F_17Jul2018_v1_Run2016F_SkimDileptonJERC_July2019_LowJetPt/190802_102522/TOTAL.root";
  if(samplename=="DoubleMuon2016G") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2016G_17Jul2018_v1_Run2016G_SkimDileptonJERC_July2019_LowJetPt/190802_102645/TOTAL.root";
  if(samplename=="DoubleMuon2016H") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2016H_17Jul2018_v1_Run2016H_SkimDileptonJERC_July2019_LowJetPt/190802_102813/TOTAL.root";
  
  if(samplename=="DoubleMuon2017B") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2017B_31Mar2018_v1_Run2017B_SkimDileptonJERC_July2019_LowJetPt/190802_091127/TOTAL.root";
  if(samplename=="DoubleMuon2017C") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2017C_31Mar2018_v1_Run2017C_SkimDileptonJERC_July2019_LowJetPt/190802_091248/TOTAL.root";
  if(samplename=="DoubleMuon2017D") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2017D_31Mar2018_v1_Run2017D_SkimDileptonJERC_July2019_LowJetPt/190802_091405/TOTAL.root";
  if(samplename=="DoubleMuon2017E") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2017E_31Mar2018_v1_Run2017E_SkimDileptonJERC_July2019_LowJetPt/190802_091525/TOTAL.root";
  if(samplename=="DoubleMuon2017F") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2017F_31Mar2018_v1_Run2017F_SkimDileptonJERC_July2019_LowJetPt/190802_091640/TOTAL.root";

  //UL
  if(samplename=="DoubleEGUL2017B") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2017B_09Aug2019_UL2017_v1_RunUL2017B_SkimDileptonJERC_June2020/200710_085859/TOTAL.root";
  if(samplename=="DoubleEGUL2017C") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2017C_09Aug2019_UL2017_v1_RunUL2017C_SkimDileptonJERC_June2020/200710_090004/TOTAL.root";
  if(samplename=="DoubleEGUL2017D") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2017D_09Aug2019_UL2017_v1_RunUL2017D_SkimDileptonJERC_June2020/200802_071846/TOTAL.root";
  if(samplename=="DoubleEGUL2017E") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2017E_09Aug2019_UL2017_v1_RunUL2017E_SkimDileptonJERC_June2020/200710_090212/TOTAL.root";
  if(samplename=="DoubleEGUL2017F") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleEG/Run2017F_09Aug2019_UL2017_v1_RunUL2017F_SkimDileptonJERC_June2020/200710_090315/TOTAL.root";

  if(samplename=="DoubleMuonUL2017B") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2017B_09Aug2019_UL2017_v1_RunUL2017B_SkimDileptonJERC_June2020/200710_090416/TOTAL.root";
  if(samplename=="DoubleMuonUL2017C") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2017C_09Aug2019_UL2017_v1_RunUL2017C_SkimDileptonJERC_June2020/200710_090519/TOTAL.root";
  if(samplename=="DoubleMuonUL2017D") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2017D_09Aug2019_UL2017_v1_RunUL2017D_SkimDileptonJERC_June2020/200710_090621/TOTAL.root";
  if(samplename=="DoubleMuonUL2017E") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2017E_09Aug2019_UL2017_v1_RunUL2017E_SkimDileptonJERC_June2020/200710_090728/TOTAL.root";
  if(samplename=="DoubleMuonUL2017F") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2017F_09Aug2019_UL2017_v1_RunUL2017F_SkimDileptonJERC_June2020/200710_090831/TOTAL.root"; 

  if(samplename=="DoubleMuon2018A") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2018A_17Sep2018_v2_Run2018A_SkimDileptonJERC_July2019_LowJetPt/190731_150758/TOTAL.root";
  if(samplename=="DoubleMuon2018B") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2018B_17Sep2018_v1_Run2018B_SkimDileptonJERC_July2019_LowJetPt/190731_151447/TOTAL.root";
  if(samplename=="DoubleMuon2018C") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2018C_17Sep2018_v1_Run2018C_SkimDileptonJERC_July2019_LowJetPt/190731_151619/TOTAL.root";
  if(samplename=="DoubleMuon2018D") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2018D_PromptReco_v2_Run2018D_SkimDileptonJERC_July2019_LowJetPt/190731_151751/TOTAL.root";

  if(samplename=="DoubleEG2018A") fname ="/pnfs/iihe/cms/store/user/lathomas/EGamma/Run2018A_17Sep2018_v2_Run2018A_SkimDileptonJERC_July2019_LowJetPt/190731_151920/TOTAL.root";
  if(samplename=="DoubleEG2018B") fname ="/pnfs/iihe/cms/store/user/lathomas/EGamma/Run2018B_17Sep2018_v1_Run2018B_SkimDileptonJERC_July2019_LowJetPt/190731_152048/TOTAL.root";
  if(samplename=="DoubleEG2018C") fname ="/pnfs/iihe/cms/store/user/lathomas/EGamma/Run2018C_17Sep2018_v1_Run2018C_SkimDileptonJERC_July2019_LowJetPt/190731_152212/TOTAL.root";
  if(samplename=="DoubleEG2018D") fname ="/pnfs/iihe/cms/store/user/lathomas/EGamma/Run2018D_PromptReco_v2_Run2018D_SkimDileptonJERC_July2019_LowJetPt/190731_152332/TOTAL.root";

  //if(samplename=="DoubleMuonUL2018A") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2018A_12Nov2019_UL2018_v2_RunUL2018A_SkimDileptonJERC_Sept2020/200922_233951/TOTAL.root";
  if(samplename=="DoubleMuonUL2018A") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2018A_UL2018_MiniAODv2_v1_DataUL2018A_ZJetsResiduals_May2021/211014_235158/TOTAL.root";
  if(samplename=="DoubleMuonUL2018B") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2018B_12Nov2019_UL2018_v2_RunUL2018B_SkimDileptonJERC_Sept2020/200922_234138/TOTAL.root";
  if(samplename=="DoubleMuonUL2018C") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2018C_12Nov2019_UL2018_v2_RunUL2018C_SkimDileptonJERC_Sept2020/200922_234250/TOTAL.root";
  if(samplename=="DoubleMuonUL2018D") fname ="/pnfs/iihe/cms/store/user/lathomas/DoubleMuon/Run2018D_12Nov2019_UL2018_v3_RunUL2018D_SkimDileptonJERC_Sept2020/200922_234404/TOTAL.root";
 
  if(samplename=="DoubleEGUL2018A") fname ="/pnfs/iihe/cms/store/user/lathomas/EGamma/Run2018A_12Nov2019_UL2018_v2_RunUL2018A_SkimDileptonJERC_Sept2020/200922_234519/TOTAL.root";
  if(samplename=="DoubleEGUL2018B") fname ="/pnfs/iihe/cms/store/user/lathomas/EGamma/Run2018B_12Nov2019_UL2018_v2_RunUL2018B_SkimDileptonJERC_Sept2020/200922_234630/TOTAL.root";
  if(samplename=="DoubleEGUL2018C") fname ="/pnfs/iihe/cms/store/user/lathomas/EGamma/Run2018C_12Nov2019_UL2018_v2_RunUL2018C_SkimDileptonJERC_Sept2020/200922_234741/TOTAL.root";
  if(samplename=="DoubleEGUL2018D") fname ="/pnfs/iihe/cms/store/user/lathomas/EGamma/Run2018D_12Nov2019_UL2018_v4_RunUL2018D_SkimDileptonJERC_Sept2020/200922_234853/TOTAL.root";
  
  if(samplename=="ZtoMuMuTestSoumyaMCUL2017") fname ="/user/lathomas/L1Prefiring/CMSSW_10_6_25/src/JetMETStudies/JMEAnalyzer/python/smalltest.root";
  TFile * myf = new TFile(fname,"open");
  
  //Determine whether this is MC or data
  bool isMC = true;
  if(fname.Index("DoubleMuon")>=0||fname.Index("DoubleEG")>=0 ||fname.Index("EGamma")>=0) isMC=false;
  TString residualstring = useresiduals ? "" : "_recalcresiduals";
  if(isMC) seconditeration = false; 
  //Some labels for the ouput file name (specifying what options were used)
  TString iter = seconditeration ?  "2ndit" :"";
  TString puidstr = applypuid ? "_pumvageq0_":"";
  TString outname = "output_"+samplename+"_"+residualstring+iter+"uncty"+unctytype+puidstr+"Oct2021.root";

  TFile * outf = new TFile(outname, "recreate");
  TTree * thetree = samplename=="ZtoMuMuTestSoumyaMCUL2017" || samplename == "DoubleMuonUL2018A" ||samplename =="MCUL2018ZtoMuMu" ? (TTree*) myf->Get("jmeanalyzer/tree")  : (TTree*) myf->Get("FakeElectrons/fakeTree");
  

  //Now some histos for pt balance (=pt (jet)/pt(Z)):
  //This is for JER SF 
  
  TH1D *  h_PtGenJetoverPtGenZ[Nptbins][Netabins][3];  //Pt balance at gen level
  TH1D *  h_PtRecoZoverPtGenZ[Nptbins][Netabins][3];  //Pt balance at reco level: 
  TH1D *  h_PtRecoJetoverPtRecoZ_testsample[Nptbins][Netabins][3];  //Pt balance at reco level: test sample (complementary to events used in the histo h_PtRecoZoverPtGenZ
  TH1D *  h_PtRecoJetoverPtRecoZ_NoGenMatchedJet[Nptbins][Netabins][3]; //Pt balance at reco level gen unmatched jets only (= PU jets)
  TH1D *  h_PtRecoJetoverPtRecoZ_GenMatchedJet[Nptbins][Netabins][3];//Pt balance at reco level, gen matched jets only (=hard scattering)
  TH1D *  h_PtRecoJetoverPtRecoZ_PU[Nptbins][Netabins][3]; //Pt balance for PU (using a template)  

  //Histos for response = pt(reco jet)/pt(gen jet)
  TH1D *  h_PtRecoJetoverPtGenJet[Nptbins][Netabins][3]; // Here the Nptbins is a binning on reco Zpt or reco jet pt
  TH1D *  h_PtRecoJetoverPtGenJet_GentPtBinning[Nptbins][Netabins][3]; //  Here the Nptbins is a binning on gen jet pt 
  

  
  //Jet reco pt distribution for PU events. This is used for building h_PtRecoJetoverPtRecoZ_PU.
  TH1D * h_PtRecoJetPU[Netabins];
  //Z reco pt distribution for PU events. This is used for building h_PtRecoJetoverPtRecoZ_PU.
  TH1D * h_PtRecoZPU =new TH1D("h_PtRecoZPU","",200,0,100);  h_PtRecoZPU->Sumw2();
  //Next histo is just a cross check 
  TH1D * h_PtRecoZPU_ptlt30 =new TH1D("h_PtRecoZPU_ptlt30","",200,0,100);
  h_PtRecoZPU_ptlt30->Sumw2();
   
  //NPV(=#reco vertices) validation plots. To study whether a given selection introduces a bias in NPV.
  TH1D* h_npv = new TH1D("h_npv","",100,0,100);h_npv->Sumw2();
  TH1D* h_npv_max1centraljet = new TH1D("h_npv_max1centraljet","",100,0,100);h_npv_max1centraljet->Sumw2();
  TH1D* h_npv_max1jet = new TH1D("h_npv_max1jet","",100,0,100);h_npv_max1jet->Sumw2();
  TH1D* h_npvall = new TH1D("h_npvall","",100,0,100);h_npvall->Sumw2();
  TH1D* h_npv_max1jetpt10 = new TH1D("h_npv_max1jetpt10","",100,0,100);h_npv_max1jetpt10->Sumw2();
   

  //Initializing the whole histograms
  //Looping on pt 
  for(int i = 0; i <  Nptbins ; i++){
    double ptmin = ptbins[i]; 
    double ptmax = ptbins[i+1];
    TString ptmin_str = (TString) Form("%.f",ptmin);
    TString ptmax_str = (TString) Form("%.f",ptmax);
    ptmin_str.ReplaceAll(".","p");
    ptmax_str.ReplaceAll(".","p");
    //Making several folders based on whether the binning is performed in Zpt, recojet pt or both conditions (A<Zpt<B AND A<reco jet pt < B)
    TString dirname="zpt_"+ptmin_str+"_"+ptmax_str;
    outf->mkdir(dirname);
    dirname="jetpt_"+ptmin_str+"_"+ptmax_str;
    outf->mkdir(dirname);
    dirname="jetandzpt_"+ptmin_str+"_"+ptmax_str;
    outf->mkdir(dirname);
    dirname="genjetpt_"+ptmin_str+"_"+ptmax_str;
    outf->mkdir(dirname);
    //Looping on eta
    for(int j = 0; j <  Netabins ; j++){
      double etamin = etabins[j]; 
      double etamax = etabins[j+1];
      TString etamin_str = (TString) Form("%3.3f",etamin);
      TString etamax_str = (TString) Form("%3.3f",etamax);
      etamin_str.ReplaceAll(".","p");
      etamax_str.ReplaceAll(".","p");
      etamin_str.ReplaceAll("-","min");
      etamax_str.ReplaceAll("-","min");
      TString  suffix_etaptbin = "_pt"+ptmin_str + "to"+ptmax_str + "_eta"+ etamin_str+ "to"+etamax_str ;
      if(i==0){	 h_PtRecoJetPU[j] =new TH1D("h_PtRecoJetPU_eta"+ etamin_str+ "to"+etamax_str,"",200,0,100);	h_PtRecoJetPU[j]->Sumw2();  }
      
      //Looping on the type of pt binning: Zpt, recojet pt or both conditions (A<Zpt<B AND A<reco jet pt < B)
      
      for(int k = 0 ; k<3 ; k++){
	TString pttype ;
	suffix_etaptbin = "_pt"+ptmin_str + "to"+ptmax_str + "_eta"+ etamin_str+ "to"+etamax_str ;
	if(k==0) pttype = "_zptbinning";
	if(k==1) pttype = "_jetptbinning";
	if(k==2) pttype = "_jetandzptbinning";
	suffix_etaptbin +=pttype;
	
	h_PtGenJetoverPtGenZ[i][j][k] =new TH1D("h_PtGenJetoverPtGenZ"+suffix_etaptbin,"",1000,0,5);
	h_PtRecoJetoverPtRecoZ_testsample[i][j][k] =new TH1D("h_PtRecoJetoverPtRecoZ_testsample"+suffix_etaptbin,"",1000,0,5);
	h_PtRecoJetoverPtRecoZ_PU[i][j][k] =new TH1D("h_PtRecoJetoverPtRecoZ_PU"+suffix_etaptbin,"",1000,0,5);
	h_PtRecoJetoverPtRecoZ_GenMatchedJet[i][j][k] =new TH1D("h_PtRecoJetoverPtRecoZ_GenMatchedJet"+suffix_etaptbin,"",1000,0,5);
	h_PtRecoJetoverPtRecoZ_NoGenMatchedJet[i][j][k] =new TH1D("h_PtRecoJetoverPtRecoZ_NoGenMatchedJet"+suffix_etaptbin,"",1000,0,5);
	h_PtRecoZoverPtGenZ[i][j][k] =new TH1D("h_PtRecoZoverPtGenZ"+suffix_etaptbin,"",1000,0,5);
	h_PtRecoJetoverPtGenJet[i][j][k] =new TH1D("h_PtRecoJetoverPtGenJet"+suffix_etaptbin,"",2000,0,2);
	
	
	suffix_etaptbin = "_pt"+ptmin_str + "to"+ptmax_str + "_eta"+ etamin_str+ "to"+etamax_str ;
	h_PtRecoJetoverPtGenJet_GentPtBinning[i][j][k] =new TH1D("h_PtRecoJetoverPtGenJet_GentPtBinning"+suffix_etaptbin+"_"+pubin[k],"",1000,0,5);
	
	h_PtGenJetoverPtGenZ[i][j][k]->Sumw2();
	h_PtRecoJetoverPtRecoZ_testsample[i][j][k]->Sumw2();
	h_PtRecoJetoverPtRecoZ_PU[i][j][k]->Sumw2();
	h_PtRecoJetoverPtRecoZ_GenMatchedJet[i][j][k]->Sumw2();
	h_PtRecoJetoverPtRecoZ_NoGenMatchedJet[i][j][k]->Sumw2();
	h_PtRecoZoverPtGenZ[i][j][k]->Sumw2();
	h_PtRecoJetoverPtGenJet[i][j][k]->Sumw2();
	h_PtRecoJetoverPtGenJet_GentPtBinning[i][j][k]->Sumw2();
      }
    }
  }
  
  
  Long64_t nentries = thetree->GetEntries()*fractouse;
  //When useresiduals = false, only MC JECs are applied to data (and not the residuals). In that case, the attempt is to derive both residuals JES and JERSF. 
  //For now keep using useresiduals = true and focus on JER SF on fully calibrated jets
  Init(thetree,useresiduals);
  
  //Some counters
  int ctr_nomatchedjet(0), ctr_nomatchedl(0);
  int ctr_all(0);
  
  //First event loop to get jet pt spectrum from PU and Zpt spectrum when there is no additional jet
  cout << nentries <<endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry%100000==0) cout << jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    
    if (ientry < 0) break;
    thetree->GetEntry(jentry);
    ShiftSyst(unctytype);//Shift a given value to its up/down variation
    double totweight=1;
    //PU reweighting. No reweighting for UL right now (to be added, small effect anyways)
    if( samplename.Index("UL2017")>=0 &&samplename.Index("MC")>=0 ) totweight*=1;
    else if( samplename.Index("UL2018")>=0 &&samplename.Index("MC")>=0 ) totweight*=1;
    else if( samplename.Index("2018")>=0 && samplename.Index("MC")>=0 ) totweight*=PUReweighting2018vs102X(trueNVtx);
    else if( samplename.Index("2017")>=0 && samplename.Index("MC")>=0 ) totweight*=PUReweighting2017vs94X(trueNVtx);
    else if( samplename.Index("2016")>=0 && samplename.Index("MC")>=0 ) totweight*=PUReweighting2016vs94X(trueNVtx);

    TLorentzVector p4dilepton;
    TLorentzVector p4dilepton_gen;
    p4dilepton_gen.SetPtEtaPhiM(0.,0.,0.,0.);
    
    //Look for a dilepton pair
    bool isdileptonpair= FindDileptonPair( p4dilepton, p4dilepton_gen, samplename);
    if(!isdileptonpair ) continue;
    
    
    h_nvtx->Fill(_n_PV,totweight);
    int njetspt30_eta2p4(0),njetspt25(0), nmatchedjets(0);
    int njetspt20_eta2p4(0), njetspt10_eta2p4(0);
    int njetspt30_alleta(0), njetspt10_alleta(0);
    int ncentraljetslowpt(0);
    //Compute the nb of central jets with pt>30
    vector <int> idxcandjet ;
    bool passjetsel = JetStuff( idxcandjet , njetspt30_eta2p4,  njetspt25,  nmatchedjets,  njetspt20_eta2p4 , njetspt10_eta2p4 , njetspt30_alleta,  njetspt10_alleta, ncentraljetslowpt );
    if(p4dilepton.Pt()>30&& p4dilepton.Pt()<50 &&p4dilepton.Mag()>80&& p4dilepton.Mag()<100&& idxcandjet.size()>0){
      if(njetspt10_alleta>=1)h_npvall->Fill(_n_PV,totweight);
      if(njetspt10_alleta==1) h_npv_max1jetpt10->Fill(_n_PV,totweight);
    }
    
    if(!passjetsel) continue;

    //Filling Zpt distribution
    if(njetspt20_eta2p4 ==0&& njetspt30_alleta ==0&& p4dilepton.Pt()>ptZmin)  h_PtRecoZPU->Fill(p4dilepton.Pt(),totweight);
    if(njetspt20_eta2p4 ==0&& njetspt30_alleta <=1&& p4dilepton.Pt()>ptZmin)  h_PtRecoZPU_ptlt30->Fill(p4dilepton.Pt(),totweight);
    
    for(unsigned int i = 0; i <idxcandjet.size();i++){
      int idxj = idxcandjet[i];
      double dphi = fabs(acos(cos( (*_jetPhi)[idxj] -p4dilepton.Phi() )));
      
      int etabin(-1) , ptbin (-1);
      FindEtaPtbin( (*_jetPt)[idxj]  ,(*_jetEta)[idxj], ptbin, etabin);
      bool passTMLPUID[3]; 
      passPUID(idxj, passTMLPUID);
      int mybinpt(-1), mybineta(-1);
      FindEtaPtbinforPUID( (*_jetPt)[idxj],(*_jetEta)[idxj], mybinpt,mybineta );
      //Fill jet pt distribution for PU template
      //Select events with a low pt Z, and large dphi(jet, Z)
      if(dphi<2.0&&dphi>1.0 && p4dilepton.Pt()<5&&p4dilepton.Mag()>80&& p4dilepton.Mag()<100 && etabin>=0) {
	h_PtRecoJetPU[etabin]->Fill((*_jetPt)[idxj],totweight);   
      }
      if(mybinpt>=0&&mybineta>=0){
	
	if(dphi<2.0&&dphi>1.0 && p4dilepton.Pt()<5&&p4dilepton.Mag()>80&& p4dilepton.Mag()<100 ) {
	  
	  //PUID studies (here for PU jets)
	  h_PUID_PUjetsden->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(passTMLPUID[0]) h_PUID_T_PUjets->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(passTMLPUID[1]) h_PUID_M_PUjets->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(passTMLPUID[2]) h_PUID_L_PUjets->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(mybinpt>=0&&mybineta>=0)h_PUMVAID_PUjets[mybineta][mybinpt]->Fill( (*_jetPUMVAUpdate2017)[idxj],totweight); 
	  
	  for(int i = 0; i<nbinsfornvtx;i++){ 
	    if(_n_PV>nvtxbins[i]&&_n_PV<=nvtxbins[i+1]){
	      if(passTMLPUID[0]) h_PUID_T_PUjets_nvtxbinned[i]->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	      if(passTMLPUID[1]) h_PUID_M_PUjets_nvtxbinned[i]->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	      if(passTMLPUID[2]) h_PUID_L_PUjets_nvtxbinned[i]->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	      h_PUID_PUjetsden_nvtxbinned[i]->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	      h_nvtx_nvtxbinned[i]->Fill(_n_PV,totweight);
	    }
	  }
	}
	if(isMC && (*_jetPtGen)[idxj]==0 ){	  
	  h_PUMVAID_PUgenjets[mybineta][mybinpt]->Fill( (*_jetPUMVAUpdate2017)[idxj],totweight); 
	  h_PUID_PUgenjetsden->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(passTMLPUID[0]) h_PUID_T_PUgenjets->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(passTMLPUID[1]) h_PUID_M_PUgenjets->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(passTMLPUID[2]) h_PUID_L_PUgenjets->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	}
      }
    }
  }
  

  //Second loop to get the pt balance
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry%100000==0) cout << jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    
    if (ientry < 0) break;
    thetree->GetEntry(jentry);
    ShiftSyst(unctytype);
    double totweight=1;
    //PU reweighting. No reweighting for UL right now (to be added, small effect anyways)
    if( samplename.Index("UL2017")>=0 &&samplename.Index("MC")>=0 ) totweight*=1;
    else if( samplename.Index("UL2018")>=0 &&samplename.Index("MC")>=0 ) totweight*=1;
    else if(samplename.Index("2018")>=0  && samplename.Index("MC")>=0  ) totweight*=PUReweighting2018vs102X(trueNVtx);
    else if(samplename.Index("2017")>=0  && samplename.Index("MC")>=0 ) totweight*=PUReweighting2017vs94X(trueNVtx);
    else if(samplename.Index("2016")>=0  && samplename.Index("MC")>=0 ) totweight*=PUReweighting2016vs94X(trueNVtx);
    TLorentzVector p4dilepton;
    TLorentzVector p4dilepton_gen;
    p4dilepton_gen.SetPtEtaPhiM(0.,0.,0.,0.);
    
    bool isdileptonpair = FindDileptonPair( p4dilepton, p4dilepton_gen, samplename);
    if(!isdileptonpair ) continue;
    
    int njetspt30_eta2p4(0),njetspt25(0), nmatchedjets(0);
    int njetspt20_eta2p4(0), njetspt10_eta2p4(0);
    int njetspt30_alleta(0), njetspt10_alleta(0);
    int ncentraljetslowpt(0);
    //Compute the nb of central jets with pt>30
    
    vector <int> idxcandjet ;
    if(!JetStuff( idxcandjet , njetspt30_eta2p4,  njetspt25,  nmatchedjets,  njetspt20_eta2p4 , njetspt10_eta2p4 , njetspt30_alleta,  njetspt10_alleta , ncentraljetslowpt ))continue;
    
    //Building the PU template
    //Use a similar selection as below but require 0 extra jet instead of 1 
    //Then "add" a random PU jet using the pt distribution obtained above for PU jet.
    if(njetspt20_eta2p4 ==0&& njetspt30_alleta ==0 && ncentraljetslowpt ==0) {
      int thebinptz ;
      FindPtbin(p4dilepton.Pt(), thebinptz);
      for(int ibineta = 0; ibineta <Netabins ; ibineta++){
	//Generate  5 (this is an arbitrary number) entries according to the h_PtRecoJetPU distribution
	for(int kk =0; kk <5;kk++){
	  double ptbalpu =0;
	  double randompt= h_PtRecoJetPU[ibineta]->GetRandom() ; 
	  int thebinptjet; 
	  if ( h_PtRecoJetPU[ibineta]->Integral()>0.) 	   ptbalpu =   randompt /p4dilepton.Pt();
	  if(ptbalpu>0) h_PtRecoJetoverPtRecoZ_PU[thebinptz][ibineta][0]->Fill(ptbalpu,totweight);
	  FindPtbin(randompt, thebinptjet);
	  if(ptbalpu>0) h_PtRecoJetoverPtRecoZ_PU[thebinptjet][ibineta][1]->Fill(ptbalpu,totweight);
	}
      }
    }

    
    //Now looking at events with >=1 candidate jet
    if(idxcandjet.size()<=0 ) continue;
    for(unsigned int i = 0; i <idxcandjet.size();i++){
      int idxj = idxcandjet[i];
      
      double dphi = fabs(acos(cos(  (*_jetPhi)[idxj] -p4dilepton.Phi() )));
      double dphisigned = dphi /3.1415926536 ; 
      
      if(sin(  (*_jetPhi)[idxj] -p4dilepton.Phi())<0 ) dphisigned = (2.-dphisigned);
      
      //Next lines are related to PUID studies
      bool passTMLPUID[3]; 
      passPUID(idxj, passTMLPUID);
      int mybinpt(-1),mybineta(-1) ;
      FindEtaPtbinforPUID( (*_jetPt)[idxj],(*_jetEta)[idxj], mybinpt,mybineta );
      
      if(mybinpt>=0&&mybineta>=0) {
	int ik= (   p4dilepton.Pt() /  (*_jetPt)[idxj]  > 0.5 && p4dilepton.Pt() /  (*_jetPt)[idxj] <1.5  )  ? 0 : 1;
	if(   p4dilepton.Pt() /  (*_jetPt)[idxj] >1.5) ik=2;
	
	if(_eventNb%2==0||!isMC){
	  h_dphiZjets_ALL[mybineta][mybinpt][ik]->Fill(dphisigned,totweight);
	  if(passTMLPUID[0])h_dphiZjets_T[mybineta][mybinpt][ik]->Fill(dphisigned,totweight);
	  if(passTMLPUID[1])h_dphiZjets_M[mybineta][mybinpt][ik]->Fill(dphisigned,totweight);
	  if(passTMLPUID[2])h_dphiZjets_L[mybineta][mybinpt][ik]->Fill(dphisigned,totweight);
	}
	if(_eventNb%2==1&& isMC&& (*_jetPtGen)[idxj]>0){
	  h_dphiZjets_genmatched_ALL[mybineta][mybinpt][ik]->Fill(dphisigned,totweight);
	  if(passTMLPUID[0])h_dphiZjets_genmatched_T[mybineta][mybinpt][ik]->Fill(dphisigned,totweight);
	  if(passTMLPUID[1])h_dphiZjets_genmatched_M[mybineta][mybinpt][ik]->Fill(dphisigned,totweight);
	  if(passTMLPUID[2])h_dphiZjets_genmatched_L[mybineta][mybinpt][ik]->Fill(dphisigned,totweight);
	}
	if(_eventNb%2==1&& isMC&& (*_jetPtGen)[idxj]==0){
	  h_dphiZjets_genunmatched_ALL[mybineta][mybinpt][ik]->Fill(dphisigned,totweight);
	  if(passTMLPUID[0])h_dphiZjets_genunmatched_T[mybineta][mybinpt][ik]->Fill(dphisigned,totweight);
	  if(passTMLPUID[1])h_dphiZjets_genunmatched_M[mybineta][mybinpt][ik]->Fill(dphisigned,totweight);
	  if(passTMLPUID[2])h_dphiZjets_genunmatched_L[mybineta][mybinpt][ik]->Fill(dphisigned,totweight);
	}
	
      }
	
      if((*_jetPt)[idxj] / p4dilepton.Pt() <1.3&& (*_jetPt)[idxj] / p4dilepton.Pt()>0.7){
	if(dphi>2.9){//dphi too tight could be an issue in HF? 
	  h_PUID_Realjetsden->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(passTMLPUID[0]) h_PUID_T_Realjets->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(passTMLPUID[1]) h_PUID_M_Realjets->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(passTMLPUID[2]) h_PUID_L_Realjets->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(mybinpt>=0&&mybineta>=0){
	    h_PUMVAID_Realjets[mybineta][mybinpt]->Fill( (*_jetPUMVAUpdate2017)[idxj],totweight);
	    
	    if(isMC && (*_jetPtGen)[idxj]>0){
		h_PUID_Realgenjetsden->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
		if(passTMLPUID[0]) h_PUID_T_Realgenjets->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
		if(passTMLPUID[1]) h_PUID_M_Realgenjets->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
		if(passTMLPUID[2]) h_PUID_L_Realgenjets->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
		h_PUMVAID_Realgenjets[mybineta][mybinpt]->Fill( (*_jetPUMVAUpdate2017)[idxj],totweight);
		
	    }
	  }
	}
	else if(dphi<2.0&&dphi>1.0){
	  h_PUID_RealjetsPUTempden->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(passTMLPUID[0]) h_PUID_T_RealjetsPUTemp->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(passTMLPUID[1]) h_PUID_M_RealjetsPUTemp->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(passTMLPUID[2]) h_PUID_L_RealjetsPUTemp->Fill((*_jetPt)[idxj],(*_jetEta)[idxj],totweight);
	  if(mybinpt>=0&&mybineta>=0) h_PUMVAID_RealjetsPUTemp[mybineta][mybinpt]->Fill( (*_jetPUMVAUpdate2017)[idxj],totweight);
	}
      }

      //End of PUID stuff
      
      
      //Now JER SF stuff
      int thebinpt(-1),thebineta(-1);
      //Loop over the three binnings
      //k == 0: binning using Zpt, k == 1 binning using jet pt, k == 2, binning using A < Zpt < B AND A < jetpt < B
      for(int k=0;k<3;k++){
	thebinpt = -1;thebineta=-1;
	// Check the eta/pt bin index associated to a given pt/eta
	if(k==0) FindEtaPtbin(p4dilepton.Pt(),(*_jetEta)[idxj], thebinpt,thebineta);
	if(k>=1) FindEtaPtbin( (*_jetPt)[idxj]  ,(*_jetEta)[idxj], thebinpt,thebineta);
	if(k==2&& (p4dilepton.Pt()< ptbins[thebinpt] || p4dilepton.Pt()> ptbins[thebinpt+1] ) )continue;
	
	if(thebinpt<0|| thebineta<0) continue;
	
	//Z and jets should be back to back
	if(dphi<3.0) continue;
	
	h_npv->Fill(_n_PV,totweight);
	
	if(njetspt20_eta2p4>1) continue; //At most 1 jet with pt >20 GeV in tracker acceptance (|eta|<2.4) 
	if(k==0)h_npv_max1centraljet->Fill(_n_PV,totweight);
	if(njetspt30_alleta>1) continue; //At most 1 jet with pt >30 GeV in full |eta| range 
	if(k==0)h_npv_max1jet->Fill(_n_PV,totweight);

	//In MC, split events according to their event nb. Half used to build the distribution used in the fit, half for comparison with data.
	if(_eventNb%2==0||!isMC) h_PtRecoJetoverPtRecoZ_testsample[thebinpt][thebineta][k]->Fill((*_jetPt)[idxj] / p4dilepton.Pt() ,totweight); //For comparison with data
	if(_eventNb%2==1 && isMC){//For building the fit template
	  //Fill the various pt balance histos
	    if( (*_jetPtGen)[idxj] !=0&& p4dilepton_gen.Pt() !=0) h_PtRecoJetoverPtRecoZ_GenMatchedJet[thebinpt][thebineta][k]->Fill((*_jetPt)[idxj] / p4dilepton.Pt(),totweight );
	    if( (*_jetPtGen)[idxj] ==0&& p4dilepton_gen.Pt() !=0) h_PtRecoJetoverPtRecoZ_NoGenMatchedJet[thebinpt][thebineta][k]->Fill((*_jetPt)[idxj] / p4dilepton.Pt() ,totweight);
	    if( (*_jetPtGen)[idxj] !=0&& p4dilepton_gen.Pt() !=0) h_PtGenJetoverPtGenZ[thebinpt][thebineta][k]->Fill((*_jetPtGen)[idxj] / p4dilepton_gen.Pt() ,totweight);
	    if(  p4dilepton_gen.Pt() !=0 )h_PtRecoZoverPtGenZ[thebinpt][thebineta][k]->Fill(p4dilepton.Pt() / p4dilepton_gen.Pt() ,totweight);
	    if( (*_jetPtGen)[idxj] !=0 )  h_PtRecoJetoverPtGenJet[thebinpt][thebineta][k]->Fill( (*_jetPt)[idxj]/ (*_jetPtGen)[idxj] ,totweight);
	}

      }
      if(isMC&& (*_jetPtGen)[idxj] !=0 ){//gen matched jet
	FindEtaPtbin( (*_jetPtGen)[idxj]  ,(*_jetEta)[idxj], thebinpt,thebineta);
	int ipu  ;
	if(_n_PV<20) ipu =0;
	else if(_n_PV<=40) ipu =1;
	else ipu =2;
	  if(_eventNb%2==1) h_PtRecoJetoverPtGenJet_GentPtBinning[thebinpt][thebineta][ipu]->Fill( (*_jetPt)[idxj]/ (*_jetPtGen)[idxj] ,totweight);
      }
      
      //Some counters, not important
      ctr_all++;
      if(isMC&& (*_jetPtGen)[idxj] >=0 ) ctr_nomatchedjet ++;
      if(isMC&& p4dilepton_gen.Pt() !=0 ) ctr_nomatchedl ++;
      
    }
      
  }
   
  //Now just writing the histos into a file
   outf->cd();
   
   for(int j = 0; j <  Netabins ; j++)  h_PtRecoJetPU[j]->Write();
   h_PtRecoZPU->Write();
   h_PtRecoZPU_ptlt30->Write();
   
   h_npv->Write();
   h_npv_max1centraljet->Write();
   h_npv_max1jet->Write();
   h_npvall->Write();
   h_npv_max1jetpt10->Write();
   
   for(int i = 0; i <  Nptbins ; i++){
     double ptmin = ptbins[i];
     double ptmax = ptbins[i+1];
     TString ptmin_str = (TString) Form("%.f",ptmin);
     TString ptmax_str = (TString) Form("%.f",ptmax);
     ptmin_str.ReplaceAll(".","p");
     ptmax_str.ReplaceAll(".","p");
     for(int k = 0 ; k<4 ; k++){
       TString dirname;
       if(k==0) dirname ="zpt_"+ptmin_str+"_"+ptmax_str;
       if(k==1) dirname ="jetpt_"+ptmin_str+"_"+ptmax_str;
       if(k==2) dirname ="jetandzpt_"+ptmin_str+"_"+ptmax_str;
       if(k==3)  dirname ="genjetpt_"+ptmin_str+"_"+ptmax_str;
       outf->cd(dirname);
       for(int j = 0; j <  Netabins ; j++){
	 if(k<3){
	   
	   h_PtRecoJetoverPtRecoZ_testsample[i][j][k]->Write();
	   h_PtRecoJetoverPtRecoZ_PU[i][j][k]->Write();
	   if(isMC) {
	     h_PtRecoJetoverPtRecoZ_GenMatchedJet[i][j][k]->Write();
	     h_PtRecoJetoverPtRecoZ_NoGenMatchedJet[i][j][k]->Write();
	     h_PtRecoZoverPtGenZ[i][j][k]->Write();
	     h_PtRecoJetoverPtGenJet[i][j][k]->Write();
	     h_PtGenJetoverPtGenZ[i][j][k]->Write(); 
	   }
	 }	   
	 if(k==3&&isMC) {
	   for(int l=0; l<3 ;l++) h_PtRecoJetoverPtGenJet_GentPtBinning[i][j][l]->Write();
	 }
       }
     }
   }
   
   outf->cd();
   //All what follows is related to PUID stuff
   for(int i = 0; i < nbinseta_forPUID ; i++){
     for(int j = 0; j < nbinspt_forPUID ; j++){
       if(isMC)h_PUMVAID_PUjets[i][j]->Write();
       if(isMC)h_PUMVAID_PUgenjets[i][j]->Write();
       h_PUMVAID_Realjets[i][j]->Write();
       h_PUMVAID_RealjetsPUTemp[i][j]->Write();
       h_PUMVAID_Realgenjets[i][j]->Write();
       
       for(int k = 0; k<3;k++){
	 
	 if(k==2){
	   h_dphiZjets_ALL[i][j][k]->Add( h_dphiZjets_ALL[i][j][0]);  h_dphiZjets_ALL[i][j][k]->Add( h_dphiZjets_ALL[i][j][1]); 
	   h_dphiZjets_T[i][j][k]->Add( h_dphiZjets_T[i][j][0]);  h_dphiZjets_T[i][j][k]->Add( h_dphiZjets_T[i][j][1]); 
	   h_dphiZjets_M[i][j][k]->Add( h_dphiZjets_M[i][j][0]);  h_dphiZjets_M[i][j][k]->Add( h_dphiZjets_M[i][j][1]); 
	   h_dphiZjets_L[i][j][k]->Add( h_dphiZjets_L[i][j][0]);  h_dphiZjets_L[i][j][k]->Add( h_dphiZjets_L[i][j][1]); 
	 }
	 h_dphiZjets_ALL[i][j][k]->Write();
	 h_dphiZjets_T[i][j][k]->Write();
	 h_dphiZjets_M[i][j][k]->Write();
	 h_dphiZjets_L[i][j][k]->Write();
	 
	 if(isMC){
	   if(k==2){
	     h_dphiZjets_genmatched_ALL[i][j][k]->Add( h_dphiZjets_genmatched_ALL[i][j][0]);  h_dphiZjets_genmatched_ALL[i][j][k]->Add( h_dphiZjets_genmatched_ALL[i][j][1]); 
	     h_dphiZjets_genmatched_T[i][j][k]->Add( h_dphiZjets_genmatched_T[i][j][0]);  h_dphiZjets_genmatched_T[i][j][k]->Add( h_dphiZjets_genmatched_T[i][j][1]); 
	     h_dphiZjets_genmatched_M[i][j][k]->Add( h_dphiZjets_genmatched_M[i][j][0]);  h_dphiZjets_genmatched_M[i][j][k]->Add( h_dphiZjets_genmatched_M[i][j][1]); 
	     h_dphiZjets_genmatched_L[i][j][k]->Add( h_dphiZjets_genmatched_L[i][j][0]);  h_dphiZjets_genmatched_L[i][j][k]->Add( h_dphiZjets_genmatched_L[i][j][1]); 
	     
	     h_dphiZjets_genunmatched_ALL[i][j][k]->Add( h_dphiZjets_genunmatched_ALL[i][j][0]);  h_dphiZjets_genunmatched_ALL[i][j][k]->Add( h_dphiZjets_genunmatched_ALL[i][j][1]); 
	     h_dphiZjets_genunmatched_T[i][j][k]->Add( h_dphiZjets_genunmatched_T[i][j][0]);  h_dphiZjets_genunmatched_T[i][j][k]->Add( h_dphiZjets_genunmatched_T[i][j][1]); 
	     h_dphiZjets_genunmatched_M[i][j][k]->Add( h_dphiZjets_genunmatched_M[i][j][0]);  h_dphiZjets_genunmatched_M[i][j][k]->Add( h_dphiZjets_genunmatched_M[i][j][1]); 
	     h_dphiZjets_genunmatched_L[i][j][k]->Add( h_dphiZjets_genunmatched_L[i][j][0]);  h_dphiZjets_genunmatched_L[i][j][k]->Add( h_dphiZjets_genunmatched_L[i][j][1]); 
	     
	   }
	   
	   h_dphiZjets_genmatched_ALL[i][j][k]->Write();
	   h_dphiZjets_genmatched_T[i][j][k]->Write();
	   h_dphiZjets_genmatched_M[i][j][k]->Write();
	   h_dphiZjets_genmatched_L[i][j][k]->Write();
	   
	   h_dphiZjets_genunmatched_ALL[i][j][k]->Write();
	   h_dphiZjets_genunmatched_T[i][j][k]->Write();
	   h_dphiZjets_genunmatched_M[i][j][k]->Write();
	   h_dphiZjets_genunmatched_L[i][j][k]->Write();
	 }
       }
     }
   }
   
   
   h_PUID_T_PUjets->Write();
   h_PUID_M_PUjets->Write();
   h_PUID_L_PUjets->Write();
   h_PUID_PUjetsden->Write();
   h_PUID_T_Realjets->Write();
   h_PUID_M_Realjets->Write();
   h_PUID_L_Realjets->Write();
   h_PUID_Realjetsden->Write();
   
   h_PUID_T_RealjetsPUTemp->Write();
   h_PUID_M_RealjetsPUTemp->Write();
   h_PUID_L_RealjetsPUTemp->Write();
   h_PUID_RealjetsPUTempden->Write();
   
   if(isMC){
     h_PUID_T_PUgenjets->Write();
     h_PUID_M_PUgenjets->Write();
     h_PUID_L_PUgenjets->Write();
     h_PUID_PUgenjetsden->Write();
     h_PUID_T_Realgenjets->Write();
     h_PUID_M_Realgenjets->Write();
     h_PUID_L_Realgenjets->Write();
     h_PUID_Realgenjetsden->Write();
   }
   
   
   h_nvtx->Write();
   for(int i = 0; i<nbinsfornvtx;i++){ 
     h_PUID_T_PUjets_nvtxbinned[i]->Write();
     h_PUID_M_PUjets_nvtxbinned[i]->Write();
     h_PUID_L_PUjets_nvtxbinned[i]->Write();
     h_PUID_PUjetsden_nvtxbinned[i]->Write();
     h_nvtx_nvtxbinned[i]->Write();
   }
   
   
   outf->Close();
   myf->Close();
   
   
   cout<< "ALL,  matched jet/lepton " <<ctr_all<<"," << ctr_nomatchedjet <<", "<< ctr_nomatchedl <<endl;
   
}


TLorentzVector templateforJERC::GenLP4(int irecol){
  double drMAX = 0.3;
  TLorentzVector lgenptp4 ;
  lgenptp4.SetPtEtaPhiM(0.,0.,0.,0.);
  int idx = -1;
  if(_lgenPt==0)   return lgenptp4 ;
  for(unsigned int igen =0; igen<(*_lgenPt).size(); igen++){
    double deta = (*_lEta)[irecol] - (*_lgenEta)[igen];
    double dphi = fabs(acos(cos ((*_lPhi)[irecol] - (*_lgenPhi)[igen])));
    double dr =sqrt(deta*deta+dphi*dphi);
    double ptbal  = fabs( (*_lPt)[irecol] - (*_lgenPt)[igen]  )/ ( (*_lPt)[irecol] + (*_lgenPt)[igen]) ;
    if(dr< drMAX &&ptbal<0.3 ) {
      drMAX=dr;
      idx = igen;
      double ml = ( fabs(  (*_lgenpdgId)[igen] ) == 13) ? 0.1:0.; 
      lgenptp4.SetPtEtaPhiM( (*_lgenPt)[igen], (*_lgenEta)[igen] , (*_lgenPhi)[igen],ml);
    }
  }
  return lgenptp4 ; 
}

void templateforJERC::FindEtaPtbin(double jetpt, double jeteta, int & thebinpt , int & thebineta){

  thebinpt = -1;
  for(int i = 0; i<= Nptbins ;i++){
    if(jetpt<ptbins[i] ){ thebinpt = i-1; break ;}
  }
  thebineta = -1;
  for(int i = 0; i<= Netabins ;i++){
    if(jeteta<etabins[i] ){ thebineta = i-1; break ;}
  }
      
}

void templateforJERC::FindPtbin(double jetpt, int & thebinpt ){

  thebinpt = -1;
  for(int i = 0; i<= Nptbins ;i++){
    if(jetpt<ptbins[i] ){ thebinpt = i-1; break ;}
  }
}



void templateforJERC::FindEtaPtbinforPUID(double jetpt, double jeteta, int & thebinpt , int & thebineta){
  thebinpt = -1;
  for(int i = 0; i<= nbinspt_forPUID ;i++){

    if(jetpt<binpt_forPUID[i] ){ thebinpt = i-1; break ;}
  }
  thebineta = -1;
  for(int i = 0; i<= nbinseta_forPUID ;i++){
    if(jeteta<bineta_forPUID[i] ){ thebineta = i-1; break ;}
  }
      
}


bool templateforJERC::FindDileptonPair(TLorentzVector & p4dil_reco, TLorentzVector & p4dil_gen, TString samplename, bool applyptcut){
  
  //  if(trueNVtx<55) return false;
  //  if(_n_PV>=40||_n_PV<20) return false;
  p4dil_gen.SetPtEtaPhiM(0.,0.,0.,0.);
  int ndileptonpair(0), fldileptonpair(-1);  
  for(unsigned int i = 0; i< (*_lPt).size() ; i++){
    if( abs((*_lpdgId )[i])  != 11 &&  abs((*_lpdgId )[i])  != 13)continue;
    if( !(*_lPassTightID )[i] ||  (*_lPt)[i]<15  ) continue;
    for(unsigned int j = 0; j< i ; j++){
      if( !(*_lPassTightID )[j]||  (*_lPt)[j]<20  ) continue;
      if( abs((*_lpdgId)[j]) ==11&&  (*_lPt)[j]<25  ) continue;
      
      if( (*_lpdgId)[i] != -(*_lpdgId)[j])continue;
      
      TLorentzVector l1,l2;
      double mlepton = ( abs((*_lpdgId)[i]) ==13) ? 0.1: 0.;
      l1.SetPtEtaPhiM( (*_lPt)[i], (*_lEta)[i], (*_lPhi)[i],mlepton);
      l2.SetPtEtaPhiM( (*_lPt)[j], (*_lEta)[j], (*_lPhi)[j],mlepton);
      ndileptonpair++;
      p4dil_reco = l1+l2;
      fldileptonpair = (abs((*_lpdgId)[i]) -11 )/2   ; //0 for ee, 1 for mumu
      TLorentzVector l1gen = GenLP4(i);
      TLorentzVector l2gen = GenLP4(j);
      p4dil_gen = l1gen+l2gen;
    }
  }

  if(ndileptonpair!=1) return false;
  if(!PassTrigger(fldileptonpair,samplename) ) return false;
  //  if(fldileptonpair ==1 && !HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8) return false;
  // if( fldileptonpair ==0 && ! HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL) return false;
  
  if(p4dil_reco.Mag()<70|| p4dil_reco.Mag()>110) return false;
  //  if(p4dil_reco.Mag()<85|| p4dil_reco.Mag()>95) return false;
  if(p4dil_reco.Pt()< ptZmin &&applyptcut) return false;
  
  //  if(_rho >10) return false;
  return true;
}


bool templateforJERC::JetStuff(vector <int>& idxjets, int & njetspt30_eta2p4, int & njetspt25, int & nmatchedjets, int & njetspt20_eta2p4 , int & njetspt10_eta2p4 , int &njetspt30_alleta, int & njetspt10_alleta , int & ncentraljetslowpt){
  
  ncentraljetslowpt = 0;
  int njetspt20_alleta = 0;
  for (unsigned int i = 0; i <(*_jetPt).size(); i++){
    if(seconditeration){
      int thebin =  hjecs->FindBin((*_jetPt) [i],fabs( (*_jetEta)[i]) );
      if( (*_jetPt)[i] <30 ) thebin =  hjecs->FindBin(31.,fabs( (*_jetEta)[i]) ); 
      double thecorr = hjecs->GetBinContent(thebin);

      //if((*_jetPt) [i]>25 )cout << "Orig jetpt " <<     (*_jetPt) [i]<< ", " <<(*_jetEta)[i]<<endl;
      (*_jetPt) [i] =   (*_jetPt) [i] / (1+thecorr);
      //      if((*_jetPt) [i]>25 )cout << "Recorr jetpt " <<     (*_jetPt) [i]<<endl;
    }

    if(!IsGoodJet(i) )continue;
    //    if((*_jetPt)[i]>20 &&fabs((*_jetEta)[i])<2.4 ) cout << "Jet pt/: " << (*_jetPt)[i]<<", " << (*_jetPtGen)[i]<<", "<< (*_jetPUMVAUpdate2017)[i] <<endl;
    if( (*_jetPt)[i]>20 ) njetspt20_alleta ++;
    if(fabs((*_jetEta)[i])<2.5&&(*_jetPt)[i]>12 &&(*_jetPt)[i]<ptbins[0]&& (*_jetPUMVAUpdate2017)[i] >0) ncentraljetslowpt++;
    if( (*_jetPt)[i]>30 ) njetspt30_alleta++;
    if(fabs((*_jetEta)[i])<2.4&&(*_jetPt)[i]>30 ) njetspt30_eta2p4++;
    if(fabs((*_jetEta)[i])<2.4&&(*_jetPt)[i]>20 ) njetspt20_eta2p4++;
    if(fabs((*_jetEta)[i])<2.4&&(*_jetPt)[i]>10 ) njetspt10_eta2p4++;
    if((*_jetPt)[i]>10 ) njetspt10_alleta++;
  }
  
  
  //This loop ensures there are no additional jet with pt >30 (all eta) ; pt>20 (|eta|<2.4), 12<pt<20 (|eta|<2.4), pass PUID)
  for (unsigned int i = 0; i <(*_jetPt).size(); i++){
    if(!IsGoodJet(i) )continue;
    int thenjetspt30_alleta = njetspt30_alleta;
    int thenjetspt20_eta2p4 = njetspt20_eta2p4;
    int thencentraljetslowpt  = ncentraljetslowpt;
    if((*_jetPt)[i]>30) thenjetspt30_alleta --; 
    if(fabs((*_jetEta)[i])<2.4&&(*_jetPt)[i]>20  )thenjetspt20_eta2p4--;
    if(fabs((*_jetEta)[i])<2.5&&(*_jetPt)[i]>12 &&(*_jetPt)[i]<ptbins[0]&& (*_jetPUMVAUpdate2017)[i] >0)thencentraljetslowpt--; 
    if( thenjetspt20_eta2p4 ==0&& thenjetspt30_alleta==0 &&thencentraljetslowpt ==0)  idxjets.push_back(i);
    
  }

  if(njetspt20_eta2p4>=2) return false;
  if(njetspt30_alleta>=2) return false;


  if(ncentraljetslowpt>=2) return false;  
  //  the next line is wrong !!


  //Checks only !!! 
  //  if( njetspt20_alleta>=2) return false;

  //  if(  idxjets.size() >=2 && njetspt10_alleta >=2) return false;
  return true;
}


bool templateforJERC::IsGoodJet(int ij){
  if(!(*_jetPassID)[ij] )return false;
  if((*_jetPt)[ij]< ThePtMin)return false;
  if(fabs((*_jetEta)[ij])<2.5&&(*_jetPt)[ij]<50 &&(*_jetPUMVAUpdate2017)[ij]<0. &&applypuid) return false;
  bool lmatch(false);
  for(unsigned int j = 0; j< (*_lPt).size() ; j++){
    if( !(*_lPassLooseID )[j]|| (*_lPt)[j]<5 ) continue;
    double deta = (*_jetEta)[ij] - (*_lEta)[j];
    double dphi = fabs(acos(cos ((*_jetPhi)[ij] - (*_lPhi)[j])));
    double dr =sqrt(deta*deta+dphi*dphi);
    if(dr<0.3) {lmatch =true; break;}
  }
  if(lmatch) return false;
  return true;
}

bool templateforJERC::PassTrigger(int fldileptonpair, TString samplename) {

  bool iselechannel = samplename.Index("ZtoEE")>=0 || samplename.Index("DoubleEG")>=0;
  bool ismuchannel = samplename.Index("ZtoMuMu")>=0 || samplename.Index("DoubleMuon")>=0;
  if(fldileptonpair ==0 && !iselechannel) return false;
  if(fldileptonpair ==1 && !ismuchannel) return false;
  if(fldileptonpair ==1 &&(samplename.Index("2016")>=0) ) return HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL ||HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ||HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ; 
  
  else if(fldileptonpair ==1 &&samplename.Index("2017B")>=0) return true;//I have an issue with trigger in this era.
  else if(fldileptonpair ==1 &&samplename.Index("2017")>=0) return HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
  else if(fldileptonpair ==1 &&samplename.Index("2018")>=0) return HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
  else if(fldileptonpair ==0 &&samplename.Index("2016")>=0) return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;    
  else if(fldileptonpair ==0 &&samplename.Index("2017")>=0) return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;    
  else if(fldileptonpair ==0 &&samplename.Index("2018")>=0) return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;    

  return false;
}




void templateforJERC::passPUID(int j, bool (&idresult)[3]){//TLM
  if((*_jetPt)[j]>50){
    idresult[0] = true;
    idresult[1] = true;
    idresult[2] = true;
    return;
  }
  else  if ((*_jetPt)[j]>30){
    idresult[0] = false;
    idresult[1] = false;
    idresult[2] = false;
    double abseta = fabs((*_jetEta)[j]); 
    if(abseta<2.5){
      if((*_jetPUMVAUpdate2017)[j]> 0.86)idresult[0] = true; 
      if((*_jetPUMVAUpdate2017)[j]> 0.61)idresult[1] = true; 
      if((*_jetPUMVAUpdate2017)[j]>-0.89)idresult[2] = true; 
    }
    else if(abseta<2.75){
      if((*_jetPUMVAUpdate2017)[j]> -0.1)idresult[0] = true; 
      if((*_jetPUMVAUpdate2017)[j]> -0.35)idresult[1] = true; 
      if((*_jetPUMVAUpdate2017)[j]> -0.52)idresult[2] = true; 
    }
    else if(abseta<3.0){
      if((*_jetPUMVAUpdate2017)[j]> -0.05)idresult[0] = true; 
      if((*_jetPUMVAUpdate2017)[j]> -0.23)idresult[1] = true; 
      if((*_jetPUMVAUpdate2017)[j]> -0.38)idresult[2] = true; 
    }
    else if(abseta<5.0){
      if((*_jetPUMVAUpdate2017)[j]> -0.01)idresult[0] = true; 
      if((*_jetPUMVAUpdate2017)[j]> -0.17)idresult[1] = true; 
      if((*_jetPUMVAUpdate2017)[j]> -0.30)idresult[2] = true; 
    }
  }
  else  if ((*_jetPt)[j]<=30){
    idresult[0] = false;
    idresult[1] = false;
    idresult[2] = false;
    double abseta = fabs((*_jetEta)[j]); 
    if(abseta<2.5){
      if((*_jetPUMVAUpdate2017)[j]> 0.69)idresult[0] = true; 
      if((*_jetPUMVAUpdate2017)[j]> 0.18)idresult[1] = true; 
      if((*_jetPUMVAUpdate2017)[j]>-0.97)idresult[2] = true; 
    }
    else if(abseta<2.75){
      if((*_jetPUMVAUpdate2017)[j]> -0.35)idresult[0] = true; 
      if((*_jetPUMVAUpdate2017)[j]> -0.55)idresult[1] = true; 
      if((*_jetPUMVAUpdate2017)[j]> -0.68)idresult[2] = true; 
    }
    else if(abseta<3.0){
      if((*_jetPUMVAUpdate2017)[j]> -0.26)idresult[0] = true; 
      if((*_jetPUMVAUpdate2017)[j]> -0.42)idresult[1] = true; 
      if((*_jetPUMVAUpdate2017)[j]> -0.53)idresult[2] = true; 
    }
    else if(abseta<5.0){
      if((*_jetPUMVAUpdate2017)[j]> -0.21)idresult[0] = true; 
      if((*_jetPUMVAUpdate2017)[j]> -0.36)idresult[1] = true; 
      if((*_jetPUMVAUpdate2017)[j]> -0.47)idresult[2] = true; 
    }
  }
  return;
}

void templateforJERC::ShiftSyst(TString unctytype){

  if(unctytype=="central" || unctytype=="") return;
  

  for (unsigned int i = 0; i <(*_jetPt).size(); i++){
    //cout <<"before " << (*_jetPt)[i]<<endl; 
    if(unctytype=="up")(*_jetPt)[i] = (*_jetPt)[i]*(1+ (*_jetJECuncty)[i]); 
    if(unctytype=="down")(*_jetPt)[i] = (*_jetPt)[i]*(1- (*_jetJECuncty)[i]); 
    if((*_jetPt)[i]<0) (*_jetPt)[i] = 1;
    //    cout << "after " <<(*_jetPt)[i]<<endl; 
  }
  return ;



}
