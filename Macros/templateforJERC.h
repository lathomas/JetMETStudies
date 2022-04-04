//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 13 16:15:10 2021 by ROOT version 6.12/07
// from TTree tree/tree
// found on file: /user/lathomas/L1Prefiring/CMSSW_10_6_25/src/JetMETStudies/JMEAnalyzer/python/smalltest.root
//////////////////////////////////////////////////////////

#ifndef templateforJERC_h
#define templateforJERC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "TLorentzVector.h"

class templateforJERC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *_jetEta;
   vector<float>   *_jetPhi;
   vector<float>   *_jetPt;
   vector<float>   *_jetRawPt;
   vector<float>   *_jetArea;
   vector<bool>    *_jetPassID;
   vector<bool>    *_jetLeptonPhotonCleaned;
   vector<float>   *_jetPtNoL2L3Res;
   vector<float>   *_jetPUMVAUpdate2017;
   vector<float>   *_jetPUMVAUpdate2018;
   vector<float>   *_puppijetEta;
   vector<float>   *_puppijetPhi;
   vector<float>   *_puppijetPt;
   vector<float>   *_puppijetRawPt;
   vector<bool>    *_puppijetLeptonPhotonCleaned;
   vector<bool>    *_puppijetPassID;
   vector<float>   *_puppijetPtNoL2L3Res;
   vector<float>   *_jetPtGen;
   vector<float>   *_jetEtaGen;
   vector<float>   *_jetPhiGen;
   vector<float>   *_jetPtGenWithNu;
   vector<float>   *_puppijetPtGen;
   vector<float>   *_puppijetPtGenWithNu;
   Int_t           trueNVtx;
   Float_t         _genmet;
   Float_t         _genmet_phi;
   vector<float>   *_jetJECuncty;
   vector<float>   *_puppijetJECuncty;
   Float_t         _met;
   Float_t         _met_phi;
   Float_t         _puppimet;
   Float_t         _puppimet_phi;
   Float_t         _rawmet;
   Float_t         _rawmet_phi;
   Float_t         _puppirawmet;
   Float_t         _puppirawmet_phi;
   Float_t         _rawchsmet;
   Float_t         _rawchsmet_phi;
   Float_t         _chsmet;
   Float_t         _chsmet_phi;
   Int_t           _n_PV;
   Float_t         _rho;
   ULong64_t       _eventNb;
   ULong64_t       _runNb;
   ULong64_t       _lumiBlock;
   Float_t         _ptll;
   Float_t         _phill;
   vector<float>   *_lEta;
   vector<float>   *_lPhi;
   vector<float>   *_lPt;
   vector<float>   *_lPtcorr;
   vector<int>     *_lpdgId;
   vector<bool>    *_lPassTightID;
   vector<bool>    *_lPassLooseID;
   Int_t           _nEles;
   Int_t           _nMus;
   Float_t         _ptll_gen;
   Float_t         _phill_gen;
   vector<float>   *_lgenPt;
   vector<float>   *_lgenEta;
   vector<float>   *_lgenPhi;
   vector<int>     *_lgenpdgId;
   Bool_t          HLT_Ele35_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf;
   Bool_t          HLT_Ele27_WPTight_Gsf;
   Bool_t          HLT_IsoMu27;
   Bool_t          HLT_IsoMu24;
   Bool_t          HLT_IsoTkMu24;
   Bool_t          HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Double_t        _weight;

   // List of branches
   TBranch        *b__jetEta;   //!
   TBranch        *b__jetPhi;   //!
   TBranch        *b__jetPt;   //!
   TBranch        *b__jetRawPt;   //!
   TBranch        *b__jetArea;   //!
   TBranch        *b__jetPassID;   //!
   TBranch        *b__jetLeptonPhotonCleaned;   //!
   TBranch        *b__jetPtNoL2L3Res;   //!
   TBranch        *b__jetPUMVAUpdate2017;   //!
   TBranch        *b__jetPUMVAUpdate2018;   //!
   TBranch        *b__puppijetEta;   //!
   TBranch        *b__puppijetPhi;   //!
   TBranch        *b__puppijetPt;   //!
   TBranch        *b__puppijetRawPt;   //!
   TBranch        *b__puppijetLeptonPhotonCleaned;   //!
   TBranch        *b__puppijetPassID;   //!
   TBranch        *b__puppijetPtNoL2L3Res;   //!
   TBranch        *b__jetPtGen;   //!
   TBranch        *b__jetEtaGen;   //!
   TBranch        *b__jetPhiGen;   //!
   TBranch        *b__jetPtGenWithNu;   //!
   TBranch        *b__puppijetPtGen;   //!
   TBranch        *b__puppijetPtGenWithNu;   //!
   TBranch        *b_trueNVtx;   //!
   TBranch        *b__genmet;   //!
   TBranch        *b__genmet_phi;   //!
   TBranch        *b__jetJECuncty;   //!
   TBranch        *b__puppijetJECuncty;   //!
   TBranch        *b__met;   //!
   TBranch        *b__met_phi;   //!
   TBranch        *b__puppimet;   //!
   TBranch        *b__puppimet_phi;   //!
   TBranch        *b__rawmet;   //!
   TBranch        *b__rawmet_phi;   //!
   TBranch        *b__puppirawmet;   //!
   TBranch        *b__puppirawmet_phi;   //!
   TBranch        *b__rawchsmet;   //!
   TBranch        *b__rawchsmet_phi;   //!
   TBranch        *b__chsmet;   //!
   TBranch        *b__chsmet_phi;   //!
   TBranch        *b__n_PV;   //!
   TBranch        *b__rho;   //!
   TBranch        *b__eventNb;   //!
   TBranch        *b__runNb;   //!
   TBranch        *b__lumiBlock;   //!
   TBranch        *b__ptll;   //!
   TBranch        *b__phill;   //!
   TBranch        *b__lEta;   //!
   TBranch        *b__lPhi;   //!
   TBranch        *b__lPt;   //!
   TBranch        *b__lPtcorr;   //!
   TBranch        *b__lpdgId;   //!
   TBranch        *b__lPassTightID;   //!
   TBranch        *b__lPassLooseID;   //!
   TBranch        *b__nEles;   //!
   TBranch        *b__nMus;   //!
   TBranch        *b__ptll_gen;   //!
   TBranch        *b__phill_gen;   //!
   TBranch        *b__lgenPt;   //!
   TBranch        *b__lgenEta;   //!
   TBranch        *b__lgenPhi;   //!
   TBranch        *b__lgenpdgId;   //!
   TBranch        *b_HLT_Ele35_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele27_WPTight_Gsf;   //!
   TBranch        *b_HLT_IsoMu27;   //!
   TBranch        *b_HLT_IsoMu24;   //!
   TBranch        *b_HLT_IsoTkMu24;   //!
   TBranch        *b_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b__weight;
   templateforJERC(TTree *tree=0);
   virtual ~templateforJERC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree,bool useresiduals=true);
   virtual void     Loop(TString samplename, bool useresiduals, TString unctytype="central");
   virtual void passPUID(int j, bool (&idresult)[3]);
   virtual void FindEtaPtbin(double jetpt, double jeteta, int & thebinpt , int & thebineta);
   virtual void FindEtaPtbinforPUID(double jetpt, double jeteta, int & thebinpt , int & thebineta);
   virtual void FindPtbin(double jetpt, int & thebinpt );
   virtual bool PassTrigger(int fldileptonpair, TString samplename) ;
   virtual bool IsGoodJet(int ij);
   virtual bool JetStuff(vector <int>&  idxjets, int & njetspt30_eta2p4, int & njetspt25, int & nmatchedjets, int & njetspt20_eta2p4 , int & njetspt10_eta2p4 , int &njetspt30_alleta, int & njetspt10_alleta , int & ncentraljetslowpt); 
   virtual bool FindDileptonPair(TLorentzVector & p4dil_reco, TLorentzVector & p4dil_gen, TString samplename, bool applyptcut=true);
   
   virtual TLorentzVector GenLP4(int irecol);
   virtual void ShiftSyst(TString unctytype);
   virtual void GetRecoilProjections(double TheMET, double TheMETphi, double zpt, double zphi, double &upar, double &uperp, TString hname);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};

#endif

#ifdef templateforJERC_cxx
templateforJERC::templateforJERC(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   Init(tree);
}

templateforJERC::~templateforJERC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t templateforJERC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t templateforJERC::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void templateforJERC::Init(TTree *tree, bool useresiduals)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   _jetEta = 0;
   _jetPhi = 0;
   _jetPt = 0;
   _jetRawPt = 0;
   _jetArea = 0;
   _jetPassID = 0;
   _jetLeptonPhotonCleaned = 0;
   _jetPtNoL2L3Res = 0;
   _jetPUMVAUpdate2017 = 0;
   _jetPUMVAUpdate2018 = 0;
   _puppijetEta = 0;
   _puppijetPhi = 0;
   _puppijetPt = 0;
   _puppijetRawPt = 0;
   _puppijetLeptonPhotonCleaned = 0;
   _puppijetPassID = 0;
   _puppijetPtNoL2L3Res = 0;
   _jetPtGen = 0;
   _jetEtaGen = 0;
   _jetPhiGen = 0;
   _jetPtGenWithNu = 0;
   _puppijetPtGen = 0;
   _puppijetPtGenWithNu = 0;
   _jetJECuncty = 0;
   _puppijetJECuncty = 0;
   _lEta = 0;
   _lPhi = 0;
   _lPt = 0;
   _lPtcorr = 0;
   _lpdgId = 0;
   _lPassTightID = 0;
   _lPassLooseID = 0;
   _lgenPt = 0;
   _lgenEta = 0;
   _lgenPhi = 0;
   _lgenpdgId = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("_jetEta", &_jetEta, &b__jetEta);
   fChain->SetBranchAddress("_jetPhi", &_jetPhi, &b__jetPhi);
   if(useresiduals)fChain->SetBranchAddress("_jetPt", &_jetPt, &b__jetPt);
   else fChain->SetBranchAddress("_jetPtNoL2L3Res", &_jetPt, &b__jetPt);
   fChain->SetBranchAddress("_jetRawPt", &_jetRawPt, &b__jetRawPt);
   fChain->SetBranchAddress("_jetArea", &_jetArea, &b__jetArea);
   fChain->SetBranchAddress("_jetPassID", &_jetPassID, &b__jetPassID);
   fChain->SetBranchAddress("_jetLeptonPhotonCleaned", &_jetLeptonPhotonCleaned, &b__jetLeptonPhotonCleaned);
   fChain->SetBranchAddress("_jetPtNoL2L3Res", &_jetPtNoL2L3Res, &b__jetPtNoL2L3Res);
   fChain->SetBranchAddress("_jetPUMVAUpdate2017", &_jetPUMVAUpdate2017, &b__jetPUMVAUpdate2017);
   fChain->SetBranchAddress("_jetPUMVAUpdate2018", &_jetPUMVAUpdate2018, &b__jetPUMVAUpdate2018);
   fChain->SetBranchAddress("_puppijetEta", &_puppijetEta, &b__puppijetEta);
   fChain->SetBranchAddress("_puppijetPhi", &_puppijetPhi, &b__puppijetPhi);
   fChain->SetBranchAddress("_puppijetPt", &_puppijetPt, &b__puppijetPt);
   fChain->SetBranchAddress("_puppijetRawPt", &_puppijetRawPt, &b__puppijetRawPt);
   fChain->SetBranchAddress("_puppijetLeptonPhotonCleaned", &_puppijetLeptonPhotonCleaned, &b__puppijetLeptonPhotonCleaned);
   fChain->SetBranchAddress("_puppijetPassID", &_puppijetPassID, &b__puppijetPassID);
   fChain->SetBranchAddress("_puppijetPtNoL2L3Res", &_puppijetPtNoL2L3Res, &b__puppijetPtNoL2L3Res);
   fChain->SetBranchAddress("_jetPtGen", &_jetPtGen, &b__jetPtGen);
   fChain->SetBranchAddress("_jetEtaGen", &_jetEtaGen, &b__jetEtaGen);
   fChain->SetBranchAddress("_jetPhiGen", &_jetPhiGen, &b__jetPhiGen);
   fChain->SetBranchAddress("_jetPtGenWithNu", &_jetPtGenWithNu, &b__jetPtGenWithNu);
   fChain->SetBranchAddress("_puppijetPtGen", &_puppijetPtGen, &b__puppijetPtGen);
   fChain->SetBranchAddress("_puppijetPtGenWithNu", &_puppijetPtGenWithNu, &b__puppijetPtGenWithNu);
   fChain->SetBranchAddress("trueNVtx", &trueNVtx, &b_trueNVtx);
   fChain->SetBranchAddress("_genmet", &_genmet, &b__genmet);
   fChain->SetBranchAddress("_genmet_phi", &_genmet_phi, &b__genmet_phi);
   fChain->SetBranchAddress("_jetJECuncty", &_jetJECuncty, &b__jetJECuncty);
   fChain->SetBranchAddress("_puppijetJECuncty", &_puppijetJECuncty, &b__puppijetJECuncty);
   fChain->SetBranchAddress("_met", &_met, &b__met);
   fChain->SetBranchAddress("_met_phi", &_met_phi, &b__met_phi);
   fChain->SetBranchAddress("_puppimet", &_puppimet, &b__puppimet);
   fChain->SetBranchAddress("_puppimet_phi", &_puppimet_phi, &b__puppimet_phi);
   fChain->SetBranchAddress("_rawmet", &_rawmet, &b__rawmet);
   fChain->SetBranchAddress("_rawmet_phi", &_rawmet_phi, &b__rawmet_phi);
   fChain->SetBranchAddress("_puppirawmet", &_puppirawmet, &b__puppirawmet);
   fChain->SetBranchAddress("_puppirawmet_phi", &_puppirawmet_phi, &b__puppirawmet_phi);
   fChain->SetBranchAddress("_rawchsmet", &_rawchsmet, &b__rawchsmet);
   fChain->SetBranchAddress("_rawchsmet_phi", &_rawchsmet_phi, &b__rawchsmet_phi);
   fChain->SetBranchAddress("_chsmet", &_chsmet, &b__chsmet);
   fChain->SetBranchAddress("_chsmet_phi", &_chsmet_phi, &b__chsmet_phi);
   fChain->SetBranchAddress("_n_PV", &_n_PV, &b__n_PV);
   fChain->SetBranchAddress("_rho", &_rho, &b__rho);
   fChain->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
   fChain->SetBranchAddress("_runNb", &_runNb, &b__runNb);
   fChain->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
   fChain->SetBranchAddress("_ptll", &_ptll, &b__ptll);
   fChain->SetBranchAddress("_phill", &_phill, &b__phill);
   fChain->SetBranchAddress("_lEta", &_lEta, &b__lEta);
   fChain->SetBranchAddress("_lPhi", &_lPhi, &b__lPhi);
   fChain->SetBranchAddress("_lPt", &_lPt, &b__lPt);
   fChain->SetBranchAddress("_lPtcorr", &_lPtcorr, &b__lPtcorr);
   fChain->SetBranchAddress("_lpdgId", &_lpdgId, &b__lpdgId);
   fChain->SetBranchAddress("_lPassTightID", &_lPassTightID, &b__lPassTightID);
   fChain->SetBranchAddress("_lPassLooseID", &_lPassLooseID, &b__lPassLooseID);
   fChain->SetBranchAddress("_nEles", &_nEles, &b__nEles);
   fChain->SetBranchAddress("_nMus", &_nMus, &b__nMus);
   fChain->SetBranchAddress("_ptll_gen", &_ptll_gen, &b__ptll_gen);
   fChain->SetBranchAddress("_phill_gen", &_phill_gen, &b__phill_gen);
   fChain->SetBranchAddress("_lgenPt", &_lgenPt, &b__lgenPt);
   fChain->SetBranchAddress("_lgenEta", &_lgenEta, &b__lgenEta);
   fChain->SetBranchAddress("_lgenPhi", &_lgenPhi, &b__lgenPhi);
   fChain->SetBranchAddress("_lgenpdgId", &_lgenpdgId, &b__lgenpdgId);
   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf, &b_HLT_Ele35_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf, &b_HLT_Ele32_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf, &b_HLT_Ele27_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu24, &b_HLT_IsoTkMu24);
   fChain->SetBranchAddress("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, &b_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("_weight", &_weight, &b__weight);
   
   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("*met*",1);
   fChain->SetBranchStatus("_lpdgId",1);
   fChain->SetBranchStatus("_lPassTightID",1);
   fChain->SetBranchStatus("_lPassLooseID",1);
   fChain->SetBranchStatus("_lPt",1);
   fChain->SetBranchStatus("_lEta",1);
   fChain->SetBranchStatus("_lPhi",1);
   if(useresiduals) fChain->SetBranchStatus("_jetPt",1);
   else fChain->SetBranchStatus("_jetPtNoL2L3Res",1);
   fChain->SetBranchStatus("_jetPUMVAUpdate2017",1);
   fChain->SetBranchStatus("_jetPhi",1);
   fChain->SetBranchStatus("_jetEta",1);
   fChain->SetBranchStatus("_jetPassID",1);
   fChain->SetBranchStatus("_jetPtGen",1);
   fChain->SetBranchStatus("_jetJECuncty",1);

   fChain->SetBranchStatus("_lgenPt",1);
   fChain->SetBranchStatus("_lgenEta",1);
   fChain->SetBranchStatus("_lgenpdgId",1);
   fChain->SetBranchStatus("_lgenPhi",1);
   fChain->SetBranchStatus("HLT_*",1);
   fChain->SetBranchStatus("_n_PV",1);
   fChain->SetBranchStatus("_rho",1);
   fChain->SetBranchStatus("_eventNb",1);
   fChain->SetBranchStatus("trueNVtx",1);
   fChain->SetBranchStatus("_weight",1);




   Notify();
}

Bool_t templateforJERC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void templateforJERC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t templateforJERC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef templateforJERC_cxx
