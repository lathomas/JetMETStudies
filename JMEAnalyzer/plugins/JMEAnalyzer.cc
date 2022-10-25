// Original Author:  Laurent Thomas
//         Created:  Fri, 26 Apr 2019 12:51:46 GMT

// system include files

#include <memory>
#include <iostream>
#include <fstream>
#include <string>

// user include files
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 


#include "JetMETStudies/JMEAnalyzer/interface/Tools.h"
#include "RecoJets/JetProducers/plugins/PileupJetIdProducer.h" 
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"

//#include "JetMETStudies/JMEAnalyzer/python/RochesterCorrections/RoccoR.h"
#include "JetMETStudies/JMEAnalyzer/interface/RoccoR.h"
const int  N_METFilters=18;
enum METFilterIndex{
  idx_Flag_goodVertices,
  idx_Flag_globalTightHalo2016Filter,
  idx_Flag_globalSuperTightHalo2016Filter,
  idx_Flag_HBHENoiseFilter,
  idx_Flag_HBHENoiseIsoFilter,
  idx_Flag_EcalDeadCellTriggerPrimitiveFilter,
  idx_Flag_BadPFMuonFilter,
  idx_Flag_BadPFMuonDzFilter,
  idx_Flag_hfNoisyHitsFilter,
  idx_Flag_BadChargedCandidateFilter,
  idx_Flag_eeBadScFilter,
  idx_Flag_ecalBadCalibFilter,
  idx_Flag_ecalLaserCorrFilter,
  idx_Flag_EcalDeadCellBoundaryEnergyFilter,
  idx_PassecalBadCalibFilter_Update,
  idx_PassecalLaserCorrFilter_Update,
  idx_PassEcalDeadCellBoundaryEnergyFilter_Update,
  idx_PassBadChargedCandidateFilter_Update
};


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
using namespace edm;
using namespace std;
using namespace reco;
using namespace tools;


class JMEAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns>  {
   public:
      explicit JMEAnalyzer(const edm::ParameterSet&);
      ~JMEAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&);

  virtual bool PassSkim();
  virtual bool GetMETFilterDecision(const edm::Event& iEvent, edm::Handle<TriggerResults> METFilterResults, TString studiedfilter);
  virtual bool GetIdxFilterDecision(int it);
  virtual TString GetIdxFilterName(int it);
  virtual void InitandClearStuff();
  virtual void CalcDileptonInfo(const int& i, const int& j, Float_t & themass, Float_t & theptll, Float_t & thepzll,  Float_t & theyll, Float_t & thephill, Float_t & thedphill, Float_t & thecosthll);
  virtual void CalcDileptonInfoGen(const int& i, const int& j, Float_t & themass, Float_t & theptll, Float_t & thepzll,  Float_t & theyll, Float_t & thephill, Float_t & thedphill, Float_t & thecosthll);
  virtual int GetRecoIdx(const reco::GenJet * genjet , vector <Float_t> recojetpt,  vector <Float_t> recojeteta, vector <Float_t> recojetphi, vector <Float_t> recojetptgen );
  virtual bool IsTauCleaned(const pat::Jet *iJ);
  virtual bool IsLeptonPhotonCleaned(const pat::Jet *iJ);
  virtual bool IsLeptonPhotonCleaned(const reco::CaloJet *iJ);
  virtual bool PassJetPreselection(const pat::Jet * iJ, double genjetpt, bool ispuppi);
  virtual bool PassJetPreselection(const reco::CaloJet * iJ, double genjetpt);
  bool PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,const pat::Muon *muonit, const edm::Event&);
  bool PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,const pat::Electron *eleit, const edm::Event&);
  bool PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,const pat::Photon *photonit, const edm::Event&);
  bool PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,const pat::Jet *jetit, const edm::Event&);


  bool PassTriggerLeg(std::string triggerlegstring, const pat::Muon *muonit, const edm::Event& theevent){return PassTriggerLeg(triggerlegstring,"Noalttrigger",muonit,theevent);};
  bool PassTriggerLeg(std::string triggerlegstring, const pat::Electron *eleit, const edm::Event& theevent){ return PassTriggerLeg(triggerlegstring,"Noalttrigger",eleit,theevent);};
  bool PassTriggerLeg(std::string triggerlegstring, const pat::Photon *photonit, const edm::Event& theevent){ return PassTriggerLeg(triggerlegstring,"Noalttrigger",photonit,theevent);};
  bool PassTriggerLeg(std::string triggerlegstring, const pat::Jet *jetit, const edm::Event& theevent){ return PassTriggerLeg(triggerlegstring,"Noalttrigger",jetit,theevent);};



 
  // ----------member data ---------------------------
  edm::EDGetTokenT<TriggerResults> metfilterspatToken_; 
  edm::EDGetTokenT<TriggerResults> metfiltersrecoToken_; 
  edm::EDGetTokenT<bool> ecalBadCalibFilterUpdateToken_;
  edm::EDGetTokenT<bool> ecalLaserCorrFilterUpdateToken_;  
  edm::EDGetTokenT<bool> ecalDeadCellBoundaryEnergyFilterUpdateToken_;
  edm::EDGetTokenT<bool> badChargedCandidateFilterUpdateToken_;
  edm::EDGetTokenT<std::vector<Vertex> > verticesToken_; 
  edm::EDGetTokenT<double> rhoJetsToken_;
  edm::EDGetTokenT<double> rhoJetsNCToken_;


  edm::EDGetTokenT<std::vector< pat::Jet> > jetToken_;
  edm::EDGetTokenT<std::vector< pat::Jet> > jetAK8Token_;
  edm::EDGetTokenT<std::vector< pat::Jet> > jetPuppiToken_;
  edm::EDGetTokenT<std::vector< pat::Jet> > jetPuppiAK8Token_;
  edm::EDGetTokenT<std::vector< reco::CaloJet> > jetCaloToken_;
  edm::EDGetTokenT<std::vector< pat::Jet> > jetnoCHSToken_;

  edm::EDGetTokenT<std::vector< reco::GenJet> > genjetToken_;
  edm::EDGetTokenT<std::vector< reco::GenJet> > genAK8jetToken_;

  edm::EDGetTokenT<edm::ValueMap<float> > pileupJetIdDiscriminantUpdateToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > pileupJetIdDiscriminantUpdate2017Token_;
  edm::EDGetTokenT<edm::ValueMap<float> > pileupJetIdDiscriminantUpdate2018Token_;
  edm::EDGetTokenT<edm::ValueMap<StoredPileupJetIdentifier> > pileupJetIdVariablesUpdateToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > qgLToken_;

  edm::EDGetTokenT<std::vector< pat::PackedCandidate>> pfcandsToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiweightsToken_;
 
  edm::EDGetTokenT<std::vector< pat::MET> > metToken_;
  edm::EDGetTokenT<std::vector< pat::MET> > puppimetToken_;
  edm::EDGetTokenT<std::vector< pat::Electron> > electronToken_;
  edm::EDGetTokenT<std::vector< pat::Muon> > muonToken_;
  edm::EDGetTokenT<std::vector< pat::Tau> > tauToken_;

  edm::EDGetTokenT<std::vector< pat::Photon> > photonToken_;

  edm::EDGetTokenT<GenParticleCollection> genpartToken_;
  edm::EDGetTokenT<GenEventInfoProduct> geninfoToken_;
  edm::EDGetTokenT<LHEEventProduct> lheEventToken_;
  edm::EDGetTokenT<LHEEventProduct> lheEventALTToken_;

  edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetAssocCHSToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetWithNuAssocCHSToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetAssocPuppiToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetWithNuAssocPuppiToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetAssocCaloToken_;


  edm::EDGetTokenT<vector<PileupSummaryInfo> > puInfoToken_;

  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigobjectToken_;
  edm::EDGetTokenT<BXVector<GlobalAlgBlk>> l1GtToken_;
  edm::EDGetTokenT<l1t::MuonBxCollection>l1MuonToken_;
  edm::EDGetTokenT<l1t::EGammaBxCollection>l1EGammaToken_;
  edm::EDGetTokenT<l1t::JetBxCollection>l1JetToken_;
  edm::EDGetTokenT<GlobalExtBlkBxCollection> UnprefirableEventToken_;

  l1t::L1TGlobalUtil* l1GtUtils_;
  InputTag  algTag_, extTag_;
  //  edm::ESGetToken<JetCorrectorParametersCollection, JetCorrectionsRecord> jecsToken_; 


  Float_t JetPtCut_;
  Float_t AK8JetPtCut_;
  Float_t ElectronPtCut_;
  string ElectronVetoWP_, ElectronTightWP_, ElectronLooseWP_;
  Float_t MuonPtCut_;
  Float_t TauPtCut_;
  string RochCorrFile_;
  Float_t PhotonPtCut_;
  string PhotonTightWP_;
  Float_t PFCandPtCut_;

  Bool_t SaveTree_, IsMC_, SavePUIDVariables_, SaveAK8Jets_, SaveCaloJets_, SavenoCHSJets_, DropUnmatchedJets_, DropBadJets_, SavePFinJets_, ApplyPhotonID_;
  string Skim_;
  Bool_t Debug_;

  Float_t _PtCutPFforMultiplicity[6]={0,0.3,0.5,1,5,10};

  //Some histos to be saved for simple checks 
  TH1F *h_PFMet, *h_PuppiMet, *h_nvtx;

  //These two histos are used to count the nb of events processed and the simulated PU distribution
  //They should be filled before any skim is applied 
  TH1D *h_Counter, *h_trueNVtx;

  //The output TTree
  TTree* outputTree;
  TTree* jetPFTree;
  
  ifstream myfile_unprefevts;
  //Variables associated to leaves of the TTree

  unsigned long _eventNb;
  unsigned long _runNb;
  unsigned long _lumiBlock;
  unsigned long _bx;

  //Nb of primary vertices
  int _n_PV;
  Float_t _LV_x,_LV_y,_LV_z;
  Float_t _LV_errx,_LV_erry,_LV_errz;
  Float_t _PUV1_x,_PUV1_y,_PUV1_z;
  int trueNVtx;

  //Rho and RhoNC;
  Float_t _rho, _rhoNC;

  //MINIAOD original MET filters decisions
  bool Flag_goodVertices;
  bool Flag_globalTightHalo2016Filter;
  bool Flag_globalSuperTightHalo2016Filter;
  bool Flag_HBHENoiseFilter;
  bool Flag_HBHENoiseIsoFilter;
  bool Flag_EcalDeadCellTriggerPrimitiveFilter;
  bool Flag_BadPFMuonFilter;
  bool Flag_BadChargedCandidateFilter;
  bool Flag_eeBadScFilter;
  bool Flag_ecalBadCalibFilter;
  bool Flag_ecalLaserCorrFilter; 
  bool Flag_EcalDeadCellBoundaryEnergyFilter;

  bool Flag_BadPFMuonDzFilter;
  bool Flag_hfNoisyHitsFilter;
  //Decision obtained rerunning the filters on top of MINIAOD
  bool PassecalBadCalibFilter_Update;
  bool PassecalLaserCorrFilter_Update;  
  bool PassEcalDeadCellBoundaryEnergyFilter_Update;
  bool PassBadChargedCandidateFilter_Update;
  bool Flag_IsUnprefirable;

  //AK4 CHS Jets 
  vector<Float_t>  _jetEta;
  vector<Float_t>  _jetPhi;
  vector<Float_t>  _jetPt;
  vector<Float_t>  _jetRawPt;
  vector<Float_t>  _jet_CHEF;
  vector<Float_t>  _jet_NHEF;
  vector<Float_t>  _jet_NEEF;
  vector<Float_t>  _jet_CEEF;
  vector<Float_t>  _jet_MUEF;
  vector <int>  _jet_CHM;
  vector <int>  _jet_NHM;
  vector <int>  _jet_PHM;
  vector <int>  _jet_NM;
  vector <Float_t>  _jetArea;
  vector <bool> _jetPassID;
  vector <bool> _jetLeptonPhotonCleaned;
  vector <bool> _jetTauCleaned;
  
  vector <Float_t>  _jethfsigmaEtaEta;
  vector <Float_t>  _jethfsigmaPhiPhi;
  vector <Int_t> _jethfcentralEtaStripSize;
  vector <Int_t> _jethfadjacentEtaStripsSize;
  
  vector <Float_t>  _jetPtGen;
  vector <Float_t>  _jetEtaGen;
  vector <Float_t>  _jetPhiGen;
  vector <Float_t>  _jetPtGenWithNu;
  vector<Float_t>  _jetJECuncty;
  vector<Float_t>  _jetPUMVA; 
  vector<Float_t>  _jetPUMVAUpdate2017;
  vector<Float_t>  _jetPUMVAUpdate2018;
  vector<Float_t>  _jetPUMVAUpdate; 
  vector<Float_t>  _jetPtNoL2L3Res;
  vector<Float_t> _jet_corrjecs;
  vector<int> _jethadronFlavour;
  vector<int> _jetpartonFlavour;
  
  vector<Float_t> _jetDeepJet_b;
  vector<Float_t> _jetParticleNet_b;
  vector<Float_t> _jetDeepJet_c;
  vector<Float_t> _jetDeepJet_uds;
  vector<Float_t> _jetDeepJet_g;
  vector<Float_t> _jetQuarkGluonLikelihood;

  vector<Float_t>  _jet_beta ;
  vector<Float_t>  _jet_dR2Mean ;
  vector<Float_t>  _jet_majW ;
  vector<Float_t>  _jet_minW ;
  vector<Float_t>  _jet_frac01 ;
  vector<Float_t>  _jet_frac02 ;
  vector<Float_t>  _jet_frac03 ;
  vector<Float_t>  _jet_frac04 ;
  vector<Float_t>  _jet_ptD ;
  vector<Float_t>  _jet_betaStar ;
  vector<Float_t>  _jet_pull ;
  vector<Float_t>  _jet_jetR ;
  vector<Float_t>  _jet_jetRchg ;
  vector<int>  _jet_nParticles ;
  vector<int>  _jet_nCharged ;

  vector<Float_t>  _ak8jetEta ;
  vector<Float_t>  _ak8jetPhi ;
  vector<Float_t>  _ak8jetPt ;
  vector<Float_t>  _ak8jetArea ;
  vector<Float_t>  _ak8jetRawPt ;
  vector<Float_t>  _ak8jetPtGen ;


  //AK4 Puppi jets
  vector<Float_t>  _puppijetEta;
  vector<Float_t>  _puppijetPhi;
  vector<Float_t>  _puppijetPt;
  vector<Float_t>  _puppijetRawPt;
  vector<Float_t>  _puppijetJECuncty;
  vector <Float_t>  _puppijetPtGen;
  vector <Float_t>  _puppijetPtGenWithNu;
  vector<bool> _puppijetPassID;
  vector<bool> _puppijetLeptonPhotonCleaned;
  vector <Float_t>  _puppijetPtNoL2L3Res;

  //AK8 Puppi jets
  vector<Float_t>  _puppiak8jetEta ;
  vector<Float_t>  _puppiak8jetPhi ;
  vector<Float_t>  _puppiak8jetPt ;
  vector<Float_t>  _puppiak8jetRawPt ;
  vector<Float_t>  _puppiak8jetPtGen ;
  vector<Float_t>  _puppiak8jet_tau1 ;
  vector<Float_t>  _puppiak8jet_tau2 ;
  vector<Float_t>  _puppiak8jet_tau3 ;

  //AK4 Calo jets
  vector<Float_t>  _calojetEta;
  vector<Float_t>  _calojetPhi;
  vector<Float_t>  _calojetPt;
  vector<Float_t>  _calojetRawPt;
  vector <Float_t>  _calojetPtGen;
  vector<bool> _calojetLeptonPhotonCleaned;
  
  //Non CHS jets 
  vector<Float_t>  _noCHSjetEta;
  vector<Float_t>  _noCHSjetPhi;
  vector<Float_t>  _noCHSjetPt;
  vector<Float_t>  _noCHSjetRawPt;
  vector<Float_t>  _noCHSjetPtGen;
  vector<bool> _noCHSjetLeptonPhotonCleaned;

  //Gen AK4 jets
  vector<Float_t>  _genjetEta;
  vector<Float_t>  _genjetPhi;
  vector<Float_t>  _genjetPt;
  vector<Int_t>  _genjet_noCHSIdx;
  vector<Int_t>  _genjet_CaloIdx;
  vector<Int_t>  _genjet_CHSIdx;
  vector<Int_t>  _genjet_PuppiIdx;

  //Gen AK8 jets
  vector<Float_t>  _genAK8jetEta;
  vector<Float_t>  _genAK8jetPhi;
  vector<Float_t>  _genAK8jetPt;
  vector<Int_t>  _genAK8jet_CHSIdx;
  vector<Int_t>  _genAK8jet_PuppiIdx;


  //Leptons
  vector<Float_t>  _lEta;
  vector<Float_t>  _lPhi;
  vector<Float_t>  _lPt;
  vector<Float_t>  _lPtcorr;
  vector<Float_t>  _lPtSC;
  
  vector<Float_t>  _ldz;
  vector<Float_t>  _ldzError;
  vector<Float_t>  _ldxy;
  vector<Float_t>  _ldxyError;
  vector<Float_t>  _l3dIP;
  vector<Float_t>  _l3dIPError;

  vector<Bool_t>  _lpassHLT_IsoMu24;
  vector<Bool_t>  _lpassHLT_Ele32_WPTight_Gsf;

  vector<Bool_t>  _lPassTightID;
  vector<Bool_t>  _lPassLooseID;
  vector<Bool_t> _lisSAMuon ;
  vector<int> _lpdgId;
  int _nEles, _nMus;

  vector<Float_t>  _tauEta;
  vector<Float_t>  _tauPhi;
  vector<Float_t>  _tauPt;
  vector<Bool_t>  _tauPassMediumID;



  //Gen leptons
  vector<Float_t>  _lgenEta;
  vector<Float_t>  _lgenPhi;
  vector<Float_t>  _lgenPt;
  vector<int> _lgenpdgId;

  //Reco Photons
  vector<Float_t>  _phEta;
  vector<Float_t>  _phPhi;
  vector<Float_t>  _phPt;
  vector<Bool_t>  _phPassIso;
  vector<Float_t>  _phPtcorr;
  vector<Bool_t>  _phPassTightID;

  //Gen photons
  vector<Float_t>  _phgenEta;
  vector<Float_t>  _phgenPhi;
  vector<Float_t>  _phgenPt;
  
  //Event variables (reco) 
  //For dilepton studies
  Float_t _mll;
  Float_t _ptll;
  Float_t _pzll;
  Float_t _yll;
  Float_t _dphill;
  Float_t _phill;
  Float_t _costhCSll;
  Int_t _nElesll ;
  //For single photon studies
  Float_t _ptgamma;
  Float_t _phigamma;

  //Event variables (gen)
  Float_t _mll_gen;
  Float_t _ptll_gen;
  Float_t _pzll_gen;
  Float_t _yll_gen;
  Float_t _dphill_gen;
  Float_t _phill_gen;
  Float_t _costhCSll_gen;
  Int_t _ngenElesll ;
  Float_t _ptgamma_gen;
  Float_t _phigamma_gen;

  

  //Some gen variables
  Float_t _genHT, _weight;

  //PF candidates
  vector <Float_t> _PFcand_pt;
  vector <Float_t> _PFcand_eta;
  vector <Float_t> _PFcand_phi;
  vector <int> _PFcand_pdgId;
  vector <int> _PFcand_fromPV;
  vector <Float_t> _PFcand_dz;
  vector <Float_t> _PFcand_dzError;
  vector <Float_t> _PFcand_hcalFraction;
  vector <int> _PFcand_PVfitidx;
  vector <Float_t> _PFcand_puppiweight;
  

  //Nb of CH in PV fit and corresponding HT, for different pt cuts
  int _n_CH_fromvtxfit[6];
  Float_t _HT_CH_fromvtxfit[6];
  vector <Float_t> _METCH_PV;
  vector <Float_t> _METPhiCH_PV;
  vector <Float_t> _SumPT2CH_PV;
  vector <Float_t> _DztoLV_PV;
  
  //Nb of electrons in vtx fit
  int _n_PFele_fromvtxfit;
  int _n_PFmu_fromvtxfit;


  //Reco MET (PFME, PUPPIMET, CHS MET, T1 first and then Raw)
  Float_t _met;
  Float_t _met_phi;
  Float_t _puppimet;
  Float_t _puppimet_phi;
  Float_t _chsmet;
  Float_t _chsmet_phi;
  Float_t _rawmet;
  Float_t _rawmet_phi;
  Float_t _puppirawmet;
  Float_t _puppirawmet_phi;
  Float_t _rawchsmet;
  Float_t _rawchsmet_phi;

  //GenMET
  Float_t _genmet;
  Float_t _genmet_phi;

  //Triggers
  bool HLT_Photon110EB_TightID_TightIso;
  bool HLT_Photon165_R9Id90_HE10_IsoM;
  bool HLT_Photon120_R9Id90_HE10_IsoM;
  bool HLT_Photon90_R9Id90_HE10_IsoM;
  bool HLT_Photon75_R9Id90_HE10_IsoM;
  bool HLT_Photon50_R9Id90_HE10_IsoM;
  bool HLT_Photon200;
  bool HLT_Photon175;
  bool HLT_DiJet110_35_Mjj650_PFMET110;
  bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
  bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF;
  bool HLT_PFMET120_PFMHT120_IDTight_PFHT60;
  bool HLT_PFMET120_PFMHT120_IDTight;
  bool HLT_PFHT1050;
  bool HLT_PFHT900;
  bool HLT_PFJet500;
  bool HLT_AK8PFJet500;
  bool HLT_Ele35_WPTight_Gsf ;
  bool HLT_Ele32_WPTight_Gsf ;
  bool HLT_Ele27_WPTight_Gsf;
  bool HLT_IsoMu27;
  bool HLT_IsoMu24;
  bool HLT_IsoTkMu24;
  bool HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  bool HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
  bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
  bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
  bool HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
  bool HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  bool HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  bool HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
  bool HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
  bool HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
  bool HLT_PixelClusters_WP1;
  bool HLT_PixelClusters_WP2;
  bool HLT_PixelClusters_WP2_split;

  bool _l1prefire;


  //These are variables for a TTree where each entry is a jet and where one stores all PF cands
  vector <Float_t> _Jet_PFcand_pt;
  vector <Float_t> _Jet_PFcand_eta;
  vector <Float_t> _Jet_PFcand_phi;
  vector <int> _Jet_PFcand_pdgId;
  vector <int> _Jet_PFcand_fromPV;
  vector <Float_t> _Jet_PFcand_dz;
  vector <Float_t> _Jet_PFcand_dzError;
  Float_t _Jet_Pt;
  Float_t _Jet_PtGen;
  Float_t _Jet_Eta;
  Float_t _Jet_EtaGen;
  Float_t _Jet_Phi;
  Float_t _Jet_PhiGen;


  //Trigger matching variables: is a reco object matched to a trigger filter
  vector < bool >hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09;
  vector < bool >hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09;
  vector < bool >hltEle27WPTightGsfTrackIsoFilter; 
  vector < bool >hltEle35noerWPTightGsfTrackIsoFilter; 
  vector < bool >hltEle32WPTightGsfTrackIsoFilter; 
  vector < bool >hltEG175HEFilter; 
  vector < bool >hltEG200HEFilter; 
  vector < bool >hltEG165R9Id90HE10IsoMTrackIsoFilter; 
  vector < bool >hltEG120R9Id90HE10IsoMTrackIsoFilter; 
  vector < bool >hltEG90R9Id90HE10IsoMTrackIsoFilter; 
  vector < bool >hltEG75R9Id90HE10IsoMTrackIsoFilter; 
  vector < bool >hltEG50R9Id90HE10IsoMTrackIsoFilter; 
  vector < bool >hltEG110EBTightIDTightIsoTrackIsoFilter;
  vector < bool >hltSinglePFJet60; 
  vector < bool >hltSinglePFJet80; 
  vector < bool >hltSinglePFJet140 ;
  vector < bool >hltSinglePFJet200; 
  vector < bool >hltSinglePFJet260; 
  vector < bool >hltSinglePFJet320; 
  vector < bool >hltSinglePFJet400; 
  vector < bool >hltSinglePFJet450; 
  vector < bool >hltSinglePFJet500; 
  vector<bool> hltHIPhoton40Eta3p1;
  vector<bool> hltEle17WPLoose1GsfTrackIsoFilterForHI;
  vector<bool> hltEle15WPLoose1GsfTrackIsoFilterForHI;
  vector<bool> hltL3fL1sMu10lqL1f0L2f10L3Filtered12;
  vector<bool> hltL3fL1sMu10lqL1f0L2f10L3Filtered15;

  bool passL1_Initial_bxmin1[512];
  bool passL1_Initial_bx0[512];
  bool passL1_Initial_bxplus1[512];

  //L1 muon
  vector <int> _L1mu_Qual;
  vector <Float_t> _L1mu_pt;
  vector <Float_t> _L1mu_eta;
  vector <Float_t> _L1mu_phi;
  vector <int> _L1mu_TFIdx;
  vector <int> _L1mu_bx;
  
  vector <int> _L1mu_dXY;
  vector <Float_t> _L1mu_upt;
  vector <int> _L1mu_charge;


  //L1 eg
  vector <Float_t> _L1eg_pt;
  vector <Float_t> _L1eg_eta;
  vector <Float_t> _L1eg_phi;
  vector <int> _L1eg_bx;
  vector <int> _L1eg_iso;

  //L1 jet
  vector <Float_t> _L1jet_pt;
  vector <Float_t> _L1jet_eta;
  vector <Float_t> _L1jet_phi;
  vector <int> _L1jet_bx;


  //Rochester correction (for muons)
  RoccoR rc; 
  //JEC uncertainties
  //  JetCorrectionUncertainty *jecUnc; 
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
JMEAnalyzer::JMEAnalyzer(const edm::ParameterSet& iConfig)
 :
  metfilterspatToken_(consumes<TriggerResults>(iConfig.getParameter<edm::InputTag>("METFiltersPAT"))),
  metfiltersrecoToken_(consumes<TriggerResults>(iConfig.getParameter<edm::InputTag>("METFiltersRECO"))),
  ecalBadCalibFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("ECALBadCalibFilterUpdate"))),
  ecalLaserCorrFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("ECALLaserCorrFilterUpdate"))),
  ecalDeadCellBoundaryEnergyFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("ECALDeadCellBoundaryEnergyFilterUpdate"))),
  badChargedCandidateFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilterUpdate"))),
  verticesToken_(consumes<std::vector<Vertex> > (iConfig.getParameter<edm::InputTag>("Vertices"))),
  rhoJetsToken_(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll",""))),
  rhoJetsNCToken_(consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral",""))),
  jetToken_(consumes< std::vector< pat::Jet> >(iConfig.getParameter<edm::InputTag>("Jets"))),
  jetAK8Token_(consumes< std::vector< pat::Jet> >(iConfig.getParameter<edm::InputTag>("JetsAK8"))),
  jetPuppiToken_(consumes< std::vector< pat::Jet> >(iConfig.getParameter<edm::InputTag>("JetsPuppi"))),
  jetPuppiAK8Token_(consumes< std::vector< pat::Jet> >(iConfig.getParameter<edm::InputTag>("JetsPuppiAK8"))),
  jetCaloToken_(consumes< std::vector< reco::CaloJet> >(iConfig.getParameter<edm::InputTag>("JetsCalo"))),
  jetnoCHSToken_(consumes< std::vector< pat::Jet> >(iConfig.getParameter<edm::InputTag>("JetsPFnoCHS"))),
  genjetToken_(consumes< std::vector< reco::GenJet> >(iConfig.getParameter<edm::InputTag>("GenJets"))),
  genAK8jetToken_(consumes< std::vector< reco::GenJet> >(iConfig.getParameter<edm::InputTag>("GenAK8Jets"))),
  pileupJetIdDiscriminantUpdateToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pileupJetIdDiscriminantUpdate"))),
  pileupJetIdDiscriminantUpdate2017Token_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pileupJetIdDiscriminantUpdate2017"))),
  pileupJetIdDiscriminantUpdate2018Token_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pileupJetIdDiscriminantUpdate2018"))),
  pileupJetIdVariablesUpdateToken_(consumes<edm::ValueMap<StoredPileupJetIdentifier> >(iConfig.getParameter<edm::InputTag>("pileupJetIdVariablesUpdate"))),
  qgLToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("QuarkGluonLikelihood"))),
  pfcandsToken_(consumes<std::vector< pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
  puppiweightsToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("PuppiWeights"))),
  metToken_(consumes<std::vector<pat::MET> > (iConfig.getParameter<edm::InputTag>("PFMet"))),
  puppimetToken_(consumes<std::vector<pat::MET> > (iConfig.getParameter<edm::InputTag>("PuppiMet"))),
  electronToken_(consumes< std::vector< pat::Electron> >(iConfig.getParameter<edm::InputTag>("Electrons"))),
  muonToken_(consumes< std::vector< pat::Muon> >(iConfig.getParameter<edm::InputTag>("Muons"))),
  tauToken_(consumes< std::vector< pat::Tau> >(iConfig.getParameter<edm::InputTag>("Taus"))),
  photonToken_(consumes< std::vector< pat::Photon> >(iConfig.getParameter<edm::InputTag>("Photons"))),
  genpartToken_(consumes<GenParticleCollection> (iConfig.getParameter<edm::InputTag>("GenParticles"))),
  geninfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("GenInfo"))),
  lheEventToken_(consumes<LHEEventProduct> ( iConfig.getParameter<InputTag>("LHELabel"))),
  lheEventALTToken_(consumes<LHEEventProduct> ( iConfig.getParameter<InputTag>("LHELabelALT"))),
  genJetAssocCHSToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("GenJetMatchCHS"))),
  genJetWithNuAssocCHSToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("GenJetWithNuMatchCHS"))),
  genJetAssocPuppiToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("GenJetMatchPuppi"))),
  genJetWithNuAssocPuppiToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("GenJetWithNuMatchPuppi"))),
  genJetAssocCaloToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("GenJetMatchCalo"))),
  puInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PULabel"))),
  trgresultsToken_(consumes<TriggerResults>(iConfig.getParameter<edm::InputTag>("Triggers"))),
  trigobjectToken_(consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("slimmedPatTrigger"))),
  l1GtToken_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("l1GtSrc"))),
  l1MuonToken_(consumes<l1t::MuonBxCollection>(edm::InputTag("gmtStage2Digis","Muon"))),
  l1EGammaToken_(consumes<l1t::EGammaBxCollection>(edm::InputTag("caloStage2Digis","EGamma"))),
  l1JetToken_(consumes<l1t::JetBxCollection>(edm::InputTag("caloStage2Digis","Jet"))),
  UnprefirableEventToken_(consumes<GlobalExtBlkBxCollection>(edm::InputTag("simGtExtUnprefireable"))),
//  jecsToken_(esConsumes(edm::ESInputTag("","AK4PFchs"))),
  JetPtCut_(iConfig.getParameter<double>("JetPtCut")),
  AK8JetPtCut_(iConfig.getParameter<double>("AK8JetPtCut")),
  ElectronPtCut_(iConfig.getParameter<double>("ElectronPtCut")),
  ElectronVetoWP_(iConfig.getParameter<string>("ElectronVetoWorkingPoint")),
  ElectronTightWP_(iConfig.getParameter<string>("ElectronTightWorkingPoint")),
  ElectronLooseWP_(iConfig.getParameter<string>("ElectronLooseWorkingPoint")),
  MuonPtCut_(iConfig.getParameter<double>("MuonPtCut")),
  TauPtCut_(iConfig.getParameter<double>("TauPtCut")),
  RochCorrFile_(iConfig.getParameter<string>("RochCorrFile")),
  PhotonPtCut_(iConfig.getParameter<double>("PhotonPtCut")),
  PhotonTightWP_(iConfig.getParameter<string>("PhotonTightWorkingPoint")),
  PFCandPtCut_(iConfig.getParameter<double>("PFCandPtCut")),
  SaveTree_(iConfig.getParameter<bool>("SaveTree")), 
  IsMC_(iConfig.getParameter<bool>("IsMC")),
  SavePUIDVariables_(iConfig.getParameter<bool>("SavePUIDVariables")),
  SaveAK8Jets_(iConfig.getParameter<bool>("SaveAK8Jets")),
  SaveCaloJets_(iConfig.getParameter<bool>("SaveCaloJets")),
  SavenoCHSJets_(iConfig.getParameter<bool>("SavenoCHSJets")),
  DropUnmatchedJets_(iConfig.getParameter<bool>("DropUnmatchedJets")),
  DropBadJets_(iConfig.getParameter<bool>("DropBadJets")),
  SavePFinJets_(iConfig.getParameter<bool>("SavePFinJets")),
  ApplyPhotonID_(iConfig.getParameter<bool>("ApplyPhotonID")),
  Skim_(iConfig.getParameter<string>("Skim")),
  Debug_(iConfig.getParameter<bool>("Debug"))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs; 
  h_nvtx  = fs->make<TH1F>("h_nvtx" , "Number of reco vertices;N_{vtx};Events"  ,    100, 0., 100.);
  h_PFMet  = fs->make<TH1F>("h_PFMet" , "Type 1 PFMET (GeV);Type 1 PFMET (GeV);Events"  ,    1000, 0., 5000.);
  h_PuppiMet  = fs->make<TH1F>("h_PuppiMet" , "PUPPI MET (GeV);PUPPI MET (GeV);Events"  ,    1000, 0., 5000.);
  h_Counter = fs->make<TH1D>("h_Counter", "Events counter", 5,0,5);
  h_trueNVtx = fs->make<TH1D>("h_trueNVtx", "Nb of generated vertices", 200,0,200);

  outputTree = fs->make<TTree>("tree","tree");
  if(SavePFinJets_) jetPFTree = fs->make<TTree>("jetPFtree","tree");

  
  rc.init(RochCorrFile_);
  algTag_ = iConfig.getParameter<edm::InputTag>("l1GtSrc") ; 
  extTag_ = iConfig.getParameter<edm::InputTag>("l1GtSrc");
  l1GtUtils_ =new l1t::L1TGlobalUtil(iConfig, consumesCollector(), *this, algTag_, extTag_, l1t::UseEventSetupIn::Event);
}


JMEAnalyzer::~JMEAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void JMEAnalyzer::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  //JES uncties: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties

  //  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  //iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);


  //The following lines are needed for picking events from a list (currently unprefirable events)
  int runnb = iRun.id().run();
  TString str_runnb =Form("%d", runnb) ; 
  myfile_unprefevts.open("UnprefireableEventList/run_"+str_runnb+".txt");
  cout << "Run is " <<str_runnb <<endl;
}


void JMEAnalyzer::endRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{

  myfile_unprefevts.close();
  //  delete jecUnc;
}

// ------------ method called for each event  ------------
void
JMEAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  if(Debug_) std::cout <<"Entering analyze" <<std::endl;
  /*
  if(auto handle = iSetup.getHandle(jecsToken_)) {
    auto const& JetCorParColl = *handle;
    JetCorrectorParameters const & JetCorPar = (JetCorParColl)["Uncertainty"];
    jecUnc = new JetCorrectionUncertainty(JetCorPar);
  }
  else jecUnc = 0;
  */


  
  InitandClearStuff();
  
  
  _runNb = iEvent.id().run();
  _eventNb = iEvent.id().event();
  _lumiBlock = iEvent.luminosityBlock();
  _bx=iEvent.bunchCrossing();

  
  //Gen info, needed pre skim !! 
  edm::Handle<GenEventInfoProduct> GenInfoHandle;
  iEvent.getByToken(geninfoToken_, GenInfoHandle);
  if(GenInfoHandle.isValid())_weight=GenInfoHandle->weight();
  else _weight = 0;
  h_Counter->Fill(0.,_weight);
  /*for(unsigned int a =0; a< (GenInfoHandle->binningValues()).size();a++){
    double thebinning = GenInfoHandle->hasBinningValues() ? (GenInfoHandle->binningValues())[a] : 0.0 ;
    }*/
  
  //Pile up info at gen level
  Handle<std::vector<PileupSummaryInfo> > puInfo;
  iEvent.getByToken(puInfoToken_, puInfo);
  if(puInfo.isValid()){
  vector<PileupSummaryInfo>::const_iterator pvi;
  for (pvi = puInfo->begin(); pvi != puInfo->end(); ++pvi) {
    if (pvi->getBunchCrossing() == 0) trueNVtx = pvi->getTrueNumInteractions();
  }
  }
  else trueNVtx = -1.;
  h_trueNVtx->Fill(trueNVtx);

  if(Debug_) std::cout <<"Entering analyze 2" <<std::endl;
  //Triggers 
  edm::Handle<TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  if( !trigResults.failedToGet() ) {
    int N_Triggers = trigResults->size();
    const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      if (trigResults.product()->accept(i_Trig)) {
	TString TrigPath =trigName.triggerName(i_Trig);
	if(TrigPath.Contains("HLT_Photon110EB_TightID_TightIso_v"))HLT_Photon110EB_TightID_TightIso =true;
	if(TrigPath.Contains("HLT_Photon165_R9Id90_HE10_IsoM_v"))HLT_Photon165_R9Id90_HE10_IsoM =true;
	if(TrigPath.Contains("HLT_Photon120_R9Id90_HE10_IsoM_v"))HLT_Photon120_R9Id90_HE10_IsoM =true;
	if(TrigPath.Contains("HLT_Photon90_R9Id90_HE10_IsoM_v"))HLT_Photon90_R9Id90_HE10_IsoM =true;
	if(TrigPath.Contains("HLT_Photon75_R9Id90_HE10_IsoM_v"))HLT_Photon75_R9Id90_HE10_IsoM =true;
	if(TrigPath.Contains("HLT_Photon50_R9Id90_HE10_IsoM_v"))HLT_Photon50_R9Id90_HE10_IsoM =true;
	if(TrigPath.Contains("HLT_Photon200_v"))HLT_Photon200 =true;
	if(TrigPath.Contains("HLT_Photon175_v"))HLT_Photon175 =true;
	if(TrigPath.Contains("HLT_DiJet110_35_Mjj650_PFMET110_v"))HLT_DiJet110_35_Mjj650_PFMET110 =true;
	if(TrigPath.Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v"))HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 =true;
	if(TrigPath.Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"))HLT_PFMETNoMu120_PFMHTNoMu120_IDTight =true;
	if(TrigPath.Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF_v"))HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF =true;
	if(TrigPath.Contains("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v"))HLT_PFMET120_PFMHT120_IDTight_PFHT60 =true;
	if(TrigPath.Contains("HLT_PFMET120_PFMHT120_IDTight_v"))HLT_PFMET120_PFMHT120_IDTight =true;
	if(TrigPath.Contains("HLT_PFHT1050_v"))HLT_PFHT1050 =true;
	if(TrigPath.Contains("HLT_PFHT900_v"))HLT_PFHT900 =true;
	if(TrigPath.Contains("HLT_PFJet500_v"))HLT_PFJet500 =true;
	if(TrigPath.Contains("HLT_AK8PFJet500_v"))HLT_AK8PFJet500 =true;
	if(TrigPath.Contains("HLT_Ele35_WPTight_Gsf_v"))HLT_Ele35_WPTight_Gsf =true;
	if(TrigPath.Contains("HLT_Ele32_WPTight_Gsf_v"))HLT_Ele32_WPTight_Gsf =true;
	if(TrigPath.Contains("HLT_Ele27_WPTight_Gsf_v"))HLT_Ele27_WPTight_Gsf =true;
	if(TrigPath.Contains("HLT_IsoMu27_v"))HLT_IsoMu27 =true;
	if(TrigPath.Contains("HLT_IsoMu24_v"))HLT_IsoMu24 =true;
	if(TrigPath.Contains("HLT_IsoTkMu24_v"))HLT_IsoTkMu24 =true;
	if(TrigPath.Contains("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"))HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ =true;
	if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"))HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ =true;
	if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"))HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL =true;
	if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"))HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ =true;
	if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"))HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 =true;
	if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL =true;
	if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"))HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ =true;
	if(TrigPath.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"))HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ =true;
	if(TrigPath.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"))HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ =true;
	if(TrigPath.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"))HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL =true;
	if(TrigPath.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"))HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL =true;

	if(TrigPath.Contains("HLT_PixelClusters_WP1_v"))HLT_PixelClusters_WP1 =true;
	if(TrigPath.Contains("HLT_PixelClusters_WP2_v"))HLT_PixelClusters_WP2 =true;
	if(TrigPath.Contains("HLT_PixelClusters_WP2_split_v"))HLT_PixelClusters_WP2_split =true;
      }
    }
  }

  if(Debug_) std::cout <<"Entering analyze 3" <<std::endl;
  //Prefiring, see: https://github.com/nsmith-/PrefireAnalysis/#usage
  edm::Handle<BXVector<GlobalAlgBlk>> l1GtHandle;
  iEvent.getByToken(l1GtToken_, l1GtHandle);


  l1GtUtils_->retrieveL1(iEvent, iSetup, l1GtToken_);
  const std::vector<std::pair<std::string, bool>> initialDecisions = l1GtUtils_->decisionsInitial();
  
  /*
  for(int i = 0; i < 20; i++){
    std::cout << "L1 algo: " << initialDecisions[i].first<<std::endl;
    }*/

  _l1prefire= false;
  if(!IsMC_)_l1prefire = l1GtHandle->begin(-1)->getFinalOR();


  for(int i =0; i <512; i++){
    if(!IsMC_) passL1_Initial_bxmin1[i]= l1GtHandle->begin(-1)->getAlgoDecisionInitial(i);
    else passL1_Initial_bxmin1[i]= false;
    passL1_Initial_bx0[i]= l1GtHandle->begin(0)->getAlgoDecisionInitial(i);
    if(!IsMC_) passL1_Initial_bxplus1[i]= l1GtHandle->begin(+1)->getAlgoDecisionInitial(i);
    else passL1_Initial_bxmin1[i]= false;
  }




  if(Debug_) std::cout <<"Entering analyze 4" <<std::endl;
  edm::Handle<l1t::MuonBxCollection> l1muoncoll;
  iEvent.getByToken(l1MuonToken_ , l1muoncoll);
  //  const  int nbx = IsMC_ ? 0:2;
  for(int i = l1muoncoll->getFirstBX() ; i<= l1muoncoll->getLastBX() ;i++){
    for( l1t::MuonBxCollection::const_iterator l1muonit= l1muoncoll->begin(i); l1muonit != l1muoncoll->end(i) ; ++l1muonit){
      if(l1muonit->pt()<1) continue;
      _L1mu_Qual.push_back( l1muonit->hwQual() );
      _L1mu_pt.push_back( l1muonit->pt() );
      _L1mu_eta.push_back( l1muonit->eta() );
      _L1mu_phi.push_back( l1muonit->phi() );
      _L1mu_bx.push_back( i);
      _L1mu_dXY.push_back( l1muonit->hwDXY() );
      _L1mu_upt.push_back( l1muonit->ptUnconstrained() );
      _L1mu_charge.push_back(l1muonit->charge());
      _L1mu_TFIdx.push_back( l1muonit->tfMuonIndex());
    }
  }
  
  if(Debug_) std::cout <<"Entering analyze 5" <<std::endl;

  edm::Handle<l1t::EGammaBxCollection> l1egcoll;
  iEvent.getByToken(l1EGammaToken_ , l1egcoll);
  //  const  int nbx = IsMC_ ? 0:2;
  for(int i = l1egcoll->getFirstBX() ; i<= l1egcoll->getLastBX() ;i++){
    for( l1t::EGammaBxCollection::const_iterator l1egit= l1egcoll->begin(i); l1egit != l1egcoll->end(i) ; ++l1egit){
      if(l1egit->pt()<5) continue;
      _L1eg_pt.push_back( l1egit->pt() );
      _L1eg_eta.push_back( l1egit->eta() );
      _L1eg_phi.push_back( l1egit->phi() );
      _L1eg_bx.push_back( i);
      _L1eg_iso.push_back( l1egit->hwIso());

      //      cout << "L1EG pt, eta, phi: "<<l1egit->pt() <<", "<<l1egit->eta() <<", " <<l1egit->phi() <<", "  <<endl;
      //cout << "L1EG hwIso, hwQual, hwEta:  "<< l1egit->hwIso() <<", " <<l1egit->hwQual() <<", " <<l1egit->hwEta() <<", " <<endl; 
    }
  }


  edm::Handle<l1t::JetBxCollection> l1jetcoll;
  iEvent.getByToken(l1JetToken_ , l1jetcoll);
  //  const  int nbx = IsMC_ ? 0:2;
  for(int i = l1jetcoll->getFirstBX() ; i<= l1jetcoll->getLastBX() ;i++){
    for( l1t::JetBxCollection::const_iterator l1jetit= l1jetcoll->begin(i); l1jetit != l1jetcoll->end(i) ; ++l1jetit){
      if(l1jetit->pt()<25) continue;
      _L1jet_pt.push_back( l1jetit->pt() );
      _L1jet_eta.push_back( l1jetit->eta() );
      _L1jet_phi.push_back( l1jetit->phi() );
      _L1jet_bx.push_back( i);
    }
  }




  //Unprefirable
  Flag_IsUnprefirable = false;
  edm::Handle<GlobalExtBlkBxCollection> handleUnprefEventResults;
  iEvent.getByToken(UnprefirableEventToken_, handleUnprefEventResults);
  if(handleUnprefEventResults.isValid()){
    if (handleUnprefEventResults->size() != 0) {
      Flag_IsUnprefirable = handleUnprefEventResults->at(0, 0).getExternalDecision(GlobalExtBlk::maxExternalConditions - 1);
    }
  }
  if(Debug_) std::cout <<"Entering analyze 6" <<std::endl;
  
  //Vertices
  // NB there are usually two available collection of vertices in MINIAOD
  edm::Handle<std::vector<Vertex> > theVertices;
  iEvent.getByToken(verticesToken_,theVertices) ;
  _n_PV = theVertices->size();
  Vertex::Point PV(0,0,0);
  if(_n_PV){ PV = theVertices->begin()->position();}


  //Storing leading vertex coordinates, as well as leading PU vertex (primarily for bad vertex choice studies)
  if(Debug_) std::cout <<"Entering analyze 7" <<std::endl;
  const Vertex* LVtx = &((*theVertices)[0]);
  _LV_x = LVtx->x();
  _LV_y = LVtx->y();
  _LV_z = LVtx->z();
  _LV_errx = LVtx->xError();
  _LV_erry = LVtx->yError();
  _LV_errz = LVtx->zError();

  _PUV1_x = 0;
  _PUV1_y = 0;
  _PUV1_z = 0;
  
  if(_n_PV>=2) {
    _PUV1_x = ((*theVertices)[1]).x();
    _PUV1_y = ((*theVertices)[1]).y();
    _PUV1_z = ((*theVertices)[1]).z();
  }

  //Storing information from all vertex to study vertex sorting
  vector <TVector2 > metPV;
  for(unsigned int i = 0;i < theVertices->size(); i++){
    const Vertex* PVtx = &((*theVertices)[i]);
    TVector2 dummyT2V(0.,0.);
    metPV.push_back(dummyT2V);
    _METCH_PV.push_back(-99.);
    _METPhiCH_PV.push_back(-99.);
    _SumPT2CH_PV.push_back(0.);
    _DztoLV_PV.push_back( PVtx->z()- LVtx->z());
  }

  //Rho
  edm::Handle<double> rhoJets;
  iEvent.getByToken(rhoJetsToken_,rhoJets);
  _rho = *rhoJets;

  //Rho Neutral Central 
  edm::Handle<double> rhoJetsNC;
  iEvent.getByToken(rhoJetsNCToken_,rhoJetsNC);
  _rhoNC = *rhoJetsNC;


  if(Debug_) std::cout <<"Entering analyze 8" <<std::endl;  
  //MET filters are stored in TriggerResults::RECO or TriggerResults::PAT . Should take the latter if it exists
  edm::Handle<TriggerResults> METFilterResults;
  iEvent.getByToken(metfilterspatToken_, METFilterResults);
  if(!(METFilterResults.isValid())) iEvent.getByToken(metfiltersrecoToken_, METFilterResults);
  
  Flag_goodVertices= GetMETFilterDecision(iEvent,METFilterResults,"Flag_goodVertices");
  if(!Flag_goodVertices && Skim_=="L1Study")return;
  Flag_globalTightHalo2016Filter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_globalTightHalo2016Filter");
  Flag_globalSuperTightHalo2016Filter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_globalSuperTightHalo2016Filter");
  Flag_HBHENoiseFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_HBHENoiseFilter");
  Flag_HBHENoiseIsoFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_HBHENoiseIsoFilter");
  Flag_EcalDeadCellTriggerPrimitiveFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_EcalDeadCellTriggerPrimitiveFilter");
  Flag_BadPFMuonFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_BadPFMuonFilter");
  Flag_BadPFMuonDzFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_BadPFMuonDzFilter");
  Flag_BadChargedCandidateFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_BadChargedCandidateFilter");
  Flag_eeBadScFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_eeBadScFilter");
  Flag_ecalBadCalibFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_ecalBadCalibFilter");
  Flag_EcalDeadCellBoundaryEnergyFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_EcalDeadCellBoundaryEnergyFilter");
  Flag_ecalLaserCorrFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_ecalLaserCorrFilter");
  Flag_hfNoisyHitsFilter = GetMETFilterDecision(iEvent,METFilterResults,"Flag_hfNoisyHitsFilter");

  //Now accessing the decisions of some filters that we reran on top of MINIAOD
  edm::Handle<bool> handle_PassecalBadCalibFilter_Update ;
  iEvent.getByToken(ecalBadCalibFilterUpdateToken_,handle_PassecalBadCalibFilter_Update);
  if(handle_PassecalBadCalibFilter_Update.isValid()) PassecalBadCalibFilter_Update =  (*handle_PassecalBadCalibFilter_Update );
  else{ 
    if(Debug_) std::cout <<"handle_PassecalBadCalibFilter_Update.isValid() =false" <<endl;
    PassecalBadCalibFilter_Update = true;
  }
  
  edm::Handle<bool> handle_PassecalLaserCorrFilter_Update ;
  iEvent.getByToken(ecalLaserCorrFilterUpdateToken_,handle_PassecalLaserCorrFilter_Update);
  if(handle_PassecalLaserCorrFilter_Update.isValid())PassecalLaserCorrFilter_Update =  (*handle_PassecalLaserCorrFilter_Update );
  else{ 
    if(Debug_) std::cout <<"handle_PassecalLaserCorrFilter_Update.isValid() =false" <<endl;
    PassecalLaserCorrFilter_Update = true;
  }

  edm::Handle<bool> handle_PassEcalDeadCellBoundaryEnergyFilter_Update;
  iEvent.getByToken(ecalDeadCellBoundaryEnergyFilterUpdateToken_,handle_PassEcalDeadCellBoundaryEnergyFilter_Update);
  if(handle_PassEcalDeadCellBoundaryEnergyFilter_Update.isValid())PassEcalDeadCellBoundaryEnergyFilter_Update =  (*handle_PassEcalDeadCellBoundaryEnergyFilter_Update );
  else{  
    if(Debug_) std::cout <<"handle_PassEcalDeadCellBoundaryEnergyFilter_Update.isValid =false" <<endl; 
    PassEcalDeadCellBoundaryEnergyFilter_Update = true;
  }

  edm::Handle<bool> handle_PassBadChargedCandidateFilter_Update;
  iEvent.getByToken(badChargedCandidateFilterUpdateToken_,handle_PassBadChargedCandidateFilter_Update);
  if(handle_PassBadChargedCandidateFilter_Update.isValid())PassBadChargedCandidateFilter_Update =  (*handle_PassBadChargedCandidateFilter_Update );
  else{  
    if(Debug_) std::cout <<"handle_PassBadChargedCandidateFilter_Update.isValid =false" <<endl; 
    PassBadChargedCandidateFilter_Update = true;
  }


  //Electrons
  edm::Handle< std::vector<pat::Electron> > thePatElectrons;
  iEvent.getByToken(electronToken_,thePatElectrons);
  for( std::vector<pat::Electron>::const_iterator electron = (*thePatElectrons).begin(); electron != (*thePatElectrons).end(); electron++ ) {
    if((&*electron)->pt() <3) continue; //Loose cut  on uncorrected pt 
    //Implementing smearing/scaling EGM corrections
    //See here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2#Applying_the_Energy_Scale_and_sm and here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaUL2016To2018
    double ptelecorr = (&*electron)->pt();
    if((&*electron)->hasUserFloat("ecalEnergyPostCorr") ){
         ptelecorr = ptelecorr * (&*electron)->userFloat("ecalTrkEnergyPostCorr") /  (&*electron)->energy() ;
    }
    
    bool passvetoid = (&*electron)->electronID(ElectronVetoWP_)&& (&*electron)->pt()>3;
    if(!passvetoid) continue;
    //Counting the number of electrons, not all of them will be stored
    _nEles++;
    if((&*electron)->pt()<ElectronPtCut_)continue;    
    _lEta.push_back((&*electron)->eta());
    _lPhi.push_back((&*electron)->phi());
    _lPt.push_back((&*electron)->pt());
    _lPtcorr.push_back(ptelecorr );
    _lPtSC.push_back( ((&*electron)->superCluster()->energy()) / cosh((&*electron)->superCluster()->eta())   );
    _lpdgId.push_back(-11*(&*electron)->charge());
    _lPassTightID.push_back( (&*electron)->electronID(ElectronTightWP_) );
    _lPassLooseID.push_back( (&*electron)->electronID(ElectronLooseWP_) );
    _lisSAMuon.push_back( false); 

    _ldz.push_back( (&*electron)->gsfTrack()->dz(PV));
    _ldzError.push_back( (&*electron)->gsfTrack()->dzError());
    _ldxy.push_back( (&*electron)->gsfTrack()->dxy(PV));
    _ldxyError.push_back( (&*electron)->gsfTrack()->dxyError());
    _l3dIP.push_back( (&*electron)->dB(pat::Electron::PV3D));
    _l3dIPError.push_back( (&*electron)->edB(pat::Electron::PV3D));

    _lpassHLT_IsoMu24.push_back(false);
    //if((&*electron)->pt()>30)std::cout << "Pt, Eta, Phi, pass Ele32 " <<  ((&*electron)->triggerObjectMatchByPath ("HLT_Ele32_WPTight_Gsf_v*",true,true) != nullptr) << ", "<<   ((&*electron)->triggerObjectMatchByPath ("HLT_Ele32_WPTight_Gsf_v",true,true) != nullptr)  <<", " <<  ((&*electron)->triggerObjectMatchByPath ("HLT_Ele35_WPTight_Gsf_v*",true,true) != nullptr)<<", " << ((&*electron)->triggerObjectMatchByPath ("HLT_Ele32_WPTight_Gsf_v*",true,false) != nullptr) <<", " << ((&*electron)->triggerObjectMatchByPath ("HLT_Ele32_WPTight_Gsf_v*",false,false) != nullptr)  <<std::endl; 
    _lpassHLT_Ele32_WPTight_Gsf.push_back( PassTriggerLeg("hltEle32WPTightGsfTrackIsoFilter",&*electron,iEvent) );

    hltEle27WPTightGsfTrackIsoFilter.push_back(PassTriggerLeg("hltEle27WPTightGsfTrackIsoFilter",&*electron,iEvent));
    hltEle35noerWPTightGsfTrackIsoFilter.push_back(PassTriggerLeg("hltEle35noerWPTightGsfTrackIsoFilter",&*electron,iEvent));
    hltEle32WPTightGsfTrackIsoFilter.push_back(PassTriggerLeg("hltEle32WPTightGsfTrackIsoFilter",&*electron,iEvent));
    hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09.push_back(false);
    hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09.push_back(false);

    hltL3fL1sMu10lqL1f0L2f10L3Filtered12.push_back(false);
    hltL3fL1sMu10lqL1f0L2f10L3Filtered15.push_back(false);
    
    hltHIPhoton40Eta3p1.push_back(PassTriggerLeg("hltHIPhoton40Eta3p1",&*electron,iEvent));
    hltEle17WPLoose1GsfTrackIsoFilterForHI.push_back(PassTriggerLeg("hltEle17WPLoose1GsfTrackIsoFilterForHI",&*electron,iEvent));
    hltEle15WPLoose1GsfTrackIsoFilterForHI.push_back(PassTriggerLeg("hltEle15WPLoose1GsfTrackIsoFilterForHI",&*electron,iEvent));  
  }
  
  //Muons
  edm::Handle< std::vector<pat::Muon> > thePatMuons;
  iEvent.getByToken(muonToken_,thePatMuons);
  for( std::vector<pat::Muon>::const_iterator muon = (*thePatMuons).begin(); muon != (*thePatMuons).end(); muon++ ) {
    if((&*muon)->pt() <0) continue; //Loose cut  on uncorrected pt 

    //Rochester corrections: https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon#Rochester_Correction
    //https://indico.cern.ch/event/926898/contributions/3897122/attachments/2052816/3441285/roccor.pdf
    double ptmuoncorr= (&*muon)->pt();

    if( !IsMC_) ptmuoncorr *= rc.kScaleDT( (&*muon)->charge(),  (&*muon)->pt(), (&*muon)->eta(),(&*muon)->phi());
    else{
      if( (&*muon)->genLepton() !=0)  ptmuoncorr *= rc.kSpreadMC( (&*muon)->charge(),  (&*muon)->pt(), (&*muon)->eta(),(&*muon)->phi(), (&*muon)->genLepton()->pt() );
      else if(! ((&*muon)->innerTrack()).isNull()) ptmuoncorr *= rc.kSmearMC( (&*muon)->charge(),  (&*muon)->pt(), (&*muon)->eta(),(&*muon)->phi(),  (&*muon)->innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());
    }
    
    bool passvetoid=  (&*muon)->passed(reco::Muon::CutBasedIdLoose)&& (&*muon)->passed(reco::Muon::PFIsoVeryLoose)&&(&*muon)->pt()>10 ;  
    passvetoid=  (&*muon)->isStandAloneMuon() || (&*muon)->pt()>5;
    if(!passvetoid) continue;
    //Counting the number of muons, not all of them will be stored 
    _nMus++;
    if((&*muon)->pt()<MuonPtCut_)continue;
    //std::cout << "Muon pt, eta, phi, pass IsoMu24" << (&*muon)->pt()<<", " << (&*muon)->eta()<<", " << (&*muon)->phi() <<", "<<(&*muon)->triggered("HLT_IsoMu24_v*")  <<std::endl;
    _lEta.push_back((&*muon)->eta());
    _lPhi.push_back((&*muon)->phi());
    _lPt.push_back((&*muon)->pt());
    _lPtcorr.push_back( ptmuoncorr );
    _lPtSC.push_back( -1 );
    _lpdgId.push_back(-13*(&*muon)->charge());
    _lPassTightID.push_back(  (&*muon)->passed(reco::Muon::CutBasedIdMediumPrompt )&& (&*muon)->passed(reco::Muon::PFIsoTight ) );
    _lPassLooseID.push_back(  (&*muon)->passed(reco::Muon::CutBasedIdLoose )&& (&*muon)->passed(reco::Muon::PFIsoLoose ) );
    _lisSAMuon.push_back( (&*muon)->isStandAloneMuon());
    if( !((&*muon)->innerTrack()).isNull()){
      _ldz.push_back( (&*muon)->innerTrack()->dz(PV));
      _ldzError.push_back( (&*muon)->innerTrack()->dzError());
      _ldxy.push_back( (&*muon)->innerTrack()->dxy(PV));
      _ldxyError.push_back( (&*muon)->innerTrack()->dxyError());
      _l3dIP.push_back( (&*muon)->dB(pat::Muon::PV3D));
      _l3dIPError.push_back((&*muon)->edB(pat::Muon::PV3D));
    }
    else{
      _ldz.push_back( 0.);
      _ldzError.push_back(0.);
      _ldxy.push_back(0.);
      _ldxyError.push_back(0.);
      _l3dIP.push_back(0.);
      _l3dIPError.push_back(0.);
    }
    _lpassHLT_IsoMu24.push_back((&*muon)->triggered("HLT_IsoMu24_v*"));
    _lpassHLT_Ele32_WPTight_Gsf.push_back(false);

    hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09.push_back(PassTriggerLeg("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09","hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07",&*muon,iEvent));
    hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09.push_back(PassTriggerLeg("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09","hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07",&*muon,iEvent));
    hltEle27WPTightGsfTrackIsoFilter.push_back(false);
    hltEle35noerWPTightGsfTrackIsoFilter.push_back(false);
    hltEle32WPTightGsfTrackIsoFilter.push_back(false);

    hltL3fL1sMu10lqL1f0L2f10L3Filtered12.push_back(PassTriggerLeg("hltL3fL1sMu10lqL1f0L2f10L3Filtered12",&*muon,iEvent));
    hltL3fL1sMu10lqL1f0L2f10L3Filtered15.push_back(PassTriggerLeg("hltL3fL1sMu10lqL1f0L2f10L3Filtered15",&*muon,iEvent));
    hltHIPhoton40Eta3p1.push_back(false);
    hltEle17WPLoose1GsfTrackIsoFilterForHI.push_back(false);
    hltEle15WPLoose1GsfTrackIsoFilterForHI.push_back(false);

  }




  edm::Handle< std::vector<pat::Tau> > thePatTaus;
  iEvent.getByToken(tauToken_,thePatTaus);
  for( std::vector<pat::Tau>::const_iterator tau = (*thePatTaus).begin(); tau != (*thePatTaus).end(); tau++ ) {
    if((&*tau)->pt()<TauPtCut_)continue;
    _tauEta.push_back((&*tau)->eta());
    _tauPhi.push_back((&*tau)->phi());
    _tauPt.push_back((&*tau)->pt());
    _tauPassMediumID.push_back((&*tau)->tauID("byMediumDeepTau2017v2p1VSjet"));
  }


  //Photons
  //These two variables store the pt/phi of the photon in monophoton events
  _ptgamma=0;
  _phigamma=0;
  edm::Handle< std::vector<pat::Photon> > thePatPhotons;
  iEvent.getByToken(photonToken_,thePatPhotons);
  for( std::vector<pat::Photon>::const_iterator photon = (*thePatPhotons).begin(); photon != (*thePatPhotons).end(); photon++ ) {
    //    cout << "Photon pt, eta, phi " << (&*photon)->pt()<< ", " <<(&*photon)->eta()<<", "<<(&*photon)->phi()<<endl;
    if((&*photon)->pt() <5) continue; //Loose cut  on uncorrected pt 
    double ptphotoncorr = (&*photon)->pt(); //This is (possibly) the smeared/scaled pt for data ! 
    //Implementing smearing/scaling EGM corrections
    //See here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2#Applying_the_Energy_Scale_and_sm and here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaUL2016To2018
    //For |eta|>2.5 the corrections are not adapted and one should instead pick the uncorrected value:
    if((&*photon)->hasUserFloat("ecalEnergyPostCorr") && fabs((&*photon)->eta())<2.5){
         ptphotoncorr = ptphotoncorr * (&*photon)->userFloat("ecalEnergyPostCorr") /  (&*photon)->energy() ;
    }
    
    if(ptphotoncorr <PhotonPtCut_)continue;
    bool passtightid = (&*photon)->photonID(PhotonTightWP_) && (&*photon)->passElectronVeto()&& !((&*photon)->hasPixelSeed()  ) ;
    if(!passtightid&& ApplyPhotonID_) continue;
    //Barrel photons + R9 cut to ensure a precise energy measurement
    if(_ptgamma==0&& passtightid  && fabs((&*photon)->eta())<1.4442&& (&*photon)->r9()>0.9 ) {
      _ptgamma=(&*photon)->pt();
      _phigamma=(&*photon)->phi();
    }
    _phEta.push_back((&*photon)->eta());
    _phPhi.push_back((&*photon)->phi());
    _phPt.push_back( (&*photon)->pt());
    _phPtcorr.push_back( ptphotoncorr);
    _phPassTightID.push_back(passtightid);
    
    
    bool passmediumiso =false;
    // bool passmediumid =false; 

    double eta=(&*photon)->superCluster()->eta();
    double pt=(&*photon)->pt();
    double chIso = (&*photon)->chargedHadronIso();
    double nhIso = (&*photon)->neutralHadronIso();
    double gIso = (&*photon)->photonIso();
    bool barrel = (fabs(eta) <= 1.479);
    bool endcap = (!barrel && fabs(eta) < 32.5);

    //double sigmaIetaIeta = (&*photon)->full5x5_sigmaIetaIeta();
    //double hoe = (&*photon)->hadTowOverEm();

    double chArea(0.),nhArea(0.),gArea(0.);
    if(fabs(eta)<1.0){ chArea = 0.0385;nhArea= 0.0636 ;gArea=0.1240;}
    else if(fabs(eta)<1.479){ chArea = 0.0468 ;nhArea=0.1103 ;gArea=0.1093;}
    else if(fabs(eta)<2.0){ chArea = 0.0435 ;nhArea=0.0759 ;gArea=0.0631;}
    else if(fabs(eta) < 2.2){ chArea = 0.0378 ;nhArea=0.0236 ;gArea=0.0779;}
    else if(fabs(eta) < 2.3){ chArea = 0.0338 ;nhArea=0.0151 ;gArea=0.0999;}
    else if(fabs(eta) < 2.4){ chArea = 0.0314 ;nhArea=0.00007 ;gArea=0.1155;}
    else { chArea =0.0269 ;nhArea=0.0132 ;gArea=0.1373 ;}

    /*    if ( barrel
             && hoe < 0.035
	     && sigmaIetaIeta < 0.0103)  passmediumid = true;*/
    if ( barrel
	 && TMath::Max(chIso-chArea*_rho,0.0) < 1.416
	 && TMath::Max(nhIso-nhArea*_rho,0.0) < 2.491 + 0.0126*pt + 0.000026*pt*pt
	 && TMath::Max(gIso-gArea*_rho, 0.0) < 2.952 + 0.0035*pt ) passmediumiso = true;
    /*    if ( endcap
             && hoe <  0.027
	     && sigmaIetaIeta < 0.0271) passmediumid = true;*/
    if ( endcap
	 && TMath::Max(chIso-chArea*_rho,0.0) < 1.012
	 && TMath::Max(nhIso-nhArea*_rho,0.0) <  9.131 + 0.0119* pt+ 0.000025*pt*pt
	 && TMath::Max(gIso-gArea*_rho, 0.0) < 4.095 + 0.0040*pt ) passmediumiso = true;
    
    _phPassIso.push_back(passmediumiso);

    hltEG175HEFilter.push_back(PassTriggerLeg("hltEG175HEFilter",&*photon,iEvent));
    hltEG200HEFilter.push_back(PassTriggerLeg("hltEG200HEFilter",&*photon,iEvent));
    hltEG165R9Id90HE10IsoMTrackIsoFilter.push_back(PassTriggerLeg("hltEG165R9Id90HE10IsoMTrackIsoFilter",&*photon,iEvent));
    hltEG120R9Id90HE10IsoMTrackIsoFilter.push_back(PassTriggerLeg("hltEG120R9Id90HE10IsoMTrackIsoFilter",&*photon,iEvent));
    hltEG90R9Id90HE10IsoMTrackIsoFilter.push_back(PassTriggerLeg("hltEG90R9Id90HE10IsoMTrackIsoFilter",&*photon,iEvent));
    hltEG75R9Id90HE10IsoMTrackIsoFilter.push_back(PassTriggerLeg("hltEG75R9Id90HE10IsoMTrackIsoFilter",&*photon,iEvent));
    hltEG50R9Id90HE10IsoMTrackIsoFilter.push_back(PassTriggerLeg("hltEG50R9Id90HE10IsoMTrackIsoFilter",&*photon,iEvent));
    hltEG110EBTightIDTightIsoTrackIsoFilter.push_back(PassTriggerLeg("hltEG110EBTightIDTightIsoTrackIsoFilter",&*photon,iEvent));
  }
  
  
  //Compute dilepton variables
  Float_t mll(0),ptll(0),pzll(0),yll(0),phill(0),dphill(0),costhCSll(0);
  for(unsigned int i = 0; i < _lPt.size(); i++){
    if(_lPt.size() !=2) continue;
    if(_lPt[i]<5) continue;
      if(!_lPassTightID[i])  continue;
      if(fabs(_lpdgId[i]) !=11 && fabs(_lpdgId[i])!=13 ) continue;
      for(unsigned int j = 0; j < i; j++){
        if(_lPt[j]<5) continue;
        if(!_lPassTightID[j])  continue;
        if(fabs(_lpdgId[j]) !=11 && fabs(_lpdgId[j])!=13 ) continue;
        //if( _lpdgId[i] != -_lpdgId[j]  ) continue;
        CalcDileptonInfo(i,j, mll,ptll,pzll,yll,phill,dphill,costhCSll);
	_mll= mll; _ptll=ptll; _pzll=pzll; _yll=yll; _phill=phill; _dphill=dphill; _costhCSll=costhCSll;
	if(fabs(_lpdgId[j]) ==11 && fabs(_lpdgId[j])  ==11) _nElesll =2 ;
	else if(fabs(_lpdgId[j]) ==13 && fabs(_lpdgId[j])  ==13) _nElesll =0 ;
	else _nElesll = 1;
      }
  }
    

  //  if(!PassSkim()) return;

    
  //Jets
  
  edm::Handle< std::vector< pat::Jet> > theJets;
  iEvent.getByToken(jetToken_,theJets );

  //Value map for the recalculated PU ID
  edm::Handle<edm::ValueMap<float> > pileupJetIdDiscriminantUpdate;
  //Value map to access the recalculated variables entering PU ID
  edm::Handle<edm::ValueMap<StoredPileupJetIdentifier> > pileupJetIdVariablesUpdate;
  edm::Handle<edm::ValueMap<float> > pileupJetIdDiscriminantUpdate2017;
  edm::Handle<edm::ValueMap<float> > pileupJetIdDiscriminantUpdate2018;

  //Value map for Quark/gluon likelihood
  edm::Handle<edm::ValueMap<float> > quarkgluonlikelihood;

  //  Accessing the matching with an updated gen collection
  edm::Handle<edm::Association<reco::GenJetCollection>> genJetMatch;
  iEvent.getByToken(genJetAssocCHSToken_, genJetMatch);
   
  edm::Handle<edm::Association<reco::GenJetCollection>> genJetWithNuMatch;
  iEvent.getByToken(genJetWithNuAssocCHSToken_, genJetWithNuMatch);

  Float_t leadjetpt (0.);
  bool useupdategenjets = true;
  if(theJets.isValid()){
    for( std::vector<pat::Jet>::const_iterator jet = (*theJets).begin(); jet != (*theJets).end(); jet++ ) {
      edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>( theJets , jet -theJets->begin()));
      //As of today only gen jets with pt >8 are saved in MINIAOD, see: https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/slimming/slimmedGenJets_cfi.py
      //The gen jets are defined excluding neutrinos. 
      //In order to decrease this pt cut and/or include neutrinos, one needs to recluster gen jets. 
      //Boolean to use the updated gen jet collection or the default one. Should probably be made configurable at some point.
      
      //      if(abs((&*jet)->eta())>3 ) cout <<"LOL pt, eta, phi"<<  (&*jet)->pt() <<", " <<  (&*jet)->eta() <<", " << (&*jet)->phi() <<endl;

      const reco::GenJet * updatedgenjet =0;
      const reco::GenJet * updatedgenjetwithnu =0; 
      if(genJetMatch.isValid() && genJetWithNuMatch.isValid() && useupdategenjets){
	updatedgenjet = ( (*genJetMatch)[jetRef].isNonnull() && (*genJetMatch)[jetRef].isAvailable()) ? &*(*genJetMatch)[jetRef] : 0;
	updatedgenjetwithnu = ( (*genJetWithNuMatch)[jetRef].isNonnull() && (*genJetWithNuMatch)[jetRef].isAvailable()) ? &*(*genJetWithNuMatch)[jetRef] : 0;
      }
 
      const reco::GenJet * genjet = updatedgenjet != 0 ? updatedgenjet : (&*jet) ->genJet()  ;
      Float_t jetptgen(-99.), jetetagen(-99.),jetphigen(-99.);
      Float_t jetptgenwithnu(-99.);
      if( genjet !=0 ){
	jetptgen= genjet->pt() ;
	jetetagen= genjet->eta() ;
	jetphigen= genjet->phi() ;
      }
      
      if(!PassJetPreselection((&*jet) , jetptgen , false) ) continue;
	 
      bool passid = PassJetID(  (&*jet) ,"2018");
      bool isleptoncleaned =IsLeptonPhotonCleaned ((&*jet) );
      if( DropBadJets_ && (!passid || !isleptoncleaned )  ) continue;//Drop bad jets (mostly leptons).
      if( genjet ==0   && DropUnmatchedJets_ && (&*jet)->pt()<50 ) continue;//Drop genunmatched jets (mostly PU). Keep those with pt>50 as these probably require special attention.

      if((&*jet)->pt() >leadjetpt) leadjetpt = (&*jet)->pt();
      _jetEta.push_back((&*jet)->eta());
      _jetPhi.push_back((&*jet)->phi());
      _jetPt.push_back((&*jet)->pt());
      _jet_CHEF.push_back((&*jet)->chargedHadronEnergyFraction());
      _jet_NHEF.push_back((&*jet)->neutralHadronEnergyFraction() );
      _jet_NEEF.push_back((&*jet)->neutralEmEnergyFraction() );
      _jet_CEEF.push_back((&*jet)->chargedEmEnergyFraction() );
      _jet_MUEF.push_back((&*jet)->muonEnergyFraction() );
      _jet_CHM.push_back((&*jet)->chargedMultiplicity());
      _jet_NHM.push_back((&*jet)->neutralHadronMultiplicity());
      _jet_PHM.push_back((&*jet)->photonMultiplicity());
      _jet_NM.push_back((&*jet)->neutralMultiplicity());
      _jetArea.push_back((&*jet)->jetArea());
      _jetPassID.push_back(passid);
      _jetLeptonPhotonCleaned.push_back(isleptoncleaned);
      _jetTauCleaned.push_back( IsTauCleaned((&*jet)));

      //Accessing the default PU ID stored in MINIAOD https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
      if((&*jet)->hasUserFloat("hfJetShowerShape:sigmaEtaEta")) _jethfsigmaEtaEta.push_back((&*jet)->userFloat("hfJetShowerShape:sigmaEtaEta"));
      if((&*jet)->hasUserFloat("hfJetShowerShape:sigmaPhiPhi")) _jethfsigmaPhiPhi.push_back((&*jet)->userFloat("hfJetShowerShape:sigmaPhiPhi"));
      if((&*jet)->hasUserInt("hfJetShowerShape:centralEtaStripSize")) _jethfcentralEtaStripSize.push_back((&*jet)->userInt("hfJetShowerShape:centralEtaStripSize"));
      if((&*jet)->hasUserInt("hfJetShowerShape:adjacentEtaStripsSize")) _jethfadjacentEtaStripsSize.push_back((&*jet)->userInt("hfJetShowerShape:adjacentEtaStripsSize"));
      if((&*jet)->hasUserFloat("pileupJetId:fullDiscriminant") )_jetPUMVA.push_back( (&*jet)->userFloat("pileupJetId:fullDiscriminant") );
      else _jetPUMVA.push_back(-99.);
      //Accessing the recomputed PU ID. This must be done with a value map. 
      iEvent.getByToken(pileupJetIdDiscriminantUpdateToken_,pileupJetIdDiscriminantUpdate);
      if(pileupJetIdDiscriminantUpdate.isValid()) _jetPUMVAUpdate.push_back((*pileupJetIdDiscriminantUpdate)[jetRef] );
      else  _jetPUMVAUpdate.push_back(-1 );
      iEvent.getByToken(pileupJetIdDiscriminantUpdate2017Token_,pileupJetIdDiscriminantUpdate2017);
      if(pileupJetIdDiscriminantUpdate2017.isValid()) _jetPUMVAUpdate2017.push_back((*pileupJetIdDiscriminantUpdate2017)[jetRef] );
      else  _jetPUMVAUpdate2017.push_back(-1 );
      iEvent.getByToken(pileupJetIdDiscriminantUpdate2018Token_,pileupJetIdDiscriminantUpdate2018);
      if(pileupJetIdDiscriminantUpdate2018.isValid()) _jetPUMVAUpdate2018.push_back((*pileupJetIdDiscriminantUpdate2018)[jetRef] );
      else  _jetPUMVAUpdate2018.push_back(-1 );
      
      //Accessing the recomputed input variables to the PUID BDT
      iEvent.getByToken(pileupJetIdVariablesUpdateToken_,pileupJetIdVariablesUpdate);
      if(pileupJetIdVariablesUpdate.isValid()){
	StoredPileupJetIdentifier pujetidentifier = (*pileupJetIdVariablesUpdate)[jetRef] ;
	
	_jet_beta.push_back(pujetidentifier.beta());
	_jet_dR2Mean.push_back(pujetidentifier.dR2Mean());
	_jet_majW.push_back(pujetidentifier.majW());
	_jet_minW.push_back(pujetidentifier.minW());
	_jet_frac01.push_back(pujetidentifier.frac01());
	_jet_frac02.push_back(pujetidentifier.frac02());
	_jet_frac03.push_back(pujetidentifier.frac03());
	_jet_frac04.push_back(pujetidentifier.frac04());
	_jet_ptD.push_back(pujetidentifier.ptD());
	_jet_betaStar.push_back(pujetidentifier.betaStar());
	_jet_pull.push_back(pujetidentifier.pull());
	_jet_jetR.push_back(pujetidentifier.jetR());
	_jet_jetRchg.push_back(pujetidentifier.jetRchg());
	_jet_nParticles.push_back(pujetidentifier.nParticles());
	_jet_nCharged.push_back(pujetidentifier.nCharged());
	
      }
      else if(Debug_) cout << "PUID variables are not valid"<<endl;

      // Parton flavour (gen level)
      _jethadronFlavour.push_back((&*jet)->hadronFlavour());  
      _jetpartonFlavour.push_back((&*jet)->partonFlavour());   
      
      //Flavour tagging (reco)
      //Deep Jet https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
      _jetDeepJet_b.push_back(  (&*jet)->bDiscriminator("pfDeepFlavourJetTags:probb")+ (&*jet)->bDiscriminator("pfDeepFlavourJetTags:probbb") + (&*jet)->bDiscriminator("pfDeepFlavourJetTags:problepb") );
      _jetDeepJet_c.push_back( (&*jet)->bDiscriminator("pfDeepFlavourJetTags:probc") );
      _jetDeepJet_uds.push_back( (&*jet)->bDiscriminator("pfDeepFlavourJetTags:probuds")  );
      _jetDeepJet_g.push_back(  (&*jet)->bDiscriminator("pfDeepFlavourJetTags:probg")  );
      _jetParticleNet_b.push_back(  (&*jet)->bDiscriminator("pfParticleNetAK4JetTags:probb")+ (&*jet)->bDiscriminator("pfParticleNetAK4JetTags:probbb"));

      //Quark Gluon likelihood  https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuarkGluonLikelihood
      iEvent.getByToken(qgLToken_, quarkgluonlikelihood);
      if(quarkgluonlikelihood.isValid() )_jetQuarkGluonLikelihood.push_back( (*quarkgluonlikelihood)[jetRef] );
      else _jetQuarkGluonLikelihood.push_back( -1.);
      

      _jetRawPt.push_back( (&*jet)->correctedP4("Uncorrected").Pt() );
      _jetPtNoL2L3Res.push_back( (&*jet)->correctedP4("L3Absolute") .Pt() ); 
      _jet_corrjecs.push_back((&*jet)->pt() / (&*jet)->correctedP4("Uncorrected").Pt() );

      /*jecUnc->setJetEta((&*jet)->eta());
      jecUnc->setJetPt((&*jet)->pt());
      _jetJECuncty.push_back( jecUnc->getUncertainty(true) );
      */

      if(updatedgenjetwithnu !=0) jetptgenwithnu = updatedgenjetwithnu->pt() ;
      _jetPtGen.push_back(jetptgen);
      _jetEtaGen.push_back(jetetagen);
      _jetPhiGen.push_back(jetphigen);
      _jetPtGenWithNu.push_back(jetptgenwithnu);
    
      //Looping over PF candidates
      for (unsigned i = 0; i < jet->numberOfSourceCandidatePtrs(); ++i) {
	const reco::Candidate* icand = jet->sourceCandidatePtr(i).get();
	const pat::PackedCandidate* lPack = dynamic_cast<const pat::PackedCandidate*>(icand);
	_Jet_PFcand_pt.push_back(lPack->pt());
	_Jet_PFcand_eta.push_back(lPack->eta());
	_Jet_PFcand_phi.push_back(lPack->phi());
	_Jet_PFcand_pdgId.push_back(lPack->pdgId());
	_Jet_PFcand_fromPV.push_back(lPack->fromPV());
	_Jet_PFcand_dz.push_back(lPack->dz());
	if(lPack-> hasTrackDetails())   _Jet_PFcand_dzError.push_back(lPack->dzError());
	else  _Jet_PFcand_dzError.push_back(-99);
      }
      
      _Jet_Pt = (&*jet)->pt();
      _Jet_Eta = (&*jet)->eta();
      _Jet_Phi = (&*jet)->phi();
      _Jet_PtGen = jetptgen;
      _Jet_EtaGen = jetetagen;
      _Jet_PhiGen = jetphigen;
      
      if(SavePFinJets_)jetPFTree->Fill();
      _Jet_PFcand_pt.clear();
      _Jet_PFcand_eta.clear();
      _Jet_PFcand_phi.clear();
      _Jet_PFcand_pdgId.clear();
      _Jet_PFcand_fromPV.clear();
      _Jet_PFcand_dz.clear();
      _Jet_PFcand_dzError.clear();
      
     
      hltSinglePFJet60.push_back(PassTriggerLeg("hltSinglePFJet60",&*jet,iEvent));
      hltSinglePFJet80.push_back(PassTriggerLeg("hltSinglePFJet80",&*jet,iEvent));
      hltSinglePFJet140.push_back(PassTriggerLeg("hltSinglePFJet140",&*jet,iEvent));
      hltSinglePFJet200.push_back(PassTriggerLeg("hltSinglePFJet200",&*jet,iEvent));
      hltSinglePFJet260.push_back(PassTriggerLeg("hltSinglePFJet260",&*jet,iEvent));
      hltSinglePFJet320.push_back(PassTriggerLeg("hltSinglePFJet320",&*jet,iEvent));
      hltSinglePFJet400.push_back(PassTriggerLeg("hltSinglePFJet400",&*jet,iEvent));
      hltSinglePFJet450.push_back(PassTriggerLeg("hltSinglePFJet450",&*jet,iEvent));
      hltSinglePFJet500.push_back(PassTriggerLeg("hltSinglePFJet500",&*jet,iEvent));
 
    }
  }
  else if(Debug_){cout << "Invalid jet collection"<<endl;}
  


  //Now PUPPI jets
  
  edm::Handle< std::vector< pat::Jet> > thePuppiJets;
  iEvent.getByToken(jetPuppiToken_,thePuppiJets );

  //  Accessing the matching with an updated gen collection
  edm::Handle<edm::Association<reco::GenJetCollection>> genJetMatchPuppi;
  iEvent.getByToken(genJetAssocPuppiToken_, genJetMatchPuppi);
   
  edm::Handle<edm::Association<reco::GenJetCollection>> genJetWithNuMatchPuppi;
  iEvent.getByToken(genJetWithNuAssocPuppiToken_, genJetWithNuMatchPuppi);

  if(thePuppiJets.isValid()){
    for( std::vector<pat::Jet>::const_iterator jet = (*thePuppiJets).begin(); jet != (*thePuppiJets).end(); jet++ ) {
      edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>( thePuppiJets , jet -thePuppiJets->begin()));
      
      const reco::GenJet * updatedgenjet =0;
      const reco::GenJet * updatedgenjetwithnu =0;
      if(genJetMatchPuppi.isValid() && genJetWithNuMatchPuppi.isValid() && useupdategenjets){
        updatedgenjet = ( (*genJetMatchPuppi)[jetRef].isNonnull() && (*genJetMatchPuppi)[jetRef].isAvailable()) ? &*(*genJetMatchPuppi)[jetRef] : 0;
        updatedgenjetwithnu = ( (*genJetWithNuMatchPuppi)[jetRef].isNonnull() && (*genJetWithNuMatchPuppi)[jetRef].isAvailable()) ? &*(*genJetWithNuMatchPuppi)[jetRef] : 0;
      }
      const reco::GenJet * genjet = useupdategenjets?updatedgenjet: (&*jet) ->genJet()  ;

      Float_t jetptgen(-99.);//, jetetagen(-99.),jetphigen(-99.); 
      Float_t jetptgenwithnu(-99.);
      if( genjet !=0 ){
        jetptgen= genjet->pt() ;
	//jetetagen= genjet->eta() ;
	//jetphigen= genjet->phi() ;
      }
      if(updatedgenjetwithnu !=0) jetptgenwithnu = updatedgenjetwithnu->pt() ;

      if(!PassJetPreselection((&*jet) , jetptgen , true) ) continue;

      if( genjet ==0   && DropUnmatchedJets_ && (&*jet)->pt()<50 ) continue;
      bool isleptoncleaned =IsLeptonPhotonCleaned ((&*jet) );
      if( DropBadJets_ && ( !isleptoncleaned )  ) continue;
      _puppijetLeptonPhotonCleaned.push_back(isleptoncleaned);
      //Still need to update the jet id for PUPPI
      _puppijetPassID.push_back(true);
      _puppijetEta.push_back((&*jet)->eta());
      _puppijetPhi.push_back((&*jet)->phi());
      _puppijetPt.push_back((&*jet)->pt());
      _puppijetRawPt.push_back((&*jet)-> correctedP4("Uncorrected").Pt());
      _puppijetPtNoL2L3Res.push_back( (&*jet)->correctedP4("L3Absolute") .Pt() );

      /*      jecUnc->setJetEta((&*jet)->eta());
      jecUnc->setJetPt((&*jet)->pt());
      _puppijetJECuncty.push_back( jecUnc->getUncertainty(true) );
      */
      if(updatedgenjetwithnu !=0) jetptgenwithnu = updatedgenjetwithnu->pt() ;
      _puppijetPtGen.push_back(jetptgen);
      _puppijetPtGenWithNu.push_back(jetptgenwithnu);
    }
  }
  

  //AK8 Puppi jets
  edm::Handle< std::vector< pat::Jet> > thePuppiAK8Jets;
  iEvent.getByToken(jetPuppiAK8Token_,thePuppiAK8Jets );
  if(thePuppiAK8Jets.isValid()){
    for( std::vector<pat::Jet>::const_iterator jet = (*thePuppiAK8Jets).begin(); jet != (*thePuppiAK8Jets).end(); jet++ ) {
      if(Skim_=="MCJECs" &&  (&*jet) ->genJet() ==0) continue;
      else if(Skim_=="MCJECs" &&  (&*jet) ->genJet() ->pt()<10) continue;
      else if((Skim_=="ZJetsResiduals" || Skim_=="GammaJetsResiduals") &&  (&*jet)->correctedP4("Uncorrected").Pt() <10)  continue;
      else if((&*jet)->pt()<AK8JetPtCut_) continue;

      //      if(fabs((&*jet)->eta())>2.4) continue;
      _puppiak8jetEta.push_back((&*jet)->eta());
      _puppiak8jetPhi.push_back((&*jet)->phi());
      _puppiak8jetPt.push_back((&*jet)->pt());
      _puppiak8jetRawPt.push_back((&*jet)->correctedP4("Uncorrected").Pt());
      if( (&*jet) ->genJet()  !=0) _puppiak8jetPtGen.push_back(  (&*jet) ->genJet()->pt() );
      else  _puppiak8jetPtGen.push_back(  -1);
      _puppiak8jet_tau1.push_back((&*jet)->userFloat("NjettinessAK8Puppi:tau1"));
      _puppiak8jet_tau2.push_back((&*jet)->userFloat("NjettinessAK8Puppi:tau2"));
      _puppiak8jet_tau3.push_back((&*jet)->userFloat("NjettinessAK8Puppi:tau3"));
    }
  }


  //AK8  jets
  edm::Handle< std::vector< pat::Jet> > theAK8Jets;
  iEvent.getByToken(jetAK8Token_,theAK8Jets );
  if(theAK8Jets.isValid()){
    for( std::vector<pat::Jet>::const_iterator jet = (*theAK8Jets).begin(); jet != (*theAK8Jets).end(); jet++ ) {
      if(Skim_=="MCJECs" &&  (&*jet) ->genJet() ==0) continue;
      else if(Skim_=="MCJECs" &&  (&*jet) ->genJet() ->pt()<10) continue;
      else if((Skim_=="ZJetsResiduals" || Skim_=="GammaJetsResiduals") &&  (&*jet)->correctedP4("Uncorrected").Pt() <10)  continue;
      else if((&*jet)->pt()<AK8JetPtCut_) continue;

      //      if(fabs((&*jet)->eta())>2.4) continue;
      _ak8jetEta.push_back((&*jet)->eta());
      _ak8jetPhi.push_back((&*jet)->phi());
      _ak8jetPt.push_back((&*jet)->pt());
      _ak8jetArea.push_back((&*jet)->jetArea());
      _ak8jetRawPt.push_back((&*jet)->correctedP4("Uncorrected").Pt());
      if( (&*jet) ->genJet()  !=0) _ak8jetPtGen.push_back(  (&*jet) ->genJet()->pt() );
      else  _ak8jetPtGen.push_back(  -1);
    }
  }



  //Calo jets
  edm::Handle< std::vector< reco::CaloJet> > theCaloJets;
  iEvent.getByToken(jetCaloToken_,theCaloJets );
  edm::Handle<edm::Association<reco::GenJetCollection>> genJetMatchCalo;
  iEvent.getByToken(genJetAssocCaloToken_, genJetMatchCalo);

  if(theCaloJets.isValid()){
    for( std::vector<reco::CaloJet>::const_iterator jet = (*theCaloJets).begin(); jet != (*theCaloJets).end(); jet++ ) {
          
      edm::RefToBase<reco::CaloJet> jetRef(edm::Ref<reco::CaloJetCollection>( theCaloJets , jet -theCaloJets->begin()));
      const reco::GenJet * updatedgenjet =0;
      if(genJetMatchCalo.isValid() ) updatedgenjet = ( (*genJetMatchCalo)[jetRef].isNonnull() && (*genJetMatchCalo)[jetRef].isAvailable()) ? &*(*genJetMatchCalo)[jetRef] : 0 ;
      Float_t jetptgen(-99.);
      if( updatedgenjet !=0 ) jetptgen= updatedgenjet->pt() ;        
      if( updatedgenjet ==0   && DropUnmatchedJets_ && (&*jet)->pt()<50 ) continue;
      if(!PassJetPreselection((&*jet) , jetptgen ) ) continue;

      bool isleptoncleaned =IsLeptonPhotonCleaned ((&*jet) );
      _calojetLeptonPhotonCleaned.push_back(isleptoncleaned);
      _calojetPtGen.push_back(jetptgen );
      _calojetEta.push_back((&*jet)->eta());
      _calojetPhi.push_back((&*jet)->phi());
      _calojetPt.push_back((&*jet)->pt());
      _calojetRawPt.push_back(0.);

    }
  }

  //AK4 PF jets (no CHS)
  edm::Handle< std::vector< pat::Jet> > thenoCHSJets;
  iEvent.getByToken(jetnoCHSToken_,thenoCHSJets );
  if(thenoCHSJets.isValid()){
    for( std::vector<pat::Jet>::const_iterator jet = (*thenoCHSJets).begin(); jet != (*thenoCHSJets).end(); jet++ ) {
      const reco::GenJet * genjet = (&*jet) ->genJet()  ;
      if( genjet ==0   && DropUnmatchedJets_ && (&*jet)->pt()<50 ) continue;
      Float_t jetptgen(-99.);
      if( genjet !=0 ) jetptgen= genjet->pt() ;
      if(!PassJetPreselection((&*jet) , jetptgen , false) ) continue;
      
      bool isleptoncleaned =IsLeptonPhotonCleaned ((&*jet) );
      _noCHSjetLeptonPhotonCleaned.push_back(isleptoncleaned);
      _noCHSjetEta.push_back((&*jet)->eta());
      _noCHSjetPhi.push_back((&*jet)->phi());
      _noCHSjetPt.push_back((&*jet)->pt());
      _noCHSjetRawPt.push_back((&*jet)-> correctedP4("Uncorrected").Pt()   );
      _noCHSjetPtGen.push_back(jetptgen);
    }
  }

  
  //AK4 gen jets
  //Mostly needed for MC jecs: it can happen that a low pt gen jet is not matched to any reco jet. Not considering those would introduce a bias 
  edm::Handle< std::vector< reco::GenJet> > theGenJets;
  iEvent.getByToken(genjetToken_,theGenJets );
  if(theGenJets.isValid()){
    for( std::vector<reco::GenJet>::const_iterator jet = (*theGenJets).begin(); jet != (*theGenJets).end(); jet++ ) {
      if((&*jet)->pt()<5&& Skim_=="MCJECs") continue;
      _genjetEta.push_back((&*jet)->eta());
      _genjetPt.push_back((&*jet)->pt());
      _genjetPhi.push_back((&*jet)->phi());
      _genjet_CHSIdx.push_back( GetRecoIdx((&*jet),_jetPt, _jetEta,_jetPhi, _jetPtGen));
      _genjet_noCHSIdx.push_back( GetRecoIdx((&*jet),_noCHSjetPt, _noCHSjetEta,_noCHSjetPhi, _noCHSjetPtGen));
      _genjet_PuppiIdx.push_back( GetRecoIdx((&*jet),_puppijetPt, _puppijetEta,_puppijetPhi, _puppijetPtGen));
      _genjet_CaloIdx.push_back( GetRecoIdx((&*jet),_calojetPt, _calojetEta,_calojetPhi, _calojetPtGen));
      
      }
  }
  //AK8 gen jets
  edm::Handle< std::vector< reco::GenJet> > theGenAK8Jets;
  iEvent.getByToken(genAK8jetToken_,theGenAK8Jets );
  if(theGenAK8Jets.isValid()){
    for( std::vector<reco::GenJet>::const_iterator jet = (*theGenAK8Jets).begin(); jet != (*theGenAK8Jets).end(); jet++ ) {
      _genAK8jetEta.push_back((&*jet)->eta());
      _genAK8jetPt.push_back((&*jet)->phi());
      _genAK8jetPhi.push_back((&*jet)->pt());
      _genAK8jet_PuppiIdx.push_back( GetRecoIdx((&*jet),_puppiak8jetPt, _puppiak8jetEta,_puppiak8jetPhi, _puppiak8jetPtGen));
      _genAK8jet_CHSIdx.push_back( GetRecoIdx((&*jet),_ak8jetPt, _ak8jetEta,_ak8jetPhi, _ak8jetPtGen));

    }
  }




  //Type 1 PFMET
  edm::Handle< vector<pat::MET> > ThePFMET;
  iEvent.getByToken(metToken_, ThePFMET);
  const vector<pat::MET> *pfmetcol = ThePFMET.product();
  const pat::MET *pfmet;
  pfmet = &(pfmetcol->front());
  _met = pfmet->pt();
  _met_phi = pfmet->phi();
  _rawmet =  pfmet->uncorPt();
  _rawmet_phi =  pfmet->uncorPhi();

  _chsmet = pfmet->corPt(pat::MET::METCorrectionLevel::Type01);
  _chsmet_phi = pfmet->corPhi(pat::MET::METCorrectionLevel::Type01);
  _rawchsmet = pfmet->corPt(pat::MET::METCorrectionLevel::RawChs);
  _rawchsmet_phi = pfmet->corPhi(pat::MET::METCorrectionLevel::RawChs);



  //PUPPI MET
  edm::Handle< vector<pat::MET> > ThePUPPIMET;
  iEvent.getByToken(puppimetToken_, ThePUPPIMET);
  const vector<pat::MET> *puppimetcol = ThePUPPIMET.product();
  const pat::MET *puppimet;
  puppimet = &(puppimetcol->front());
  _puppimet = puppimet->pt();
  _puppimet_phi = puppimet->phi();
  _puppirawmet = puppimet->uncorPt();
  _puppirawmet_phi = puppimet->uncorPhi();
  
  
  _n_PFele_fromvtxfit=0;
  _n_PFmu_fromvtxfit=0;

  //PF candidates
  edm::Handle<std::vector< pat::PackedCandidate>> pfcands;
  iEvent.getByToken(pfcandsToken_ ,pfcands);
  for(int i = 0; i< 6; i++){
    _n_CH_fromvtxfit[i] = 0;
    _HT_CH_fromvtxfit[i] = 0;
  }
  //Accessing the puppi weights !
  edm::Handle<edm::ValueMap<float> > puppiweights;
  iEvent.getByToken(puppiweightsToken_, puppiweights);
  /*
    To map vertex idx with vertex ref
    https://github.com/cms-sw/cmssw/blob/master/DataFormats/PatCandidates/interface/PackedCandidate.h#L705-L706
    I need to check that the vertex idx used in fromPV is equivalent to the vertex position in the vertex collection
    
  */
  for( std::vector<pat::PackedCandidate>::const_iterator p = (*pfcands).begin(); p != (*pfcands).end(); p++ ) {
    int idxrefvtx= (p->vertexRef().isNonnull()) ?  p->vertexRef().key() : -1;
    
    edm::RefToBase<pat::PackedCandidate> pfcandRef (edm::Ref<pat::PackedCandidateCollection>( pfcands , p  - pfcands->begin()));
    /*
    cout << "PFCH cand pt, eta, phi, pdgid: " << p->pt() <<", "<<p->eta()<<", " << p->phi() <<", "<<p->pdgId() <<endl;
    cout << "PFCH idxrefvtx, fromPV(idxrefvtx), fromPV(0):  " << idxrefvtx<<", "<< p->fromPV(idxrefvtx)<< ", "<<p->fromPV()<<endl; 
    cout << "the puppi weight is " <<   (*puppiweights)[pfcandRef]  <<endl; 
    */
    if( idxrefvtx >=0 && p->fromPV(idxrefvtx) >=2 && fabs(p->pdgId())  == 211 ){
      TVector2 ptandphi;
      ptandphi.SetMagPhi(p->pt(),p->phi());
      metPV[idxrefvtx] -= ptandphi;
      _SumPT2CH_PV[idxrefvtx] += p->pt()*p->pt();
    }
    
    if(fabs(p->pdgId())  == 11&& p->fromPV(0)==3 &&p->pt()>10) _n_PFele_fromvtxfit ++; 
    if(fabs(p->pdgId())  == 13&& p->fromPV(0)==3 &&p->pt()>10) _n_PFmu_fromvtxfit ++;
    if(fabs(p->pdgId())  == 211&& p->fromPV(0)==3 ){
      for(int i = 0; i< 6; i++){
	if(p->pt()>_PtCutPFforMultiplicity[i]){
	  _n_CH_fromvtxfit[i] ++;
	  _HT_CH_fromvtxfit[i] += p->pt();
	}
      }
    }
        
    if(p->pt()<PFCandPtCut_)continue;
    _PFcand_pt.push_back(p->pt());
    _PFcand_eta.push_back(p->eta());
    _PFcand_phi.push_back(p->phi());
    _PFcand_pdgId.push_back(p->pdgId());
    _PFcand_fromPV.push_back(p->fromPV(0));//See https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Packed_ParticleFlow_Candidates
    _PFcand_dz.push_back(p->dz(0));
    if(p-> hasTrackDetails())  _PFcand_dzError.push_back(p->dzError());
    else _PFcand_dzError.push_back(-99);

    _PFcand_hcalFraction.push_back(p->hcalFraction());
    _PFcand_PVfitidx.push_back(idxrefvtx);
    // _PFcand_puppiweight.push_back((*puppiweights)[pfcandRef]);
    //    cout<< "reaching that point " <<endl;
    if(puppiweights.isValid())_PFcand_puppiweight.push_back((*puppiweights)[pfcandRef]);
    else _PFcand_puppiweight.push_back(0.);
  }
  //Now computing the met from CH for each PV
  for(unsigned int i = 0;i<theVertices->size() ; i++) {_METCH_PV[i] = metPV[i].Mod(); _METPhiCH_PV[i] = metPV[i].Phi(); }


  //Gen particle info
  edm::Handle<GenParticleCollection> TheGenParticles;
  iEvent.getByToken(genpartToken_, TheGenParticles);
  TLorentzVector Gen0;
  Gen0.SetPtEtaPhiE( 0, 0, 0, 0);
  _genHT = 0;
  if(TheGenParticles.isValid()){
    for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
      int id = TMath::Abs(p->pdgId());
      /*      if ((id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id== 6 || id == 21 || id == 22 ) ) cout << "p ID, status, pt, eta, phi " << 
	      p->pdgId() <<", "<<
	      p->status() <<", "<<
	      p->pt() <<", "<<
	      p->eta() <<", "<<
	      p->phi() <<
	      endl;
      */
      if ((id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 21 || id == 22 ) && (p->status() == 23&& TMath::Abs(p->mother()->pdgId())!=6 &&TMath::Abs(p->mother()->pdgId())!=24 )){
	_genHT += p->pt();
      }
      if ( (id == 12 || id == 14 || id == 16 ) && (p->status() == 1) ) {
	TLorentzVector Gen;
	Gen.SetPtEtaPhiE( p->pt(), p->eta(), p->phi(), p->energy() );
	Gen0 += Gen;
      }
      if( id ==11 && p->status() == 1 &&  (p->pt() > 0.8* ElectronPtCut_ || p->pt()>50 ) ){ 
	_lgenPt.push_back(p->pt());
	_lgenEta.push_back(p->eta());
	_lgenPhi.push_back(p->phi());
	_lgenpdgId.push_back(p->pdgId());
      }
      if( id==13 && p->status() == 1 &&  (p->pt() > 0.8* MuonPtCut_ || p->pt()>50 ) ){ 
	_lgenPt.push_back(p->pt());
	_lgenEta.push_back(p->eta());
	_lgenPhi.push_back(p->phi());
	_lgenpdgId.push_back(p->pdgId());
      }
 
      if( (id ==22) && p->status() == 1  &&  (p->pt() > 0.8* PhotonPtCut_ || p->pt()>50 ) ){ 
	_phgenPt.push_back(p->pt());
	_phgenEta.push_back(p->eta());
	_phgenPhi.push_back(p->phi());
	
      }
    }


    //Compute dilepton variables
    Float_t mll_gen(0),ptll_gen(0),pzll_gen(0),yll_gen(0),phill_gen(0),dphill_gen(0),costhCSll_gen(0);
    if(_lgenPt.size() >=2){
      
      for(unsigned int i =0; i < _lgenPt.size(); i++){
	if (fabs(_lgenpdgId[i]) !=11 && fabs(_lgenpdgId[i])!=13 ) continue;
	for(unsigned int j =0; j < i; j++){
	  if(fabs(_lgenpdgId[j]) !=11 && fabs(_lgenpdgId[j])!=13 ) continue;
	  if( _lgenpdgId[i] * _lgenpdgId[j] >0 ) continue;
	  CalcDileptonInfoGen(i,j, mll_gen,ptll_gen,pzll_gen,yll_gen,phill_gen,dphill_gen,costhCSll_gen);
	  _mll_gen= mll_gen; _ptll_gen=ptll_gen; _pzll_gen=pzll_gen; _yll_gen=yll_gen;  _phill_gen = phill_gen; _dphill_gen=dphill_gen; _costhCSll_gen=costhCSll_gen;
	  if(fabs(_lgenpdgId[j]) ==11 && fabs(_lgenpdgId[j])  ==11) _ngenElesll =2 ;
	  else if(fabs(_lgenpdgId[j]) ==13 && fabs(_lgenpdgId[j])  ==13) _ngenElesll =0 ;
	  else _ngenElesll = 1;
	  
	}
      }
    }
    

  }
  if (Gen0.E()!=0) {
    _genmet = Gen0.Pt();
    _genmet_phi = Gen0.Phi();
  } else {
    _genmet = 0;
    _genmet_phi = 0;
  }

  
  edm::Handle<LHEEventProduct> lhe_handle;
  iEvent.getByToken(lheEventToken_, lhe_handle);
  if ( !lhe_handle.isValid()) iEvent.getByToken(lheEventALTToken_, lhe_handle);
  double lheht = 0;
  if (lhe_handle.isValid()){ 
    std::vector<lhef::HEPEUP::FiveVector> lheParticles = lhe_handle->hepeup().PUP;
    ROOT::Math::PxPyPzEVector cand_;

    for (unsigned i = 0; i < lheParticles.size(); ++i) {
      cand_ = ROOT::Math::PxPyPzEVector(lheParticles[i][0],lheParticles[i][1],lheParticles[i][2],lheParticles[i][3]);
    
      int imotherfirst = lhe_handle->hepeup().MOTHUP[i].first-1;// See the warning here: https://github.com/cms-sw/cmssw/blob/master/GeneratorInterface/AlpgenInterface/src/AlpgenEventRecordFixes.cc#L4-L16
      int imotherlast = lhe_handle->hepeup().MOTHUP[i].second-1;
      bool istopdecay = fabs( lhe_handle->hepeup().IDUP[imotherfirst]) ==6 ||fabs( lhe_handle->hepeup().IDUP[imotherlast])==6 ;
      bool isZdecay = fabs( lhe_handle->hepeup().IDUP[imotherfirst]) ==23 ||fabs( lhe_handle->hepeup().IDUP[imotherlast])==23 ;
      bool isWdecay = fabs( lhe_handle->hepeup().IDUP[imotherfirst]) ==24 ||fabs( lhe_handle->hepeup().IDUP[imotherlast])==24 ;
      //      cout <<"LHE parti: pdgid, pt: "<< lhe_handle->hepeup().IDUP[i] << ", "<< cand_.Pt()<< ", "<< lhe_handle->hepeup().IDUP[imotherfirst] <<", "<<lhe_handle->hepeup().IDUP[imotherlast] <<endl;
      //cout <<"imotherfirst/last " <<imotherfirst <<", "<<imotherlast<<endl;
      if(istopdecay) continue;
      if(isZdecay) continue;
      if(isWdecay) continue;
      if(fabs( lhe_handle->hepeup().IDUP[i]) <=5|| lhe_handle->hepeup().IDUP[i] ==21 ) lheht += cand_.Pt();//That definition works to retrieve the HT in madgraph HT binned QCD 
    }
  }
  _genHT = lheht;
  //Tested with QCD/photon jets/DY with madgraphm


  //Filling trees and histos   
  if(PassSkim()){
      if(SaveTree_)outputTree->Fill();
      h_PFMet->Fill(_met);
      h_PuppiMet->Fill(_puppimet);
      h_nvtx->Fill(_n_PV);
    }
}


// ------------ method called once each job just before starting event loop  ------------
void
JMEAnalyzer::beginJob()
{



  outputTree->Branch("_tauEta",&_tauEta);
  outputTree->Branch("_tauPhi",&_tauPhi);
  outputTree->Branch("_tauPt",&_tauPt);
  outputTree->Branch("_tauPassMediumID",&_tauPassMediumID);

 
  if(Skim_.find("L1Study") !=std::string::npos){


    outputTree->Branch("_eventNb",   &_eventNb,   "_eventNb/l");
    outputTree->Branch("_runNb",     &_runNb,     "_runNb/l");
    outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");
    outputTree->Branch("_bx", &_bx, "_bx/l");
    outputTree->Branch("_n_PV", &_n_PV, "_n_PV/I");
    outputTree->Branch("_LV_x", &_LV_x, "_LV_x/F");
    outputTree->Branch("_LV_y", &_LV_y, "_LV_y/F");
    outputTree->Branch("_LV_z", &_LV_z, "_LV_z/F");

    outputTree->Branch("_l1prefire",&_l1prefire,"_l1prefire/O");



    outputTree->Branch("_lEta",&_lEta);
    outputTree->Branch("_lPhi",&_lPhi);
    outputTree->Branch("_lPt",&_lPt);
    outputTree->Branch("_lPtcorr",&_lPtcorr);
    outputTree->Branch("_lPtSC",&_lPtSC);
    outputTree->Branch("_lpdgId",&_lpdgId);
    outputTree->Branch("_lPassTightID",&_lPassTightID);
    outputTree->Branch("_lPassLooseID",&_lPassLooseID);
    outputTree->Branch("_lisSAMuon",&_lisSAMuon);
    outputTree->Branch("_nEles", &_nEles, "_nEles/I");
    outputTree->Branch("_nMus", &_nMus, "_nMus/I");
    outputTree->Branch("_ldz",&_ldz);
    outputTree->Branch("_ldzError",&_ldzError);
    outputTree->Branch("_ldxy",&_ldxy);
    outputTree->Branch("_ldxyError",&_ldxyError);
    
    outputTree->Branch("_l3dIP",&_l3dIP);
    outputTree->Branch("_l3dIPError",&_l3dIPError);
    outputTree->Branch("_lpassHLT_IsoMu24",&_lpassHLT_IsoMu24);
    outputTree->Branch("_lpassHLT_Ele32_WPTight_Gsf",&_lpassHLT_Ele32_WPTight_Gsf);
    

    outputTree->Branch("_phEta",&_phEta);
    outputTree->Branch("_phPhi",&_phPhi);
    outputTree->Branch("_phPt",&_phPt);
    outputTree->Branch("_phPtcorr",&_phPtcorr);
    outputTree->Branch("_phPassTightID",&_phPassTightID);
    outputTree->Branch("_phPassIso",&_phPassIso);

    //outputTree->Branch("passL1_Initial_bxmin1",&passL1_Initial_bxmin1,"passL1_Initial_bxmin1[512]/O");
    outputTree->Branch("passL1_Initial_bx0",&passL1_Initial_bx0,"passL1_Initial_bx0[512]/O");
    //outputTree->Branch("passL1_Initial_bxplus1",&passL1_Initial_bxplus1,"passL1_Initial_bxplus1[512]/O");
    

    outputTree->Branch("_L1mu_Qual",&_L1mu_Qual);
    outputTree->Branch("_L1mu_pt",&_L1mu_pt);
    outputTree->Branch("_L1mu_eta",&_L1mu_eta);
    outputTree->Branch("_L1mu_phi",&_L1mu_phi);
    outputTree->Branch("_L1mu_bx",&_L1mu_bx);
    outputTree->Branch("_L1mu_TFIdx",&_L1mu_TFIdx);
    outputTree->Branch("_L1mu_upt",&_L1mu_upt);
    outputTree->Branch("_L1mu_charge",&_L1mu_charge);
    outputTree->Branch("_L1mu_dXY",&_L1mu_dXY);

    
    outputTree->Branch("_ptgamma", &_ptgamma, "_ptgamma/F");
    outputTree->Branch("_phigamma", &_phigamma, "_phigamma/F");
    outputTree->Branch("_mll", &_mll, "_mll/F");
    outputTree->Branch("_ptll", &_ptll, "_ptll/F");
    outputTree->Branch("_phill", &_phill, "_phill/F");

    if(Skim_.find("ZTo") !=std::string::npos){

    outputTree->Branch("_pzll", &_pzll, "_pzll/F");
    outputTree->Branch("_yll", &_yll, "_yll/F");
    outputTree->Branch("_dphill", &_dphill, "_dphill/F");
    outputTree->Branch("_costhCSll", &_costhCSll, "_costhCSll/F");
    outputTree->Branch("_nElesll",&_nElesll,"_nElesll/I");
    }

    
    if(Skim_.find("ZToEE") !=std::string::npos){
    outputTree->Branch("_L1eg_pt",&_L1eg_pt);
    outputTree->Branch("_L1eg_eta",&_L1eg_eta);
    outputTree->Branch("_L1eg_phi",&_L1eg_phi);
    outputTree->Branch("_L1eg_bx",&_L1eg_bx);
    outputTree->Branch("_L1eg_iso",&_L1eg_iso);
    }


    if(Skim_.find("JME") !=std::string::npos){
    
    outputTree->Branch("_L1jet_pt",&_L1jet_pt);
    outputTree->Branch("_L1jet_eta",&_L1jet_eta);
    outputTree->Branch("_L1jet_phi",&_L1jet_phi);
    outputTree->Branch("_L1jet_bx",&_L1jet_bx);
    }


    // outputTree->Branch("_PUV1_x", &_PUV1_x, "_PUV1_x/F");
    //outputTree->Branch("_PUV1_y", &_PUV1_y, "_PUV1_y/F");
    // outputTree->Branch("_PUV1_z", &_PUV1_z, "_PUV1_z/F");
    //    outputTree->Branch("_n_CH_fromvtxfit",&_n_CH_fromvtxfit,"_n_CH_fromvtxfit[6]/I");
    //outputTree->Branch("_HT_CH_fromvtxfit", &_HT_CH_fromvtxfit, "_HT_CH_fromvtxfit[6]/F");
    
    outputTree->Branch("_met", &_met, "_met/F");
    outputTree->Branch("_met_phi", &_met_phi, "_met_phi/F");
    outputTree->Branch("_puppimet", &_puppimet, "_puppimet/F");
    outputTree->Branch("_puppimet_phi", &_puppimet_phi, "_puppimet_phi/F");
    
    outputTree->Branch("_rawmet", &_rawmet, "_rawmet/F");
    outputTree->Branch("_rawmet_phi", &_rawmet_phi, "_rawmet_phi/F");
    outputTree->Branch("_puppirawmet", &_puppirawmet, "_puppirawmet/F");
    outputTree->Branch("_puppirawmet_phi", &_puppirawmet_phi, "_puppirawmet_phi/F");
    
    //outputTree->Branch("_rawchsmet",&_rawchsmet,"_rawchsmet/F");
    //outputTree->Branch("_rawchsmet_phi",&_rawchsmet_phi,"_rawchsmet_phi/F");
    //outputTree->Branch("_chsmet",&_chsmet,"_chsmet/F");
    //outputTree->Branch("_chsmet_phi",&_chsmet_phi,"_chsmet_phi/F");
    
    outputTree->Branch("_jetEta",&_jetEta);
    outputTree->Branch("_jetPhi",&_jetPhi);
    outputTree->Branch("_jetPt",&_jetPt);
    //    outputTree->Branch("_jetRawPt",&_jetRawPt);
    outputTree->Branch("_jet_CHEF",&_jet_CHEF);
    outputTree->Branch("_jet_NHEF",&_jet_NHEF);
    outputTree->Branch("_jet_NEEF",&_jet_NEEF);
    outputTree->Branch("_jet_CEEF",&_jet_CEEF);
    outputTree->Branch("_jet_MUEF",&_jet_MUEF);
    /*outputTree->Branch("_jet_CHM",&_jet_CHM);
    outputTree->Branch("_jet_NHM",&_jet_NHM);
    outputTree->Branch("_jet_PHM",&_jet_PHM);
    outputTree->Branch("_jet_NM",&_jet_NM);*/
    outputTree->Branch("_jetPassID",&_jetPassID);
    outputTree->Branch("_jetLeptonPhotonCleaned",&_jetLeptonPhotonCleaned);
    outputTree->Branch("_jetTauCleaned",&_jetTauCleaned);
    
    outputTree->Branch("_jetDeepJet_b",&_jetDeepJet_b);
    /*outputTree->Branch("_jetParticleNet_b",&_jetParticleNet_b);
    outputTree->Branch("_jetDeepJet_c",&_jetDeepJet_c);
    outputTree->Branch("_jetDeepJet_uds",&_jetDeepJet_uds);
    outputTree->Branch("_jetDeepJet_g",&_jetDeepJet_g);
    outputTree->Branch("_jetQuarkGluonLikelihood",&_jetQuarkGluonLikelihood);*/
    outputTree->Branch("_jetPUMVAUpdate",&_jetPUMVAUpdate);
    outputTree->Branch("_jethfsigmaEtaEta",&_jethfsigmaEtaEta);
    outputTree->Branch("_jethfsigmaPhiPhi",&_jethfsigmaPhiPhi);
    outputTree->Branch("_jethfcentralEtaStripSize",&_jethfcentralEtaStripSize);
    outputTree->Branch("_jethfadjacentEtaStripsSize",&_jethfadjacentEtaStripsSize);

    

    outputTree->Branch("Flag_IsUnprefirable",&Flag_IsUnprefirable,"Flag_IsUnprefirable/O");
    outputTree->Branch("Flag_goodVertices",&Flag_goodVertices,"Flag_goodVertices/O");
    outputTree->Branch("Flag_globalTightHalo2016Filter",&Flag_globalTightHalo2016Filter,"Flag_globalTightHalo2016Filter/O");
    outputTree->Branch("Flag_globalSuperTightHalo2016Filter",&Flag_globalSuperTightHalo2016Filter,"Flag_globalSuperTightHalo2016Filter/O");
    outputTree->Branch("Flag_HBHENoiseFilter",&Flag_HBHENoiseFilter,"Flag_HBHENoiseFilter/O");
    outputTree->Branch("Flag_HBHENoiseIsoFilter",&Flag_HBHENoiseIsoFilter,"Flag_HBHENoiseIsoFilter/O");
    outputTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter",&Flag_EcalDeadCellTriggerPrimitiveFilter,"Flag_EcalDeadCellTriggerPrimitiveFilter/O");
    outputTree->Branch("Flag_BadPFMuonFilter",&Flag_BadPFMuonFilter,"Flag_BadPFMuonFilter/O");
    outputTree->Branch("Flag_BadPFMuonDzFilter",&Flag_BadPFMuonDzFilter,"Flag_BadPFMuonDzFilter/O");
    outputTree->Branch("Flag_hfNoisyHitsFilter",&Flag_hfNoisyHitsFilter,"Flag_hfNoisyHitsFilter/O");
    outputTree->Branch("Flag_BadChargedCandidateFilter",&Flag_BadChargedCandidateFilter,"Flag_BadChargedCandidateFilter/O");
    outputTree->Branch("Flag_eeBadScFilter",&Flag_eeBadScFilter,"Flag_eeBadScFilter/O");
    outputTree->Branch("Flag_ecalBadCalibFilter",&Flag_ecalBadCalibFilter,"Flag_ecalBadCalibFilter/O");
    outputTree->Branch("Flag_ecalLaserCorrFilter",&Flag_ecalLaserCorrFilter,"Flag_ecalLaserCorrFilter/O");
    outputTree->Branch("Flag_EcalDeadCellBoundaryEnergyFilter",&Flag_EcalDeadCellBoundaryEnergyFilter,"Flag_EcalDeadCellBoundaryEnergyFilter/O");


    outputTree->Branch("HLT_PixelClusters_WP1",&HLT_PixelClusters_WP1,"HLT_PixelClusters_WP1/O");
    outputTree->Branch("HLT_PixelClusters_WP2",&HLT_PixelClusters_WP2,"HLT_PixelClusters_WP2/O");
    outputTree->Branch("HLT_PixelClusters_WP2_split",&HLT_PixelClusters_WP2_split,"HLT_PixelClusters_WP2_split/O");
    outputTree->Branch("HLT_IsoMu24",&HLT_IsoMu24,"HLT_IsoMu24/O");
    outputTree->Branch("HLT_Ele32_WPTight_Gsf",&HLT_Ele32_WPTight_Gsf,"HLT_Ele32_WPTight_Gsf/O");


    outputTree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60/O");
    outputTree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight,"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight/O");
    outputTree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF",&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF,"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF/O");
    outputTree->Branch("HLT_PFHT1050",&HLT_PFHT1050,"HLT_PFHT1050/O");
    outputTree->Branch("HLT_DiJet110_35_Mjj650_PFMET110",&HLT_DiJet110_35_Mjj650_PFMET110,"HLT_DiJet110_35_Mjj650_PFMET110/O");
    outputTree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8/O");
    outputTree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL/O");
    outputTree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ/O");
    outputTree->Branch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL/O");
    outputTree->Branch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ/O");

    outputTree->Branch("HLT_Photon110EB_TightID_TightIso",&HLT_Photon110EB_TightID_TightIso,"HLT_Photon110EB_TightID_TightIso/O");

    outputTree->Branch("HLT_Photon50_R9Id90_HE10_IsoM",&HLT_Photon50_R9Id90_HE10_IsoM,"HLT_Photon50_R9Id90_HE10_IsoM/O");

    outputTree->Branch("hltSinglePFJet60",&hltSinglePFJet60);
    outputTree->Branch("hltSinglePFJet80",&hltSinglePFJet80);
    outputTree->Branch("hltSinglePFJet140",&hltSinglePFJet140);
    outputTree->Branch("hltSinglePFJet200",&hltSinglePFJet200);
    outputTree->Branch("hltSinglePFJet260",&hltSinglePFJet260);
    outputTree->Branch("hltSinglePFJet320",&hltSinglePFJet320);
    outputTree->Branch("hltSinglePFJet400",&hltSinglePFJet400);
    outputTree->Branch("hltSinglePFJet450",&hltSinglePFJet450);
    outputTree->Branch("hltSinglePFJet500",&hltSinglePFJet500);

    outputTree->Branch("hltEG200HEFilter",&hltEG200HEFilter );
    outputTree->Branch("hltEG110EBTightIDTightIsoTrackIsoFilter",&hltEG110EBTightIDTightIsoTrackIsoFilter);
    

    return;
  }
 

  if(Skim_=="ZJetsResiduals" || Skim_=="GammaJetsResiduals"  ){
    outputTree->Branch("_jetEta",&_jetEta);
    outputTree->Branch("_jetPhi",&_jetPhi);
    outputTree->Branch("_jetPt",&_jetPt);
    outputTree->Branch("_jetRawPt",&_jetRawPt);
    outputTree->Branch("_jetArea",&_jetArea);
    outputTree->Branch("_jetPassID",&_jetPassID);
    outputTree->Branch("_jetLeptonPhotonCleaned",&_jetLeptonPhotonCleaned);
    outputTree->Branch("_jetTauCleaned",&_jetTauCleaned);
    outputTree->Branch("_jetPtNoL2L3Res",&_jetPtNoL2L3Res);
    outputTree->Branch("_jetPUMVAUpdate2017",&_jetPUMVAUpdate2017);
    outputTree->Branch("_jetPUMVAUpdate2018",&_jetPUMVAUpdate2018);
    outputTree->Branch("_puppijetEta",&_puppijetEta);
    outputTree->Branch("_puppijetPhi",&_puppijetPhi);
    outputTree->Branch("_puppijetPt",&_puppijetPt);
    outputTree->Branch("_puppijetRawPt",&_puppijetRawPt);
    outputTree->Branch("_puppijetLeptonPhotonCleaned",&_puppijetLeptonPhotonCleaned);
    outputTree->Branch("_puppijetPassID",&_puppijetPassID);
    outputTree->Branch("_puppijetPtNoL2L3Res",&_puppijetPtNoL2L3Res);


    if(IsMC_){
      outputTree->Branch("_jetPtGen",&_jetPtGen);
      outputTree->Branch("_jetEtaGen",&_jetEtaGen);
      outputTree->Branch("_jetPhiGen",&_jetPhiGen);
      outputTree->Branch("_jetPtGenWithNu",&_jetPtGenWithNu);
      outputTree->Branch("_puppijetPtGen",&_puppijetPtGen);
      outputTree->Branch("_puppijetPtGenWithNu",&_puppijetPtGenWithNu);
      outputTree->Branch("trueNVtx", &trueNVtx,"trueNVtx/I");
      outputTree->Branch("_genmet", &_genmet, "_genmet/F");
      outputTree->Branch("_genmet_phi", &_genmet_phi, "_genmet_phi/F");
      outputTree->Branch("_jetJECuncty",&_jetJECuncty); 
      outputTree->Branch("_puppijetJECuncty",&_puppijetJECuncty);
    }
    outputTree->Branch("_met", &_met, "_met/F");
    outputTree->Branch("_met_phi", &_met_phi, "_met_phi/F");
    outputTree->Branch("_puppimet", &_puppimet, "_puppimet/F");
    outputTree->Branch("_puppimet_phi", &_puppimet_phi, "_puppimet_phi/F");
    outputTree->Branch("_rawmet", &_rawmet, "_rawmet/F");
    outputTree->Branch("_rawmet_phi", &_rawmet_phi, "_rawmet_phi/F");
    outputTree->Branch("_puppirawmet", &_puppirawmet, "_puppirawmet/F");
    outputTree->Branch("_puppirawmet_phi", &_puppirawmet_phi, "_puppirawmet_phi/F");
    outputTree->Branch("_rawchsmet",&_rawchsmet,"_rawchsmet/F");
    outputTree->Branch("_rawchsmet_phi",&_rawchsmet_phi,"_rawchsmet_phi/F");
    outputTree->Branch("_chsmet",&_chsmet,"_chsmet/F");
    outputTree->Branch("_chsmet_phi",&_chsmet_phi,"_chsmet_phi/F");
    outputTree->Branch("_n_PV", &_n_PV, "_n_PV/I");
    outputTree->Branch("_rho", &_rho, "_rho/F");
    outputTree->Branch("_eventNb",   &_eventNb,   "_eventNb/l");
    outputTree->Branch("_runNb",     &_runNb,     "_runNb/l");
    outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");

  }


  if(Skim_=="ZJetsResiduals"  ){
    outputTree->Branch("_ptll", &_ptll, "_ptll/F");
    outputTree->Branch("_phill", &_phill, "_phill/F");
    outputTree->Branch("_lEta",&_lEta);
    outputTree->Branch("_lPhi",&_lPhi);
    outputTree->Branch("_lPt",&_lPt);
    outputTree->Branch("_lPtcorr",&_lPtcorr);
    outputTree->Branch("_lpdgId",&_lpdgId);
    outputTree->Branch("_lPassTightID",&_lPassTightID);
    outputTree->Branch("_lPassLooseID",&_lPassLooseID);
    outputTree->Branch("_nEles", &_nEles, "_nEles/I");
    outputTree->Branch("_nMus", &_nMus, "_nMus/I");

    if(IsMC_){
    outputTree->Branch("_ptll_gen", &_ptll_gen, "_ptll_gen/F");
    outputTree->Branch("_phill_gen", &_phill_gen, "_phill_gen/F");
    outputTree->Branch("_lgenPt",&_lgenPt);
    outputTree->Branch("_lgenEta",&_lgenEta);
    outputTree->Branch("_lgenPhi",&_lgenPhi);
    outputTree->Branch("_lgenpdgId",&_lgenpdgId);
    }

    outputTree->Branch("HLT_Ele35_WPTight_Gsf",&HLT_Ele35_WPTight_Gsf,"HLT_Ele35_WPTight_Gsf/O");
    outputTree->Branch("HLT_Ele32_WPTight_Gsf",&HLT_Ele32_WPTight_Gsf,"HLT_Ele32_WPTight_Gsf/O");
    outputTree->Branch("HLT_Ele27_WPTight_Gsf",&HLT_Ele27_WPTight_Gsf,"HLT_Ele27_WPTight_Gsf/O");
    outputTree->Branch("HLT_IsoMu27",&HLT_IsoMu27,"HLT_IsoMu27/O");
    outputTree->Branch("HLT_IsoMu24",&HLT_IsoMu24,"HLT_IsoMu24/O");
    outputTree->Branch("HLT_IsoTkMu24",&HLT_IsoTkMu24,"HLT_IsoTkMu24/O");
    outputTree->Branch("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",&HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ,"HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ/O");
    outputTree->Branch("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",&HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ,"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ/O");
    outputTree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL/O");
    outputTree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ/O");
    outputTree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8/O");
    outputTree->Branch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL/O");
    outputTree->Branch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ/O");
  }


  if( Skim_=="GammaJetsResiduals"  ){
    outputTree->Branch("_ptgamma", &_ptgamma, "_ptgamma/F");
    outputTree->Branch("_phigamma", &_phigamma, "_phigamma/F");
    outputTree->Branch("_phEta",&_phEta);
    outputTree->Branch("_phPhi",&_phPhi);
    outputTree->Branch("_phPt",&_phPt);
    outputTree->Branch("_phPtcorr",&_phPtcorr);
    outputTree->Branch("_phPassTightID",&_phPassTightID);
    outputTree->Branch("_phPassIso",&_phPassIso);

    if(IsMC_){
    outputTree->Branch("_ptgamma_gen", &_ptgamma_gen, "_ptgamma_gen/F");
    outputTree->Branch("_phigamma_gen", &_phigamma_gen, "_phigamma_gen/F");
    outputTree->Branch("_phgenEta",&_phgenEta);
    outputTree->Branch("_phgenPhi",&_phgenPhi);
    outputTree->Branch("_phgenPt",&_phgenPt);
    }
    outputTree->Branch("HLT_Photon110EB_TightID_TightIso",&HLT_Photon110EB_TightID_TightIso,"HLT_Photon110EB_TightID_TightIso/O");
    outputTree->Branch("HLT_Photon165_R9Id90_HE10_IsoM",&HLT_Photon165_R9Id90_HE10_IsoM,"HLT_Photon165_R9Id90_HE10_IsoM/O");
    outputTree->Branch("HLT_Photon120_R9Id90_HE10_IsoM",&HLT_Photon120_R9Id90_HE10_IsoM,"HLT_Photon120_R9Id90_HE10_IsoM/O");
    outputTree->Branch("HLT_Photon90_R9Id90_HE10_IsoM",&HLT_Photon90_R9Id90_HE10_IsoM,"HLT_Photon90_R9Id90_HE10_IsoM/O");
    outputTree->Branch("HLT_Photon75_R9Id90_HE10_IsoM",&HLT_Photon75_R9Id90_HE10_IsoM,"HLT_Photon75_R9Id90_HE10_IsoM/O");
    outputTree->Branch("HLT_Photon50_R9Id90_HE10_IsoM",&HLT_Photon50_R9Id90_HE10_IsoM,"HLT_Photon50_R9Id90_HE10_IsoM/O");
    outputTree->Branch("HLT_Photon200",&HLT_Photon200,"HLT_Photon200/O");
    outputTree->Branch("HLT_Photon175",&HLT_Photon175,"HLT_Photon175/O");
  }
  if(Skim_=="ZJetsResiduals" || Skim_=="GammaJetsResiduals"  ) return;

  
  outputTree->Branch("_eventNb",   &_eventNb,   "_eventNb/l");
  outputTree->Branch("_runNb",     &_runNb,     "_runNb/l");
  outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");
  outputTree->Branch("_bx", &_bx, "_bx/l");
  outputTree->Branch("_n_PV", &_n_PV, "_n_PV/I");
  outputTree->Branch("_LV_x", &_LV_x, "_LV_x/F");
  outputTree->Branch("_LV_y", &_LV_y, "_LV_y/F");
  outputTree->Branch("_LV_z", &_LV_z, "_LV_z/F");
  outputTree->Branch("_LV_errx", &_LV_errx, "_LV_errx/F");
  outputTree->Branch("_LV_erry", &_LV_erry, "_LV_erry/F");
  outputTree->Branch("_LV_errz", &_LV_errz, "_LV_errz/F");


  outputTree->Branch("_rho", &_rho, "_rho/F");
  outputTree->Branch("_rhoNC", &_rhoNC, "_rhoNC/F");
  
  outputTree->Branch("Flag_IsUnprefirable",&Flag_IsUnprefirable,"Flag_IsUnprefirable/O");
  outputTree->Branch("Flag_goodVertices",&Flag_goodVertices,"Flag_goodVertices/O");
  outputTree->Branch("Flag_globalTightHalo2016Filter",&Flag_globalTightHalo2016Filter,"Flag_globalTightHalo2016Filter/O");
  outputTree->Branch("Flag_globalSuperTightHalo2016Filter",&Flag_globalSuperTightHalo2016Filter,"Flag_globalSuperTightHalo2016Filter/O");
  outputTree->Branch("Flag_HBHENoiseFilter",&Flag_HBHENoiseFilter,"Flag_HBHENoiseFilter/O");
  outputTree->Branch("Flag_HBHENoiseIsoFilter",&Flag_HBHENoiseIsoFilter,"Flag_HBHENoiseIsoFilter/O");
  outputTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter",&Flag_EcalDeadCellTriggerPrimitiveFilter,"Flag_EcalDeadCellTriggerPrimitiveFilter/O");
  outputTree->Branch("Flag_BadPFMuonFilter",&Flag_BadPFMuonFilter,"Flag_BadPFMuonFilter/O");
  outputTree->Branch("Flag_BadPFMuonDzFilter",&Flag_BadPFMuonDzFilter,"Flag_BadPFMuonDzFilter/O");
  outputTree->Branch("Flag_hfNoisyHitsFilter",&Flag_hfNoisyHitsFilter,"Flag_hfNoisyHitsFilter/O");
  outputTree->Branch("Flag_BadChargedCandidateFilter",&Flag_BadChargedCandidateFilter,"Flag_BadChargedCandidateFilter/O");
  outputTree->Branch("Flag_eeBadScFilter",&Flag_eeBadScFilter,"Flag_eeBadScFilter/O");
  outputTree->Branch("Flag_ecalBadCalibFilter",&Flag_ecalBadCalibFilter,"Flag_ecalBadCalibFilter/O");
  outputTree->Branch("Flag_ecalLaserCorrFilter",&Flag_ecalLaserCorrFilter,"Flag_ecalLaserCorrFilter/O");
  outputTree->Branch("Flag_EcalDeadCellBoundaryEnergyFilter",&Flag_EcalDeadCellBoundaryEnergyFilter,"Flag_EcalDeadCellBoundaryEnergyFilter/O");

  outputTree->Branch("PassecalBadCalibFilter_Update",&PassecalBadCalibFilter_Update,"PassecalBadCalibFilter_Update/O");
  outputTree->Branch("PassecalLaserCorrFilter_Update",&PassecalLaserCorrFilter_Update,"PassecalLaserCorrFilter_Update/O");
  outputTree->Branch("PassEcalDeadCellBoundaryEnergyFilter_Update",&PassEcalDeadCellBoundaryEnergyFilter_Update,"PassEcalDeadCellBoundaryEnergyFilter_Update/O");
  outputTree->Branch("PassBadChargedCandidateFilter_Update",&PassBadChargedCandidateFilter_Update,"PassBadChargedCandidateFilter_Update/O");


  outputTree->Branch("_jetEta",&_jetEta);
  outputTree->Branch("_jetPhi",&_jetPhi);
  outputTree->Branch("_jetPt",&_jetPt);
  outputTree->Branch("_jetRawPt",&_jetRawPt);
  
  if(Skim_!="MCJECs"){
  outputTree->Branch("_jet_CHEF",&_jet_CHEF);
  outputTree->Branch("_jet_NHEF",&_jet_NHEF);
  outputTree->Branch("_jet_NEEF",&_jet_NEEF);
  outputTree->Branch("_jet_CEEF",&_jet_CEEF);
  outputTree->Branch("_jet_MUEF",&_jet_MUEF);
  outputTree->Branch("_jet_CHM",&_jet_CHM);
  outputTree->Branch("_jet_NHM",&_jet_NHM);
  outputTree->Branch("_jet_PHM",&_jet_PHM);
  outputTree->Branch("_jet_NM",&_jet_NM);
  outputTree->Branch("_jethfsigmaEtaEta",&_jethfsigmaEtaEta);
  outputTree->Branch("_jethfsigmaPhiPhi",&_jethfsigmaPhiPhi);
  outputTree->Branch("_jethfcentralEtaStripSize",&_jethfcentralEtaStripSize);
  outputTree->Branch("_jethfadjacentEtaStripsSize",&_jethfadjacentEtaStripsSize);
  outputTree->Branch("_jetPtNoL2L3Res",&_jetPtNoL2L3Res);
  outputTree->Branch("_jet_corrjecs",&_jet_corrjecs);
  if(IsMC_)outputTree->Branch("_jethadronFlavour",&_jethadronFlavour);
  if(IsMC_)outputTree->Branch("_jetpartonFlavour",&_jetpartonFlavour);
  outputTree->Branch("_jetJECuncty",&_jetJECuncty);
  }
  outputTree->Branch("_jetArea",&_jetArea);
  outputTree->Branch("_jetPassID",&_jetPassID);
  outputTree->Branch("_jetLeptonPhotonCleaned",&_jetLeptonPhotonCleaned);
  outputTree->Branch("_jetTauCleaned",&_jetTauCleaned);

  if(IsMC_){
  outputTree->Branch("_jetPtGen",&_jetPtGen);
  outputTree->Branch("_jetEtaGen",&_jetEtaGen);
  outputTree->Branch("_jetPhiGen",&_jetPhiGen);
  outputTree->Branch("_jetPtGenWithNu",&_jetPtGenWithNu);
  }
  outputTree->Branch("_jetPUMVA",&_jetPUMVA);
  outputTree->Branch("_jetPUMVAUpdate",&_jetPUMVAUpdate);
  outputTree->Branch("_jetPUMVAUpdate2017",&_jetPUMVAUpdate2017);
  outputTree->Branch("_jetPUMVAUpdate2018",&_jetPUMVAUpdate2018);
  outputTree->Branch("_jetDeepJet_b",&_jetDeepJet_b);
  outputTree->Branch("_jetParticleNet_b",&_jetParticleNet_b);
  outputTree->Branch("_jetDeepJet_c",&_jetDeepJet_c);
  outputTree->Branch("_jetDeepJet_uds",&_jetDeepJet_uds);
  outputTree->Branch("_jetDeepJet_g",&_jetDeepJet_g);
  outputTree->Branch("_jetQuarkGluonLikelihood",&_jetQuarkGluonLikelihood);


  if(SavePUIDVariables_){ 
  outputTree->Branch("_jet_beta",&_jet_beta);
  outputTree->Branch("_jet_dR2Mean",&_jet_dR2Mean);
  outputTree->Branch("_jet_majW",&_jet_majW);
  outputTree->Branch("_jet_minW",&_jet_minW);
  outputTree->Branch("_jet_frac01",&_jet_frac01);
  outputTree->Branch("_jet_frac02",&_jet_frac02);
  outputTree->Branch("_jet_frac03",&_jet_frac03);
  outputTree->Branch("_jet_frac04",&_jet_frac04);
  outputTree->Branch("_jet_ptD",&_jet_ptD);
  outputTree->Branch("_jet_betaStar",&_jet_betaStar);
  outputTree->Branch("_jet_pull",&_jet_pull);
  outputTree->Branch("_jet_jetR",&_jet_jetR);
  outputTree->Branch("_jet_jetRchg",&_jet_jetRchg);
  outputTree->Branch("_jet_nParticles",&_jet_nParticles);
  outputTree->Branch("_jet_nCharged",&_jet_nCharged);
  }
  
  outputTree->Branch("_puppijetEta",&_puppijetEta);
  outputTree->Branch("_puppijetPhi",&_puppijetPhi);
  outputTree->Branch("_puppijetPt",&_puppijetPt);
  outputTree->Branch("_puppijetRawPt",&_puppijetRawPt);
  outputTree->Branch("_puppijetPtGen",&_puppijetPtGen);
  outputTree->Branch("_puppijetPtGenWithNu",&_puppijetPtGenWithNu);
  outputTree->Branch("_puppijetLeptonPhotonCleaned",&_puppijetLeptonPhotonCleaned);
  outputTree->Branch("_puppijetPassID",&_puppijetPassID);
  outputTree->Branch("_puppijetPtNoL2L3Res",&_puppijetPtNoL2L3Res);
  
  if(SaveAK8Jets_){
  outputTree->Branch("_puppiak8jetEta",&_puppiak8jetEta);
  outputTree->Branch("_puppiak8jetPhi",&_puppiak8jetPhi);
  outputTree->Branch("_puppiak8jetPt",&_puppiak8jetPt);
  outputTree->Branch("_puppiak8jetRawPt",&_puppiak8jetRawPt);
  outputTree->Branch("_puppiak8jetPtGen",&_puppiak8jetPtGen);
  outputTree->Branch("_puppiak8jet_tau1",&_puppiak8jet_tau1);
  outputTree->Branch("_puppiak8jet_tau2",&_puppiak8jet_tau2);
  outputTree->Branch("_puppiak8jet_tau3",&_puppiak8jet_tau3);
  outputTree->Branch("_ak8jetEta",&_ak8jetEta);
  outputTree->Branch("_ak8jetPhi",&_ak8jetPhi);
  outputTree->Branch("_ak8jetPt",&_ak8jetPt);
  outputTree->Branch("_ak8jetArea",&_ak8jetArea);
  outputTree->Branch("_ak8jetRawPt",&_ak8jetRawPt);
  outputTree->Branch("_ak8jetPtGen",&_ak8jetPtGen);

  }

  if(SaveCaloJets_){
    outputTree->Branch("_calojetEta",&_calojetEta);
    outputTree->Branch("_calojetPhi",&_calojetPhi);
    outputTree->Branch("_calojetPt",&_calojetPt);
    outputTree->Branch("_calojetRawPt",&_calojetRawPt);
    outputTree->Branch("_calojetPtGen",&_calojetPtGen);
    outputTree->Branch("_calojetLeptonPhotonCleaned",&_calojetLeptonPhotonCleaned);

  }
  
  if(SavenoCHSJets_){
    outputTree->Branch("_noCHSjetEta",&_noCHSjetEta);
    outputTree->Branch("_noCHSjetPhi",&_noCHSjetPhi);
    outputTree->Branch("_noCHSjetPt",&_noCHSjetPt);
    outputTree->Branch("_noCHSjetRawPt",&_noCHSjetRawPt);
    outputTree->Branch("_noCHSjetPtGen",&_noCHSjetPtGen);
    outputTree->Branch("_noCHSjetLeptonPhotonCleaned",&_noCHSjetLeptonPhotonCleaned);
  }

  
  
  if(Skim_=="MCJECs")
  {
  outputTree->Branch("_genjetEta",&_genjetEta);
  outputTree->Branch("_genjetPhi",&_genjetPhi);
  outputTree->Branch("_genjetPt",&_genjetPt);
  
  outputTree->Branch("_genjet_CHSIdx",&_genjet_CHSIdx);
  outputTree->Branch("_genjet_noCHSIdx",&_genjet_noCHSIdx);
  outputTree->Branch("_genjet_CaloIdx",&_genjet_CaloIdx);
  outputTree->Branch("_genjet_PuppiIdx",&_genjet_PuppiIdx);
 

  outputTree->Branch("_genAK8jetEta",&_genAK8jetEta);
  outputTree->Branch("_genAK8jetPhi",&_genAK8jetPhi);
  outputTree->Branch("_genAK8jetPt",&_genAK8jetPt);
  outputTree->Branch("_genAK8jet_PuppiIdx",&_genAK8jet_PuppiIdx);
  outputTree->Branch("_genAK8jet_CHSIdx",&_genAK8jet_CHSIdx);
 
  }
  

  if(PFCandPtCut_<1000){
  outputTree->Branch("_PFcand_pt",&_PFcand_pt);
  outputTree->Branch("_PFcand_eta",&_PFcand_eta);
  outputTree->Branch("_PFcand_phi",&_PFcand_phi);
  outputTree->Branch("_PFcand_pdgId",&_PFcand_pdgId);
  outputTree->Branch("_PFcand_fromPV",&_PFcand_fromPV);
  outputTree->Branch("_PFcand_dz",&_PFcand_dz);
  outputTree->Branch("_PFcand_dzError",&_PFcand_dzError);
  outputTree->Branch("_PFcand_hcalFraction",&_PFcand_hcalFraction);
  outputTree->Branch("_PFcand_PVfitidx",&_PFcand_PVfitidx);
  outputTree->Branch("_PFcand_puppiweight",&_PFcand_puppiweight);
  }




  
  if(IsMC_){
  outputTree->Branch("_genmet", &_genmet, "_genmet/F");
  outputTree->Branch("_genmet_phi", &_genmet_phi, "_genmet_phi/F");
  outputTree->Branch("trueNVtx", &trueNVtx,"trueNVtx/I");
  }
  outputTree->Branch("_met", &_met, "_met/F");
  outputTree->Branch("_met_phi", &_met_phi, "_met_phi/F");
  outputTree->Branch("_puppimet", &_puppimet, "_puppimet/F");
  outputTree->Branch("_puppimet_phi", &_puppimet_phi, "_puppimet_phi/F");

  outputTree->Branch("_rawmet", &_rawmet, "_rawmet/F");
  outputTree->Branch("_rawmet_phi", &_rawmet_phi, "_rawmet_phi/F");
  outputTree->Branch("_puppirawmet", &_puppirawmet, "_puppirawmet/F");
  outputTree->Branch("_puppirawmet_phi", &_puppirawmet_phi, "_puppirawmet_phi/F");
  
  outputTree->Branch("_rawchsmet",&_rawchsmet,"_rawchsmet/F");
  outputTree->Branch("_rawchsmet_phi",&_rawchsmet_phi,"_rawchsmet_phi/F");
  outputTree->Branch("_chsmet",&_chsmet,"_chsmet/F");
  outputTree->Branch("_chsmet_phi",&_chsmet_phi,"_chsmet_phi/F");

  if(SavePFinJets_){
    jetPFTree->Branch("_Jet_Pt", &_Jet_Pt,"_Jet_Pt/F");
    jetPFTree->Branch("_Jet_Eta", &_Jet_Eta,"_Jet_Eta/F");
    jetPFTree->Branch("_Jet_Phi", &_Jet_Phi,"_Jet_Phi/F");
    jetPFTree->Branch("_Jet_PtGen", &_Jet_PtGen,"_Jet_PtGen/F");
    jetPFTree->Branch("_Jet_EtaGen", &_Jet_EtaGen,"_Jet_EtaGen/F");
    jetPFTree->Branch("_Jet_PhiGen", &_Jet_PhiGen,"_Jet_PhiGen/F");

    jetPFTree->Branch("_PFcand_pt",&_Jet_PFcand_pt);
    jetPFTree->Branch("_PFcand_eta",&_Jet_PFcand_eta);
    jetPFTree->Branch("_PFcand_phi",&_Jet_PFcand_phi);
    jetPFTree->Branch("_PFcand_pdgId",&_Jet_PFcand_pdgId);
    jetPFTree->Branch("_PFcand_fromPV",&_Jet_PFcand_fromPV);
    jetPFTree->Branch("_PFcand_dz",&_Jet_PFcand_dz);
    jetPFTree->Branch("_PFcand_dzError",&_Jet_PFcand_dzError);
  }


  if(Skim_=="MCJECs")return;

  outputTree->Branch("_PUV1_x", &_PUV1_x, "_PUV1_x/F");
  outputTree->Branch("_PUV1_y", &_PUV1_y, "_PUV1_y/F");
  outputTree->Branch("_PUV1_z", &_PUV1_z, "_PUV1_z/F");
  outputTree->Branch("_n_CH_fromvtxfit",&_n_CH_fromvtxfit,"_n_CH_fromvtxfit[6]/I");
  outputTree->Branch("_HT_CH_fromvtxfit", &_HT_CH_fromvtxfit, "_HT_CH_fromvtxfit[6]/F");


  

  outputTree->Branch("_lEta",&_lEta);
  outputTree->Branch("_lPhi",&_lPhi);
  outputTree->Branch("_lPt",&_lPt);
  outputTree->Branch("_lPtcorr",&_lPtcorr);
  outputTree->Branch("_lpdgId",&_lpdgId);
  outputTree->Branch("_lPassTightID",&_lPassTightID);
  outputTree->Branch("_lPassLooseID",&_lPassLooseID);
  outputTree->Branch("_lisSAMuon",&_lisSAMuon);
  outputTree->Branch("_nEles", &_nEles, "_nEles/I");
  outputTree->Branch("_nMus", &_nMus, "_nMus/I");
  

  outputTree->Branch("_ldz",&_ldz);
  outputTree->Branch("_ldzError",&_ldzError);
  outputTree->Branch("_ldxy",&_ldxy);
  outputTree->Branch("_ldxyError",&_ldxyError);
   
  outputTree->Branch("_l3dIP",&_l3dIP);
  outputTree->Branch("_l3dIPError",&_l3dIPError);

  
  if(Skim_=="ZToEEorMuMu" || Skim_=="Dilepton" || Skim_=="DileptonInfo" || Skim_=="ZHFJet"){
    outputTree->Branch("_mll", &_mll, "_mll/F");
    outputTree->Branch("_ptll", &_ptll, "_ptll/F");
    outputTree->Branch("_pzll", &_pzll, "_pzll/F");
    outputTree->Branch("_yll", &_yll, "_yll/F");
    outputTree->Branch("_dphill", &_dphill, "_dphill/F");
    outputTree->Branch("_phill", &_phill, "_phill/F");
    outputTree->Branch("_costhCSll", &_costhCSll, "_costhCSll/F");
    outputTree->Branch("_nElesll",&_nElesll,"_nElesll/I");
    
    if(IsMC_){
    outputTree->Branch("_mll_gen", &_mll_gen, "_mll_gen/F");
    outputTree->Branch("_ptll_gen", &_ptll_gen, "_ptll_gen/F");
    outputTree->Branch("_pzll_gen", &_pzll_gen, "_pzll_gen/F");
    outputTree->Branch("_yll_gen", &_yll_gen, "_yll_gen/F");
    outputTree->Branch("_dphill_gen", &_dphill_gen, "_dphill_gen/F");
    outputTree->Branch("_phill_gen", &_phill_gen, "_phill_gen/F");
    outputTree->Branch("_costhCSll_gen", &_costhCSll_gen, "_costhCSll_gen/F");
    outputTree->Branch("_ngenElesll",&_ngenElesll,"_ngenElesll/I");
    }
  }
  outputTree->Branch("_n_PFele_fromvtxfit",&_n_PFele_fromvtxfit,"_n_PFele_fromvtxfit/I");
  outputTree->Branch("_n_PFmu_fromvtxfit",&_n_PFmu_fromvtxfit,"_n_PFmu_fromvtxfit/I");


  if(IsMC_){
    outputTree->Branch("_lgenEta",&_lgenEta);
    outputTree->Branch("_lgenPhi",&_lgenPhi);
    outputTree->Branch("_lgenPt",&_lgenPt);
    outputTree->Branch("_lgenpdgId",&_lgenpdgId);
    outputTree->Branch("_phgenEta",&_phgenEta);
    outputTree->Branch("_phgenPhi",&_phgenPhi);
    outputTree->Branch("_phgenPt",&_phgenPt);
    outputTree->Branch("_genHT",&_genHT,"_genHT/F");
    outputTree->Branch("_weight",&_weight,"_weight/F");
    

  }
  
  if(PhotonPtCut_<1000){
  outputTree->Branch("_phEta",&_phEta);
  outputTree->Branch("_phPhi",&_phPhi);
  outputTree->Branch("_phPt",&_phPt);
  outputTree->Branch("_phPtcorr",&_phPtcorr);
  outputTree->Branch("_phPassTightID",&_phPassTightID);
  outputTree->Branch("_phPassIso",&_phPassIso);
  }



  if(Skim_=="VtxInfo" ){
  outputTree->Branch("_METCH_PV",&_METCH_PV);
  outputTree->Branch("_METPhiCH_PV",&_METPhiCH_PV);
  outputTree->Branch("_SumPT2CH_PV",&_SumPT2CH_PV);
  outputTree->Branch("_DztoLV_PV",&_DztoLV_PV);
  }




  outputTree->Branch("HLT_Photon110EB_TightID_TightIso",&HLT_Photon110EB_TightID_TightIso,"HLT_Photon110EB_TightID_TightIso/O");
  outputTree->Branch("HLT_Photon165_R9Id90_HE10_IsoM",&HLT_Photon165_R9Id90_HE10_IsoM,"HLT_Photon165_R9Id90_HE10_IsoM/O");
  outputTree->Branch("HLT_Photon120_R9Id90_HE10_IsoM",&HLT_Photon120_R9Id90_HE10_IsoM,"HLT_Photon120_R9Id90_HE10_IsoM/O");
  outputTree->Branch("HLT_Photon90_R9Id90_HE10_IsoM",&HLT_Photon90_R9Id90_HE10_IsoM,"HLT_Photon90_R9Id90_HE10_IsoM/O");
  outputTree->Branch("HLT_Photon75_R9Id90_HE10_IsoM",&HLT_Photon75_R9Id90_HE10_IsoM,"HLT_Photon75_R9Id90_HE10_IsoM/O");
  outputTree->Branch("HLT_Photon50_R9Id90_HE10_IsoM",&HLT_Photon50_R9Id90_HE10_IsoM,"HLT_Photon50_R9Id90_HE10_IsoM/O");
  outputTree->Branch("HLT_Photon200",&HLT_Photon200,"HLT_Photon200/O");
  outputTree->Branch("HLT_Photon175",&HLT_Photon175,"HLT_Photon175/O");
  outputTree->Branch("HLT_DiJet110_35_Mjj650_PFMET110",&HLT_DiJet110_35_Mjj650_PFMET110,"HLT_DiJet110_35_Mjj650_PFMET110/O");
  outputTree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60/O");
  outputTree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight,"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight/O");
  outputTree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF",&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF,"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF/O");
  outputTree->Branch("HLT_PFMET120_PFMHT120_IDTight_PFHT60",&HLT_PFMET120_PFMHT120_IDTight_PFHT60,"HLT_PFMET120_PFMHT120_IDTight_PFHT60/O");
  outputTree->Branch("HLT_PFMET120_PFMHT120_IDTight",&HLT_PFMET120_PFMHT120_IDTight,"HLT_PFMET120_PFMHT120_IDTight/O");
  outputTree->Branch("HLT_PFHT1050",&HLT_PFHT1050,"HLT_PFHT1050/O");
  outputTree->Branch("HLT_PFHT900",&HLT_PFHT900,"HLT_PFHT900/O");
  outputTree->Branch("HLT_PFJet500",&HLT_PFJet500,"HLT_PFJet500/O");
  outputTree->Branch("HLT_AK8PFJet500",&HLT_AK8PFJet500,"HLT_AK8PFJet500/O");
  outputTree->Branch("HLT_Ele35_WPTight_Gsf",&HLT_Ele35_WPTight_Gsf,"HLT_Ele35_WPTight_Gsf/O");
  outputTree->Branch("HLT_Ele32_WPTight_Gsf",&HLT_Ele32_WPTight_Gsf,"HLT_Ele32_WPTight_Gsf/O");
  outputTree->Branch("HLT_Ele27_WPTight_Gsf",&HLT_Ele27_WPTight_Gsf,"HLT_Ele27_WPTight_Gsf/O");
  outputTree->Branch("HLT_IsoMu27",&HLT_IsoMu27,"HLT_IsoMu27/O");
  outputTree->Branch("HLT_IsoMu24",&HLT_IsoMu24,"HLT_IsoMu24/O");
  outputTree->Branch("HLT_IsoTkMu24",&HLT_IsoTkMu24,"HLT_IsoTkMu24/O");
  outputTree->Branch("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",&HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ,"HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ/O");
  outputTree->Branch("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",&HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ,"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ/O");
  outputTree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL/O");
  outputTree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ/O");
  outputTree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8/O");
  outputTree->Branch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL/O");
  outputTree->Branch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ/O");
  outputTree->Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",&HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ/O");
  outputTree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ/O");
  outputTree->Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",&HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL/O");
  outputTree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL/O");

  if(!IsMC_)outputTree->Branch("_l1prefire",&_l1prefire,"_l1prefire/O");

  if(Skim_=="L1Unprefirable" || Skim_==""){
  outputTree->Branch("_L1mu_Qual",&_L1mu_Qual);
  outputTree->Branch("_L1mu_pt",&_L1mu_pt);
  outputTree->Branch("_L1mu_eta",&_L1mu_eta);
  outputTree->Branch("_L1mu_phi",&_L1mu_phi);
  outputTree->Branch("_L1mu_bx",&_L1mu_bx);
  outputTree->Branch("_L1mu_TFIdx",&_L1mu_TFIdx);
  outputTree->Branch("_L1mu_upt",&_L1mu_upt);
  outputTree->Branch("_L1mu_charge",&_L1mu_charge);
  outputTree->Branch("_L1mu_dXY",&_L1mu_dXY);


  outputTree->Branch("_L1eg_pt",&_L1eg_pt);
  outputTree->Branch("_L1eg_eta",&_L1eg_eta);
  outputTree->Branch("_L1eg_phi",&_L1eg_phi);
  outputTree->Branch("_L1eg_bx",&_L1eg_bx);
  outputTree->Branch("_L1eg_iso",&_L1eg_iso);


  outputTree->Branch("_L1jet_pt",&_L1jet_pt);
  outputTree->Branch("_L1jet_eta",&_L1jet_eta);
  outputTree->Branch("_L1jet_phi",&_L1jet_phi);
  outputTree->Branch("_L1jet_bx",&_L1jet_bx);


  outputTree->Branch("HLT_PixelClusters_WP1",&HLT_PixelClusters_WP1,"HLT_PixelClusters_WP1/O");
  outputTree->Branch("HLT_PixelClusters_WP2",&HLT_PixelClusters_WP2,"HLT_PixelClusters_WP2/O");
  outputTree->Branch("HLT_PixelClusters_WP2_split",&HLT_PixelClusters_WP2_split,"HLT_PixelClusters_WP2_split/O");


  }
  

  outputTree->Branch("hltEle27WPTightGsfTrackIsoFilter",&hltEle27WPTightGsfTrackIsoFilter);
  outputTree->Branch("hltEle35noerWPTightGsfTrackIsoFilter",&hltEle35noerWPTightGsfTrackIsoFilter);
  outputTree->Branch("hltEle32WPTightGsfTrackIsoFilter",&hltEle32WPTightGsfTrackIsoFilter);

  outputTree->Branch("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09",&hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09);
  outputTree->Branch("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09",&hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09);

  outputTree->Branch("hltEG175HEFilter",&hltEG175HEFilter );
  outputTree->Branch("hltEG200HEFilter",&hltEG200HEFilter );
  outputTree->Branch("hltEG165R9Id90HE10IsoMTrackIsoFilter",&hltEG165R9Id90HE10IsoMTrackIsoFilter );
  outputTree->Branch("hltEG120R9Id90HE10IsoMTrackIsoFilter",&hltEG120R9Id90HE10IsoMTrackIsoFilter );
  outputTree->Branch("hltEG90R9Id90HE10IsoMTrackIsoFilter",&hltEG90R9Id90HE10IsoMTrackIsoFilter );
  outputTree->Branch("hltEG75R9Id90HE10IsoMTrackIsoFilter",&hltEG75R9Id90HE10IsoMTrackIsoFilter );
  outputTree->Branch("hltEG50R9Id90HE10IsoMTrackIsoFilter",&hltEG50R9Id90HE10IsoMTrackIsoFilter );
  outputTree->Branch("hltEG110EBTightIDTightIsoTrackIsoFilter",&hltEG110EBTightIDTightIsoTrackIsoFilter);


  outputTree->Branch("hltEle17WPLoose1GsfTrackIsoFilterForHI",&hltEle17WPLoose1GsfTrackIsoFilterForHI);
  outputTree->Branch("hltHIPhoton40Eta3p1",&hltHIPhoton40Eta3p1);
  outputTree->Branch("hltL3fL1sMu10lqL1f0L2f10L3Filtered12",&hltL3fL1sMu10lqL1f0L2f10L3Filtered12);
  outputTree->Branch("hltL3fL1sMu10lqL1f0L2f10L3Filtered15",&hltL3fL1sMu10lqL1f0L2f10L3Filtered15);
  outputTree->Branch("hltEle15WPLoose1GsfTrackIsoFilterForHI",&hltEle15WPLoose1GsfTrackIsoFilterForHI);

  outputTree->Branch("hltSinglePFJet60",&hltSinglePFJet60);
  outputTree->Branch("hltSinglePFJet80",&hltSinglePFJet80);
  outputTree->Branch("hltSinglePFJet140",&hltSinglePFJet140);
  outputTree->Branch("hltSinglePFJet200",&hltSinglePFJet200);
  outputTree->Branch("hltSinglePFJet260",&hltSinglePFJet260);
  outputTree->Branch("hltSinglePFJet320",&hltSinglePFJet320);
  outputTree->Branch("hltSinglePFJet400",&hltSinglePFJet400);
  outputTree->Branch("hltSinglePFJet450",&hltSinglePFJet450);
  outputTree->Branch("hltSinglePFJet500",&hltSinglePFJet500);



}

// ------------ method called once each job just after ending the event loop  ------------
void
JMEAnalyzer::endJob()
{
}





// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JMEAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

bool JMEAnalyzer::GetMETFilterDecision(const edm::Event& iEvent,edm::Handle<TriggerResults> METFilterResults, TString studiedfilter){
  
  if( !METFilterResults.failedToGet() ) {
    int N_MetFilters = METFilterResults->size();
    const edm::TriggerNames & metfilterName = iEvent.triggerNames(*METFilterResults);
    for( int i_Metfilter = 0; i_Metfilter < N_MetFilters; ++i_Metfilter ) {
      TString MetfilterPath =metfilterName.triggerName(i_Metfilter);
      //      cout << MetfilterPath<<endl;
      if(MetfilterPath.Index(studiedfilter) >=0)  return METFilterResults.product()->accept(i_Metfilter);

    }
  }
   return true; 
}


bool JMEAnalyzer::GetIdxFilterDecision(int it){
  if(it== idx_Flag_goodVertices)return  Flag_goodVertices;
  else if(it==  idx_Flag_globalTightHalo2016Filter)return   Flag_globalTightHalo2016Filter;
  else if(it==  idx_Flag_globalSuperTightHalo2016Filter)return   Flag_globalSuperTightHalo2016Filter;
  else if(it==  idx_Flag_HBHENoiseFilter)return   Flag_HBHENoiseFilter;
  else if(it==  idx_Flag_HBHENoiseIsoFilter)return   Flag_HBHENoiseIsoFilter;
  else if(it==  idx_Flag_EcalDeadCellTriggerPrimitiveFilter)return   Flag_EcalDeadCellTriggerPrimitiveFilter;
  else if(it==  idx_Flag_BadPFMuonFilter)return   Flag_BadPFMuonFilter;
  else if(it==  idx_Flag_BadPFMuonDzFilter)return   Flag_BadPFMuonDzFilter;
  else if(it==  idx_Flag_hfNoisyHitsFilter)return   Flag_hfNoisyHitsFilter;
  else if(it==  idx_Flag_BadChargedCandidateFilter)return   Flag_BadChargedCandidateFilter;
  else if(it==  idx_Flag_eeBadScFilter)return   Flag_eeBadScFilter;
  else if(it==  idx_Flag_ecalBadCalibFilter)return   Flag_ecalBadCalibFilter;
  else if(it==  idx_Flag_ecalLaserCorrFilter)return   Flag_ecalLaserCorrFilter;
  else if(it==  idx_Flag_EcalDeadCellBoundaryEnergyFilter)return   Flag_EcalDeadCellBoundaryEnergyFilter;
  else if(it==  idx_PassecalBadCalibFilter_Update)return   PassecalBadCalibFilter_Update;
  else if(it==  idx_PassecalLaserCorrFilter_Update)return   PassecalLaserCorrFilter_Update;
  else if(it==  idx_PassEcalDeadCellBoundaryEnergyFilter_Update)return   PassEcalDeadCellBoundaryEnergyFilter_Update;
  else if(it==  idx_PassBadChargedCandidateFilter_Update)return   PassBadChargedCandidateFilter_Update;
  else return false;
}

TString JMEAnalyzer::GetIdxFilterName(int it){
  if(it==idx_Flag_goodVertices)return "Flag_goodVertices";
  else if(it== idx_Flag_globalTightHalo2016Filter)return "Flag_globalTightHalo2016Filter";
  else if(it== idx_Flag_globalSuperTightHalo2016Filter)return "Flag_globalSuperTightHalo2016Filter";
  else if(it== idx_Flag_HBHENoiseFilter)return "Flag_HBHENoiseFilter";
  else if(it== idx_Flag_HBHENoiseIsoFilter)return "Flag_HBHENoiseIsoFilter";
  else if(it== idx_Flag_EcalDeadCellTriggerPrimitiveFilter)return "Flag_EcalDeadCellTriggerPrimitiveFilter";
  else if(it== idx_Flag_BadPFMuonFilter)return "Flag_BadPFMuonFilter";
  else if(it== idx_Flag_BadPFMuonDzFilter)return "Flag_BadPFMuonDzFilter";
  else if(it== idx_Flag_hfNoisyHitsFilter)return "Flag_hfNoisyHitsFilter";
  else if(it== idx_Flag_BadChargedCandidateFilter)return "Flag_BadChargedCandidateFilter";
  else if(it== idx_Flag_eeBadScFilter)return "Flag_eeBadScFilter";
  else if(it== idx_Flag_ecalBadCalibFilter)return "Flag_ecalBadCalibFilter";
  else if(it== idx_Flag_ecalLaserCorrFilter)return "Flag_ecalLaserCorrFilter";
  else if(it== idx_Flag_EcalDeadCellBoundaryEnergyFilter)return "Flag_EcalDeadCellBoundaryEnergyFilter";
  else if(it== idx_PassecalBadCalibFilter_Update)return "PassecalBadCalibFilter_Update";
  else if(it== idx_PassecalLaserCorrFilter_Update)return "PassecalLaserCorrFilter_Update";
  else if(it== idx_PassEcalDeadCellBoundaryEnergyFilter_Update )return "PassEcalDeadCellBoundaryEnergyFilter_Update";
  else if(it== idx_PassBadChargedCandidateFilter_Update )return "PassBadChargedCandidateFilter_Update";
  else return "";
}


void JMEAnalyzer::InitandClearStuff(){
  _mll=0;
  _ptll=0;
  _pzll=0;
  _yll=0;
  _dphill=0;
  _phill=0;
  _costhCSll=0;
  _nElesll =0;
  
  _mll_gen=0;
  _ptll_gen=0;
  _pzll_gen=0;
  _yll_gen=0;
  _dphill_gen=0;
  _costhCSll_gen=0;
  _phill_gen=0;
  _ngenElesll =0;
  
  if(!IsMC_) DropUnmatchedJets_=false;

  Flag_goodVertices=false;
  Flag_globalTightHalo2016Filter=false;
  Flag_globalSuperTightHalo2016Filter=false;
  Flag_HBHENoiseFilter=false;
  Flag_HBHENoiseIsoFilter=false;
  Flag_EcalDeadCellTriggerPrimitiveFilter=false;
  Flag_BadPFMuonFilter=false;
  Flag_BadPFMuonDzFilter=false;
  Flag_hfNoisyHitsFilter=false;
  Flag_BadChargedCandidateFilter=false;
  Flag_eeBadScFilter=false;
  Flag_ecalBadCalibFilter=false;
  Flag_ecalLaserCorrFilter=false;
  Flag_EcalDeadCellBoundaryEnergyFilter=false;
  PassecalBadCalibFilter_Update=false;
  PassecalLaserCorrFilter_Update=false;
  PassEcalDeadCellBoundaryEnergyFilter_Update=false;
  PassBadChargedCandidateFilter_Update=false;
  
  _jetEta.clear();
  _jetPhi.clear();
  _jetPt.clear();
  _jetRawPt.clear();
  _jet_CHEF.clear();
  _jet_NHEF.clear();
  _jet_NEEF.clear();
  _jet_CEEF.clear();
  _jet_MUEF.clear();
  _jet_CHM.clear();
  _jet_NHM.clear();
  _jet_PHM.clear();
  _jet_NM.clear();
  _jetArea.clear();
  _jetPassID.clear();
  _jetLeptonPhotonCleaned.clear();
  _jetTauCleaned.clear();
  _jethfsigmaEtaEta.clear();
  _jethfsigmaPhiPhi.clear();
  _jethfcentralEtaStripSize.clear();
  _jethfadjacentEtaStripsSize.clear();
  _jetPtGen.clear();
  _jetEtaGen.clear();
  _jetPhiGen.clear();
  _jetPtGenWithNu.clear();
  _jetJECuncty.clear();
  _jetPUMVA.clear();
  _jetPUMVAUpdate.clear();
  _jetPUMVAUpdate2017.clear();
  _jetPUMVAUpdate2018.clear();

  _jetPtNoL2L3Res.clear();
  _jet_corrjecs.clear();
  _jetpartonFlavour.clear();
  _jethadronFlavour.clear();
  _jetDeepJet_b.clear();
  _jetParticleNet_b.clear();
  _jetDeepJet_c.clear();
  _jetDeepJet_uds.clear();
  _jetDeepJet_g.clear();
  _jetQuarkGluonLikelihood.clear();




  _jet_beta.clear();
  _jet_dR2Mean.clear();
  _jet_majW.clear();
  _jet_minW.clear();
  _jet_frac01.clear();
  _jet_frac02.clear();
  _jet_frac03.clear();
  _jet_frac04.clear();
  _jet_ptD.clear();
  _jet_betaStar.clear();
  _jet_pull.clear();
  _jet_jetR.clear();
  _jet_jetRchg.clear();
  _jet_nParticles.clear();
  _jet_nCharged.clear();

  _puppijetEta.clear();
  _puppijetPhi.clear();
  _puppijetPt.clear(); 
  _puppijetRawPt.clear(); 
  _puppijetPtGen.clear();
  _puppijetPtGenWithNu.clear();
  _puppijetLeptonPhotonCleaned.clear();
  _puppijetPassID.clear();
  _puppijetJECuncty.clear();
  _puppijetPtNoL2L3Res.clear();

  _puppiak8jetEta.clear();
  _puppiak8jetPhi.clear();
  _puppiak8jetPt.clear();
  _puppiak8jetRawPt.clear();
  _puppiak8jetPtGen.clear();
  _puppiak8jet_tau1.clear();
  _puppiak8jet_tau2.clear();
  _puppiak8jet_tau3.clear();

  _ak8jetEta.clear();
  _ak8jetPhi.clear();
  _ak8jetPt.clear();
  _ak8jetArea.clear();
  _ak8jetRawPt.clear();
  _ak8jetPtGen.clear();

  _calojetEta.clear();
  _calojetPhi.clear();
  _calojetPt.clear();
  _calojetRawPt.clear();
  _calojetPtGen.clear();
  _calojetLeptonPhotonCleaned.clear();

  _noCHSjetEta.clear();
  _noCHSjetPhi.clear();
  _noCHSjetPt.clear();
  _noCHSjetRawPt.clear();
  _noCHSjetPtGen.clear();
  _noCHSjetLeptonPhotonCleaned.clear();

  
  _genjetEta.clear();
  _genjetPhi.clear();
  _genjetPt.clear();
  _genjet_CHSIdx.clear();
  _genjet_noCHSIdx.clear();
  _genjet_CaloIdx.clear();
  _genjet_PuppiIdx.clear();


  _genAK8jetEta.clear();
  _genAK8jetPhi.clear();
  _genAK8jetPt.clear();
  _genAK8jet_CHSIdx.clear();
  _genAK8jet_PuppiIdx.clear();


  _tauEta.clear();
  _tauPhi.clear();
  _tauPt.clear();
  _tauPassMediumID.clear();

  _lEta.clear();
  _lPhi.clear();
  _lPt.clear();
  _lPtcorr.clear();
  _lPtSC.clear();
  _lpdgId.clear();
  _lPassTightID.clear();
  _lPassLooseID.clear();
  _lisSAMuon.clear();

  _ldz.clear();
  _ldzError.clear();
  _ldxy.clear();
  _ldxyError.clear();
  _l3dIP.clear();
  _l3dIPError.clear();
  
  _lpassHLT_IsoMu24.clear();
  _lpassHLT_Ele32_WPTight_Gsf.clear();

  _nEles=0;
  _nMus=0;
  

  _lgenEta.clear();
  _lgenPhi.clear();
  _lgenPt.clear();
  _lgenpdgId.clear();

  _phEta.clear();
  _phPhi.clear();
  _phPt.clear();
  _phPtcorr.clear();
  _phPassTightID.clear();
  _phPassIso.clear();
  _phgenEta.clear();
  _phgenPhi.clear();
  _phgenPt.clear();


  _PFcand_pt.clear();
  _PFcand_eta.clear();
  _PFcand_phi.clear();
  _PFcand_pdgId.clear();
  _PFcand_fromPV.clear();
  _PFcand_dz.clear();
  _PFcand_dzError.clear();
  _PFcand_hcalFraction.clear();
  _PFcand_PVfitidx.clear();
  _PFcand_puppiweight.clear();


  _Jet_PFcand_pt.clear();
  _Jet_PFcand_eta.clear();
  _Jet_PFcand_phi.clear();
  _Jet_PFcand_pdgId.clear();
  _Jet_PFcand_fromPV.clear();
  _Jet_PFcand_dz.clear();
  _Jet_PFcand_dzError.clear();
  


  _METCH_PV.clear();
  _METPhiCH_PV.clear();
  _SumPT2CH_PV.clear();
  _DztoLV_PV.clear();

  
  _L1mu_Qual.clear();
  _L1mu_pt.clear();
  _L1mu_eta.clear();
  _L1mu_phi.clear();
  _L1mu_bx.clear();
  _L1mu_TFIdx.clear();
  _L1mu_dXY.clear();
  _L1mu_upt.clear();
  _L1mu_charge.clear();
  _L1mu_TFIdx.clear();




  _L1eg_pt.clear();
  _L1eg_eta.clear();
  _L1eg_phi.clear();
  _L1eg_bx.clear();
  _L1eg_iso.clear();

  _L1jet_pt.clear();
  _L1jet_eta.clear();
  _L1jet_phi.clear();
  _L1jet_bx.clear();


  HLT_Photon110EB_TightID_TightIso=false;
  HLT_Photon165_R9Id90_HE10_IsoM=false;
  HLT_Photon120_R9Id90_HE10_IsoM=false;
  HLT_Photon90_R9Id90_HE10_IsoM=false;
  HLT_Photon75_R9Id90_HE10_IsoM=false;
  HLT_Photon50_R9Id90_HE10_IsoM=false;
  HLT_Photon200=false;
  HLT_Photon175=false;
  HLT_DiJet110_35_Mjj650_PFMET110=false;
  HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60=false;
  HLT_PFMETNoMu120_PFMHTNoMu120_IDTight=false;
  HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF=false;
  HLT_PFMET120_PFMHT120_IDTight_PFHT60=false;
  HLT_PFMET120_PFMHT120_IDTight=false;
  HLT_PFHT1050=false;
  HLT_PFHT900=false;
  HLT_PFJet500=false;
  HLT_AK8PFJet500=false;
  HLT_Ele35_WPTight_Gsf =false;
  HLT_Ele32_WPTight_Gsf =false;
  HLT_Ele27_WPTight_Gsf=false;
  HLT_IsoMu27=false;
  HLT_IsoMu24=false;
  HLT_IsoTkMu24=false;
  HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ=false;
  HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ=false;
  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL=false;
  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ=false;
  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8=false;
  HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL=false;
  HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ=false;
  HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ=false;
  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ=false;
  HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL=false;
  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL=false;  


  
  hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09.clear();
  hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09.clear();
  hltEle27WPTightGsfTrackIsoFilter.clear(); 
  hltEle35noerWPTightGsfTrackIsoFilter.clear(); 
  hltEle32WPTightGsfTrackIsoFilter.clear(); 
  hltEG175HEFilter.clear(); 
  hltEG200HEFilter.clear(); 
  hltEG165R9Id90HE10IsoMTrackIsoFilter.clear(); 
  hltEG120R9Id90HE10IsoMTrackIsoFilter.clear(); 
  hltEG90R9Id90HE10IsoMTrackIsoFilter.clear(); 
  hltEG75R9Id90HE10IsoMTrackIsoFilter.clear(); 
  hltEG50R9Id90HE10IsoMTrackIsoFilter.clear(); 
  hltEG110EBTightIDTightIsoTrackIsoFilter.clear();
  hltSinglePFJet60.clear(); 
  hltSinglePFJet80.clear(); 
  hltSinglePFJet140.clear();
  hltSinglePFJet200.clear(); 
  hltSinglePFJet260.clear(); 
  hltSinglePFJet320.clear(); 
  hltSinglePFJet400.clear(); 
  hltSinglePFJet450.clear(); 
  hltSinglePFJet500.clear(); 


  hltEle17WPLoose1GsfTrackIsoFilterForHI.clear();
  hltHIPhoton40Eta3p1.clear();
  hltL3fL1sMu10lqL1f0L2f10L3Filtered12.clear();
  hltL3fL1sMu10lqL1f0L2f10L3Filtered15.clear();
  hltEle15WPLoose1GsfTrackIsoFilterForHI.clear();

}

bool JMEAnalyzer::PassSkim(){
  
  if(Skim_=="L1Study") return Flag_goodVertices;

  if((Skim_=="L1Study_ZToMuMu" || Skim_=="L1Study_SingleMuforJME") && !HLT_IsoMu24) return false;
  if(Skim_=="L1Study_ZToEE"&&!HLT_Ele32_WPTight_Gsf) return false;
  if(Skim_=="L1Study_SinglePhotonforJME"&& !HLT_Photon50_R9Id90_HE10_IsoM&& !HLT_Photon110EB_TightID_TightIso) return false;
  
  if(Skim_=="ZToEEorMuMu" || Skim_=="Dilepton" || Skim_=="L1Study_ZToMuMu" || Skim_=="L1Study_ZToEE"){


    if(_mll>20&&Skim_=="Dilepton")  return true;
    if(_mll>70&&_mll<110) return true;
    return false;
  }
  else if(Skim_=="Photon" || Skim_=="L1Study_SinglePhotonforJME"){
    
    int ngoodphotons =0;
    for(unsigned int i = 0; i < _phPt.size(); i++){
      if(_phPt[i]<20) continue;
      if(fabs(_phEta[i])>1.4442) continue;
      ngoodphotons++;
    }
    if( ngoodphotons>0) return true;
    return false;
  }
  else if(Skim_=="MET100"&&_met<100) return false;
  else if(Skim_=="MET100"&&_met>100) return true;
  else if(Skim_=="FourLeptons" &&_lPt.size()<4 ) return false;
  else if(Skim_=="FourLeptons" &&_lPt.size() >=4) {
    int nl= 0;
    for(unsigned int i = 0; i < _lPt.size(); i++){
      if(!_lPassLooseID[i])  continue;
      nl++;
    }
    if(nl>=4){
      //      cout << "found " <<endl;
      return true;
    }
    else return false;
  }
  else if(Skim_=="HighHT") return (HLT_PFHT1050 || HLT_PFHT900 || HLT_PFJet500 || HLT_AK8PFJet500) ; 
  else if(Skim_=="L1Unprefirable" ){
    std::string str_run = std::to_string(_runNb);
    std::string str_lumi = std::to_string(_lumiBlock);
    std::string str_event = std::to_string(_eventNb);
    string line;

    if(_runNb == 305064 && _lumiBlock == 128 && _eventNb==202006565) cout << _runNb<<", " <<_lumiBlock <<", " <<_eventNb <<endl;
    myfile_unprefevts.clear();
    myfile_unprefevts.seekg(0);
    std::string stringtosearch = "run:"+str_run+":"+str_lumi+":"+str_event;
    //    cout << stringtosearch<<endl;
    //if(_runNb == 305248&& _lumiBlock == 640 && _eventNb== 1141243194)  cout <<" _runNb == 305248&& _lumiBlock == 640 && _eventNb== 1141243194 "<<endl;
    if (myfile_unprefevts.is_open()) {
      
      while ( getline (myfile_unprefevts,line) ) {
	//cout << "line in file is " <<endl;
	//cout << line<<endl;
	
	//	if(_runNb == 305248&& _lumiBlock == 640 && _eventNb== 1141243194)  cout <<"in line " <<line <<endl;
	if(line.find(stringtosearch) != string::npos) {
	  myfile_unprefevts.clear();
	  //cout << line<<endl;
	  return true; 
	}
      }
      
    }
    //    else cout << "myfile_unprefevts.is_open() = false " << endl;
    return false;   
  }

  else if(Skim_=="ZJetsResiduals" || Skim_=="GammaJetsResiduals"  ){
    
    if(Skim_=="ZJetsResiduals" && (_ptll==0|| _mll<70 || _mll>110 || _phPt.size()>0)) return false;
    if(Skim_=="GammaJetsResiduals"  && (_ptgamma<50 || _phPt.size()>1  ||_lPt.size()>0 ) ) return false; 
    

    int ntotjets(0), ntotpuppijets(0);
    for(unsigned int i = 0; i < _jetPt.size(); i++){
      if( !_jetPassID[i] ) continue; 
      if( !_jetLeptonPhotonCleaned[i] ) continue; 
      if(_jetPtNoL2L3Res[i] <30) continue;
      if(_jetPtNoL2L3Res[i] <40&&abs(_jetEta[i])>2.4 ) continue;

      ntotjets++;
    }

    for(unsigned int i = 0; i < _puppijetPt.size(); i++){
      if( !_puppijetPassID[i] ) continue;
      if( !_puppijetLeptonPhotonCleaned[i] ) continue;
      if(_puppijetPtNoL2L3Res[i] <30) continue;
      if(_puppijetPtNoL2L3Res[i] <40&&abs(_puppijetEta[i])>2.4 ) continue;
      
      ntotpuppijets++;
    }

    if(ntotjets>1&& ntotpuppijets>1) return false;
  }
  else if(Skim_=="HFJet" || Skim_=="PhotonHFJet" || Skim_=="ZHFJet"){
    int nhfjet=0;
    for(unsigned int i = 0; i < _jetPt.size(); i++){
      if(fabs(_jetEta[i])>2.95&&_jetPt[i]>70)nhfjet++;
    }
  
    if ( Skim_=="HFJet" && nhfjet >0 ) return true;
    if ( Skim_=="ZHFJet" && nhfjet >0 &&_ptll>50&&_mll>70&&_mll<110 && (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_IsoMu27 || HLT_IsoMu24 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL )) return true;
    int ngoodphotons =0;
    for(unsigned int i = 0; i < _phPt.size(); i++){
      if(_phPt[i]<50) continue;
      if(fabs(_phEta[i])>1.4442) continue;
      ngoodphotons++;
    }
    if ( Skim_=="PhotonHFJet"  && nhfjet>0&& _met<50&& ngoodphotons >0 && (  HLT_Photon110EB_TightID_TightIso  ||HLT_Photon165_R9Id90_HE10_IsoM  ||HLT_Photon120_R9Id90_HE10_IsoM  ||HLT_Photon90_R9Id90_HE10_IsoM || HLT_Photon75_R9Id90_HE10_IsoM  ||HLT_Photon50_R9Id90_HE10_IsoM  ||HLT_Photon200  ||HLT_Photon175) ) return true;
    return false;

  }

  return true; 
}

void
JMEAnalyzer::CalcDileptonInfo(const int& i, const int& j, Float_t & themass, Float_t & theptll, Float_t & thepzll,  Float_t & theyll, Float_t & thephill, Float_t & thedphill, Float_t & thecosthll){

  TLorentzVector lep1;
  TLorentzVector lep2;

  Float_t et1 = (_lPt)[i];
  Float_t et2 = (_lPt)[j];
  lep1.SetPtEtaPhiE(et1, (_lEta)[i], (_lPhi)[i], (et1 * cosh((_lEta)[i])));
  lep2.SetPtEtaPhiE(et2, (_lEta)[j], (_lPhi)[j], (et2 * cosh((_lEta)[j])));

  themass = (lep1+lep2).Mag() ;
  theptll = (lep1+lep2).Pt() ;
  thedphill = fabs(acos(cos(lep1.Phi()- lep2.Phi())));
  thephill =  (lep1+lep2).Phi();
  thepzll = (lep1+lep2).Pz() ;
  theyll= (lep1+lep2).Rapidity();

  Float_t leppz = ((_lpdgId)[i]>0)?   lep1.Pz(): lep2.Pz();
  Float_t posipz = ((_lpdgId)[i]<0)?  lep1.Pz(): lep2.Pz();
  Float_t lepenergy = ((_lpdgId)[i]>0)? lep1.E(): lep2.E();
  Float_t posienergy = ((_lpdgId)[i]<0)? lep1.E(): lep2.E();

  thecosthll = thepzll/fabs(thepzll)*2/themass/sqrt(themass*themass+theptll*theptll)*(leppz*posienergy-posipz*lepenergy) ;


}



void
JMEAnalyzer::CalcDileptonInfoGen(const int& i, const int& j, Float_t & themass, Float_t & theptll, Float_t & thepzll,  Float_t & theyll,Float_t & thephill,  Float_t & thedphill, Float_t & thecosthll){

  TLorentzVector lep1;
  TLorentzVector lep2;

  Float_t et1 = (_lgenPt)[i];
  Float_t et2 = (_lgenPt)[j];
  lep1.SetPtEtaPhiE(et1, (_lgenEta)[i], (_lgenPhi)[i], (et1 * cosh((_lgenEta)[i])));
  lep2.SetPtEtaPhiE(et2, (_lgenEta)[j], (_lgenPhi)[j], (et2 * cosh((_lgenEta)[j])));

  themass = (lep1+lep2).Mag() ;
  theptll = (lep1+lep2).Pt() ;
  thedphill = fabs(acos(cos(lep1.Phi()- lep2.Phi())));
  thephill = (lep1+lep2).Phi() ;
  thepzll = (lep1+lep2).Pz() ;
  theyll= (lep1+lep2).Rapidity();

  Float_t leppz = ((_lgenpdgId)[i]>0)?   lep1.Pz(): lep2.Pz();
  Float_t posipz = ((_lgenpdgId)[i]<0)?  lep1.Pz(): lep2.Pz();
  Float_t lepenergy = ((_lgenpdgId)[i]>0)? lep1.E(): lep2.E();
  Float_t posienergy = ((_lgenpdgId)[i]<0)? lep1.E(): lep2.E();

  thecosthll = thepzll/fabs(thepzll)*2/themass/sqrt(themass*themass+theptll*theptll)*(leppz*posienergy-posipz*lepenergy) ;


}








bool JMEAnalyzer::PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,  const pat::Muon *muonit, const edm::Event& iEvent ){

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::InputTag triggerBits_("TriggerResults","","HLT");
  edm::InputTag  triggerObjects_("slimmedPatTrigger");
  iEvent.getByToken(trigobjectToken_, triggerObjects);
  iEvent.getByToken(trgresultsToken_, triggerBits);
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    obj.unpackNamesAndLabels(iEvent,*triggerBits);                                                                                                                                                           
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
      string myfillabl=obj.filterLabels()[h];
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& deltaR(muonit->eta(),muonit->phi(), obj.eta(),obj.phi())<0.4 ) return true;
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& myfillabl.find("hltL1s")!=string::npos &&deltaR(muonit->eta(),muonit->phi(), obj.eta(),obj.phi())<0.5 ) return true; 
    }
  }
  return false;
}


bool JMEAnalyzer::PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,const pat::Electron *eleit, const edm::Event& iEvent ){

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::InputTag triggerBits_("TriggerResults","","HLT");
  edm::InputTag  triggerObjects_("slimmedPatTrigger");
  iEvent.getByToken(trigobjectToken_, triggerObjects);

  iEvent.getByToken(trgresultsToken_, triggerBits);

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    obj.unpackNamesAndLabels(iEvent,*triggerBits);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){

      string myfillabl=obj.filterLabels()[h];
                                                                                                                                                                    

      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& deltaR(eleit->eta(),eleit->phi(), obj.eta(),obj.phi())<0.4 ) return true;
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& myfillabl.find("L1sL1")!=string::npos &&deltaR(eleit->eta(),eleit->phi(), obj.eta(),
																				  obj.phi())<0.5 ) return true; 
      if(triggerlegstring =="hltDoubleEle8Mass8Filter" && myfillabl.find(triggerlegstring)!=string::npos )return true;

 
    }
  }
  return false;
}


bool JMEAnalyzer::PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,const pat::Photon *photonit, const edm::Event& iEvent ){

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::InputTag triggerBits_("TriggerResults","","HLT");
  edm::InputTag  triggerObjects_("slimmedPatTrigger");
  iEvent.getByToken(trigobjectToken_, triggerObjects);

  iEvent.getByToken(trgresultsToken_, triggerBits);

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    obj.unpackNamesAndLabels(iEvent,*triggerBits);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){

      string myfillabl=obj.filterLabels()[h];
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& deltaR(photonit->eta(),photonit->phi(), obj.eta(),obj.phi())<0.4 ) return true;
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& myfillabl.find("L1sL1")!=string::npos &&deltaR(photonit->eta(),photonit->phi(), 
																				  obj.eta(),obj.phi())<0.5 ) return true; 

    }
  }
  return false;
}



bool JMEAnalyzer::PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,const pat::Jet *jetit, const edm::Event& iEvent ){

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::InputTag triggerBits_("TriggerResults","","HLT");
  edm::InputTag  triggerObjects_("slimmedPatTrigger");
  iEvent.getByToken(trigobjectToken_, triggerObjects);

  iEvent.getByToken(trgresultsToken_, triggerBits);

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {                                                                               
    obj.unpackNamesAndLabels(iEvent,*triggerBits);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){

      string myfillabl=obj.filterLabels()[h];
 
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& deltaR(jetit->eta(),jetit->phi(), obj.eta(),obj.phi())<0.4 ) return true;
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& myfillabl.find("L1sL1")!=string::npos &&deltaR(jetit->eta(),jetit->phi(), 
																				  obj.eta(),obj.phi())<0.5 ) return true;
    }
  }
  return false;
}



int  JMEAnalyzer::GetRecoIdx(const reco::GenJet * genjet , vector <Float_t> recojetpt,  vector <Float_t> recojeteta, vector <Float_t> recojetphi, vector <Float_t> recojetptgen ){
  Float_t etagen = genjet->eta();
  Float_t phigen = genjet->phi();
  Float_t ptgen = genjet->pt();
  
  for(unsigned int i = 0; i<  recojeteta.size() ; i++){
    if(abs( recojetptgen[i] - ptgen)>0.1 )continue;// Requiring the two to be exactly equal might not work due to rounding issues? 
    Float_t dr = deltaR(etagen,phigen, recojeteta[i] ,  recojetphi[i] );
    if(dr>0.3) continue; 
    return i;
  }
  return -1;
}
      

bool JMEAnalyzer::IsTauCleaned(const pat::Jet *iJ){
  for(unsigned int i = 0 ; i < _tauPt.size(); i++){
    if(!_tauPassMediumID[i])continue;
    double dr = deltaR( _tauEta[i],_tauPhi[i],iJ->eta(),iJ->phi());
    if(dr<0.3) return false;
  }
  return true;
}


bool JMEAnalyzer::IsLeptonPhotonCleaned(const pat::Jet *iJ){
  for(unsigned int i = 0 ; i < _lPt.size(); i++){
    double dr = deltaR( _lEta[i],_lPhi[i],iJ->eta(),iJ->phi()); 
    if(dr<0.3) return false; 
  }
  for(unsigned int i = 0 ; i < _phPt.size(); i++){
    double dr = deltaR( _phEta[i],_phPhi[i],iJ->eta(),iJ->phi()); 
    if(dr<0.3) return false; 
  }
  return true;
}

bool JMEAnalyzer::IsLeptonPhotonCleaned(const reco::CaloJet *iJ){
  for(unsigned int i = 0 ; i < _lPt.size(); i++){
    double dr = deltaR( _lEta[i],_lPhi[i],iJ->eta(),iJ->phi()); 
    if(dr<0.3) return false; 
  }
  for(unsigned int i = 0 ; i < _phPt.size(); i++){
    double dr = deltaR( _phEta[i],_phPhi[i],iJ->eta(),iJ->phi()); 
    if(dr<0.3) return false; 
  }
  return true;
}

bool JMEAnalyzer::PassJetPreselection(const pat::Jet * iJ, double genjetpt, bool ispuppi){
  //These values should be chosen small enough that the resulting cut is always looser than the final cut that will be applied on the corrected pt after including L2L3res
  double ptcut = 10;
  double ptcut_foractivity_central =  15;
  double ptcut_foractivity_fwd =  20;
  if(Skim_=="MCJECs" &&  genjetpt<5) return false;
  else if(Skim_=="ZJetsResiduals" ){
    //There are 3 kinds of jets we care of: 
    //a) Those with large dphi(Z,j) for ptZ>10 => used for balance measurement
    //b) Those with low dphi(Z,j) for ptZ<5 => used to build the PU template
    //c) All other jets are used to assess the jet activity. For those we only care about large pt. 
    if(iJ->correctedP4("L3Absolute").Pt()<ptcut)return false;
    double dphi_jet_ll = fabs(deltaPhi(_phill, iJ->phi()));
    bool  passdphicond =  (_ptll>10&& dphi_jet_ll > 2.9 ) || ( _ptll<5 && dphi_jet_ll<2. && dphi_jet_ll>1.0);
    if( !passdphicond && iJ->correctedP4("L3Absolute").Pt()  <   ptcut_foractivity_central ) return false;
    if( !passdphicond && iJ->correctedP4("L3Absolute").Pt()  <   ptcut_foractivity_fwd && fabs(iJ->eta())>2.4) return false;
  }
  else if(Skim_=="GammaJetsResiduals") {
    //For gamma + jets, item b) above doesn't apply 
    if(iJ->correctedP4("L3Absolute").Pt()<ptcut)return false;
    double dphi_jet_gamma = fabs(deltaPhi(_phigamma, iJ->phi()));
    bool  passdphicond =  (_ptgamma>50&& dphi_jet_gamma > 2.9 );
    if( !passdphicond && iJ->correctedP4("L3Absolute").Pt()  <   ptcut_foractivity_central ) return false;
    if( !passdphicond && iJ->correctedP4("L3Absolute").Pt()  <   ptcut_foractivity_fwd && fabs(iJ->eta())>2.4) return false;
  }
  else if(  iJ->pt()< JetPtCut_) return false;
  return true;
}


bool JMEAnalyzer::PassJetPreselection(const reco::CaloJet * iJ, double genjetpt){
  if(Skim_=="MCJECs" &&  genjetpt<5) return false;
  else if(  iJ->pt()< JetPtCut_) return false;
  return true;

}
//define this as a plug-in
DEFINE_FWK_MODULE(JMEAnalyzer);



/*
  - prescale for photon triggers      
  - add jer variation?
  - update jet id for puppi
  - add AK8 CHS
-Add dz and isolation for leptons? 


*/
/* Then the residual procedure is as follows:
   a) Fit ptbalance to extract JES and JER for both data and MC => One can derive an approximate residual corr as 1+JES(MC)-JER(data) for a given ref pt bin
   b) Convert the ref pt bin into a <jetpt> bin 
   c) Fit the Residual vs <jetpt> for various eta bins
   d) Apply the residual to data and go back to a) 

*/
