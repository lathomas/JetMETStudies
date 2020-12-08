// Original Author:  Laurent Thomas
//         Created:  Fri, 26 Apr 2019 12:51:46 GMT
//
//


// system include files
#include <memory>


// user include files
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 


#include "JetMETStudies/JMEAnalyzer/interface/Tools.h"
#include "RecoJets/JetProducers/plugins/PileupJetIdProducer.h" 

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
const int  N_METFilters=16;
enum METFilterIndex{
  idx_Flag_goodVertices,
  idx_Flag_globalTightHalo2016Filter,
  idx_Flag_globalSuperTightHalo2016Filter,
  idx_Flag_HBHENoiseFilter,
  idx_Flag_HBHENoiseIsoFilter,
  idx_Flag_EcalDeadCellTriggerPrimitiveFilter,
  idx_Flag_BadPFMuonFilter,
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
  edm::EDGetTokenT<std::vector< pat::Jet> > jetPuppiToken_;
  edm::EDGetTokenT<std::vector< pat::Jet> > jetPuppiAK8Token_;
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
  edm::EDGetTokenT<std::vector< pat::Photon> > photonToken_;

  edm::EDGetTokenT<GenParticleCollection> genpartToken_;
  edm::EDGetTokenT<GenEventInfoProduct> geninfoToken_;
  edm::EDGetTokenT<LHEEventProduct> lheEventToken_;
  edm::EDGetTokenT<LHEEventProduct> lheEventALTToken_;

  edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetAssocCHSToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetWithNuAssocCHSToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetAssocPuppiToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetWithNuAssocPuppiToken_;

  edm::EDGetTokenT<vector<PileupSummaryInfo> > puInfoToken_;

  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;
  edm::EDGetTokenT<BXVector<GlobalAlgBlk>> l1GtToken_;

  Float_t JetPtCut_;
  Float_t AK8JetPtCut_;
  Float_t ElectronPtCut_;
  string ElectronVetoWP_, ElectronTightWP_;
  Float_t MuonPtCut_;
  string RochCorrFile_;
  Float_t PhotonPtCut_;
  string PhotonTightWP_;
  Float_t PFCandPtCut_;

  Bool_t SaveTree_, IsMC_, SavePUIDVariables_, SaveAK8Jets_, DropUnmatchedJets_, DropBadJets_, SavePFinJets_, ApplyPhotonID_;
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
  
  //Variables associated to leaves of the TTree

  unsigned long _eventNb;
  unsigned long _runNb;
  unsigned long _lumiBlock;
  unsigned long _bx;

  //Nb of primary vertices
  int _n_PV;
  Float_t _LV_z;
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

  //Decision obtained rerunning the filters on top of MINIAOD
  bool PassecalBadCalibFilter_Update;
  bool PassecalLaserCorrFilter_Update;  
  bool PassEcalDeadCellBoundaryEnergyFilter_Update;
  bool PassBadChargedCandidateFilter_Update;

  //Jets 
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

  //Puppi jet
  vector<Float_t>  _puppijetEta;
  vector<Float_t>  _puppijetPhi;
  vector<Float_t>  _puppijetPt;
  vector <Float_t>  _puppijetPtGen;
  vector <Float_t>  _puppijetPtGenWithNu;

  vector<Float_t>  _puppiak8jetEta ;
  vector<Float_t>  _puppiak8jetPt ;
  vector<Float_t>  _puppiak8jet_tau1 ;
  vector<Float_t>  _puppiak8jet_tau2 ;
  vector<Float_t>  _puppiak8jet_tau3 ;


  //Leptons
  vector<Float_t>  _lEta;
  vector<Float_t>  _lPhi;
  vector<Float_t>  _lPt;
  vector<Float_t>  _lPtcorr;
  vector<Float_t>  _lPassTightID;
  vector<int> _lpdgId;
  int _nEles, _nMus;

  vector<Float_t>  _lgenEta;
  vector<Float_t>  _lgenPhi;
  vector<Float_t>  _lgenPt;
  vector<int> _lgenpdgId;
  
  Float_t _mll;
  Float_t _ptll;
  Float_t _pzll;
  Float_t _yll;
  Float_t _dphill;
  Float_t _phill;
  Float_t _costhCSll;
  Int_t _nElesll ;

  Float_t _mll_gen;
  Float_t _ptll_gen;
  Float_t _pzll_gen;
  Float_t _yll_gen;
  Float_t _dphill_gen;
  Float_t _phill_gen;
  Float_t _costhCSll_gen;
  
  Int_t _ngenElesll ;
  
  //Nb of electrons in vtx fit
  int _n_PFele_fromvtxfit;
  int _n_PFmu_fromvtxfit;

  //Photons
  vector<Float_t>  _phEta;
  vector<Float_t>  _phPhi;
  vector<Float_t>  _phPt;
  vector<Float_t>  _phPtcorr;
  
  vector<Float_t>  _phgenEta;
  vector<Float_t>  _phgenPhi;
  vector<Float_t>  _phgenPt;

  Float_t _genHT, _weight;

  //PF candidates
  vector <Float_t> _PFcand_pt;
  vector <Float_t> _PFcand_eta;
  vector <Float_t> _PFcand_phi;
  vector <int> _PFcand_pdgId;
  vector <int> _PFcand_fromPV;
  vector <Float_t> _PFcand_dz;
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


  //MET
  Float_t _met;
  Float_t _met_phi;
  Float_t _puppimet;
  Float_t _puppimet_phi;

  Float_t _rawmet;
  Float_t _rawmet_phi;
  Float_t _puppirawmet;
  Float_t _puppirawmet_phi;


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
  bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
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
  bool _l1prefire;


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

  RoccoR rc; 
  JetCorrectionUncertainty *jecUnc; 
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
  jetPuppiToken_(consumes< std::vector< pat::Jet> >(iConfig.getParameter<edm::InputTag>("JetsPuppi"))),
  jetPuppiAK8Token_(consumes< std::vector< pat::Jet> >(iConfig.getParameter<edm::InputTag>("JetsPuppiAK8"))),
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
  photonToken_(consumes< std::vector< pat::Photon> >(iConfig.getParameter<edm::InputTag>("Photons"))),
  genpartToken_(consumes<GenParticleCollection> (iConfig.getParameter<edm::InputTag>("GenParticles"))),
  geninfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("GenInfo"))),
  lheEventToken_(consumes<LHEEventProduct> ( iConfig.getParameter<InputTag>("LHELabel"))),
  lheEventALTToken_(consumes<LHEEventProduct> ( iConfig.getParameter<InputTag>("LHELabelALT"))),
  genJetAssocCHSToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("GenJetMatchCHS"))),
  genJetWithNuAssocCHSToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("GenJetWithNuMatchCHS"))),
  genJetAssocPuppiToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("GenJetMatchPuppi"))),
  genJetWithNuAssocPuppiToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("GenJetWithNuMatchPuppi"))),
  puInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PULabel"))),
  trgresultsToken_(consumes<TriggerResults>(iConfig.getParameter<edm::InputTag>("Triggers"))),
  l1GtToken_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("l1GtSrc"))),
  JetPtCut_(iConfig.getParameter<double>("JetPtCut")),
  AK8JetPtCut_(iConfig.getParameter<double>("AK8JetPtCut")),
  ElectronPtCut_(iConfig.getParameter<double>("ElectronPtCut")),
  ElectronVetoWP_(iConfig.getParameter<string>("ElectronVetoWorkingPoint")),
  ElectronTightWP_(iConfig.getParameter<string>("ElectronTightWorkingPoint")),
  MuonPtCut_(iConfig.getParameter<double>("MuonPtCut")),
  RochCorrFile_(iConfig.getParameter<string>("RochCorrFile")),
  PhotonPtCut_(iConfig.getParameter<double>("PhotonPtCut")),
  PhotonTightWP_(iConfig.getParameter<string>("PhotonTightWorkingPoint")),
  PFCandPtCut_(iConfig.getParameter<double>("PFCandPtCut")),
  SaveTree_(iConfig.getParameter<bool>("SaveTree")), 
  IsMC_(iConfig.getParameter<bool>("IsMC")),
  SavePUIDVariables_(iConfig.getParameter<bool>("SavePUIDVariables")),
  SaveAK8Jets_(iConfig.getParameter<bool>("SaveAK8Jets")),
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

  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);

  if(JetCorParColl.isValid()){
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  jecUnc = new JetCorrectionUncertainty(JetCorPar);
  }
  else jecUnc = 0;

}


void JMEAnalyzer::endRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{

  delete jecUnc;
}

// ------------ method called for each event  ------------
void
JMEAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
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
	if(TrigPath.Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v"))HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 =true;
	if(TrigPath.Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"))HLT_PFMETNoMu120_PFMHTNoMu120_IDTight =true;
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
      }
    }
  }

  //Prefiring, see: https://github.com/nsmith-/PrefireAnalysis/#usage
  edm::Handle<BXVector<GlobalAlgBlk>> l1GtHandle;
  iEvent.getByToken(l1GtToken_, l1GtHandle);
  _l1prefire= false;
  if(!IsMC_)_l1prefire = l1GtHandle->begin(-1)->getFinalOR();



  
  //Vertices
  edm::Handle<std::vector<Vertex> > theVertices;
  iEvent.getByToken(verticesToken_,theVertices) ;
  _n_PV = theVertices->size();

  const Vertex* LVtx = &((*theVertices)[0]);
  _LV_z = LVtx->z();
  
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


  
  //MET filters are stored in TriggerResults::RECO or TriggerResults::PAT . Should take the latter if it exists
  edm::Handle<TriggerResults> METFilterResults;
  iEvent.getByToken(metfilterspatToken_, METFilterResults);
  if(!(METFilterResults.isValid())) iEvent.getByToken(metfiltersrecoToken_, METFilterResults);
  
  Flag_goodVertices= GetMETFilterDecision(iEvent,METFilterResults,"Flag_goodVertices");
  Flag_globalTightHalo2016Filter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_globalTightHalo2016Filter");
  Flag_globalSuperTightHalo2016Filter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_globalSuperTightHalo2016Filter");
  Flag_HBHENoiseFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_HBHENoiseFilter");
  Flag_HBHENoiseIsoFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_HBHENoiseIsoFilter");
  Flag_EcalDeadCellTriggerPrimitiveFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_EcalDeadCellTriggerPrimitiveFilter");
  Flag_BadPFMuonFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_BadPFMuonFilter");
  Flag_BadChargedCandidateFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_BadChargedCandidateFilter");
  Flag_eeBadScFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_eeBadScFilter");
  Flag_ecalBadCalibFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_ecalBadCalibFilter");
  Flag_EcalDeadCellBoundaryEnergyFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_EcalDeadCellBoundaryEnergyFilter");
  Flag_ecalLaserCorrFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_ecalLaserCorrFilter");
  

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



  edm::Handle< std::vector<pat::Electron> > thePatElectrons;
  iEvent.getByToken(electronToken_,thePatElectrons);
  for( std::vector<pat::Electron>::const_iterator electron = (*thePatElectrons).begin(); electron != (*thePatElectrons).end(); electron++ ) {
    if((&*electron)->pt() <5) continue; //Loose cut  on uncorrected pt 
    //Implementing smearing/scaling EGM corrections
    //See here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2#Applying_the_Energy_Scale_and_sm and here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaUL2016To2018
    double ptelecorr = (&*electron)->pt();
    if((&*electron)->hasUserFloat("ecalEnergyPostCorr") ){
         ptelecorr = ptelecorr * (&*electron)->userFloat("ecalTrkEnergyPostCorr") /  (&*electron)->energy() ;
	 //   cout << "Electron energy, Precorr, PostCorr, pt*cosh(eta) " << (&*electron)->energy() <<", "<<(&*electron)->userFloat("ecalEnergyPreCorr") <<", "<< (&*electron)->userFloat("ecalEnergyPostCorr") << ", "<<(&*electron)->pt() *cosh((&*electron)->eta()) <<endl; 
    }
    
    bool passvetoid = (&*electron)->electronID(ElectronVetoWP_)&& (&*electron)->pt()>10;
    if(!passvetoid) continue;
    _nEles++;
    if((&*electron)->pt()<ElectronPtCut_)continue;    
    _lEta.push_back((&*electron)->eta());
    _lPhi.push_back((&*electron)->phi());
    _lPt.push_back((&*electron)->pt());
    _lPtcorr.push_back(ptelecorr );
    _lpdgId.push_back(-11*(&*electron)->charge());
    _lPassTightID.push_back( (&*electron)->electronID(ElectronTightWP_) );
    
  }

  edm::Handle< std::vector<pat::Muon> > thePatMuons;
  iEvent.getByToken(muonToken_,thePatMuons);
  for( std::vector<pat::Muon>::const_iterator muon = (*thePatMuons).begin(); muon != (*thePatMuons).end(); muon++ ) {
    if((&*muon)->pt() <5) continue; //Loose cut  on uncorrected pt 

    //Rochester corrections: https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon#Rochester_Correction
    //https://indico.cern.ch/event/926898/contributions/3897122/attachments/2052816/3441285/roccor.pdf
    double ptmuoncorr= (&*muon)->pt();

    if( !IsMC_) ptmuoncorr *= rc.kScaleDT( (&*muon)->charge(),  (&*muon)->pt(), (&*muon)->eta(),(&*muon)->phi());
    else{
      if( (&*muon)->genLepton() !=0)  ptmuoncorr *= rc.kSpreadMC( (&*muon)->charge(),  (&*muon)->pt(), (&*muon)->eta(),(&*muon)->phi(), (&*muon)->genLepton()->pt() );
      else if(! ((&*muon)->innerTrack()).isNull()) ptmuoncorr *= rc.kSmearMC( (&*muon)->charge(),  (&*muon)->pt(), (&*muon)->eta(),(&*muon)->phi(),  (&*muon)->innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());
    }

    bool passvetoid=  (&*muon)->passed(reco::Muon::CutBasedIdLoose)&& (&*muon)->passed(reco::Muon::PFIsoVeryLoose)&&(&*muon)->pt()>10 ;  
    if(!passvetoid) continue;
    _nMus++;
    if((&*muon)->pt()<MuonPtCut_)continue;
    _lEta.push_back((&*muon)->eta());
    _lPhi.push_back((&*muon)->phi());
    _lPt.push_back((&*muon)->pt());
    _lPtcorr.push_back( ptmuoncorr );
    _lpdgId.push_back(-13*(&*muon)->charge());
    _lPassTightID.push_back(  (&*muon)->passed(reco::Muon::CutBasedIdMediumPrompt )&& (&*muon)->passed(reco::Muon::PFIsoTight ) );
  }

  edm::Handle< std::vector<pat::Photon> > thePatPhotons;
  iEvent.getByToken(photonToken_,thePatPhotons);
  for( std::vector<pat::Photon>::const_iterator photon = (*thePatPhotons).begin(); photon != (*thePatPhotons).end(); photon++ ) {
    if((&*photon)->pt() <10) continue; //Loose cut  on uncorrected pt 
    double ptphotoncorr = (&*photon)->pt(); //This is (possibly) the smeared/scaled pt for data ! 
    //Implementing smearing/scaling EGM corrections
    //See here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2#Applying_the_Energy_Scale_and_sm and here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaUL2016To2018
    //For |eta|>2.5 the corrections are not adapted and one should instead pick the uncorrected value:
    if((&*photon)->hasUserFloat("ecalEnergyPostCorr") && fabs((&*photon)->eta())<2.5){
         ptphotoncorr = ptphotoncorr * (&*photon)->userFloat("ecalEnergyPostCorr") /  (&*photon)->energy() ;
    }
        
    if(ptphotoncorr <PhotonPtCut_)continue;
    bool passtightid = (&*photon)->photonID(PhotonTightWP_) && (&*photon)->passElectronVeto()&& !((&*photon)->hasPixelSeed()  ) &&fabs((&*photon)->eta())<1.4442&& (&*photon)->r9()>0.9 ; 
    if(!passtightid&& ApplyPhotonID_) continue;
    _phEta.push_back((&*photon)->eta());
    _phPhi.push_back((&*photon)->phi());
    _phPt.push_back( (&*photon)->pt());
    _phPtcorr.push_back( ptphotoncorr);
    
  }
  
  
  //Compute dilepton variables
  Float_t mll(0),ptll(0),pzll(0),yll(0),phill(0),dphill(0),costhCSll(0);
  for(unsigned int i = 0; i < _lPt.size(); i++){
    if(_lPt.size() !=2) continue;
    if(_lPt[i]<20) continue;
      if(!_lPassTightID[i])  continue;
      if(fabs(_lpdgId[i]) !=11 && fabs(_lpdgId[i])!=13 ) continue;
      for(unsigned int j = 0; j < i; j++){
        if(_lPt[j]<20) continue;
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
    

  if(!PassSkim()) return;

    
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

      const reco::GenJet * updatedgenjet =0;
      const reco::GenJet * updatedgenjetwithnu =0; 
      if(genJetMatch.isValid() && genJetWithNuMatch.isValid() && useupdategenjets){
	updatedgenjet = ( (*genJetMatch)[jetRef].isNonnull() && (*genJetMatch)[jetRef].isAvailable()) ? &*(*genJetMatch)[jetRef] : 0;
	updatedgenjetwithnu = ( (*genJetWithNuMatch)[jetRef].isNonnull() && (*genJetWithNuMatch)[jetRef].isAvailable()) ? &*(*genJetWithNuMatch)[jetRef] : 0;
      }
 
      const reco::GenJet * genjet = useupdategenjets?updatedgenjet: (&*jet) ->genJet()  ;
      
      if((&*jet)->pt() >leadjetpt) leadjetpt = (&*jet)->pt();
      if((&*jet)->pt()<JetPtCut_) continue;
      bool passid = PassJetID(  (&*jet) ,"2018");
      if( DropBadJets_ && !passid  ) continue;//Drop bad jets (mostly leptons).

      if( genjet ==0   && DropUnmatchedJets_ && (&*jet)->pt()<50 ) continue;//Drop genunmatched jets (mostly PU). Keep those with pt>50 as these probably require special attention.
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
      //Accessing the default PU ID stored in MINIAOD https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
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

      //Quark Gluon likelihood  https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuarkGluonLikelihood
      iEvent.getByToken(qgLToken_, quarkgluonlikelihood);
      if(quarkgluonlikelihood.isValid() )_jetQuarkGluonLikelihood.push_back( (*quarkgluonlikelihood)[jetRef] );
      else _jetQuarkGluonLikelihood.push_back( -1.);
      

      _jetRawPt.push_back( (&*jet)->correctedP4("Uncorrected").Pt() );
      _jetPtNoL2L3Res.push_back( (&*jet)->correctedP4("L3Absolute") .Pt() ); 
      _jet_corrjecs.push_back((&*jet)->pt() / (&*jet)->correctedP4("Uncorrected").Pt() );
      //Accessing uncertainties
      //      if(jecUnc ==0) _jetJECuncty.push_back( -99.);
      //else{
	jecUnc->setJetEta((&*jet)->eta());
	jecUnc->setJetPt((&*jet)->pt());
	_jetJECuncty.push_back( jecUnc->getUncertainty(true) );
	//}
      Float_t jetptgen(-99.), jetetagen(-99.),jetphigen(-99.);
      Float_t jetptgenwithnu(-99.);
      if( genjet !=0 ){
	jetptgen= genjet->pt() ;
	jetetagen= genjet->eta() ;
	jetphigen= genjet->phi() ;
      }
      if(updatedgenjetwithnu !=0) jetptgenwithnu = updatedgenjetwithnu->pt() ;
      _jetPtGen.push_back(jetptgen);
      _jetEtaGen.push_back(jetetagen);
      _jetPhiGen.push_back(jetphigen);
      _jetPtGenWithNu.push_back(jetptgenwithnu);
    
      for (unsigned i = 0; i < jet->numberOfSourceCandidatePtrs(); ++i) {
	//const pat::PackedCandidate * icand = jet->sourceCandidatePtr(i).get();
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
      if((&*jet)->pt()<JetPtCut_) continue;
      _puppijetEta.push_back((&*jet)->eta());
      _puppijetPhi.push_back((&*jet)->phi());
      _puppijetPt.push_back((&*jet)->pt());
      Float_t jetptgen(-99.);//, jetetagen(-99.),jetphigen(-99.);
      Float_t jetptgenwithnu(-99.);
      
      if( genjet !=0 ){
        jetptgen= genjet->pt() ;
	//        jetetagen= genjet->eta() ;
        //jetphigen= genjet->phi() ;
      }
      if(updatedgenjetwithnu !=0) jetptgenwithnu = updatedgenjetwithnu->pt() ;
      _puppijetPtGen.push_back(jetptgen);
      _puppijetPtGenWithNu.push_back(jetptgenwithnu);
    }
  }
  

  //AK8                                                                                                                                                                                                                     
  edm::Handle< std::vector< pat::Jet> > thePuppiAK8Jets;
  iEvent.getByToken(jetPuppiAK8Token_,thePuppiAK8Jets );
  if(thePuppiAK8Jets.isValid()){
    for( std::vector<pat::Jet>::const_iterator jet = (*thePuppiAK8Jets).begin(); jet != (*thePuppiAK8Jets).end(); jet++ ) {
      if((&*jet)->pt()<AK8JetPtCut_) continue;
      if(fabs((&*jet)->eta())>2.4) continue;
      _puppiak8jetEta.push_back((&*jet)->eta());
      _puppiak8jetPt.push_back((&*jet)->pt());
      _puppiak8jet_tau1.push_back((&*jet)->userFloat("NjettinessAK8Puppi:tau1"));
      _puppiak8jet_tau2.push_back((&*jet)->userFloat("NjettinessAK8Puppi:tau2"));
      _puppiak8jet_tau3.push_back((&*jet)->userFloat("NjettinessAK8Puppi:tau3"));
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
  //  for(pat::PackedCandidateCollection::const_reverse_iterator p = pfcands->rbegin() ; p != pfcands->rend() ; p++ ) {
  for( std::vector<pat::PackedCandidate>::const_iterator p = (*pfcands).begin(); p != (*pfcands).end(); p++ ) {

    
    int idxrefvtx= (p->vertexRef().isNonnull()) ?  p->vertexRef().key() : -1;
    //cout << "idxrefvtx, vtx size  " << idxrefvtx << ", " << _n_PV <<endl;
    edm::RefToBase<pat::PackedCandidate> pfcandRef (edm::Ref<pat::PackedCandidateCollection>( pfcands , p  - pfcands->begin()));

    //    edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>( thePuppiJets , jet -thePuppiJets->begin()));
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
    _PFcand_hcalFraction.push_back(p->hcalFraction());
    _PFcand_PVfitidx.push_back(idxrefvtx);
    _PFcand_puppiweight.push_back((*puppiweights)[pfcandRef]);
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

  outputTree->Branch("_eventNb",   &_eventNb,   "_eventNb/l");
  outputTree->Branch("_runNb",     &_runNb,     "_runNb/l");
  outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");
  outputTree->Branch("_bx", &_bx, "_bx/l");
  outputTree->Branch("_n_PV", &_n_PV, "_n_PV/I");
  outputTree->Branch("_LV_z", &_LV_z, "_LV_z/f");
  outputTree->Branch("_rho", &_rho, "_rho/f");
  outputTree->Branch("_rhoNC", &_rhoNC, "_rhoNC/f");
  

  outputTree->Branch("Flag_goodVertices",&Flag_goodVertices,"Flag_goodVertices/O");
  outputTree->Branch("Flag_globalTightHalo2016Filter",&Flag_globalTightHalo2016Filter,"Flag_globalTightHalo2016Filter/O");
  outputTree->Branch("Flag_globalSuperTightHalo2016Filter",&Flag_globalSuperTightHalo2016Filter,"Flag_globalSuperTightHalo2016Filter/O");
  outputTree->Branch("Flag_HBHENoiseFilter",&Flag_HBHENoiseFilter,"Flag_HBHENoiseFilter/O");
  outputTree->Branch("Flag_HBHENoiseIsoFilter",&Flag_HBHENoiseIsoFilter,"Flag_HBHENoiseIsoFilter/O");
  outputTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter",&Flag_EcalDeadCellTriggerPrimitiveFilter,"Flag_EcalDeadCellTriggerPrimitiveFilter/O");
  outputTree->Branch("Flag_BadPFMuonFilter",&Flag_BadPFMuonFilter,"Flag_BadPFMuonFilter/O");
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
  outputTree->Branch("_jet_CHEF",&_jet_CHEF);
  outputTree->Branch("_jet_NHEF",&_jet_NHEF);
  outputTree->Branch("_jet_NEEF",&_jet_NEEF);
  outputTree->Branch("_jet_CEEF",&_jet_CEEF);
  outputTree->Branch("_jet_MUEF",&_jet_MUEF);
  outputTree->Branch("_jet_CHM",&_jet_CHM);
  outputTree->Branch("_jet_NHM",&_jet_NHM);
  outputTree->Branch("_jet_PHM",&_jet_PHM);
  outputTree->Branch("_jet_NM",&_jet_NM);
  outputTree->Branch("_jetArea",&_jetArea);
  outputTree->Branch("_jetPassID",&_jetPassID);

  outputTree->Branch("_jetPtGen",&_jetPtGen);
  outputTree->Branch("_jetEtaGen",&_jetEtaGen);
  outputTree->Branch("_jetPhiGen",&_jetPhiGen);
  outputTree->Branch("_jetPtGenWithNu",&_jetPtGenWithNu);
  outputTree->Branch("_jetJECuncty",&_jetJECuncty);
  outputTree->Branch("_jetPUMVA",&_jetPUMVA);
  outputTree->Branch("_jetPUMVAUpdate",&_jetPUMVAUpdate);
  outputTree->Branch("_jetPUMVAUpdate2017",&_jetPUMVAUpdate2017);
  outputTree->Branch("_jetPUMVAUpdate2018",&_jetPUMVAUpdate2018);
  outputTree->Branch("_jetPtNoL2L3Res",&_jetPtNoL2L3Res);
  outputTree->Branch("_jet_corrjecs",&_jet_corrjecs);
  outputTree->Branch("_jethadronFlavour",&_jethadronFlavour);
  outputTree->Branch("_jetpartonFlavour",&_jetpartonFlavour);
  outputTree->Branch("_jetDeepJet_b",&_jetDeepJet_b);
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
  outputTree->Branch("_puppijetPtGen",&_puppijetPtGen);
  outputTree->Branch("_puppijetPtGenWithNu",&_puppijetPtGenWithNu);

  if(SaveAK8Jets_){
  outputTree->Branch("_puppiak8jetEta",&_puppiak8jetEta);
  outputTree->Branch("_puppiak8jetPt",&_puppiak8jetPt);
  outputTree->Branch("_puppiak8jet_tau1",&_puppiak8jet_tau1);
  outputTree->Branch("_puppiak8jet_tau2",&_puppiak8jet_tau2);
  outputTree->Branch("_puppiak8jet_tau3",&_puppiak8jet_tau3);
  }

  
  outputTree->Branch("_lEta",&_lEta);
  outputTree->Branch("_lPhi",&_lPhi);
  outputTree->Branch("_lPt",&_lPt);
  outputTree->Branch("_lPtcorr",&_lPtcorr);
  outputTree->Branch("_lpdgId",&_lpdgId);
  outputTree->Branch("_lPassTightID",&_lPassTightID);
  outputTree->Branch("_nEles", &_nEles, "_nEles/I");
  outputTree->Branch("_nMus", &_nMus, "_nMus/I");
  
  if(Skim_=="ZToEEorMuMu" || Skim_=="Dilepton" || Skim_=="DileptonInfo" ){
    outputTree->Branch("_mll", &_mll, "_mll/f");
    outputTree->Branch("_ptll", &_ptll, "_ptll/f");
    outputTree->Branch("_pzll", &_pzll, "_pzll/f");
    outputTree->Branch("_yll", &_yll, "_yll/f");
    outputTree->Branch("_dphill", &_dphill, "_dphill/f");
    outputTree->Branch("_phill", &_phill, "_phill/f");
    outputTree->Branch("_costhCSll", &_costhCSll, "_costhCSll/f");
    outputTree->Branch("_nElesll",&_nElesll,"_nElesll/I");
    
    outputTree->Branch("_mll_gen", &_mll_gen, "_mll_gen/f");
    outputTree->Branch("_ptll_gen", &_ptll_gen, "_ptll_gen/f");
    outputTree->Branch("_pzll_gen", &_pzll_gen, "_pzll_gen/f");
    outputTree->Branch("_yll_gen", &_yll_gen, "_yll_gen/f");
    outputTree->Branch("_dphill_gen", &_dphill_gen, "_dphill_gen/f");
    outputTree->Branch("_phill_gen", &_phill_gen, "_phill_gen/f");
    outputTree->Branch("_costhCSll_gen", &_costhCSll_gen, "_costhCSll_gen/f");
    outputTree->Branch("_ngenElesll",&_ngenElesll,"_ngenElesll/I");
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
    outputTree->Branch("_genHT",&_genHT,"_genHT/f");
    outputTree->Branch("_weight",&_weight,"_weight/f");
    

  }
  
  if(PhotonPtCut_<1000){
  outputTree->Branch("_phEta",&_phEta);
  outputTree->Branch("_phPhi",&_phPhi);
  outputTree->Branch("_phPt",&_phPt);
  outputTree->Branch("_phPtcorr",&_phPtcorr);
  }
  if(PFCandPtCut_<1000){
  outputTree->Branch("_PFcand_pt",&_PFcand_pt);
  outputTree->Branch("_PFcand_eta",&_PFcand_eta);
  outputTree->Branch("_PFcand_phi",&_PFcand_phi);
  outputTree->Branch("_PFcand_pdgId",&_PFcand_pdgId);
  outputTree->Branch("_PFcand_fromPV",&_PFcand_fromPV);
  outputTree->Branch("_PFcand_dz",&_PFcand_dz);
  outputTree->Branch("_PFcand_hcalFraction",&_PFcand_hcalFraction);
  outputTree->Branch("_PFcand_PVfitidx",&_PFcand_PVfitidx);
  outputTree->Branch("_PFcand_puppiweight",&_PFcand_puppiweight);
  outputTree->Branch("_n_CH_fromvtxfit",&_n_CH_fromvtxfit,"_n_CH_fromvtxfit[6]/I");
  outputTree->Branch("_HT_CH_fromvtxfit", &_HT_CH_fromvtxfit, "_HT_CH_fromvtxfit[6]/f");
  }
  
  if(Skim_=="VtxInfo" ){
  outputTree->Branch("_METCH_PV",&_METCH_PV);
  outputTree->Branch("_METPhiCH_PV",&_METPhiCH_PV);
  outputTree->Branch("_SumPT2CH_PV",&_SumPT2CH_PV);
  outputTree->Branch("_DztoLV_PV",&_DztoLV_PV);
  }

  if(IsMC_){
  outputTree->Branch("_genmet", &_genmet, "_genmet/f");
  outputTree->Branch("_genmet_phi", &_genmet_phi, "_genmet_phi/f");
  outputTree->Branch("trueNVtx", &trueNVtx,"trueNVtx/I");
  }
  outputTree->Branch("_met", &_met, "_met/f");
  outputTree->Branch("_met_phi", &_met_phi, "_met_phi/f");
  outputTree->Branch("_puppimet", &_puppimet, "_puppimet/f");
  outputTree->Branch("_puppimet_phi", &_puppimet_phi, "_puppimet_phi/f");

  outputTree->Branch("_rawmet", &_rawmet, "_rawmet/f");
  outputTree->Branch("_rawmet_phi", &_rawmet_phi, "_rawmet_phi/f");
  outputTree->Branch("_puppirawmet", &_puppirawmet, "_puppirawmet/f");
  outputTree->Branch("_puppirawmet_phi", &_puppirawmet_phi, "_puppirawmet_phi/f");
  

  outputTree->Branch("HLT_Photon110EB_TightID_TightIso",&HLT_Photon110EB_TightID_TightIso,"HLT_Photon110EB_TightID_TightIso/O");
  outputTree->Branch("HLT_Photon165_R9Id90_HE10_IsoM",&HLT_Photon165_R9Id90_HE10_IsoM,"HLT_Photon165_R9Id90_HE10_IsoM/O");
  outputTree->Branch("HLT_Photon120_R9Id90_HE10_IsoM",&HLT_Photon120_R9Id90_HE10_IsoM,"HLT_Photon120_R9Id90_HE10_IsoM/O");
  outputTree->Branch("HLT_Photon90_R9Id90_HE10_IsoM",&HLT_Photon90_R9Id90_HE10_IsoM,"HLT_Photon90_R9Id90_HE10_IsoM/O");
  outputTree->Branch("HLT_Photon75_R9Id90_HE10_IsoM",&HLT_Photon75_R9Id90_HE10_IsoM,"HLT_Photon75_R9Id90_HE10_IsoM/O");
  outputTree->Branch("HLT_Photon50_R9Id90_HE10_IsoM",&HLT_Photon50_R9Id90_HE10_IsoM,"HLT_Photon50_R9Id90_HE10_IsoM/O");
  outputTree->Branch("HLT_Photon200",&HLT_Photon200,"HLT_Photon200/O");
  outputTree->Branch("HLT_Photon175",&HLT_Photon175,"HLT_Photon175/O");
  outputTree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60/O");
  outputTree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight,"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight/O");
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

  if(SavePFinJets_){
    jetPFTree->Branch("_Jet_Pt", &_Jet_Pt,"_Jet_Pt/f");
    jetPFTree->Branch("_Jet_Eta", &_Jet_Eta,"_Jet_Eta/f");
    jetPFTree->Branch("_Jet_Phi", &_Jet_Phi,"_Jet_Phi/f");
    jetPFTree->Branch("_Jet_PtGen", &_Jet_PtGen,"_Jet_PtGen/f");
    jetPFTree->Branch("_Jet_EtaGen", &_Jet_EtaGen,"_Jet_EtaGen/f");
    jetPFTree->Branch("_Jet_PhiGen", &_Jet_PhiGen,"_Jet_PhiGen/f");

    jetPFTree->Branch("_PFcand_pt",&_Jet_PFcand_pt);
    jetPFTree->Branch("_PFcand_eta",&_Jet_PFcand_eta);
    jetPFTree->Branch("_PFcand_phi",&_Jet_PFcand_phi);
    jetPFTree->Branch("_PFcand_pdgId",&_Jet_PFcand_pdgId);
    jetPFTree->Branch("_PFcand_fromPV",&_Jet_PFcand_fromPV);
    jetPFTree->Branch("_PFcand_dz",&_Jet_PFcand_dz);
    jetPFTree->Branch("_PFcand_dzError",&_Jet_PFcand_dzError);
  }

  
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
  _puppijetPtGen.clear();
  _puppijetPtGenWithNu.clear();

  _puppiak8jetEta.clear();
  _puppiak8jetPt.clear();
  _puppiak8jet_tau1.clear();
  _puppiak8jet_tau2.clear();
  _puppiak8jet_tau3.clear();


  _lEta.clear();
  _lPhi.clear();
  _lPt.clear();
  _lPtcorr.clear();
  _lpdgId.clear();
  _lPassTightID.clear();
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

  _phgenEta.clear();
  _phgenPhi.clear();
  _phgenPt.clear();


  _PFcand_pt.clear();
  _PFcand_eta.clear();
  _PFcand_phi.clear();
  _PFcand_pdgId.clear();
  _PFcand_fromPV.clear();
  _PFcand_dz.clear();
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




  HLT_Photon110EB_TightID_TightIso=false;
  HLT_Photon165_R9Id90_HE10_IsoM=false;
  HLT_Photon120_R9Id90_HE10_IsoM=false;
  HLT_Photon90_R9Id90_HE10_IsoM=false;
  HLT_Photon75_R9Id90_HE10_IsoM=false;
  HLT_Photon50_R9Id90_HE10_IsoM=false;
  HLT_Photon200=false;
  HLT_Photon175=false;
  HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60=false;
  HLT_PFMETNoMu120_PFMHTNoMu120_IDTight=false;
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
}

bool JMEAnalyzer::PassSkim(){
  
  if(Skim_=="ZToEEorMuMu" || Skim_=="Dilepton"){


    if(_mll>20&&Skim_=="Dilepton")  return true;
    if(_mll>70&&_mll<110) return true;
    return false;
  }
  else if(Skim_=="Photon"){
    
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
  else if(Skim_=="HighHT") return (HLT_PFHT1050 || HLT_PFHT900 || HLT_PFJet500 || HLT_AK8PFJet500) ; 
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




//define this as a plug-in
DEFINE_FWK_MODULE(JMEAnalyzer);
