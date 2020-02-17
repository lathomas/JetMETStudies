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

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 


#include "JetMETStudies/JMEAnalyzer/interface/Tools.h"
 

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"

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


class JMEAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit JMEAnalyzer(const edm::ParameterSet&);
      ~JMEAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual bool GetMETFilterDecision(const edm::Event& iEvent, edm::Handle<TriggerResults> METFilterResults, TString studiedfilter);
  virtual bool GetIdxFilterDecision(int it);
  virtual TString GetIdxFilterName(int it);
  virtual void InitandClearStuff();
  
 
  // ----------member data ---------------------------
  edm::EDGetTokenT<TriggerResults> metfilterspatToken_; 
  edm::EDGetTokenT<TriggerResults> metfiltersrecoToken_; 
  edm::EDGetTokenT<bool> ecalBadCalibFilterUpdateToken_;
  edm::EDGetTokenT<bool> ecalLaserCorrFilterUpdateToken_;  
  edm::EDGetTokenT<bool> ecalDeadCellBoundaryEnergyFilterUpdateToken_;
  edm::EDGetTokenT<bool> badChargedCandidateFilterUpdateToken_;
  edm::EDGetTokenT<std::vector<Vertex> > verticesToken_; 

  edm::EDGetTokenT<std::vector< pat::Jet> > jetToken_;
  edm::EDGetTokenT<std::vector< pat::Jet> > jetPuppiToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_;
  
  edm::EDGetTokenT<std::vector< pat::MET> > metToken_;
  edm::EDGetTokenT<std::vector< pat::MET> > puppimetToken_;
  edm::EDGetTokenT<std::vector< pat::Electron> > electronToken_;
  edm::EDGetTokenT<std::vector< pat::Muon> > muonToken_;
  edm::EDGetTokenT<std::vector< pat::Photon> > photonToken_;

  edm::EDGetTokenT<GenParticleCollection> genpartToken_;
  edm::EDGetTokenT<vector<PileupSummaryInfo> > puInfoToken_;

  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;

  Float_t JetPtCut_, ElectronPtCut_, MuonPtCut_, PhotonPtCut_, PFCandPtCut_;
  Bool_t SaveTree_, IsMC_;

  //Some histos to be saved for simple checks 
  TH1F *h_PFMet, *h_PuppiMet, *h_nvtx;
  //The output TTree
  TTree* outputTree;

  //Variables associated to leaves of the TTree

  unsigned long _eventNb;
  unsigned long _runNb;
  unsigned long _lumiBlock;
  unsigned long _bx;

  //Nb of primary vertices
  int _n_PV;
  int trueNVtx;
  
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
  vector<double>  _jetEta;
  vector<double>  _jetPhi;
  vector<double>  _jetPt;
  vector<double>  _jetRawPt;
  vector<double>  _jet_CHEF;
  vector<double>  _jet_NHEF;
  vector<double>  _jet_NEEF;
  vector<double>  _jet_CEEF;
  vector<double>  _jet_MUEF;
  vector <int>  _jet_CHM;
  vector <int>  _jet_NHM;
  vector <int>  _jet_PHM;
  vector <int>  _jet_NM;
  vector <bool> _jetPassID;
  vector <double>  _jetPtGen;
  vector <double>  _jetEtaGen;
  vector <double>  _jetPhiGen;
  vector<double>  _jetJECuncty;
  vector<double>  _jetPUMVA; 
  vector<double>  _jetPtNoL2L3Res;
  vector<double> _jet_corrjecs;
  vector<int> _jethadronFlavour;
  vector<int> _jetpartonFlavour;

  //Leptons
  vector<double>  _lEta;
  vector<double>  _lPhi;
  vector<double>  _lPt;
  vector<int> _lpdgId;

  vector<double>  _lgenEta;
  vector<double>  _lgenPhi;
  vector<double>  _lgenPt;
  vector<int> _lgenpdgId;


  //Photons
  vector<double>  _phEta;
  vector<double>  _phPhi;
  vector<double>  _phPt;

  vector<double>  _phgenEta;
  vector<double>  _phgenPhi;
  vector<double>  _phgenPt;


  //PF candidates
  vector <double> _PFcand_pt;
  vector <double> _PFcand_eta;
  vector <double> _PFcand_phi;
  vector <int> _PFcand_pdgId;
  vector <int> _PFcand_fromPV;


  //MET
  double _met;
  double _met_phi;
  double _puppimet;
  double _puppimet_phi;

  double _genmet;
  double _genmet_phi;

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
  jetToken_(consumes< std::vector< pat::Jet> >(iConfig.getParameter<edm::InputTag>("Jets"))),
  jetPuppiToken_(consumes< std::vector< pat::Jet> >(iConfig.getParameter<edm::InputTag>("JetsPuppi"))),
  pfcandsToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("PFCandCollection"))),
  metToken_(consumes<std::vector<pat::MET> > (iConfig.getParameter<edm::InputTag>("PFMet"))),
  puppimetToken_(consumes<std::vector<pat::MET> > (iConfig.getParameter<edm::InputTag>("PuppiMet"))),
  electronToken_(consumes< std::vector< pat::Electron> >(iConfig.getParameter<edm::InputTag>("Electrons"))),
  muonToken_(consumes< std::vector< pat::Muon> >(iConfig.getParameter<edm::InputTag>("Muons"))),
  photonToken_(consumes< std::vector< pat::Photon> >(iConfig.getParameter<edm::InputTag>("Photons"))),
  genpartToken_(consumes<GenParticleCollection> (iConfig.getParameter<edm::InputTag>("GenParticles"))),
  puInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PULabel"))),
  trgresultsToken_(consumes<TriggerResults>(iConfig.getParameter<edm::InputTag>("Triggers"))),
  JetPtCut_(iConfig.getParameter<double>("JetPtCut")),
  ElectronPtCut_(iConfig.getParameter<double>("ElectronPtCut")),
  MuonPtCut_(iConfig.getParameter<double>("MuonPtCut")),
  PhotonPtCut_(iConfig.getParameter<double>("PhotonPtCut")),
  PFCandPtCut_(iConfig.getParameter<double>("PFCandPtCut")),
  SaveTree_(iConfig.getParameter<bool>("SaveTree"))
  //  IsMC_(iConfig.getParameter<bool>("IsMC"))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs; 
  h_nvtx  = fs->make<TH1F>("h_nvtx" , "Number of reco vertices;N_{vtx};Events"  ,    100, 0., 100.);
  h_PFMet  = fs->make<TH1F>("h_PFMet" , "Type 1 PFMET (GeV);Type 1 PFMET (GeV);Events"  ,    1000, 0., 5000.);
  h_PuppiMet  = fs->make<TH1F>("h_PuppiMet" , "PUPPI MET (GeV);PUPPI MET (GeV);Events"  ,    1000, 0., 5000.);

  outputTree = fs->make<TTree>("tree","tree");

}


JMEAnalyzer::~JMEAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
JMEAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  InitandClearStuff();
  
  
  _runNb = iEvent.id().run();
  _eventNb = iEvent.id().event();
  _lumiBlock = iEvent.luminosityBlock();
  _bx=iEvent.bunchCrossing();
  
  //Vertices
  edm::Handle<std::vector<Vertex> > theVertices;
  iEvent.getByToken(verticesToken_,theVertices) ;
  _n_PV = theVertices->size();

  
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
  else PassecalBadCalibFilter_Update = true;
  
  edm::Handle<bool> handle_PassecalLaserCorrFilter_Update ;
  iEvent.getByToken(ecalLaserCorrFilterUpdateToken_,handle_PassecalLaserCorrFilter_Update);
  if(handle_PassecalLaserCorrFilter_Update.isValid())PassecalLaserCorrFilter_Update =  (*handle_PassecalLaserCorrFilter_Update );
  else PassecalLaserCorrFilter_Update = true;
  
  edm::Handle<bool> handle_PassEcalDeadCellBoundaryEnergyFilter_Update;
  iEvent.getByToken(ecalDeadCellBoundaryEnergyFilterUpdateToken_,handle_PassEcalDeadCellBoundaryEnergyFilter_Update);
  if(handle_PassEcalDeadCellBoundaryEnergyFilter_Update.isValid())PassEcalDeadCellBoundaryEnergyFilter_Update =  (*handle_PassEcalDeadCellBoundaryEnergyFilter_Update );
  else{  std::cout <<"handle_PassEcalDeadCellBoundaryEnergyFilter_Update.isValid =false" <<endl; PassEcalDeadCellBoundaryEnergyFilter_Update = true;}

  edm::Handle<bool> handle_PassBadChargedCandidateFilter_Update;
  iEvent.getByToken(badChargedCandidateFilterUpdateToken_,handle_PassBadChargedCandidateFilter_Update);
  if(handle_PassBadChargedCandidateFilter_Update.isValid())PassBadChargedCandidateFilter_Update =  (*handle_PassBadChargedCandidateFilter_Update );
  else{  std::cout <<"handle_PassBadChargedCandidateFilter_Update.isValid =false" <<endl; PassBadChargedCandidateFilter_Update = true;}

    
  //Jets
  edm::Handle< std::vector< pat::Jet> > theJets;
  iEvent.getByToken(jetToken_,theJets );
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

  double leadjetpt (0.);
  for( std::vector<pat::Jet>::const_iterator jet = (*theJets).begin(); jet != (*theJets).end(); jet++ ) {
    if((&*jet)->pt() >leadjetpt) leadjetpt = (&*jet)->pt();
    if((&*jet)->pt()<JetPtCut_) continue;

    
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
    _jetPUMVA.push_back( (&*jet)->userFloat("pileupJetId:fullDiscriminant") );
    
    _jethadronFlavour.push_back((&*jet)->hadronFlavour());  
    _jetpartonFlavour.push_back((&*jet)->partonFlavour());   

    
    bool passid = PassJetID(  (&*jet) ,"2018");
    _jetPassID.push_back(passid);
    _jetRawPt.push_back( (&*jet)->correctedP4("Uncorrected").Pt() );
    _jetPtNoL2L3Res.push_back( (&*jet)->correctedP4("L3Absolute") .Pt() ); 
    _jet_corrjecs.push_back((&*jet)->pt() / (&*jet)->correctedP4("Uncorrected").Pt() );
    //Accessing uncertainties
    jecUnc->setJetEta((&*jet)->eta());
    jecUnc->setJetPt((&*jet)->pt());
    _jetJECuncty.push_back( jecUnc->getUncertainty(true) );
    
    double jetptgen(-99.), jetetagen(-99.),jetphigen(-99.);
    if( (&*jet) ->genJet() !=0 ){
      jetptgen= (&*jet)->genJet()->pt() ;
      jetetagen= (&*jet)->genJet()->eta() ;
      jetphigen= (&*jet)->genJet()->phi() ;
    }
    _jetPtGen.push_back(jetptgen);
    _jetEtaGen.push_back(jetetagen);
    _jetPhiGen.push_back(jetphigen);

  }
  
  
  edm::Handle< std::vector<pat::Electron> > thePatElectrons;
  iEvent.getByToken(electronToken_,thePatElectrons);
  for( std::vector<pat::Electron>::const_iterator electron = (*thePatElectrons).begin(); electron != (*thePatElectrons).end(); electron++ ) {
    if((&*electron)->pt()<ElectronPtCut_)continue;
    bool passid = (&*electron)->electronID("mvaEleID-Fall17-iso-V2-wp90");
    if(!passid) continue;
    
    _lEta.push_back((&*electron)->eta());
    _lPhi.push_back((&*electron)->phi());
    _lPt.push_back((&*electron)->pt());
    _lpdgId.push_back(-11*(&*electron)->charge());

  }

  edm::Handle< std::vector<pat::Muon> > thePatMuons;
  iEvent.getByToken(muonToken_,thePatMuons);
  for( std::vector<pat::Muon>::const_iterator muon = (*thePatMuons).begin(); muon != (*thePatMuons).end(); muon++ ) {
    if((&*muon)->pt()<MuonPtCut_)continue;
    bool passid=  (&*muon)->passed(reco::Muon::CutBasedIdMediumPrompt )&& (&*muon)->passed(reco::Muon::PFIsoTight ); 
    if(!passid) continue;
    
    _lEta.push_back((&*muon)->eta());
    _lPhi.push_back((&*muon)->phi());
    _lPt.push_back((&*muon)->pt());
    _lpdgId.push_back(-13*(&*muon)->charge());
  }

  edm::Handle< std::vector<pat::Photon> > thePatPhotons;
  iEvent.getByToken(photonToken_,thePatPhotons);
  for( std::vector<pat::Photon>::const_iterator photon = (*thePatPhotons).begin(); photon != (*thePatPhotons).end(); photon++ ) {
    if((&*photon)->pt()<PhotonPtCut_)continue;
    bool passid = (&*photon)->photonID("mvaPhoID-RunIIFall17-v1p1-wp80") && (&*photon)->passElectronVeto()&& !((&*photon)->hasPixelSeed()) ; 
    if(!passid) continue;
        
    _phEta.push_back((&*photon)->eta());
    _phPhi.push_back((&*photon)->phi());
    _phPt.push_back((&*photon)->pt());
  }
  

   
  //Type 1 PFMET
  edm::Handle< vector<pat::MET> > ThePFMET;
  iEvent.getByToken(metToken_, ThePFMET);
  const vector<pat::MET> *pfmetcol = ThePFMET.product();
  const pat::MET *pfmet;
  pfmet = &(pfmetcol->front());
  _met = pfmet->pt();
  _met_phi = pfmet->phi();

  //PUPPI MET
  edm::Handle< vector<pat::MET> > ThePUPPIMET;
  iEvent.getByToken(puppimetToken_, ThePUPPIMET);
  const vector<pat::MET> *puppimetcol = ThePUPPIMET.product();
  const pat::MET *puppimet;
  puppimet = &(puppimetcol->front());
  _puppimet = puppimet->pt();
  _puppimet_phi = puppimet->phi();
  
  //PF candidates
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfcandsToken_ ,pfcands);
  for(pat::PackedCandidateCollection::const_reverse_iterator p = pfcands->rbegin() ; p != pfcands->rend() ; p++ ) {
    if(p->pt()<PFCandPtCut_)continue;
    _PFcand_pt.push_back(p->pt());
    _PFcand_eta.push_back(p->eta());
    _PFcand_phi.push_back(p->phi());
    _PFcand_pdgId.push_back(p->pdgId());
    _PFcand_fromPV.push_back(p->fromPV(0));//See https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Packed_ParticleFlow_Candidates
  }
  


  //Gen particle info
  edm::Handle<GenParticleCollection> TheGenParticles;
  iEvent.getByToken(genpartToken_, TheGenParticles);
  TLorentzVector Gen0;
  Gen0.SetPtEtaPhiE( 0, 0, 0, 0);
  if(TheGenParticles.isValid()){
    for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
      int id = TMath::Abs(p->pdgId());
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
  }
  if (Gen0.E()!=0) {
    _genmet = Gen0.Pt();
    _genmet_phi = Gen0.Phi();
  } else {
    _genmet = 0;
    _genmet_phi = 0;
  }

  Handle<std::vector<PileupSummaryInfo> > puInfo;
  iEvent.getByToken(puInfoToken_, puInfo);
  if(puInfo.isValid()){
  vector<PileupSummaryInfo>::const_iterator pvi;
  for (pvi = puInfo->begin(); pvi != puInfo->end(); ++pvi) {
    if (pvi->getBunchCrossing() == 0) trueNVtx = pvi->getTrueNumInteractions();
  }
  }
  else trueNVtx = -1.;


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
  
  
  //Filling trees and histos   
  if(SaveTree_)outputTree->Fill();
  h_PFMet->Fill(_met);
  h_PuppiMet->Fill(_puppimet);
  h_nvtx->Fill(_n_PV);
  
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
  outputTree->Branch("_jetPassID",&_jetPassID);

  outputTree->Branch("_jetPtGen",&_jetPtGen);
  outputTree->Branch("_jetEtaGen",&_jetEtaGen);
  outputTree->Branch("_jetPhiGen",&_jetPhiGen);
  outputTree->Branch("_jetJECuncty",&_jetJECuncty);
  outputTree->Branch("_jetPUMVA",&_jetPUMVA);
  outputTree->Branch("_jetPtNoL2L3Res",&_jetPtNoL2L3Res);
  outputTree->Branch("_jet_corrjecs",&_jet_corrjecs);
  outputTree->Branch("_jethadronFlavour",&_jethadronFlavour);
  outputTree->Branch("_jetpartonFlavour",&_jetpartonFlavour);

  
  outputTree->Branch("_lEta",&_lEta);
  outputTree->Branch("_lPhi",&_lPhi);
  outputTree->Branch("_lPt",&_lPt);
  outputTree->Branch("_lpdgId",&_lpdgId);

  outputTree->Branch("_lgenEta",&_lgenEta);
  outputTree->Branch("_lgenPhi",&_lgenPhi);
  outputTree->Branch("_lgenPt",&_lgenPt);
  outputTree->Branch("_lgenpdgId",&_lgenpdgId);


  outputTree->Branch("_phEta",&_phEta);
  outputTree->Branch("_phPhi",&_phPhi);
  outputTree->Branch("_phPt",&_phPt);
  
  outputTree->Branch("_PFcand_pt",&_PFcand_pt);
  outputTree->Branch("_PFcand_eta",&_PFcand_eta);
  outputTree->Branch("_PFcand_phi",&_PFcand_phi);
  outputTree->Branch("_PFcand_pdgId",&_PFcand_pdgId);
  outputTree->Branch("_PFcand_fromPV",&_PFcand_fromPV);

  //  if(IsMC_){
  outputTree->Branch("_genmet", &_genmet, "_genmet/D");
  outputTree->Branch("_genmet_phi", &_genmet_phi, "_genmet_phi/D");
  outputTree->Branch("trueNVtx", &trueNVtx,"trueNVtx/I");
  //}
  outputTree->Branch("_met", &_met, "_met/D");
  outputTree->Branch("_met_phi", &_met_phi, "_met_phi/D");
  outputTree->Branch("_puppimet", &_puppimet, "_puppimet/D");
  outputTree->Branch("_puppimet_phi", &_puppimet_phi, "_puppimet_phi/D");
  

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
  _jetPassID.clear();
  _jetPtGen.clear();
  _jetEtaGen.clear();
  _jetPhiGen.clear();
  _jetJECuncty.clear();
  _jetPUMVA.clear();
  _jetPtNoL2L3Res.clear();
  _jet_corrjecs.clear();
  _jetpartonFlavour.clear();
  _jethadronFlavour.clear();

  _lEta.clear();
  _lPhi.clear();
  _lPt.clear();
  _lpdgId.clear();

  _lgenEta.clear();
  _lgenPhi.clear();
  _lgenPt.clear();
  _lgenpdgId.clear();

  _phEta.clear();
  _phPhi.clear();
  _phPt.clear();

  _phgenEta.clear();
  _phgenPhi.clear();
  _phgenPt.clear();


  _PFcand_pt.clear();
  _PFcand_eta.clear();
  _PFcand_phi.clear();
  _PFcand_pdgId.clear();
  _PFcand_fromPV.clear();



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

//define this as a plug-in
DEFINE_FWK_MODULE(JMEAnalyzer);
