#include "JetMETStudies/JMEAnalyzer/interface/Tools.h"
using namespace edm;

bool tools::PassJetID(const pat::Jet *iJ , string runera){
  
  bool vetol = false;
  
  if( (runera.find("2017") != string::npos || runera.find("2018") != string::npos )){

    if( TMath::Abs( iJ->eta() ) < 2.7 ){
      if( iJ->neutralHadronEnergyFraction() >= 0.9 ) return false;
      if( iJ->neutralEmEnergyFraction() >= 0.9 ) return false;
      if( (iJ->chargedMultiplicity() + iJ->neutralMultiplicity() ) < 2 ) return false;
      if( iJ->muonEnergyFraction()>=0.8 && vetol) return false;
      if( TMath::Abs( iJ->eta() ) < 2.4&&iJ->chargedMultiplicity() == 0 ) return false;
      if( TMath::Abs( iJ->eta() ) < 2.4&&iJ->chargedHadronEnergyFraction() <= 0. ) return false;
      if( TMath::Abs( iJ->eta() ) < 2.4&&iJ->chargedEmEnergyFraction() >= 0.8 &&vetol ) return false;
    }
    if( TMath::Abs( iJ->eta() ) > 2.7 &&TMath::Abs( iJ->eta() ) <  3.0 ){
      if( iJ->neutralEmEnergyFraction() >= 0.99 ) return false;
      if( iJ->neutralEmEnergyFraction() <= 0.02 ) return false;
      if( iJ->neutralMultiplicity()  <=2)  return false;
    }
    if( TMath::Abs( iJ->eta() ) >=3.0 ){
      if( iJ->neutralEmEnergyFraction() >= 0.9 ) return false;
      if( iJ->neutralHadronEnergyFraction() <= 0.02 ) return false;
      if(iJ->neutralMultiplicity()  <=10)  return false;
    }
  }

  if( (runera.find("2016") != string::npos  )){
    if( TMath::Abs( iJ->eta() ) < 2.7 ){
      if( iJ->neutralHadronEnergyFraction() >= 0.99 ) return false;
      if( iJ->neutralEmEnergyFraction() >= 0.99 ) return false;
      if( (iJ->chargedMultiplicity() + iJ->neutralMultiplicity() ) < 2 ) return false;
      if( iJ->muonEnergyFraction()>=0.8 && vetol) return false;
      if( TMath::Abs( iJ->eta() ) < 2.4&&iJ->chargedMultiplicity() == 0 ) return false;
      if( TMath::Abs( iJ->eta() ) < 2.4&&iJ->chargedHadronEnergyFraction() <= 0. ) return false;
      if( TMath::Abs( iJ->eta() ) < 2.4&&iJ->chargedEmEnergyFraction() >= 0.99 ) return false;
      if( TMath::Abs( iJ->eta() ) < 2.4&&iJ->chargedEmEnergyFraction() >= 0.8 &&vetol ) return false;
    }
    if( TMath::Abs( iJ->eta() ) > 2.7 &&TMath::Abs( iJ->eta() ) <  3.0 ){
      if(iJ->chargedHadronEnergyFraction()>=0.98 ) return false;
      if( iJ->neutralEmEnergyFraction() <= 0.01 ) return false;
      if(iJ->neutralMultiplicity()  <=2)  return false;
    }
    if( TMath::Abs( iJ->eta() ) >=3.0 ){
      if( iJ->neutralEmEnergyFraction() >= 0.9 ) return false;
      if( iJ->neutralHadronEnergyFraction() <= 0.02 ) return false;
      if(iJ->neutralMultiplicity()  <=10)  return false;
    }
  }
  return true;
}
