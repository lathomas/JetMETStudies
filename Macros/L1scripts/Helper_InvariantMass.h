#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
 
using namespace ROOT;
using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;



bool PassDiMumass4to56(ROOT::VecOps::RVec<float>pt, ROOT::VecOps::RVec<float>eta, ROOT::VecOps::RVec<float>phi, ROOT::VecOps::RVec<int>charge ){

  for(unsigned int i = 0; i<pt.size(); i++){
    for(unsigned int j = 0; j<pt.size(); j++){
      if(charge[i]==charge[j])continue;
      TLorentzVector mu1, mu2 ; 
      mu1.SetPtEtaPhiM(pt[i], eta[i], phi[i], 0.1);
      mu2.SetPtEtaPhiM(pt[j], eta[j], phi[j], 0.1);
      float mass = (mu1+mu2).Mag();
      if(mass>4||mass<56) return true;
    }
   
  }

 return false;
}


bool PassDiJet140_70_Mjj900(ROOT::VecOps::RVec<float>pt, ROOT::VecOps::RVec<float>eta, ROOT::VecOps::RVec<float>phi){

  for(unsigned int i = 0; i<pt.size(); i++){
    if(pt[i]<140)continue;
    for(unsigned int j = 0; j<pt.size(); j++){
      if(eta[i]*eta[j]>0) continue;
      if(pt[j]<70)continue;
      TLorentzVector jet1, jet2;
      jet1.SetPtEtaPhiM(pt[i], eta[i], phi[i], 0.);
      jet2.SetPtEtaPhiM(pt[j], eta[j], phi[j], 0.);
      float mass = (jet1+jet2).Mag();
      if(mass>900) return true;
    }
  }
  return false;
}



vector<float> HighestMjj(ROOT::VecOps::RVec<float>pt, ROOT::VecOps::RVec<float>eta, ROOT::VecOps::RVec<float>phi){

  vector<float> result;
  result.push_back(0);
  result.push_back(0);
  result.push_back(0);
  result.push_back(0);
  for(unsigned int i = 0; i<pt.size(); i++){
    if(pt[i]<20.)continue;
    for(unsigned int j = 0; j<pt.size(); j++){
      if(eta[i]*eta[j]>0) continue;
      if(pt[j]<20.)continue;
      TLorentzVector jet1, jet2;
      jet1.SetPtEtaPhiM(pt[i], eta[i], phi[i], 0.);
      jet2.SetPtEtaPhiM(pt[j], eta[j], phi[j], 0.);
      float mass = (jet1+jet2).Mag();
      if(mass>result[0]){
	result[0] = mass;
	result[3] = abs(eta[i]-eta[j]);
	if(pt[j]> pt[i]){ result[1] = pt[j] ; result[2] = pt[i] ;}
	else {  result[1] = pt[i] ; result[2] = pt[j] ; }
      }
    }
  }
  return result;
}
