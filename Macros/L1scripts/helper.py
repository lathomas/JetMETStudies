import ROOT
from array import array
jetEtaBins = [0., 1.3, 2.5, 3., 3.5, 4., 5.]
egEtaBins = [0., 1.479, 2.5]
muEtaBins = [0., 0.83, 1.24, 2.4]


ht_bins = array('f', [ i*10 for i in range(50) ] + [ 500+ i*20 for i in range(25) ] + [1000 + i*50 for i in range(10)] +[1500,1600,1700,1800,2000,2500,3000])
leptonpt_bins = array('f',[ i for i in range(50) ] + [ 50+2*i for i in range(10) ] + [ 70+3*i for i in range(10) ] + [100+10*i for i in range(10) ] + [200, 250, 300, 400, 500])
jetmetpt_bins = array('f',[ i*5 for i in range(50) ] +  [250+10*i for i in range(25) ]  + [500+20*i for i in range(10) ] + [700, 800, 900, 1000, 1200, 1500, 2000 ])
#String printing stuff for a few events
stringToPrint = '''
if(EventsToPrint <100) {

cout << "*********New Event********"<<endl;
cout << "_runNb " << _runNb<<endl;

for(unsigned int i = 0;i< (_lPt).size();i++ ){
cout << "Lepton Pt, Eta, Phi: " << (_lPt)[i]<<", "<<(_lEta)[i]<<", "<<(_lPhi)[i]<<endl;
cout << "Lepton pdgId: " << (_lpdgId)[i]<<endl;
}
for(unsigned int i = 0;i< (_phPt).size();i++ ){
cout << "Photon Pt, Eta, Phi: " << (_phPt)[i]<<", "<<(_phEta)[i]<<", "<<(_phPhi)[i]<<endl;
}
for(unsigned int i = 0;i< (_jetPt).size();i++ ){
cout << "jet Pt, Eta, Phi: " << (_jetPt)[i]<<", "<<(_jetEta)[i]<<", "<<(_jetPhi)[i]<<endl;
cout << "jet PassID, MUEF, CEEF, CHEF: " << (_jetPassID)[i]<<", "<<(_jet_MUEF)[i]<<", "<<(_jet_CEEF)[i]<<", "<<  (_jet_CHEF)[i]<<endl;
}

for(unsigned int i = 0;i< (_L1eg_pt).size();i++ ){
cout << "L1 EG Pt, Eta, Phi: " << (_L1eg_pt)[i]<<", "<<(_L1eg_eta)[i]<<", "<<(_L1eg_phi)[i]<<endl;
}

cout << "probe_L1PtoverRecoPt probe_Pt sizes: "<<probe_L1PtoverRecoPt.size()<<", " <<probe_Pt.size()<<endl;
for(unsigned int i = 0;i< probe_Pt.size();i++ ){
cout << "probe_Pt, probe_L1PtoverRecoPt: "<<probe_L1PtoverRecoPt[i]<<", "<<probe_Pt[i]<<endl;
}

cout << "response, denominator_pt sizes: " << response.size()<<", " <<denominator_pt.size()<<endl;
for(unsigned int i = 0;i<response.size() ;i++ ){
cout <<"response, denominator_pt "<<response[i]<<", "<<denominator_pt[i]<<endl; 
} 
cout <<endl;

EventsToPrint++;
} 
return true;
'''










stringToPrintHF = '''
if(EventsToPrint <100) {

cout << "*********New Event********"<<endl;
cout << "_runNb " << _runNb<<endl;
cout << "met, met_phi " << _met <<", "<<_met_phi<<endl;
for(unsigned int i = 0;i< (_lPt).size();i++ ){
cout << "Lepton Pt, Eta, Phi: " << (_lPt)[i]<<", "<<(_lEta)[i]<<", "<<(_lPhi)[i]<<endl;
cout << "Lepton pdgId: " << (_lpdgId)[i]<<endl;
}
for(unsigned int i = 0;i< (_phPt).size();i++ ){
cout << "Photon Pt, Eta, Phi: " << (_phPt)[i]<<", "<<(_phEta)[i]<<", "<<(_phPhi)[i]<<endl;
}
for(unsigned int i = 0;i< (_jetPt).size();i++ ){
cout << "jet Pt, Eta, Phi: " << (_jetPt)[i]<<", "<<(_jetEta)[i]<<", "<<(_jetPhi)[i]<<endl;
cout << "jet PassID, MUEF, CEEF, CHEF, NHEF: " << (_jetPassID)[i]<<", "<<(_jet_MUEF)[i]<<", "<<(_jet_CEEF)[i]<<", "<<  (_jet_CHEF)[i]<<", "<<  (_jet_NHEF)[i]<<endl;
cout << "jethfsigmaEtaEta jethfsigmaPhiPhi jethfcentralEtaStripSize "<< (_jethfsigmaEtaEta)[i]<<", "<<  (_jethfsigmaPhiPhi)[i]<<", "<<(_jethfcentralEtaStripSize)[i]<<endl;
}


for(unsigned int i = 0;i< (_L1jet_pt).size();i++ ){
cout << "L1 JET Pt, Eta, Phi, Bx: " << (_L1jet_pt)[i]<<", "<<(_L1jet_eta)[i]<<", "<<(_L1jet_phi)[i]<<", " << (_L1jet_bx)[i]<<endl;
}



cout <<endl;

EventsToPrint++;
} 
return true;
'''






def SinglePhotonSelection(df):
    '''
    Select events with exactly one photon with pT>20 GeV.
    The event must pass a photon trigger. 
    '''
    df = df.Filter('_met<50')
    df = df.Filter('HLT_Photon110EB_TightID_TightIso')
    df = df.Define('photonsptgt20','_phPt>20')
    df = df.Filter('Sum(photonsptgt20)==1','=1 photon with p_{T}>20 GeV')
    
    df = df.Define('isRefPhoton','_phPassTightID&&_phPassIso&&_phPt>115&&abs(_phEta)<1.479')
    df = df.Filter('Sum(isRefPhoton)==1','Photon has p_{T}>115 GeV, passes tight ID and is in EB')
    
    df = df.Define('cleanPh_Pt','_phPt[isRefPhoton]')
    df = df.Define('cleanPh_Eta','_phEta[isRefPhoton]')
    df = df.Define('cleanPh_Phi','_phPhi[isRefPhoton]')
    
    df = df.Define('ref_Pt','cleanPh_Pt[0]')
    df = df.Define('ref_Phi','cleanPh_Phi[0]')
    
    return df
    
    
    
    
def MuonJet_MuonSelection(df):
    '''
    Select events with >= 1 muon with pT>25 GeV.
    The event must pass a single muon trigger. 
    '''
    df = df.Filter('HLT_IsoMu24')
    df = df.Define('goodmuonPt25','_lPt>25&&abs(_lpdgId)==13&&_lPassTightID')
    df = df.Filter('Sum(goodmuonPt25)>=1','>=1 muon with p_{T}>25 GeV')
    df = df.Define('badmuonPt10','_lPt>10&&abs(_lpdgId)==13&&_lPassTightID==0')
    df = df.Filter('Sum(badmuonPt10)==0','No bad quality muon')
    return df
    
    

def ZEE_EleSelection(df):
    df = df.Filter('HLT_Ele32_WPTight_Gsf')

    df = df.Define('isTag','_lPt>35&&abs(_lpdgId)==11&&_lPassTightID&&_lpassHLT_Ele32_WPTight_Gsf')
    df = df.Filter('Sum(isTag)>0')
    df = df.Define('isProbe','_lPt>5&&abs(_lpdgId)==11&&_lPassTightID&& (Sum(isTag)>=2|| isTag==0)')
    df = df.Filter('_mll>80&&_mll<100')

    df = df.Define('probe_Pt','_lPt[isProbe]')
    df = df.Define('probe_Eta','_lEta[isProbe]')
    df = df.Define('probe_Phi','_lPhi[isProbe]')
    
    return df


def ZMuMu_MuSelection(df):
    df = df.Filter('HLT_IsoMu24')

    df = df.Define('isTag','_lPt>25&&abs(_lpdgId)==13&&_lPassTightID&&_lpassHLT_IsoMu24')
    df = df.Filter('Sum(isTag)>0')
    df = df.Define('isProbe','_lPt>5&&abs(_lpdgId)==13&&_lPassTightID&& (Sum(isTag)>=2|| isTag==0)')
    df = df.Filter('_mll>80&&_mll<100')

    df = df.Define('probe_Pt','_lPt[isProbe]')
    df = df.Define('probe_Eta','_lEta[isProbe]')
    df = df.Define('probe_Phi','_lPhi[isProbe]')
    
    return df



    

def makehistosforturnons_inprobeetaranges(df, histos, etavarname, phivarname, ptvarname, responsevarname, etabins, l1varname, l1thresholds, prefix, binning, l1thresholdforeffvsrunnb, offlinethresholdforeffvsrunnb):
    '''Make histos for turnons vs pt (1D histos for numerator and denominator) in ranges of eta
    Also look at response vs run number (2D histo) '''
    for i in range(len(etabins)-1):
        str_bineta = "eta{}to{}".format(etabins[i],etabins[i+1]).replace(".","p")
        #Define columns corresponding to pt and response for the selected eta range 
        df_etarange = df.Define('inEtaRange','abs({})>={}'.format(etavarname, etabins[i])+'&&abs({})<{}'.format(etavarname, etabins[i+1]))
        df_etarange = df_etarange.Define('denominator_pt',ptvarname+'[inEtaRange]')
        df_etarange = df_etarange.Define('response',responsevarname+'[inEtaRange]')
        df_etarange = df_etarange.Define('runnb',"return ROOT::VecOps::RVec<int>(response.size(), _runNb);")

        #Response vs pt and vs runnb (2d)
        histos[prefix+str_bineta] = df_etarange.Histo1D(ROOT.RDF.TH1DModel('h_{}_{}'.format(prefix, str_bineta), '', len(binning)-1, binning), 'denominator_pt')
        histos[prefix+str_bineta+'_ResponseVsPt'] = df_etarange.Histo2D(ROOT.RDF.TH2DModel('h_ResponseVsPt_{}_{}'.format(prefix, str_bineta), '', 200, 0, 200, 100, 0, 2), 'denominator_pt', 'response')
        histos[prefix+str_bineta+'_ResponseVsRunNb'] = df_etarange.Histo2D(ROOT.RDF.TH2DModel('h_ResponseVsRunNb_{}_{}'.format(prefix, str_bineta), '', 5000, 355000, 360000, 100, 0 , 2), 'runnb', 'response')

        if i ==1 and prefix == 'EGNonIso_plots':
            df_etarange = df_etarange.Filter(stringToPrint)


        #Numerator/denominator for plateau eff vs runnb
        df_etarange = df_etarange.Define('inplateau','{}>={}&&inEtaRange'.format(ptvarname,offlinethresholdforeffvsrunnb))
        df_etarange = df_etarange.Define('N_inplateau','Sum(inplateau)')
        df_etarange = df_etarange.Define('runnb_inplateau',"return ROOT::VecOps::RVec<int>(N_inplateau, _runNb);")
        df_etarange = df_etarange.Define('inplateaupassL1','inplateau && {}>={}'.format(l1varname,l1thresholdforeffvsrunnb))
        df_etarange = df_etarange.Define('N_inplateaupassL1','Sum(inplateaupassL1)')
        df_etarange = df_etarange.Define('runnb_inplateaupassL1',"return ROOT::VecOps::RVec<int>(N_inplateaupassL1, _runNb);")
        histos[prefix+"_plateaueffvsrunnb_numerator_"+str_bineta] = df_etarange.Histo1D(ROOT.RDF.TH1DModel('h_PlateauEffVsRunNb_Numerator_{}_{}'.format(prefix, str_bineta), '', 5000, 355000, 360000),'runnb_inplateaupassL1')
        histos[prefix+"_plateaueffvsrunnb_denominator_"+str_bineta] = df_etarange.Histo1D(ROOT.RDF.TH1DModel('h_PlateauEffVsRunNb_Denominator_{}_{}'.format(prefix, str_bineta), '', 5000, 355000, 360000),'runnb_inplateau')

        for ipt in l1thresholds:
            str_binetapt = "eta{}to{}_l1thrgeq{}".format(etabins[i], etabins[i+1],ipt).replace(".","p")
            df_loc = df_etarange.Define('passL1Cond','{}>={}'.format(l1varname, ipt))
            df_loc = df_loc.Define('numerator_pt',ptvarname+'[inEtaRange&&passL1Cond]')
            histos[prefix+str_binetapt] = df_loc.Histo1D(ROOT.RDF.TH1DModel('h_{}_{}'.format(prefix, str_binetapt), '', len(binning)-1, binning), 'numerator_pt')

        

    return df #, histos

    

def ZEE_Plots(df):
    histos = {}
    
    label = ['EGNonIso','EGLooseIso', 'EGTightIso']
    df_eg = [None, None, None]
    for i in range(3): 
        
        if i ==0:
            df_eg[i] = df.Define('probe_idxL1jet','FindL1ObjIdx(_L1eg_eta, _L1eg_phi, probe_Eta, probe_Phi)')
        if i ==1:
            df_eg[i] = df.Define('probe_idxL1jet','FindL1ObjIdx(_L1eg_eta, _L1eg_phi, probe_Eta, probe_Phi, _L1eg_iso, 2)')
        if i ==2:
            df_eg[i] = df.Define('probe_idxL1jet','FindL1ObjIdx(_L1eg_eta, _L1eg_phi, probe_Eta, probe_Phi, _L1eg_iso, 3)')

        df_eg[i] = df_eg[i].Define('probe_L1Pt','GetVal(probe_idxL1jet, _L1eg_pt)')
        df_eg[i] = df_eg[i].Define('probe_L1Bx','GetVal(probe_idxL1jet, _L1eg_bx)')
        df_eg[i] = df_eg[i].Define('probe_L1PtoverRecoPt','probe_L1Pt/probe_Pt')
        
        df_eg[i] = makehistosforturnons_inprobeetaranges(df_eg[i], histos, etavarname='probe_Eta', phivarname='probe_Phi', ptvarname='probe_Pt', responsevarname='probe_L1PtoverRecoPt', etabins=egEtaBins, l1varname='probe_L1Pt', l1thresholds=[5,10,15,20,25,30,35], prefix=label[i]+"_plots", binning=leptonpt_bins, l1thresholdforeffvsrunnb = 30, offlinethresholdforeffvsrunnb = 35 )
        
        
        df_eg[i] = df_eg[i].Define('probePt30_Eta','probe_Eta[probe_Pt>30]')
        df_eg[i] = df_eg[i].Define('probePt30_Phi','probe_Phi[probe_Pt>30]')
        df_eg[i] = df_eg[i].Define('probePt30PassL1EG25_Eta','probe_Eta[probe_Pt>30&&probe_L1Pt>25]')
        df_eg[i] = df_eg[i].Define('probePt30PassL1EG25_Phi','probe_Phi[probe_Pt>30&&probe_L1Pt>25]')
        histos['h_EG25_EtaPhi_Numerator'+label[i]] = df_eg[i].Histo2D(ROOT.RDF.TH2DModel('h_EG25_EtaPhi_Numerator', '', 100, -5,5, 100, -3.1416, 3.1416), 'probePt30PassL1EG25_Eta', 'probePt30PassL1EG25_Phi')
        histos['h_EG25_EtaPhi_Denominator'+label[i]] = df_eg[i].Histo2D(ROOT.RDF.TH2DModel('h_EG25_EtaPhi_Denominator', '', 100, -5,5, 100, -3.1416, 3.1416), 'probePt30_Eta', 'probePt30_Phi')

        
        if i ==0:
            df_eg[i] = df_eg[i].Define('probeL1EG15to26Bxmin1_Eta','probe_Eta[probe_Pt>12&&probe_Pt<23&&probe_L1Pt>15&&probe_L1Pt<=26&&probe_L1Bx==-1]')
            df_eg[i] = df_eg[i].Define('probeL1EG15to26Bxmin1_Phi','probe_Phi[probe_Pt>12&&probe_Pt<23&&probe_L1Pt>15&&probe_L1Pt<=26&&probe_L1Bx==-1]')
            df_eg[i] = df_eg[i].Define('probeL1EG15to26Bx0_Eta','probe_Eta[probe_Pt>12&&probe_Pt<23&&probe_L1Pt>15&&probe_L1Pt<=26&&probe_L1Bx==0]')
            df_eg[i] = df_eg[i].Define('probeL1EG15to26Bx0_Phi','probe_Phi[probe_Pt>12&&probe_Pt<23&&probe_L1Pt>15&&probe_L1Pt<=26&&probe_L1Bx==0]')
            df_eg[i] = df_eg[i].Define('probeL1EG15to26Bxplus1_Eta','probe_Eta[probe_Pt>12&&probe_Pt<23&&probe_L1Pt>15&&probe_L1Pt<=26&&probe_L1Bx==1]')
            df_eg[i] = df_eg[i].Define('probeL1EG15to26Bxplus1_Phi','probe_Phi[probe_Pt>12&&probe_Pt<23&&probe_L1Pt>15&&probe_L1Pt<=26&&probe_L1Bx==1]')
            
            
            histos['L1EG15to26_bxmin1_etaphi'] = df_eg[i].Histo2D(ROOT.RDF.TH2DModel('L1EG15to26_bxmin1_etaphi', '', 100, -5,5, 100, -3.1416, 3.1416), 'probeL1EG15to26Bxmin1_Eta', 'probeL1EG15to26Bxmin1_Phi')
            histos['L1EG15to26_bx0_etaphi'] = df_eg[i].Histo2D(ROOT.RDF.TH2DModel('L1EG15to26_bx0_etaphi', '', 100, -5,5, 100, -3.1416, 3.1416), 'probeL1EG15to26Bx0_Eta', 'probeL1EG15to26Bx0_Phi')
            histos['L1EG15to26_bxplus1_etaphi'] = df_eg[i].Histo2D(ROOT.RDF.TH2DModel('L1EG15to26_bxplus1_etaphi', '', 100, -5,5, 100, -3.1416, 3.1416), 'probeL1EG15to26Bxplus1_Eta', 'probeL1EG15to26Bxplus1_Phi')

    return df, histos
    


def ZMuMu_Plots(df):
                                                                                          
    histos = {}
    label = ['AllQual', 'Qual8', 'Qual12']
    df_mu = [None, None, None]
    
    for i in range(3):
        if i ==0:
            df_mu[i] = df.Define('probe_idxL1jet','FindL1ObjIdx(_L1mu_eta, _L1mu_phi, probe_Eta, probe_Phi)')
        if i ==1:
            df_mu[i] = df.Define('probe_idxL1jet','FindL1ObjIdx(_L1mu_eta, _L1mu_phi, probe_Eta, probe_Phi, _L1mu_Qual, 8)')
        if i ==2:
            df_mu[i] = df.Define('probe_idxL1jet','FindL1ObjIdx(_L1mu_eta, _L1mu_phi, probe_Eta, probe_Phi, _L1mu_Qual, 12)')

            

        df_mu[i] = df_mu[i].Define('probe_L1Pt','GetVal(probe_idxL1jet, _L1mu_pt)')
        df_mu[i] = df_mu[i].Define('probe_L1Bx','GetVal(probe_idxL1jet, _L1mu_bx)')
        df_mu[i] = df_mu[i].Define('probe_L1Qual','GetVal(probe_idxL1jet, _L1mu_Qual)')
        df_mu[i] = df_mu[i].Define('probe_L1PtoverRecoPt','probe_L1Pt/probe_Pt')
        
        
        df_mu[i] = makehistosforturnons_inprobeetaranges(df_mu[i], histos, etavarname='probe_Eta', phivarname='probe_Phi', ptvarname='probe_Pt', responsevarname='probe_L1PtoverRecoPt', etabins=muEtaBins, l1varname='probe_L1Pt', l1thresholds=[5,10,15,20,22,25],  prefix=label[i]+"_plots" , binning = leptonpt_bins, l1thresholdforeffvsrunnb = 22, offlinethresholdforeffvsrunnb = 27)
        
    
        df_mu[i] = df_mu[i].Define('probePt30_Eta','probe_Eta[probe_Pt>30]')
        df_mu[i] = df_mu[i].Define('probePt30_Phi','probe_Phi[probe_Pt>30]')
        df_mu[i] = df_mu[i].Define('probePt30PassL1Mu22_Eta','probe_Eta[probe_Pt>30&&probe_L1Pt>22]')
        df_mu[i] = df_mu[i].Define('probePt30PassL1Mu22_Phi','probe_Phi[probe_Pt>30&&probe_L1Pt>22]')
        histos['h_Mu22_EtaPhi_Denominator'+label[i]] = df_mu[i].Histo2D(ROOT.RDF.TH2DModel('h_Mu22_EtaPhi_Denominator', '', 100, -5,5, 100, -3.1416, 3.1416), 'probePt30_Eta', 'probePt30_Phi')
        histos['h_Mu22_EtaPhi_Numerator'+label[i]] = df_mu[i].Histo2D(ROOT.RDF.TH2DModel('h_Mu22_EtaPhi_Numerator', '', 100, -5,5, 100, -3.1416, 3.1416), 'probePt30PassL1Mu22_Eta', 'probePt30PassL1Mu22_Phi')
        
        
        if i== 2:
            df_mu[i] = df_mu[i].Define('probeL1Mu10to21Bxmin1_Eta','probe_Eta[probe_Pt>8&&probe_Pt<25&&probe_L1Pt>10&&probe_L1Pt<=21&&probe_L1Bx==-1&&probe_L1Qual>=12]')
            df_mu[i] = df_mu[i].Define('probeL1Mu10to21Bxmin1_Phi','probe_Phi[probe_Pt>8&&probe_Pt<25&&probe_L1Pt>10&&probe_L1Pt<=21&&probe_L1Bx==-1&&probe_L1Qual>=12]')
            df_mu[i] = df_mu[i].Define('probeL1Mu10to21Bx0_Eta','probe_Eta[probe_Pt>8&&probe_Pt<25&&probe_L1Pt>10&&probe_L1Pt<=21&&probe_L1Bx==0&&probe_L1Qual>=12]')
            df_mu[i] = df_mu[i].Define('probeL1Mu10to21Bx0_Phi','probe_Phi[probe_Pt>8&&probe_Pt<25&&probe_L1Pt>10&&probe_L1Pt<=21&&probe_L1Bx==0&&probe_L1Qual>=12]')
            df_mu[i] = df_mu[i].Define('probeL1Mu10to21Bxplus1_Eta','probe_Eta[probe_Pt>8&&probe_Pt<25&&probe_L1Pt>10&&probe_L1Pt<=21&&probe_L1Bx==1&&probe_L1Qual>=12]')
            df_mu[i] = df_mu[i].Define('probeL1Mu10to21Bxplus1_Phi','probe_Phi[probe_Pt>8&&probe_Pt<25&&probe_L1Pt>10&&probe_L1Pt<=21&&probe_L1Bx==1&&probe_L1Qual>=12]')
            
            
            histos['L1Mu10to21_bxmin1_etaphi'] = df_mu[i].Histo2D(ROOT.RDF.TH2DModel('L1Mu10to21_bxmin1_etaphi', '', 100, -5,5, 100, -3.1416, 3.1416), 'probeL1Mu10to21Bxmin1_Eta', 'probeL1Mu10to21Bxmin1_Phi')
            histos['L1Mu10to21_bx0_etaphi'] = df_mu[i].Histo2D(ROOT.RDF.TH2DModel('L1Mu10to21_bx0_etaphi', '', 100, -5,5, 100, -3.1416, 3.1416), 'probeL1Mu10to21Bx0_Eta', 'probeL1Mu10to21Bx0_Phi')
            histos['L1Mu10to21_bxplus1_etaphi'] = df_mu[i].Histo2D(ROOT.RDF.TH2DModel('L1Mu10to21_bxplus1_etaphi', '', 100, -5,5, 100, -3.1416, 3.1416), 'probeL1Mu10to21Bxplus1_Eta', 'probeL1Mu10to21Bxplus1_Phi')

    return df, histos
    




def CleanJets(df):
    #List of cleaned jets (noise cleaning + lepton/photon overlap removal)
    df = df.Define('isCleanJet','_jetPassID&&_jetLeptonPhotonCleaned&&_jetPt>30')
    df = df.Define('cleanJet_Pt','_jetPt[isCleanJet]')
    df = df.Define('cleanJet_Eta','_jetEta[isCleanJet]')
    df = df.Define('cleanJet_Phi','_jetPhi[isCleanJet]')
    
    df = df.Filter('Sum(isCleanJet)>=1','>=1 clean jet with p_{T}>30 GeV')

    return df



def EtSum(df):
    histos = {}
    #HT=scalar pt sum of all jets with pt>30 and |eta|<2.5 
    df = df.Define('iscentraljet','cleanJet_Pt>30&&abs(cleanJet_Eta)<2.5')
    #df = df.Filter('Sum(iscentraljet)>0')
    df = df.Define('HT','Sum(cleanJet_Pt[cleanJet_Pt>30&&abs(cleanJet_Eta)<2.5])')

    df = df.Define('muons_px','Sum(_lPt[abs(_lpdgId)==13]*cos(_lPhi[abs(_lpdgId)==13]))')
    df = df.Define('muons_py','Sum(_lPt[abs(_lpdgId)==13]*sin(_lPhi[abs(_lpdgId)==13]))')
    df = df.Define('metnomu_x','_met*cos(_met_phi)+muons_px')
    df = df.Define('metnomu_y','_met*sin(_met_phi)+muons_py')
    df = df.Define('MetNoMu','sqrt(metnomu_x*metnomu_x+metnomu_y*metnomu_y)')
    df = df.Define('L1_ETMHF80','passL1_Initial_bx0[419]')
    df = df.Define('L1_ETMHF100','passL1_Initial_bx0[421]')

    histos['h_MetNoMu_Denominator'] = df.Histo1D(ROOT.RDF.TH1DModel('h_MetNoMu_Denominator', '', len(jetmetpt_bins)-1, array('d',jetmetpt_bins)), 'MetNoMu') 
    
    dfmetl1 = df.Filter('passL1_Initial_bx0[419]')
    histos['L1_ETMHF80'] = dfmetl1.Histo1D(ROOT.RDF.TH1DModel('h_MetNoMu_ETMHF80', '', len(jetmetpt_bins)-1, array('d',jetmetpt_bins)), 'MetNoMu')
    dfmetl1 = df.Filter('passL1_Initial_bx0[421]')
    histos['L1_ETMHF100'] = dfmetl1.Histo1D(ROOT.RDF.TH1DModel('h_MetNoMu_ETMHF100', '', len(jetmetpt_bins)-1, array('d',jetmetpt_bins)), 'MetNoMu')

    histos['HLT_PFMETNoMu120_PFMHTNoMu120_IDTight'] =  df.Filter('HLT_PFMETNoMu120_PFMHTNoMu120_IDTight').Histo1D(ROOT.RDF.TH1DModel('h_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight', '', len(jetmetpt_bins)-1, array('d',jetmetpt_bins)), 'MetNoMu')
    histos['HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60'] =  df.Filter('HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60').Histo1D(ROOT.RDF.TH1DModel('h_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60', '', len(jetmetpt_bins)-1, array('d',jetmetpt_bins)), 'MetNoMu')
    

    histos['h_HT_Denominator'] = df.Filter('_met<50').Histo1D(ROOT.RDF.TH1DModel('h_HT_Denominator', '', len(ht_bins)-1, array('d',ht_bins)), 'HT') 
    histos['L1_HTT200er'] = df.Filter('passL1_Initial_bx0[400]').Filter('_met<50').Histo1D(ROOT.RDF.TH1DModel('h_HT_L1_HTT200er', '', len(ht_bins)-1, array('d',ht_bins)), 'HT')  
    histos['L1_HTT280er'] = df.Filter('passL1_Initial_bx0[402]').Filter('_met<50').Histo1D(ROOT.RDF.TH1DModel('h_HT_L1_HTT280er', '', len(ht_bins)-1, array('d',ht_bins)), 'HT')
    histos['L1_HTT360er'] = df.Filter('passL1_Initial_bx0[404]').Filter('_met<50').Histo1D(ROOT.RDF.TH1DModel('h_HT_L1_HTT360er', '', len(ht_bins)-1, array('d',ht_bins)), 'HT')
    histos['HLT_PFHT1050'] =  df.Filter('HLT_PFHT1050').Filter('_met<50').Histo1D(ROOT.RDF.TH1DModel('h_HLT_PFHT1050', '', len(ht_bins)-1, array('d',ht_bins)), 'HT')

    return df, histos

def AnalyzeCleanJets(df, JetRecoPtCut, L1JetPtCut):    
    histos = {}
    #Find L1 jets matched to the offline jet
    df = df.Define('cleanJet_idxL1jet','FindL1ObjIdx(_L1jet_eta, _L1jet_phi, cleanJet_Eta, cleanJet_Phi)')
    df = df.Define('cleanJet_L1Pt','GetVal(cleanJet_idxL1jet,_L1jet_pt)')
    df = df.Define('cleanJet_L1Bx','GetVal(cleanJet_idxL1jet,_L1jet_bx)')
    df = df.Define('cleanJet_L1PtoverRecoPt','cleanJet_L1Pt/cleanJet_Pt')
    #Now some plotting (turn ons for now)
    L1PtCuts = [30., 40., 60., 80., 100., 120., 140., 160., 180., 200.]


    df = makehistosforturnons_inprobeetaranges(df, histos, etavarname='cleanJet_Eta', phivarname='cleanJet_Phi', ptvarname='cleanJet_Pt', responsevarname='cleanJet_L1PtoverRecoPt', etabins=jetEtaBins, l1varname='cleanJet_L1Pt', l1thresholds=L1PtCuts, prefix="Jet_plots", binning=jetmetpt_bins, l1thresholdforeffvsrunnb = L1JetPtCut, offlinethresholdforeffvsrunnb = JetRecoPtCut )

    '''
    for i in range(len(jetEtaBins)-1):
        str_bineta = "eta{}to{}".format(jetEtaBins[i],jetEtaBins[i+1]).replace(".","p")
        df_JetsBinnedInEta = df.Define('cleanJet_EtaRestricted','abs(cleanJet_Eta)>={}&&abs(cleanJet_Eta)<{}'.format(jetEtaBins[i],jetEtaBins[i+1]))
        df_JetsBinnedInEta = df_JetsBinnedInEta.Define('cleanJetPt_EtaRestricted','cleanJet_Pt[cleanJet_EtaRestricted]')
        histos[str_bineta] = df_JetsBinnedInEta.Histo1D(ROOT.RDF.TH1DModel('h_jetpt_{}'.format(str_bineta), '', len(jetmetpt_bins)-1, array('d',jetmetpt_bins)), 'cleanJetPt_EtaRestricted')

        df_JetsBinnedInEta = df_JetsBinnedInEta.Define('cleanJetL1Pt_EtaRestricted','cleanJet_L1Pt[cleanJet_EtaRestricted]')
        df_JetsBinnedInEta = df_JetsBinnedInEta.Define('cleanJetL1OverRecoPt_EtaRestricted','cleanJetL1Pt_EtaRestricted/cleanJetPt_EtaRestricted')
        df_JetsBinnedInEta = df_JetsBinnedInEta.Define('runnb',"return ROOT::VecOps::RVec<int>(cleanJetL1OverRecoPt_EtaRestricted.size(), _runNb);")
        histos[str_bineta+'_ResponseVsRunNb'] = df_JetsBinnedInEta.Histo2D(ROOT.RDF.TH2DModel('h_JetResponseVsRunNb_{}'.format(str_bineta), '', 5000, 355000, 360000, 100, 0 , 2), 'runnb', 'cleanJetL1OverRecoPt_EtaRestricted')
    
        for ipt in L1PtCuts:
            str_binetapt = "eta{}to{}_L1ptgt{}".format(jetEtaBins[i],jetEtaBins[i+1],ipt).replace(".","p")
            df_JetsBinnedInEtaPt = df_JetsBinnedInEta.Define('cleanJet_L1Cond','cleanJet_L1Pt>={}'.format(ipt))
            df_JetsBinnedInEtaPt = df_JetsBinnedInEtaPt.Define('cleanJetPt_Numerator','cleanJet_Pt[cleanJet_EtaRestricted&&cleanJet_L1Cond]')
            histos[str_binetapt] = df_JetsBinnedInEtaPt.Histo1D(ROOT.RDF.TH1DModel('h_jetpt_{}'.format(str_binetapt), '', len(jetmetpt_bins)-1, array('d',jetmetpt_bins)), 'cleanJetPt_Numerator')
    '''        

    
    df = df.Define('cleanHighPtJet_Eta','cleanJet_Eta[cleanJet_Pt>{}]'.format(JetRecoPtCut))
    df = df.Define('cleanHighPtJet_Phi','cleanJet_Phi[cleanJet_Pt>{}]'.format(JetRecoPtCut))
    df = df.Define('cleanHighPtJet_Eta_PassL1Jet','cleanJet_Eta[cleanJet_L1Pt>={}&&cleanJet_Pt>{}]'.format(L1JetPtCut, JetRecoPtCut))
    df = df.Define('cleanHighPtJet_Phi_PassL1Jet','cleanJet_Phi[cleanJet_L1Pt>={}&&cleanJet_Pt>{}]'.format(L1JetPtCut, JetRecoPtCut))
    
    histos["L1JetvsEtaPhi_Numerator"] = df.Histo2D(ROOT.RDF.TH2DModel('h_L1Jet{}vsEtaPhi_Numerator'.format(int(L1JetPtCut)), '', 100,-5,5,100,-3.1416,3.1416), 'cleanHighPtJet_Eta_PassL1Jet','cleanHighPtJet_Phi_PassL1Jet')
    histos["L1JetvsEtaPhi_EtaRestricted"] = df.Histo2D(ROOT.RDF.TH2DModel('h_L1Jet{}vsEtaPhi_EtaRestricted'.format(int(L1JetPtCut)), '', 100,-5,5,100,-3.1416,3.1416), 'cleanHighPtJet_Eta','cleanHighPtJet_Phi')
    


    df = df.Define('probeL1Jet100to150Bxmin1_Eta','cleanJet_Eta[cleanJet_Pt>90&&cleanJet_Pt<160&&cleanJet_L1Pt>100&&cleanJet_L1Pt<=150&&cleanJet_L1Bx==-1]')
    df = df.Define('probeL1Jet100to150Bxmin1_Phi','cleanJet_Phi[cleanJet_Pt>90&&cleanJet_Pt<160&&cleanJet_L1Pt>100&&cleanJet_L1Pt<=150&&cleanJet_L1Bx==-1]')
    df = df.Define('probeL1Jet100to150Bx0_Eta','cleanJet_Eta[cleanJet_Pt>90&&cleanJet_Pt<160&&cleanJet_L1Pt>100&&cleanJet_L1Pt<=150&&cleanJet_L1Bx==0]')
    df = df.Define('probeL1Jet100to150Bx0_Phi','cleanJet_Phi[cleanJet_Pt>90&&cleanJet_Pt<160&&cleanJet_L1Pt>100&&cleanJet_L1Pt<=150&&cleanJet_L1Bx==0]')
    df = df.Define('probeL1Jet100to150Bxplus1_Eta','cleanJet_Eta[cleanJet_Pt>90&&cleanJet_Pt<160&&cleanJet_L1Pt>100&&cleanJet_L1Pt<=150&&cleanJet_L1Bx==1]')
    df = df.Define('probeL1Jet100to150Bxplus1_Phi','cleanJet_Phi[cleanJet_Pt>90&&cleanJet_Pt<160&&cleanJet_L1Pt>100&&cleanJet_L1Pt<=150&&cleanJet_L1Bx==1]')




    histos['L1Jet100to150_bxmin1_etaphi'] = df.Histo2D(ROOT.RDF.TH2DModel('L1Jet100to150_bxmin1_etaphi', '', 100, -5,5, 100, -3.1416, 3.1416), 'probeL1Jet100to150Bxmin1_Eta', 'probeL1Jet100to150Bxmin1_Phi')
    histos['L1Jet100to150_bx0_etaphi'] = df.Histo2D(ROOT.RDF.TH2DModel('L1Jet100to150_bx0_etaphi', '', 100, -5,5, 100, -3.1416, 3.1416), 'probeL1Jet100to150Bx0_Eta', 'probeL1Jet100to150Bx0_Phi')
    histos['L1Jet100to150_bxplus1_etaphi'] = df.Histo2D(ROOT.RDF.TH2DModel('L1Jet100to150_bxplus1_etaphi', '', 100, -5,5, 100, -3.1416, 3.1416), 'probeL1Jet100to150Bxplus1_Eta', 'probeL1Jet100to150Bxplus1_Phi')



    return df, histos



def HFNoiseStudy(df):
    histos={}
    #debugging what happens in the HF 
    dfhfnoise = df.Define('isHFJet','cleanJet_Pt>30&&((cleanJet_Eta>3.0&&cleanJet_Eta<5)||(cleanJet_Eta>-5&&cleanJet_Eta<-3.))' )
    dfhfnoise = dfhfnoise.Define('nHFJets','Sum(isHFJet)')
    dfhfnoise = dfhfnoise.Define('isHFJetPt150','cleanJet_Pt>150&&((cleanJet_Eta>3.0&&cleanJet_Eta<5)||(cleanJet_Eta>-5&&cleanJet_Eta<-3.))' )
    dfhfnoise = dfhfnoise.Filter('Sum(isHFJetPt150)>0')
    dfhfnoise = dfhfnoise.Define('passL1Pt80','cleanJet_L1Pt>80')
    
    dfhfnoise = dfhfnoise.Define('HighPtHFJet_Pt','cleanJet_Pt[isHFJetPt150]')
    dfhfnoise = dfhfnoise.Define('HighPtHFJet_Eta','cleanJet_Eta[isHFJetPt150]')
    
    

    suffix = ['failL1Jet80', 'passL1Jet80']
    df_passvsfailL1 = [dfhfnoise.Filter('Sum(passL1Pt80&&isHFJetPt150)==0'), dfhfnoise.Filter('Sum(passL1Pt80&&isHFJetPt150)>=1')]
    
    for i, s in enumerate(suffix):
        if i == 0:
            df_passvsfailL1[i] = df_passvsfailL1[i].Filter(stringToPrintHF)
        
        histos[s+'_photonpt'] = df_passvsfailL1[i].Histo1D(ROOT.RDF.TH1DModel(s+'_photonpt', '', 100, 0, 1000), 'ref_Pt')
        histos[s+'_nhfjets'] = df_passvsfailL1[i].Histo1D(ROOT.RDF.TH1DModel(s+'_nhfjets', '', 100, 0, 100), 'nHFJets')
        histos[s+'_npv'] = df_passvsfailL1[i].Histo1D(ROOT.RDF.TH1DModel(s+'_npv', '', 100, 0, 100), '_n_PV')
        histos[s+'_runnb'] = df_passvsfailL1[i].Histo1D(ROOT.RDF.TH1DModel(s+'_runnb', '', 5000, 355000, 360000), '_runNb')
        histos[s+'_ptbalance'] = df_passvsfailL1[i].Histo1D(ROOT.RDF.TH1DModel(s+'_ptbalance', '', 100,0,2), 'ptbalance')
        histos[s+'_ptbalanceL1'] = df_passvsfailL1[i].Histo1D(ROOT.RDF.TH1DModel(s+'_ptbalanceL1', '', 100,0,2), 'ptbalanceL1')

    return df, histos
    
def PtBalanceSelection(df):
    '''
    Compute pt balance = pt(jet)/pt(ref)
    ref can be a photon or a Z.
    '''
    
    #Back to back condition
    df = df.Filter('abs(acos(cos(ref_Phi-cleanJet_Phi[0])))>2.9','DeltaPhi(ph,jet)>2.9')

    #Compute Pt balance = pt(jet)/pt(ref) => here ref is a photon
    #Reco first
    df = df.Define('ptbalance','cleanJet_Pt[0]/ref_Pt')
    df = df.Define('ptbalanceL1','_L1jet_pt[cleanJet_idxL1jet[0]]/ref_Pt')
    df = df.Define('probe_Eta','cleanJet_Eta[0]') 
    df = df.Define('probe_Phi','cleanJet_Phi[0]')
    return df

def AnalyzePtBalance(df):
    histos = {}
    df_JetsBinnedInEta ={}
    histos['L1JetvsEtaPhi'] = df.Histo3D(ROOT.RDF.TH3DModel('h_L1PtBalanceVsEtaPhi', 'ptbalanceL1', 100, -5, 5, 100, -3.1416, 3.1416, 100, 0, 2), 'probe_Eta','probe_Phi','ptbalanceL1')
    for i in range(len(jetEtaBins)-1):
        str_bineta = "eta{}to{}".format(jetEtaBins[i],jetEtaBins[i+1]).replace(".","p")
        df_JetsBinnedInEta[str_bineta] = df.Filter('abs(cleanJet_Eta[0])>={}&&abs(cleanJet_Eta[0])<{}'.format(jetEtaBins[i],jetEtaBins[i+1]))
        histos['RecoJetvsRunNb'+str_bineta] = df_JetsBinnedInEta[str_bineta].Histo2D(ROOT.RDF.TH2DModel('h_PtBalanceVsRunNb_{}'.format(str_bineta), 'ptbalance', 5000, 355000,360000, 100,0,2), '_runNb','ptbalance')
        histos['L1JetvsRunNb'+str_bineta] = df_JetsBinnedInEta[str_bineta].Histo2D(ROOT.RDF.TH2DModel('h_L1PtBalanceVsRunNb_{}'.format(str_bineta), 'ptbalanceL1', 5000, 355000,360000, 100,0,2), '_runNb','ptbalanceL1')
        histos['L1JetvsPU'+str_bineta] = df_JetsBinnedInEta[str_bineta].Histo2D(ROOT.RDF.TH2DModel('h_L1PtBalanceVsPU_{}'.format(str_bineta), 'ptbalanceL1', 100, 0, 100, 100, 0, 2), '_n_PV','ptbalanceL1')

        
    return df, histos



    
