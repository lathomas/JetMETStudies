# Standard importts
import os,sys,socket,argparse
import os
import ROOT
import math
from array import array
import numpy as np

ROOT.gROOT.SetBatch(True)
ReducedBinning = False
correctforMCPU = False

'''
Still to do 
merge bins at high pt high eta
add - eta side

check impact with tight id
'''
# RooFit
ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libRooFitCore.so")
ROOT.gROOT.SetStyle("Plain") # Not sure this is needed
ROOT.gSystem.SetIncludePath( "-I$ROOFITSYS/include/" )

def ConvFit( ptbalance, pu, pli, jer, output, pt, eta, binning, isData=False):

    print "Performing a fit using convoluted templates and a pdf to get the mean and the width for the JER" 
    # Declare the observable mean, and import the histogram to a RooDataHist

    tmp_sigma    = ptbalance.GetRMS()
    ptfloat = float(pt.split("to")[0])
    if binning == "zpt":
        variable_name = "p_{T}^{Z}"
    elif binning == "jetandzpt":
        variable_name = "p_{T}^{Z}"
    else:
        variable_name = "p_{T}^{jet}"
        
    balance      = ROOT.RooRealVar("balance","p_{T}_{Jet}/p_{T}_{Z}",0., 3.) ;
    dh_ptbalance = ROOT.RooDataHist("dh_ptbalance"  ,"dh_ptbalance"  ,ROOT.RooArgList(balance),ROOT.RooFit.Import(ptbalance)) ;

    # plot the data hist with error from sum of weighted events
    frame        = balance.frame(ROOT.RooFit.Title("ptbalance"))
    if isData:
        dh_ptbalance.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    else:
        dh_ptbalance.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2)) ;
     
    # Now import the templates
    dh_pu  = ROOT.RooDataHist("dh_pu",  "dh_pu" , ROOT.RooArgList(balance),ROOT.RooFit.Import(pu)) ;
    dh_pli = ROOT.RooDataHist("dh_pli", "dh_pli", ROOT.RooArgList(balance),ROOT.RooFit.Import(pli)) ;
    dh_jer = ROOT.RooDataHist("dh_jer", "dh_jer", ROOT.RooArgList(balance),ROOT.RooFit.Import(jer)) ;
    
    #Convert them to pdf
    pdf_pu  = ROOT.RooHistPdf("pdf_pu" , "pdf_pu" , ROOT.RooArgSet(balance), dh_pu);
    pdf_pli = ROOT.RooHistPdf("pdf_pli", "pdf_pli", ROOT.RooArgSet(balance), dh_pli);
    pdf_jer = ROOT.RooHistPdf("pdf_jer", "pdf_jer", ROOT.RooArgSet(balance), dh_jer);
    
    # create a simple gaussian pdf
    #gauss_mean  = ROOT.RooRealVar("mean","mean",0) #use this line if you don't want the mean (=scale) to float
    gauss_mean  = ROOT.RooRealVar("mean","mean",0.001,-0.1,0.1)
    gauss_sigma = ROOT.RooRealVar("sigma jer","sigma gauss",tmp_sigma,0.01,0.5)
    gauss       = ROOT.RooGaussian("gauss","gauss",balance,gauss_mean,gauss_sigma) 
    
    #Finale template for signal : PLI (gen level) convoluted with gauss
    #right now don't know how to do a 3 component convolution so only taking pli and gauss
    tmpxg = ROOT.RooFFTConvPdf("tmpxg","templates x gauss" ,balance,pdf_pli, gauss)
    tmpxg.setBufferFraction(0.7)

    #Alternatively can also convolute with a logN
    logn_mean = ROOT.RooRealVar("lognmean","lognmean",0.,10.)
    logn_sigma = ROOT.RooRealVar("lognsigma","lognsigma",0.,10.)
    logn = ROOT.RooLognormal("logn","logn",balance,logn_mean,logn_sigma)
    tmpxlogn = ROOT.RooFFTConvPdf("tmpxlogn","templates x logn" ,balance,pdf_pli, logn)

    #To let the fractions be calculated on the fly, need to convert them to extended pdfs
    nentries  = ptbalance.Integral()
    npu       = ROOT.RooRealVar("npu","npu",0,nentries)
    n         = ROOT.RooRealVar("n","n",0,nentries)
    
    extpdf_pu    = ROOT.RooExtendPdf("extpdf_pu"   , "extpdf_pu"   , pdf_pu   , npu)
    extpdf_pli   = ROOT.RooExtendPdf("extpdf_pli"  , "extpdf_pli"  , pdf_pli  , n)
    extpdf_tmpxg = ROOT.RooExtendPdf("extpdf_tmpxg", "extpdf_tmpxg", tmpxg, n)
    extpdf_tmpxlogn = ROOT.RooExtendPdf("extpdf_tmpxlogn", "extpdf_tmpxlogn", tmpxlogn, n)


    #Try background template from analytical functions
    a = ROOT.RooRealVar("a","a",0,nentries)
    b = ROOT.RooRealVar("b","b",0,100)
    c = ROOT.RooRealVar("c","c",0)
    d = ROOT.RooRealVar("d","d",0,100)
    gp = ROOT.RooGenericPdf("gp","Generic PDF","a*exp(-b*(balance-c))*pow((balance),d)", ROOT.RooArgList(balance,a,b,c,d))
    extpdf_gp = ROOT.RooExtendPdf("extpdf_gp","extpdf_gp",gp,npu)
    

    # Second gaussian to smear PU template
    gausspu_mean  = ROOT.RooRealVar("meanpu","mean pu",0.0)
    gausspu_sigma = ROOT.RooRealVar("sigmapu jer","sigma gauss pu",0.05,0.0,0.15)
    gausspu       = ROOT.RooGaussian("gausspu","gauss",balance,gausspu_mean,gausspu_sigma)
    
    tmpxgpu = ROOT.RooFFTConvPdf("tmpxgpu","templates x gauss pu" ,balance,pdf_pu, gausspu)
    tmpxgpu.setBufferFraction(0.7)
    extpdf_tmpxgpu = ROOT.RooExtendPdf("extpdf_tmpxgpu", "extpdf_tmpxgpu", tmpxgpu, npu)
    tmpfit = ROOT.RooAddPdf("tmpfit","tmpfit",ROOT.RooArgList(extpdf_tmpxg,extpdf_pu))
    
    #Depending on the pt binning, you may want to estimate the PU contribution differently. 
    #For now, always use the PU template smeared with a gaussian
    if binning == "jetpt":
        if pt is "20to25" or pt is "25to30" or  pt is "30to40" or pt is "40to50" :
            #tmpfit = ROOT.RooAddPdf("tmpfit","tmpfit",ROOT.RooArgList(extpdf_tmpxg))
            tmpfit = ROOT.RooAddPdf("tmpfit","tmpfit",ROOT.RooArgList(extpdf_pu,extpdf_tmpxg))
            #tmpfit = ROOT.RooAddPdf("tmpfit","tmpfit",ROOT.RooArgList(extpdf_pu,extpdf_tmpxlogn))
   
        else :
            tmpfit = ROOT.RooAddPdf("tmpfit","tmpfit",ROOT.RooArgList(extpdf_pu,extpdf_tmpxg)) 
            #tmpfit = ROOT.RooAddPdf("tmpfit","tmpfit",ROOT.RooArgList(extpdf_pu,extpdf_tmpxg))
            #tmpfit = ROOT.RooAddPdf("tmpfit","tmpfit",ROOT.RooArgList(extpdf_pu,extpdf_tmpxlogn))
            #tmpfit = ROOT.RooAddPdf("tmpfit","tmpfit",ROOT.RooArgList(extpdf_tmpxlogn)) 



    minX = 0.6
    maxX = 2.
    minxcorfact = 1
    ptfloat = float(pt.split("to")[0])
    minX =  12/ptfloat                                                                                                                                                                
    
    if binning == "zpt":
        if ptfloat < 100 : 
            maxX = 2
        else :
            maxX=  2
        '''
        if pt is "20to25" :
            minX, maxX = 0.6*minxcorfact, 2.
        elif pt is "25to30" :
            minX, maxX = 0.5*minxcorfact, 2.
        elif pt is "30to40" :
            minX, maxX = 0.4*minxcorfact, 2.
        elif pt is "40to50" :
            minX, maxX = 0.3*minxcorfact, 2.
        elif pt is "50to60" :
            minX, maxX = 0.25*minxcorfact, 2.
        elif pt is "60to85" :
            minX, maxX = 0.2*minxcorfact, 2.
        elif pt is "85to105" :
            minX, maxX = 0.15*minxcorfact, 2.
        elif pt is "105to130" :
            minX, maxX = 0.15*minxcorfact, 2.
        elif pt is "130to175" :
            minX, maxX = 0.15*minxcorfact, 2.
        #           gauss_sigma.setVal(0.11)
        #          gauss_sigma.setRange(0.05,0.13)
        #          gauss_mean.setVal(0.02)
        #         gauss_mean.setRange(-0.02,0.06)
        else :
        minX, maxX = 0.15, 2.
        '''
    if binning == "jetpt":
        minX = 0.4
        maxX = 2.5
        '''
        if pt is "20to25" :
            minX, maxX = 0.4*minxcorfact, 2.
        elif pt is "25to30" :
            minX, maxX = 0.4*minxcorfact, 2.
        elif pt is "30to40" :
            minX, maxX = 0.4*minxcorfact, 2.
        elif pt is "40to50" :
            minX, maxX = 0.6*minxcorfact, 2.
        elif pt is "50to60" :
            minX, maxX = 0.6*minxcorfact, 2.
        elif pt is "60to85" :
            minX, maxX = 0.6*minxcorfact, 2.
        elif pt is "85to105" :
            minX, maxX = 0.6*minxcorfact, 1.6
        elif pt is "105to130" :
            minX, maxX = 0.6*minxcorfact, 1.5
        elif pt is "130to175" :
            minX, maxX = 0.6*minxcorfact, 1.4
        else : 
            minX, maxX = 0.6*minxcorfact, 1.4
        '''
    tmpfit.fitTo(dh_ptbalance,ROOT.RooFit.Save(),ROOT.RooFit.SumW2Error(False) ,ROOT.RooFit.Range(minX,maxX))
   
    tmpfit.plotOn(frame)


    tmpfit.plotOn(frame, ROOT.RooFit.Components("extpdf_pu"),ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineColor(8))
    tmpfit.plotOn(frame, ROOT.RooFit.Components("extpdf_tmpxg"),ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineColor(2))
    #tmpfit.plotOn(frame, ROOT.RooFit.Components("extpdf_pli"),ROOT.RooFit.LineStyle(2))

    argset_fit = ROOT.RooArgSet(gauss_mean,gauss_sigma)
    #tmpxg.paramOn(frame,ROOT.RooFit.Format("NELU",ROOT.RooFit.AutoPrecision(1)),ROOT.RooFit.Layout(0.4,0.89,0.9)) 
    frame.SetMaximum(frame.GetMaximum()*1.25)
    frame.SetMinimum(0)

    # add chi2 info
    chi2_text = ROOT.TPaveText(0.15,0.72,0.15,0.88,"NBNDC")
    chi2_text.SetTextAlign(11)
    chi2_text.AddText("#chi^{2} fit = %s" %round(frame.chiSquare(6),2))
    chi2_text.AddText("#sigma_{gauss} "+"= {} #pm {}".format(round(gauss_sigma.getVal(),3), round(gauss_sigma.getError(),3)))
    chi2_text.AddText("#mu_{gauss} "+"= {} #pm {}".format(round(gauss_mean.getVal(),3), round(gauss_mean.getError(),3)) )
    #chi2_text.AddText("rms_{jer} "+"= {} #pm {}".format(round(jer.GetRMS(),3), round(jer.GetRMSError(),3)))
    chi2_text.SetTextSize(0.03)
    chi2_text.SetTextColor(2)
    chi2_text.SetShadowColor(0)
    chi2_text.SetFillColor(0)
    chi2_text.SetLineColor(0)
    frame.addObject(chi2_text)

    cfit = ROOT.TCanvas("cfit","cfit",600,600)
    cfit.SetLogx(False)
    frame.Draw()

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.3*cfit.GetTopMargin())
    latex2.SetTextFont(42)
    latex2.SetTextAlign(31) # align right                                                     

    if isData:        
        latex2.DrawLatex(0.90, 0.93,pt.split("to")[0]+" GeV < "+variable_name+" < "+pt.split("to")[1]+" GeV, "+eta.split("to")[0].replace("p",".")+" < |#eta_{jet}| < "+eta.split("to")[1].replace("p",".")+ ", Data")
    else:
        latex2.DrawLatex(0.90, 0.93,pt.split("to")[0]+" GeV < "+variable_name+" < "+pt.split("to")[1]+" GeV, "+eta.split("to")[0].replace("p",".")+" < |#eta_{jet}| < "+eta.split("to")[1].replace("p",".")+ ", MC")
    latex2.Draw("same")

    frame.Print()
 
    legend = ROOT.TLegend(0.60,0.75,0.88,0.88)
    legend.AddEntry(frame.findObject("tmpfit_Norm[balance]_Comp[extpdf_pu]"),"pile up","l")
    legend.AddEntry(frame.findObject("tmpfit_Norm[balance]_Comp[extpdf_tmpxg]"),"pli #otimes gauss","l")
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw("same")
 
    fit_filename = "fit_"+pt+"_"+eta
    #if not os.path.exists(fit_plot_directory): os.makedirs(fit_plot_directory)
    if isData:
        cfit.SaveAs(os.path.join(output, fit_filename+"_Data.pdf"))
        cfit.SaveAs(os.path.join(output, fit_filename+"_Data.png"))
    else:
        cfit.SaveAs(os.path.join(output, fit_filename+".pdf"))
        cfit.SaveAs(os.path.join(output, fit_filename+".png"))
    del cfit

    return gauss_sigma.getVal(), gauss_sigma.getError(), gauss_mean.getVal(), gauss_mean.getError(), frame.chiSquare(6)




def GetHisto (inputfile, histoprefix, histosuffix,  rebinfactor, theptbinning, ptbin, etabin, etaNegbin):
    '''
    This function returns a histogram taken from the 'inputfile' file that corresponds to a given pt/eta range
    The various arguments are: 
    histoprefix: string, prefix of the histo
    histosuffix: string, suffix of the histo
    rebinfactor: int, rebin factor
    theptbinning: string, what pt binning is used (zpt, jetpt)
    ptbin: string, ptrange considered
    etabin: string, eta range
    etaNegbin: string, eta range on the negative side
    '''
    #inputfile = ROOT.TFile(fname,"READ")
    
    folder            = theptbinning+"_"+ptbin.replace("to","_")

    ptbinlowedge = float(ptbin.split("to")[0]) 
    ptbinhighedge = float(ptbin.split("to")[1])

    etabinlowedge = float(etabin.replace("p",".").replace("min","-").split("to")[0]) 
    etabinhighedge = float(etabin.replace("p",".").replace("min","-").split("to")[1])


    #The histograms in the input file have a very fine eta binning. 
    #When ReducedBinning is true, a wider eta binning is used and the original histos are merged together.
    if ReducedBinning :
        #These correspond to the original eta bins
        _etaNeg= ["min5p191tomin4p013","min4p013tomin3p489","min3p489tomin3p139","min3p139tomin2p964","min2p964tomin2p853","min2p853tomin2p650","min2p650tomin2p500","min2p500tomin2p322","min2p322tomin2p043","min2p043tomin1p930","min1p930tomin1p740","min1p740tomin1p305","min1p305tomin1p131","min1p131tomin0p783","min0p783tomin0p522","min0p522to0p000"]      
        _eta= ["0p000to0p522","0p522to0p783","0p783to1p131","1p131to1p305","1p305to1p740","1p740to1p930","1p930to2p043","2p043to2p322","2p322to2p500","2p500to2p650","2p650to2p853","2p853to2p964","2p964to3p139","3p139to3p489","3p489to4p013","4p013to5p191"]                      
        _pt    = ["20to22","22to24","24to26","26to28","28to30","30to35","35to40","40to50","50to60","60to85","85to105","105to130","130to175","175to230","230to300","300to400"]
        histo = None

        for i in range(0, len(_pt)):
            localptbin_lowedge = float(_pt[i].split("to")[0])
            localptbin_highedge = float(_pt[i].split("to")[1])
            if localptbin_highedge <= ptbinlowedge : 
                continue 
            if localptbin_lowedge >= ptbinhighedge :
                continue 
            for j in range(0,len(_eta)):
                localetabin_lowedge = float(_eta[j].replace("p",".").split("to")[0])
                localetabin_highedge = float(_eta[j].replace("p",".").split("to")[1])
                if localetabin_highedge <= etabinlowedge :
                    continue
                if localetabin_lowedge >= etabinhighedge :
                    continue
                folder = theptbinning+"_"+_pt[i].replace("to","_")
                hname = folder +"/"+histoprefix +  "_pt"+_pt[i]+"_eta"+_eta[j]+histosuffix
                print hname 
                histotemp = inputfile.Get(hname).Clone()
                if histo == None :
                    histo = histotemp.Clone()
                    histo.SetName(histoprefix+ptbin+etabin+histosuffix)
                else :
                    histo.Add(histotemp)
            #Now negative eta side
            for j in range(0,len(_etaNeg)):
                localetabin_lowedge = float(_etaNeg[j].replace("p",".").replace("min","-").split("to")[0])
                localetabin_highedge = float(_etaNeg[j].replace("p",".").replace("min","-").split("to")[1])
                if localetabin_highedge <= etabinlowedge :
                    continue
                if localetabin_lowedge >= etabinhighedge :
                    continue
                folder = theptbinning+"_"+_pt[i].replace("to","_")
                hname = folder +"/"+histoprefix +  "_pt"+_pt[i]+"_eta"+_etaNeg[j]+histosuffix
                print hname 
                histotemp = inputfile.Get(hname).Clone()
                if histo == None :
                    histo = histotemp.Clone()
                    histo.SetName(histoprefix+ptbin+etabin+histosuffix)
                else :
                    histo.Add(histotemp)
         

                    #if hnamenegside != "":
                #histoNeg =  inputfile.Get(hnamenegside).Clone()
                #histo.Add(histoNeg)
                
        histo.Rebin(rebinfactor)
                
        return histo

    else :
        hname = histoprefix +  "_pt"+ptbin+"_eta"+etabin+histosuffix 
        hnamenegside = histoprefix +  "_pt"+ptbin+"_eta"+etaNegbin+histosuffix
        print hname," ", hnamenegside
        histo = inputfile.Get(hname).Clone()
        if hnamenegside != "":
            histoNeg =  inputfile.Get(hnamenegside).Clone()
            histo.Add(histoNeg)
            histo.Rebin(rebinfactor)
        return histo




def main():

    parser = argparse.ArgumentParser(description='Extract resolution')
    parser.add_argument("-o", "--output", dest="output", help="output folder name", type=str)
    parser.add_argument("-y", "--year",   dest="year",   help="data year", type=str)
    parser.add_argument("-c", "--channel",dest="channel",help="MuMu or EE", type=str)
    parser.add_argument("-b", "--binning",dest="binning",help="zpt or jetpt", type=str)
    parser.add_argument("-r", "--dorebin",dest="dorebin",help="true for rebinning, false by default", type=bool)
    args = parser.parse_args()    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetMarkerSize(0.5)
    ROOT.gStyle.SetOptLogx()

    #Rebin factor for the ptbalance distribution
    rebinfactor = 5
    laurentfolder = "/afs/cern.ch/work/l/lathomas/JERStudies/ForSoumya/"

    data_filename = ""
    mc_filename = ""

    if args.year == "UL2018":
        if args.channel == "MuMu":
            data_filename = laurentfolder+"output_DoubleMuonUL2018A_unctyOct2021.root"
            mc_filename    = laurentfolder+"output_MCUL2018ZtoMuMu_unctyOct2021.root"
        elif args.channel == "EE":
            data_filename = laurentfolder+""
            mc_filename    = laurentfolder+""
        else:
            data_filename = ""
            mc_filename   = ""

    else:
        print("Year does not exist")

    print "OPENING FILE", data_filename, "AND", mc_filename, args.binning
    f_data = ROOT.TFile(data_filename,"READ")
    f_mc   = ROOT.TFile(mc_filename,"READ")
    


    #Binning 
    _etaNeg= []
    _eta = []
    
    if ReducedBinning :
        _etaNeg = ["min5p191tomin2p964","min2p964tomin2p500","min2p500tomin1p305","min1p305to0p000"]
        _eta = ["0p000to1p305","1p305to2p500","2p500to2p964","2p964to5p191"]
        _pt    = ["20to22","22to24","24to26","26to28","28to30","30to35","35to40","40to50","50to60","60to85","85to105","105to130","130to175","175to230","230to300","300to400"]
    
    else: 
        _etaNeg= ["min5p191tomin4p013","min4p013tomin3p489","min3p489tomin3p139","min3p139tomin2p964","min2p964tomin2p853","min2p853tomin2p650","min2p650tomin2p500","min2p500tomin2p322","min2p322tomin2p043","min2p043tomin1p930","min1p930tomin1p740","min1p740tomin1p305","min1p305tomin1p131","min1p131tomin0p783","min0p783tomin0p522","min0p522to0p000"]      
        _eta= ["0p000to0p522","0p522to0p783","0p783to1p131","1p131to1p305","1p305to1p740","1p740to1p930","1p930to2p043","2p043to2p322","2p322to2p500","2p500to2p650","2p650to2p853","2p853to2p964","2p964to3p139","3p139to3p489","3p489to4p013","4p013to5p191"]                      
        #_pt    = ["20to25","25to30","30to40","40to50","50to60","60to85","85to105","105to130","130to175","175to230","230to300","300to400"]
        _pt    = ["20to22","22to24","24to26","26to28","28to30","30to35","35to40","40to50","50to60","60to85","85to105","105to130","130to175","175to230","230to300","300to400"]

    #Wide eta binning used to compute correction to the PU template
    _etabinsinBarrel = ["min1p305tomin1p131","min1p131tomin0p783","min0p783tomin0p522","min0p522to0p000","0p000to0p522","0p522to0p783","0p783to1p131","1p131to1p305"]
    _etabinsinEndcap1 = ["min2p500tomin2p322","min2p322tomin2p043","min2p043tomin1p930","min1p930tomin1p740","min1p740tomin1p305","1p305to1p740","1p740to1p930","1p930to2p043","2p043to2p322","2p322to2p500"]
    _etabinsinEndcap2 = ["min2p964tomin2p853","min2p853tomin2p650","min2p650tomin2p500","2p500to2p650","2p650to2p853","2p853to2p964"]
    _etabinsinHF = ["min5p191tomin4p013","min4p013tomin3p489","min3p489tomin3p139","min3p139tomin2p964","2p964to3p139","3p139to3p489","3p489to4p013","4p013to5p191"]


    ### first define bining
    xbins = []
    ybins = []
    for ptBin in _pt:
        xbins.append(float(ptBin.split("to")[0]))                                
        xbins.append(float(ptBin.split("to")[1]))                                
    for etaBin in _eta:
        ybins.append(float((etaBin.split("to")[0].replace("p",".")).replace("min","-")))               
        ybins.append(float((etaBin.split("to")[1].replace("p",".")).replace("min","-")))
    
    xbins.sort()
    ybins.sort()

    ## Transform to numpy array for ROOT 
    xbinsT = np.array(xbins)
    ybinsT = np.array(ybins)
    ## Just in case there are duplicates
    xbinsTab = np.unique(xbinsT)
    ybinsTab = np.unique(ybinsT)

    #Declaring histos    
    htitle = '#sigma_{jer}'
    hname  = 'jer'
    
    #JER (sigma of the Gaussian)  
    hdata = ROOT.TH2F(hname,htitle+" Data",xbinsTab.size-1,xbinsTab,ybinsTab.size-1,ybinsTab)
    hmc   = ROOT.TH2F(hname,htitle+" MC",xbinsTab.size-1,xbinsTab,ybinsTab.size-1,ybinsTab)
    #JES (mu of the Gaussian)
    hjec_data = ROOT.TH2F("hjec_data","#mu Data",xbinsTab.size-1,xbinsTab,ybinsTab.size-1,ybinsTab)
    hjec_mc   = ROOT.TH2F("hjec_mc","#mu MC",xbinsTab.size-1,xbinsTab,ybinsTab.size-1,ybinsTab)

    #Residuals (mu(data) - mu(MC)) 
    hresidual  = ROOT.TH2F("residual","residuals (#mu_{data} - #mu_{mc})",xbinsTab.size-1,xbinsTab,ybinsTab.size-1,ybinsTab)

    #JER and JES from a direct fit to pt(reco)/pt(gen)
    #sigma, mean from a Gaussian fit
    hsigmagauss_mctruth   = ROOT.TH2F("hsigmagauss_mctruth" ,"Gaussian sigma from pt(reco)/pt(gen)",xbinsTab.size-1,xbinsTab,ybinsTab.size-1,ybinsTab)
    hmeangauss_mctruth   = ROOT.TH2F("hmeangauss_mctruth" ,"Gaussian mean from pt(reco)/pt(gen) -1 ",xbinsTab.size-1,xbinsTab,ybinsTab.size-1,ybinsTab)
    #RMS, mean from the distribution directly
    hrms_mctruth   = ROOT.TH2F("hrms_mctruth" ,"RMS from pt(reco)/pt(gen)",xbinsTab.size-1,xbinsTab,ybinsTab.size-1,ybinsTab)
    hmean_mctruth   = ROOT.TH2F("hmean_mctruth" ,"< pt(reco)/pt(gen) - 1>",xbinsTab.size-1,xbinsTab,ybinsTab.size-1,ybinsTab)
    
    #Chi2/Nddof
    hchi2_mc = ROOT.TH2F("hchi2_mc","chi2 MC",xbinsTab.size-1,xbinsTab,ybinsTab.size-1,ybinsTab)
    hchi2_data   = ROOT.TH2F("hchi2_data","chi2 Data",xbinsTab.size-1,xbinsTab,ybinsTab.size-1,ybinsTab)
    
    hdata.Sumw2()
    hmc.Sumw2()
    hjec_data.Sumw2()
    hjec_mc.Sumw2()
    hresidual.Sumw2()
    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat("4.2f")
    

    
    #Loop over all pt bins
    for i in range(0,len(_pt)):
        print _pt[i]
        
        #next few lines needed to compute MC corrections to PU template (from PU gen)
        folder            = args.binning+"_"+_pt[i].replace("to","_")+"/"
        #_histos_PUtrueMC => this is the pt balance distribution obtained from jets not matched to a gen jet
        _histos_PUtrueMC = [ 
            (f_mc.Get("zpt_40_50/h_PtRecoJetoverPtRecoZ_NoGenMatchedJet_pt40to50_eta0p000to0p522_zptbinning")).Clone() ,
            (f_mc.Get("zpt_40_50/h_PtRecoJetoverPtRecoZ_NoGenMatchedJet_pt40to50_eta0p000to0p522_zptbinning")).Clone() ,
            (f_mc.Get("zpt_40_50/h_PtRecoJetoverPtRecoZ_NoGenMatchedJet_pt40to50_eta0p000to0p522_zptbinning")).Clone() ,
            (f_mc.Get("zpt_40_50/h_PtRecoJetoverPtRecoZ_NoGenMatchedJet_pt40to50_eta0p000to0p522_zptbinning")).Clone() 
        ]
        #_histos_PUtemplateMC => this is the PU template computed in the same way as in data
        _histos_PUtemplateMC = [ 
            (f_mc.Get("zpt_40_50/h_PtRecoJetoverPtRecoZ_PU_pt40to50_eta0p000to0p522_zptbinning")).Clone(),
            (f_mc.Get("zpt_40_50/h_PtRecoJetoverPtRecoZ_PU_pt40to50_eta0p000to0p522_zptbinning")).Clone(),
            (f_mc.Get("zpt_40_50/h_PtRecoJetoverPtRecoZ_PU_pt40to50_eta0p000to0p522_zptbinning")).Clone(),
            (f_mc.Get("zpt_40_50/h_PtRecoJetoverPtRecoZ_PU_pt40to50_eta0p000to0p522_zptbinning")).Clone()
        ]
        
        
        suffixDetRegion = ["HB","HE1","HE2","HF"]
        print "_histos_PUtrueMC entries", _histos_PUtrueMC[0].GetEntries()
        for k in range(0,len(_histos_PUtrueMC)):
            _histos_PUtrueMC[k].SetName("_histos_PUtrueMC_"+_pt[i]+suffixDetRegion[k])
            _histos_PUtrueMC[k].Reset()
            _histos_PUtemplateMC[k].SetName("_histos_PUtemplateMC_"+_pt[i]+suffixDetRegion[k])
            _histos_PUtemplateMC[k].Reset()

    
        #Get the histos
        if ReducedBinning:
            for j in range(0,len(_eta)):
                h_ptbal_unmatched = GetHisto (f_mc,
                                              "h_PtRecoJetoverPtRecoZ_NoGenMatchedJet","_"+args.binning+"binning",
                                              rebinfactor,args.binning,_pt[i], _eta[j],_etaNeg[len(_eta)-1-j])
                
                h_ptbal_pu        = GetHisto (f_mc,
                                              "h_PtRecoJetoverPtRecoZ_PU","_"+args.binning+"binning",
                                              rebinfactor,args.binning,_pt[i], _eta[j],_etaNeg[len(_eta)-1-j])
                
                if _eta[j] == "0p000to1p305":
                    _histos_PUtrueMC[0].Add(h_ptbal_unmatched)
                    _histos_PUtemplateMC[0].Add(h_ptbal_pu)
                elif _eta[j] == "1p305to2p500":
                    _histos_PUtrueMC[1].Add(h_ptbal_unmatched)
                    _histos_PUtemplateMC[1].Add(h_ptbal_pu)
                elif _eta[j] == "2p500to2p964":
                    _histos_PUtrueMC[2].Add(h_ptbal_unmatched)
                    _histos_PUtemplateMC[2].Add(h_ptbal_pu)
                elif _eta[j] == "2p964to5p191":
                    _histos_PUtrueMC[3].Add(h_ptbal_unmatched)
                    _histos_PUtemplateMC[3].Add(h_ptbal_pu)

        else: 
            for j in range(0,len(_eta)):
                print(folder+"h_PtRecoJetoverPtRecoZ_PU_pt"+_pt[i]+"_eta"+_eta[j]+"_"+args.binning+"binning")
                h_ptbal_pu        = (f_mc.Get(folder+"h_PtRecoJetoverPtRecoZ_PU_pt"+_pt[i]+"_eta"+_eta[j]+"_"+args.binning+"binning")).Clone()
                h_ptbal_unmatched = (f_mc.Get(folder+"h_PtRecoJetoverPtRecoZ_NoGenMatchedJet_pt"+_pt[i]+"_eta"+_eta[j]+"_"+args.binning+"binning")).Clone()
                hNeg_ptbal_pu        = (f_mc.Get(folder+"h_PtRecoJetoverPtRecoZ_PU_pt"+_pt[i]+"_eta"+_etaNeg[len(_eta)-1-j]+"_"+args.binning+"binning")).Clone()
                hNeg_ptbal_unmatched = (f_mc.Get(folder+"h_PtRecoJetoverPtRecoZ_NoGenMatchedJet_pt"+_pt[i]+"_eta"+_etaNeg[len(_eta)-1-j]+"_"+args.binning+"binning")).Clone()
                h_ptbal_pu.Add(hNeg_ptbal_pu)
                h_ptbal_unmatched.Add(hNeg_ptbal_unmatched)
                for etabinstr in _etabinsinBarrel:
                    if _eta[j] == etabinstr : 
                        _histos_PUtrueMC[0].Add(h_ptbal_unmatched)
                        _histos_PUtemplateMC[0].Add(h_ptbal_pu)
                        
                for etabinstr in _etabinsinEndcap1:
                    if _eta[j] == etabinstr :
                        _histos_PUtrueMC[1].Add(h_ptbal_unmatched)
                        _histos_PUtemplateMC[1].Add(h_ptbal_pu)

                for etabinstr in _etabinsinEndcap2:
                    if _eta[j] == etabinstr :
                        _histos_PUtrueMC[2].Add(h_ptbal_unmatched)
                        _histos_PUtemplateMC[2].Add(h_ptbal_pu)
                        
                for etabinstr in _etabinsinHF:
                    if _eta[j] == etabinstr :
                        _histos_PUtrueMC[3].Add(h_ptbal_unmatched)
                        _histos_PUtemplateMC[3].Add(h_ptbal_pu)
        
        thefitfn = [None]* 4
        for k in range(0,len(_histos_PUtrueMC)):
            print "entries is " 
            print _histos_PUtrueMC[k].GetEntries()
            print _histos_PUtemplateMC[k].GetEntries()
            print _histos_PUtrueMC[k].GetNbinsX()
            _histos_PUtrueMC[k].Rebin(10)
            _histos_PUtemplateMC[k].Rebin(10)
            thefitfn[k] = ROOT.TF1("myf","expo",0.2,2.)
            thefitfn[k].SetParLimits(1,-10,-0.0001)

            if _histos_PUtrueMC[k].GetEntries() <= 10 or _histos_PUtemplateMC[k].GetEntries() <=10 : 
                thefitfn[k].SetParameter(0,0)
                thefitfn[k].SetParameter(1,0)
                print "This is zero "
                continue

            _histos_PUtrueMC[k].Scale(1./_histos_PUtrueMC[k].Integral())
            _histos_PUtemplateMC[k].Scale(1./_histos_PUtemplateMC[k].Integral())
            _histos_PUtrueMC[k].Divide(_histos_PUtemplateMC[k])
            _histos_PUtrueMC[k].Fit(thefitfn[k], "","")
            c_putemplatecorr = ROOT.TCanvas("c_putemplatecorr","c_putemplatecorr",600,600)
            c_putemplatecorr.SetLogx(False)
            _histos_PUtrueMC[k].GetYaxis().SetRangeUser(0.,1.5)
            _histos_PUtrueMC[k].Draw()
            thefitfn[k].Draw("same")
            if k == 0:
                c_putemplatecorr.SaveAs(os.path.join(args.output, "PUtemplatecorr"+_pt[i]+"_barrel.pdf"))
                c_putemplatecorr.SaveAs(os.path.join(args.output, "PUtemplatecorr"+_pt[i]+"_barrel.png"))
            if k == 1:
                c_putemplatecorr.SaveAs(os.path.join(args.output, "PUtemplatecorr"+_pt[i]+"_endcap1.pdf"))
                c_putemplatecorr.SaveAs(os.path.join(args.output, "PUtemplatecorr"+_pt[i]+"_endcap1.png"))
            if k == 2:
                c_putemplatecorr.SaveAs(os.path.join(args.output, "PUtemplatecorr"+_pt[i]+"_endcap2.pdf"))
                c_putemplatecorr.SaveAs(os.path.join(args.output, "PUtemplatecorr"+_pt[i]+"_endcap2.png"))
            if k == 3:
                c_putemplatecorr.SaveAs(os.path.join(args.output, "PUtemplatecorr"+_pt[i]+"_hf.pdf"))
                c_putemplatecorr.SaveAs(os.path.join(args.output, "PUtemplatecorr"+_pt[i]+"_hf.png"))
            
            #c_putemplatecorr.SaveAs("PUtemplatecorr.root")
            print "ok" 
        
           
        for j in range(0,len(_eta)):
            
            #folder            = args.binning+"_"+_pt[i].replace("to","_") 

            
            h_ptbal           = GetHisto (f_mc,
                                          folder+"h_PtRecoJetoverPtRecoZ_testsample","_"+args.binning+"binning",
                                          rebinfactor,args.binning,_pt[i], _eta[j],_etaNeg[len(_eta)-1-j])
                                          #correctPUtemplate)
            
            h_ptbal_pu        = GetHisto (f_mc,
                                          folder+"h_PtRecoJetoverPtRecoZ_PU","_"+args.binning+"binning",
                                          rebinfactor,args.binning,_pt[i], _eta[j],_etaNeg[len(_eta)-1-j])
                                          #correctPUtemplate)

            h_ptbal_matched   = GetHisto (f_mc,
                                          folder+"h_PtRecoJetoverPtRecoZ_GenMatchedJet","_"+args.binning+"binning", 
                                          rebinfactor,args.binning,_pt[i], _eta[j],_etaNeg[len(_eta)-1-j])
                                          #correctPUtemplate)

            h_ptbal_unmatched = GetHisto (f_mc,
                                          folder+"h_PtRecoJetoverPtRecoZ_NoGenMatchedJet","_"+args.binning+"binning",
                                          rebinfactor,args.binning,_pt[i], _eta[j],_etaNeg[len(_eta)-1-j])
                                          #correctPUtemplate)
                                          
            h_pli             = GetHisto (f_mc,
                                          folder+"h_PtGenJetoverPtGenZ","_"+args.binning+"binning",
                                          rebinfactor,args.binning,_pt[i], _eta[j],_etaNeg[len(_eta)-1-j])
                                          #correctPUtemplate)
            
            #JER from ptreco/ptgen (all PU bins together)
            foldergen = "genjetpt_"+_pt[i].replace("to","_")+"/"
            h_jer         = GetHisto (f_mc,
                                      foldergen+"h_PtRecoJetoverPtGenJet_GentPtBinning","_pu20to40",
                                      rebinfactor,"genjetpt",_pt[i], _eta[j],_etaNeg[len(_eta)-1-j])
                                      #correctPUtemplate)

            h_jer_2         = GetHisto (f_mc,
                                      foldergen+"h_PtRecoJetoverPtGenJet_GentPtBinning","_pult20",
                                      rebinfactor,"genjetpt",_pt[i], _eta[j],_etaNeg[len(_eta)-1-j])
                                      #correctPUtemplate)

            h_jer_3         = GetHisto (f_mc,
                                      foldergen+"h_PtRecoJetoverPtGenJet_GentPtBinning","_pugt40",
                                      rebinfactor,"genjetpt",_pt[i], _eta[j],_etaNeg[len(_eta)-1-j])
                                      #correctPUtemplate)
            h_jer.Add(h_jer_2)
            h_jer.Add(h_jer_3)
            
            h_ptbaldata       = GetHisto (f_data,
                                          folder+"h_PtRecoJetoverPtRecoZ_testsample","_"+args.binning+"binning",
                                          rebinfactor,args.binning,_pt[i], _eta[j],_etaNeg[len(_eta)-1-j])
                                          #correctPUtemplate)

            h_ptbal_pudata    = GetHisto (f_data,
                                          folder+"h_PtRecoJetoverPtRecoZ_PU","_"+args.binning+"binning",
                                          rebinfactor,args.binning,_pt[i], _eta[j],_etaNeg[len(_eta)-1-j])
                                          #correctPUtemplate)


            

            if correctforMCPU and args.binning=="zpt":
                for etabinstr in _etabinsinBarrel:
                    if _eta[j] == etabinstr or _eta[j] =="0p000to1p305":
                        h_ptbal_pu.Multiply(thefitfn[0])
                        h_ptbal_pudata.Multiply(thefitfn[0])
                        
                for etabinstr in _etabinsinEndcap1:
                    if _eta[j] == etabinstr or _eta[j] =="1p305to2p500":
                        h_ptbal_pu.Multiply(thefitfn[1])
                        h_ptbal_pudata.Multiply(thefitfn[1])

                for etabinstr in _etabinsinEndcap2:
                    if _eta[j] == etabinstr or _eta[j] =="2p500to2p964":
                        h_ptbal_pu.Multiply(thefitfn[2])
                        h_ptbal_pudata.Multiply(thefitfn[2])

                for etabinstr in _etabinsinHF:
                    if _eta[j] == etabinstr or _eta[j] =="2p964to5p191":
                        h_ptbal_pu.Multiply(thefitfn[3])
                        h_ptbal_pudata.Multiply(thefitfn[3])
            

            sigma, sigma_error, mean, mean_error, chi2 = ConvFit(h_ptbal, h_ptbal_pu, h_pli, h_jer, args.output,_pt[i], _eta[j], args.binning)
            sigma_data, sigma_error_data, mean_data, mean_error_data, chi2_data = ConvFit(h_ptbaldata, h_ptbal_pudata, h_pli, h_jer, args.output,_pt[i], _eta[j], args.binning,isData=True)

            #Some ideas to improve the fit: test various templates, only look at low PU events, change fit range, ...
            #sigma_data, sigma_error_data, mean_data, mean_error_data, chi2_data = ConvFit(h_ptbaldata, h_ptbal_unmatched, h_pli, h_jer, args.output,_pt[i], _eta[j], args.binning,isData=True)
            #sigma, sigma_error, mean, mean_error = ConvFit(h_ptbal, h_ptbal_unmatched, h_pli, h_jer, args.output,_pt[i], _eta[j], args.binning)
            #sigma_data, sigma_error_data, mean_data, mean_error_data = ConvFit(h_ptbaldata, h_ptbal_unmatched, h_pli, h_jer, args.output,_pt[i], _eta[j], args.binning,isData=True)
            #Using PU template from MC (h_ptbal_pu instead of h_ptbal_pudata) !!!!
            #sigma_data, sigma_error_data, mean_data, mean_error_data = ConvFit(h_ptbaldata, h_ptbal_pu, h_pli, h_jer, args.output,_pt[i], _eta[j], args.binning,isData=True)

            thefitresponse = ROOT.TF1("myf","gaus",0.6,1.5)
            if h_jer.Integral()>0 :
                h_jer.Scale(1./h_jer.Integral())
            h_jer.Fit(thefitresponse,"R","")
            c_genresol = ROOT.TCanvas("c_genresol","c_genresol",600,600)
            c_genresol.SetLogx(False)
            h_jer.Draw()
            thefitresponse.Draw("same")
            c_genresol.SaveAs(os.path.join(args.output, "genresol"+_pt[i]+ "_"+_eta[j]+".pdf"))
            c_genresol.SaveAs(os.path.join(args.output, "genresol"+_pt[i]+ "_"+_eta[j]+".png"))
            
            hsigmagauss_mctruth.SetBinContent(i+1,j+1, round(float(thefitresponse.GetParameter(2)),4))
            hsigmagauss_mctruth.SetBinError(i+1,j+1, round(float(thefitresponse.GetParError(2)),4))
            hmeangauss_mctruth.SetBinContent(i+1,j+1, round(float(thefitresponse.GetParameter(1))-1,4))
            hmeangauss_mctruth.SetBinError(i+1,j+1, round(float(thefitresponse.GetParError(1)),4))

            hmean_mctruth.SetBinContent(i+1,j+1, round(float(h_jer.GetMean())-1,4))
            hmean_mctruth.SetBinError(i+1,j+1, round(float(h_jer.GetMeanError()),4))
            hrms_mctruth.SetBinContent(i+1,j+1, round(float(h_jer.GetRMS()),4))
            hrms_mctruth.SetBinError(i+1,j+1, round(float(h_jer.GetRMSError()),4))

            #if sigma<0.001:
            #    print "Something went very WRONG! Taking RMS of the h_res"
            #    sigma = h_jer.GetRMS()
            #    sigma_error = h_jer.GetRMSError()
                    
            hmc.SetBinContent(i+1,j+1, round(float(sigma/(1+mean)),4))
            hmc.SetBinError  (i+1,j+1, round(float(sigma_error/(1+mean)),4))

            hdata.SetBinContent(i+1,j+1, round(float(sigma_data/(1+mean_data)),4))
            hdata.SetBinError  (i+1,j+1, round(float(sigma_error_data/(1+mean_data)),4))
    
            hjec_mc.SetBinContent(i+1,j+1, round(float(mean),4))
            hjec_mc.SetBinError(i+1,j+1, round(float(mean_error),4))

            hjec_data.SetBinContent(i+1,j+1, round(float(mean_data),4))
            hjec_data.SetBinError(i+1,j+1, round(float(mean_error_data),4))
            
            hresidual.SetBinContent(i+1,j+1, round(float(mean_data)-float(mean),4))
            #hresidual.SetBinError(i+1,j+1, round(float(mean_error_data)-float(mean_error),4))
            hresidual.SetBinError(i+1,j+1, round(math.sqrt(float(mean_error_data)*float(mean_error_data)+float(mean_error)*float(mean_error)),4))
            
            hchi2_data.SetBinContent(i+1,j+1, round(float(chi2_data),2))
            hchi2_mc.SetBinContent(i+1,j+1, round(float(chi2),2))
            
    if args.binning == "zpt":
        xname = "Z p_{T} [GeV]"
    else:
        xname = "Jet p_{T} [GeV]"

    hmc.GetYaxis().SetTitle("Jet #eta")    
    hmc.GetXaxis().SetTitle(xname)
    hmc.GetXaxis().SetMoreLogLabels()

    hdata.GetYaxis().SetTitle("Jet #eta")
    hdata.GetXaxis().SetTitle(xname)
    hdata.GetXaxis().SetMoreLogLabels()

    hresidual.GetYaxis().SetTitle("Jet #eta")
    hresidual.GetXaxis().SetTitle(xname)
    hresidual.GetXaxis().SetMoreLogLabels()

    hjec_mc.GetYaxis().SetTitle("Jet #eta")
    hjec_mc.GetXaxis().SetTitle(xname)
    hjec_mc.GetXaxis().SetMoreLogLabels()

    hjec_data.GetYaxis().SetTitle("Jet #eta")
    hjec_data.GetXaxis().SetTitle(xname)
    hjec_data.GetXaxis().SetMoreLogLabels()
    
    hsigmagauss_mctruth.GetYaxis().SetTitle("Jet #eta")
    hsigmagauss_mctruth.GetXaxis().SetTitle(xname)
    hsigmagauss_mctruth.GetXaxis().SetMoreLogLabels()
    hmeangauss_mctruth.GetYaxis().SetTitle("Jet #eta")
    hmeangauss_mctruth.GetXaxis().SetTitle(xname)
    hmeangauss_mctruth.GetXaxis().SetMoreLogLabels()

    hrms_mctruth.GetYaxis().SetTitle("Jet #eta")
    hrms_mctruth.GetXaxis().SetTitle(xname)
    hrms_mctruth.GetXaxis().SetMoreLogLabels()
    hmean_mctruth.GetYaxis().SetTitle("Jet #eta")
    hmean_mctruth.GetXaxis().SetTitle(xname)
    hmean_mctruth.GetXaxis().SetMoreLogLabels()

    hchi2_mc .GetYaxis().SetTitle("Jet #eta")
    hchi2_mc .GetXaxis().SetTitle(xname)
    hchi2_mc .GetXaxis().SetMoreLogLabels()
    hchi2_data.GetYaxis().SetTitle("Jet #eta")
    hchi2_data.GetXaxis().SetTitle(xname)
    hchi2_data.GetXaxis().SetMoreLogLabels()


    c2 = ROOT.TCanvas("c2","c2",600,600)
    hmc.SetMinimum(0.0)
    hmc.SetMaximum(0.5)
    hmc.Draw("colztexterr")
    c2.SaveAs(os.path.join(args.output, "h2_jer_mc"+".pdf"))
    c2.SaveAs(os.path.join(args.output, "h2_jer_mc"+".png"))

    c3 = ROOT.TCanvas("c3","c3",600,600)
    hdata.SetMinimum(0.0)
    hdata.SetMaximum(0.5)
    hdata.Draw("colztexterr")
    c3.SaveAs(os.path.join(args.output, "h2_jer_data"+".pdf"))
    c3.SaveAs(os.path.join(args.output, "h2_jer_data"+".png"))

    c4 = ROOT.TCanvas("c4","c4",600,600)
    hdata.Divide(hmc)
    hdata.SetNameTitle("jersf","#sigma_{jer} SF")
    hdata.SetMaximum(1.5)
    hdata.SetMinimum(0.5)
    hdata.Draw("colztexterr")
    c4.SaveAs(os.path.join(args.output, "h2_jer_sf"+".pdf"))
    c4.SaveAs(os.path.join(args.output, "h2_jer_sf"+".png"))
    c4.SaveAs(os.path.join(args.output, "h2_jer_sf"+".root"))
    hdata.SaveAs(os.path.join(args.output, "th2_jer_sf"+".root"))

    c5 = ROOT.TCanvas("c5","c5",600,600)
    hresidual.SetMinimum(-0.1)
    hresidual.SetMaximum(0.1)
    hresidual.Draw("colztexterr")
    c5.SaveAs(os.path.join(args.output, "h2_residual"+".pdf"))
    c5.SaveAs(os.path.join(args.output, "h2_residual"+".png"))
    c5.SaveAs(os.path.join(args.output, "h2_residual"+".root"))
    hresidual.SaveAs(os.path.join(args.output, "th2_residual"+".root"))
    c6 = ROOT.TCanvas("c6","c6",600,600)
    hjec_mc.SetMinimum(-0.05)
    hjec_mc.SetMaximum(0.15)
    hjec_mc.Draw("colztexterr")
    c6.SaveAs(os.path.join(args.output, "h2_jec_mc"+".pdf"))
    c6.SaveAs(os.path.join(args.output, "h2_jec_mc"+".png"))

    c7 = ROOT.TCanvas("c7","c7",600,600)
    hjec_data.SetMinimum(-0.05)
    hjec_data.SetMaximum(0.15)
    hjec_data.Draw("colztexterr")
    c7.SaveAs(os.path.join(args.output, "h2_jec_data"+".pdf"))
    c7.SaveAs(os.path.join(args.output, "h2_jec_data"+".png"))
    
    c8 = ROOT.TCanvas("c8","c8",600,600)
    hmean_mctruth.SetMinimum(-0.05)
    hmean_mctruth.SetMaximum(0.15)
    hmean_mctruth.Draw("colztexterr")
    c8.SaveAs(os.path.join(args.output, "h2_mean_mctruth"+".pdf"))
    c8.SaveAs(os.path.join(args.output, "h2_mean_mctruth"+".png"))

    c9 = ROOT.TCanvas("c9","c9",600,600)
    hrms_mctruth.SetMinimum(0)
    hrms_mctruth.SetMaximum(0.5)
    hrms_mctruth.Draw("colztexterr")
    c9.SaveAs(os.path.join(args.output, "h2_rms_mctruth"+".pdf"))
    c9.SaveAs(os.path.join(args.output, "h2_rms_mctruth"+".png"))


    c10 = ROOT.TCanvas("c10","c10",600,600)
    hmeangauss_mctruth.SetMinimum(-0.05)
    hmeangauss_mctruth.SetMaximum(0.15)
    hmeangauss_mctruth.Draw("colztexterr")
    c10.SaveAs(os.path.join(args.output, "h2_meangauss_mctruth"+".pdf"))
    c10.SaveAs(os.path.join(args.output, "h2_meangauss_mctruth"+".png"))

    c11 = ROOT.TCanvas("c11","c11",600,600)
    hsigmagauss_mctruth.SetMinimum(0)
    hsigmagauss_mctruth.SetMaximum(0.5)
    hsigmagauss_mctruth.Draw("colztexterr")
    c11.SaveAs(os.path.join(args.output, "h2_sigmagauss_mctruth"+".pdf"))
    c11.SaveAs(os.path.join(args.output, "h2_sigmagauss_mctruth"+".png"))
 
    c12 = ROOT.TCanvas("c12","c12",600,600)
    hchi2_mc.SetMinimum(0.1)
    hchi2_mc.SetMaximum(50)
    hchi2_mc.Draw("colztexterr")
    c12.SetLogz()
    c12.SaveAs(os.path.join(args.output, "h2_chi2_mc"+".pdf"))
    c12.SaveAs(os.path.join(args.output, "h2_chi2_mc"+".png"))
    
    c13 = ROOT.TCanvas("c13","c13",600,600)
    hchi2_data.SetMinimum(0.1)
    hchi2_data.SetMaximum(50)
    hchi2_data.Draw("colztexterr")
    c13.SetLogz()
    c13.SaveAs(os.path.join(args.output, "h2_chi2_data"+".pdf"))
    c13.SaveAs(os.path.join(args.output, "h2_chi2_data"+".png"))

    del c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13
    
if __name__ == '__main__':
    main()
