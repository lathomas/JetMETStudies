import ROOT
import os
import sys
import argparse

colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kOrange, ROOT.kMagenta, ROOT.kGreen+2, ROOT.kGray+1, ROOT.kCyan+2, ROOT.kYellow+2, ROOT.kOrange+2]
dirname = 'plotL1Run3_prov2/'


def main():
    parser = argparse.ArgumentParser(
        description='''Plotter''',
        usage='use "%(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-t", "--type", dest="type", help="Type of plots (hists, eff)", type=str, default='')
    parser.add_argument("-i", "--input", dest="inputFiles", help="Input file", nargs='+', type=str, default='')
    parser.add_argument("--den", dest="den", help="Name(s) of the histo denominator(s). There can be 1 or the same number as numerator histos ", nargs='+', type=str, default='')
    parser.add_argument("--num", dest="num", help="Name(s) of the histo numerator(s)", nargs='+', type=str, default='')
    parser.add_argument("--xtitle", dest="xtitle", help="X axis title", type=str, default='p_{T} (GeV)')
    parser.add_argument("--ytitle", dest="ytitle", help="Y axis title", type=str, default='Number of entries')
    parser.add_argument("--ztitle", dest="ztitle", help="Z axis title", type=str, default='Number of entries')
    parser.add_argument("--legend", dest="legendlabels", help="Legend labels", nargs='+', type=str, default='')
    parser.add_argument("--extralabel", dest="extralabel", help="Extra label", type=str, default='')
    parser.add_argument("--setlogx", dest="setlogx", help="Set log x", type=bool, default=False)
    parser.add_argument("--plotname", dest="plotname", help="Name of the plot", type=str, default='plot')
    parser.add_argument("--h2d", dest="h2d", help="2D histo (for profile etc)", nargs='+',type=str)
    parser.add_argument("--axisranges", dest="axisranges", help="Axis ranges [xmin, xmax, ymin, ymax, zmin, zmax]", nargs='+', type=float, default=[])
    parser.add_argument("--addnumtoden", dest="addnumtoden", help="Add numerator histo to denominator (because it only contains complementary events e.g. failing probes)",type=bool, default=False)
    parser.add_argument("--saveplot", dest="saveplot", help="Save plots or not",type=bool, default = False)
    parser.add_argument("--interactive", dest="interactive", help="Run in interactive mode (keep plot drawn)", type=bool, default=False)
    parser.add_argument("--suffix_files", dest="suffix_files", help="Input files suffix", nargs='+', type=str, default='')
        
    args = parser.parse_args()
    
    if args.interactive == False:
        ROOT.gROOT.SetBatch(1)
        
    inputFiles = []
    for i in args.inputFiles:
        inputFiles.append(ROOT.TFile(i,"read"))


    
    if args.type=='efficiency':
        h_dens = [] 
        h_nums = []
        for inputFile in inputFiles:
            print(inputFile.GetName())
            if len(args.den) !=1 and len(args.den)!=len(args.num):
                print("Numerator has {} histos while denominator has {}. Exiting.".format(len(h_nums),len(h_dens)))
                return 
            for i in range(len(args.num)):
                h_nums.append(inputFile.Get(args.num[i]).Clone())
                if len(args.den)==1:
                    h_dens.append(inputFile.Get(args.den[0]).Clone())
                else:
                    h_dens.append(inputFile.Get(args.den[i]).Clone())
                h_nums[i].SetName(h_nums[i].GetName()+"_{}".format(i))
                h_dens[i].SetName(h_dens[i].GetName()+"_{}".format(i))
        effs = compute_eff(h_dens, h_nums, args.addnumtoden)
        drawplots(effs, legendlabels = args.legendlabels, xtitle=args.xtitle, ytitle=args.ytitle, ztitle=args.ztitle, extralabel=args.extralabel, setlogx=args.setlogx, plotname=args.plotname, axisranges=args.axisranges, saveplot = args.saveplot, interactive=args.interactive, suffix_files = args.suffix_files)
            
    if args.type=='profilex_fromh2':
        h2ds = []
        for i in args.h2d:
            h2ds.append(inputFiles[0].Get(i).Clone())
        profiles = compute_profilex(h2ds)
        drawplots(profiles, legendlabels = args.legendlabels, xtitle=args.xtitle, ytitle=args.ytitle, ztitle=args.ztitle, extralabel=args.extralabel, setlogx=args.setlogx, plotname=args.plotname, axisranges=args.axisranges, saveplot = args.saveplot, interactive=args.interactive)

    if args.type=='resolvsx':
        h2ds = []
        for inputFile in inputFiles:
            for i in range(len(args.h2d)):
                h2ds.append(inputFile.Get(args.h2d[i]).Clone())
                h2ds[i].SetName(h2ds[i].GetName()+"_{}".format(i))
        hresponse, hresol = compute_ResolutionvsX(h2ds)
        drawplots(hresponse, legendlabels = args.legendlabels, xtitle=args.xtitle, ytitle='#mu'+args.ytitle, ztitle=args.ztitle, extralabel=args.extralabel, setlogx=args.setlogx, plotname='mu_'+args.plotname, axisranges=args.axisranges, saveplot = args.saveplot, interactive=args.interactive, suffix_files = args.suffix_files)
        axisranges = args.axisranges
        axisranges[2] = 0
        axisranges[3] = 0.5
        drawplots(hresol, legendlabels = args.legendlabels, xtitle=args.xtitle, ytitle='#sigma_{scale corr.}'+args.ytitle, ztitle=args.ztitle, extralabel=args.extralabel, setlogx=args.setlogx, plotname='resol_'+args.plotname, axisranges=axisranges, saveplot = args.saveplot, interactive=args.interactive)

        
        

def canvas():
    c = ROOT.TCanvas("c_eff","c_eff",700,600)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.15)
    return c


def drawplots(objs, legendlabels, xtitle='', ytitle='', ztitle='',  extralabel='', setlogx=False, plotname='plot', axisranges=[], saveplot=False, interactive=False, suffix_files = []):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTextFont(42)
    c = canvas()

    labelsize = len(legendlabels)
    for suffix in suffix_files:
        for j in range(labelsize):
            legendlabels.append(legendlabels[j] + suffix)

    if len(suffix_files)>0:
        del legendlabels[0:labelsize]
    
    if len(legendlabels) != len(objs):
        print("Some histos are missing a legend. ")
        legendlabels=[]
        for i in range(len(objs)):
            legendlabels.append('')
            

    
    legend = ROOT.TLegend(0.6,0.12,0.85,0.12+0.04*len(objs),"","mlNDC")
    if plotname.find('resol')>=0:
        legend.SetY1(0.5)
        legend.SetY2(0.5+0.04*len(objs))
        
    legend.SetTextFont(42)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
        
    objs[0].SetTitle(";{};{};{}".format(xtitle,ytitle,ztitle))
    if objs[0].GetDimension () == 1:
        objs[0].Draw()
    else:
        objs[0].Draw("zcol")

    if setlogx:
        c.SetLogx()

    c.Update()
    drawnobject = None
    if type(objs[0])==ROOT.TEfficiency and objs[0].GetDimension () == 1:
        drawnobject = objs[0].GetPaintedGraph()
    elif type(objs[0])==ROOT.TEfficiency and objs[0].GetDimension () == 2:
        drawnobject = objs[0].GetPaintedHistogram()
    else:
        drawnobject = objs[0]
        
    drawnobject.GetXaxis().SetMoreLogLabels()
    drawnobject.GetXaxis().SetNoExponent()
    drawnobject.GetXaxis().SetTitleSize(0.04)
    drawnobject.GetYaxis().SetTitleSize(0.04)
    if objs[0].GetDimension() >=2:
        drawnobject.GetZaxis().SetTitleSize(0.04)
    if len(axisranges)>=2:
        drawnobject.GetXaxis().SetRangeUser(axisranges[0],axisranges[1])
    if len(axisranges)>=4:
        drawnobject.GetYaxis().SetRangeUser(axisranges[2],axisranges[3])
    if len(axisranges)>=6:
        drawnobject.GetZaxis().SetRangeUser(axisranges[4],axisranges[5])
    if len(axisranges)<4 and objs[0].GetDimension() == 1 and type(objs[0]) == ROOT.TEfficiency : 
        drawnobject.GetYaxis().SetRangeUser(0., 1.3) 
        
    print(type(objs))
    for i, h in enumerate(objs):
        if legendlabels[i] != '':
            legend.AddEntry(h,legendlabels[i],"lep")
        if i>0:
            objs[i].Draw("same")


    label_cms = ROOT.TLatex()
    label_cms.SetTextSize(0.05)
    label_cms.DrawLatexNDC(0.15, 0.92, "#bf{CMS} #it{Internal}")
    label_lumi = ROOT.TLatex()
    label_lumi.SetTextSize(0.04)
    label_lumi.DrawLatexNDC(0.5, 0.92, "#sqrt{s} = 13.6 TeV, L_{int} #approx 9 fb^{-1}")

    label_extra = ROOT.TPaveLabel(0.17,0.75,0.9,0.9,"#color[2]{"+extralabel+"}","brNDC")
    label_extra.SetTextFont(43)
    label_extra.SetTextAlign(12)
    label_extra.SetTextSize(18)
    label_extra.SetFillStyle(0)
    label_extra.SetBorderSize(0)
    label_extra.Draw()

    
    legend.Draw()
    if interactive:
        input()
    if saveplot: 
        c.SaveAs(dirname+'/'+plotname+'.png')
        c.SaveAs(dirname+'/'+plotname+'.pdf')
    

def compute_eff(hdens, hnums, addnumtoden):
    effs = []

    for i in range(len(hnums)):
        if addnumtoden:
            hdens[i].Add(hnums[i])
        eff = ROOT.TEfficiency(hnums[i], hdens[i])
        eff.SetFillStyle(0)
        eff.SetLineColor(colors[i])
        eff.SetMarkerStyle(20)
        eff.SetMarkerSize(0.5)
        eff.SetMarkerColor(colors[i])
        effs.append(eff)

    return effs

def setstyle(h, i):
    h.SetFillStyle(0)
    h.SetLineColor(colors[i])
    h.SetMarkerStyle(20)
    h.SetMarkerSize(0.5)
    h.SetMarkerColor(colors[i])
    return h


def compute_profilex(h2d):
    profiles = []
    for i, h in enumerate(h2d):
        profile = h.ProfileX().Clone()
        profile.SetFillStyle(0)
        profile.SetLineColor(colors[i])
        profile.SetMarkerStyle(20)
        profile.SetMarkerSize(0.5)
        profile.SetMarkerColor(colors[i])
        profiles.append(profile)
    
    return profiles

def compute_ResolutionvsX(h2d):
    histos_response, histos_resol = [], []
    for ctr, h in enumerate(h2d):
        h_responsevsX = h.ProjectionX().Clone()
        h_resolvsX = h.ProjectionX().Clone()
        setstyle(h_responsevsX, ctr)
        setstyle(h_resolvsX, ctr)
        for i in range(h.GetNbinsX()+1):
            proj = h.ProjectionY("py",i,i).Clone()
           
            f_gaus = ROOT.TF1("f_gaus","gaus")
            proj.Fit(f_gaus)
            if f_gaus.GetParameter(1) >0 and f_gaus.GetParError(2)/f_gaus.GetParameter(1)<0.03:
                h_responsevsX.SetBinContent(i,f_gaus.GetParameter(1))
                h_responsevsX.SetBinError(i,f_gaus.GetParError(1))
                h_resolvsX.SetBinContent(i,f_gaus.GetParameter(2)/f_gaus.GetParameter(1))
                h_resolvsX.SetBinError(i,f_gaus.GetParError(2)/f_gaus.GetParameter(1))
            else:
                h_responsevsX.SetBinContent(i,0)
                h_responsevsX.SetBinError(i,0)
                h_resolvsX.SetBinContent(i,0)
                h_resolvsX.SetBinError(i,0)

        histos_resol.append(h_resolvsX)
        histos_response.append(h_responsevsX)
    return histos_response, histos_resol
if __name__ == '__main__':
    main()


