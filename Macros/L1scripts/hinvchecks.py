from datetime import datetime
import ROOT
import json
import os
import sys
import argparse



leadingpt = [60, 70, 80, 90, 100, 110, 120, 130]
trailingpt = [30, 35, 40, 50, 60, 70]
deta = [2., 2.5, 3., 3.5, 4.0, 4.5, 5]
mjj = [400., 500., 550., 600., 650., 700.]
dphi = [0.8, 1.6, 2.4, 3.2]

#In case you want to load an helper for C++ functions
ROOT.gInterpreter.Declare('#include "Helper.h"')
ROOT.gInterpreter.Declare('#include "Helper_InvariantMass.h"')
#Importing stuff from other python files
from helper import * 

def main():
    ###Arguments 
    parser = argparse.ArgumentParser(
        description='''
        H invisible checks
        ''',
        usage='use "%(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--maxEvents", dest="maxEvents", help="Maximum number of events to analyze. Default=-1 i.e. run on all events.", type=int, default=-1)
    parser.add_argument("-i", "--input", dest="inputfile", help="Input file", type=str, default='')
    parser.add_argument("-o", "--output", dest="outputfile", help="Output file", type=str, default='output.root')
    parser.add_argument("-s", "--study", dest="study", help="Type of study (rate or efficiency)", type=str, default='efficiency')
    parser.add_argument("-f", "--format", dest="dataformat", help="Input format (NANO or CustomNtuple)", type=str, default='CustomNtuple')
    args = parser.parse_args() 

    
    ###Define the RDataFrame from the input tree
    inputfile = args.inputfile
    if args.dataformat == 'NANO' and args.inputfile == '':
        inputfile = '/user/lathomas/PortingCodeToGitIIHECMS/CMSSW_12_4_8/src/GenericTreeProducerMINIAOD/Ntuplizer/python/outputvbfhinv.root'

    df = None 
    if args.dataformat == 'NANO':
        df = ROOT.RDataFrame('Events', inputfile)
    #    elif args.study == 'rate':
    else:
        df = ROOT.RDataFrame('ntuplizer/tree', inputfile)
    nEvents = df.Count().GetValue()
    print('There are {} events'.format(nEvents))
    
    #Max events to run on 
    maxEvents = min(nEvents, args.maxEvents) if args.maxEvents >=0 else nEvents
    df = df.Range(0, maxEvents)
    #Next line to monitor event loop progress
    df = df.Filter('if(tdfentry_ %100000 == 0) {cout << "Event is  " << tdfentry_ << endl;} return true;')
    

    out = ROOT.TFile(args.outputfile, "recreate")
    ####The sequence of filters/column definition starts here

    
    if args.dataformat == 'NANO': 
        df = df.Define("mjj","(HighestMjj(GenJet_pt, GenJet_eta, GenJet_phi))[0]")
        df = df.Define("leading_pt","(HighestMjj(GenJet_pt, GenJet_eta, GenJet_phi))[1]")
        df = df.Define("subleading_pt","(HighestMjj(GenJet_pt, GenJet_eta, GenJet_phi))[2]")
        df = df.Define("deta","(HighestMjj(GenJet_pt, GenJet_eta, GenJet_phi))[3]")
        df = df.Define("dphi","(HighestMjj(GenJet_pt, GenJet_eta, GenJet_phi))[4]")
        
    else:
        df = df.Define("mjj","(HighestMjj_L1(_L1jet_pt, _L1jet_eta, _L1jet_phi, _L1jet_bx))[0]")
        df = df.Define("leading_pt","(HighestMjj_L1(_L1jet_pt, _L1jet_eta, _L1jet_phi, _L1jet_bx))[1]")
        df = df.Define("subleading_pt","(HighestMjj_L1(_L1jet_pt, _L1jet_eta, _L1jet_phi, _L1jet_bx))[2]")
        df = df.Define("deta","(HighestMjj_L1(_L1jet_pt, _L1jet_eta, _L1jet_phi, _L1jet_bx))[3]")
        df = df.Define("dphi","(HighestMjj_L1(_L1jet_pt, _L1jet_eta, _L1jet_phi, _L1jet_bx))[4]")
        

        for h in dphi: 
            
            for i in leadingpt:
                for j in trailingpt:
                    if j > i:
                        continue
                    for k in deta:
                        df = df.Define("passL1_DoubleJetPt{}_{}_dEta{}_dPhi{}".format(i,j,k,h).replace(".","p"),"L1SeedDoubleJetEtaMin(_L1jet_pt, _L1jet_eta, _L1jet_phi, _L1jet_bx, {}, {}, {}, {})".format(i, j, k, h))
                    for k in mjj:
                        df = df.Define("passL1_DoubleJetPt{}_{}_Mjj{}_dPhi{}".format(i,j,k,h).replace(".","p"),"L1SeedDoubleJetMassMin(_L1jet_pt, _L1jet_eta, _L1jet_phi, _L1jet_bx, {}, {}, {}, {})".format(i, j, k, h))
                    
    histos = {}
    if args.dataformat == 'NANO':
        histos['met'] = df.Histo1D(ROOT.RDF.TH1DModel('hmet', '', 5000, 0, 5000), 'GenMET_pt')

    #df = df.Filter('passL1_DoubleJetPt100_30_Mjj500p0')
    histos['mjj'] = df.Histo1D(ROOT.RDF.TH1DModel('hmjj', '', 5000, 0, 5000), 'mjj')
    histos['leadingpt'] = df.Histo1D(ROOT.RDF.TH1DModel('hleadingpt', '', 1000, 0, 1000), 'leading_pt')
    histos['deta'] = df.Histo1D(ROOT.RDF.TH1DModel('hdeta', '', 1000, 0, 10), 'deta')
    histos['dphi'] = df.Histo1D(ROOT.RDF.TH1DModel('hdphi', '', 100, 0, 3.1416), 'dphi')
    histos['subleadingpt'] = df.Histo1D(ROOT.RDF.TH1DModel('hsubleadingpt', '', 1000, 0, 1000), 'subleading_pt')

    if args.dataformat != 'NANO':
        for h in dphi:
            for i in leadingpt:
                for j in trailingpt:
                    print(i, j)
                    print('i> j', i > j)
                    if j > i :
                        continue
                    print('passing i, j', i, j)
                    for k in deta:
                        print('k', k)
                        print('passL1_DoubleJetPt{}_{}_dEta{}_dPhi{}'.format(i,j,k,h).replace(".","p"))
                        histos['passL1_DoubleJetPt{}_{}_dEta{}_dPhi{}'.format(i,j,k,h).replace(".","p")] = df.Histo1D(ROOT.RDF.TH1DModel('passL1_DoubleJetPt{}_{}_dEta{}_dPhi{}'.format(i,j,k,h).replace(".","p"), '', 2, 0, 2), 'passL1_DoubleJetPt{}_{}_dEta{}_dPhi{}'.format(i,j,k,h).replace(".","p"))
                    for k in mjj:
                        histos['passL1_DoubleJetPt{}_{}_Mkk{}_dPhi{}'.format(i,j,k,h).replace(".","p")] = df.Histo1D(ROOT.RDF.TH1DModel('passL1_DoubleJetPt{}_{}_Mjj{}_dPhi{}'.format(i,j,k,h).replace(".","p"), '', 2, 0, 2), 'passL1_DoubleJetPt{}_{}_Mjj{}_dPhi{}'.format(i,j,k,h).replace(".","p"))
    
    for i in histos:
        histos[i].Write()
    

if __name__ == '__main__':
    main()
