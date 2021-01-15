####################################################
####################################################

import os
import sys
import numpy as np

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from Plot_Helper import LoadNtuples, MakeStack, CreateCanvas, DrawOnCanv, SaveCanvPic, MakeLumiLabel, MakeCMSDASLabel, ScaleSignal, MakeRatioPlot, MakeLegend
from Analyzer_Helper import MuPair, DimuCandidates, IsBtagLoose, IsBtagMed, IsBtagTight, CountYield
import Analyzer_Configs as AC
import Plot_Configs     as PC

from Analyzer_ALP import PIso2D, plot2D, plot2D_CONT

import CMS_lumi, tdrstyle

from xgboost import XGBClassifier
import pickle

from array import array
import math
import copy

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
parser.add_argument('-C', '--CR', dest='CR', action='store_true', default=False, help='make control region')
args = parser.parse_args()

mass = 'massIndependent'



gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

mva = False
if mva:
    name = mass + '_mva'
else:
    name = 'ShapeSys'

if args.year == '2016':
    file_out = 'plots_16'
    out_name = "ALP_plot_data16_"+name+".root"
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_massindependent_2016.pkl"
    mvaCut = 0.8675
elif args.year == '2017':
    file_out = 'plots_17'
    out_name = "ALP_plot_data17_"+name+".root"
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_massindependent_2017.pkl"
    mvaCut = 0.8571
elif args.year == '2018':
    file_out = 'plots_18'
    out_name = "ALP_plot_data18_"+name+".root"
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_massindependent_2018.pkl"
    mvaCut = 0.9204
else:
    print "do not include at 2016/2017/2018"
    exit(0)

file_plot = file_out + "/plot_"+name
# load the model from disk
model = pickle.load(open(BDT_filename, 'rb'))

#def DrawShapeSys():
#    return 0

def main():

    if not os.path.exists(file_out):
        os.makedirs(file_out)

    if not os.path.exists(file_plot):
        os.makedirs(file_plot)

    out_file = TFile( file_out + '/' + out_name , "RECREATE")

    analyzer_cfg = AC.Analyzer_Config('inclusive', mass,args.year)
    #analyzer_cfg.samp_names = ["DYJetsToLL", "data"]
    #analyzer_cfg.sample_loc = analyzer_cfg.sample_loc + '/CR'
    print "using samples: "
    print analyzer_cfg.samp_names
    print "getting ntuples from:" + analyzer_cfg.sample_loc
    #analyzer_cfg.Print_Config()
    ntuples = LoadNtuples(analyzer_cfg)


    rootfile = "/publicfs/cms/user/wangzebing/ALP/Plot/plots_17/ALP_plot_data17_massIndependent_mva_CR.root"
    rootfile_sys = "/publicfs/cms/user/wangzebing/ALP/Plot/plots_17/ALP_plot_data17_sysShift.root"
    files = TFile(rootfile)
    files_sys = TFile(rootfile_sys)

    #var_names = ['Z_m', 'H_m', 'ALP_m','pho1Pt', 'pho1eta', 'pho1phi', 'pho1R9', 'pho1IetaIeta', 'pho1IetaIeta55','pho1PIso_noCorr' ,'pho2Pt', 'pho2eta', 'pho2phi', 'pho2R9', 'pho2IetaIeta', 'pho2IetaIeta55','pho2PIso_noCorr','ALP_calculatedPhotonIso', 'var_dR_Za', 'var_dR_g1g2', 'var_dR_g1Z', 'var_PtaOverMh', 'var_Pta', 'var_MhMZ', 'H_pt', 'var_PtaOverMa', 'var_MhMa']
    var_names = ['Z_m']
    Hist_data = {}
    Hist_MC = {}
    Hist_up = {}
    Hist_down = {}
    for var_name in var_names:

        canv = TCanvas(var_name, var_name, 600,750)
        canv.SetBottomMargin(0.01)
        canv.cd()

        Hist_data[var_name] = files.raw_plots.Get(var_name + '_data')
        Hist_MC[var_name] = files.raw_plots.Get(var_name + '_DYJetsToLL')

        h_data = copy.deepcopy(Hist_data[var_name])
        Ntot_Data = h_data.Integral()
        h_data.Sumw2()
        print "Total data: " + str(Ntot_Data)

        h_MC = copy.deepcopy(Hist_MC[var_name])
        Ntot_MC = h_MC.Integral()
        print "Total MC: " + str(Ntot_MC)
        scale_MC = Ntot_Data*1.0/Ntot_MC
        h_MC.Sumw2()
        h_MC.Scale(scale_MC)
        print "After scale: " + str(h_MC.Integral())

        sys = 0.5
        Hist_up[var_name] = files_sys.plots_up.Get('UP_' + var_name + '_DYJetsToLL')
        Hist_down[var_name] = files_sys.plots_down.Get('Down_' + var_name + '_DYJetsToLL')

        h_MCP = copy.deepcopy(Hist_up[var_name])
        h_MCM = copy.deepcopy(Hist_down[var_name])

        h_MCP.Scale(scale_MC)
        h_MCM.Scale(scale_MC)

        print "Scale up: " + str(h_MCP.Integral())
        print "Scale down: " + str(h_MCM.Integral())

        BinTotal = int(h_data.GetXaxis().GetNbins())
        BinXHig = h_data.GetXaxis().GetXmin()
        BinXLow = h_data.GetXaxis().GetXmax()

        #print BinTotal
        #print BinXHig
        #print BinXLow

        h_MCUP = TH1F("h_MCUP","",BinTotal,BinXLow,BinXHig)
        h_MCDN = TH1F("h_MCDN","",BinTotal,BinXLow,BinXHig)
        h_MCUPNorm = TH1F("h_MCUPNorm","",BinTotal,BinXLow,BinXHig)
        h_MCDNNorm = TH1F("h_MCDNNorm","",BinTotal,BinXLow,BinXHig)

        YMean = array('f',BinTotal*[0.])
        YMeanNorm = array('f',BinTotal*[0.])
        XMean = array('f',BinTotal*[0.])
        X_ErrH = array('f',BinTotal*[0.])
        X_ErrL = array('f',BinTotal*[0.])
        MCSys_ErrH = array('f',BinTotal*[0.])
        MCSys_ErrL = array('f',BinTotal*[0.])
        MCSys_RelErrH = array('f',BinTotal*[0.])
        MCSys_RelErrL = array('f',BinTotal*[0.])
        WidthBin=h_data.GetBinWidth(1)
        Chi2 = 0.

        for ibin in range(BinTotal):
            Nd = h_data.GetBinContent(ibin+1)
            Nm = h_MC.GetBinContent(ibin+1)
            NmErr = h_MC.GetBinError(ibin+1)

            Chi2 += math.fabs(NmErr)>1e-9 and (Nm-Nd)*(Nm-Nd)*1.0/(NmErr*NmErr) or 0.0

            YMean[ibin] = Nm
            YMeanNorm[ibin] = 1.0
            XMean[ibin] = BinXLow+ (ibin+0.5)*WidthBin
            X_ErrH[ibin] = 0.5*WidthBin
            X_ErrL[ibin] = 0.5*WidthBin
            dNmp = h_MCP.GetBinContent(ibin+1) - Nm
            dNmm = h_MCM.GetBinContent(ibin+1) - Nm
            dplus = dNmp>dNmm and dNmp or dNmm
            if(dplus<0): dplus *= -1
            dminus = dNmp<dNmm and dNmp or dNmm
            if(dminus>0):
                dminus *= -1
                dminus *= -1

            MCSys_ErrH[ibin] = dplus
            MCSys_ErrL[ibin] = dminus

            h_MCUP.SetBinContent(ibin+1, Nm+dplus)
            h_MCDN.SetBinContent(ibin+1, Nm-dminus)
            MCSys_RelErrH[ibin] = Nm>0 and dplus/Nm or 0.
            MCSys_RelErrL[ibin] = Nm>0 and dminus/Nm or 0.
            h_MCUPNorm.SetBinContent(ibin+1, 1.+MCSys_RelErrH[ibin])
            h_MCDNNorm.SetBinContent(ibin+1, 1.-MCSys_RelErrL[ibin])

        gr_MCSys = TGraphAsymmErrors(BinTotal, XMean, YMean, X_ErrL, X_ErrH, MCSys_ErrL, MCSys_ErrH)
        gr_MCSysNorm = TGraphAsymmErrors(BinTotal, XMean, YMeanNorm, X_ErrL, X_ErrH, MCSys_RelErrL, MCSys_RelErrH)

        print gr_MCSys

        gr_MCSys.SetFillColor(kRed-10)
        gr_MCSys.SetFillStyle(3001)
        gr_MCSysNorm.SetFillColor(kRed-10)
        gr_MCSysNorm.SetFillStyle(3001)

        h_MCUP.SetLineColor(2)
        h_MCUP.SetLineStyle(1)
        h_MCUP.SetLineWidth(2)
        h_MCDN.SetLineColor(2)
        h_MCDN.SetLineStyle(1)
        h_MCDN.SetLineWidth(2)
        h_MCUPNorm.SetLineColor(2)
        h_MCUPNorm.SetLineStyle(1)
        h_MCUPNorm.SetLineWidth(2)
        h_MCDNNorm.SetLineColor(2)
        h_MCDNNorm.SetLineStyle(1)
        h_MCDNNorm.SetLineWidth(2)


        maxY=max(h_data.GetMaximum(),h_MC.GetMaximum())
        maxY = max(maxY, h_MCUP.GetMaximum())
        maxY = max(maxY, h_MCDN.GetMaximum())
        minY=min(h_data.GetMinimum(),h_MC.GetMinimum())
        h_data.GetYaxis().SetRangeUser(0.95*minY, 1.05*maxY)
        minY = min(minY, h_MCUP.GetMinimum())
        minY = min(minY, h_MCDN.GetMinimum())

        h_data.SetLineColor(1)
        h_data.SetFillColor(0)
        h_data.SetLineStyle(1)
        h_data.SetLineWidth(2)

        PreTitleY = "Events / {0} ".format(WidthBin)
        TitleY =  PreTitleY + "GeV"
        h_data.GetYaxis().SetTitle(TitleY)

        h_data.SetTitleSize(0.05,"X")
        h_data.SetTitleSize(0.05,"Y")
        h_data.SetTitleOffset(1.1, "Y")
        h_data.SetMarkerColor(kBlack)
        h_data.SetMarkerSize(1.0)
        h_data.SetMarkerStyle(20)

        h_MC.SetFillColor(8)
        h_MC.SetMarkerStyle(0)
        h_MC.SetLineColor(8)
        h_MC.SetLineStyle(1)
        h_MC.SetLineWidth(2)

        legend = TLegend(0.5,0.8,0.95,0.92)
        legend.AddEntry(h_data,"Data","pe")
        legend.AddEntry(h_MC,"Z#rightarrow#mu#mu#gamma","f")
        legend.AddEntry(gr_MCSys,"MC syst.","f")

        #prepare 2 pads

        h_data.GetXaxis().SetLabelColor(0)
        h_data.SetNdivisions(510 ,"X")

        #upper_pad = TPad("upperpad_"+var_name, "upperpad_"+var_name, 0,0.2, 1,1)
        #upper_pad.SetBottomMargin(0.05)
        #upper_pad.Draw()
        #upper_pad.cd()

        h_data.Draw("PE1")
        h_MC.Draw("hist,same")
        gr_MCSys.Draw("same2")
        #h_MCUP.Draw("same")
        #h_MCDN.Draw("same")
        legend.Draw("same")
        h_data.Draw("samePE1")
        h_data.Draw("Axissame")

        #canv.SetBottomMargin(0.01)
        #canv.cd()

        #upper_pad = TPad("upperpad_"+var_name, "upperpad_"+var_name, 0,0.2, 1,1)
        #upper_pad.SetBottomMargin(0.05)
        #upper_pad.Draw()
        #upper_pad.cd()

        #canv.Draw()
        PrintInfor1="#bf{CMS} #it{} #it{Preliminary}"
        PrintInfor2="35.9 fb^{-1} (13TeV)"
        tex = TLatex(0.16,0.94, PrintInfor1)
        tex.SetNDC()
        tex.SetTextFont(42)
        tex.SetTextSize(0.045)
        tex.SetLineWidth(2)
        tex.Draw()

        tex = TLatex(0.70,0.94, PrintInfor2)
        tex.SetNDC()
        tex.SetTextFont(42)
        tex.SetTextSize(0.045)
        tex.SetLineWidth(2)
        tex.Draw()


        Line1 = TLine(h_data.GetBinLowEdge(1),1,h_data.GetBinLowEdge(h_data.GetNbinsX())+ h_data.GetBinWidth(h_data.GetNbinsX()),1)
        Line1.SetLineColor(1)
        Line1.SetLineWidth(2)
        Line1.SetLineStyle(4)

        canv.SaveAs('test.png')







    '''
    #var_names = ['Z_m', 'H_m', 'ALP_m', 'pho1IetaIeta', 'pho1IetaIeta55', 'pho1PIso', 'var_dR_Za', 'var_dR_g1g2', 'var_dR_g1Z', 'var_PtaOverMa', 'var_PtaOverMh', 'var_Pta', 'var_MhMa', 'var_MhMZ', 'ALP_calculatedPhotonIso']
    var_names = ['Z_m', 'H_m', 'ALP_m','pho1Pt', 'pho1eta', 'pho1phi', 'pho1R9', 'pho1IetaIeta', 'pho1IetaIeta55','pho1PIso_noCorr' ,'pho2Pt', 'pho2eta', 'pho2phi', 'pho2R9', 'pho2IetaIeta', 'pho2IetaIeta55','pho2PIso_noCorr','ALP_calculatedPhotonIso', 'var_dR_Za', 'var_dR_g1g2', 'var_dR_g1Z', 'var_PtaOverMh', 'var_Pta', 'var_MhMZ', 'H_pt', 'var_PtaOverMa', 'var_MhMa']
    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)

    ### declare histograms

    histos = {}

    for var_name in var_names:
        histos[var_name] = {}
    for sample in analyzer_cfg.samp_names:
        histos['pho1Pt'][sample]    = TH1F('pho1Pt'    + '_' + sample, 'pho1Pt'    + '_' + sample, 50,  5., 50.)
        histos['pho1eta'][sample]    = TH1F('pho1eta'    + '_' + sample, 'pho1eta'    + '_' + sample, 100,  -3., 3.)
        histos['pho1phi'][sample]    = TH1F('pho1phi'    + '_' + sample, 'pho1phi'    + '_' + sample, 100,  -4., 4.)
        histos['pho1R9'][sample]    = TH1F('pho1R9'    + '_' + sample, 'pho1R9'    + '_' + sample, 100,  0., 1.)
        histos['pho1IetaIeta'][sample]    = TH1F('pho1IetaIeta'    + '_' + sample, 'pho1IetaIeta'    + '_' + sample, 100,  0., 0.07)
        histos['pho1IetaIeta55'][sample]    = TH1F('pho1IetaIeta55'    + '_' + sample, 'pho1IetaIeta55'    + '_' + sample, 100,  0., 0.07)
        histos['pho1PIso_noCorr'][sample]    = TH1F('pho1PIso_noCorr'    + '_' + sample, 'pho1PIso_noCorr'    + '_' + sample, 100, 0., 40.)
        histos['pho2Pt'][sample]    = TH1F('pho2Pt'    + '_' + sample, 'pho2Pt'    + '_' + sample, 50,  5., 30.)
        histos['pho2eta'][sample]    = TH1F('pho2eta'    + '_' + sample, 'pho2eta'    + '_' + sample, 100,  -3., 3.)
        histos['pho2phi'][sample]    = TH1F('pho2phi'    + '_' + sample, 'pho2phi'    + '_' + sample, 100,  -4., 4.)
        histos['pho2R9'][sample]    = TH1F('pho2R9'    + '_' + sample, 'pho2R9'    + '_' + sample, 100,  0., 1.)
        histos['pho2IetaIeta'][sample]    = TH1F('pho2IetaIeta'    + '_' + sample, 'pho2IetaIeta'    + '_' + sample, 100,  0., 0.07)
        histos['pho2IetaIeta55'][sample]    = TH1F('pho2IetaIeta55'    + '_' + sample, 'pho2IetaIeta55'    + '_' + sample, 100,  0., 0.07)
        histos['pho2PIso_noCorr'][sample]    = TH1F('pho2PIso_noCorr'    + '_' + sample, 'pho2PIso_noCorr'    + '_' + sample, 100, 0., 40.)
        histos['Z_m'][sample]    = TH1F('Z_m'    + '_' + sample, 'Z_m'    + '_' + sample, 100,  50., 130.)
        histos['H_m'][sample]    = TH1F('H_m'    + '_' + sample, 'H_m'    + '_' + sample, 20,  118., 140.)
        histos['H_pt'][sample]    = TH1F('H_pt'    + '_' + sample, 'H_pt'    + '_' + sample, 100,  0., 200.)
        histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 50, 0., 40.)
        histos['var_dR_g1g2'][sample] = TH1F('var_dR_g1g2' + '_' + sample, 'var_dR_g1g2' + '_' + sample, 50, 0., 5)
        histos['var_PtaOverMa'][sample] = TH1F('var_PtaOverMa' + '_' + sample, 'var_PtaOverMa' + '_' + sample, 100, 0., 100.)
        histos['var_dR_Za'][sample] = TH1F('var_dR_Za' + '_' + sample, 'var_dR_Za' + '_' + sample, 50, 0., 7.)
        histos['var_dR_g1Z'][sample] = TH1F('var_dR_g1Z' + '_' + sample, 'var_dR_g1Z' + '_' + sample, 50, 0., 8.)
        histos['var_PtaOverMh'][sample] = TH1F('var_PtaOverMh' + '_' + sample, 'var_PtaOverMh' + '_' + sample, 100, 0., 1.)
        histos['var_Pta'][sample] = TH1F('var_Pta' + '_' + sample, 'var_Pta' + '_' + sample, 50, 0., 60.)
        histos['var_MhMa'][sample] = TH1F('var_MhMa' + '_' + sample, 'var_MhMa' + '_' + sample, 50, 100., 200.)
        histos['var_MhMZ'][sample] = TH1F('var_MhMZ' + '_' + sample, 'var_MhMZ' + '_' + sample, 50, 120., 250.)
        histos['ALP_calculatedPhotonIso'][sample] = TH1F('ALP_calculatedPhotonIso' + '_' + sample, 'ALP_calculatedPhotonIso' + '_' + sample, 100, 0., 5.)



    ### loop over samples and events
    for sample in analyzer_cfg.samp_names:
        ntup = ntuples[sample] # just a short name
        print '\n\nOn sample: %s' %sample
        print 'total events: %d' %ntup.GetEntries()

        for iEvt in range( ntup.GetEntries() ):
            ntup.GetEvent(iEvt)
            if (iEvt == 5000): break


            if (iEvt % 100000 == 1):
                print "looking at event %d" %iEvt


            weight = ntup.factor

            if (ntup.H_m > -90):

                #if ntup.dR_pho < 0.02: continue
                #if not ntup.passEleVeto: continue
                if not ntup.passChaHadIso: continue
                if not ntup.passNeuHadIso: continue
                if not ntup.passdR_gl: continue
                if not ntup.passHOverE: continue
                if ntup.H_m>130. or ntup.H_m<118.: continue
                #if ntup.H_m<130. and ntup.H_m>118.: continue

                #if abs(ntup.l1_id)!=13: continue

                #print 'pho1: pt:' + str(ntup.pho1Pt) + ', eta: ' + str(ntup.pho1eta) + ', SFs: ' + str(pho1_SFs)
                if mva:
                    #MVA_list = [ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_PtaOverMh, ntup.var_MhMZ, ntup.var_dR_g1g2]
                    MVA_list = [ntup.pho1Pt, ntup.pho1eta, ntup.pho1phi, ntup.pho1R9, ntup.pho1IetaIeta55 ,ntup.pho2Pt, ntup.pho2eta, ntup.pho2phi, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, ntup.var_Pta, ntup.var_MhMZ, ntup.H_pt ]
                    MVA_value = model.predict_proba(MVA_list)[:, 1]
                    if MVA_value < mvaCut:continue

                #histos['Z_befor50'][sample].Fill( ntup.Z_befor50, weight )
                #if ntup.pho1Pt<15: continue
                histos['pho1Pt'][sample].Fill( ntup.pho1Pt, weight )
                histos['pho1eta'][sample].Fill( ntup.pho1eta, weight )
                histos['pho1phi'][sample].Fill( ntup.pho1phi, weight )
                histos['pho1R9'][sample].Fill( ntup.pho1R9, weight )
                histos['pho1IetaIeta'][sample].Fill( ntup.pho1IetaIeta, weight )
                histos['pho1IetaIeta55'][sample].Fill( ntup.pho1IetaIeta55, weight )
                histos['pho1PIso_noCorr'][sample].Fill( ntup.pho1PIso_noCorr, weight )
                histos['pho2Pt'][sample].Fill( ntup.pho2Pt, weight )
                histos['pho2eta'][sample].Fill( ntup.pho2eta, weight )
                histos['pho2phi'][sample].Fill( ntup.pho2phi, weight )
                histos['pho2R9'][sample].Fill( ntup.pho2R9, weight )
                histos['pho2IetaIeta'][sample].Fill( ntup.pho2IetaIeta, weight )
                histos['pho2IetaIeta55'][sample].Fill( ntup.pho2IetaIeta55, weight )
                histos['pho2PIso_noCorr'][sample].Fill( ntup.pho2PIso_noCorr, weight )


                histos['H_m'][sample].Fill( ntup.H_m, weight )
                histos['H_pt'][sample].Fill( ntup.H_pt, weight )
                histos['ALP_m'][sample].Fill( ntup.ALP_m, weight )
                histos['Z_m'][sample].Fill( ntup.Z_m, weight )

                histos['var_dR_Za'][sample].Fill( ntup.var_dR_Za, weight )
                histos['var_dR_g1g2'][sample].Fill( ntup.var_dR_g1g2, weight )
                histos['var_dR_g1Z'][sample].Fill( ntup.var_dR_g1Z, weight )
                histos['var_PtaOverMa'][sample].Fill( ntup.var_PtaOverMa, weight )
                histos['var_PtaOverMh'][sample].Fill( ntup.var_PtaOverMh, weight )
                histos['var_Pta'][sample].Fill( ntup.var_Pta, weight )
                histos['var_MhMa'][sample].Fill( ntup.var_MhMa, weight )
                histos['var_MhMZ'][sample].Fill( ntup.var_MhMZ, weight )
                histos['ALP_calculatedPhotonIso'][sample].Fill( ntup.ALP_calculatedPhotonIso, weight )



        ## End of for iEvt in range( ntup.GetEntries() )
    ## End of for sample in analyzer_cfg.samp_names

    ### save raw histograms
    raw_dir = out_file.mkdir('raw_plots')
    raw_dir.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.samp_names:
            plot_cfg.SetHistStyles(histos[var_name][sample], sample)
            histos[var_name][sample].Write()

    ### save stack plots and make ratio plots
    out_file.cd()
    lumi_label = MakeLumiLabel(plot_cfg.lumi)
    cms_label  = MakeCMSDASLabel()

    scaled_sig = {}
    for var_name in var_names:
        stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
        #scaled_sig = 0
        #ratio_plot = 0
        for sample in analyzer_cfg.sig_names:
            scaled_sig[sample] = ScaleSignal(plot_cfg, stacks[sample], histos[var_name][sample], var_name)
        ratio_plot = MakeRatioPlot(histos[var_name]['data'], stacks['all'].GetStack().Last(), var_name)
        legend = MakeLegend(plot_cfg, histos[var_name], scaled_sig)

        canv = CreateCanvas(var_name)
        DrawShapeSys(canv, var_name, stacks, histos[var_name])

        canv.Write()
        SaveCanvPic(canv, file_plot, var_name)
    '''
    print '\n\n'
    #CountYield(analyzer_cfg, histos['ALP_m'])
    out_file.Close()

    print 'Done'


main()
