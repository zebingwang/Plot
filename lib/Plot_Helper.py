import os
import Plot_Configs as PC
import Analyzer_Configs as AC
from ROOT import *
import ROOT
import numpy as np

import CMS_lumi
from array import array

import math
#####################################################################

def LoadNtuples(ana_cfg):
    ntuples = {}
    for sample in ana_cfg.samp_names:
        ntuples[sample] = TChain("passedEvents","chain_" + sample)
        ntuples[sample]. Add(ana_cfg.sample_loc + '/ALP_%s.root' %sample)
    return ntuples


def MakeStack(histos, ana_cfg, var_name):
    stacks = {}
    stacks['all']  = THStack("h_stack_all_"+var_name, "all_"+var_name)
    stacks['sig']  = THStack("h_stack_sig_"+var_name, "sig_"+var_name)
    stacks['bkg']  = THStack("h_stack_bkg_"+var_name, "bkg_"+var_name)

    for sample in ana_cfg.samp_names:
        stacks[sample] = THStack("h_stack_"+sample+"_"+var_name, sample+"_"+var_name)

    for sample in ana_cfg.bkg_names:
        stacks['bkg'].Add(histos[sample])
        stacks['all'].Add(histos[sample])

    for sample in ana_cfg.sig_names:
        stacks['sig'].Add(histos[sample])
        stacks[sample].Add(histos[sample])

    return stacks

def CreateCanvas(canv_name):
    canv = TCanvas(canv_name, canv_name, 600,750)
    return canv

def MakeLumiLabel(lumi):
    tex = TLatex()
    tex.SetTextSize(0.035)
    tex.SetTextAlign(31)
    tex.DrawLatexNDC(0.90, 0.91, '%s fb^{-1} (13 TeV)' %lumi)
    return tex

def MakeCMSDASLabel():
    #tex = TLatex()
    #tex.SetTextSize(0.03)
    #tex.DrawLatexNDC(0.15, 0.85, '#scale[1.5]{CMSDAS} H To Z + ALP')
    #return tex

    onTop=False
    text='#bf{CMS} #scale[0.75]{#it{Simulation Preliminary}  H#rightarrow#gamma#gamma}'
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.05)
    latex.DrawLatex(0.1, 0.85 if not onTop else 0.93, text)
    return latex

def ScaleSignal(plt_cfg, stack_sig, hist_data, var_name):
    sig_hist = hist_data
    sig_hist.SetLineWidth(3)
    sig_hist.SetFillStyle(0)

    sig_hist.GetXaxis().SetTitle(var_name)
    sig_hist.GetXaxis().SetTitleSize(0.5)
    sig_hist.GetYaxis().SetTitle('Events / %.2f' %sig_hist.GetBinWidth(1))
    return sig_hist

def MakeRatioPlot(h_data, h_MC, var_name):
    ratio_plot = TGraphAsymmErrors()
    ratio_plot.Divide(h_data, h_MC, "pois")
    ratio_plot.SetName("ratiograph_" + var_name)
    ratio_plot.SetMinimum(0.4)
    ratio_plot.SetMaximum(1.6)
    ratio_plot.SetMarkerStyle(20)

    ratio_plot.GetXaxis().SetLimits( h_data.GetXaxis().GetXmin(), h_data.GetXaxis().GetXmax() )
    ratio_plot.GetXaxis().SetLabelSize(0.15)
    ratio_plot.GetXaxis().SetTitle(var_name)
    ratio_plot.GetXaxis().SetTitleSize(0.15)
    ratio_plot.GetXaxis().SetTitleOffset(0.8)
    ratio_plot.GetXaxis().SetTickLength(0.15)

    ratio_plot.GetYaxis().SetNdivisions(505)
    ratio_plot.GetYaxis().SetLabelSize(0.1)
    ratio_plot.GetYaxis().SetTitle("data/MC")
    ratio_plot.GetYaxis().SetTitleSize(0.15)
    ratio_plot.GetYaxis().SetTitleOffset(0.2)

    return ratio_plot

def MakeLegend(plt_cfg, histos, scaled_signal):
    legend = TLegend(0.5,0.8,0.95,0.92)
    legend.SetNColumns(3)
    legend.AddEntry(histos["data"], "data")
    for sample in plt_cfg.ana_cfg.sig_names:
        legend.AddEntry(scaled_signal[sample], sample )

    for sample in plt_cfg.ana_cfg.bkg_names:
        legend.AddEntry(histos[sample], sample )

    return legend

def Get_StatUnc(hist):
    TH1.Sumw2

    BinTotal = hist.GetNbinsX()
    WidthBin = hist.GetBinWidth(1)

    XMean=[]
    YMean=[]
    YMeanNorm=[]
    X_ErrH=[]
    X_ErrL=[]
    Y_ErrH=[]
    Y_ErrL=[]
    YNorm_ErrH=[]
    YNorm_ErrL=[]

    xaxis = hist.GetXaxis()

    graph = TGraphAsymmErrors()

    for iBin in range(1, BinTotal+1):
        N = hist.GetBinContent(iBin)
        NErr = hist.GetBinError(iBin)
        XMean.append(xaxis.GetBinCenter(iBin))
        YMean.append(N)
        YMeanNorm.append(1.0)
        X_ErrH.append(0.5*WidthBin)
        X_ErrL.append(0.5*WidthBin)
        Y_ErrH.append(NErr)
        Y_ErrL.append(NErr)

        if N>0:
            errNorm = NErr/N
        else:
            errNorm = 0.
        YNorm_ErrH.append(errNorm)
        YNorm_ErrL.append(errNorm)


    #print YMean
    #print Y_ErrH

    #print 'stat:'
    #print 'high:'
    #print Y_ErrH
    #print Y_ErrL

    graph = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMean), np.array(X_ErrL), np.array(X_ErrH), np.array(Y_ErrL), np.array(Y_ErrH))
    graph_norm = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMeanNorm), np.array(X_ErrL), np.array(X_ErrH), np.array(YNorm_ErrH), np.array(YNorm_ErrL))

    return graph, graph_norm


def Get_SysUnc(hist, hist_up, hist_dn):
    TH1.Sumw2

    BinTotal = hist.GetNbinsX()
    WidthBin = hist.GetBinWidth(1)

    XMean=[]
    YMean=[]
    YMeanNorm=[]
    X_ErrH=[]
    X_ErrL=[]
    YSys_ErrH=[]
    YSys_ErrL=[]
    YSysNorm_ErrH=[]
    YSysNorm_ErrL=[]

    xaxis = hist.GetXaxis()

    graph_sys = TGraphAsymmErrors()

    for iBin in range(1, BinTotal+1):
        N = hist.GetBinContent(iBin)
        NErr = hist.GetBinError(iBin)

        dNp = hist_up.GetBinContent(iBin) - N
        dNm = hist_dn.GetBinContent(iBin) - N

        dplus = dNp if (dNp>dNm) else dNm
        if dplus < 0: dplus = -1.0*dplus
        dminus = dNp if (dNp<dNm) else dNm
        if dminus > 0: dminus = -1.0*dminus
        dminus = -1.0*dminus

        XMean.append(xaxis.GetBinCenter(iBin))
        YMean.append(N)
        YMeanNorm.append(1.0)
        X_ErrH.append(0.5*WidthBin)
        X_ErrL.append(0.5*WidthBin)
        YSys_ErrH.append(dplus)
        YSys_ErrL.append(dminus)

        if N>0:
            errNormH = dplus/N
            errNormL = dminus/N
        else:
            errNormH = 0.
            errNormL = 0.
        YSysNorm_ErrH.append(errNormH)
        YSysNorm_ErrL.append(errNormL)

    graph_sys = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMean), np.array(X_ErrL), np.array(X_ErrH), np.array(YSys_ErrH), np.array(YSys_ErrL))
    #graph_sys_norm = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMeanNorm), np.array(X_ErrL), np.array(X_ErrH), np.array(YSysNorm_ErrH), np.array(YSysNorm_ErrL))

    #print "high:"
    #print YSys_ErrH
    #print "low:"
    #print YSys_ErrL

    return YSys_ErrH, YSys_ErrL, YSysNorm_ErrH, YSysNorm_ErrL

def Total_Unc(hist_norm, hist_sys, sys_names, analyzer_cfg):
    TH1.Sumw2
    hist_norm = hist_norm.GetStack().Last()

    stacks_sys = {}
    for sys in sys_names:
        stacks_sys['bkgSys_'+sys]  = THStack("h_stack_bkgSys_"+sys, "bkgSys_"+sys)

    
    for sys in sys_names:
        for sample in analyzer_cfg.bkg_names:
            stacks_sys['bkgSys_'+sys].Add(hist_sys[sample][sys])

    #print stacks_sys

    YSys_ErrH = {}
    YSys_ErrL = {}
    YSysNorm_ErrH = {}
    YSysNorm_ErrL = {}
    for i in range(len(sys_names)/2):
        #print stacks_sys['bkgSys_'+sys_names[2*i]].GetStack().Last()
        YSys_ErrH[i], YSys_ErrL[i], YSysNorm_ErrH[i], YSysNorm_ErrL[i] = Get_SysUnc(hist_norm, stacks_sys['bkgSys_'+sys_names[2*i]].GetStack().Last(), stacks_sys['bkgSys_'+sys_names[2*i+1]].GetStack().Last())


    BinTotal = hist_norm.GetNbinsX()
    WidthBin = hist_norm.GetBinWidth(1)

    XMean=[]
    YMean=[]
    YMeanNorm=[]
    X_ErrH=[]
    X_ErrL=[]
    Y_ErrH=[]
    Y_ErrL=[]
    YNorm_ErrH=[]
    YNorm_ErrL=[]

    xaxis_total = hist_norm.GetXaxis()

    graph_Total = TGraphAsymmErrors()
    graph_norm_Total = TGraphAsymmErrors()

    for iBin in range(1, BinTotal+1):
        N = hist_norm.GetBinContent(iBin)
        NErr = hist_norm.GetBinError(iBin)
        XMean.append(xaxis_total.GetBinCenter(iBin))
        YMean.append(N)
        YMeanNorm.append(1.0)
        X_ErrH.append(0.5*WidthBin)
        X_ErrL.append(0.5*WidthBin)
        
        NErrH = 0.0
        NErrL = 0.0
        for i in range(len(sys_names)/2):
            NErrH = NErrH + YSys_ErrH[i][iBin-1] * YSys_ErrH[i][iBin-1] 
            NErrL = NErrL + YSys_ErrL[i][iBin-1] * YSys_ErrL[i][iBin-1]

        NErrH = np.sqrt(NErrH + NErr*NErr)
        NErrL = np.sqrt(NErrL + NErr*NErr)

        Y_ErrH.append(NErrH)
        Y_ErrL.append(NErrL)

        if N>0:
            errNormH = NErrH/N
            errNormL = NErrL/N
        else:
            errNormH = 0.
            errNormL = 0.
        YNorm_ErrH.append(errNormH)
        YNorm_ErrL.append(errNormL)

    graph_Total = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMean), np.array(X_ErrL), np.array(X_ErrH), np.array(Y_ErrL), np.array(Y_ErrH))
    graph_norm_Total = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMeanNorm), np.array(X_ErrL), np.array(X_ErrH), np.array(YNorm_ErrH), np.array(YNorm_ErrL))

    return [graph_Total, graph_norm_Total]

def Draw_unc(graph, color):
    graph.SetFillColor(color)
    graph.SetFillStyle(3001)
    graph.SetLineColor(color)
    #graph.Draw("SAME2")
    graph.Draw("SAME2")

def DrawOnCanv(canv, var_name, plt_cfg, stacks, histos, scaled_sig, ratio_plot, legend, lumi_label, cms_label, total_unc, logY):
    #print var_name
    canv.SetBottomMargin(0.01)
    canv.cd()

    upper_pad = TPad("upperpad_"+var_name, "upperpad_"+var_name, 0,0.2, 1,1)
    upper_pad.SetBottomMargin(0.05)
    upper_pad.Draw()
    upper_pad.cd()

    if logY:
        upper_pad.SetLogy()
        stacks['all'].SetMinimum(1e-3)
        stacks['all'].SetMaximum(1e5)

        histos['data'].SetMinimum(1e-3)
        histos['data'].SetMaximum(1e5)

    if histos['data'].GetMaximum() > stacks['all'].GetMaximum():
        h_max = histos['data'].GetMaximum()
    else:
        h_max = stacks['all'].GetMaximum()
    if h_max < stacks['sig'].GetMaximum():
        h_max = stacks['sig'].GetMaximum()

    histos['data'].SetMaximum(h_max*1.3)
    stacks['all'].SetMaximum(h_max*1.3)

    histos['data'].Draw('PE')
    histos['data'].GetXaxis().SetLabelSize(0)
    histos['data'].GetYaxis().SetTitle('Events / (%.2f GeV)' %histos['data'].GetBinWidth(1))
    histos['data'].GetYaxis().SetTitleSize(0.045)
    stacks['all'].Draw('HISTSAME')

    for sample in plt_cfg.ana_cfg.sig_names:
        scaled_sig[sample].Draw('HISTSAME')

    histos['data'].SetMarkerStyle(20)
    #histos['data'].Draw('SAMEPE')
    #histos['data'].Draw("Axissame")
    
    ### Draw the uncertainties
    global stat_err, stat_err_norm
    stat_err,  stat_err_norm= Get_StatUnc(stacks['bkg'].GetStack().Last())
    #Draw_unc(stat_err, kGray+10)
    #legend.AddEntry(stat_err,"stat.","f")
    Draw_unc(total_unc[0], kGray+10)
    Draw_unc(stat_err, kRed-10)

    histos['data'].Draw('SAMEPE')
    histos['data'].Draw("Axissame")

    legend.AddEntry(stat_err,"stat.","f")
    legend.AddEntry(total_unc[0],"sys.","f")
    
    legend.Draw()

    # CMS style
    CMS_lumi.cmsText = "CMS"
    #CMS_lumi.extraText = "Preliminary"
    CMS_lumi.extraText = "Private"
    CMS_lumi.cmsTextSize = 0.95
    CMS_lumi.outOfFrame = True
    CMS_lumi.CMS_lumi(canv,4,11,plt_cfg.year)

    canv.cd()
    lower_pad = TPad("lowerpad_"+var_name, "lowerpad_"+var_name, 0, 0.01, 1,0.22)
    lower_pad.SetTopMargin(0.)
    lower_pad.SetBottomMargin(0.3)
    lower_pad.SetGridy()
    lower_pad.Draw()
    lower_pad.cd()
    ratio_plot.Draw("APSAME")

    Draw_unc(total_unc[1], kGray+10)
    Draw_unc(stat_err_norm, kRed-10)
    #ratio_plot.Draw("SAME")
    

def SaveCanvPic(canv, save_dir, save_name):
    canv.cd()
    #canv.SaveAs(save_dir + '/' + save_name + '.pdf')
    canv.SaveAs(save_dir + '/' + save_name + '.png')
