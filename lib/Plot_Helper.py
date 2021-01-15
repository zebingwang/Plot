import os
import Plot_Configs as PC
import Analyzer_Configs as AC
from ROOT import *
import ROOT

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
    stacks['all']  = THStack("h_stack_"+var_name, var_name)
    stacks['sig']  = THStack("h_stack_"+var_name, var_name)
    stacks['bkg']  = THStack("h_stack_"+var_name, var_name)

    for sample in ana_cfg.samp_names:
        stacks[sample] = THStack("h_stack_"+var_name, var_name)

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
    ratio_plot.SetMinimum(0.2)
    ratio_plot.SetMaximum(3.4)
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

def DrawOnCanv(canv, var_name, plt_cfg, stacks, histos, scaled_sig, ratio_plot, legend, lumi_label, cms_label, logY):

    canv.SetBottomMargin(0.01)
    canv.cd()

    upper_pad = TPad("upperpad_"+var_name, "upperpad_"+var_name, 0,0.2, 1,1)
    upper_pad.SetBottomMargin(0.05)
    upper_pad.Draw()
    upper_pad.cd()

    if logY:
        upper_pad.SetLogy()
        stacks['all'].SetMinimum(1e-3)
        stacks['all'].SetMaximum(1e8)

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
    histos['data'].Draw('SAMEPE')
    histos['data'].Draw("Axissame")

    legend.Draw()
    #cms_label.DrawLatexNDC(0.1, 0.95, '#scale[1.5]{CMS} H To Z + ALP')
    #cms_label.Draw('same')
    #lumi_label.DrawLatexNDC(0.90, 0.91, '%s fb^{-1} (13 TeV)' %plt_cfg.lumi)
    #lumi_label.Draw('same')
    # CMS style
    CMS_lumi.cmsText = "CMS"
    #CMS_lumi.extraText = "Preliminary"
    CMS_lumi.extraText = "Private"
    CMS_lumi.cmsTextSize = 0.95
    CMS_lumi.outOfFrame = True
    CMS_lumi.CMS_lumi(canv,4,11)

    canv.cd()
    lower_pad = TPad("lowerpad_"+var_name, "lowerpad_"+var_name, 0, 0.01, 1,0.22)
    lower_pad.SetTopMargin(0.)
    lower_pad.SetBottomMargin(0.3)
    lower_pad.SetGridy()
    lower_pad.Draw()
    lower_pad.cd()
    ratio_plot.Draw("AP")


def SaveCanvPic(canv, save_dir, save_name):
    canv.cd()
    #canv.SaveAs(save_dir + '/' + save_name + '.pdf')
    canv.SaveAs(save_dir + '/' + save_name + '.png')
