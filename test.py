from ROOT import *
import numpy as np

def Get_Unc(hist):

    BinTotal = hist.GetNbinsX()
    WidthBin = hist.GetBinWidth(1)

    XMean=[]
    YMean=[]
    X_ErrH=[]
    X_ErrL=[]
    Y_ErrH=[]
    Y_ErrL=[]

    xaxis = hist.GetXaxis()

    graph = TGraphAsymmErrors()

    for iBin in range(1, BinTotal+1):
        N = hist.GetBinContent(iBin)
        NErr = hist.GetBinError(iBin)
        XMean.append(xaxis.GetBinCenter(iBin))
        YMean.append(N)
        X_ErrH.append(0.5*WidthBin)
        X_ErrL.append(0.5*WidthBin)
        Y_ErrH.append(NErr)
        Y_ErrL.append(NErr)

    #print YMean
    #print Y_ErrH

    graph = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMean), np.array(X_ErrL), np.array(X_ErrH), np.array(Y_ErrL), np.array(Y_ErrH))

    return graph

file = TFile("plots_run2/ALP_plot_run2_param_SR.root")
canvas = TCanvas('a','a',1000,800)
hist = file.raw_plots.Get("Z_m_DYJetsToLL")

hist.Draw("HIST")

g = Get_Unc(hist)

g.SetFillColor(kRed)
g.SetFillStyle(3001)
g.SetLineColor(kRed)

g.Draw("SAME2")

canvas.SaveAs("canvas.png")

