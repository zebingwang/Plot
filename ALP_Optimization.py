####################################################
####################################################

import os
import sys
import numpy as np
from array import array
import math

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
import Analyzer_Configs as AC
#import Plot_Configs     as PC

#from Analyzer_ALP import PIso2D, plot2D, plot2D_CONT

import CMS_lumi, tdrstyle

from xgboost import XGBClassifier
import pickle

import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

#import argparse
from optparse import OptionParser

#parser = argparse.ArgumentParser(description="A simple ttree plotter")
#parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
#parser.add_option("-c", "--nCats",   dest="nCats",    default=5,    type="int",    help="nCats"  )
#parser.add_argument('-b', '--blind', dest='blind', action='store_true', default=False, help='Blind signal region')

parser = OptionParser()
parser.add_option("-y", "--Year", dest="year", default='run2', type="str", help="which year's datasetes")
parser.add_option("-c", "--nCats",   dest="nCats",    default=5,    type="int",    help="nCats"  )
parser.add_option('-o', "--outDir", dest='outDir', default="./optimize", type="string", help="outDir")
parser.add_option('-p', '--plot', dest='plot', action='store_true', default=False, help='Plot?')
parser.add_option('-s', '--sigma', dest='sigma', action='store_true', default=False, help='mass region?')
parser.add_option('--sigVSscore', dest='sigVSscore', action='store_true', default=False, help='mass region?')
parser.add_option('--doOpt', dest='doOpt', action='store_true', default=False, help='whether optimize the category?')
parser.add_option('--ele', dest='ele', action='store_true', default=False, help='electron channel?')
parser.add_option('--mu', dest='mu', action='store_true', default=False, help='muon channel?')
(options, args) = parser.parse_args()



gROOT.SetBatch(True)
tdrstyle.setTDRStyle()


if options.year == '2016':
    file_out = 'plots_16'
    name_SR = "ALP_plot_data16_param_SR.root"
    name_CR = "ALP_plot_data16_param_CR.root"
elif options.year == '2017':
    file_out = 'plots_17'
    name_SR = "ALP_plot_data17_param_SR.root"
    name_CR = "ALP_plot_data17_param_CR.root"
elif options.year == '2018':
    file_out = 'plots_18'
    name_SR = "ALP_plot_data18_param_SR.root"
    name_CR = "ALP_plot_data18_param_CR.root"
elif options.year == 'run2':
    file_out = 'plots_run2UL'
    if options.ele:
        name_SR = "ALP_plot_run2_UL_SR_ele.root"
        name_CR = "ALP_plot_run2_UL_CR_ele.root"
    elif options.mu:
        name_SR = "ALP_plot_run2_UL_SR_mu.root"
        name_CR = "ALP_plot_run2_UL_CR_mu.root"
    else:
        name_SR = "ALP_plot_run2_UL_SR.root"
        name_CR = "ALP_plot_run2_UL_CR.root"

        #name_SR = "ALP_plot_run2_UL_onlyM10_SR.root"
        #name_CR = "ALP_plot_run2_UL_onlyM10_CR.root"
        
        #name_SR = "ALP_plot_run2_UL_Ns_10_SR.root"
        #name_CR = "ALP_plot_run2_UL_Ns_10_CR.root"
    
else:
    print "do not include at 2016/2017/2018"
    exit(0)


def SetgStyle():

    gStyle.SetFrameFillColor(0)
    gStyle.SetStatColor(0)
    gStyle.SetOptStat(0)
    gStyle.SetTitleFillColor(0)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetFrameBorderMode(0)
    gStyle.SetPadColor(kWhite)
    gStyle.SetCanvasColor(kWhite)
    
    
    gStyle.SetCanvasDefH(600) #Height of canvas
    gStyle.SetCanvasDefW(600) #Width of canvas
    gStyle.SetCanvasDefX(0)   #POsition on screen
    gStyle.SetCanvasDefY(0)

    
    gStyle.SetPadLeftMargin(0.13)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetPadTopMargin(0.085)
    gStyle.SetPadBottomMargin(0.12)
    
    # For hgg axis titles:
    gStyle.SetTitleColor(1, "XYZ")
    gStyle.SetTitleFont(42, "XYZ")
    gStyle.SetTitleSize(0.04, "XYZ")
    gStyle.SetTitleXOffset(0.95)#//0.9)
    gStyle.SetTitleYOffset(1.15)# // => 1.15 if exponents
    
    # For hgg axis labels:
    gStyle.SetLabelColor(1, "XYZ")
    gStyle.SetLabelFont(42, "XYZ")
    gStyle.SetLabelOffset(0.007, "XYZ")
    gStyle.SetLabelSize(0.04, "XYZ")
    
    # Legends
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendFillColor(kWhite)
    gStyle.SetLegendFont(42)
    
    gStyle.SetFillColor(10)
    # Nothing for now
    gStyle.SetTextFont(42)
    gStyle.SetTextSize(0.03)


def open_file(file_name):
    in_file = TFile( file_name , "r")
    return in_file

def get_hist(file, mvaVal_name):
    hist = file.raw_plots.Get(mvaVal_name)
    return hist

def hist2graph(hist, mva_low = 0.1):
    
    Nbins = hist.GetNbinsX()
    bin_x_low = hist.FindBin(mva_low)
    xaxis = hist.GetXaxis()

    bin_x_Center=[]
    bin_y=[]
    for x in range(Nbins+1):
        # remove BDT score less than mva_low
        if x < bin_x_low: continue
        bin_x_Center.append(xaxis.GetBinCenter(x))
        bin_y.append(hist.GetBinContent(x))

    graph = TGraph(Nbins-bin_x_low+1, np.array(bin_x_Center), np.array(bin_y))

    return graph

def smooth(graph, hist, mva_low = 0.1):
    h_smooth = TH1F(hist.GetName()+"_smooth",hist.GetName()+"_smooth",hist.GetNbinsX(),0.0,1.0)
    h_smooth_up = TH1F(hist.GetName()+"_smooth_up",hist.GetName()+"_smooth_up",hist.GetNbinsX(),0.0,1.0)
    h_smooth_dn = TH1F(hist.GetName()+"_smooth_dn",hist.GetName()+"_smooth_dn",hist.GetNbinsX(),0.0,1.0)

    smoother = TGraphSmooth()
    g_smooth = smoother.SmoothSuper(graph)

    xaxis = hist.GetXaxis()
    x = array('d', [0])
    y = array('d', [0])

    for i in range(1, hist.GetNbinsX()+1):

        h_x = xaxis.GetBinCenter(i)

        if i < hist.FindBin(mva_low):
            y[0] = 0.0
        else:
            g_smooth.GetPoint(i,x,y)

        h_smooth.SetBinContent(i,y[0])

    for i in range(1, hist.GetNbinsX()+1):

        if i < hist.FindBin(mva_low):
            y[0] = 0.0
        else:
            y = h_smooth.GetBinContent(i)

        #print y
        #print np.sqrt(y)

        if y>=0.:
            h_smooth_up.SetBinContent(i,y+np.sqrt(y))
            if (y-np.sqrt(y))>0.: h_smooth_dn.SetBinContent(i,y-np.sqrt(y))
            else: h_smooth_dn.SetBinContent(i,0.)
        else:
            h_smooth_up.SetBinContent(i,0.)
            h_smooth_dn.SetBinContent(i,0.)


            

    return [h_smooth, h_smooth_up, h_smooth_dn]

def compare(hist, hist_smooth, ma):
    canv = TCanvas("cc", "cc", 650, 600)
    canv.cd()
    canv.SetLogy()
    SetgStyle()
    '''
    canv.SetRightMargin(0.1)
    canv.SetLeftMargin(0.13)
    canv.SetTopMargin(0.085)
    canv.SetBottomMargin(0.12)
    '''

    hist.SetMinimum(1e-1)
    hist.SetMaximum(1e4)
    hist_smooth[0].SetMinimum(1e-1)
    hist_smooth[0].SetMaximum(1e4)

    hist.SetFillColor(10)
    hist.GetXaxis().SetTitle('BDT Output')
    hist.GetYaxis().SetTitle('Events')

    hist.GetYaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetLabelColor(1)
    hist.GetYaxis().SetLabelFont(42)
    hist.GetYaxis().SetLabelOffset(0.007)
    hist.GetYaxis().SetTitleColor(1)
    hist.GetYaxis().SetTitleFont(42)
    hist.GetYaxis().SetTitleSize(0.05)

    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetXaxis().SetLabelColor(1)
    hist.GetXaxis().SetLabelFont(42)
    hist.GetXaxis().SetLabelOffset(0.007)
    hist.GetXaxis().SetTitleColor(1)
    hist.GetXaxis().SetTitleFont(42)
    hist.GetXaxis().SetTitleSize(0.05)


    hist.Draw("HIST")



    hist_smooth[1].SetLineColor(kGreen)
    hist_smooth[1].SetLineWidth(2)
    hist_smooth[1].Draw("SAME")

    hist_smooth[2].SetLineColor(kBlue)
    hist_smooth[2].SetLineWidth(2)
    hist_smooth[2].Draw("SAME")

    hist_smooth[0].SetLineColor(kRed)
    hist_smooth[0].SetLineWidth(2)
    hist_smooth[0].Draw("SAME")

    #canv.cd()
    global legend
    legend = TLegend(0.6,0.65,0.88,0.88)
    #legend.AddEntry(hist, "Background")
    legend.AddEntry(hist_smooth[1], "Smoothing + 1#sigma", "l")
    legend.AddEntry(hist_smooth[0], "Smoothing", "l")
    legend.AddEntry(hist_smooth[2], "Smoothing - 1#sigma", "l")
    legend.Draw("SAME")

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextColor(kBlack)
    latex.SetTextFont(42)
    latex.SetTextSize(0.042)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.45,0.8, r"m_{a} = %s GeV" % (ma.lstrip("M")))


    canv.Update()
    #print legend

    #cmsText     = "CMS";
    #cmsTextFont   = 61
    #latex = rt.TLatex()

    

    return canv


def GetHist(file, var):

    hist = get_hist(file, var)
    
    graph = hist2graph(hist, 0.01)
    hist_smooth = smooth(graph, hist, 0.01)

    return hist, hist_smooth

def computeSignificance(s,b,d_noSmooth):
  significance = -999.
  if b>0 and s>0.: significance = (2*(s+b)*math.log(1+(s/b))) - 2*s
  if significance>0. and d_noSmooth>=8. and b>0.: return np.sqrt((2*(s+b)*math.log(1+(s/b))) - 2*s)
  else: return -999.

def sumSignificance(partition, h_sig_SR, h_bkg_SR,h_data_SB_noSmooth):
  sum = 0.
  for pair in partition:
    s = 2*h_sig_SR.Integral(pair[0],pair[1])
    b = h_bkg_SR.Integral(pair[0],pair[1])
    d_noSmooth = h_data_SB_noSmooth.Integral(pair[0],pair[1])
    significance = computeSignificance(s,b,d_noSmooth)
    #print h_sig_SR.GetBinCenter(pair[0])-h_bdt_signal_SR.GetBinWidth(pair[0])/2.,significance,b
    if significance>0.: sum += significance*significance
    else: return -999.
  return np.sqrt(sum)

def getResults(h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB, nCats, nBins):

    significance_final = -999.
    partition_final = []
    sig_all = {}

    #1 categories
    if nCats == 1:

        for i in range(1,nBins+1):
            partition = [[i,nBins]]
            significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
            #print h_bdt_signal_SR.GetBinCenter(partition[0][0])-h_bdt_signal_SR.GetBinWidth(partition[0][0])/2,"1. --->",significance,h_bdt_data_SB.Integral(partition[0][0],nBins)
            #print h_bdt_signal_SR.GetBinCenter(partition[0][0])-h_bdt_signal_SR.GetBinWidth(partition[0][0])/2,"1. --->",significance
            if significance>significance_final:
                significance_final = significance
                partition_final = partition
            sig_all[i] = significance
        output = nCats," - Best category: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2,"1. --->",significance_final,"signal total: ",h_bdt_signal_SR.Integral(1,nBins+1),"events:",h_bdt_signal_SR.GetEntries(),"signal cut",h_bdt_signal_SR.Integral(partition_final[0][0],nBins),"smoothed background: ",h_bdt_datamix_SR_weighted_smooth.Integral(partition_final[0][0],nBins),"data: ",h_bdt_data_SB.Integral(partition_final[0][0],nBins),"NEXT BIN: ","smoothed background: ",h_bdt_datamix_SR_weighted_smooth.Integral(partition_final[0][0]+1,nBins),"data: ",h_bdt_data_SB.Integral(partition_final[0][0]+1,nBins)
        print ' '.join(map(str,output))

    #2 categories
    elif nCats == 2:

        for i in range(1,nBins+1):
            for j in range(i+1,nBins+1):
                partition = [[1,i],[j,nBins]]
                if abs(i-j)==1:
                    #print partition
                    significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
                    if significance>significance_final:
                        significance_final = significance
                        partition_final = partition

        output = nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2,h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2,"1. --->",significance_final
        print ' '.join(map(str,output))

    #3 categories
    elif nCats == 3:

        for i in range(1,nBins+1):
            for j in range(i+1,nBins+1):
                for k in range(j+1,nBins+1):
                    partition = [[1,i],[j,k-1],[k,nBins]]
                    if abs(i-j)==1:
                        #print partition
                        significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
                        if significance>significance_final:
                            significance_final = significance
                            partition_final = partition

        output = nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, "1. --->",significance_final
        print ' '.join(map(str,output))

    #4 categories
    elif nCats == 4:

        for i in range(1,nBins+1):
            for j in range(i+1,nBins+1):
                for k in range(j+1,nBins+1):
                    for d in range(k+1,nBins+1):
                        partition = [[1,i],[j,k-1],[k,d-1],[d,nBins]]
                        if abs(i-j)==1:
                            #print partition
                            significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
                            if significance>significance_final:
                                significance_final = significance
                                partition_final = partition

        output = nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[3][0])-h_bdt_signal_SR.GetBinWidth(partition_final[3][0])/2, "1. --->",significance_final
        print ' '.join(map(str,output))

    #5 categories
    elif nCats == 5:
        for i in range(1,nBins+1):
            for j in range(i+1,nBins+1):
                for k in range(j+1,nBins+1):
                    for d in range(k+1,nBins+1):
                        for f in range(d+1,nBins+1):
                            partition = [[1,i],[j,k-1],[k,d-1],[d,f-1],[f,nBins]]
                            if abs(i-j)==1:
                            #print partition
                                significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
                                if significance>significance_final:
                                    significance_final = significance
                                    partition_final = partition

        output = nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[3][0])-h_bdt_signal_SR.GetBinWidth(partition_final[3][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[4][0])-h_bdt_signal_SR.GetBinWidth(partition_final[4][0])/2, "1. --->",significance_final
        print ' '.join(map(str,output))

    else:
        print "Number of categories not supported, choose: 1, 2, 3, 4, 5, 6, 7, 8 or 9!"
        sys.exit()

    
    return partition_final, significance_final, output, sig_all



def main():

    #analyzer_cfg.sig_names
    analyzer_cfg = AC.Analyzer_Config('inclusive',options.year)

    file = {}
    file['SR'] = open_file(file_out + '/' + name_SR)
    file['CR'] = open_file(file_out + '/' + name_CR)


    signal_region = ['1sigma', '1P5sigma', '2sigma', '3sigma','all']
    hist_SR = {}
    hist_SR_smooth = {}

    hist_CR = {}
    hist_CR_smooth = {}

    hist_signal = {}

    for r in signal_region:
        hist_SR[r] = {}
        hist_SR_smooth[r] = {}

        hist_signal[r] = {}

    y = {}
    y_list = ['significance', 'boundary', 'Nsignal', 'Nbackground']
    for l in y_list:
        y[l] = {}

    for sig in analyzer_cfg.sig_names:

        #hist_signal[sig] = get_hist(file['SR'], "mvaVal_"+sig+"_"+sig)

        ###### smoothing his
        #hist_SR[sig], hist_SR_smooth[sig] = GetHist(file['SR'], "mvaVal_"+sig+"_DYJetsToLL")
        #hist_SR_1sigma[sig], hist_SR_1sigma_smooth[sig] = GetHist(file['SR'], "mvaVal_"+sig+"_DYJetsToLL")
        #hist_SR_1P5sigma[sig], hist_SR_1P5sigma_smooth[sig] = GetHist(file['SR'], "mvaVal_"+sig+"_DYJetsToLL")
        #hist_SR_2sigma[sig], hist_SR_2sigma_smooth[sig] = GetHist(file['SR'], "mvaVal_"+sig+"_DYJetsToLL")
        #hist_SR_3sigma[sig], hist_SR_3sigma_smooth[sig] = GetHist(file['SR'], "mvaVal_"+sig+"_DYJetsToLL")
        for r in signal_region:
            if r == 'all':
                hist_signal[r][sig] = get_hist(file['SR'], "mvaVal_"+sig+"_"+sig)
                hist_SR[r][sig], hist_SR_smooth[r][sig] = GetHist(file['SR'], "mvaVal_"+sig+"_DYJetsToLL")
            else:
                hist_signal[r][sig] = get_hist(file['SR'], "mvaVal_"+r+"_"+sig+"_"+sig)
                hist_SR[r][sig], hist_SR_smooth[r][sig] = GetHist(file['SR'], "mvaVal_"+r+"_"+sig+"_DYJetsToLL")
        hist_CR[sig], hist_CR_smooth[sig] = GetHist(file['CR'], "mvaVal_"+sig+"_data")

        if options.plot:
            for r in signal_region:
                canv_SR = compare(hist_SR[r][sig], hist_SR_smooth[r][sig],sig)
                canv_SR.SaveAs(options.outDir+"/"+sig+"_"+r+"_DYJetsToLL_SR.png")

            canv_CR = compare(hist_CR[sig], hist_CR_smooth[sig],sig)
            canv_CR.SaveAs(options.outDir+"/"+sig+"_data_CR.png")

        if options.doOpt:
            ###### category optimization
            nCats = options.nCats
            nBins = hist_signal['all'][sig].GetNbinsX()

            partition_final = {}
            significance_final = {}
            output = {}
            significance_all = {}

            partition_final_up = {}
            significance_final_up = {}
            output_up = {}
            significance_all_up = {}
            partition_final_dn = {}
            significance_final_dn = {}
            output_dn = {}
            significance_all_dn = {}

            for r in signal_region:

                partition_final[r], significance_final[r], output[r], significance_all[r] = getResults(hist_signal[r][sig], hist_SR_smooth[r][sig][0], hist_CR[sig], nCats, nBins)
                partition_final_up[r], significance_final_up[r], output_up[r], significance_all_up[r] = getResults(hist_signal[r][sig], hist_SR_smooth[r][sig][1], hist_CR[sig], nCats, nBins)
                partition_final_dn[r], significance_final_dn[r], output_dn[r], significance_all_dn[r] = getResults(hist_signal[r][sig], hist_SR_smooth[r][sig][2], hist_CR[sig], nCats, nBins)
                #partition_final[r], significance_final[r], output[r], significance_all[r] = getResults(hist_signal[r][sig], hist_SR_smooth[r][sig][1], hist_CR[sig], nCats, nBins)
                #partition_final[r], significance_final[r], output[r], significance_all[r] = getResults(hist_signal[r][sig], hist_SR_smooth[r][sig][2], hist_CR[sig], nCats, nBins)

                #print "nSig: ", 

                out_file_name = options.outDir + '/categorize_' + r + '_' + sig + '.txt'

                if not os.path.exists(options.outDir):
                    os.makedirs(options.outDir)

                outfile = open(out_file_name, 'a')
                outfile.write(' '.join(map(str,output[r]))+'\n')

                outfile_up = open(options.outDir + '/categorize_' + r + '_' + sig + '_up.txt', 'a')
                outfile_up.write(' '.join(map(str,output_up[r]))+'\n')

                outfile_dn = open(options.outDir + '/categorize_' + r + '_' + sig + '_dn.txt', 'a')
                outfile_dn.write(' '.join(map(str,output_dn[r]))+'\n')

        
        if options.sigVSscore:
            
            if nCats == 1:
                
                w = hist_signal['all']['M1'].GetBinWidth(1)
                '''
                plt.plot([(i-1)*w for i in significance_all['all'].keys()], significance_all['all'].values(), "o-", markersize=2., label='Significance')
                #plt.plot([(i-1)*w for i in significance_all['all'].keys()], significance_all_up['all'].values(), "o-", c='red')
                #plt.plot([(i-1)*w for i in significance_all['all'].keys()], significance_all_dn['all'].values(), "o-", c='green')
                plt.fill_between([(i-1)*w for i in significance_all['all'].keys()], significance_all_up['all'].values(), significance_all_dn['all'].values(), alpha=0.2, color='red',label=r'Significance$\pm 1\sigma$')
                plt.ylim([0.0, 150.0])
                plt.xlim([0.5, 1.0])
                plt.xlabel('boundary_'+sig)
                plt.ylabel('significance')
                plt.legend(loc='best')
                plt.plot([(partition_final['all'][0][0]-1)*w, (partition_final['all'][0][0]-1)*w], [0., 150], c='black', linestyle='--')
                plt.grid()
                plt.savefig(options.outDir+ '/sigVSscore_' + sig + '.png')
                plt.close('all')
                '''

                canv_sigVSscore = TCanvas("cc1", "cc1", 650, 600)
                canv_sigVSscore.cd()
                

                num_arr = len(significance_all['all'].values())
                x_arr = array('d')
                y_arr = array('d')
                y_arr_err = array('d')
                for i in significance_all['all'].keys():
                    x_arr.append((i-1)*w)
                    y_arr.append(significance_all['all'][i])
                gr = TGraph( num_arr, x_arr, y_arr)
                gr_err = TGraph( 2 * num_arr)

                for i in range(num_arr):
                    gr_err.SetPoint(i, x_arr[i], significance_all_up['all'].values()[i])
                    gr_err.SetPoint(2*num_arr-1-i, x_arr[i], significance_all_dn['all'].values()[i])
                
                gr_err.SetFillColorAlpha(2,0.3)
                gr_err.GetHistogram().SetMaximum(100.0)
                gr_err.GetHistogram().SetMinimum(0.0)
                gr_err.GetYaxis().SetTitle( 'Significance (AMS)' )
                gr_err.GetXaxis().SetTitle( 'BDT select boundary' )
                gr_err.GetYaxis().SetLabelSize(0.04)
                gr_err.GetYaxis().SetLabelColor(1)
                gr_err.GetYaxis().SetLabelFont(42)
                gr_err.GetYaxis().SetLabelOffset(0.007)
                gr_err.GetYaxis().SetTitleColor(1)
                gr_err.GetYaxis().SetTitleFont(42)
                gr_err.GetYaxis().SetTitleSize(0.05)

                gr_err.GetXaxis().SetLabelSize(0.04)
                gr_err.GetXaxis().SetLabelColor(1)
                gr_err.GetXaxis().SetLabelFont(42)
                gr_err.GetXaxis().SetLabelOffset(0.007)
                gr_err.GetXaxis().SetTitleColor(1)
                gr_err.GetXaxis().SetTitleFont(42)
                gr_err.GetXaxis().SetTitleSize(0.05)

                SetgStyle()
                gr.SetMarkerSize(0.5)

                gr_err.Draw('AF')
                gr.Draw('LP')
                Line = TLine((partition_final['all'][0][0]-1)*w, 0, (partition_final['all'][0][0]-1)*w, 100)
                #SetLineStyle
                Line.Draw("SAME")

                canv_sigVSscore.SaveAs(options.outDir+ '/sigVSscore_' + sig + '.png')

        #make sigma plots
        if nCats == 1:
            for l in y_list:
                y[l][sig] = []
            for r in signal_region:
                
                s = 2*hist_signal['all'][sig].Integral(partition_final[r][0][0],nBins)
                b = hist_SR_smooth['all'][sig][0].Integral(partition_final[r][0][0],nBins)

                if b > 0.0:
                    y['boundary'][sig].append(hist_signal[r][sig].GetBinCenter(partition_final[r][0][0])-hist_signal[r][sig].GetBinWidth(partition_final[r][0][0])/2)
                    #y['significance'][sig].append(significance_final[r])

                    y['significance'][sig].append(np.sqrt((2*(s+b)*math.log(1+(s/b))) - 2*s))
                    y['Nsignal'][sig].append(hist_signal['all'][sig].Integral(partition_final[r][0][0],nBins))
                    y['Nbackground'][sig].append(hist_SR_smooth['all'][sig][0].Integral(partition_final[r][0][0],nBins))
                else:
                    b = hist_SR_smooth['all'][sig][0].Integral(partition_final['all'][0][0],nBins)

                    y['boundary'][sig].append(hist_signal['all'][sig].GetBinCenter(partition_final['all'][0][0])-hist_signal['all'][sig].GetBinWidth(partition_final['all'][0][0])/2)

                    y['significance'][sig].append(np.sqrt((2*(s+b)*math.log(1+(s/b))) - 2*s))
                    y['Nsignal'][sig].append(hist_signal['all'][sig].Integral(partition_final['all'][0][0],nBins))
                    y['Nbackground'][sig].append(hist_SR_smooth['all'][sig][0].Integral(partition_final['all'][0][0],nBins))



    if options.sigma and nCats == 1:
        for l in y_list:
            colors = ['red','orange','blue','purple','darkturquoise', 'yellow', 'pink', 'green', 'black', 'gray', 'violet', 'gold', 'cyan', 'brown']

            plt.xlabel('signal mass region')
            plt.ylabel(l)

            for sig in analyzer_cfg.sig_names:

                plt.plot(signal_region, y[l][sig], c=colors[analyzer_cfg.sig_names.index(sig)], label=sig)


            plt.grid()
            plt.legend()
            plt.savefig(options.outDir+ '/' + l + '.png')
            plt.close('all')


main()
