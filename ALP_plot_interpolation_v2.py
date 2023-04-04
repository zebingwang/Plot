####################################################
####################################################

import os
import sys
import numpy as np

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from Plot_Helper import LoadNtuples, MakeStack, CreateCanvas, DrawOnCanv, SaveCanvPic, MakeLumiLabel, MakeCMSDASLabel, ScaleSignal, MakeRatioPlot, MakeLegend, Total_Unc
from Analyzer_Helper import MuPair, DimuCandidates, IsBtagLoose, IsBtagMed, IsBtagTight, CountYield
import Analyzer_Configs as AC
import Plot_Configs     as PC

import CMS_lumi, tdrstyle

from xgboost import XGBClassifier
import pickle
import copy 
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy import interpolate

from array import array


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
    gStyle.SetTitleSize(0.05, "XYZ")
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

   



def main():

    mass_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
    m_x = array('d')
    for m in mass_list:
        m_x.append(m)
    mass_interp = range(1,31)
    years = ['16', '16APV', '17', '18']
    
    path_basic_num = '/publicfs/cms/user/wangzebing/ALP/Analysis_code/fit/fit_run2_UL/'
    path_basic_dem = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/'
    eff = {}
    # v1
    #resolution = {'16':{'ele':[2.95, 2.47, 2.12, 2.03, 2.36, 2.27, 2.30, 2.49, 2.22, 2.47, 2.40, 2.69, 2.27, 2.11], 
    #                     'mu':[2.39, 1.45, 1.58, 1.55, 1.50, 1.61, 1.51, 1.55, 1.48, 1.51, 1.40, 1.43, 1.55, 1.43]}, 
    #           '16APV':{'ele':[3.77, 2.37, 2.32, 2.32, 2.20, 2.26, 2.29, 2.19, 2.13, 2.11, 2.11, 2.12, 2.01, 2.20], 
    #                     'mu':[2.49, 1.58, 1.69, 1.54, 1.76, 1.76, 1.57, 1.69, 1.65, 1.47, 1.41, 1.37, 1.44, 1.40]}, 
    #              '17':{'ele':[3.54, 2.38, 2.36, 2.40, 2.36, 2.37, 2.46, 2.21, 2.49, 2.33, 2.15, 2.25, 2.11, 2.34], 
    #                     'mu':[2.60, 1.45, 1.59, 1.64, 1.72, 1.42, 1.52, 1.46, 1.64, 1.59, 1.61, 1.39, 1.38, 1.46]}, 
    #              '18':{'ele':[3.89, 2.31, 2.14, 2.46, 2.29, 2.22, 2.45, 2.21, 2.35, 2.20, 2.10, 2.21, 2.36, 2.15], 
    #                     'mu':[3.22, 1.64, 1.51, 1.48, 1.68, 1.45, 1.68, 1.56, 1.69, 1.52, 1.38, 1.40, 1.34, 1.46]}}
    #v2
    resolution = {'16':{'ele':[2.90, 2.41, 2.40, 2.16, 2.25, 2.23, 2.35, 2.47, 2.18, 2.43, 2.23, 2.33, 2.12, 2.27], 
                         'mu':[2.28, 1.52, 1.64, 1.52, 1.62, 1.50, 1.49, 1.55, 1.67, 1.46, 1.40, 1.46, 1.55, 1.41]}, 
               '16APV':{'ele':[3.38, 2.43, 2.32, 2.36, 2.12, 2.49, 2.38, 2.14, 2.21, 2.00, 1.97, 2.03, 2.22, 2.23], 
                         'mu':[2.61, 1.66, 1.48, 1.56, 1.64, 1.61, 1.69, 1.63, 1.54, 1.45, 1.52, 1.53, 1.40, 1.50]}, 
                  '17':{'ele':[4.22, 2.46, 2.38, 2.18, 2.53, 2.53, 2.46, 2.33, 2.34, 2.20, 2.24, 2.17, 1.93, 2.34], 
                         'mu':[2.09, 1.45, 1.59, 1.61, 1.73, 1.56, 1.50, 1.57, 1.54, 1.56, 1.57, 1.58, 1.47, 1.50]}, 
                  '18':{'ele':[3.26, 2.32, 2.31, 2.46, 2.37, 2.35, 2.42, 2.22, 2.42, 2.33, 2.38, 2.29, 2.25, 2.05], 
                         'mu':[3.07, 1.69, 1.45, 1.48, 1.72, 1.49, 1.68, 1.70, 1.69, 1.62, 1.52, 1.50, 1.46, 1.37]}}
    

    for year in years:
        eff[year] = {}
        for channel in ['ele','mu']:
            eff[year][channel] = array('d')
            for mass in mass_list:
                file_name_num = path_basic_num + 'ALP_data_sig_Am' + str(mass) + '_Hm125_' + year + '_workspace_'+channel+'.root'
                file_name_dem = path_basic_dem + year + '/' + 'ALP_M' + str(mass) + '.root'

                file_num = TFile(file_name_num)
                file_dem = TFile(file_name_dem)
                
                #print file_name
                num_final = file_num.CMS_hza_workspace.data("ggh_125_13TeV_cat0").sumEntries()
                dem_final = file_dem.nEvents_weight_ntuple.GetBinContent(1) + file_dem.nEvents_weight_ntuple.GetBinContent(2)

                #print 'year: ', year, 'channel: ', channel, 'mass: ', mass, 'dem: ', dem_final, 'num: ', num_final

                eff[year][channel].append(num_final/dem_final*100.)

            #print "{0} ALP mass:{1}, num_final:{2}, eff:{3}".format(year, mass, num_final, eff)
    f = open('./interpolation/eff.txt','w')


    for channel in ['ele','mu']:
        
        canvas = TCanvas(channel,channel,650,600)
        canvas.cd()
        SetgStyle()
        canvas.SetRightMargin(0.1)
        canvas.SetLeftMargin(0.13)
        canvas.SetTopMargin(0.085)
        canvas.SetBottomMargin(0.12)
        colors = [2, 5, 3, 4]
        legend = TLegend(0.6,0.65,0.88,0.88)
        #gStyle.SetLegendFont(42)
        #gStyle.SetLegendTextSize(0.03)
        #legend.SetBorderSize(0)
        

        gr = {}
        marks = [20, 23, 21, 22]

        eff_16_mean = array('d')
        for i in range(len(mass_list)):
            eff_16_mean.append((eff['16'][channel][i] + eff['16APV'][channel][i])/2.0)

        for year in years:

            if '16APV' == year: continue

            if '16' in year:
                gr[year] = TGraph( len(mass_list), m_x, eff_16_mean )
            else:
                gr[year] = TGraph( len(mass_list), m_x, eff[year][channel] )
            gr[year].SetLineColor( colors[years.index(year)] )
            gr[year].SetLineWidth( 2 )
            gr[year].SetMarkerColor( colors[years.index(year)] )
            gr[year].SetMarkerStyle( marks[years.index(year)] )
            #gr[year].SetMarkerSize(2)
            gr[year].SetTitle('')
            gr[year].GetHistogram().SetMaximum(2.5)
            gr[year].GetHistogram().SetMinimum(0.2)
            gr[year].GetXaxis().SetTitle( 'm_{a} (GeV)' )
            #gr[year].GetXaxis().SetLabelSize(0.045)
            #gr[year].GetXaxis().SetTitleSize(0.045)
            gr[year].GetYaxis().SetTitle( 'Efficiency #times Acceptance (%)' )

            legend.AddEntry(gr[year],"20"+year,"lp")

            
        
        gr['16'].Draw("ALP")
        gr['17'].Draw("LP")
        gr['18'].Draw("LP")

        # CMS style
        CMS_lumi.cmsText = "CMS"
        CMS_lumi.extraText = "Preliminary"
        CMS_lumi.cmsTextSize = 0.67
        CMS_lumi.lumiTextSize = 0.6
        CMS_lumi.lumiText_posY = 0.0075
        CMS_lumi.lumiText_posX = 0.0
        CMS_lumi.CMSText_posX = -0.03
        CMS_lumi.outOfFrame = True
        CMS_lumi.CMS_lumi(canvas,4,0,'run2')

        legend.Draw('SAME')

        canvas.SaveAs('./interpolation/eff_'+channel+'.png')
        canvas.Close()
        



        


    print '\n\n'

    print 'Done'


main()
