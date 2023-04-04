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
    resolution = {'16':{'ele':[2.69, 2.41, 2.27, 2.20, 2.46, 2.15, 2.35, 2.47, 2.26, 2.31, 2.44, 2.34, 2.21, 2.19], 
                         'mu':[2.26, 1.51, 1.60, 1.48, 1.46, 1.62, 1.62, 1.67, 1.61, 1.62, 1.48, 1.48, 1.51, 1.51]}, 
               '16APV':{'ele':[3.05, 2.31, 2.32, 2.35, 2.38, 2.19, 2.46, 2.19, 2.32, 2.23, 2.02, 2.13, 2.18, 2.21], 
                         'mu':[2.42, 1.68, 1.58, 1.63, 1.71, 1.67, 1.69, 1.62, 1.67, 1.57, 1.36, 1.49, 1.51, 1.32]}, 
                  '17':{'ele':[4.08, 2.24, 2.15, 2.18, 2.31, 2.26, 2.41, 2.36, 2.20, 2.18, 2.27, 2.07, 2.24, 2.11], 
                         'mu':[2.29, 1.50, 1.39, 1.63, 1.59, 1.47, 1.50, 1.56, 1.66, 1.52, 1.59, 1.32, 1.47, 1.36]}, 
                  '18':{'ele':[3.26, 2.38, 2.24, 2.23, 2.28, 2.17, 2.42, 2.13, 2.51, 2.13, 2.38, 2.19, 2.09, 1.97], 
                         'mu':[3.09, 1.58, 1.62, 1.44, 1.53, 1.40, 1.64, 1.60, 1.53, 1.61, 1.52, 1.48, 1.48, 1.47]}}
    


    for channel in ['ele','mu']:
        
        canvas = TCanvas(channel,channel,650,600)
        canvas.cd()
        SetgStyle()
        canvas.SetRightMargin(0.1)
        canvas.SetLeftMargin(0.13)
        canvas.SetTopMargin(0.085)
        canvas.SetBottomMargin(0.12)
        colors = [2, 6, 3, 4]
        legend = TLegend(0.6,0.65,0.88,0.88)

        

        gr = {}

        resolu = {}

        for year in years:
            if '16APV' == year: continue

            resolu[year] = array('d')

            if '16' in year:
                for i in range(len(mass_list)):
                    resolu[year].append((resolution['16'][channel][i]+resolution['16APV'][channel][i])/2.0)
            else:
                for r in resolution[year][channel]:
                    resolu[year].append(r)

        marks = [20, 23, 21, 22]

        for year in years:
            
            if '16APV' == year: continue

            gr[year] = TGraph( len(mass_list), m_x, resolu[year] )
            gr[year].SetLineColor( colors[years.index(year)] )
            gr[year].SetLineWidth( 2 )
            gr[year].SetMarkerColor( colors[years.index(year)] )
            gr[year].SetMarkerStyle( marks[years.index(year)] )
            #gr[year].SetMarkerSize(2)
            gr[year].SetTitle('')
            gr[year].GetHistogram().SetMaximum(5.2)
            gr[year].GetHistogram().SetMinimum(1.1)
            gr[year].GetXaxis().SetTitle( 'm_{a} (GeV)' )
            #gr[year].GetXaxis().SetLabelSize(0.045)
            #gr[year].GetXaxis().SetTitleSize(0.045)
            gr[year].GetYaxis().SetTitle( 'Resolution' )

            legend.AddEntry(gr[year],"20"+year,"lp")

            
        
        gr['16'].Draw("ALP")
        #gr['16APV'].Draw("CP")
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
        
        canvas.SaveAs('./interpolation/resolution_'+channel+'.png')
        canvas.Close()
        



        


    print '\n\n'

    print 'Done'


main()
