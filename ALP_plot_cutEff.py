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

def getVariableHistsEventsNumber_weight(Tree,varName,cut):

    #Make a canvas to do the work on
    canvas = TCanvas(varName,varName,1000,800)

    #Extract the relevant variable from the trees
    Tree.Draw("{0}>>tree{0}".format(varName),"factor*pho1SFs*pho2SFs*({0})".format(cut))
    Hist = gDirectory.Get("tree{0}".format(varName))

    canvas.Clear()

    return Hist.Integral()

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
    gStyle.SetLegendTextSize(0.03)
    gStyle.SetFillColor(10)
    # Nothing for now
    gStyle.SetTextFont(42)
    gStyle.SetTextSize(0.03)


def main():

    mass_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
    m_x = array('d')
    for m in mass_list:
        m_x.append(m)

    years = ['16', '16APV', '17', '18']
    #years = ['16APV', '18']
    #years = ['run2']
    
    path_basic = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/'
    eff = {}
    num = {}

    #cuts_lable = ["PreLeptonSelection", "PreLeptonSelection + At least 2 photons", "PreLeptonSelection + At least 2 photons + PhotonSelection"]
    cuts_lable = ["Trigger + Z candidate selection", "#splitline{Trigger + Z candidate selection}{+ ALP candidate selection(At least 2 photons)}", "#splitline{Trigger + Z candidate selection}{#splitline{+ ALP candidate selection(At least 2 photons)}{+ Higgs candidate selection}}"]
    cuts =["lepton selection", "1", "passChaHadIso&&passNeuHadIso&&passHOverE&&passdR_gl&&(H_m>95&&H_m<180)"]

    #cuts_lable = ["Trigger + Z candidate selection", "Trigger + Z candidate selection + ALP candidate selection(At least 2 photons)", "Trigger + Z candidate selection + ALP candidate selection(At least 2 photons) \n + passChaHadIso", "Trigger + Z candidate selection + ALP candidate selection(At least 2 photons) \n + passChaHadIso+passNeuHadIso", "Trigger + Z candidate selection + ALP candidate selection(At least 2 photons)\n + passChaHadIso+passNeuHadIso+passHOverE", "Trigger + Z candidate selection + ALP candidate selection(At least 2 photons) \n + Za candidate selection"]
    #cuts =["lepton selection", "1", "passChaHadIso", "passChaHadIso&&passNeuHadIso", "passChaHadIso&&passNeuHadIso&&passHOverE", "passChaHadIso&&passNeuHadIso&&passHOverE&&passdR_gl&&(H_m>95&&H_m<180)"]
    
    
    #cuts_lable = ["2 photons", "passChaHadIso","passNeuHadIso","passHOverE","passdR_gl","H_m"]
    #cuts =["1", "passChaHadIso","passNeuHadIso","passHOverE","passdR_gl"]

    '''
    cuts_lable = ["2 leptons", "isolation", "lep tight", "m_{ll} > 50GeV"]
    cuts = cuts_lable
    '''

    #cuts =["passChaHadIso&&passNeuHadIso&&passHOverE&&passdR_gl&&(H_m>110&&H_m<180)"]

    for year in years:
        eff[year] = {}
        num[year] = {}
        #for channel in ['ele','mu']:
        for channel in ['ele']:
            eff[year][channel] = {}
            num[year][channel] = {}
            dem_final = {}
            for mass in mass_list:
                file_name = path_basic + year + '/' + 'ALP_M' + str(mass) + '.root'

                file = TFile(file_name)
                filesTree = file.Get("passedEvents")

                weight_value = file.Events_weight.GetBinContent(1)
                
                #print file_name
                dem_final[mass] = file.nEvents_weight_ntuple.GetBinContent(1) + file.nEvents_weight_ntuple.GetBinContent(2)
                #dem_final[mass] = file.nEvents_ntuple.GetEntries() * weight_value
                num[year][channel][mass] = {}


                cut = ""
                for i in range(len(cuts)):
                    eff[year][channel][cuts[i]] = array('d')
                    
                    if i == 0:
                        num[year][channel][mass][cuts[i]] = file.Z_50.GetEntries()*weight_value
                    else:
                        cut = cut + "&" + cuts[i]
                        cut = cut.lstrip("&")
                        num[year][channel][mass][cuts[i]] = getVariableHistsEventsNumber_weight(filesTree, "H_m", cut)

                    #print "year: " + str(year) + " mass:" + str(mass) + "dem: " + str(dem_final[mass])+", "+str(file.nEvents_weight_ntuple.GetBinContent(1) + file.nEvents_weight_ntuple.GetBinContent(2)) + " , num: " +str(num[year][channel][mass][cuts[i]])

            print "year:",year,", channel:", channel, ", dem:", dem_final

            for i in range(len(cuts)):
                for mass in mass_list:
                    eff[year][channel][cuts[i]].append(num[year][channel][mass][cuts[i]]/dem_final[mass])

            #print "{0} ALP mass:{1}, num_final:{2}, eff:{3}".format(year, mass, num_final, eff)

    #for channel in ['ele','mu']:
    for channel in ['ele']:
        for year in years:
            print year, eff[year][channel]
            canvas = TCanvas(channel+year,channel+year,650,600)
            canvas.cd()
            SetgStyle()
            canvas.SetRightMargin(0.1)
            canvas.SetLeftMargin(0.13)
            canvas.SetTopMargin(0.085)
            canvas.SetBottomMargin(0.12)
            colors = [2, 3, 4, 6]
            legend = TLegend(0.2,0.5,0.88,0.90)

            marks = [20, 23, 21, 22]

            gr = {}
            for cut in cuts:
                gr[cut] = TGraph( len(mass_list), m_x, eff[year][channel][cut] )
                gr[cut].SetLineColor( colors[cuts.index(cut)] )
                gr[cut].SetLineWidth( 2 )
                gr[cut].SetMarkerColor( colors[cuts.index(cut)] )
                gr[cut].SetMarkerStyle( marks[cuts.index(cut)] )
                #gr[year].SetMarkerSize(2)
                gr[cut].SetTitle('')
                gr[cut].GetHistogram().SetMaximum(0.8)
                gr[cut].GetHistogram().SetMinimum(0.0)
                gr[cut].GetXaxis().SetTitle( 'm_{a} (GeV)' )
                #gr[year].GetXaxis().SetLabelSize(0.045)
                #gr[year].GetXaxis().SetTitleSize(0.045)
                gr[cut].GetYaxis().SetTitle( 'Selection efficiency' )

                legend.AddEntry(gr[cut],cuts_lable[cuts.index(cut)],"ep")

                if cuts.index(cut) == 0:
                    gr[cut].Draw("ALP")
                else:
                    gr[cut].Draw("LP")
                
            legend.Draw('SAME')

            # CMS style
            #CMS_lumi.cmsText = "CMS"
            #CMS_lumi.extraText = "Simulation Preliminary"
            CMS_lumi.cmsText = ""
            CMS_lumi.extraText = ""
            CMS_lumi.cmsTextSize = 0.67
            CMS_lumi.lumiTextSize = 0.6
            CMS_lumi.lumiText_posY = 0.0075
            CMS_lumi.lumiText_posX = 0.0
            CMS_lumi.CMSText_posX = -0.03
            CMS_lumi.extraText_posX = 0.08
            CMS_lumi.outOfFrame = True
            CMS_lumi.CMS_lumi(canvas,5,0,'run2')

            canvas.SaveAs('./interpolation/cuteff_'+year+'.png')
            canvas.SaveAs('./interpolation/cuteff_'+year+'.pdf')
            canvas.Close()
        


    print '\n\n'

    print 'Done'


main()
