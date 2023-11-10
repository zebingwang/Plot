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

   
def GetSys(list_eff, list_effSys, sys_names, mass_list):
    
    sys_up = {}
    sys_dn = {}

    sys_total_up = array('d')
    sys_total_dn = array('d')
    for m in range(len(mass_list)):
        norm = list_eff[m]
        sys_up[m] = array('d')
        sys_dn[m] = array('d')

        for s in sys_names:
            if 'up' in s:
                if list_effSys[s][m]-norm >0.0:
                    sys_up[m].append(list_effSys[s][m]-norm)
                else:
                    sys_dn[m].append(norm-list_effSys[s][m])
            else:
                if norm - list_effSys[s][m] >0.0:
                    sys_dn[m].append(norm-list_effSys[s][m])
                else:
                    sys_up[m].append(list_effSys[s][m]-norm)

        sys_total_up.append(np.sqrt(sys_up[m][0]*sys_up[m][0]+sys_up[m][1]*sys_up[m][1]+sys_up[m][2]*sys_up[m][2]))
        sys_total_dn.append(np.sqrt(sys_dn[m][0]*sys_dn[m][0]+sys_dn[m][1]*sys_dn[m][1]+sys_dn[m][2]*sys_dn[m][2]))

    return sys_total_up, sys_total_dn
    
def main():

    mass_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
    #mass_list = [1]
    #years = ['18']
    m_x = array('d')
    for m in mass_list:
        m_x.append(m)
    mass_interp = range(1,31)
    years = ['16', '16APV', '17', '18']
    sys_names = ['CMS_eff_g_up','CMS_eff_g_dn','CMS_pileup_up','CMS_pileup_dn','CMS_eff_lep_up','CMS_eff_lep_dn']
    
    #path_basic_num = '/publicfs/cms/user/wangzebing/ALP/Analysis_code/fit/fit_run2_UL_eff/'
    path_basic_num = '/publicfs/cms/user/wangzebing/ALP/Analysis_code/fit/fit_run2_UL_eff_beforBDT/'
    path_basic_dem = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/'
    eff = {}
    eff_sys = {}
    stat_dn = {}
    stat_up = {}
    eff_sys_total_up = {}
    eff_sys_total_dn = {}

    print "start"
    for year in years:
        eff[year] = {}
        eff_sys[year] = {}
        stat_dn[year] = {}
        stat_up[year] = {}
        eff_sys_total_up[year] = {}
        eff_sys_total_dn[year] = {}
        for channel in ['ele','mu']:
            eff[year][channel] = array('d')
            eff_sys[year][channel] = {}
            stat_dn[year][channel] = array('d')
            stat_up[year][channel] = array('d')
            eff_sys_total_up[year][channel] = array('d')
            eff_sys_total_dn[year][channel] = array('d')
            for s in sys_names:
                eff_sys[year][channel][s] = array('d')
            for mass in mass_list:
                file_name_num = path_basic_num + 'ALP_data_sig_Am' + str(mass) + '_Hm125_' + year + '_workspace_'+channel+'.root'
                file_name_dem = path_basic_dem + year + '/' + 'ALP_M' + str(mass) + '.root'

                file_num = TFile(file_name_num)
                file_dem = TFile(file_name_dem)
                
                #print file_name
                signal_scale = 1.0
                num_final_sys = {}
                num_final = file_num.CMS_hza_workspace.data("ggh_125_13TeV_cat0").sumEntries() * signal_scale
                dem_final = ( file_dem.nEvents_weight_ntuple.GetBinContent(1) + file_dem.nEvents_weight_ntuple.GetBinContent(2) ) * signal_scale
                num_final_entry = file_num.CMS_hza_workspace.data("ggh_125_13TeV_cat0_norm").sumEntries() * signal_scale
                dem_final_entry = ( file_dem.nEvents_ntuple.GetBinContent(1) + file_dem.nEvents_ntuple.GetBinContent(2) ) * signal_scale

                #print 'year: ', year, 'channel: ', channel, 'mass: ', mass, 'dem: ', dem_final, 'num: ', num_final, 'eff: ', num_final/dem_final*100., 'nEvents_entry_num: ', num_final_entry, 'nEvents_entry_dem', dem_final_entry

                eff[year][channel].append(num_final/dem_final*100.)
                #file_num.CMS_hza_workspace.Print()

                #stat_dn[year][channel].append(num_final/dem_final*100.0 - TEfficiency.ClopperPearson(dem_final, num_final, 0.683,  0)*100.0)
                #stat_up[year][channel].append(TEfficiency.ClopperPearson(dem_final, num_final, 0.683,  1)*100.0 - num_final/dem_final*100.0)

                stat_dn[year][channel].append(num_final/dem_final*100.0 - TEfficiency.ClopperPearson(dem_final_entry, num_final_entry, 0.683,  0)*100.0)
                stat_up[year][channel].append(TEfficiency.ClopperPearson(dem_final_entry, num_final_entry, 0.683,  1)*100.0 - num_final/dem_final*100.0)

                for s in sys_names:
                    #print "[[Debug]]", "ggh_125_13TeV_cat0_{0}".format(s)
                    num_final_sys[s] = file_num.CMS_hza_workspace.data("ggh_125_13TeV_cat0_{0}".format(s)).sumEntries()*signal_scale

                    eff_sys[year][channel][s].append(num_final_sys[s]/dem_final*100.)
            
            sys_up, sys_dn = GetSys(eff[year][channel], eff_sys[year][channel], sys_names, mass_list)

            for m in range(len(stat_up[year][channel])):

                eff_sys_total_up[year][channel].append(np.sqrt(stat_up[year][channel][m]*stat_up[year][channel][m]+sys_up[m]*sys_up[m]))
                eff_sys_total_dn[year][channel].append(np.sqrt(stat_dn[year][channel][m]*stat_dn[year][channel][m]+sys_dn[m]*sys_dn[m]))
            
            #print 'stat: ', stat_up[year][channel], 'sys: ', sys_up

            #for  s in sys_names:
            #    print "{0}, {1}, sys:{2}, eff:{3}, nom_eff:{4}".format(year, channel, s, eff_sys[year][channel][s], eff[year][channel])
            #print "{0} ALP mass:{1}, num_final:{2}, eff:{3}".format(year, mass, num_final, eff)

    

    for channel in ['ele','mu']:
        
        canvas = TCanvas(channel,channel,650,600)
        canvas.cd()
        SetgStyle()
        canvas.SetRightMargin(0.1)
        canvas.SetLeftMargin(0.13)
        canvas.SetTopMargin(0.085)
        canvas.SetBottomMargin(0.12)
        colors = [2, 6, 209, 4]
        legend = TLegend(0.6,0.65,0.88,0.88)
        #gStyle.SetLegendFont(42)
        #gStyle.SetLegendTextSize(0.03)
        #legend.SetBorderSize(0)
        

        gr = {}
        marks = [20, 23, 21, 22]

        eff_16_mean = array('d')
        err_x = array('d')
        for i in range(len(mass_list)):
            eff_16_mean.append((eff['16'][channel][i] + eff['16APV'][channel][i])/2.0)
            err_x.append(0.0)

        for year in years:

            if '16APV' == year: continue

            if '16' in year:
                gr[year] = TGraphAsymmErrors( len(mass_list), m_x, eff_16_mean, err_x, err_x, eff_sys_total_dn[year][channel], eff_sys_total_up[year][channel])
            else:
                gr[year] = TGraphAsymmErrors( len(mass_list), m_x, eff[year][channel], err_x, err_x, eff_sys_total_dn[year][channel], eff_sys_total_up[year][channel])
            gr[year].SetLineColor( colors[years.index(year)] )
            gr[year].SetLineWidth( 2 )
            gr[year].SetMarkerColor( colors[years.index(year)] )
            gr[year].SetMarkerStyle( marks[years.index(year)] )
            #gr[year].SetMarkerSize(2)
            gr[year].SetTitle('')
            gr[year].GetHistogram().SetMaximum(4.0)
            gr[year].GetHistogram().SetMinimum(0.2)
            gr[year].GetXaxis().SetTitle( 'm_{a} (GeV)' )
            #gr[year].GetXaxis().SetLabelSize(0.045)
            #gr[year].GetXaxis().SetTitleSize(0.045)
            gr[year].GetYaxis().SetTitle( 'Efficiency #times Acceptance (%)' )

            legend.AddEntry(gr[year],"20"+year,"ep")
            
            gr[year].SetFillColor(colors[years.index(year)])
            gr[year].SetFillStyle(3001)

            
        
        gr['16'].Draw("ALP")
        gr['17'].Draw("LP")
        gr['18'].Draw("LP")

        # CMS style
        CMS_lumi.cmsText = "CMS"
        CMS_lumi.extraText = "Simulation Supplementary"
        #CMS_lumi.extraText = "Simulation"
        #CMS_lumi.cmsText = ""
        #CMS_lumi.extraText = ""
        CMS_lumi.cmsTextSize = 0.67
        CMS_lumi.lumiTextSize = 0.6
        CMS_lumi.lumiText_posY = 0.0075
        CMS_lumi.lumiText_posX = 0.0
        CMS_lumi.CMSText_posX = -0.03
        CMS_lumi.extraText_posX = 0.11
        CMS_lumi.outOfFrame = True
        CMS_lumi.CMS_lumi(canvas,5,0,'run2')

        latex = TLatex()
        latex.SetNDC()
        latex.SetTextColor(kBlack)
        latex.SetTextFont(42)
        latex.SetTextSize(0.045)
        latex.SetTextAlign(31)
        if channel=='ele':
            latex.DrawLatex(0.5,0.85, r"H \rightarrow Za \rightarrow ee + 2\gamma")
        else:
            latex.DrawLatex(0.5,0.85, r"H \rightarrow Za \rightarrow \mu\mu + 2\gamma")

        legend.Draw('SAME')

        canvas.SaveAs('./interpolation/eff_'+channel+'_beforeBDT.png')
        canvas.SaveAs('./interpolation/eff_'+channel+'_beforeBDT.pdf')
        canvas.Close()

        for m in range(len(mass_list)):
            print "{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}".format(mass_list[m],eff_16_mean[m],eff_sys_total_dn['16'][channel][m],eff['17'][channel][m],eff_sys_total_dn['17'][channel][m],eff['18'][channel][m],eff_sys_total_dn['18'][channel][m])
        
      
    

    print '\n\n'

    print 'Done'


main()
