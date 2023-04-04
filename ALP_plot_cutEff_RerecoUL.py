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


def getVariableHistsEventsNumber_weight(Tree,varName,cut):

    #Make a canvas to do the work on
    canvas = TCanvas(varName,varName,1000,800)

    #Extract the relevant variable from the trees
    Tree.Draw("{0}>>tree{0}".format(varName),"factor*pho1SFs*pho2SFs*({0})".format(cut))
    Hist = gDirectory.Get("tree{0}".format(varName))

    canvas.Clear()

    return Hist.Integral()


def main():

    mass_list = [1,5,15,30]
    
    #years = ['16', '16APV', '17', '18']
    years = ['UL/17', 'Rereco/17/massInde']
    
    path_basic = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/'
    eff = {}
    num = {}

    #cuts_lable = ["HLT", "lep_selection", "2 photons", "pho_selection"]
    #cuts =["HLT", "lepton selection", "1", "passChaHadIso&&passNeuHadIso&&passHOverE&&passdR_gl&&(H_m>110&&H_m<180)"]
    cuts_lable = ["2 photons", "passChaHadIso","passNeuHadIso","passHOverE","passdR_gl"]
    cuts =["1", "passChaHadIso","passNeuHadIso","passHOverE","passdR_gl"]

    '''
    cuts_lable = ["2 leptons", "isolation", "lep tight", "m_{ll} > 50GeV"]
    cuts = cuts_lable
    '''

    #cuts =["passChaHadIso&&passNeuHadIso&&passHOverE&&passdR_gl&&(H_m>110&&H_m<180)"]

    for year in years:
        eff[year] = {}
        num[year] = {}
        for channel in ['ele','mu']:
            eff[year][channel] = {}
            num[year][channel] = {}
            dem_final = {}
            for mass in mass_list:
                file_name = path_basic + year + '/' + 'ALP_M' + str(mass) + '.root'

                file = TFile(file_name)
                filesTree = file.Get("passedEvents")
                
                
                weight_value = file.Events_weight.GetBinContent(1)
                
                
                if 'Rereco' in year:
                    dem_final[mass] = file.nEvents_ntuple.GetEntries() * weight_value *1.5
                else:
                    dem_final[mass] = file.nEvents_ntuple.GetEntries() * weight_value 
                    #dem_final[mass] = file.nEvents_weight_ntuple.GetBinContent(1) + file.nEvents_weight_ntuple.GetBinContent(2)
                
                num[year][channel][mass] = {}


                cut = ""
                for i in range(len(cuts)):
                    eff[year][channel][cuts[i]] = []
                    
                    if i == 0:
                        num[year][channel][mass][cuts[i]] = file.nEvents_trig.GetEntries()*weight_value
                    elif i == 1:
                        num[year][channel][mass][cuts[i]] = file.Z_50.GetEntries()*weight_value
                    else:
                        cut = cut + "&" + cuts[i]
                        cut = cut.lstrip("&")
                        num[year][channel][mass][cuts[i]] = getVariableHistsEventsNumber_weight(filesTree, "H_m", cut)
                    
                    '''
                    if i ==0:
                        num[year][channel][mass][cuts[i]] = (file.Z_e_nocut.GetEntries() + file.Z_mu_nocut.GetEntries())*weight_value
                    elif i == 1:
                        num[year][channel][mass][cuts[i]] = (file.Z_e_lIso.GetEntries() + file.Z_mu_lIso.GetEntries())*weight_value
                    elif i == 2:
                        num[year][channel][mass][cuts[i]] = (file.Z_e_lIso_lTight.GetEntries() + file.Z_mu_lIso_lTight.GetEntries())*weight_value
                    else:
                        num[year][channel][mass][cuts[i]] = file.Z_50.GetEntries()*weight_value
                    '''

                #print "year: " + str(year) + " mass:" + str(mass) + "dem: " + str(file.nEvents_ntuple.GetEntries() * weight_value)+", "+str(file.nEvents_weight_ntuple.GetBinContent(1) + file.nEvents_weight_ntuple.GetBinContent(2)) + " , num: " +str(num[year][channel][mass][cuts[0]])

            for i in range(len(cuts)):
                for mass in mass_list:
                    eff[year][channel][cuts[i]].append(num[year][channel][mass][cuts[i]]/dem_final[mass])

            #print "{0} ALP mass:{1}, num_final:{2}, eff:{3}".format(year, mass, num_final, eff)

    for channel in ['ele','mu']:
        #colors = ['red','orange','blue', 'green']
        colors = ['orange','green']
        shapes = ['o-', '^-', 'x-', 's-', '>-', 'd-', '<-']
        plt.xlim(xmax=40.0,xmin=0.0)
        #plt.ylim(ymax=1.0,ymin=0)
        #plt.yscale('symlog')
        #plt.yscale('log')
        plt.xlabel('m(a) GeV')
        plt.ylabel('$Efficiency\\times Acceptance$')


        for year in years:
            for cut in cuts:
                plt.plot(mass_list, eff[year][channel][cut], shapes[cuts.index(cut)], c=colors[years.index(year)], label=year+'_'+cuts_lable[cuts.index(cut)])
        
        plt.grid()
        plt.legend(fontsize=8)
        plt.savefig('./interpolation/cuteff_RerecoUL_'+channel+'.png')
        plt.close('all')
        


    print '\n\n'

    print 'Done'


main()
