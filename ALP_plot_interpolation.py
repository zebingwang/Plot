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

from Analyzer_ALP import PIso2D, plot2D, plot2D_CONT

import CMS_lumi, tdrstyle

from xgboost import XGBClassifier
import pickle
import copy 
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-y", "--Year", dest="year", default="run2", help="which year's datasetes")
parser.add_argument('-C', '--CR', dest='CR', action='store_true', default=False, help='make control region')
parser.add_argument('-S', '--SR', dest='SR', action='store_true', default=False, help='make signal region')
parser.add_argument('-m', '--mva', dest='mva', action='store_true', default=False, help='use mva or not')
parser.add_argument('--cut', dest='cut', action='store_true', default=False, help='apply mva')
parser.add_argument("--mA", dest="mA", default="M5", help="ALP mass")
parser.add_argument('-b', '--blind', dest='blind', action='store_true', default=False, help='Blind signal region')
parser.add_argument('-ln', '--ln', dest='ln', action='store_true', default=False, help='log plot?')
args = parser.parse_args()



def main():

    mass_list = [1, 5, 15, 30]
    years = ['16', '17', '18']
    
    path_basic = '/publicfs/cms/user/wangzebing/ALP/Analysis_code/fit/fit_run2_param/'
    num_total = {1:97000.0, 5:97000.0, 15:95000.0, 30:94000.0}
    eff = {}

    for year in years:
        eff[year] = []

        for mass in mass_list:
            file_name = path_basic + 'ALP_data_sig_Am' + str(mass) + '_' + year + '_workspace.root'

            file = TFile(file_name)
            #print file_name
            num_final = file.CMS_hza_workspace.data("data_125_13TeV_cat0").numEntries()

            eff[year].append(num_final/num_total[mass])

            #print "{0} ALP mass:{1}, num_final:{2}, eff:{3}".format(year, mass, num_final, eff)

    colors = ['red','orange','blue']
    plt.xlim(xmax=35.0,xmin=0.0)
    #plt.ylim(ymax=1.0,ymin=0)
    plt.xlabel('m(a) GeV')
    plt.ylabel('$Efficiency\\times Acceptance$')
    area = np.pi * 4**2

    #print eff

    for year in years:
        plt.plot(mass_list, eff[year], "o-", c=colors[years.index(year)], label=year)
    #plt.plot(mass, eff, s=area)
    plt.grid()
    plt.legend()
    plt.savefig('test.png')
    plt.close('all')


    print '\n\n'

    print 'Done'


main()
