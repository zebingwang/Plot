####################################################
####################################################

import os
import sys
import numpy as np

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from Plot_Helper import LoadNtuples, MakeStack, CreateCanvas, DrawOnCanv, SaveCanvPic, MakeLumiLabel, MakeCMSDASLabel, ScaleSignal, MakeRatioPlot, MakeLegend
from Analyzer_Helper import MuPair, DimuCandidates, IsBtagLoose, IsBtagMed, IsBtagTight, CountYield
import Analyzer_Configs as AC
import Plot_Configs     as PC

from Analyzer_ALP import CountCutFlow, CountCutFlow_less, CountCutFlow_mva, CountCutFlow_mva_less, PIso2D, plot2D

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
args = parser.parse_args()


mva = False
if mva:
    name = mass + '_mva'
else:
    name = 'SFs'

if args.year == '2016':
    file_out = 'plots_16'
    out_name = "ALP_plot_data16_"+name+".root"
elif args.year == '2017':
    file_out = 'plots_17'
    out_name = "ALP_plot_data17_"+name+".root"
elif args.year == '2018':
    file_out = 'plots_18'
    out_name = "ALP_plot_data18_"+name+".root"
else:
    print "do not include at 2016/2017/2018"
    exit(0)

file_plot = file_out + "/plot_"+name

gROOT.SetBatch(True)
#gStyle.SetOptStat(1011111)

def main():

    #mass = 'M30'
    mass = 'massIndependent'

    analyzer_cfg = AC.Analyzer_Config('inclusive', mass, args.year)
    analyzer_cfg.Print_Config()

    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)

    CountCutFlow(analyzer_cfg, plot_cfg.lumi, file_out + '/')
    #CountCutFlow_less(analyzer_cfg, plot_cfg.lumi, 'plots/')
    #CountCutFlow_mva(analyzer_cfg, plot_cfg.lumi, 'plots/')
    #CountCutFlow_mva_less(analyzer_cfg, plot_cfg.lumi, 'plots/')

main()
