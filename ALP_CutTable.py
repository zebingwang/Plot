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

from Analyzer_ALP import CountCutFlow, CountCutFlow_less, CountCutFlow_mva, CountCutFlow_mva_less, PIso2D

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
args = parser.parse_args()

file_out = 'plots_run2UL'

gROOT.SetBatch(True)
#gStyle.SetOptStat(1011111)

def main():

    #mass = 'M30'
    #mass = 'massIndependent'
    version = 'UL'

    analyzer_cfg = AC.Analyzer_Config('inclusive', args.year,"","")
    analyzer_cfg.Print_Config()

    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)

    CountCutFlow(analyzer_cfg, plot_cfg.lumi, file_out + '/', version)

main()
