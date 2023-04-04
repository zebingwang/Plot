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

from Analyzer_ALP import PIso2D, plot2D_mgg, plot2D_mllgg, plot2D_CONT

import CMS_lumi, tdrstyle

from xgboost import XGBClassifier
import pickle
import copy 

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
parser.add_argument('--ele', dest='ele', action='store_true', default=False, help='electron channel?')
parser.add_argument('--mu', dest='mu', action='store_true', default=False, help='muon channel?')
args = parser.parse_args()

version = 'UL'



gROOT.SetBatch(True)
#tdrstyle.setTDRStyle()

mva = args.mva
if args.CR:
    name = version + '_CR'
elif args.SR:
    name = version + '_SR'
elif mva:
    name = version + '_mva'
else:
    name = version

if args.cut: name = name + '_cut_' + args.mA

if args.ele:
    name = name + '_ele'
if args.mu:
    name = name + '_mu'


if args.year == '2016':
    file_out = 'plots_16'
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2016.pkl"
    #mvaCut = 0.8675
elif args.year == '2017':
    file_out = 'plots_17'
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2017.pkl"
    #mvaCut = 0.8365
elif args.year == '2018':
    file_out = 'plots_18'
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2018.pkl"
    #mvaCut = 0.9766
elif args.year == 'run2':
    file_out = 'plots_run2UL'

    if args.ele:
        BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/UL/model_ALP_BDT_param_ele.pkl"
    elif args.mu:
        BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/UL/model_ALP_BDT_param_mu.pkl"
    else:
        BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/UL/model_ALP_BDT_param.pkl"

else:
    print "do not include at 2016/2017/2018"
    exit(0)

file_plot = file_out + "/plot_"+name
# load the model from disk
model = pickle.load(open(BDT_filename, 'rb'))

def main():

    if not os.path.exists(file_out):
        os.makedirs(file_out)

    if not os.path.exists(file_plot):
        os.makedirs(file_plot)


    analyzer_cfg = AC.Analyzer_Config('inclusive', args.year)

    analyzer_cfg.Print_Config()
    ntuples = LoadNtuples(analyzer_cfg)


    #### systematic uncertainties ######
    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)

    ### declare histograms

    histos = {}
    histos_sys = {}

    mass_list = {'M1':1.0, 'M2':2.0, 'M3':3.0, 'M4':4.0, 'M5':5.0, 'M6':6.0, 'M7':7.0, 'M8':8.0, 'M9':9.0, 'M10':10.0, 'M15':15.0, 'M20':20.0, 'M25':25.0, 'M30':30.0}
    histos['mvaVal_bkg']    = TH2F('mvaVal_bkg', 'mvaVal_bkg', 35,  0., 35., 100, 0., 1)#2D

    histos['mvaVal_bkg_mgg'] = {}
    histos['mvaVal_bkg_mllgg'] = {}
    for ALP_mass in mass_list:
        histos['mvaVal_bkg_mgg'][ALP_mass] = TH2F('mvaVal_bkg_mgg_'+ALP_mass, 'mvaVal_bkg_mgg_'+ALP_mass, 35,  0., 35., 100, 0., 1)#2D
        histos['mvaVal_bkg_mllgg'][ALP_mass] = TH2F('mvaVal_bkg_mllgg'+ALP_mass, 'mvaVal_bkg_mllgg'+ALP_mass, 100,  115., 180., 100, 0., 1)#2D

   
    ### loop over samples and events
    
    for sample in analyzer_cfg.bkg_names:
        ntup = ntuples[sample] # just a short name
        print '\n\nOn sample: %s' %sample
        print 'total events: %d' %ntup.GetEntries()

        for iEvt in range( ntup.GetEntries() ):
            ntup.GetEvent(iEvt)
            #if (iEvt == 50): break


            if (iEvt % 100000 == 1):
                print "looking at event %d" %iEvt


            if args.ele:
                if abs(ntup.l1_id) == 13: 
                    continue
            if args.mu:
                if abs(ntup.l1_id) == 11: 
                    continue


            weight = ntup.factor * ntup.pho1SFs * ntup.pho2SFs

            if (ntup.H_m > -90):
                if not ntup.passChaHadIso: continue
                if not ntup.passNeuHadIso: continue
                if not ntup.passdR_gl: continue
                if not ntup.passHOverE: continue
                if ntup.H_m>180. or ntup.H_m<110.: continue

                if  args.CR:
                    if ntup.H_m<135. and ntup.H_m>115.: continue
                elif  args.SR:
                    if ntup.H_m>135. or ntup.H_m<115.: continue


                MVA_value = {}
                if mva:
                    
                    for ALP_mass in mass_list:

                        param = (ntup.ALP_m - mass_list[ALP_mass])/ntup.H_m

                        MVA_list = [ntup.pho1Pt, ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr ,ntup.pho2Pt, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.pho2PIso_noCorr,ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, ntup.var_PtaOverMh, ntup.H_pt, param]
                        MVA_value[ALP_mass] = model.predict_proba(MVA_list)[:, 1]


                    if args.cut:
                        if MVA_value[args.mA] < mvaCut[args.mA]: continue

                    for ALP_mass in mass_list:
                        histos['mvaVal_bkg'].Fill(mass_list[ALP_mass], MVA_value[ALP_mass], weight)
                        histos['mvaVal_bkg_mgg'][ALP_mass].Fill(ntup.ALP_m, MVA_value[ALP_mass], weight)
                        histos['mvaVal_bkg_mllgg'][ALP_mass].Fill(ntup.H_m, MVA_value[ALP_mass], weight)


        ## End of for iEvt in range( ntup.GetEntries() )
    ## End of for sample in analyzer_cfg.samp_names

    gStyle.SetPalette(kGreyScale)
    TColor.InvertPalette()
    
    canv2D = CreateCanvas("mvaVal_bkg")
    plot2D_mgg(canv2D, histos["mvaVal_bkg"])
    SaveCanvPic(canv2D, file_plot, "mvaVal_bkg")

    canv2D_mgg = {}
    canv2D_mllgg = {}
    for ALP_mass in mass_list:
        canv2D_mgg[ALP_mass] = CreateCanvas("mvaVal_bkg_mgg_"+ALP_mass)
        plot2D_mgg(canv2D_mgg[ALP_mass], histos['mvaVal_bkg_mgg'][ALP_mass])
        SaveCanvPic(canv2D_mgg[ALP_mass], file_plot, "mvaVal_bkg_mgg_"+ALP_mass)

        canv2D_mllgg[ALP_mass] = CreateCanvas("mvaVal_bkg_mllgg_"+ALP_mass)
        plot2D_mllgg(canv2D_mllgg[ALP_mass], histos['mvaVal_bkg_mllgg'][ALP_mass])
        SaveCanvPic(canv2D_mllgg[ALP_mass], file_plot, "mvaVal_bkg_mllgg_"+ALP_mass)

    print histos["mvaVal_bkg"].GetCorrelationFactor()
    print '\n\n'
    #CountYield(analyzer_cfg, histos['ALP_m'])

    print 'Done'


main()
