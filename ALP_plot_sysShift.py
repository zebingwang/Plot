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

from Analyzer_ALP import PIso2D, plot2D, plot2D_CONT

import CMS_lumi, tdrstyle

from xgboost import XGBClassifier
import pickle

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
parser.add_argument('-C', '--CR', dest='CR', action='store_true', default=False, help='make control region')
parser.add_argument("-s", "--sys", dest="sys", default="0.0", help="systematic shift")
args = parser.parse_args()

mass = 'massIndependent'



gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

mva = True

name = 'sysShift'

if args.year == '2016':
    file_out = 'plots_16'
    out_name = "ALP_plot_data16_"+name+".root"
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_massindependent_2016.pkl"
    mvaCut = 0.8675
elif args.year == '2017':
    file_out = 'plots_17'
    out_name = "ALP_plot_data17_"+name+".root"
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_massindependent_2017.pkl"
    mvaCut = 0.8571
elif args.year == '2018':
    file_out = 'plots_18'
    out_name = "ALP_plot_data18_"+name+".root"
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_massindependent_2018.pkl"
    mvaCut = 0.9204
else:
    print "do not include at 2016/2017/2018"
    exit(0)

file_plot = file_out + "/plot_"+name
# load the model from disk
model = pickle.load(open(BDT_filename, 'rb'))

def main():

    if not os.path.exists(file_out):
        os.makedirs(file_out)

    out_file = TFile( file_out + '/' + out_name , "RECREATE")

    analyzer_cfg = AC.Analyzer_Config('inclusive', mass,args.year)

    ############ CR plots#########
    #if args.CR:
        #analyzer_cfg.sample_loc = analyzer_cfg.sample_loc.replace('massInde','massInde/CR')
        #analyzer_cfg.sig_names  = ['']
    ##############################
    analyzer_cfg.Print_Config()
    ntuples = LoadNtuples(analyzer_cfg)


    #var_names = ['Z_m', 'H_m', 'ALP_m', 'pho1IetaIeta', 'pho1IetaIeta55', 'pho1PIso', 'var_dR_Za', 'var_dR_g1g2', 'var_dR_g1Z', 'var_PtaOverMa', 'var_PtaOverMh', 'var_Pta', 'var_MhMa', 'var_MhMZ', 'ALP_calculatedPhotonIso']
    var_names = ['Z_m', 'H_m', 'ALP_m','pho1Pt', 'pho1eta', 'pho1phi', 'pho1R9', 'pho1IetaIeta', 'pho1IetaIeta55','pho1PIso_noCorr' ,'pho2Pt', 'pho2eta', 'pho2phi', 'pho2R9', 'pho2IetaIeta', 'pho2IetaIeta55','pho2PIso_noCorr','ALP_calculatedPhotonIso', 'var_dR_Za', 'var_dR_g1g2', 'var_dR_g1Z', 'var_PtaOverMh', 'var_Pta', 'var_MhMZ', 'H_pt', 'var_PtaOverMa', 'var_MhMa']
    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)

    ### declare histograms

    histos_up = {}
    histos_down = {}

    for var_name in var_names:
        histos_up[var_name] = {}
        histos_down[var_name] = {}
    for sample in analyzer_cfg.bkg_names:
        histos_up['pho1Pt'][sample]    = TH1F('UP_pho1Pt'    + '_' + sample, 'pho1Pt'    + '_' + sample, 50,  5., 50.)
        histos_up['pho1eta'][sample]    = TH1F('UP_pho1eta'    + '_' + sample, 'pho1eta'    + '_' + sample, 100,  -3., 3.)
        histos_up['pho1phi'][sample]    = TH1F('UP_pho1phi'    + '_' + sample, 'pho1phi'    + '_' + sample, 100,  -4., 4.)
        histos_up['pho1R9'][sample]    = TH1F('UP_pho1R9'    + '_' + sample, 'pho1R9'    + '_' + sample, 100,  0., 1.)
        histos_up['pho1IetaIeta'][sample]    = TH1F('UP_pho1IetaIeta'    + '_' + sample, 'pho1IetaIeta'    + '_' + sample, 100,  0., 0.07)
        histos_up['pho1IetaIeta55'][sample]    = TH1F('UP_pho1IetaIeta55'    + '_' + sample, 'pho1IetaIeta55'    + '_' + sample, 100,  0., 0.07)
        histos_up['pho1PIso_noCorr'][sample]    = TH1F('UP_pho1PIso_noCorr'    + '_' + sample, 'pho1PIso_noCorr'    + '_' + sample, 100, 0., 40.)
        histos_up['pho2Pt'][sample]    = TH1F('UP_pho2Pt'    + '_' + sample, 'pho2Pt'    + '_' + sample, 50,  5., 30.)
        histos_up['pho2eta'][sample]    = TH1F('UP_pho2eta'    + '_' + sample, 'pho2eta'    + '_' + sample, 100,  -3., 3.)
        histos_up['pho2phi'][sample]    = TH1F('UP_pho2phi'    + '_' + sample, 'pho2phi'    + '_' + sample, 100,  -4., 4.)
        histos_up['pho2R9'][sample]    = TH1F('UP_pho2R9'    + '_' + sample, 'pho2R9'    + '_' + sample, 100,  0., 1.)
        histos_up['pho2IetaIeta'][sample]    = TH1F('UP_pho2IetaIeta'    + '_' + sample, 'pho2IetaIeta'    + '_' + sample, 100,  0., 0.07)
        histos_up['pho2IetaIeta55'][sample]    = TH1F('UP_pho2IetaIeta55'    + '_' + sample, 'pho2IetaIeta55'    + '_' + sample, 100,  0., 0.07)
        histos_up['pho2PIso_noCorr'][sample]    = TH1F('UP_pho2PIso_noCorr'    + '_' + sample, 'pho2PIso_noCorr'    + '_' + sample, 100, 0., 40.)
        histos_up['Z_m'][sample]    = TH1F('UP_Z_m'    + '_' + sample, 'Z_m'    + '_' + sample, 100,  50., 130.)
        histos_up['H_m'][sample]    = TH1F('UP_H_m'    + '_' + sample, 'H_m'    + '_' + sample, 20,  118., 140.)
        histos_up['H_pt'][sample]    = TH1F('UP_H_pt'    + '_' + sample, 'H_pt'    + '_' + sample, 100,  0., 200.)
        histos_up['ALP_m'][sample] = TH1F('UP_ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 50, 0., 40.)
        histos_up['var_dR_g1g2'][sample] = TH1F('UP_var_dR_g1g2' + '_' + sample, 'var_dR_g1g2' + '_' + sample, 50, 0., 5)
        histos_up['var_PtaOverMa'][sample] = TH1F('UP_var_PtaOverMa' + '_' + sample, 'var_PtaOverMa' + '_' + sample, 100, 0., 100.)
        histos_up['var_dR_Za'][sample] = TH1F('UP_var_dR_Za' + '_' + sample, 'var_dR_Za' + '_' + sample, 50, 0., 7.)
        histos_up['var_dR_g1Z'][sample] = TH1F('UP_var_dR_g1Z' + '_' + sample, 'var_dR_g1Z' + '_' + sample, 50, 0., 8.)
        histos_up['var_PtaOverMh'][sample] = TH1F('UP_var_PtaOverMh' + '_' + sample, 'var_PtaOverMh' + '_' + sample, 100, 0., 1.)
        histos_up['var_Pta'][sample] = TH1F('UP_var_Pta' + '_' + sample, 'var_Pta' + '_' + sample, 50, 0., 60.)
        histos_up['var_MhMa'][sample] = TH1F('UP_var_MhMa' + '_' + sample, 'var_MhMa' + '_' + sample, 50, 100., 200.)
        histos_up['var_MhMZ'][sample] = TH1F('UP_var_MhMZ' + '_' + sample, 'var_MhMZ' + '_' + sample, 50, 120., 250.)
        histos_up['ALP_calculatedPhotonIso'][sample] = TH1F('UP_ALP_calculatedPhotonIso' + '_' + sample, 'ALP_calculatedPhotonIso' + '_' + sample, 100, 0., 5.)

        histos_down['pho1Pt'][sample]    = TH1F('Down_pho1Pt'    + '_' + sample, 'pho1Pt'    + '_' + sample, 50,  5., 50.)
        histos_down['pho1eta'][sample]    = TH1F('Down_pho1eta'    + '_' + sample, 'pho1eta'    + '_' + sample, 100,  -3., 3.)
        histos_down['pho1phi'][sample]    = TH1F('Down_pho1phi'    + '_' + sample, 'pho1phi'    + '_' + sample, 100,  -4., 4.)
        histos_down['pho1R9'][sample]    = TH1F('Down_pho1R9'    + '_' + sample, 'pho1R9'    + '_' + sample, 100,  0., 1.)
        histos_down['pho1IetaIeta'][sample]    = TH1F('Down_pho1IetaIeta'    + '_' + sample, 'pho1IetaIeta'    + '_' + sample, 100,  0., 0.07)
        histos_down['pho1IetaIeta55'][sample]    = TH1F('Down_pho1IetaIeta55'    + '_' + sample, 'pho1IetaIeta55'    + '_' + sample, 100,  0., 0.07)
        histos_down['pho1PIso_noCorr'][sample]    = TH1F('Down_pho1PIso_noCorr'    + '_' + sample, 'pho1PIso_noCorr'    + '_' + sample, 100, 0., 40.)
        histos_down['pho2Pt'][sample]    = TH1F('Down_pho2Pt'    + '_' + sample, 'pho2Pt'    + '_' + sample, 50,  5., 30.)
        histos_down['pho2eta'][sample]    = TH1F('Down_pho2eta'    + '_' + sample, 'pho2eta'    + '_' + sample, 100,  -3., 3.)
        histos_down['pho2phi'][sample]    = TH1F('Down_pho2phi'    + '_' + sample, 'pho2phi'    + '_' + sample, 100,  -4., 4.)
        histos_down['pho2R9'][sample]    = TH1F('Down_pho2R9'    + '_' + sample, 'pho2R9'    + '_' + sample, 100,  0., 1.)
        histos_down['pho2IetaIeta'][sample]    = TH1F('Down_pho2IetaIeta'    + '_' + sample, 'pho2IetaIeta'    + '_' + sample, 100,  0., 0.07)
        histos_down['pho2IetaIeta55'][sample]    = TH1F('Down_pho2IetaIeta55'    + '_' + sample, 'pho2IetaIeta55'    + '_' + sample, 100,  0., 0.07)
        histos_down['pho2PIso_noCorr'][sample]    = TH1F('Down_pho2PIso_noCorr'    + '_' + sample, 'pho2PIso_noCorr'    + '_' + sample, 100, 0., 40.)
        histos_down['Z_m'][sample]    = TH1F('Down_Z_m'    + '_' + sample, 'Z_m'    + '_' + sample, 100,  50., 130.)
        histos_down['H_m'][sample]    = TH1F('Down_H_m'    + '_' + sample, 'H_m'    + '_' + sample, 20,  118., 140.)
        histos_down['H_pt'][sample]    = TH1F('Down_H_pt'    + '_' + sample, 'H_pt'    + '_' + sample, 100,  0., 200.)
        histos_down['ALP_m'][sample] = TH1F('Down_ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 50, 0., 40.)
        histos_down['var_dR_g1g2'][sample] = TH1F('Down_var_dR_g1g2' + '_' + sample, 'var_dR_g1g2' + '_' + sample, 50, 0., 5)
        histos_down['var_PtaOverMa'][sample] = TH1F('Down_var_PtaOverMa' + '_' + sample, 'var_PtaOverMa' + '_' + sample, 100, 0., 100.)
        histos_down['var_dR_Za'][sample] = TH1F('Down_var_dR_Za' + '_' + sample, 'var_dR_Za' + '_' + sample, 50, 0., 7.)
        histos_down['var_dR_g1Z'][sample] = TH1F('Down_var_dR_g1Z' + '_' + sample, 'var_dR_g1Z' + '_' + sample, 50, 0., 8.)
        histos_down['var_PtaOverMh'][sample] = TH1F('Down_var_PtaOverMh' + '_' + sample, 'var_PtaOverMh' + '_' + sample, 100, 0., 1.)
        histos_down['var_Pta'][sample] = TH1F('Down_var_Pta' + '_' + sample, 'var_Pta' + '_' + sample, 50, 0., 60.)
        histos_down['var_MhMa'][sample] = TH1F('Down_var_MhMa' + '_' + sample, 'var_MhMa' + '_' + sample, 50, 100., 200.)
        histos_down['var_MhMZ'][sample] = TH1F('Down_var_MhMZ' + '_' + sample, 'var_MhMZ' + '_' + sample, 50, 120., 250.)
        histos_down['ALP_calculatedPhotonIso'][sample] = TH1F('Down_ALP_calculatedPhotonIso' + '_' + sample, 'ALP_calculatedPhotonIso' + '_' + sample, 100, 0., 5.)



    ### loop over samples and events
    for sample in analyzer_cfg.bkg_names:
        ntup = ntuples[sample] # just a short name
        print '\n\nOn sample: %s' %sample
        print 'total events: %d' %ntup.GetEntries()

        for iEvt in range( ntup.GetEntries() ):
            ntup.GetEvent(iEvt)
            #if (iEvt == 50000): break


            if (iEvt % 100000 == 1):
                print "looking at event %d" %iEvt


            weight = ntup.factor

            if (ntup.H_m > -90):

                #if ntup.dR_pho < 0.02: continue
                #if not ntup.passEleVeto: continue
                if not ntup.passChaHadIso: continue
                if not ntup.passNeuHadIso: continue
                if not ntup.passdR_gl: continue
                if not ntup.passHOverE: continue
                #if ntup.H_m>130. or ntup.H_m<118.: continue
                if ntup.H_m<130. and ntup.H_m>118.: continue

                #if abs(ntup.l1_id)!=13: continue

                #print 'pho1: pt:' + str(ntup.pho1Pt) + ', eta: ' + str(ntup.pho1eta) + ', SFs: ' + str(pho1_SFs)
                if mva:
                    #MVA_list = [ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_PtaOverMh, ntup.var_MhMZ, ntup.var_dR_g1g2]
                    MVA_list = [ntup.pho1Pt, ntup.pho1eta, ntup.pho1phi, ntup.pho1R9, ntup.pho1IetaIeta55 ,ntup.pho2Pt, ntup.pho2eta, ntup.pho2phi, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, ntup.var_Pta, ntup.var_MhMZ, ntup.H_pt ]
                    MVA_value = model.predict_proba(MVA_list)[:, 1]
                    if MVA_value < mvaCut:continue

                #histos['Z_befor50'][sample].Fill( ntup.Z_befor50, weight )
                #if ntup.pho1Pt<15: continue

                shift_up = 1. + float(args.sys)
                shift_down = 1. - float(args.sys)

                histos_up['pho1Pt'][sample].Fill( ntup.pho1Pt*shift_up, weight )
                histos_up['pho1eta'][sample].Fill( ntup.pho1eta*shift_up, weight )
                histos_up['pho1phi'][sample].Fill( ntup.pho1phi*shift_up, weight )
                histos_up['pho1R9'][sample].Fill( ntup.pho1R9*shift_up, weight )
                histos_up['pho1IetaIeta'][sample].Fill( ntup.pho1IetaIeta*shift_up, weight )
                histos_up['pho1IetaIeta55'][sample].Fill( ntup.pho1IetaIeta55*shift_up, weight )
                histos_up['pho1PIso_noCorr'][sample].Fill( ntup.pho1PIso_noCorr*shift_up, weight )
                histos_up['pho2Pt'][sample].Fill( ntup.pho2Pt*shift_up, weight )
                histos_up['pho2eta'][sample].Fill( ntup.pho2eta*shift_up, weight )
                histos_up['pho2phi'][sample].Fill( ntup.pho2phi*shift_up, weight )
                histos_up['pho2R9'][sample].Fill( ntup.pho2R9*shift_up, weight )
                histos_up['pho2IetaIeta'][sample].Fill( ntup.pho2IetaIeta*shift_up, weight )
                histos_up['pho2IetaIeta55'][sample].Fill( ntup.pho2IetaIeta55*shift_up, weight )
                histos_up['pho2PIso_noCorr'][sample].Fill( ntup.pho2PIso_noCorr*shift_up, weight )
                histos_up['H_m'][sample].Fill( ntup.H_m*shift_up, weight )
                histos_up['H_pt'][sample].Fill( ntup.H_pt*shift_up, weight )
                histos_up['ALP_m'][sample].Fill( ntup.ALP_m*shift_up, weight )
                histos_up['Z_m'][sample].Fill( ntup.Z_m*shift_up, weight )
                histos_up['var_dR_Za'][sample].Fill( ntup.var_dR_Za*shift_up, weight )
                histos_up['var_dR_g1g2'][sample].Fill( ntup.var_dR_g1g2*shift_up, weight )
                histos_up['var_dR_g1Z'][sample].Fill( ntup.var_dR_g1Z*shift_up, weight )
                histos_up['var_PtaOverMa'][sample].Fill( ntup.var_PtaOverMa*shift_up, weight )
                histos_up['var_PtaOverMh'][sample].Fill( ntup.var_PtaOverMh*shift_up, weight )
                histos_up['var_Pta'][sample].Fill( ntup.var_Pta*shift_up, weight )
                histos_up['var_MhMa'][sample].Fill( ntup.var_MhMa*shift_up, weight )
                histos_up['var_MhMZ'][sample].Fill( ntup.var_MhMZ*shift_up, weight )
                histos_up['ALP_calculatedPhotonIso'][sample].Fill( ntup.ALP_calculatedPhotonIso*shift_up, weight )

                histos_down['pho1Pt'][sample].Fill( ntup.pho1Pt*shift_down, weight )
                histos_down['pho1eta'][sample].Fill( ntup.pho1eta*shift_down, weight )
                histos_down['pho1phi'][sample].Fill( ntup.pho1phi*shift_down, weight )
                histos_down['pho1R9'][sample].Fill( ntup.pho1R9*shift_down, weight )
                histos_down['pho1IetaIeta'][sample].Fill( ntup.pho1IetaIeta*shift_down, weight )
                histos_down['pho1IetaIeta55'][sample].Fill( ntup.pho1IetaIeta55*shift_down, weight )
                histos_down['pho1PIso_noCorr'][sample].Fill( ntup.pho1PIso_noCorr*shift_down, weight )
                histos_down['pho2Pt'][sample].Fill( ntup.pho2Pt*shift_down, weight )
                histos_down['pho2eta'][sample].Fill( ntup.pho2eta*shift_down, weight )
                histos_down['pho2phi'][sample].Fill( ntup.pho2phi*shift_down, weight )
                histos_down['pho2R9'][sample].Fill( ntup.pho2R9*shift_down, weight )
                histos_down['pho2IetaIeta'][sample].Fill( ntup.pho2IetaIeta*shift_down, weight )
                histos_down['pho2IetaIeta55'][sample].Fill( ntup.pho2IetaIeta55*shift_down, weight )
                histos_down['pho2PIso_noCorr'][sample].Fill( ntup.pho2PIso_noCorr*shift_down, weight )
                histos_down['H_m'][sample].Fill( ntup.H_m*shift_down, weight )
                histos_down['H_pt'][sample].Fill( ntup.H_pt*shift_down, weight )
                histos_down['ALP_m'][sample].Fill( ntup.ALP_m*shift_down, weight )
                histos_down['Z_m'][sample].Fill( ntup.Z_m*shift_down, weight )
                histos_down['var_dR_Za'][sample].Fill( ntup.var_dR_Za*shift_down, weight )
                histos_down['var_dR_g1g2'][sample].Fill( ntup.var_dR_g1g2*shift_down, weight )
                histos_down['var_dR_g1Z'][sample].Fill( ntup.var_dR_g1Z*shift_down, weight )
                histos_down['var_PtaOverMa'][sample].Fill( ntup.var_PtaOverMa*shift_down, weight )
                histos_down['var_PtaOverMh'][sample].Fill( ntup.var_PtaOverMh*shift_down, weight )
                histos_down['var_Pta'][sample].Fill( ntup.var_Pta*shift_down, weight )
                histos_down['var_MhMa'][sample].Fill( ntup.var_MhMa*shift_down, weight )
                histos_down['var_MhMZ'][sample].Fill( ntup.var_MhMZ*shift_down, weight )
                histos_down['ALP_calculatedPhotonIso'][sample].Fill( ntup.ALP_calculatedPhotonIso*shift_down, weight )



        ## End of for iEvt in range( ntup.GetEntries() )
    ## End of for sample in analyzer_cfg.samp_names

    ### save raw histograms
    raw_dir_up = out_file.mkdir('plots_up')
    raw_dir_up.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.bkg_names:
            plot_cfg.SetHistStyles(histos_up[var_name][sample], sample)
            histos_up[var_name][sample].Write()


    out_file.cd()
    raw_dir_down = out_file.mkdir('plots_down')
    raw_dir_down.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.bkg_names:
            plot_cfg.SetHistStyles(histos_down[var_name][sample], sample)
            histos_down[var_name][sample].Write()

    ### save stack plots and make ratio plots
    out_file.cd()



    print '\n\n'
    #CountYield(analyzer_cfg, histos['ALP_m'])
    out_file.Close()

    print 'Done'


main()
