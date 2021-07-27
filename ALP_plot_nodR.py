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

#from xgboost import XGBClassifier
#import pickle

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
parser.add_argument('-C', '--CR', dest='CR', action='store_true', default=False, help='make control region')
parser.add_argument('-m', '--mva', dest='mva', action='store_true', default=False, help='use mva or not')
parser.add_argument("--blind", dest="blind", default="store_true", default=False, help="Blind signal region?")
args = parser.parse_args()

mass = 'massIndependent'



gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

mva = args.mva
if args.CR:
    name = mass + '_CR'
elif mva:
    name = mass + '_mva'
else:
    name = mass+ '_SFs'

if args.year == '2016':
    file_out = 'plots_16'
    out_name = "ALP_plot_data16_"+name+".root"
    #BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_massindependent_2016.pkl"
    #mvaCut = 0.8675
elif args.year == '2017':
    file_out = 'plots_17'
    out_name = "ALP_plot_data17_"+name+".root"
    #BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_massindependent_2017.pkl"
    #mvaCut = 0.8365
elif args.year == '2018':
    file_out = 'plots_18'
    out_name = "ALP_plot_data18_"+name+".root"
    #BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_massindependent_2018.pkl"
    #mvaCut = 0.8923
else:
    print "do not include at 2016/2017/2018"
    exit(0)

file_plot = file_out + "/plot_"+name
# load the model from disk
#model = pickle.load(open(BDT_filename, 'rb'))

def main():

    if not os.path.exists(file_out):
        os.makedirs(file_out)

    if not os.path.exists(file_plot):
        os.makedirs(file_plot)

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
    var_names = ['Z_m', 'H_m', 'ALP_m','pho1Pt', 'pho1eta', 'pho1phi', 'pho1R9', 'pho1IetaIeta', 'pho1IetaIeta55','pho1PIso_noCorr' ,'pho2Pt', 'pho2eta', 'pho2phi', 'pho2R9', 'pho2IetaIeta', 'pho2IetaIeta55','pho2PIso_noCorr','ALP_calculatedPhotonIso', 'var_dR_Za', 'var_dR_g1g2', 'var_dR_g1Z', 'var_PtaOverMh', 'var_Pta', 'var_MhMZ', 'H_pt', 'var_PtaOverMa', 'var_MhMa', 'mvaVal']
    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)

    ### declare histograms

    histos = {}

    for var_name in var_names:
        histos[var_name] = {}
    for sample in analyzer_cfg.samp_names:
        histos['pho1Pt'][sample]    = TH1F('pho1Pt'    + '_' + sample, 'pho1Pt'    + '_' + sample, 50,  5., 50.)
        histos['pho1eta'][sample]    = TH1F('pho1eta'    + '_' + sample, 'pho1eta'    + '_' + sample, 100,  -3., 3.)
        histos['pho1phi'][sample]    = TH1F('pho1phi'    + '_' + sample, 'pho1phi'    + '_' + sample, 100,  -4., 4.)
        histos['pho1R9'][sample]    = TH1F('pho1R9'    + '_' + sample, 'pho1R9'    + '_' + sample, 100,  0., 1.)
        histos['pho1IetaIeta'][sample]    = TH1F('pho1IetaIeta'    + '_' + sample, 'pho1IetaIeta'    + '_' + sample, 100,  0., 0.07)
        histos['pho1IetaIeta55'][sample]    = TH1F('pho1IetaIeta55'    + '_' + sample, 'pho1IetaIeta55'    + '_' + sample, 100,  0., 0.07)
        histos['pho1PIso_noCorr'][sample]    = TH1F('pho1PIso_noCorr'    + '_' + sample, 'pho1PIso_noCorr'    + '_' + sample, 100, 0., 40.)
        histos['pho2Pt'][sample]    = TH1F('pho2Pt'    + '_' + sample, 'pho2Pt'    + '_' + sample, 50,  5., 30.)
        histos['pho2eta'][sample]    = TH1F('pho2eta'    + '_' + sample, 'pho2eta'    + '_' + sample, 100,  -3., 3.)
        histos['pho2phi'][sample]    = TH1F('pho2phi'    + '_' + sample, 'pho2phi'    + '_' + sample, 100,  -4., 4.)
        histos['pho2R9'][sample]    = TH1F('pho2R9'    + '_' + sample, 'pho2R9'    + '_' + sample, 100,  0., 1.)
        histos['pho2IetaIeta'][sample]    = TH1F('pho2IetaIeta'    + '_' + sample, 'pho2IetaIeta'    + '_' + sample, 100,  0., 0.07)
        histos['pho2IetaIeta55'][sample]    = TH1F('pho2IetaIeta55'    + '_' + sample, 'pho2IetaIeta55'    + '_' + sample, 100,  0., 0.07)
        histos['pho2PIso_noCorr'][sample]    = TH1F('pho2PIso_noCorr'    + '_' + sample, 'pho2PIso_noCorr'    + '_' + sample, 100, 0., 40.)
        histos['Z_m'][sample]    = TH1F('Z_m'    + '_' + sample, 'Z_m'    + '_' + sample, 100,  50., 130.)
        histos['H_m'][sample]    = TH1F('H_m'    + '_' + sample, 'H_m'    + '_' + sample, 20,  118., 140.)
        histos['H_pt'][sample]    = TH1F('H_pt'    + '_' + sample, 'H_pt'    + '_' + sample, 100,  0., 200.)
        histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 50, 0., 40.)
        histos['var_dR_g1g2'][sample] = TH1F('var_dR_g1g2' + '_' + sample, 'var_dR_g1g2' + '_' + sample, 50, 0., 5)
        histos['var_PtaOverMa'][sample] = TH1F('var_PtaOverMa' + '_' + sample, 'var_PtaOverMa' + '_' + sample, 100, 0., 100.)
        histos['var_dR_Za'][sample] = TH1F('var_dR_Za' + '_' + sample, 'var_dR_Za' + '_' + sample, 50, 0., 7.)
        histos['var_dR_g1Z'][sample] = TH1F('var_dR_g1Z' + '_' + sample, 'var_dR_g1Z' + '_' + sample, 50, 0., 8.)
        histos['var_PtaOverMh'][sample] = TH1F('var_PtaOverMh' + '_' + sample, 'var_PtaOverMh' + '_' + sample, 100, 0., 1.)
        histos['var_Pta'][sample] = TH1F('var_Pta' + '_' + sample, 'var_Pta' + '_' + sample, 50, 0., 60.)
        histos['var_MhMa'][sample] = TH1F('var_MhMa' + '_' + sample, 'var_MhMa' + '_' + sample, 50, 100., 200.)
        histos['var_MhMZ'][sample] = TH1F('var_MhMZ' + '_' + sample, 'var_MhMZ' + '_' + sample, 50, 120., 250.)
        histos['ALP_calculatedPhotonIso'][sample] = TH1F('ALP_calculatedPhotonIso' + '_' + sample, 'ALP_calculatedPhotonIso' + '_' + sample, 100, 0., 5.)

        histos['mvaVal'][sample]    = TH1F('mvaVal'    + '_' + sample, 'mvaVal'    + '_' + sample, 50,  -0.1, 1.1)


    ### loop over samples and events
    for sample in analyzer_cfg.samp_names:
        ntup = ntuples[sample] # just a short name
        print '\n\nOn sample: %s' %sample
        print 'total events: %d' %ntup.GetEntries()

        for iEvt in range( ntup.GetEntries() ):
            ntup.GetEvent(iEvt)
            #if (iEvt == 1): break


            if (iEvt % 100000 == 1):
                print "looking at event %d" %iEvt


            weight = ntup.factor * ntup.pho1SFs * ntup.pho2SFs

            if (ntup.H_m > -90):

                #if ntup.dR_pho < 0.02: continue
                #if not ntup.passEleVeto: continue
                if not ntup.passChaHadIso: continue
                if not ntup.passNeuHadIso: continue
                if not ntup.passdR_gl: continue
                if not ntup.passHOverE: continue
                if not args.CR:
                    if ntup.H_m>130. or ntup.H_m<118.: continue
                else:
                    if ntup.H_m<130. and ntup.H_m>118.: continue

                #if abs(ntup.l1_id)!=13: continue

                #print 'pho1: pt:' + str(ntup.pho1Pt) + ', eta: ' + str(ntup.pho1eta) + ', SFs: ' + str(pho1_SFs)
                if mva:
                    #MVA_list = [ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_PtaOverMh, ntup.var_MhMZ, ntup.var_dR_g1g2]

                    #MVA_list = [ntup.pho1Pt, ntup.pho1eta, ntup.pho1phi, ntup.pho1R9, ntup.pho1IetaIeta55 ,ntup.pho2Pt, ntup.pho2eta, ntup.pho2phi, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, ntup.var_Pta, ntup.var_MhMZ, ntup.H_pt ]
                    #MVA_value = model.predict_proba(MVA_list)[:, 1]
                    #if MVA_value < mvaCut:continue
                    if ntup.passBDT < 0.5: continue


                #histos['Z_befor50'][sample].Fill( ntup.Z_befor50, weight )
                #if ntup.pho1Pt<15: continue
                histos['pho1Pt'][sample].Fill( ntup.pho1Pt, weight )
                histos['pho1eta'][sample].Fill( ntup.pho1eta, weight )
                histos['pho1phi'][sample].Fill( ntup.pho1phi, weight )
                histos['pho1R9'][sample].Fill( ntup.pho1R9, weight )
                histos['pho1IetaIeta'][sample].Fill( ntup.pho1IetaIeta, weight )
                histos['pho1IetaIeta55'][sample].Fill( ntup.pho1IetaIeta55, weight )
                histos['pho1PIso_noCorr'][sample].Fill( ntup.pho1PIso_noCorr, weight )
                histos['pho2Pt'][sample].Fill( ntup.pho2Pt, weight )
                histos['pho2eta'][sample].Fill( ntup.pho2eta, weight )
                histos['pho2phi'][sample].Fill( ntup.pho2phi, weight )
                histos['pho2R9'][sample].Fill( ntup.pho2R9, weight )
                histos['pho2IetaIeta'][sample].Fill( ntup.pho2IetaIeta, weight )
                histos['pho2IetaIeta55'][sample].Fill( ntup.pho2IetaIeta55, weight )
                histos['pho2PIso_noCorr'][sample].Fill( ntup.pho2PIso_noCorr, weight )


                histos['H_m'][sample].Fill( ntup.H_m, weight )
                histos['H_pt'][sample].Fill( ntup.H_pt, weight )
                histos['ALP_m'][sample].Fill( ntup.ALP_m, weight )
                histos['Z_m'][sample].Fill( ntup.Z_m, weight )

                histos['var_dR_Za'][sample].Fill( ntup.var_dR_Za, weight )
                histos['var_dR_g1g2'][sample].Fill( ntup.var_dR_g1g2, weight )
                histos['var_dR_g1Z'][sample].Fill( ntup.var_dR_g1Z, weight )
                histos['var_PtaOverMa'][sample].Fill( ntup.var_PtaOverMa, weight )
                histos['var_PtaOverMh'][sample].Fill( ntup.var_PtaOverMh, weight )
                histos['var_Pta'][sample].Fill( ntup.var_Pta, weight )
                histos['var_MhMa'][sample].Fill( ntup.var_MhMa, weight )
                histos['var_MhMZ'][sample].Fill( ntup.var_MhMZ, weight )
                histos['ALP_calculatedPhotonIso'][sample].Fill( ntup.ALP_calculatedPhotonIso, weight )

                if args.blind:
                    if sample == 'data' and ntup.passBDT > 0.5:
                        continue
                    else:
                        histos['mvaVal'][sample].Fill( ntup.Val_BDT, weight )
                else:
                    histos['mvaVal'][sample].Fill( ntup.Val_BDT, weight )



        ## End of for iEvt in range( ntup.GetEntries() )
    ## End of for sample in analyzer_cfg.samp_names

    ### save raw histograms
    raw_dir = out_file.mkdir('raw_plots')
    raw_dir.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.samp_names:
            plot_cfg.SetHistStyles(histos[var_name][sample], sample)
            histos[var_name][sample].Write()

    ### save stack plots and make ratio plots
    out_file.cd()
    lumi_label = MakeLumiLabel(plot_cfg.lumi)
    cms_label  = MakeCMSDASLabel()

    scaled_sig = {}
    for var_name in var_names:
        stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
        #scaled_sig = 0
        #ratio_plot = 0
        for sample in analyzer_cfg.sig_names:
            scaled_sig[sample] = ScaleSignal(plot_cfg, stacks[sample], histos[var_name][sample], var_name)
        ratio_plot = MakeRatioPlot(histos[var_name]['data'], stacks['all'].GetStack().Last(), var_name)
        legend = MakeLegend(plot_cfg, histos[var_name], scaled_sig)

        canv = CreateCanvas(var_name)
        DrawOnCanv(canv, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label,logY=False)

        canv.Write()
        SaveCanvPic(canv, file_plot, var_name)


    for var_name in var_names:
        stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
        #scaled_sig = 0
        #ratio_plot = 0
        for sample in analyzer_cfg.sig_names:
            scaled_sig[sample] = ScaleSignal(plot_cfg, stacks[sample], histos[var_name][sample], var_name)
        ratio_plot = MakeRatioPlot(histos[var_name]['data'], stacks['all'].GetStack().Last(), var_name)
        legend = MakeLegend(plot_cfg, histos[var_name], scaled_sig)

        canv = CreateCanvas(var_name)
        DrawOnCanv(canv, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label, logY=True)

        canv.Write()
        SaveCanvPic(canv, file_plot, var_name+'_log')


    print '\n\n'
    #CountYield(analyzer_cfg, histos['ALP_m'])
    out_file.Close()

    #for sample in analyzer_cfg.samp_names:
        #print sample + '\t\t' + str(histos['mvaVal'][sample].Integral())

    print 'Done'


main()
