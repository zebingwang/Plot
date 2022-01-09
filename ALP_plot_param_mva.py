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
parser.add_argument('-m', '--mva', dest='mva', action='store_true', default=False, help='use mva or not')
parser.add_argument('-b', '--blind', dest='blind', action='store_true', default=False, help='Blind signal region')
args = parser.parse_args()

mass = 'param'



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
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2018.pkl"
    mvaCut = 0.9766
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
    var_names = ['Z_m', 'H_m', 'ALP_m', 'H_pt', 'mvaVal']
    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)

    ### declare histograms

    histos = {}

    for var_name in var_names:
        histos[var_name] = {}
    for sample in analyzer_cfg.samp_names:
        histos['Z_m'][sample]    = TH1F('Z_m'    + '_' + sample, 'Z_m'    + '_' + sample, 100,  50., 130.)
        histos['H_m'][sample]    = TH1F('H_m'    + '_' + sample, 'H_m'    + '_' + sample, 20,  118., 140.)
        histos['H_pt'][sample]    = TH1F('H_pt'    + '_' + sample, 'H_pt'    + '_' + sample, 100,  0., 200.)
        histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 50, 0., 40.)

        histos['mvaVal'][sample]    = TH1F('mvaVal'    + '_' + sample, 'mvaVal'    + '_' + sample, 50,  -0.1, 1.1)


    ### loop over samples and events
    ntup_all = {}
    mass_list = {'M5':5.0, 'M15':15.0, 'M30':30.0}
    for sample in analyzer_cfg.sig_names:
        ntup_all[sample] = ntuples[sample] # just a short name
        ntup_all['data'] = ntuples['data']

        for sOrb in [sample,'data']:
            ntup = ntup_all[sOrb]
            print '\n\nOn sample: %s' %sOrb
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
                        if ntup.H_m>180. or ntup.H_m<130.: continue
                    else:
                        if ntup.H_m<130. and ntup.H_m>118.: continue

                    #if abs(ntup.l1_id)!=13: continue

                    #print 'pho1: pt:' + str(ntup.pho1Pt) + ', eta: ' + str(ntup.pho1eta) + ', SFs: ' + str(pho1_SFs)

                    param = (ntup.ALP_m - mass_list[sample])/ntup.H_m

                    if mva:
                        MVA_list = [ntup.pho1Pt, ntup.pho1eta, ntup.pho1phi, ntup.pho1R9, ntup.pho1IetaIeta55 ,ntup.pho2Pt, ntup.pho2eta, ntup.pho2phi, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, ntup.var_Pta, ntup.var_MhMZ, ntup.H_pt ,param]
                        MVA_value = model.predict_proba(MVA_list)[:, 1]
                        if MVA_value < mvaCut:continue


                    histos['H_m'][sOrb].Fill( ntup.H_m, weight )
                    histos['H_pt'][sOrb].Fill( ntup.H_pt, weight )
                    histos['ALP_m'][sOrb].Fill( ntup.ALP_m, weight )
                    histos['Z_m'][sOrb].Fill( ntup.Z_m, weight )


                    if args.blind:
                        if sOrb == 'data' and ntup.passBDT > 0.5:
                            continue
                        else:
                            histos['mvaVal'][sOrb].Fill( ntup.Val_BDT, weight )
                    else:
                        histos['mvaVal'][sOrb].Fill( ntup.Val_BDT, weight )



            ## End of for iEvt in range( ntup.GetEntries() )

        ### save raw histograms
        raw_dir = out_file.mkdir('raw_plots_'+sample)
        raw_dir.cd()
        for var_name in var_names:
            for sam in [sample,'data']:
                plot_cfg.SetHistStyles(histos[var_name][sam], sam)
                histos[var_name][sam].Write()



    print '\n\n'
    out_file.Close()


    print 'Done'


main()
