####################################################
####################################################

import os
import sys
import numpy as np

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from Plot_Helper import LoadNtuples, MakeStack, CreateCanvas, DrawOnCanv, SaveCanvPic, MakeLumiLabel, MakeCMSDASLabel, ScaleSignal, MakeRatioPlot, MakeLegend
from Analyzer_Helper import MuPair, DimuCandidates, IsBtagLoose, IsBtagMed, IsBtagTight
import Analyzer_Configs as AC
import Plot_Configs     as PC

from Analyzer_ALP import PIso2D, plot2D, plot2D_CONT, CountNormalizationYield

import CMS_lumi, tdrstyle

#from xgboost import XGBClassifier
#import pickle

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
parser.add_argument('-C', '--CR', dest='CR', action='store_true', default=False, help='make control region')
parser.add_argument('-m', '--mva', dest='mva', action='store_true', default=False, help='use mva or not')
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
    name = mass

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

    analyzer_cfg = AC.Analyzer_Config('inclusive', mass,args.year)

    ############ CR plots#########
    #if args.CR:
        #analyzer_cfg.sample_loc = analyzer_cfg.sample_loc.replace('massInde','massInde/CR')
        #analyzer_cfg.sig_names  = ['']
    ##############################
    analyzer_cfg.Print_Config()
    ntuples = LoadNtuples(analyzer_cfg)


    #sys_names = ['pho_norm','pho_SFs_Sys_up','pho_SFs_Sys_dn','pho_PUWeight_Sys_up','pho_PUWeight_Sys_dn','lep_dataMC_up','lep_dataMC_dn']
    sys_names = ['pho_norm','pho_SFs_Sys_up','pho_SFs_Sys_dn','pho_PUWeight_Sys_up','pho_PUWeight_Sys_dn','lep_dataMC_up','lep_dataMC_dn']
    var_names = ['ALP_m']
    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)

    ### declare histograms

    histos = {}

    for sample in analyzer_cfg.sig_names:
        histos[sample] = {}
        for sys_name in sys_names:
            histos[sample][sys_name] = TH1F(sys_name + 'ALP_m' + '_' + sample, sys_name + 'ALP_m' + '_' + sample, 50, 0., 40.)


    for sys_name in sys_names:
        print '##############'
        print 'Systematics : ' +  sys_name
        print '##############'
        ### loop over samples and events
        for sample in analyzer_cfg.sig_names:
            ntup = ntuples[sample] # just a short name
            print '\n\nOn sample: %s' %sample
            print 'total events: %d' %ntup.GetEntries()

            for iEvt in range( ntup.GetEntries() ):
                ntup.GetEvent(iEvt)
                #if (iEvt == 1): break


                if (iEvt % 100000 == 1):
                    print "looking at event %d" %iEvt



                if sys_name =='pho_norm':
                    weight = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs
                elif sys_name =='pho_SFs_Sys_up':
                    #weight = ntup.factor * (ntup.pho1SFs+ntup.pho1SFs_sys) * (ntup.pho2SFs+ntup.pho2SFs_sys)
                    weight = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * (ntup.pho1SFs+ntup.pho1SFs_sys_replaced) * (ntup.pho2SFs+ntup.pho2SFs_sys_replaced)
                elif sys_name =='pho_SFs_Sys_dn':
                    #weight = ntup.factor * (ntup.pho1SFs-ntup.pho1SFs_sys) * (ntup.pho2SFs-ntup.pho2SFs_sys)
                    weight = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * (ntup.pho1SFs-ntup.pho1SFs_sys_replaced) * (ntup.pho2SFs-ntup.pho2SFs_sys_replaced)
                elif sys_name =='pho_PUWeight_Sys_up':
                    weight = ntup.event_genWeight * ntup.event_pileupWeightUp * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs
                elif sys_name =='pho_PUWeight_Sys_dn':
                    weight = ntup.event_genWeight * ntup.event_pileupWeightDn * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs
                elif sys_name =='lep_dataMC_up':
                    weight = ntup.event_genWeight * ntup.event_pileupWeight * (ntup.l1_dataMC+ntup.l1_dataMCErr) * (ntup.l2_dataMC+ntup.l2_dataMCErr) * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs
                elif sys_name =='lep_dataMC_dn':
                    weight = ntup.event_genWeight * ntup.event_pileupWeight * (ntup.l1_dataMC-ntup.l1_dataMCErr) * (ntup.l2_dataMC-ntup.l2_dataMCErr) * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs


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

                    histos[sample][sys_name].Fill( ntup.ALP_m, weight )



    print '\n\n'
    CountNormalizationYield(analyzer_cfg, histos, sys_names)

    print 'Done'


main()
