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

from Analyzer_ALP import PIso2D, CountNormalizationYield

import CMS_lumi, tdrstyle

from xgboost import XGBClassifier
import pickle

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
parser.add_argument('-C', '--CR', dest='CR', action='store_true', default=False, help='make control region')
parser.add_argument('-m', '--mva', dest='mva', action='store_true', default=False, help='use mva or not')
parser.add_argument('--ele', dest='ele', action='store_true', default=False, help='electron channel?')
parser.add_argument('--mu', dest='mu', action='store_true', default=False, help='muon channel?')
args = parser.parse_args()

gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

mva = args.mva


BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/UL/model_ALP_BDT_param.pkl"
mvaCuts = {'M1':0.955, 'M2':0.98, 'M3':0.985, 'M4':0.98, 'M5':0.985, 'M6':0.99, 'M7':0.985, 'M8':0.99, 'M9':0.99, 'M10':0.99, 'M15':0.99, 'M20':0.99, 'M25':0.985, 'M30':0.98}

# load the model from disk
model = pickle.load(open(BDT_filename, 'rb'))

def main():

    analyzer_cfg = AC.Analyzer_Config('inclusive', args.year)

    ############ CR plots#########
    #if args.CR:
        #analyzer_cfg.sample_loc = analyzer_cfg.sample_loc.replace('massInde','massInde/CR')
        #analyzer_cfg.sig_names  = ['']
    ##############################
    analyzer_cfg.Print_Config()
    ntuples = LoadNtuples(analyzer_cfg)


    sys_names = ['pho_norm','CMS_IDeff_g_up','CMS_IDeff_g_dn','CMS_pileup_up','CMS_pileup_dn','CMS_IDeff_lep_up','CMS_IDeff_lep_dn']
    #sys_names = ['pho_norm','CMS_IDeff_g_up','CMS_IDeff_g_dn','CMS_pileup_up','CMS_pileup_dn','CMS_IDeff_lep_up','CMS_IDeff_lep_dn','pho_SFs_dR0','pho_SFs_dR0P1','pho_SFs_dR0P15','pho_SFs_dR0P2','pho_SFs_dR0P25','pho_SFs_dR0P3']
    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)

    ### declare histograms

    histos = {}

    for sample in analyzer_cfg.sig_names:
        histos[sample] = {}
        for sys_name in sys_names:
            histos[sample][sys_name] = TH1F(sys_name + '_' + sample, sys_name + '_' + sample, 50, 110., 180.)


    mass_list = {'M1':1.0, 'M2':2.0, 'M3':3.0, 'M4':4.0, 'M5':5.0, 'M6':6.0, 'M7':7.0, 'M8':8.0, 'M9':9.0, 'M10':10.0, 'M15':15.0, 'M20':20.0, 'M25':25.0, 'M30':30.0}
    
    ### loop over samples and events
    for sample in analyzer_cfg.sig_names:
        ntup = ntuples[sample] # just a short name
        print '\n\nOn sample: %s' %sample
        print 'total events: %d' %ntup.GetEntries()

        for iEvt in range( ntup.GetEntries() ):
            ntup.GetEvent(iEvt)
            #if (iEvt == 1): break

            if args.ele:
                if abs(ntup.l1_id) == 13: 
                    continue
            if args.mu:
                if abs(ntup.l1_id) == 11: 
                    continue


            if (iEvt % 100000 == 1):
                print "looking at event %d" %iEvt

            if (ntup.H_m > -90):

                #if ntup.dR_pho < 0.02: continue
                #if not ntup.passEleVeto: continue
                if not ntup.passChaHadIso: continue
                if not ntup.passNeuHadIso: continue
                if not ntup.passdR_gl: continue
                if not ntup.passHOverE: continue

                if ntup.H_m>180. or ntup.H_m<110.: continue


                #if abs(ntup.l1_id)!=13: continue

                #print 'pho1: pt:' + str(ntup.pho1Pt) + ', eta: ' + str(ntup.pho1eta) + ', SFs: ' + str(pho1_SFs)
                if mva:

                    param = (ntup.ALP_m - mass_list[sample])/ntup.H_m

                    MVA_list = [ntup.pho1Pt, ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr ,ntup.pho2Pt, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.pho2PIso_noCorr,ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, ntup.var_PtaOverMh, ntup.H_pt, param]
                    MVA_value = model.predict_proba(MVA_list)[:, 1]

                    if MVA_value < mvaCuts[sample]: continue

                for sys_name in sys_names:

                    if sys_name =='pho_norm':
                        weight = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs
                    elif sys_name =='CMS_IDeff_g_up':
                        #weight = ntup.factor * (ntup.pho1SFs+ntup.pho1SFs_sys) * (ntup.pho2SFs+ntup.pho2SFs_sys)
                        weight = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * (ntup.pho1SFs+ntup.pho1SFs_sys) * (ntup.pho2SFs+ntup.pho2SFs_sys)
                    elif sys_name =='CMS_IDeff_g_dn':
                        #weight = ntup.factor * (ntup.pho1SFs-ntup.pho1SFs_sys) * (ntup.pho2SFs-ntup.pho2SFs_sys)
                        weight = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * (ntup.pho1SFs-ntup.pho1SFs_sys) * (ntup.pho2SFs-ntup.pho2SFs_sys)
                    elif sys_name =='CMS_pileup_up':
                        weight = ntup.event_genWeight * ntup.event_pileupWeightUp * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs
                    elif sys_name =='CMS_pileup_dn':
                        weight = ntup.event_genWeight * ntup.event_pileupWeightDn * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs
                    elif sys_name =='CMS_IDeff_lep_up':
                        weight = ntup.event_genWeight * ntup.event_pileupWeight * (ntup.l1_dataMC+ntup.l1_dataMCErr) * (ntup.l2_dataMC+ntup.l2_dataMCErr) * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs
                    elif sys_name =='CMS_IDeff_lep_dn':
                        weight = ntup.event_genWeight * ntup.event_pileupWeight * (ntup.l1_dataMC-ntup.l1_dataMCErr) * (ntup.l2_dataMC-ntup.l2_dataMCErr) * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs
                    elif sys_name =='pdf_Higgs_gg_up':
                        weight = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs*ntup.event_pdfweight_up
                    elif sys_name =='pdf_Higgs_gg_dn':
                        weight = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs*ntup.event_pdfweight_dn
                    elif sys_name =='QCDscale_ggH_up':
                        weight = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs*ntup.event_qcdweight_up
                    elif sys_name =='QCDscale_ggH_dn':
                        weight = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs*ntup.event_qcdweight_dn
                    else:
                        print "Sys name do not exist!"



                    histos[sample][sys_name].Fill( ntup.H_m, weight )



    print '\n\n'
    if args.ele:
        CountNormalizationYield(analyzer_cfg, histos, sys_names, 'ele', args.year)   
    elif args.mu:
        CountNormalizationYield(analyzer_cfg, histos, sys_names, 'mu', args.year)
    else:
        CountNormalizationYield(analyzer_cfg, histos, sys_names, 'combine', args.year)

    print 'Done'


main()
