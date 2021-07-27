####################################################
####################################################

import os
import sys
import numpy as np

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
import ROOT
from Plot_Helper import LoadNtuples, MakeStack, CreateCanvas, DrawOnCanv, SaveCanvPic, MakeLumiLabel, MakeCMSDASLabel, ScaleSignal, MakeRatioPlot, MakeLegend
from Analyzer_Helper import MuPair, DimuCandidates, IsBtagLoose, IsBtagMed, IsBtagTight, CountYield
import Analyzer_Configs as AC
import Plot_Configs     as PC

from Analyzer_ALP import PIso2D, plot2D, plot2D_CONT, CountYield

import CMS_lumi, tdrstyle
import multiprocessing

from xgboost import XGBClassifier
import pickle

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
parser.add_argument('-C', '--CR', dest='CR', action='store_true', default=False, help='make control region')
parser.add_argument('-m', '--mva', dest='mva', action='store_true', default=False, help='use mva or not')
parser.add_argument("-s", "--sysName", dest="sysName", default="norm", help="name of the systematics")
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
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_massindependent_2016.pkl"
    mvaCut = 0.9381
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

def Cal_mvaVal(*mvalist):
    return model.predict_proba(list(mvalist))[:, 1]


def main():

    if not os.path.exists(file_out):
        os.makedirs(file_out)

    if not os.path.exists(file_plot):
        os.makedirs(file_plot)

    out_file = TFile( file_out + '/' + out_name , "RECREATE")

    analyzer_cfg = AC.Analyzer_Config('inclusive', mass,args.year)
    analyzer_cfg.Print_Config()
    ntuples = LoadNtuples(analyzer_cfg)



    var_names = ['mvaVal']
    #sys_names = ['ShowerShape', 'pho_scale', 'pho_smear', 'lep_scale', 'lep_smear']
    sys_names = []
    sys_names.append(args.sysName)

    MVA_list = {}

    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)

    ### declare histograms

    histos = {}

    for sample in analyzer_cfg.samp_names:
        histos[sample] = {}
        histos[sample]['norm'] = TH1F('norm_mvaVal'    + '_' + sample, 'norm_mvaVal'    + '_' + sample, 50,  -0.1, 1.1)
        for sys_name in sys_names:
            if sys_name == 'norm': continue
            histos[sample][sys_name] = {}
            histos[sample][sys_name]['up'] = TH1F(sys_name + '_up_' + 'mvaVal'    + '_' + sample, sys_name + '_up_' + 'mvaVal'    + '_' + sample, 50,  -0.1, 1.1)
            histos[sample][sys_name]['dn'] = TH1F(sys_name + '_dn_' + 'mvaVal'    + '_' + sample, sys_name + '_dn_' + 'mvaVal'    + '_' + sample, 50,  -0.1, 1.1)


    ### loop over samples and events
    for sample in analyzer_cfg.samp_names:
        ntup = ntuples[sample] # just a short name
        print '\n\nOn sample: %s' %sample
        print 'total events: %d' %ntup.GetEntries()

        for iEvt in range( ntup.GetEntries() ):
            ntup.GetEvent(iEvt)
            #if (iEvt == 10000): break


            if (iEvt % 100000 == 1):
                print "looking at event %d" %iEvt


            weight = ntup.factor

            if (ntup.H_m > -90):

                if not ntup.passChaHadIso: continue
                if not ntup.passNeuHadIso: continue
                if not ntup.passdR_gl: continue
                if not ntup.passHOverE: continue
                if not args.CR:
                    if ntup.H_m>130. or ntup.H_m<118.: continue
                else:
                    if ntup.H_m<130. and ntup.H_m>118.: continue



                pho1_norm = ROOT.TLorentzVector()
                pho2_norm = ROOT.TLorentzVector()
                pho1_scaleup = ROOT.TLorentzVector()
                pho2_scaleup = ROOT.TLorentzVector()
                pho1_scaledn = ROOT.TLorentzVector()
                pho2_scaledn = ROOT.TLorentzVector()
                pho1_smearup = ROOT.TLorentzVector()
                pho2_smearup = ROOT.TLorentzVector()
                pho1_smeardn = ROOT.TLorentzVector()
                pho2_smeardn = ROOT.TLorentzVector()
                l1_norm = ROOT.TLorentzVector()
                l2_norm = ROOT.TLorentzVector()
                l1_scaleup = ROOT.TLorentzVector()
                l2_scaleup = ROOT.TLorentzVector()
                l1_scaledn = ROOT.TLorentzVector()
                l2_scaledn = ROOT.TLorentzVector()
                l1_smearup = ROOT.TLorentzVector()
                l2_smearup = ROOT.TLorentzVector()
                l1_smeardn = ROOT.TLorentzVector()
                l2_smeardn = ROOT.TLorentzVector()

                pho1_norm.SetPtEtaPhiM(ntup.pho1Pt, ntup.pho1eta, ntup.pho1phi, 0.)
                pho1_scaleup.SetPtEtaPhiM(ntup.pho1scaleup, ntup.pho1eta, ntup.pho1phi, 0.)
                pho1_scaledn.SetPtEtaPhiM(ntup.pho1scaledn, ntup.pho1eta, ntup.pho1phi, 0.)
                pho1_smearup.SetPtEtaPhiM(ntup.pho1smearup, ntup.pho1eta, ntup.pho1phi, 0.)
                pho1_smeardn.SetPtEtaPhiM(ntup.pho1smeardn, ntup.pho1eta, ntup.pho1phi, 0.)
                pho2_norm.SetPtEtaPhiM(ntup.pho2Pt, ntup.pho2eta, ntup.pho2phi, 0.)
                pho2_scaleup.SetPtEtaPhiM(ntup.pho2scaleup, ntup.pho2eta, ntup.pho2phi, 0.)
                pho2_scaledn.SetPtEtaPhiM(ntup.pho2scaledn, ntup.pho2eta, ntup.pho2phi, 0.)
                pho2_smearup.SetPtEtaPhiM(ntup.pho2smearup, ntup.pho2eta, ntup.pho2phi, 0.)
                pho2_smeardn.SetPtEtaPhiM(ntup.pho2smeardn, ntup.pho2eta, ntup.pho2phi, 0.)

                l1_norm.SetPtEtaPhiM(ntup.l1_pt, ntup.l1_eta, ntup.l1_phi, ntup.l1_mass)
                l1_scaleup.SetPtEtaPhiM(ntup.l1_scaleup, ntup.l1_eta, ntup.l1_phi, ntup.l1_mass)
                l1_scaledn.SetPtEtaPhiM(ntup.l1_scaledn, ntup.l1_eta, ntup.l1_phi, ntup.l1_mass)
                l1_smearup.SetPtEtaPhiM(ntup.l1_smearup, ntup.l1_eta, ntup.l1_phi, ntup.l1_mass)
                l1_smeardn.SetPtEtaPhiM(ntup.l1_smeardn, ntup.l1_eta, ntup.l1_phi, ntup.l1_mass)
                l2_norm.SetPtEtaPhiM(ntup.l2_pt, ntup.l2_eta, ntup.l2_phi, ntup.l2_mass)
                l2_scaleup.SetPtEtaPhiM(ntup.l2_scaleup, ntup.l2_eta, ntup.l2_phi, ntup.l2_mass)
                l2_scaledn.SetPtEtaPhiM(ntup.l2_scaledn, ntup.l2_eta, ntup.l2_phi, ntup.l2_mass)
                l2_smearup.SetPtEtaPhiM(ntup.l2_smearup, ntup.l2_eta, ntup.l2_phi, ntup.l2_mass)
                l2_smeardn.SetPtEtaPhiM(ntup.l2_smeardn, ntup.l2_eta, ntup.l2_phi, ntup.l2_mass)

                Pta_norm = ROOT.TLorentzVector()
                Pta_pho_scale_up = ROOT.TLorentzVector()
                Pta_pho_scale_dn = ROOT.TLorentzVector()
                Pta_pho_smear_up = ROOT.TLorentzVector()
                Pta_pho_smear_dn = ROOT.TLorentzVector()

                Pta_norm = (pho1_norm + pho2_norm).Pt()
                Pta_pho_scale_up = (pho1_scaleup + pho2_scaleup).Pt()
                Pta_pho_scale_dn = (pho1_scaledn + pho2_scaledn).Pt()
                Pta_pho_smear_up = (pho1_smearup + pho2_smearup).Pt()
                Pta_pho_smear_dn = (pho1_smeardn + pho2_smeardn).Pt()


                MhMZ_norm = (pho1_norm + pho2_norm + l1_norm + l2_norm).M() + (l1_norm + l2_norm).M()
                MhMZ_pho_scale_up = (pho1_scaleup + pho2_scaleup + l1_norm + l2_norm).M() + (l1_norm + l2_norm).M()
                MhMZ_pho_scale_dn = (pho1_scaledn + pho2_scaledn + l1_norm + l2_norm).M() + (l1_norm + l2_norm).M()
                MhMZ_pho_smear_up = (pho1_smearup + pho2_smearup + l1_norm + l2_norm).M() + (l1_norm + l2_norm).M()
                MhMZ_pho_smear_dn = (pho1_smeardn + pho2_smeardn + l1_norm + l2_norm).M() + (l1_norm + l2_norm).M()
                MhMZ_lep_scale_up = (pho1_norm + pho2_norm + l1_scaleup + l2_scaleup).M() + (l1_scaleup + l2_scaleup).M()
                MhMZ_lep_scale_dn = (pho1_norm + pho2_norm + l1_scaledn + l2_scaledn).M() + (l1_scaledn + l2_scaledn).M()
                MhMZ_lep_smear_up = (pho1_norm + pho2_norm + l1_smearup + l2_smearup).M() + (l1_smearup + l2_smearup).M()
                MhMZ_lep_smear_dn = (pho1_norm + pho2_norm + l1_smeardn + l2_smeardn).M() + (l1_smeardn + l2_smeardn).M()


                H_pt_norm = (pho1_norm + pho2_norm + l1_norm + l2_norm).Pt()
                H_pt_pho_scale_up = (pho1_scaleup + pho2_scaleup + l1_norm + l2_norm).Pt()
                H_pt_pho_scale_dn = (pho1_scaledn + pho2_scaledn + l1_norm + l2_norm).Pt()
                H_pt_pho_smear_up = (pho1_smearup + pho2_smearup + l1_norm + l2_norm).Pt()
                H_pt_pho_smear_dn = (pho1_smeardn + pho2_smeardn + l1_norm + l2_norm).Pt()
                H_pt_lep_scale_up = (pho1_norm + pho2_norm + l1_scaleup + l2_scaleup).Pt()
                H_pt_lep_scale_dn = (pho1_norm + pho2_norm + l1_scaledn + l2_scaledn).Pt()
                H_pt_lep_smear_up = (pho1_norm + pho2_norm + l1_smearup + l2_smearup).Pt()
                H_pt_lep_smear_dn = (pho1_norm + pho2_norm + l1_smeardn + l2_smeardn).Pt()

                MVA_list = {}
                MVA_list['norm'] = [ntup.pho1Pt, ntup.pho1eta, ntup.pho1phi, ntup.pho1R9, ntup.pho1IetaIeta55 ,ntup.pho2Pt, ntup.pho2eta, ntup.pho2phi, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, ntup.var_Pta, ntup.var_MhMZ, ntup.H_pt ]

                MVA_list['ShowerShape'] = {}
                MVA_list['ShowerShape']['up'] = [ntup.pho1Pt, ntup.pho1eta, ntup.pho1phi, ntup.pho1R9 * (1.0+ntup.pho1ShowerShapeSys), ntup.pho1IetaIeta55 * (1.0+ntup.pho1ShowerShapeSys) ,ntup.pho2Pt, ntup.pho2eta, ntup.pho2phi, ntup.pho2R9 * (1.0+ntup.pho2ShowerShapeSys), ntup.pho2IetaIeta55 * (1.0+ntup.pho2ShowerShapeSys),ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, ntup.var_Pta, ntup.var_MhMZ, ntup.H_pt ]
                MVA_list['ShowerShape']['dn'] = [ntup.pho1Pt, ntup.pho1eta, ntup.pho1phi, ntup.pho1R9 * (1.0-ntup.pho1ShowerShapeSys), ntup.pho1IetaIeta55 * (1.0-ntup.pho1ShowerShapeSys) ,ntup.pho2Pt, ntup.pho2eta, ntup.pho2phi, ntup.pho2R9 * (1.0-ntup.pho2ShowerShapeSys), ntup.pho2IetaIeta55 * (1.0-ntup.pho2ShowerShapeSys),ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, ntup.var_Pta, ntup.var_MhMZ, ntup.H_pt ]

                MVA_list['pho_scale'] = {}
                MVA_list['pho_scale']['up'] = [pho1_scaleup.Pt(), pho1_scaleup.Eta(), pho1_scaleup.Phi(), ntup.pho1R9, ntup.pho1IetaIeta55, pho2_scaleup.Pt(), pho2_scaleup.Eta(), pho2_scaleup.Phi(), ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, Pta_pho_scale_up, MhMZ_pho_scale_up, H_pt_pho_scale_up ]
                MVA_list['pho_scale']['dn'] = [pho1_scaledn.Pt(), pho1_scaledn.Eta(), pho1_scaledn.Phi(), ntup.pho1R9, ntup.pho1IetaIeta55, pho2_scaledn.Pt(), pho2_scaledn.Eta(), pho2_scaledn.Phi(), ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, Pta_pho_scale_dn, MhMZ_pho_scale_dn, H_pt_pho_scale_dn ]
                MVA_list['pho_smear'] = {}
                MVA_list['pho_smear']['up'] = [pho1_smearup.Pt(), pho1_smearup.Eta(), pho1_smearup.Phi(), ntup.pho1R9, ntup.pho1IetaIeta55, pho2_smearup.Pt(), pho2_smearup.Eta(), pho2_smearup.Phi(), ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, Pta_pho_smear_up, MhMZ_pho_smear_up, H_pt_pho_smear_up ]
                MVA_list['pho_smear']['dn'] = [pho1_smeardn.Pt(), pho1_smeardn.Eta(), pho1_smeardn.Phi(), ntup.pho1R9, ntup.pho1IetaIeta55, pho2_smeardn.Pt(), pho2_smeardn.Eta(), pho2_smeardn.Phi(), ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, Pta_pho_smear_dn, MhMZ_pho_smear_dn, H_pt_pho_smear_dn ]

                MVA_list['lep_scale'] = {}
                MVA_list['lep_scale']['up'] = [pho1_norm.Pt(), pho1_norm.Eta(), pho1_norm.Phi(), ntup.pho1R9, ntup.pho1IetaIeta55, pho2_norm.Pt(), pho2_norm.Eta(), pho2_norm.Phi(), ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, Pta_norm, MhMZ_lep_scale_up, H_pt_lep_scale_up ]
                MVA_list['lep_scale']['dn'] = [pho1_norm.Pt(), pho1_norm.Eta(), pho1_norm.Phi(), ntup.pho1R9, ntup.pho1IetaIeta55, pho2_norm.Pt(), pho2_norm.Eta(), pho2_norm.Phi(), ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, Pta_norm, MhMZ_lep_scale_dn, H_pt_lep_scale_dn ]
                MVA_list['lep_smear'] = {}
                MVA_list['lep_smear']['up'] = [pho1_norm.Pt(), pho1_norm.Eta(), pho1_norm.Phi(), ntup.pho1R9, ntup.pho1IetaIeta55, pho2_norm.Pt(), pho2_norm.Eta(), pho2_norm.Phi(), ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, Pta_norm, MhMZ_lep_smear_up, H_pt_lep_smear_up ]
                MVA_list['lep_smear']['dn'] = [pho1_norm.Pt(), pho1_norm.Eta(), pho1_norm.Phi(), ntup.pho1R9, ntup.pho1IetaIeta55, pho2_norm.Pt(), pho2_norm.Eta(), pho2_norm.Phi(), ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, Pta_norm, MhMZ_lep_smear_dn, H_pt_lep_smear_dn ]


                #MVA_value = model.predict_proba(MVA_list)[:, 1]
                #if MVA_value < mvaCut:continue
                #histos['mvaVal'][sample].Fill( MVA_value, weight )

                #mva_v = [MVA_list['norm'], MVA_list['ShowerShape']['up'], MVA_list['ShowerShape']['dn'], MVA_list['pho_scale']['up'], MVA_list['pho_scale']['dn'], MVA_list['pho_smear']['up'], MVA_list['pho_smear']['dn'], MVA_list['lep_scale']['up'], MVA_list['lep_scale']['dn'], MVA_list['lep_smear']['up'], MVA_list['lep_smear']['dn']]
                #results = pool.map(Cal_mvaVal, mva_v)
                #print results
                #pool.close()
                #pool.join()

                '''
                val_norm = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['norm']))
                val_ShowerShape_up = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['ShowerShape']['up']))
                val_ShowerShape_dn = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['ShowerShape']['dn']))
                val_pho_scale_up = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['pho_scale']['up']))
                val_pho_scale_dn = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['pho_scale']['dn']))
                val_pho_smear_up = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['pho_smear']['up']))
                val_pho_smear_dn = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['pho_smear']['dn']))
                val_lep_scale_up = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['lep_scale']['up']))
                val_lep_scale_dn = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['lep_scale']['dn']))
                val_lep_smear_up = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['lep_smear']['up']))
                val_lep_smear_dn = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['lep_smear']['dn']))

                val_ShowerShape_up.start()
                val_ShowerShape_dn.start()
                val_pho_scale_up.start()
                val_pho_scale_dn.start()
                val_pho_smear_up.start()
                val_pho_smear_dn.start()
                val_lep_scale_up.start()
                val_lep_scale_dn.start()
                val_lep_smear_up.start()
                val_lep_smear_dn.start()

                val_ShowerShape_up.join()
                val_ShowerShape_dn.join()
                val_pho_scale_up.join()
                val_pho_scale_dn.join()
                val_pho_smear_up.join()
                val_pho_smear_dn.join()
                val_lep_scale_up.join()
                val_lep_scale_dn.join()
                val_lep_smear_up.join()
                val_lep_smear_dn.join()

                print val_ShowerShape_up #+ 'zzz' + val_ShowerShape_dn + 'zzz' + val_pho_scale_up + 'zzz' + val_pho_scale_dn
                '''

                ['ShowerShape', 'pho_scale', 'pho_smear', 'lep_scale', 'lep_smear']
                if args.sysName == 'ShowerShape':
                    val_ShowerShape_up = model.predict_proba(MVA_list['ShowerShape']['up'])[:, 1]
                    val_ShowerShape_dn = model.predict_proba(MVA_list['ShowerShape']['dn'])[:, 1]

                    histos[sample]['ShowerShape']['up'].Fill(val_ShowerShape_up, weight)
                    histos[sample]['ShowerShape']['dn'].Fill(val_ShowerShape_dn, weight)
                elif args.sysName == 'pho_scale':
                    val_pho_scale_up = model.predict_proba(MVA_list['pho_scale']['up'])[:, 1]
                    val_pho_scale_dn = model.predict_proba(MVA_list['pho_scale']['dn'])[:, 1]

                    histos[sample]['pho_scale']['up'].Fill(val_pho_scale_up , weight)
                    histos[sample]['pho_scale']['dn'].Fill(val_pho_scale_dn , weight)
                elif args.sysName == 'pho_smear':
                    val_pho_smear_up = model.predict_proba(MVA_list['pho_smear']['up'])[:, 1]
                    val_pho_smear_dn = model.predict_proba(MVA_list['pho_smear']['dn'])[:, 1]

                    histos[sample]['pho_smear']['up'].Fill(val_pho_smear_up , weight)
                    histos[sample]['pho_smear']['dn'].Fill(val_pho_smear_dn , weight)
                elif args.sysName == 'lep_scale':
                    val_lep_scale_up = model.predict_proba(MVA_list['lep_scale']['up'])[:, 1]
                    val_lep_scale_dn = model.predict_proba(MVA_list['lep_scale']['dn'])[:, 1]

                    histos[sample]['lep_scale']['up'].Fill(val_lep_scale_up , weight)
                    histos[sample]['lep_scale']['dn'].Fill(val_lep_scale_dn , weight)
                elif args.sysName == 'lep_smear':
                    val_lep_smear_up = model.predict_proba(MVA_list['lep_smear']['up'])[:, 1]
                    val_lep_smear_dn = model.predict_proba(MVA_list['lep_smear']['dn'])[:, 1]

                    histos[sample]['lep_smear']['up'].Fill(val_lep_smear_up , weight)
                    histos[sample]['lep_smear']['dn'].Fill(val_lep_smear_dn , weight)
                else:
                    val_norm = model.predict_proba(MVA_list['norm'])[:, 1]
                    histos[sample]['norm'].Fill(val_norm, weight)








        ## End of for iEvt in range( ntup.GetEntries() )
    ## End of for sample in analyzer_cfg.samp_names






    print '\n\n'
    CountYield(analyzer_cfg, histos, sys_names[0])
    out_file.Close()

    print 'Done'


main()
