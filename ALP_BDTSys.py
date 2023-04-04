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

from Analyzer_ALP import PIso2D, CountYield, CountBDTSys

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
parser.add_argument('--ele', dest='ele', action='store_true', default=False, help='electron channel?')
parser.add_argument('--mu', dest='mu', action='store_true', default=False, help='muon channel?')
args = parser.parse_args()



gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

mva = args.mva

BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/UL/model_ALP_BDT_param.pkl"
mvaCuts = {'M1':0.955, 'M2':0.98, 'M3':0.985, 'M4':0.98, 'M5':0.985, 'M6':0.99, 'M7':0.985, 'M8':0.99, 'M9':0.99, 'M10':0.99, 'M15':0.99, 'M20':0.99, 'M25':0.985, 'M30':0.98}
mass_list = {'M1':1.0, 'M2':2.0, 'M3':3.0, 'M4':4.0, 'M5':5.0, 'M6':6.0, 'M7':7.0, 'M8':8.0, 'M9':9.0, 'M10':10.0, 'M15':15.0, 'M20':20.0, 'M25':25.0, 'M30':30.0}
# load the model from disk
model = pickle.load(open(BDT_filename, 'rb'))

def Cal_mvaVal(*mvalist):
    return model.predict_proba(list(mvalist))[:, 1]


def main():

    analyzer_cfg = AC.Analyzer_Config('inclusive',args.year)
    analyzer_cfg.Print_Config()
    ntuples = LoadNtuples(analyzer_cfg)


    #sys_names = ['norm', 'ShowerShape', 'CMS_scale_g', 'CMS_smear_g', 'CMS_scale_lep', 'CMS_smear_lep']
    sys_names = ['norm', 'ShowerShape', 'CMS_scale_g', 'CMS_smear_g', 'CMS_scale_lep', 'CMS_smear_lep', 'CMS_R9_g', 'CMS_SigmaIEtaIEta_g', 'CMS_PhotonIso_g', 'CMS_ALPIso_g']#bing
    #sys_names = ['norm', 'CMS_R9_g']#bing
    #sys_names = ['norm']
    #sys_names.append(args.sysName)

    MVA_list = {}

    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)

    ### declare histograms

    histos = {}

    for sample in analyzer_cfg.sig_names:
        histos[sample] = {}
        histos[sample]['norm'] = TH1F('mvaVal_norm'    + '_' + sample, 'mvaVal_norm'    + '_' + sample, 50,  -0.1, 1.1)
        for sys_name in sys_names:
            if sys_name == 'norm': continue
            histos[sample][sys_name] = {}
            histos[sample][sys_name]['up'] = TH1F('mvaVal_' + sys_name + '_up'    + '_' + sample, 'mvaVal_' + sys_name + '_up'    + '_' + sample, 100,  -0.1, 1.1)
            histos[sample][sys_name]['dn'] = TH1F('mvaVal_' + sys_name + '_dn'    + '_' + sample, 'mvaVal_' + sys_name + '_dn'    + '_' + sample, 100,  -0.1, 1.1)

    ### loop over samples and events
    for sys_name in sys_names:
        print '##############'
        print 'Systematics : ' +  sys_name
        print '##############'
        for sample in analyzer_cfg.sig_names:
            ntup = ntuples[sample] # just a short name
            print '\n\nOn sample: %s' %sample
            print 'total events: %d' %ntup.GetEntries()

            for iEvt in range( ntup.GetEntries() ):
                ntup.GetEvent(iEvt)
                #if (iEvt == 500): break


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
                    pho1_ShowerShapeup = ROOT.TLorentzVector()
                    pho2_ShowerShapeup = ROOT.TLorentzVector()
                    pho1_ShowerShapedn = ROOT.TLorentzVector()
                    pho2_ShowerShapedn = ROOT.TLorentzVector()
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
                    pho1_f_scaleup = ntup.pho1scaleup/pho1_norm.E()
                    pho1_f_scaledn = ntup.pho1scaledn/pho1_norm.E()
                    pho1_scaleup.SetPxPyPzE(pho1_norm.Px()*pho1_f_scaleup,pho1_norm.Py()*pho1_f_scaleup,pho1_norm.Pz()*pho1_f_scaleup,ntup.pho1scaleup)
                    pho1_scaledn.SetPxPyPzE(pho1_norm.Px()*pho1_f_scaledn,pho1_norm.Py()*pho1_f_scaledn,pho1_norm.Pz()*pho1_f_scaledn,ntup.pho1scaledn)
                    pho1_f_smearup = ntup.pho1smearup/pho1_norm.E()
                    pho1_f_smeardn = ntup.pho1smeardn/pho1_norm.E()
                    pho1_smearup.SetPxPyPzE(pho1_norm.Px()*pho1_f_smearup,pho1_norm.Py()*pho1_f_smearup,pho1_norm.Pz()*pho1_f_smearup,ntup.pho1smearup)
                    pho1_smeardn.SetPxPyPzE(pho1_norm.Px()*pho1_f_smeardn,pho1_norm.Py()*pho1_f_smeardn,pho1_norm.Pz()*pho1_f_smeardn,ntup.pho1smeardn)

                    pho1_ShowerShapeup.SetPxPyPzE(pho1_norm.Px()*(1+ntup.pho1ShowerShapeSys),pho1_norm.Py()*(1+ntup.pho1ShowerShapeSys),pho1_norm.Pz()*(1+ntup.pho1ShowerShapeSys),pho1_norm.E()*(1+ntup.pho1ShowerShapeSys))
                    pho1_ShowerShapedn.SetPxPyPzE(pho1_norm.Px()*(1-ntup.pho1ShowerShapeSys),pho1_norm.Py()*(1-ntup.pho1ShowerShapeSys),pho1_norm.Pz()*(1-ntup.pho1ShowerShapeSys),pho1_norm.E()*(1-ntup.pho1ShowerShapeSys))

                    pho2_norm.SetPtEtaPhiM(ntup.pho2Pt, ntup.pho2eta, ntup.pho2phi, 0.)
                    pho2_f_scaleup = ntup.pho2scaleup/pho2_norm.E()
                    pho2_f_scaledn = ntup.pho2scaledn/pho2_norm.E()
                    pho2_scaleup.SetPxPyPzE(pho2_norm.Px()*pho2_f_scaleup,pho2_norm.Py()*pho2_f_scaleup,pho2_norm.Pz()*pho2_f_scaleup,ntup.pho2scaleup)
                    pho2_scaledn.SetPxPyPzE(pho2_norm.Px()*pho2_f_scaledn,pho2_norm.Py()*pho2_f_scaledn,pho2_norm.Pz()*pho2_f_scaledn,ntup.pho2scaledn)
                    pho2_f_smearup = ntup.pho2smearup/pho2_norm.E()
                    pho2_f_smeardn = ntup.pho2smeardn/pho2_norm.E()
                    pho2_smearup.SetPxPyPzE(pho2_norm.Px()*pho2_f_smearup,pho2_norm.Py()*pho2_f_smearup,pho2_norm.Pz()*pho2_f_smearup,ntup.pho2smearup)
                    pho2_smeardn.SetPxPyPzE(pho2_norm.Px()*pho2_f_smeardn,pho2_norm.Py()*pho2_f_smeardn,pho2_norm.Pz()*pho2_f_smeardn,ntup.pho2smeardn)

                    pho2_ShowerShapeup.SetPxPyPzE(pho2_norm.Px()*(1+ntup.pho2ShowerShapeSys),pho2_norm.Py()*(1+ntup.pho2ShowerShapeSys),pho2_norm.Pz()*(1+ntup.pho2ShowerShapeSys),pho2_norm.E()*(1+ntup.pho2ShowerShapeSys))
                    pho2_ShowerShapedn.SetPxPyPzE(pho2_norm.Px()*(1-ntup.pho2ShowerShapeSys),pho2_norm.Py()*(1-ntup.pho2ShowerShapeSys),pho2_norm.Pz()*(1-ntup.pho2ShowerShapeSys),pho2_norm.E()*(1-ntup.pho2ShowerShapeSys))


                    if abs(ntup.l1_id)==11:
                        l1_norm.SetPtEtaPhiM(ntup.l1_pt, ntup.l1_eta, ntup.l1_phi, ntup.l1_mass)
                        l1_f_scaleup = ntup.l1_scaleup/l1_norm.E()
                        l1_f_scaledn = ntup.l1_scaledn/l1_norm.E()
                        l1_scaleup.SetPxPyPzE(l1_norm.Px()*l1_f_scaleup, l1_norm.Py()*l1_f_scaleup, l1_norm.Pz()*l1_f_scaleup, ntup.l1_scaleup)
                        l1_scaledn.SetPxPyPzE(l1_norm.Px()*l1_f_scaledn, l1_norm.Py()*l1_f_scaledn, l1_norm.Pz()*l1_f_scaledn, ntup.l1_scaledn)
                        l1_f_smearup = ntup.l1_smearup/l1_norm.E()
                        l1_f_smeardn = ntup.l1_smeardn/l1_norm.E()
                        l1_smearup.SetPxPyPzE(l1_norm.Px()*l1_f_smearup, l1_norm.Py()*l1_f_smearup, l1_norm.Pz()*l1_f_smearup, ntup.l1_smearup)
                        l1_smeardn.SetPxPyPzE(l1_norm.Px()*l1_f_smeardn, l1_norm.Py()*l1_f_smeardn, l1_norm.Pz()*l1_f_smeardn, ntup.l1_smeardn)

                        l2_norm.SetPtEtaPhiM(ntup.l2_pt, ntup.l2_eta, ntup.l2_phi, ntup.l2_mass)
                        l2_f_scaleup = ntup.l2_scaleup/l2_norm.E()
                        l2_f_scaledn = ntup.l2_scaledn/l2_norm.E()
                        l2_scaleup.SetPxPyPzE(l2_norm.Px()*l2_f_scaleup, l2_norm.Py()*l2_f_scaleup, l2_norm.Pz()*l2_f_scaleup, ntup.l2_scaleup)
                        l2_scaledn.SetPxPyPzE(l2_norm.Px()*l2_f_scaledn, l2_norm.Py()*l2_f_scaledn, l2_norm.Pz()*l2_f_scaledn, ntup.l2_scaledn)
                        l2_f_smearup = ntup.l2_smearup/l2_norm.E()
                        l2_f_smeardn = ntup.l2_smeardn/l2_norm.E()
                        l2_smearup.SetPxPyPzE(l2_norm.Px()*l2_f_smearup, l2_norm.Py()*l2_f_smearup, l2_norm.Pz()*l2_f_smearup, ntup.l2_smearup)
                        l2_smeardn.SetPxPyPzE(l2_norm.Px()*l2_f_smeardn, l2_norm.Py()*l2_f_smeardn, l2_norm.Pz()*l2_f_smeardn, ntup.l2_smeardn)
                    elif abs(ntup.l1_id)==13:
                        l1_norm.SetPtEtaPhiM(ntup.l1_pt, ntup.l1_eta, ntup.l1_phi, ntup.l1_mass)
                        l1_scaleup.SetPxPyPzE(l1_norm.Px()*ntup.l1_scaleup, l1_norm.Py()*ntup.l1_scaleup, l1_norm.Pz()*ntup.l1_scaleup, l1_norm.E()*ntup.l1_scaleup)
                        l1_scaledn.SetPxPyPzE(l1_norm.Px()*ntup.l1_scaledn, l1_norm.Py()*ntup.l1_scaledn, l1_norm.Pz()*ntup.l1_scaledn, l1_norm.E()*ntup.l1_scaledn)
                        l1_smearup.SetPxPyPzE(l1_norm.Px()*ntup.l1_smearup, l1_norm.Py()*ntup.l1_smearup, l1_norm.Pz()*ntup.l1_smearup, l1_norm.E()*ntup.l1_smearup)
                        l1_smeardn.SetPxPyPzE(l1_norm.Px()*ntup.l1_smeardn, l1_norm.Py()*ntup.l1_smeardn, l1_norm.Pz()*ntup.l1_smeardn, l1_norm.E()*ntup.l1_smeardn)

                        l2_norm.SetPtEtaPhiM(ntup.l2_pt, ntup.l2_eta, ntup.l2_phi, ntup.l2_mass)
                        l2_scaleup.SetPxPyPzE(l2_norm.Px()*ntup.l2_scaleup, l2_norm.Py()*ntup.l2_scaleup, l2_norm.Pz()*ntup.l2_scaleup, l2_norm.E()*ntup.l2_scaleup)
                        l2_scaledn.SetPxPyPzE(l2_norm.Px()*ntup.l2_scaledn, l2_norm.Py()*ntup.l2_scaledn, l2_norm.Pz()*ntup.l2_scaledn, l2_norm.E()*ntup.l2_scaledn)
                        l2_smearup.SetPxPyPzE(l2_norm.Px()*ntup.l2_smearup, l2_norm.Py()*ntup.l2_smearup, l2_norm.Pz()*ntup.l2_smearup, l2_norm.E()*ntup.l2_smearup)
                        l2_smeardn.SetPxPyPzE(l2_norm.Px()*ntup.l2_smeardn, l2_norm.Py()*ntup.l2_smeardn, l2_norm.Pz()*ntup.l2_smeardn, l2_norm.E()*ntup.l2_smeardn)


                    PtaOverMh_norm = (pho1_norm + pho2_norm).Pt()/(pho1_norm + pho2_norm + l1_norm + l2_norm).M()
                    PtaOverMh_pho_scale_up = (pho1_scaleup + pho2_scaleup).Pt()/(pho1_scaleup + pho2_scaleup + l1_norm + l2_norm).M()
                    PtaOverMh_pho_scale_dn = (pho1_scaledn + pho2_scaledn).Pt()/(pho1_scaledn + pho2_scaledn + l1_norm + l2_norm).M()
                    PtaOverMh_pho_smear_up = (pho1_smearup + pho2_smearup).Pt()/(pho1_smearup + pho2_smearup + l1_norm + l2_norm).M()
                    PtaOverMh_pho_smear_dn = (pho1_smeardn + pho2_smeardn).Pt()/(pho1_smeardn + pho2_smeardn + l1_norm + l2_norm).M()
                    PtaOverMh_pho_ShowerShape_up = (pho1_ShowerShapeup + pho2_ShowerShapeup).Pt()/(pho1_ShowerShapeup + pho2_ShowerShapeup + l1_norm + l2_norm).M()
                    PtaOverMh_pho_ShowerShape_dn = (pho1_ShowerShapedn + pho2_ShowerShapedn).Pt()/(pho1_ShowerShapedn + pho2_ShowerShapedn + l1_norm + l2_norm).M()
                    PtaOverMh_lep_scale_up = (pho1_norm + pho2_norm).Pt()/(pho1_norm + pho2_norm + l1_scaleup + l2_scaleup).M()
                    PtaOverMh_lep_scale_dn = (pho1_norm + pho2_norm).Pt()/(pho1_norm + pho2_norm + l1_scaledn + l2_scaledn).M()
                    PtaOverMh_lep_smear_up = (pho1_norm + pho2_norm).Pt()/(pho1_norm + pho2_norm + l1_smearup + l2_smearup).M()
                    PtaOverMh_lep_smear_dn = (pho1_norm + pho2_norm).Pt()/(pho1_norm + pho2_norm + l1_smeardn + l2_smeardn).M()

                    H_pt_norm = (pho1_norm + pho2_norm + l1_norm + l2_norm).Pt()
                    H_pt_pho_scale_up = (pho1_scaleup + pho2_scaleup + l1_norm + l2_norm).Pt()
                    H_pt_pho_scale_dn = (pho1_scaledn + pho2_scaledn + l1_norm + l2_norm).Pt()
                    H_pt_pho_smear_up = (pho1_smearup + pho2_smearup + l1_norm + l2_norm).Pt()
                    H_pt_pho_smear_dn = (pho1_smeardn + pho2_smeardn + l1_norm + l2_norm).Pt()
                    H_pt_pho_ShowerShape_up = (pho1_ShowerShapeup + pho2_ShowerShapeup + l1_norm + l2_norm).Pt()
                    H_pt_pho_ShowerShape_dn = (pho1_ShowerShapedn + pho2_ShowerShapedn + l1_norm + l2_norm).Pt()
                    H_pt_lep_scale_up = (pho1_norm + pho2_norm + l1_scaleup + l2_scaleup).Pt()
                    H_pt_lep_scale_dn = (pho1_norm + pho2_norm + l1_scaledn + l2_scaledn).Pt()
                    H_pt_lep_smear_up = (pho1_norm + pho2_norm + l1_smearup + l2_smearup).Pt()
                    H_pt_lep_smear_dn = (pho1_norm + pho2_norm + l1_smeardn + l2_smeardn).Pt()

                    
                    param_norm = ((pho1_norm + pho2_norm).M()-mass_list[sample])/(pho1_norm + pho2_norm + l1_norm + l2_norm).M()
                    param_pho_scale_up = ((pho1_scaleup + pho2_scaleup).M()-mass_list[sample])/(pho1_scaleup + pho2_scaleup + l1_norm + l2_norm).M()
                    param_pho_scale_dn = ((pho1_scaledn + pho2_scaledn).M()-mass_list[sample])/(pho1_scaledn + pho2_scaledn + l1_norm + l2_norm).M()
                    param_pho_smear_up = ((pho1_smearup + pho2_smearup).M()-mass_list[sample])/(pho1_smearup + pho2_smearup + l1_norm + l2_norm).M()
                    param_pho_smear_dn = ((pho1_smeardn + pho2_smeardn).M()-mass_list[sample])/(pho1_smeardn + pho2_smeardn + l1_norm + l2_norm).M()
                    param_pho_ShowerShape_up = ((pho1_ShowerShapeup + pho2_ShowerShapeup).M()-mass_list[sample])/(pho1_ShowerShapeup + pho2_ShowerShapeup + l1_norm + l2_norm).M()
                    param_pho_ShowerShape_dn = ((pho1_ShowerShapedn + pho2_ShowerShapedn).M()-mass_list[sample])/(pho1_ShowerShapedn + pho2_ShowerShapedn + l1_norm + l2_norm).M()
                    param_lep_scale_up = ((pho1_norm + pho2_norm).M()-mass_list[sample])/(pho1_norm + pho2_norm + l1_scaleup + l2_scaleup).M()
                    param_lep_scale_dn = ((pho1_norm + pho2_norm).M()-mass_list[sample])/(pho1_norm + pho2_norm + l1_scaledn + l2_scaledn).M()
                    param_lep_smear_up = ((pho1_norm + pho2_norm).M()-mass_list[sample])/(pho1_norm + pho2_norm + l1_smearup + l2_smearup).M()
                    param_lep_smear_dn = ((pho1_norm + pho2_norm).M()-mass_list[sample])/(pho1_norm + pho2_norm + l1_smeardn + l2_smeardn).M()

                    MVA_list = {}
                    MVA_list['norm'] = [ntup.pho1Pt, ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr ,ntup.pho2Pt, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.pho2PIso_noCorr,ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, ntup.var_PtaOverMh, ntup.H_pt, param_norm]
                    #MVA_list['norm'] = [ntup.pho1Pt, ntup.pho1R9, ntup.pho1IetaIeta55 ,ntup.pho2Pt,ntup.pho2R9, ntup.pho2IetaIeta55,ntup.ALP_calculatedPhotonIso, ntup.var_dR_g1Z, ntup.var_Pta, ntup.var_MhMZ, ntup.H_pt ]

                    MVA_list['ShowerShape'] = {}
                    MVA_list['ShowerShape']['up'] = [pho1_ShowerShapeup.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_ShowerShapeup.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_pho_ShowerShape_up, H_pt_pho_ShowerShape_up, param_pho_ShowerShape_up]
                    MVA_list['ShowerShape']['dn'] = [pho1_ShowerShapedn.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_ShowerShapedn.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_pho_ShowerShape_dn, H_pt_pho_ShowerShape_dn, param_pho_ShowerShape_dn]

                    MVA_list['CMS_scale_g'] = {}
                    MVA_list['CMS_scale_g']['up'] = [pho1_scaleup.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_scaleup.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_pho_scale_up, H_pt_pho_scale_up, param_pho_scale_up]
                    MVA_list['CMS_scale_g']['dn'] = [pho1_scaledn.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_scaledn.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_pho_scale_dn, H_pt_pho_scale_dn, param_pho_scale_dn]
                    MVA_list['CMS_smear_g'] = {}
                    MVA_list['CMS_smear_g']['up'] = [pho1_smearup.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_smearup.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_pho_smear_up, H_pt_pho_smear_up, param_pho_smear_up]
                    MVA_list['CMS_smear_g']['dn'] = [pho1_smeardn.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_smeardn.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_pho_smear_dn, H_pt_pho_smear_dn, param_pho_smear_dn]

                    MVA_list['CMS_scale_lep'] = {}
                    MVA_list['CMS_scale_lep']['up'] = [pho1_norm.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_norm.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_norm, H_pt_lep_scale_up, param_lep_scale_up]
                    MVA_list['CMS_scale_lep']['dn'] = [pho1_norm.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_norm.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_norm, H_pt_lep_scale_dn, param_lep_scale_dn]
                    MVA_list['CMS_smear_lep'] = {}
                    MVA_list['CMS_smear_lep']['up'] = [pho1_norm.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_norm.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_norm, H_pt_lep_smear_up, param_lep_smear_up]
                    MVA_list['CMS_smear_lep']['dn'] = [pho1_norm.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_norm.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_norm, H_pt_lep_smear_dn, param_lep_smear_dn]

                    scale_value = 0.005
                    MVA_list['CMS_R9_g'] = {}
                    MVA_list['CMS_R9_g']['up'] = [pho1_norm.Pt(), ntup.pho1R9*(1.0+scale_value), ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_norm.Pt(), ntup.pho2R9*(1.0+scale_value), ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_norm, H_pt_lep_smear_up, param_lep_smear_up]
                    MVA_list['CMS_R9_g']['dn'] = [pho1_norm.Pt(), ntup.pho1R9*(1.0-scale_value), ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_norm.Pt(), ntup.pho2R9*(1.0-scale_value), ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_norm, H_pt_lep_smear_dn, param_lep_smear_dn]

                    MVA_list['CMS_SigmaIEtaIEta_g'] = {}
                    MVA_list['CMS_SigmaIEtaIEta_g']['up'] = [pho1_norm.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55*(1.0+0.01), ntup.pho1PIso_noCorr, pho2_norm.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55*(1.0+0.01), ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_norm, H_pt_lep_smear_up, param_lep_smear_up]
                    MVA_list['CMS_SigmaIEtaIEta_g']['dn'] = [pho1_norm.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55*(1.0-0.01), ntup.pho1PIso_noCorr, pho2_norm.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55*(1.0-0.01), ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_norm, H_pt_lep_smear_dn, param_lep_smear_dn]

                    MVA_list['CMS_PhotonIso_g'] = {}
                    MVA_list['CMS_PhotonIso_g']['up'] = [pho1_norm.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr*(1.0+0.01), pho2_norm.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr*(1.0+0.01), ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_norm, H_pt_lep_smear_up, param_lep_smear_up]
                    MVA_list['CMS_PhotonIso_g']['dn'] = [pho1_norm.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr*(1.0-0.01), pho2_norm.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr*(1.0-0.01), ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_norm, H_pt_lep_smear_dn, param_lep_smear_dn]

                    MVA_list['CMS_ALPIso_g'] = {}
                    MVA_list['CMS_ALPIso_g']['up'] = [pho1_norm.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_norm.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso*(1.0+0.05), ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_norm, H_pt_lep_smear_up, param_lep_smear_up]
                    MVA_list['CMS_ALPIso_g']['dn'] = [pho1_norm.Pt(), ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, pho2_norm.Pt(), ntup.pho2R9, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso*(1.0-0.05), ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, PtaOverMh_norm, H_pt_lep_smear_dn, param_lep_smear_dn]


                    #MVA_value = model.predict_proba(MVA_list)[:, 1]
                    #if MVA_value < mvaCut:continue
                    #histos['mvaVal'][sample].Fill( MVA_value, weight )

                    #mva_v = [MVA_list['norm'], MVA_list['ShowerShape']['up'], MVA_list['ShowerShape']['dn'], MVA_list['CMS_scale_g']['up'], MVA_list['CMS_scale_g']['dn'], MVA_list['CMS_smear_g']['up'], MVA_list['CMS_smear_g']['dn'], MVA_list['CMS_scale_lep']['up'], MVA_list['CMS_scale_lep']['dn'], MVA_list['CMS_smear_lep']['up'], MVA_list['CMS_smear_lep']['dn']]
                    #results = pool.map(Cal_mvaVal, mva_v)
                    #print results
                    #pool.close()
                    #pool.join()

                    '''
                    val_norm = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['norm']))
                    val_ShowerShape_up = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['ShowerShape']['up']))
                    val_ShowerShape_dn = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['ShowerShape']['dn']))
                    val_pho_scale_up = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['CMS_scale_g']['up']))
                    val_pho_scale_dn = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['CMS_scale_g']['dn']))
                    val_pho_smear_up = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['CMS_smear_g']['up']))
                    val_pho_smear_dn = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['CMS_smear_g']['dn']))
                    val_lep_scale_up = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['CMS_scale_lep']['up']))
                    val_lep_scale_dn = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['CMS_scale_lep']['dn']))
                    val_lep_smear_up = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['CMS_smear_lep']['up']))
                    val_lep_smear_dn = multiprocessing.Process(target=Cal_mvaVal,args=(MVA_list['CMS_smear_lep']['dn']))

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

                    #['ShowerShape', 'CMS_scale_g', 'CMS_smear_g', 'CMS_scale_lep', 'CMS_smear_lep']
                    if sys_name == 'ShowerShape':
                        val_ShowerShape_up = model.predict_proba(MVA_list['ShowerShape']['up'])[:, 1]
                        val_ShowerShape_dn = model.predict_proba(MVA_list['ShowerShape']['dn'])[:, 1]

                        if val_ShowerShape_up > mvaCuts[sample]:
                            histos[sample]['ShowerShape']['up'].Fill(val_ShowerShape_up, weight)
                        if val_ShowerShape_dn > mvaCuts[sample]:
                            histos[sample]['ShowerShape']['dn'].Fill(val_ShowerShape_dn, weight)
                    elif sys_name == 'CMS_scale_g':
                        val_pho_scale_up = model.predict_proba(MVA_list['CMS_scale_g']['up'])[:, 1]
                        val_pho_scale_dn = model.predict_proba(MVA_list['CMS_scale_g']['dn'])[:, 1]

                        if val_pho_scale_up > mvaCuts[sample]:
                            histos[sample]['CMS_scale_g']['up'].Fill(val_pho_scale_up , weight)
                        if val_pho_scale_dn > mvaCuts[sample]:
                            histos[sample]['CMS_scale_g']['dn'].Fill(val_pho_scale_dn , weight)
                    elif sys_name == 'CMS_smear_g':
                        val_pho_smear_up = model.predict_proba(MVA_list['CMS_smear_g']['up'])[:, 1]
                        val_pho_smear_dn = model.predict_proba(MVA_list['CMS_smear_g']['dn'])[:, 1]

                        if val_pho_smear_up > mvaCuts[sample]:
                            histos[sample]['CMS_smear_g']['up'].Fill(val_pho_smear_up , weight)
                        if val_pho_smear_dn > mvaCuts[sample]:
                            histos[sample]['CMS_smear_g']['dn'].Fill(val_pho_smear_dn , weight)
                    elif sys_name == 'CMS_scale_lep':
                        val_lep_scale_up = model.predict_proba(MVA_list['CMS_scale_lep']['up'])[:, 1]
                        val_lep_scale_dn = model.predict_proba(MVA_list['CMS_scale_lep']['dn'])[:, 1]

                        if val_lep_scale_up > mvaCuts[sample]:
                            histos[sample]['CMS_scale_lep']['up'].Fill(val_lep_scale_up , weight)
                        if val_lep_scale_dn > mvaCuts[sample]:
                            histos[sample]['CMS_scale_lep']['dn'].Fill(val_lep_scale_dn , weight)
                    elif sys_name == 'CMS_smear_lep':
                        val_lep_smear_up = model.predict_proba(MVA_list['CMS_smear_lep']['up'])[:, 1]
                        val_lep_smear_dn = model.predict_proba(MVA_list['CMS_smear_lep']['dn'])[:, 1]

                        if val_lep_smear_up > mvaCuts[sample]:
                            histos[sample]['CMS_smear_lep']['up'].Fill(val_lep_smear_up , weight)
                        if val_lep_smear_dn > mvaCuts[sample]:
                            histos[sample]['CMS_smear_lep']['dn'].Fill(val_lep_smear_dn , weight)
                    elif sys_name == 'CMS_R9_g':
                        val_CMS_R9_g_up = model.predict_proba(MVA_list['CMS_R9_g']['up'])[:, 1]
                        val_CMS_R9_g_dn = model.predict_proba(MVA_list['CMS_R9_g']['dn'])[:, 1]

                        if val_CMS_R9_g_up > mvaCuts[sample]:
                            histos[sample]['CMS_R9_g']['up'].Fill(val_CMS_R9_g_up , weight)
                        if val_CMS_R9_g_dn > mvaCuts[sample]:
                            histos[sample]['CMS_R9_g']['dn'].Fill(val_CMS_R9_g_dn , weight)
                    elif sys_name == 'CMS_SigmaIEtaIEta_g':
                        val_CMS_SigmaIEtaIEta_g_up = model.predict_proba(MVA_list['CMS_SigmaIEtaIEta_g']['up'])[:, 1]
                        val_CMS_SigmaIEtaIEta_g_dn = model.predict_proba(MVA_list['CMS_SigmaIEtaIEta_g']['dn'])[:, 1]

                        if val_CMS_SigmaIEtaIEta_g_up > mvaCuts[sample]:
                            histos[sample]['CMS_SigmaIEtaIEta_g']['up'].Fill(val_CMS_SigmaIEtaIEta_g_up , weight)
                        if val_CMS_SigmaIEtaIEta_g_dn > mvaCuts[sample]:
                            histos[sample]['CMS_SigmaIEtaIEta_g']['dn'].Fill(val_CMS_SigmaIEtaIEta_g_dn , weight)
                    elif sys_name == 'CMS_PhotonIso_g':
                        val_CMS_PhotonIso_g_up = model.predict_proba(MVA_list['CMS_PhotonIso_g']['up'])[:, 1]
                        val_CMS_PhotonIso_g_dn = model.predict_proba(MVA_list['CMS_PhotonIso_g']['dn'])[:, 1]

                        if val_CMS_PhotonIso_g_up > mvaCuts[sample]:
                            histos[sample]['CMS_PhotonIso_g']['up'].Fill(val_CMS_PhotonIso_g_up , weight)
                        if val_CMS_PhotonIso_g_dn > mvaCuts[sample]:
                            histos[sample]['CMS_PhotonIso_g']['dn'].Fill(val_CMS_PhotonIso_g_dn , weight)
                    elif sys_name == 'CMS_ALPIso_g':
                        val_CMS_ALPIso_g_up = model.predict_proba(MVA_list['CMS_ALPIso_g']['up'])[:, 1]
                        val_CMS_ALPIso_g_dn = model.predict_proba(MVA_list['CMS_ALPIso_g']['dn'])[:, 1]

                        if val_CMS_ALPIso_g_up > mvaCuts[sample]:
                            histos[sample]['CMS_ALPIso_g']['up'].Fill(val_CMS_ALPIso_g_up , weight)
                        if val_CMS_ALPIso_g_dn > mvaCuts[sample]:
                            histos[sample]['CMS_ALPIso_g']['dn'].Fill(val_CMS_ALPIso_g_dn , weight)
                    else:
                        val_norm = model.predict_proba(MVA_list['norm'])[:, 1]

                        if val_norm > mvaCuts[sample]:
                            histos[sample]['norm'].Fill(val_norm, weight)








            ## End of for iEvt in range( ntup.GetEntries() )
        ## End of for sample in analyzer_cfg.samp_names






    print '\n\n'
    #CountYield(analyzer_cfg, histos, sys_names[0])
    if args.ele:
        CountBDTSys(analyzer_cfg, histos, sys_names, "ele", args.year)
    elif args.mu:
        CountBDTSys(analyzer_cfg, histos, sys_names, "mu", args.year)
    else:
        CountBDTSys(analyzer_cfg, histos, sys_names, "combine", args.year)

    print 'Done'


main()
