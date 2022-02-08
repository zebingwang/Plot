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

from Analyzer_ALP import PIso2D

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
parser.add_argument("--cutVal", dest="cutVal", type=float, default=0.0, help="mva cut value")
parser.add_argument("--mA", dest="mA", default="M5", help="ALP mass")
parser.add_argument('-b', '--blind', dest='blind', action='store_true', default=False, help='Blind signal region')
parser.add_argument('-ln', '--ln', dest='ln', action='store_true', default=False, help='log plot?')
parser.add_argument('--ele', dest='ele', action='store_true', default=False, help='electron channel?')
args = parser.parse_args()

mass = 'param'



gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

mva = args.mva
if args.CR:
    name = mass + '_CR'
elif args.SR:
    name = mass + '_SR'
elif mva:
    name = mass + '_mva'
else:
    name = mass

if args.cut: name = name + '_cut_' + str(args.cutVal) + '_' + args.mA
'''
if args.ele:
    name = name + '_ele'
else:
    name = name + '_mu'
'''

if args.year == '2016':
    file_out = 'plots_16'
    out_name = "ALP_plot_data16_"+name+".root"
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2016.pkl"
    #mvaCut = 0.8675
elif args.year == '2017':
    file_out = 'plots_17'
    out_name = "ALP_plot_data17_"+name+".root"
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2017.pkl"
    #mvaCut = 0.8365
elif args.year == '2018':
    file_out = 'plots_18'
    out_name = "ALP_plot_data18_"+name+".root"
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2018.pkl"
    #mvaCut = 0.9766
elif args.year == 'run2':
    file_out = 'plots_run2'
    out_name = "ALP_plot_run2_"+name+".root"
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_runII.pkl"
    mvaCut = {'M1':0.94, 'M5':0.945, 'M15':0.97, 'M30':0.97}
    #sig_1sigma = {'M1':2.67, 'M5':1.81, 'M15':1.83, 'M30':1.85}
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
    var_names = ['Z_m', 'H_m', 'ALP_m','pho1Pt', 'pho1eta', 'pho1phi', 'pho1R9', 'pho1IetaIeta', 'pho1IetaIeta55','pho1PIso_noCorr' ,'pho2Pt', 'pho2eta', 'pho2phi', 'pho2R9', 'pho2IetaIeta', 'pho2IetaIeta55','pho2PIso_noCorr','ALP_calculatedPhotonIso', 'var_dR_Za', 'var_dR_g1g2', 'var_dR_g1Z', 'var_PtaOverMh', 'var_Pta', 'var_MhMZ', 'H_pt', 'var_PtaOverMa', 'var_MhMa', 'mvaVal_M1', 'mvaVal_M5', 'mvaVal_M15', 'mvaVal_M30']
    var_mva = []
    for ALP_mass in analyzer_cfg.sig_names:
        for r in ['1sigma', '1P5sigma', '2sigma', '3sigma']:
            var_mva.append('mvaVal_'+r+'_'+ALP_mass)
    var_names = var_names + var_mva

    # get mass region
    sigma_low = {}
    sigma_hig = {}
    for sample in analyzer_cfg.sig_names:
        files = TFile(analyzer_cfg.sample_loc + '/ALP_%s.root' %sample)
        filesTree = files.Get("passedEvents")
        filesTree.Draw("{0}>>tree{0}".format("H_m"),"factor*pho1SFs_dR0P15*pho1SFs_dR0P15*({0})".format("H_m>-90&&passChaHadIso&&passNeuHadIso&&passdR_gl&&passHOverE&&H_m>115&&H_m<135"))
        Hist = gDirectory.Get("tree{0}".format("H_m"))
        sigma_bin = 0
        for i in range(Hist.GetNbinsX()/2):
            if Hist.Integral(50-i,51+i)/Hist.Integral()>0.683: 
                sigma_bin = i
                break
        sigma_low[sample] = Hist.GetBinCenter(50-sigma_bin)-125.0-Hist.GetBinWidth(50-sigma_bin)/2.0
        sigma_hig[sample] = Hist.GetBinCenter(51+sigma_bin)-125.0+Hist.GetBinWidth(51+sigma_bin)/2.0

    #### systematic uncertainties ######
    sys_names = ['pho_SFs_Sys_up','pho_SFs_Sys_dn','pho_PUWeight_Sys_up','pho_PUWeight_Sys_dn','lep_dataMC_up','lep_dataMC_dn']
    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)

    ### declare histograms

    histos = {}
    histos_sys = {}

    #histos['mvaVal_bkg']    = TH2F('mvaVal_bkg', 'mvaVal_bkg', 35,  0., 35., 100, 0., 1)#2D

    for var_name in var_names:
        histos[var_name] = {}
        
    for sample in analyzer_cfg.samp_names:
        histos['pho1Pt'][sample]    = TH1F('pho1Pt'    + '_' + sample, 'pho1Pt'    + '_' + sample, 25,  5., 50.)
        histos['pho1eta'][sample]    = TH1F('pho1eta'    + '_' + sample, 'pho1eta'    + '_' + sample, 50,  -3., 3.)
        histos['pho1phi'][sample]    = TH1F('pho1phi'    + '_' + sample, 'pho1phi'    + '_' + sample, 50,  -4., 4.)
        histos['pho1R9'][sample]    = TH1F('pho1R9'    + '_' + sample, 'pho1R9'    + '_' + sample, 50,  0., 1.)
        histos['pho1IetaIeta'][sample]    = TH1F('pho1IetaIeta'    + '_' + sample, 'pho1IetaIeta'    + '_' + sample, 50,  0., 0.07)
        histos['pho1IetaIeta55'][sample]    = TH1F('pho1IetaIeta55'    + '_' + sample, 'pho1IetaIeta55'    + '_' + sample, 50,  0., 0.07)
        histos['pho1PIso_noCorr'][sample]    = TH1F('pho1PIso_noCorr'    + '_' + sample, 'pho1PIso_noCorr'    + '_' + sample, 50, 0., 40.)
        histos['pho2Pt'][sample]    = TH1F('pho2Pt'    + '_' + sample, 'pho2Pt'    + '_' + sample, 25,  5., 30.)
        histos['pho2eta'][sample]    = TH1F('pho2eta'    + '_' + sample, 'pho2eta'    + '_' + sample, 50,  -3., 3.)
        histos['pho2phi'][sample]    = TH1F('pho2phi'    + '_' + sample, 'pho2phi'    + '_' + sample, 50,  -4., 4.)
        histos['pho2R9'][sample]    = TH1F('pho2R9'    + '_' + sample, 'pho2R9'    + '_' + sample, 50,  0., 1.)
        histos['pho2IetaIeta'][sample]    = TH1F('pho2IetaIeta'    + '_' + sample, 'pho2IetaIeta'    + '_' + sample, 50,  0., 0.07)
        histos['pho2IetaIeta55'][sample]    = TH1F('pho2IetaIeta55'    + '_' + sample, 'pho2IetaIeta55'    + '_' + sample, 50,  0., 0.07)
        histos['pho2PIso_noCorr'][sample]    = TH1F('pho2PIso_noCorr'    + '_' + sample, 'pho2PIso_noCorr'    + '_' + sample, 50, 0., 40.)
        histos['Z_m'][sample]    = TH1F('Z_m'    + '_' + sample, 'Z_m'    + '_' + sample, 50,  50., 130.)
        histos['H_m'][sample]    = TH1F('H_m'    + '_' + sample, 'H_m'    + '_' + sample, 28,  110., 180.)
        histos['H_pt'][sample]    = TH1F('H_pt'    + '_' + sample, 'H_pt'    + '_' + sample, 50,  0., 200.)
        histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 25, 0., 40.)
        histos['var_dR_g1g2'][sample] = TH1F('var_dR_g1g2' + '_' + sample, 'var_dR_g1g2' + '_' + sample, 25, 0., 5)
        histos['var_PtaOverMa'][sample] = TH1F('var_PtaOverMa' + '_' + sample, 'var_PtaOverMa' + '_' + sample, 50, 0., 100.)
        histos['var_dR_Za'][sample] = TH1F('var_dR_Za' + '_' + sample, 'var_dR_Za' + '_' + sample, 25, 0., 7.)
        histos['var_dR_g1Z'][sample] = TH1F('var_dR_g1Z' + '_' + sample, 'var_dR_g1Z' + '_' + sample, 25, 0., 8.)
        histos['var_PtaOverMh'][sample] = TH1F('var_PtaOverMh' + '_' + sample, 'var_PtaOverMh' + '_' + sample, 50, 0., 1.)
        histos['var_Pta'][sample] = TH1F('var_Pta' + '_' + sample, 'var_Pta' + '_' + sample, 25, 0., 60.)
        histos['var_MhMa'][sample] = TH1F('var_MhMa' + '_' + sample, 'var_MhMa' + '_' + sample, 25, 100., 200.)
        histos['var_MhMZ'][sample] = TH1F('var_MhMZ' + '_' + sample, 'var_MhMZ' + '_' + sample, 25, 120., 250.)
        histos['ALP_calculatedPhotonIso'][sample] = TH1F('ALP_calculatedPhotonIso' + '_' + sample, 'ALP_calculatedPhotonIso' + '_' + sample, 50, 0., 5.)

        for ALP_mass in analyzer_cfg.sig_names:
            mva_bin = np.hstack((np.arange(0.0, 0.94, 0.02),np.arange(0.94,1.0,0.005)))
            histos['mvaVal_'+ALP_mass][sample]    = TH1F('mvaVal_'+ALP_mass    + '_' + sample, 'mvaVal_'+ALP_mass    + '_' + sample, mva_bin.shape[0]-1,  mva_bin)
            histos['mvaVal_1sigma_'+ALP_mass][sample]    = TH1F('mvaVal_1sigma_'+ALP_mass    + '_' + sample, 'mvaVal_1sigma_'+ALP_mass    + '_' + sample, 200,  0.0, 1.0)
            histos['mvaVal_1P5sigma_'+ALP_mass][sample]    = TH1F('mvaVal_1P5sigma_'+ALP_mass    + '_' + sample, 'mvaVal_1P5sigma_'+ALP_mass    + '_' + sample, 200,  0.0, 1.0)
            histos['mvaVal_2sigma_'+ALP_mass][sample]    = TH1F('mvaVal_2sigma_'+ALP_mass    + '_' + sample, 'mvaVal_2sigma_'+ALP_mass    + '_' + sample, 200,  0.0, 1.0)
            histos['mvaVal_3sigma_'+ALP_mass][sample]    = TH1F('mvaVal_3sigma_'+ALP_mass    + '_' + sample, 'mvaVal_3sigma_'+ALP_mass    + '_' + sample, 200,  0.0, 1.0)



   

    for var_name in var_names:
        histos_sys[var_name] = {}
        for sample in analyzer_cfg.samp_names:
            histos_sys[var_name][sample] = {}
            for sys in sys_names:
                histos_sys[var_name][sample][sys] = copy.deepcopy(histos[var_name][sample])
                histos_sys[var_name][sample][sys].SetNameTitle(var_name+'_'+sample+'_'+sys, var_name+'_'+sample+'_'+sys)

    ### loop over samples and events
    mass_list = {'M1':1.0, 'M5':5.0, 'M15':15.0, 'M30':30.0}
    for sample in analyzer_cfg.samp_names:
        ntup = ntuples[sample] # just a short name
        print '\n\nOn sample: %s' %sample
        print 'total events: %d' %ntup.GetEntries()

        for iEvt in range( ntup.GetEntries() ):
            ntup.GetEvent(iEvt)
            #if (iEvt == 1000): break


            if (iEvt % 100000 == 1):
                print "looking at event %d" %iEvt

            '''
            if args.ele:
                if abs(ntup.l1_id) == 13: 
                    continue
            else:
                if abs(ntup.l1_id) == 11: 
                    continue
            '''

            weight = ntup.factor * ntup.pho1SFs_dR0P15 * ntup.pho1SFs_dR0P15

            if (ntup.H_m > -90):
                #if ntup.dR_pho < 0.02: continue
                #if not ntup.passEleVeto: continue
                if not ntup.passChaHadIso: continue
                if not ntup.passNeuHadIso: continue
                if not ntup.passdR_gl: continue
                if not ntup.passHOverE: continue
                if ntup.H_m>180. or ntup.H_m<110.: continue

                if  args.CR:
                    if ntup.H_m<135. and ntup.H_m>115.: continue
                elif  args.SR:
                    if ntup.H_m>135. or ntup.H_m<115.: continue


                #if abs(ntup.l1_id)!=13: continue

                #print 'pho1: pt:' + str(ntup.pho1Pt) + ', eta: ' + str(ntup.pho1eta) + ', SFs: ' + str(pho1_SFs)
                MVA_value = {}
                if mva:
                    
                    for ALP_mass in mass_list:

                        param = (ntup.ALP_m - mass_list[ALP_mass])/ntup.H_m

                        MVA_list = [ntup.pho1Pt, ntup.pho1eta, ntup.pho1phi, ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr ,ntup.pho2Pt, ntup.pho2eta, ntup.pho2phi, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.pho2PIso_noCorr,ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, ntup.var_PtaOverMh, ntup.H_pt, param]
                        MVA_value[ALP_mass] = model.predict_proba(MVA_list)[:, 1]


                    if args.cut:
                        #if MVA_value[args.mA] < mvaCut[args.mA]: continue
                        if MVA_value[args.mA] < args.cutVal: continue

                    for ALP_mass in mass_list:
                        if args.blind:
                            if sample == 'data' and MVA_value[args.mA] > mvaCut[args.mA]:
                                continue
                            else:
                                histos['mvaVal_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                        else:
                            histos['mvaVal_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            
                            if ntup.H_m<(125.+sigma_hig[ALP_mass]) and ntup.H_m>(125.+sigma_low[ALP_mass]): histos['mvaVal_1sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            if ntup.H_m<(125.+sigma_hig[ALP_mass]*1.5) and ntup.H_m>(125.+sigma_low[ALP_mass]*1.5): histos['mvaVal_1P5sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            if ntup.H_m<(125.+sigma_hig[ALP_mass]*2.) and ntup.H_m>(125.+sigma_low[ALP_mass]*2.): histos['mvaVal_2sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            if ntup.H_m<(125.+sigma_hig[ALP_mass]*3.) and ntup.H_m>(125.+sigma_low[ALP_mass]*3.): histos['mvaVal_3sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            #if sample == 'DYJetsToLL':
                                #histos['mvaVal_bkg'].Fill(mass_list[ALP_mass], MVA_value[ALP_mass], weight)
                    


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

                if args.blind:
                    if not (sample == 'data' and (ntup.H_m<135. and ntup.H_m>115.)): 
                        histos['H_m'][sample].Fill( ntup.H_m, weight )
                else:        
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

                var_map = {'Z_m':ntup.Z_m, 'H_m':ntup.H_m, 'ALP_m':ntup.ALP_m,'pho1Pt':ntup.pho1Pt, 'pho1eta':ntup.pho1eta, 'pho1phi':ntup.pho1phi, 'pho1R9':ntup.pho1R9, 'pho1IetaIeta':ntup.pho1IetaIeta, 'pho1IetaIeta55':ntup.pho1IetaIeta55,'pho1PIso_noCorr':ntup.pho1PIso_noCorr ,'pho2Pt':ntup.pho2Pt, 'pho2eta':ntup.pho2eta, 'pho2phi':ntup.pho2phi, 'pho2R9':ntup.pho2R9, 'pho2IetaIeta':ntup.pho2IetaIeta, 'pho2IetaIeta55':ntup.pho2IetaIeta55,'pho2PIso_noCorr':ntup.pho2PIso_noCorr,'ALP_calculatedPhotonIso':ntup.ALP_calculatedPhotonIso, 'var_dR_Za':ntup.var_dR_Za, 'var_dR_g1g2':ntup.var_dR_g1g2, 'var_dR_g1Z':ntup.var_dR_g1Z, 'var_PtaOverMh':ntup.var_PtaOverMh, 'var_Pta':ntup.var_Pta, 'var_MhMZ':ntup.var_MhMZ, 'H_pt':ntup.H_pt, 'var_PtaOverMa':ntup.var_PtaOverMa, 'var_MhMa':ntup.var_MhMa, 'mvaVal_M1':MVA_value['M1'], 'mvaVal_M5':MVA_value['M5'], 'mvaVal_M15':MVA_value['M15'], 'mvaVal_M30':MVA_value['M30']}
                var_map_mva = {}
                for ALP_mass in analyzer_cfg.sig_names:
                    for r in ['1sigma', '1P5sigma', '2sigma', '3sigma']:
                        var_map_mva['mvaVal_'+r+'_'+ALP_mass] = MVA_value[ALP_mass]
                        
                var_map.update(var_map_mva)

                for sys_name in sys_names:

                    if sys_name =='pho_SFs_Sys_up':
                        #weight = ntup.factor * (ntup.pho1SFs+ntup.pho1SFs_sys) * (ntup.pho2SFs+ntup.pho2SFs_sys)
                        weight_sys = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * (ntup.pho1SFs_dR0P15+ntup.pho1SFs_sys_replaced) * (ntup.pho2SFs_dR0P15+ntup.pho2SFs_sys_replaced)
                    elif sys_name =='pho_SFs_Sys_dn':
                        #weight = ntup.factor * (ntup.pho1SFs-ntup.pho1SFs_sys) * (ntup.pho2SFs-ntup.pho2SFs_sys)
                        weight_sys = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * (ntup.pho1SFs_dR0P15-ntup.pho1SFs_sys_replaced) * (ntup.pho2SFs_dR0P15-ntup.pho2SFs_sys_replaced)
                    elif sys_name =='pho_PUWeight_Sys_up':
                        weight_sys = ntup.event_genWeight * ntup.event_pileupWeightUp * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs_dR0P15 * ntup.pho2SFs_dR0P15
                    elif sys_name =='pho_PUWeight_Sys_dn':
                        weight_sys = ntup.event_genWeight * ntup.event_pileupWeightDn * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs_dR0P15 * ntup.pho2SFs_dR0P15
                    elif sys_name =='lep_dataMC_up':
                        weight_sys = ntup.event_genWeight * ntup.event_pileupWeight * (ntup.l1_dataMC+ntup.l1_dataMCErr) * (ntup.l2_dataMC+ntup.l2_dataMCErr) * ntup.event_weight * ntup.pho1SFs_dR0P15 * ntup.pho2SFs_dR0P15
                    elif sys_name =='lep_dataMC_dn':
                        weight_sys = ntup.event_genWeight * ntup.event_pileupWeight * (ntup.l1_dataMC-ntup.l1_dataMCErr) * (ntup.l2_dataMC-ntup.l2_dataMCErr) * ntup.event_weight * ntup.pho1SFs_dR0P15 * ntup.pho2SFs_dR0P15

                    for var in var_names:
                        #print histos_sys[var].keys()
                        #if var in var_mva: continue
                        histos_sys[var][sample][sys_name].Fill(var_map[var], weight_sys)

                



        ## End of for iEvt in range( ntup.GetEntries() )
    ## End of for sample in analyzer_cfg.samp_names

    ### save raw histograms
    raw_dir = out_file.mkdir('raw_plots')
    raw_dir.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.samp_names:
            plot_cfg.SetHistStyles(histos[var_name][sample], sample)
            histos[var_name][sample].Write()

    #### save sys histograms
    sys_dir = out_file.mkdir('sys_dir')
    sys_dir.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.samp_names:
            for sys in sys_names:
                plot_cfg.SetHistStyles(histos_sys[var_name][sample][sys], sample)
                histos_sys[var_name][sample][sys].Write()

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

        #### uncertainty graph
        total_unc = Total_Unc(stacks['bkg'], histos_sys[var_name], sys_names, analyzer_cfg)

        if args.ln:
            canv_log = CreateCanvas(var_name+'_log')
            DrawOnCanv(canv_log, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label, total_unc, logY=True)

            canv_log.Write()
            SaveCanvPic(canv_log, file_plot, var_name+'_log')
        else:
            canv = CreateCanvas(var_name)
            DrawOnCanv(canv, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label, total_unc, logY=False)

            canv.Write()
            SaveCanvPic(canv, file_plot, var_name)

    '''
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
    '''

    
    #canv2D = CreateCanvas("mvaVal_bkg")
    #plot2D(canv2D, histos["mvaVal_bkg"])
    #canv2D.Write()
    #SaveCanvPic(canv2D, file_plot, "mvaVal_bkg")

    print '\n\n'
    #CountYield(analyzer_cfg, histos['ALP_m'])
    out_file.Close()

    print 'Done'


main()
