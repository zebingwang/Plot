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

from xgboost import XGBClassifier
import pickle

mass = 'M15'

if mass == 'M1':
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/model_ALP_M1.pkl"
    mvaCut = 0.49
elif mass == 'M5':
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/model_ALP_M5.pkl"
    mvaCut = 0.858
elif mass == 'M15':
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/model_ALP_M15.pkl"
    mvaCut = 0.8632
else:
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/model_ALP_M30.pkl"
    mvaCut = 0.7482

# load the model from disk
model = pickle.load(open(BDT_filename, 'rb'))

gROOT.SetBatch(True)
gStyle.SetOptStat(0)

def main():

    isEB = 1

    mva = False
    if mva:
        name = mass + '_mva'
    else:
        name = mass

    out_name = "ALP_plot_data17_"+name+".root"

    if not os.path.exists('plots'):
        os.makedirs('plots')
    out_file = TFile( "plots/" + out_name , "RECREATE")

    analyzer_cfg = AC.Analyzer_Config('inclusive', mass)
    analyzer_cfg.Print_Config()
    ntuples = LoadNtuples(analyzer_cfg)


    var_names = ['Z_m', 'H_m', 'ALP_m', 'pho1IetaIeta', 'pho1IetaIeta55', 'pho1PIso', 'var_dR_Za', 'var_dR_g1g2', 'var_dR_g1Z', 'var_PtaOverMa', 'var_PtaOverMh', 'var_Pta', 'var_MhMa', 'var_MhMZ', 'ALP_calculatedPhotonIso']
    plot_cfg = PC.Plot_Config(analyzer_cfg)
    #var_names = ['Z_befor50']

    ### declare histograms

    histos = {}

    histos['phoPIso2D'] = {}
    histos['phoPIso2D']['sig']    = TH2F('PIsoVSPt_sig', 'PIsoVSPt_sig', 100,  -10., 10., 100, 0., 30)#2D
    histos['phoPIso2D']['sig_noCorr']    = TH2F('PIsoVSPt_signoCorr', 'PIsoVSPt_signoCorr', 100,  0., 25., 100, 0., 30)#2D
    histos['phoPIso2D']['data']    = TH2F('PIsoVSPt_data', 'PIsoVSPt_data', 100,  -10., 10., 100, 0., 30)#2D
    histos['phoPIso2D']['data_noCorr']    = TH2F('PIsoVSPt_datanoCorr', 'PIsoVSPt_datanoCorr', 100,  0., 25., 100, 0., 30)#2D

    histos['pho1PIso2D'] = {}
    histos['pho1PIso2D']['sig']    = TH2F('pho1PIsoVSPt_sig', 'pho1PIsoVSPt_sig', 100,  -10., 10., 100, 0., 30)#2D
    histos['pho1PIso2D']['sig_noCorr']    = TH2F('pho1PIsoVSPt_signoCorr', 'pho1PIsoVSPt_signoCorr', 100,  0., 25., 100, 0., 30)#2D
    histos['pho1PIso2D']['data']    = TH2F('pho1PIsoVSPt_data', 'pho1PIsoVSPt_data', 100,  -10., 10., 100, 0., 30)#2D
    histos['pho1PIso2D']['data_noCorr']    = TH2F('pho1PIsoVSPt_datanoCorr', 'pho1PIsoVSPt_datanoCorr', 100,  0., 25., 100, 0., 30)#2D
    histos['pho2PIso2D'] = {}
    histos['pho2PIso2D']['sig']    = TH2F('pho2PIsoVSPt_sig', 'pho2PIsoVSPt_sig', 100,  -10., 10., 100, 0., 30)#2D
    histos['pho2PIso2D']['sig_noCorr']    = TH2F('pho2PIsoVSPt_signoCorr', 'pho2PIsoVSPt_signoCorr', 100,  0., 25., 100, 0., 30)#2D
    histos['pho2PIso2D']['data']    = TH2F('pho2PIsoVSPt_data', 'pho2PIsoVSPt_data', 100,  -10., 10., 100, 0., 30)#2D
    histos['pho2PIso2D']['data_noCorr']    = TH2F('pho2PIsoVSPt_datanoCorr', 'pho2PIsoVSPt_datanoCorr', 100,  0., 25., 100, 0., 30)#2D

    histos['phoIetaIeta'] = {}
    histos['phoIetaIeta']["sig"]    = TH2F('phoIetaIeta_sig', 'phoIetaIeta_sig', 100,  0., 0.06, 100, 0., 0.06)#2D
    histos['phoIetaIeta']["data"]    = TH2F('phoIetaIeta_data', 'phoIetaIeta_data', 100,  0., 0.06, 100, 0., 0.06)#2D
    histos['phoIetaIeta55'] = {}
    histos['phoIetaIeta55']["sig"]    = TH2F('phoIetaIeta55_sig', 'phoIetaIeta55_sig', 100,  0., 0.06, 100, 0., 0.06)#2D
    histos['phoIetaIeta55']["data"]    = TH2F('phoIetaIeta55_data', 'phoIetaIeta55_data', 100,  0., 0.06, 100, 0., 0.06)#2D

    histos['phoHOE'] = {}
    histos['phoHOE']["sig"]    = TH2F('phoHOE_sig', 'phoHOE_sig', 100,  0.005, 0.1, 100, 0.005, 0.1)#2D
    histos['phoHOE']["data"]    = TH2F('phoHOE_data', 'phoHOE_data', 100,  0.005, 0.1, 100, 0.005, 0.1)#2D

    histos['phoPIso'] = {}
    histos['phoPIso']['sig']    = TH2F('phoPIso_sig', 'phoPIso_sig', 100,  -10., 10., 100, -10., 10)#2D
    histos['phoPIso']['data']    = TH2F('phoPIso_data', 'phoPIso_data', 100,  -10., 10., 100, -10., 10)#2D
    histos['phoPIso']['sig_noCorr']    = TH2F('phoPIso_signoCorr', 'phoPIso_signoCorr', 100,  0., 30., 100, 0., 30)#2D
    histos['phoPIso']['data_noCorr']    = TH2F('phoPIso_datanoCorr', 'phoPIso_datanoCorr', 100,  0., 30., 100, 0., 30)#2D

    histos['MhVsMhMZ'] = {}
    histos['MhVsMhMZ']['sig']    = TH2F('MhVsMhMZ_sig', 'MhVsMhMZ_sig', 50,  180., 250., 50, 100., 160)#2D
    histos['MhVsMhMZ']['data']    = TH2F('MhVsMhMZ_data', 'MhVsMhMZ_data', 50,  180., 250., 50, 100., 160)#2D

    for var_name in var_names:
        histos[var_name] = {}
    for sample in analyzer_cfg.samp_names:
        #x = [0.2*i for i in range(15)]
        #x.append(3.0)
        #x.append(5.0)
        #x_bin = np.array(x)
        #histos['Z_befor50'][sample]    = TH1F('Z_m'    + '_' + sample, 'Z_m'    + '_' + sample, 100,  50., 130.)


        histos['Z_m'][sample]    = TH1F('Z_m'    + '_' + sample, 'Z_m'    + '_' + sample, 100,  50., 130.)
        histos['H_m'][sample]    = TH1F('H_m'    + '_' + sample, 'H_m'    + '_' + sample, 30,  100., 180.)
        #histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, len(x)-1, x_bin)
        if mass == 'M1':
            histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 20, 0., 10.)
            histos['var_dR_g1g2'][sample] = TH1F('var_dR_g1g2' + '_' + sample, 'var_dR_g1g2' + '_' + sample, 50, 0., 0.15)
            histos['var_PtaOverMa'][sample] = TH1F('var_PtaOverMa' + '_' + sample, 'var_PtaOverMa' + '_' + sample, 100, 0., 200.)
        elif mass == 'M5':
            histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 20, 0., 10.)
            histos['var_dR_g1g2'][sample] = TH1F('var_dR_g1g2' + '_' + sample, 'var_dR_g1g2' + '_' + sample, 50, 0., 0.55)
            histos['var_PtaOverMa'][sample] = TH1F('var_PtaOverMa' + '_' + sample, 'var_PtaOverMa' + '_' + sample, 50, 0., 30.)
        elif mass == 'M15':
            histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 50, 0., 50.)
            histos['var_dR_g1g2'][sample] = TH1F('var_dR_g1g2' + '_' + sample, 'var_dR_g1g2' + '_' + sample, 50, 0., 2)
            histos['var_PtaOverMa'][sample] = TH1F('var_PtaOverMa' + '_' + sample, 'var_PtaOverMa' + '_' + sample, 50, 0., 12.)
        else:
            histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 50, 0., 70.)
            histos['var_dR_g1g2'][sample] = TH1F('var_dR_g1g2' + '_' + sample, 'var_dR_g1g2' + '_' + sample, 50, 0., 6)
            histos['var_PtaOverMa'][sample] = TH1F('var_PtaOverMa' + '_' + sample, 'var_PtaOverMa' + '_' + sample, 50, 0., 12.)
        #histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 50, 0., 50.)
        histos['var_dR_Za'][sample] = TH1F('var_dR_Za' + '_' + sample, 'var_dR_Za' + '_' + sample, 50, 0., 7.)

        histos['var_dR_g1Z'][sample] = TH1F('var_dR_g1Z' + '_' + sample, 'var_dR_g1Z' + '_' + sample, 50, 0., 8.)

        histos['var_PtaOverMh'][sample] = TH1F('var_PtaOverMh' + '_' + sample, 'var_PtaOverMh' + '_' + sample, 100, 0., 1.)
        histos['var_Pta'][sample] = TH1F('var_Pta' + '_' + sample, 'var_Pta' + '_' + sample, 50, 20., 60.)
        histos['var_MhMa'][sample] = TH1F('var_MhMa' + '_' + sample, 'var_MhMa' + '_' + sample, 50, 50., 400.)
        histos['var_MhMZ'][sample] = TH1F('var_MhMZ' + '_' + sample, 'var_MhMZ' + '_' + sample, 50, 100., 500.)
        histos['ALP_calculatedPhotonIso'][sample] = TH1F('ALP_calculatedPhotonIso' + '_' + sample, 'ALP_calculatedPhotonIso' + '_' + sample, 100, 0., 5.)

        histos['pho1IetaIeta'][sample]    = TH1F('pho1IetaIeta'    + '_' + sample, 'pho1IetaIeta'    + '_' + sample, 100,  0., 0.07)
        histos['pho1IetaIeta55'][sample]    = TH1F('pho1IetaIeta55'    + '_' + sample, 'pho1IetaIeta55'    + '_' + sample, 100,  0., 0.07)
        histos['pho1PIso'][sample]    = TH1F('pho1PIso'    + '_' + sample, 'pho1PIso'    + '_' + sample, 100, 0., 40.)


    ### loop over samples and events
    for sample in analyzer_cfg.samp_names:
        ntup = ntuples[sample] # just a short name
        print '\n\nOn sample: %s' %sample
        print 'total events: %d' %ntup.GetEntries()

        for iEvt in range( ntup.GetEntries() ):
            ntup.GetEvent(iEvt)
            #if (iEvt == 50): break


            if (iEvt % 100000 == 1):
                print "looking at event %d" %iEvt


            weight = ntup.factor
            #weight = ntup.event_weight

            #if (ntup.Z_m > -90):
                #histos['Z_m'][sample].Fill( ntup.Z_m, weight )
                #histos['H_m'][sample].Fill( ntup.H_m, weight )
                #histos['ALP_m'][sample].Fill( ntup.ALP_m, weight )
                #if (ntup.ALP_m > 5.0):
                    #histos['ALP_m'][sample].AddBinContent(len(x)-1,weight)



            if (ntup.Z_dR > -90):



                #if ntup.H_pho_veto>130. or ntup.H_pho_veto<118.: continue


                if mass == 'M1':
                    cutdR_gg = ntup.dR_pho < 0.15 ## 1GeV
                elif mass == 'M5':
                    cutdR_gg = 0.1 < ntup.dR_pho < 0.5 ## 5GeV
                elif mass == 'M15':
                    cutdR_gg = 0.2 < ntup.dR_pho < 2 ## 15 GeV
                else:
                    cutdR_gg = 1. < ntup.dR_pho < 3.2  ## 30 GeV

                if not cutdR_gg: continue
                if not ntup.passHOverE: continue
                #if abs(ntup.l1_id)!=13: continue

                if mva:
                    MVA_list = [ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr, ntup.pho2IetaIeta55, ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_PtaOverMa, ntup.var_PtaOverMh, ntup.var_MhMZ]
                    MVA_value = model.predict_proba(MVA_list)[:, 1]
                    if MVA_value < mvaCut:continue

                #histos['Z_befor50'][sample].Fill( ntup.Z_befor50, weight )
                #if ntup.pho1Pt<15: continue

                histos['H_m'][sample].Fill( ntup.H_pho_veto, weight )
                histos['ALP_m'][sample].Fill( ntup.ALP_pho_veto, weight )
                histos['Z_m'][sample].Fill( ntup.Z_pho_veto, weight )
                histos['var_dR_Za'][sample].Fill( ntup.var_dR_Za, weight )
                histos['var_dR_g1g2'][sample].Fill( ntup.var_dR_g1g2, weight )
                histos['var_dR_g1Z'][sample].Fill( ntup.var_dR_g1Z, weight )
                histos['var_PtaOverMa'][sample].Fill( ntup.var_PtaOverMa, weight )
                histos['var_PtaOverMh'][sample].Fill( ntup.var_PtaOverMh, weight )
                histos['var_Pta'][sample].Fill( ntup.var_Pta, weight )
                histos['var_MhMa'][sample].Fill( ntup.var_MhMa, weight )
                histos['var_MhMZ'][sample].Fill( ntup.var_MhMZ, weight )
                histos['ALP_calculatedPhotonIso'][sample].Fill( ntup.ALP_calculatedPhotonIso, weight )

                histos['pho1IetaIeta'][sample].Fill( ntup.pho1IetaIeta, weight )
                histos['pho1IetaIeta55'][sample].Fill( ntup.pho1IetaIeta55, weight )
                histos['pho1PIso'][sample].Fill( ntup.pho1PIso_noCorr, weight )

                histos['pho1IetaIeta'][sample].Fill( ntup.pho2IetaIeta, weight )
                histos['pho1IetaIeta55'][sample].Fill( ntup.pho2IetaIeta55, weight )
                histos['pho1PIso'][sample].Fill( ntup.pho2PIso_noCorr, weight )

                if sample == 'M1':
                    histos['phoIetaIeta']["sig"].Fill(ntup.pho1IetaIeta, ntup.pho2IetaIeta, weight)
                    histos['phoIetaIeta55']["sig"].Fill(ntup.pho1IetaIeta55, ntup.pho2IetaIeta55, weight)
                    histos['phoHOE']["sig"].Fill(ntup.pho1HOE, ntup.pho2HOE, weight)
                    histos['phoPIso']['sig'].Fill(ntup.pho1PIso, ntup.pho2PIso, weight)
                    histos['phoPIso']['sig_noCorr'].Fill(ntup.pho1PIso_noCorr, ntup.pho2PIso_noCorr, weight)

                    histos['phoPIso2D']['sig'].Fill(ntup.pho1PIso, ntup.pho1Pt, weight)
                    histos['phoPIso2D']['sig_noCorr'].Fill(ntup.pho1PIso_noCorr, ntup.pho1Pt, weight)
                    histos['pho1PIso2D']['sig'].Fill(ntup.pho1PIso, ntup.pho1Pt, weight)
                    histos['pho1PIso2D']['sig_noCorr'].Fill(ntup.pho1PIso_noCorr, ntup.pho1Pt, weight)

                    histos['phoPIso2D']['sig'].Fill(ntup.pho2PIso, ntup.pho2Pt, weight)
                    histos['phoPIso2D']['sig_noCorr'].Fill(ntup.pho2PIso_noCorr, ntup.pho2Pt, weight)
                    histos['pho2PIso2D']['sig'].Fill(ntup.pho2PIso, ntup.pho2Pt, weight)
                    histos['pho2PIso2D']['sig_noCorr'].Fill(ntup.pho2PIso_noCorr, ntup.pho2Pt, weight)

                    histos['MhVsMhMZ']['sig'].Fill(ntup.var_MhMZ, ntup.H_HOE, weight)
                if sample == 'data':
                    histos['phoIetaIeta']["data"].Fill(ntup.pho1IetaIeta, ntup.pho2IetaIeta, weight)
                    histos['phoIetaIeta55']["data"].Fill(ntup.pho1IetaIeta55, ntup.pho2IetaIeta55, weight)
                    histos['phoHOE']["data"].Fill(ntup.pho1HOE, ntup.pho2HOE, weight)
                    histos['phoPIso']['data'].Fill(ntup.pho1PIso, ntup.pho2PIso, weight)
                    histos['phoPIso']['data_noCorr'].Fill(ntup.pho1PIso_noCorr, ntup.pho2PIso_noCorr, weight)

                    histos['phoPIso2D']['data'].Fill(ntup.pho1PIso, ntup.pho1Pt, weight)
                    histos['phoPIso2D']['data_noCorr'].Fill(ntup.pho1PIso_noCorr, ntup.pho1Pt, weight)
                    histos['pho1PIso2D']['data'].Fill(ntup.pho1PIso, ntup.pho1Pt, weight)
                    histos['pho1PIso2D']['data_noCorr'].Fill(ntup.pho1PIso_noCorr, ntup.pho1Pt, weight)

                    histos['phoPIso2D']['data'].Fill(ntup.pho2PIso, ntup.pho2Pt, weight)
                    histos['phoPIso2D']['data_noCorr'].Fill(ntup.pho2PIso_noCorr, ntup.pho2Pt, weight)
                    histos['pho2PIso2D']['data'].Fill(ntup.pho2PIso, ntup.pho2Pt, weight)
                    histos['pho2PIso2D']['data_noCorr'].Fill(ntup.pho2PIso_noCorr, ntup.pho2Pt, weight)

                    histos['MhVsMhMZ']['data'].Fill(ntup.var_MhMZ, ntup.H_HOE, weight)






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

    for var_name in var_names:
        stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
        scaled_sig = 0
        #ratio_plot = 0
        scaled_sig = ScaleSignal(plot_cfg, stacks['sig'], histos[var_name], var_name)
        ratio_plot = MakeRatioPlot(histos[var_name]['data'], stacks['all'].GetStack().Last(), var_name)
        legend = MakeLegend(plot_cfg, histos[var_name], scaled_sig)

        canv = CreateCanvas(var_name)
        DrawOnCanv(canv, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label, logY=False)

        canv.Write()
        SaveCanvPic(canv, "plots/plot_"+name, var_name)

    for var_name in var_names:
        stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
        scaled_sig = 0
        #ratio_plot = 0
        scaled_sig = ScaleSignal(plot_cfg, stacks['sig'], histos[var_name], var_name)
        ratio_plot = MakeRatioPlot(histos[var_name]['data'], stacks['all'].GetStack().Last(), var_name)
        legend = MakeLegend(plot_cfg, histos[var_name], scaled_sig)

        canv = CreateCanvas(var_name)
        DrawOnCanv(canv, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label, logY=True)

        canv.Write()
        SaveCanvPic(canv, "plots/plot_"+name, var_name+'_log')



    '''
    canv2D = CreateCanvas("MhVsMhMZ")
    plot2D(canv2D, histos['MhVsMhMZ']['sig'], histos['MhVsMhMZ']['data'])
    canv2D.Write()

    canv2D_CONT = CreateCanvas("MhVsMhMZ_CONT")
    plot2D_CONT(canv2D_CONT, histos['MhVsMhMZ']['sig'], histos['MhVsMhMZ']['data'])
    canv2D_CONT.Write()


    canv2D_noCorr = CreateCanvas("PIso2D_noCorr")
    PIso2D(canv2D_noCorr, histos['phoPIso2D']['sig_noCorr'], histos['phoPIso2D']['data_noCorr'], isEB)
    canv2D_noCorr.Write()

    canv2D_pho1 = CreateCanvas("pho1PIso2D")
    PIso2D(canv2D_pho1, histos['pho1PIso2D']['sig'], histos['pho1PIso2D']['data'], isEB)
    canv2D_pho1.Write()

    canv2D_pho2 = CreateCanvas("pho2PIso2D")
    PIso2D(canv2D_pho2, histos['pho2PIso2D']['sig'], histos['pho2PIso2D']['data'], isEB)
    canv2D_pho2.Write()

    canv2D_pho1_noCorr = CreateCanvas("pho1PIso2D_noCorr")
    PIso2D(canv2D_pho1_noCorr, histos['pho1PIso2D']['sig_noCorr'], histos['pho1PIso2D']['data_noCorr'], isEB)
    canv2D_pho1_noCorr.Write()

    canv2D_pho2_noCorr = CreateCanvas("pho2PIso2D_noCorr")
    PIso2D(canv2D_pho2_noCorr, histos['pho2PIso2D']['sig_noCorr'], histos['pho2PIso2D']['data_noCorr'], isEB)
    canv2D_pho2_noCorr.Write()

    canv2D_IeIe = CreateCanvas("IeIe")
    plot2D(canv2D_IeIe, histos['phoIetaIeta']["sig"], histos['phoIetaIeta']["data"])
    canv2D_IeIe.Write()

    canv2D_IeIe55 = CreateCanvas("IeIe55")
    plot2D(canv2D_IeIe55, histos['phoIetaIeta55']["sig"], histos['phoIetaIeta55']["data"])
    canv2D_IeIe55.Write()

    canv2D_HOE = CreateCanvas("HOE")
    plot2D(canv2D_HOE, histos['phoHOE']["sig"], histos['phoHOE']["data"])
    canv2D_HOE.Write()

    #canv2D_PIso = CreateCanvas("PIso")
    #plot2D(canv2D_PIso, histos['phoPIso']['sig'], histos['phoPIso']['data'])
    #canv2D_PIso.Write()

    canv2D_PIso_noCorr = CreateCanvas("PIso_noCorr")
    plot2D(canv2D_PIso_noCorr, histos['phoPIso']['sig_noCorr'], histos['phoPIso']['data_noCorr'])
    #canv2D_PIso_noCorr.Write()
    '''


    print '\n\n'
    #CountYield(analyzer_cfg, histos['ALP_m'])
    out_file.Close()

    print 'Done'


main()
