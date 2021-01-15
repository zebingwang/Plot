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

from Analyzer_ALP import CountCutFlow, CountCutFlow_less, CountCutFlow_mva, CountCutFlow_mva_less

gROOT.SetBatch(True)

def main():
    out_name = "ALP_plot_data17.root"
    if not os.path.exists('plots'):
        os.makedirs('plots')
    out_file = TFile( "plots/" + out_name , "RECREATE")

    analyzer_cfg = AC.Analyzer_Config('inclusive')
    analyzer_cfg.Print_Config()
    ntuples = LoadNtuples(analyzer_cfg)


    var_names = ['Z_m', 'H_m', 'ALP_m', 'pho1IetaIeta', 'pho1IetaIeta55', 'pho1PIso']
    plot_cfg = PC.Plot_Config(analyzer_cfg)

    ### declare histograms
    histos = {}
    for var_name in var_names:
        histos[var_name] = {}
    for sample in analyzer_cfg.samp_names:
	x = [0.2*i for i in range(15)]
	x.append(3.0)
	x.append(5.0)
	x.append(11.0)
	x.append(13.0)
	x.append(15.0)
	x_bin = np.array(x)
    	histos['Z_m'][sample]    = TH1F('Z_m'    + '_' + sample, 'Z_m'    + '_' + sample, 100,  50., 130.)
    	histos['H_m'][sample]    = TH1F('H_m'    + '_' + sample, 'H_m'    + '_' + sample, 100,  100., 180.)
    	histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, len(x)-1, x_bin)

    ### loop over samples and events
    for sample in analyzer_cfg.samp_names:
	ntup = ntuples[sample] # just a short name
	print '\n\nOn sample: %s' %sample
        print 'total events: %d' %ntup.GetEntries()

        for iEvt in range( ntup.GetEntries() ):
    	    ntup.GetEvent(iEvt)
    	    #if (iEvt % 50 != 1): continue

	    if (ntup.Z_m < 0): continue

    	    if (iEvt % 100000 == 1):
    	        print "looking at event %d" %iEvt

	    if sample == 'data':
		weight = ntup.event_weight
	    else:
		weight = ntup.event_weight
	    if sample == 'DYJetsToLL':
		if ntup.pho1_matche_PdgId == 22 and ntup.pho1_matchedR < 0.3: continue
		if ntup.pho2_matche_PdgId == 22 and ntup.pho2_matchedR < 0.3: continue
    	    histos['Z_m'][sample] .Fill( ntup.Z_m, weight )
    	    histos['H_m'][sample] .Fill( ntup.H_m, weight )
    	    histos['ALP_m'][sample].Fill( ntup.ALP_m, weight )
	    if (ntup.ALP_m > 15.0):
		histos['ALP_m'][sample].AddBinContent(len(x)-1,weight)
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
        #scaled_sig = 0
        #ratio_plot = 0
    	scaled_sig = ScaleSignal(plot_cfg, stacks['sig'], var_name)
    	ratio_plot = MakeRatioPlot(histos[var_name]['data'], stacks['all'].GetStack().Last(), var_name)
    	legend = MakeLegend(plot_cfg, histos[var_name], scaled_sig)

     	canv = CreateCanvas(var_name)
    	DrawOnCanv(canv, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label)
    	canv.Write()

    print '\n\n'
    #CountYield(analyzer_cfg, histos['ALP_m'])
    out_file.Close()

    CountCutFlow(analyzer_cfg, plot_cfg.lumi, 'plots/')
    #CountCutFlow_less(analyzer_cfg, plot_cfg.lumi, 'plots/')
    #CountCutFlow_mva(analyzer_cfg, plot_cfg.lumi, 'plots/')
    #CountCutFlow_mva_less(analyzer_cfg, plot_cfg.lumi, 'plots/')

main()
