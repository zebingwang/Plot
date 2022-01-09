import os
import sys
from ROOT import *

basic_path = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/17/massInde/'
massList = ['M1', 'M5', 'M15', 'M30']
color = {'M1':kRed, 'M5':kYellow, 'M15':kGreen, 'M30':kOrange}

global histos_dR
histos_dR = {}

canv = TCanvas("cc", "cc", 200,10,600,400)
canv.cd()
legend = TLegend(0.7,0.8,0.95,0.92)

for mass in massList:
    files = TFile(basic_path + 'ALP_' + mass + '.root')
    filesTree = files.Get("passedEvents")

    canvas = TCanvas('var_dR_g1g2'+mass,'var_dR_g1g2'+mass,1000,800)
    #Extract the relevant variable from the trees
    filesTree.Draw("{0}>>tree{1}".format('var_dR_g1g2', mass))
    histos_dR[mass] = gDirectory.Get("tree{0}".format(mass))
    #canvas.Clear()
    print histos_dR[mass]
    print histos_dR[mass].Integral()

    
    histos_dR[mass].SetLineColor(color[mass])
    histos_dR[mass].SetLineWidth(2)
    histos_dR[mass].Draw()

    #legend.AddEntry(histos_dR[mass], mass)

#legend.Draw("SAME")
#canv.Update()

canv.SaveAs("test.png")

'''
canv = TCanvas("cc", "cc", 200,10,600,400)
canv.cd()

#hist.SetFillColor(10)
#hist.GetXaxis().SetTitle('mvaVal')
#hist.GetYaxis().SetTitle('Events / (%.2f GeV)' %hist.GetBinWidth(1))
print histos_dR

#histos_dR['M1'].Draw("HIST")

legend = TLegend(0.7,0.8,0.95,0.92)

for mass in massList:

    histos_dR[mass].SetLineColor(color[mass])
    histos_dR[mass].SetLineWidth(2)
    histos_dR[mass].Draw("SAME")

    legend.AddEntry(histos_dR[mass], mass)

legend.Draw("SAME")
canv.Update()

canv.saveAs("test.png")
'''
