
import os
from Plot_Helper import MakeStack
import Plot_Configs as PC
import Analyzer_Configs as AC
from ROOT import *

#####################################################################
class MuPair:
    def __init__(self, vec, idx1, idx2):
	self.dimu_vec = vec
	self.mass     = vec.M()
	self.pt       = vec.Pt()
	self.mu1_idx  = idx1
	self.mu2_idx  = idx2

def DimuCandidates(ntup):
    muPairs = []
    mu1_vec = TLorentzVector(0,0,0,0)
    mu2_vec = TLorentzVector(0,0,0,0)
    for mu_idx1 in range(ntup.nMuons):
	    mu1_vec.SetPtEtaPhiM(ntup.muon_pt[mu_idx1], ntup.muon_eta[mu_idx1], ntup.muon_phi[mu_idx1], 0.105658367)
	    for mu_idx2 in range(ntup.nMuons):
                if ntup.muon_charge[mu_idx2] != ntup.muon_charge[mu_idx1]:
        		mu2_vec.SetPtEtaPhiM(ntup.muon_pt[mu_idx2], ntup.muon_eta[mu_idx2], ntup.muon_phi[mu_idx2], 0.105658367)

        		muPair = MuPair(mu1_vec + mu2_vec, mu_idx1, mu_idx2)
        	 	muPairs.append(muPair)
    return muPairs



def IsBtagLoose(year, deepCSV):
    if year == '2016' and deepCSV > 0.2217: return True
    if year == '2017' and deepCSV > 0.1522: return True
    if year == '2018' and deepCSV > 0.1241: return True
    return False


def IsBtagMed(year, deepCSV):
    if year == '2016' and deepCSV > 0.6321: return True
    if year == '2017' and deepCSV > 0.4941: return True
    if year == '2018' and deepCSV > 0.4184: return True
    return False

def IsBtagTight(year, deepCSV):
    if year == '2016' and deepCSV > 0.8953: return True
    if year == '2017' and deepCSV > 0.8001: return True
    if year == '2018' and deepCSV > 0.7527: return True
    return False


def CountYield(ana_cfg, histos): # take dimu_mass histos as input
    bin_110 = histos['data'].FindBin(110 + 0.01)
    bin_150 = histos['data'].FindBin(150 - 0.01)
    bin_120 = histos['data'].FindBin(120 + 0.01)
    bin_130 = histos['data'].FindBin(130 - 0.01)
    bin_123 = histos['data'].FindBin(123 + 0.01)
    bin_127 = histos['data'].FindBin(127 - 0.01)

    stacks = MakeStack(histos, ana_cfg, 'dimuon_mass')
    sum_sig = stacks['sig'].GetStack().Last()
    sum_bkg = stacks['bkg'].GetStack().Last()

    print 'Printing event yield: '
    print 'sample \t\t in [110, 150] \t\t in [120, 130] \t\t in [123, 127]'
    for sample in ana_cfg.sig_names + ana_cfg.bkg_names:
	print '%s \t\t %f \t\t %f \t\t %f' %(sample, histos[sample].Integral(bin_110, bin_150), histos[sample].Integral(bin_120, bin_130), histos[sample].Integral(bin_123, bin_127))
    print '================================='
    print 'sum_sig   \t %f \t\t %f \t\t %f' %(sum_sig.Integral(bin_110, bin_150), sum_sig.Integral(bin_120, bin_130), sum_sig.Integral(bin_123, bin_127))
    print 'sum_bkg   \t %f \t\t %f \t\t %f' %(sum_bkg.Integral(bin_110, bin_150), sum_bkg.Integral(bin_120, bin_130), sum_bkg.Integral(bin_123, bin_127))
