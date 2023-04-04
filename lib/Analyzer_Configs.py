import os
import sys
from ROOT import *

class Analyzer_Config:
    def __init__(self, channel, year):
        self.channel    = channel
        self.year       = year
        self.sample_loc = 'NONE'
        self.sig_names  = []
        self.bkg_names  = []
        self.samp_names = []

        self.Config_Analyzer()

    def Config_Analyzer(self):
        if self.channel == 'inclusive' or self.channel == 'ggH' or self.channel == 'VBF' or self.channel == 'WH_3l' or self.channel == 'ZH_4l' or self.channel == 'ttH_lep':

            if self.year == '2016':
                self.sample_loc = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/16'
            elif self.year == '-2016':
                self.sample_loc = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/16APV'
            elif self.year == '2017':
                self.sample_loc = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/17'
            elif self.year == '2018':
                self.sample_loc = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/18'
            elif self.year == 'run2Rereco':
                self.sample_loc = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/Rereco/run2'
            elif self.year == 'run2':
                self.sample_loc = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/run2'
            else:
                print 'do not included at 2016/2017/2018!'
                exit(0)

            if self.year == 'run2Rereco':
                self.sig_names  = ['M1', 'M5', 'M15', 'M30']
            else:    
                self.sig_names  = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M15', 'M20', 'M25', 'M30']

            self.bkg_names  = [ 'DYJetsToLL']
            self.samp_names = self.bkg_names + self.sig_names + ['data']
        else:
            print "channel is invalid: channel = %s" %self.channel
            sys.exit()

    def Print_Config(self):
        print 'Running analysis in channel: %s' %self.channel
        print 'getting ntuples from: %s' %self.sample_loc
        print 'using signals: '
        print self.sig_names
        print 'using backgrounds:'
        print self.bkg_names