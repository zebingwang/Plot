import os
import sys
from ROOT import *

class Analyzer_Config:
    def __init__(self, channel, mass, year):
        self.channel    = channel
        self.mass       = mass
        self.year       = year
        self.sample_loc = 'NONE'
        self.sig_names  = []
        self.bkg_names  = []
        self.samp_names = []

        self.Config_Analyzer()

    def Config_Analyzer(self):
        if self.channel == 'inclusive' or self.channel == 'ggH' or self.channel == 'VBF' or self.channel == 'WH_3l' or self.channel == 'ZH_4l' or self.channel == 'ttH_lep':

            if self.year == '2016':
                self.sample_loc = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/16/massInde'
            elif self.year == '2017':
                self.sample_loc = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/17/massInde'
            elif self.year == '2018':
                self.sample_loc = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/18/massInde'
            elif self.year == 'run2':
                self.sample_loc = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/run2'
            else:
                print 'do not included at 2016/2017/2018!'
                exit(0)

            self.sig_names  = ['M1', 'M5', 'M15', 'M30']

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