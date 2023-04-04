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

import CMS_lumi, tdrstyle

from xgboost import XGBClassifier
import pickle
import copy 
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy import interpolate



def main():

    mass_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
    mass_interp = range(1,31)
    years = ['16', '16APV', '17', '18']
    
    path_basic_num = '/publicfs/cms/user/wangzebing/ALP/Analysis_code/fit/fit_run2_UL/'
    path_basic_dem = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/'
    eff = {}
    # v1
    #resolution = {'16':{'ele':[2.95, 2.47, 2.12, 2.03, 2.36, 2.27, 2.30, 2.49, 2.22, 2.47, 2.40, 2.69, 2.27, 2.11], 
    #                     'mu':[2.39, 1.45, 1.58, 1.55, 1.50, 1.61, 1.51, 1.55, 1.48, 1.51, 1.40, 1.43, 1.55, 1.43]}, 
    #           '16APV':{'ele':[3.77, 2.37, 2.32, 2.32, 2.20, 2.26, 2.29, 2.19, 2.13, 2.11, 2.11, 2.12, 2.01, 2.20], 
    #                     'mu':[2.49, 1.58, 1.69, 1.54, 1.76, 1.76, 1.57, 1.69, 1.65, 1.47, 1.41, 1.37, 1.44, 1.40]}, 
    #              '17':{'ele':[3.54, 2.38, 2.36, 2.40, 2.36, 2.37, 2.46, 2.21, 2.49, 2.33, 2.15, 2.25, 2.11, 2.34], 
    #                     'mu':[2.60, 1.45, 1.59, 1.64, 1.72, 1.42, 1.52, 1.46, 1.64, 1.59, 1.61, 1.39, 1.38, 1.46]}, 
    #              '18':{'ele':[3.89, 2.31, 2.14, 2.46, 2.29, 2.22, 2.45, 2.21, 2.35, 2.20, 2.10, 2.21, 2.36, 2.15], 
    #                     'mu':[3.22, 1.64, 1.51, 1.48, 1.68, 1.45, 1.68, 1.56, 1.69, 1.52, 1.38, 1.40, 1.34, 1.46]}}
    #v2
    resolution = {'16':{'ele':[2.90, 2.41, 2.40, 2.16, 2.25, 2.23, 2.35, 2.47, 2.18, 2.43, 2.23, 2.33, 2.12, 2.27], 
                         'mu':[2.28, 1.52, 1.64, 1.52, 1.62, 1.50, 1.49, 1.55, 1.67, 1.46, 1.40, 1.46, 1.55, 1.41]}, 
               '16APV':{'ele':[3.38, 2.43, 2.32, 2.36, 2.12, 2.49, 2.38, 2.14, 2.21, 2.00, 1.97, 2.03, 2.22, 2.23], 
                         'mu':[2.61, 1.66, 1.48, 1.56, 1.64, 1.61, 1.69, 1.63, 1.54, 1.45, 1.52, 1.53, 1.40, 1.50]}, 
                  '17':{'ele':[4.22, 2.46, 2.38, 2.18, 2.53, 2.53, 2.46, 2.33, 2.34, 2.20, 2.24, 2.17, 1.93, 2.34], 
                         'mu':[2.09, 1.45, 1.59, 1.61, 1.73, 1.56, 1.50, 1.57, 1.54, 1.56, 1.57, 1.58, 1.47, 1.50]}, 
                  '18':{'ele':[3.26, 2.32, 2.31, 2.46, 2.37, 2.35, 2.42, 2.22, 2.42, 2.33, 2.38, 2.29, 2.25, 2.05], 
                         'mu':[3.07, 1.69, 1.45, 1.48, 1.72, 1.49, 1.68, 1.70, 1.69, 1.62, 1.52, 1.50, 1.46, 1.37]}}
    

    for year in years:
        eff[year] = {}
        for channel in ['ele','mu']:
            eff[year][channel] = []
            for mass in mass_list:
                file_name_num = path_basic_num + 'ALP_data_sig_Am' + str(mass) + '_Hm125_' + year + '_workspace_'+channel+'.root'
                file_name_dem = path_basic_dem + year + '/' + 'ALP_M' + str(mass) + '.root'

                file_num = TFile(file_name_num)
                file_dem = TFile(file_name_dem)
                
                #print file_name
                num_final = file_num.CMS_hza_workspace.data("ggh_125_13TeV_cat0").sumEntries()
                dem_final = file_dem.nEvents_weight_ntuple.GetBinContent(1) + file_dem.nEvents_weight_ntuple.GetBinContent(2)

                print 'year: ', year, 'channel: ', channel, 'mass: ', mass, 'dem: ', dem_final, 'num: ', num_final

                eff[year][channel].append(num_final/dem_final)

            #print "{0} ALP mass:{1}, num_final:{2}, eff:{3}".format(year, mass, num_final, eff)
    f = open('./interpolation/eff.txt','a')

    for channel in ['ele','mu']:
        colors = ['red','orange','blue', 'green']
        plt.xlim(xmax=35.0,xmin=0.0)
        #plt.ylim(ymax=1.0,ymin=0)
        plt.xlabel('m(a) GeV')
        plt.ylabel('$Efficiency\\times Acceptance$')

        for year in years:
            plt.plot(mass_list, eff[year][channel], "o-", c=colors[years.index(year)], label=year)
        
        plt.grid()
        plt.legend()
        plt.savefig('./interpolation/eff_'+channel+'.png')
        plt.close('all')

        ### resolution
        plt.xlim(xmax=35.0,xmin=0.0)
        #plt.ylim(ymax=1.0,ymin=0)
        plt.xlabel('m(a) GeV')
        plt.ylabel('$Signal Resolution$')
    
        for year in years:
            plt.plot(mass_list, resolution[year][channel], "o-", c=colors[years.index(year)], label=year)
     
        plt.grid()
        plt.legend()
        plt.savefig('./interpolation/resolution_'+channel+'.png')
        plt.close('all')

        ### resolution interpolation
        plt.xlim(xmax=35.0,xmin=0.0)
        #plt.ylim(ymax=1.0,ymin=0)
        plt.xlabel('m(a) GeV')
        plt.ylabel('$Efficiency\\times Acceptance$')
    
        for year in years:
            f_interp = interpolate.interp1d(mass_list, eff[year][channel], kind='cubic')
            eff_interp = f_interp(mass_interp)
            plt.plot(mass_interp, eff_interp, "o-", c=colors[years.index(year)], label=year+'_interp')

            rate=[]
            for m_interp in mass_interp:
                if m_interp in mass_list:
                    rate.append(1.0)
                else:
                    for m_nom in mass_list:
                        if m_nom<10: continue
                        if abs(m_interp - m_nom) <= 2:
                            rate.append(1.0+(eff_interp[m_interp-1] - eff_interp[m_nom-1])/eff_interp[m_nom-1])
                            break
                        else:
                            continue
            #f.write('year:'+year+' channel:'+channel+' '+str(rate)+'\n')
            f.write(year+'_'+channel+' '+str(rate)+'\n')


        plt.grid()
        plt.legend()
        plt.savefig('./interpolation/eff_'+channel+'_interp.png')
        plt.close('all')

        


    print '\n\n'

    print 'Done'


main()
