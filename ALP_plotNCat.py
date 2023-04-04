import numpy as np
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

def main():
    
    mass=['M1','M2','M3','M4','M5','M6','M7','M8','M9','M10','M15','M20','M25','M30']

    basic_path = '/publicfs/cms/user/wangzebing/ALP/Plot/optimize_run2UL'



    sig = {}

    for m in mass:
        sig[m] = []
        file = open(basic_path+'/categorize_all_'+m+'.txt', 'r')
        for line in file:
             #sig[m].append(float('%.1f'% float(line.split('--->')[1].split()[0])))
             sig[m].append(format(float(line.split('--->')[1].split()[0]), '.3f'))

    
    x = [1, 2]

    for m in sig:
        plt.xlim(xmax=2.1,xmin=0.9)
        #plt.ylim(ymax=1.0,ymin=0)
        plt.xlabel('number of category ({0:})'.format(m))
        plt.ylabel('AMS')
        area = np.pi * 4**2

        plt.plot(x, sig[m], "o-")
        plt.xticks([1,2],[1,2],rotation=0)
        plt.grid()
        plt.savefig('./BiasPlot_UL/sig_nCat_' + m + '.png')
        plt.close('all')
    

main()