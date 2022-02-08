import numpy as np
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

def main():
    
    mass=['M1','M5','M15','M30']

    sig = {}
    sig['M1'] = [26.18, 26.67]
    sig['M5'] = [79.40, 79.45]
    sig['M15'] = [72.71, 72.80]
    sig['M30'] = [69.30, 69.72]

    x = [1, 2]

    for m in sig:
        plt.xlim(xmax=2.1,xmin=0.9)
        #plt.ylim(ymax=1.0,ymin=0)
        plt.xlabel('number of category')
        plt.ylabel('AMS')
        area = np.pi * 4**2

        plt.plot(x, sig[m], "o-")
        plt.xticks([1,2],[1,2],rotation=0)
        plt.grid()
        plt.savefig('./BiasPlot/sig_nCat_' + m + '.png')
        plt.close('all')


main()