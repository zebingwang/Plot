import numpy as np
import matplotlib.pyplot as plt

def main():
    
    mass=['M1','M5','M15','M30']

    sig = {}
    sig['M1'] = [0.0, 26.18, 26.67]
    sig['M5'] = [0.0, 79.40, 79.45]
    sig['M15'] = [0.0, 72.71, 72.80]
    sig['M30'] = [0.0, 69.30, 69.72]

    x = [0.0, 1.0, 2.0]

    for m in sig:
        plt.xlim(xmax=2.1,xmin=-0.1)
        #plt.ylim(ymax=1.0,ymin=0)
        plt.xlabel('number of category')
        plt.ylabel('AMS')
        area = np.pi * 4**2

        plt.plot(x, sig[m], s=area)
        plt.grid()
        plt.savefig('./BiasPlot/sig_nCat_' + m + '.png')
        plt.close('all')


main()