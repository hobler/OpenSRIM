import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(os.path.abspath(__file__))
print("path:", path)
sys.path.append(os.path.join(path, '..', '..', 'scripts'))
from read_output import read_his2

sys.path.append(os.path.join(path, '..', 'PLOT'))
from plotMODULE import titleFromINP


def main():
    datapath='b2d_pbc'
    x, z, hist, text = read_his2(os.path.join(path, datapath+'.his2'))
    
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    
    log10_hmax = np.log10(np.max(hist))
    levels = np.linspace(log10_hmax-4, log10_hmax, 17)
    print('levels = ', levels)
    plt.contour(x, -z, np.log10(hist), levels, linestyles='solid')

    if 'OLD_PATH' in os.environ:
        oldpath = os.environ['OLD_PATH']
        if not os.path.exists(oldpath):
            warn = 'WARNING: Oldpath ' + oldpath + ' does not exist.'
            oldpath = ''
        else:
            warn = ''
    else:
        warn = 'WARNING: Environment variable OLD_PATH not set.'
        oldpath = ''

    if oldpath:
        x, z, hist, text = read_his2(os.path.join(oldpath, datapath+'.his2'))
    
    #print 'max(hist)=', np.max(hist)
    plt.contour(x, -z, np.log10(hist), levels, linestyles='dashed')

    titleFromINP(os.path.join(path, datapath+'.inp'), ax)
    if warn:
        print(warn)
    plt.show()
    
main()
