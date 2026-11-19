import numpy as np
import matplotlib.pyplot as plt
import os, sys
path = os.path.dirname(os.path.abspath(__file__))
print("path:", path)
sys.path.append(os.path.join(path, '..', '..', 'scripts'))
from read_output import read_his2

sys.path.append(os.path.join(path, '..', 'PLOT'))
from plotMODULE import titleFromINP


def main():
    fig, axs = plt.subplots(1, 3)

    for i, ax in enumerate(axs):
        datapath=f'bspot{i+1}'
        x, z, hist, text = read_his2(os.path.join(path, datapath+'.his2'))
        
        ax.set_aspect('equal')
        
        log10_hmax = np.log10(np.max(hist))
        if i == 0:
            levels = np.linspace(log10_hmax-4, log10_hmax, 9)
        print('levels = ', levels)
        print(x.shape, z.shape, hist.shape)
        ax.contour(x, -z, np.log10(hist), levels, linestyles='solid')
    
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

        x, z, hist, text = read_his2(os.path.join(oldpath, datapath+'.his2'))
    
        #print 'max(hist)=', np.max(hist)
        ax.contour(x, -z, np.log10(hist), levels, linestyles='dashed')
    
        ax.set_ylim(-1000, 0)
        ax.set_xlim(-600, 600)
        titleFromINP(os.path.join(path, datapath+'.inp'), ax)
        if warn:
            print(warn)

    plt.show()
    
main()