import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(os.path.abspath(__file__))
print("path:", path)
basename = os.path.basename(path)
sys.path.append(os.path.join(path, '..', 'PLOT'))
path += '/'
from plotMODULE import addtolegend, plotFromTEM, plotFromINP


def main():
    datapath = 'krdyn.tem'
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
    fig = plt.figure(figsize= (8.3,11.7)) 
    mainplot = fig.add_subplot(2,1,1)
    subplot  = fig.add_subplot(2,1,2)

    mainplot.set_yscale('log')
    addtolegend('Current Data',mainplot,'-')
    addtolegend('Reference',mainplot,'--')

    plotFromTEM(path+datapath, mainplot, 1, ext='cell', record=2)
    if oldpath:
        plotFromTEM(oldpath+datapath, mainplot, 1, '', '--', ext='cell', record=2)
        plotFromTEM(path+datapath, subplot, 1, '', '-', oldpath+datapath,
                    record=2, ext='cell')
    plotFromINP(path+'krdyn_10e15_stat.inp', mainplot, 1, '', '-.') 

    subplot.axis((0.0, 200.0, -10, 10))
    mainplot.axis((0.0, 200.0, 1e19, 1e23))
    mainplot.set_ylabel(r'Ion Distribution [cm$^{-3}$]')
    mainplot.set_xlabel(r'Depth [$\AA$]')
    subplot.set_xlabel(r'Depth [$\AA$]')
    #plt.savefig( sys.argv[0][:-3]+'.png', dpi=150)
    if warn:
        print(warn)
    plt.show()

main()
