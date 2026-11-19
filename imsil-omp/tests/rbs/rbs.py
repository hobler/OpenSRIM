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
    datapath = 'rbs.tem'
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

    fig = plt.figure(figsize=(8.3,11.7)) 
    mainplot = fig.add_subplot(2,1,1)
    subplot  = fig.add_subplot(2,1,2)   

    mainplot.set_yscale('linear')
    addtolegend('Current Data', subplot, '-')
    addtolegend('Reference', subplot, '--')
    plotFromTEM(path+datapath, mainplot, 1, ext='rbs')
    plotFromINP('rbs5.inp', subplot, 1, ext='rbs')
    if oldpath:
        plotFromTEM(oldpath+datapath, mainplot, 1, '', '--', ext='rbs')
        plotFromINP(oldpath+'rbs5.inp', subplot, 1, '', '--', ext='rbs')
        subplot.set_title('')

    mainplot.axis((1.35e6, 1.72e6, 0., 40.))
    subplot.axis((1.35e6, 1.72e6, 0., 200.))
    mainplot.set_ylabel('RBS signal [a.u.]')
    mainplot.set_xlabel('Energy [eV]')
    subplot.set_xlabel('Energy [eV]')
    #plt.savefig( sys.argv[0][:-3]+'.png', dpi=150)
    if warn:
        print(warn)
    plt.show()

main()
