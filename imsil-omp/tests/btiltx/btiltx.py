import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(os.path.abspath(__file__))
print("path:", path)
basename = os.path.basename(path)
sys.path.append(os.path.join(path, '..', 'PLOT'))
path += '/'
from plotMODULE import addtolegend, initShortenLabels, plotFromINP, PlotOLDData


def main():
    datapath = 'btiltx%s.inp'
    pathletters = ['1', '2']
    subplotindices = [1, 1, 0, 0, 0]
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

    fig = plt.figure(figsize=(16.5,11.7)) 
    mainplot = fig.add_subplot(1,2,1)
    subplot = list()
    subplot.append(fig.add_subplot(2,2,2))
    subplot.append(fig.add_subplot(2,2,4))

    mainplot.set_yscale('log')
    addtolegend('Current Data', mainplot, '-')
    addtolegend('Reference', mainplot, '--')
    initShortenLabels(path+datapath%pathletters[-1])
    for pathletter, subplotindex in zip(pathletters, subplotindices):
        plotFromINP(path+datapath%pathletter, mainplot, 1)
        if oldpath:
            plotFromINP(oldpath+datapath%pathletter, mainplot, 1, '', '--')
            plotFromINP(path+datapath%pathletter, subplot[subplotindex], 1, '', '-', 
                        oldpath+datapath%pathletter)

    mainplot.legend()
    subplot[0].axis((0.0, 5000.0, -10.0, 10.0))
    subplot[1].axis((0.0, 5000.0, -10.0, 10.0))
    mainplot.axis((0.0,5000.0,1e15,1e19))
    mainplot.set_ylabel(r'Concentration [cm$^{-3}$]')
    mainplot.set_xlabel(r'Depth [$\AA$]')
    subplot[0].set_xlabel(r'Depth [$\AA$]')
    subplot[1].set_xlabel(r'Depth [$\AA$]')
    #plt.savefig( sys.argv[0][:-3]+'.png', dpi=150)
    if warn:
        print(warn)
    plt.show()

main()
