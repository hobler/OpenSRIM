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
    datapath = 'btilt1%s.inp'
    pathletters = ['a','b','c','d']
    subplotindices = [1, 1, 0, 0]
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
    exppath = ''

    fig = plt.figure(figsize= (16.5,11.7)) 
    mainplot = fig.add_subplot(1,2,1)
    subplot  = list()
    subplot.append(fig.add_subplot(2,2,2))
    subplot.append(fig.add_subplot(2,2,4))

    mainplot.set_yscale('log')
    addtolegend('Current Data', mainplot, '-')
    addtolegend('2003 Data', mainplot, '--')
    initShortenLabels(path+datapath%pathletters[-1])
    for pathletter, subplotindex in zip(pathletters, subplotindices):
        plotFromINP(path+datapath%pathletter, mainplot, 1)
        if oldpath:
            plotFromINP(oldpath+datapath%pathletter, mainplot, 1, '', '--' )
            plotFromINP(path+datapath%pathletter, subplot[subplotindex], 1, '', '-', 
                        oldpath+datapath%pathletter)
    PlotOLDData(exppath+'hlb8161.i95w', mainplot, 1, '<100> Si - EXP', '-.', 'k')
    PlotOLDData(exppath+'hlb815a.i95w', mainplot, 1, '', '-.', 'k')
    PlotOLDData(exppath+'hob814a.i95w', mainplot, 1, '', '-.', 'k')
    PlotOLDData(exppath+'hlb813z.i95w', mainplot, 1, '', '-.', 'k')

    subplot[0].axis((0.0, 4000.0, -10, 10))
    subplot[1].axis((0.0, 4000.0, -10, 10))
    mainplot.axis((0.0,4000.0,1e16,1e22))
    mainplot.set_ylabel(r'Concentration [cm$^{-3}$]')
    mainplot.set_xlabel(r'Depth [$\AA$]')
    subplot[0].set_xlabel(r'Depth [$\AA$]')
    subplot[1].set_xlabel(r'Depth [$\AA$]')
    #plt.savefig( sys.argv[0][:-3]+'.png', dpi=150)
    if warn:
        print(warn)
    plt.show()

main()
