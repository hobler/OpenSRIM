import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(os.path.abspath(__file__))
print("path:", path)
basename = os.path.basename(path)
sys.path.append(os.path.join(path, '..', 'PLOT'))
path += '/'
from plotMODULE import addtolegend, initShortenLabels, plotFromINP


def main():
    datapath = 'multilayer%s.inp'
    pathletters = ['1d1', '1d2', '1d3', '2d1', '2d2']
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
    initShortenLabels(path+datapath%pathletters[-1])
    for pathletter in pathletters:
        plotFromINP(path+datapath%pathletter, mainplot, 1, '', '-' )
        if oldpath:
            plotFromINP(oldpath+datapath%pathletter, mainplot, 1, '', '--' )
            plotFromINP(path+datapath%pathletter, subplot, 1, '', '-', 
                        oldpath+datapath%pathletter)

    subplot.axis((0.0, 3000.0, -10.0, 10.0))
    mainplot.axis((0.0, 3000.0, 1e16, 1e20))
    mainplot.set_ylabel(r'Concentration [cm$^{-3}$]')
    mainplot.set_xlabel(r'Depth [$\AA$]')
    subplot.set_xlabel(r'Depth [$\AA$]')
    #plt.savefig( sys.argv[0][:-3]+'.png', dpi=150)
    if warn:
        print(warn)
    plt.show()

main()
