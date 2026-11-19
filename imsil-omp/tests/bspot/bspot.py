import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(os.path.abspath(__file__))
print("path:", path)
basename = os.path.basename(path)
sys.path.append(os.path.join(path, '..', 'PLOT'))
path += '/'
from plotMODULE import addtolegend, plotFromTEM


def main():
    datapath = 'bspot.tem'
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

    fig = plt.figure(figsize= (16.5,11.7)) 
    mainplot = list()
    mainplot.append(fig.add_subplot(2,2,1))
    mainplot.append(fig.add_subplot(2,2,3))
    subplot = list()
    subplot.append(fig.add_subplot(2,2,2))
    subplot.append(fig.add_subplot(2,2,4))

    mainplot[0].set_yscale('log')
    mainplot[1].set_yscale('log')
    addtolegend('Current Data',mainplot[0],'-')
    addtolegend('Reference',mainplot[0],'--')
    plotFromTEM(path+datapath, mainplot[0], 2, '', '-' )
    plotFromTEM(path+datapath, mainplot[1], 1, '', '-')
    if oldpath:
        plotFromTEM(oldpath+datapath, mainplot[0], 2, '', '--' )
        plotFromTEM(oldpath+datapath, mainplot[1], 1, '', '--')
        mainplot[1].set_title('')
    
        plotFromTEM(path+datapath, subplot[0], 2, '', '-', oldpath+datapath)
        plotFromTEM(path+datapath, subplot[1], 1, '', '-', oldpath+datapath)

    mainplot[0].axis((0.0, 1000.0, 1e19, 1e23))
    mainplot[1].axis((0.0, 1000.0, 1e19, 1e23))
    mainplot[0].set_ylabel(r'Damage Concentration [cm$^{-3}$]')
    mainplot[1].set_ylabel(r'Ion Concentration [cm$^{-3}$]')
    mainplot[0].set_xlabel(r'Depth [$\AA$]')
    subplot[0].axis((0.0, 1000.0, -10, 10))
    subplot[1].axis((0.0, 1000.0, -10, 10))
    subplot[0].set_xlabel(r'Depth [$\AA$]')
    #plt.savefig( sys.argv[0][:-3]+'.png', dpi=150)
    if warn:
        print(warn)
    plt.show()

main() 