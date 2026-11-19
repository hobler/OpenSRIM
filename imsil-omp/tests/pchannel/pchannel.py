import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(os.path.abspath(__file__))
print("path:", path)
basename = os.path.basename(path)
sys.path.append(os.path.join(path, '..', 'PLOT'))
path += '/'
from plotMODULE import addtolegend, plotFromTEM, PlotOLDData


def main():
    datapath = 'pchannel.tem'
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

    plotFromTEM(path+datapath, mainplot, 1)
    plotFromTEM(path+datapath, subplot, 1, '', '--',oldpath+datapath )
    plotFromTEM(oldpath+datapath, mainplot, 1, '', '--')

    PlotOLDData('1tpr3e.xii94d', mainplot, 1, '', ':', 'k')
    PlotOLDData('2tpr1e.xii94d', mainplot, 1, '', ':', 'k')

    subplot.axis((0.0, 8000.0, -10, 10))
    mainplot.axis((0.0,8000.0,1e16,1e19))
    mainplot.set_ylabel(r'Concentration [cm$^{-3}$]')
    mainplot.set_xlabel(r'Depth [$\AA$]')
    subplot.set_xlabel(r'Depth [$\AA$]')
    #plt.savefig( sys.argv[0][:-3]+'.png', dpi=150)
    if warn:
        print(warn)
    plt.show()

main()
