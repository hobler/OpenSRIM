import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(os.path.abspath(__file__))
print("path:", path)
basename = os.path.basename(path)
sys.path.append(os.path.join(path, '..', 'PLOT'))
path += '/'

from plotMODULE import plotYieldFromTEM, addtolegend


def main():
    datapath = 'arsputter.tem'
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
    mainplot = fig.add_subplot(1,2,1)
    subplot = fig.add_subplot(1,2,2)

    plotYieldFromTEM(path+datapath, mainplot, 2, 'sputter yield', '-')
    plotYieldFromTEM(path+datapath, subplot, 1, 'backscatter yield', '-')

    if oldpath:
        plotYieldFromTEM(oldpath+datapath, mainplot, 2, '', '--')
        plotYieldFromTEM(oldpath+datapath, subplot, 1, '', '--')

    mainplot.axis((-90.0, 90.0, 0, 12))
    subplot.axis((-90.0, 90.0, 0., 1.))
    ylabel_m = r'Yield'
    xlabel_m = r'Incidence Angle [$^{\circ}$]'
    mainplot.set_ylabel(ylabel_m)
    mainplot.set_xlabel(xlabel_m)
    subplot.set_ylabel(ylabel_m)
    subplot.set_xlabel(xlabel_m)

    subplot.set_title('')

    addtolegend('Current Data', mainplot, '-')
    addtolegend('Reference', mainplot, '--')
    addtolegend('Current Data', subplot, '-')
    addtolegend('Reference', subplot, '--')

    mainplot.legend(loc='best')
    subplot.legend(loc='best')

    #plt.savefig('arsputtery.png', dpi=150)
    if warn:
        print(warn)
    plt.show()

main()
