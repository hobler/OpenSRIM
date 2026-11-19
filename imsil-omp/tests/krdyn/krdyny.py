import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(os.path.abspath(__file__))
print("path:", path)
basename = os.path.basename(path)
sys.path.append(os.path.join(path, '..', 'PLOT'))
path += '/'
scripts_path = os.path.join( os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts')
sys.path.append(scripts_path)
import read_output

from plotMODULE import plotYieldFromTEM, addtolegend


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
    
    fig = plt.figure(figsize= (16.5,11.7))
    mainplot = fig.add_subplot(1,2,1)
    subplot = fig.add_subplot(1,2,2)

    plotYieldFromTEM(path+datapath, mainplot, 3, 'Si sputtered', '-', 'r')
    plotYieldFromTEM(path+datapath, mainplot, 2, 'Kr sputtered', '-', 'g')
    plotYieldFromTEM(path+datapath, mainplot, 1, 'Kr backscattered', '-', 'b')
    plotYieldFromTEM(path+datapath, subplot, 1, 'Kr backscattered', '-', 'b')
    if oldpath:
        plotYieldFromTEM(oldpath+datapath, mainplot, 3, '', '--', 'r')
        plotYieldFromTEM(oldpath+datapath, mainplot, 2, '', '--', 'g')
        plotYieldFromTEM(oldpath+datapath, mainplot, 1, '', '--', 'b')
        plotYieldFromTEM(oldpath+datapath, subplot, 1, '', '--', 'c')

    yield1 = read_output.read_out('krdyn_10e15_stat.out', ('Yield moments', 'Backscattered', 'mean'), 1)
    yield2 = read_output.read_out('krdyn_10e15_stat.out', ('Yield moments', 'Backscattered', 'mean'), 2)
    yield3 = read_output.read_out('krdyn_10e15_stat.out', ('Yield moments', 'Backscattered', 'mean'), 3)
    mainplot.plot((0,1e16), (yield1, yield1), 'k:', label='static')
    mainplot.plot((0,1e16), (yield2, yield2), 'k:')
    mainplot.plot((0,1e16), (yield3, yield3), 'k:')
    
#    mainplot.axis((-90.0, 90.0, 0, 12))
#    subplot.axis((-90.0, 90.0, 0., 1.))
    ylabel_m = r'Yield'
    xlabel_m = r'Dose [cm$^{-2}$]'
    mainplot.set_ylabel(ylabel_m)
    mainplot.set_xlabel(xlabel_m)
    subplot.set_ylabel('Backscattering ' + ylabel_m)
    subplot.set_xlabel(xlabel_m)
#    mainplot.set_xscale('log')
#    subplot.set_xscale('log')

    subplot.set_title('')

    addtolegend('Current Data', mainplot, '-')
    if oldpath:
        addtolegend('Reference', mainplot, '--')
    addtolegend('Current Data', subplot, '-')
    if oldpath:
        addtolegend('Reference', subplot, '--')

    mainplot.legend(loc='best')
    subplot.legend(loc='best')

    #plt.savefig('arsputtery.png', dpi=150)
    if warn:
        print(warn)
    plt.show()

main()
