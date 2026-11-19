import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(os.path.abspath(__file__))
print("path:", path)
basename = os.path.basename(path)
sys.path.append(os.path.join(path, '..', 'PLOT'))
path += '/'

from plotMODULE import readHISTandPlot, plotERROR, titleFromINP


def main():
    pathletters = ['0', '10', '20', '30', '40', '50', '60', '70', '75',
                   '80', '85', '87', '88', '89']
    
    root = path + 'arsputter%s'
    ext = '.hisab'
    chunk = 1
    if not os.path.exists(root%0 + ext):
        ext = '.hisb'
        chunk = 2
    
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
    oldroot = oldpath + 'arsputter%s'
    oldext = '.hisab'
    oldchunk = 1
    if not os.path.exists(oldroot%0 + oldext):
        oldext = '.hisb'
        oldchunk = 2
    
    fig = plt.figure(figsize= (16.5,11.7))
    mainplot = fig.add_subplot(1,2,1)
    subplot = fig.add_subplot(1,2,2)

    for pathletter in pathletters:
        HISB = readHISTandPlot(root%pathletter+ext, mainplot, 1, chunk, 
                               label=pathletter+' tilt', linestyle='-')
#        addtolegend('tilt='+pathletter, mainplot, '-')
        if oldpath:
            HISBold = readHISTandPlot(oldroot%pathletter+oldext, mainplot, 1, oldchunk, 
                                      '', '--')
            plotERROR(HISB, HISBold, subplot)

    titleFromINP(root%0+'.inp', mainplot)
    mainplot.set_title(mainplot.get_title().rstrip(', '+pathletters[-1]+' tilt')+
                       '\n Backscattering angular distribution')
    mainplot.set_yscale('log')
    mainplot.axis((-90.0, 90.0, 1e-5, 1.))
#    mainplot.axis((-90.0, 90.0, 0, 0.12))
    subplot.axis((-90.0, 90.0, -10., 10.))
    ylabel_m = r'Distribution [$^{\circ^{-1}}$]'
    xlabel_m = r'Angle [$^{\circ}$]'
    mainplot.set_ylabel(ylabel_m)
    mainplot.set_xlabel(xlabel_m)

    mainplot.legend(loc='best')
#    addtolegend('Current Data', mainplot, '-')
#    addtolegend('2009 Data', mainplot, '--')

    #plt.savefig('arsputterb.png', dpi=150)
    if warn:
        print(warn)
    plt.show()

main()
