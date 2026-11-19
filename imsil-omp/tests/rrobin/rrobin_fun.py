import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(__file__)
print("path:", path)
sys.path.append(os.path.join(path, '..', 'PLOT'))
path += '/'

from plotMODULE import readHISTandPlot, plotERROR, addtolegend, titleFromINP


def rrobin_fun(datapath):

    root = path + datapath
    extAAT = '.hisaat'
    extAAB = '.hisaab'
    extET = '.hiset'
    extEB = '.hiseb'
    chunk = 1
    if not os.path.exists(root + extAAT):
        extAAT = '.hist'
        extAAB = '.hisb'
        extET = '.hist'
        extEB = '.hisb'
        if 'eclipse' in path:
            chunk = 2
        else:
            chunk = 3
    
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
        
    oldroot = oldpath + datapath
    oldextAAT = '.hisaat'
    oldextAAB = '.hisaab'
    oldextET = '.hiset'
    oldextEB = '.hiseb'
    oldchunk = 1
    if not os.path.exists(oldroot + oldextAAT):
        oldextAAT = '.hist'
        oldextAAB = '.hisb'
        oldextET = '.hist'
        oldextEB = '.hisb'
        if 'eclipse' in oldpath:
            oldchunk = 2
        else:
            oldchunk = 3

    fig = plt.figure(figsize= (16.5,11.7))
    mainplot = list()
    mainplot.append(fig.add_subplot(2,2,1))
    mainplot.append(fig.add_subplot(2,2,3))
    subplot = list()
    subplot.append(fig.add_subplot(2,2,2))
    subplot.append(fig.add_subplot(2,2,4))
   
    mainplot[0].set_yscale('log')
    mainplot[1].set_yscale('log')
    
    rrobinHIST = readHISTandPlot(root+extAAT, mainplot[0], 1, chunk, 'transmitted', '-')
    rrobinHISB = readHISTandPlot(root+extAAB, mainplot[0], 1, chunk, 'backscattered')
    erobinHIST = readHISTandPlot(root+extET, mainplot[1], 1, 1, '', '-', 'auto', 1, 1)
    erobinHISB = readHISTandPlot(root+extEB, mainplot[1], 1, 1, '', '-', 'auto', 1, 1)

    if oldpath:
        rrobinHISTold = readHISTandPlot(oldroot+oldextAAT, mainplot[0], 1, oldchunk, '', '--')
        rrobinHISBold = readHISTandPlot(oldroot+oldextAAB, mainplot[0], 1, oldchunk,'', '--')
        # this assumes that the old energy distributions are still in keV
        erobinHISTold = readHISTandPlot(oldroot+oldextET, mainplot[1], 1, 1, '', '--', 'auto') 
        erobinHISBold = readHISTandPlot(oldroot+oldextEB, mainplot[1], 1, 1, '', '--', 'auto') 
        plotERROR(rrobinHIST, rrobinHISTold, subplot[0])
        plotERROR(rrobinHISB, rrobinHISBold, subplot[0])
        plotERROR(erobinHIST, erobinHISTold, subplot[1])
        plotERROR(erobinHISB, erobinHISBold, subplot[1])

    addtolegend('Current Data', mainplot[0], '-')
    if oldpath:
        addtolegend('Reference', mainplot[0], '--')
    titleFromINP(root + '.inp', mainplot[0])
    subplot[0].axis((0.0, 90.0, -10, 10))
    subplot[1].axis((0.0, 200.0, -10, 10))
    mainplot[0].axis((0.0, 90.0, 1e-4, 1e-1))
    mainplot[1].axis((0.0, 200.0, 1e-4, 1e-1))
    ylabel_m = r'Distribution [$^{\circ^{-1}}$]'
    ylabel_s = r'Distribution [eV$^{-1}$]'
    xlabel_m = r'Angle [$^{\circ}$]'
    mainplot[0].set_ylabel(ylabel_m)
    mainplot[0].set_xlabel(xlabel_m)
    mainplot[1].set_ylabel(ylabel_s)
    mainplot[1].set_xlabel('Energy [eV]')
    #plt.savefig(datapath+'.png', dpi=150)
    if warn:
        print(warn)
    plt.show()

