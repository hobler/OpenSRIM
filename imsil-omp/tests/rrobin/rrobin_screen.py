import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(__file__)
sys.path.append(os.path.join(path, '..', 'PLOT'))
#path += '/'

from plotMODULE import readHISTandPlot, plotERROR, addtolegend, titleFromINP


def main():
    root = os.path.join(path, 'rrobin_screen_ZBLuniv')
    extAAT = '.hisaat'
    extAAB = '.hisaab'
    extET = '.hiset'
    extEB = '.hiseb'
    chunk = 1
    if not os.path.exists(root + extAAT):
        extAAT = '.hist'      # old file extensions
        extAAB = '.hisb'
        extET = '.hist'
        extEB = '.hisb'
        chunk = 3

    root_ZBL = os.path.join(path, 'rrobin_screen_ZBLspec')
    extAAT_ZBL = '.hisaat'
    extAAB_ZBL = '.hisaab'
    extET_ZBL = '.hiset'
    extEB_ZBL = '.hiseb'

    root_DMol = os.path.join(path, 'rrobin_screen_DMol')
    extAAT_DMol = '.hisaat'
    extAAB_DMol = '.hisaab'
    extET_DMol = '.hiset'
    extEB_DMol = '.hiseb'
        
    oldroot  = os.path.join(path, 'rrobin' )
    oldextAAT = '.hisaat'
    oldextAAB = '.hisaab'
    oldextET = '.hiset'
    oldextEB = '.hiseb'
    oldchunk = 1
    if not os.path.exists(oldroot + oldextAAT):
        oldextAAT = '.hist'      # old file extensions
        oldextAAB = '.hisb'
        oldextET = '.hist'
        oldextEB = '.hisb'
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
    rrobinHIST = readHISTandPlot(root+extAAT, mainplot[0], 1, chunk, 'transmitted', '-', 'C0')
    rrobinHIST_ZBL = readHISTandPlot(root_ZBL+extAAT_ZBL, mainplot[0], 1, 1, '', ':', 'C1')
    rrobinHIST_DMol = readHISTandPlot(root_DMol+extAAT_DMol, mainplot[0], 1, 1, '', '-.', 'C2')
    rrobinHISTold =readHISTandPlot(oldroot+oldextAAT, mainplot[0], 1, oldchunk,'', '--', 'C3')
    rrobinHISB = readHISTandPlot(root+extAAB, mainplot[0], 1, chunk, 'backscattered', '-', 'C4')
    rrobinHISB_ZBL = readHISTandPlot(root_ZBL+extAAB_ZBL, mainplot[0], 1, 1, '', ':', 'C5')
    rrobinHISB_DMol = readHISTandPlot(root_DMol+extAAB_DMol, mainplot[0], 1, 1, '', '-.', 'C6')
    rrobinHISBold = readHISTandPlot(oldroot+oldextAAB, mainplot[0], 1, oldchunk,'', '--', 'C7')
    erobinHIST = readHISTandPlot(root+extET, mainplot[1], 1, 1, '', '-', 'C0')
    erobinHIST_ZBL = readHISTandPlot(root_ZBL+extET_ZBL, mainplot[1], 1, 1, '', ':', 'C1')
    erobinHIST_DMol = readHISTandPlot(root_DMol+extET_DMol, mainplot[1], 1, 1, '', '-.', 'C2')
    erobinHISTold = readHISTandPlot(oldroot+oldextET, mainplot[1], 1, 1, '', '--', 'C3') 
    erobinHISB = readHISTandPlot(root+extEB, mainplot[1], 1, 1, '', '-', 'C4')
    erobinHISB_ZBL = readHISTandPlot(root_ZBL+extEB_ZBL, mainplot[1], 1, 1, '', ':', 'C5')
    erobinHISB_DMol = readHISTandPlot(root_DMol+extEB_DMol, mainplot[1], 1, 1, '', ':', 'C6')
    erobinHISBold = readHISTandPlot(oldroot+oldextEB, mainplot[1], 1, 1, '', '--', 'C7') 
    plotERROR(rrobinHIST, rrobinHISTold, subplot[0])
    plotERROR(rrobinHISB, rrobinHISBold, subplot[0])
    plotERROR(erobinHIST, erobinHISTold, subplot[1])
    plotERROR(erobinHISB, erobinHISBold, subplot[1])

    addtolegend('rrobin_screen_ZBLuniv', mainplot[0], '-')
    addtolegend('rrobin_screen_ZBLspec', mainplot[0], ':')
    addtolegend('rrobin_screen_DMol', mainplot[0], '-.')
    addtolegend('rrobin', mainplot[0], '--')
    titleFromINP(root+'.inp', mainplot[0])
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
    #plt.savefig(path+'.png', dpi=150)
    plt.show()

main()
