import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(__file__)
sys.path.append(os.path.join(path, '..', 'PLOT'))
#path += '/'

from plotMODULE import readHISTandPlot, plotERROR, addtolegend, titleFromINP


def main():
    root_Legendre = os.path.join(path, 'rrobin_quad_legendre')
    fnameAAT_Legendre = root_Legendre + '.hisaat'
    fnameAAB_Legendre = root_Legendre + '.hisaab'
    fnameET_Legendre = root_Legendre + '.hiset'
    fnameEB_Legendre = root_Legendre + '.hiseb'
    chunkLegendre = 1
    if not os.path.exists(fnameAAT_Legendre):
        fnameAAT_Legendre = root_Legendre + '.hist'      # old file fnameensions
        fnameAAB_Legendre = root_Legendre + '.hisb'
        fnameET_Legendre = root_Legendre + '.hist'
        fnameEB_Legendre = root_Legendre + '.hisb'
        chunkLegendre = 3

    root_Mehler = os.path.join(path, 'rrobin_quad_mehler')
    fnameAAT_Mehler = root_Mehler + '.hisaat'
    fnameAAB_Mehler = root_Mehler + '.hisaab'
    fnameET_Mehler = root_Mehler + '.hiset'
    fnameEB_Mehler = root_Mehler + '.hiseb'

    root  = os.path.join(path, 'rrobin' )
    fnameAAT = root + '.hisaat'
    fnameAAB = root + '.hisaab'
    fnameET = root + '.hiset'
    fnameEB = root + '.hiseb'
    chunk = 1
    if not os.path.exists(fnameAAT):
        fnameAAT = root + '.hist'      # old file extensions
        fnameAAB = root + '.hisb'
        fnameET = root + '.hist'
        fnameEB = root + '.hisb'
        chunk = 3

    fig = plt.figure(figsize= (16.5,11.7))
    mainplot = list()
    mainplot.append(fig.add_subplot(2,2,1))
    mainplot.append(fig.add_subplot(2,2,3))
    subplot = list()
    subplot.append(fig.add_subplot(2,2,2))
    subplot.append(fig.add_subplot(2,2,4))
   
    mainplot[0].set_yscale('log')
    mainplot[1].set_yscale('log')
    rrobinHISTLegendre = readHISTandPlot(fnameAAT_Legendre, mainplot[0], 1, chunkLegendre, 
                                         'transmitted', '-', 'C0')
    rrobinHISTMehler = readHISTandPlot(fnameAAT_Mehler, mainplot[0], 1, 1, '', ':', 'C1')
    rrobinHIST = readHISTandPlot(fnameAAT, mainplot[0], 1, chunk,'', '--', 'C3')
    rrobinHISBLegendre = readHISTandPlot(fnameAAB_Legendre, mainplot[0], 1, chunkLegendre, 
                                         'backscattered', '-', 'C4')
    rrobinHISBMehler = readHISTandPlot(fnameAAB_Mehler, mainplot[0], 1, 1, '', ':', 'C5')
    rrobinHISB = readHISTandPlot(fnameAAB, mainplot[0], 1, chunk,'', '--', 'C7')
    erobinHISTLegendre = readHISTandPlot(fnameET_Legendre, mainplot[1], 1, 1, '', '-', 'C0')
    erobinHISTMehler = readHISTandPlot(fnameET_Mehler, mainplot[1], 1, 1, '', ':', 'C1')
    erobinHIST = readHISTandPlot(fnameET, mainplot[1], 1, 1, '', '--', 'C3') 
    erobinHISBLegendre = readHISTandPlot(fnameEB_Legendre, mainplot[1], 1, 1, '', '-', 'C4')
    erobinHISBMehler = readHISTandPlot(fnameEB_Mehler, mainplot[1], 1, 1, '', ':', 'C5')
    erobinHISB = readHISTandPlot(fnameEB, mainplot[1], 1, 1, '', '--', 'C7') 
    plotERROR(rrobinHISTLegendre, rrobinHIST, subplot[0])
    plotERROR(rrobinHISBLegendre, rrobinHISB, subplot[0])
    plotERROR(erobinHISTLegendre, erobinHIST, subplot[1])
    plotERROR(erobinHISBLegendre, erobinHISB, subplot[1])
    plotERROR(rrobinHISTMehler, rrobinHIST, subplot[0])
    plotERROR(rrobinHISBMehler, rrobinHISB, subplot[0])
    plotERROR(erobinHISTMehler, erobinHIST, subplot[1])
    plotERROR(erobinHISBMehler, erobinHISB, subplot[1])

    addtolegend('rrobin_quad_legendre', mainplot[0], '-')
    addtolegend('rrobin_quad_mehler', mainplot[0], ':')
    addtolegend('rrobin', mainplot[0], '--')
    titleFromINP(root_Legendre+'.inp', mainplot[0])
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
