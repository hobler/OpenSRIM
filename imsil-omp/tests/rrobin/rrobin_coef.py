import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(__file__)
sys.path.append(os.path.join(path, '..', 'PLOT'))
#path += '/'

from plotMODULE import readHISTandPlot, plotERROR, addtolegend, titleFromINP


def main():
    root_ZBLuniv = os.path.join(path, 'rrobin_coef_ZBLuniv')
    fnameAAT_ZBLuniv = root_ZBLuniv + '.hisaat'
    fnameAAB_ZBLuniv = root_ZBLuniv + '.hisaab'
    fnameET_ZBLuniv = root_ZBLuniv + '.hiset'
    fnameEB_ZBLuniv = root_ZBLuniv + '.hiseb'

    root_ZBLspec = os.path.join(path, 'rrobin_coef_ZBLspec')
    fnameAAT_ZBLspec = root_ZBLspec + '.hisaat'
    fnameAAB_ZBLspec = root_ZBLspec + '.hisaab'
    fnameET_ZBLspec = root_ZBLspec + '.hiset'
    fnameEB_ZBLspec = root_ZBLspec + '.hiseb'

    root_NLH = os.path.join(path, 'rrobin_coef_NLH')
    fnameAAT_NLH = root_NLH + '.hisaat'
    fnameAAB_NLH = root_NLH + '.hisaab'
    fnameET_NLH = root_NLH + '.hiset'
    fnameEB_NLH = root_NLH + '.hiseb'
        
    root_NLHlin = os.path.join(path, 'rrobin_coef_NLHlin' )
    fnameAAT_NLHlin = root_NLHlin + '.hisaat'
    fnameAAB_NLHlin = root_NLHlin + '.hisaab'
    fnameET_NLHlin = root_NLHlin + '.hiset'
    fnameEB_NLHlin = root_NLHlin + '.hiseb'

    fig = plt.figure(figsize= (16.5,11.7))
    mainplot = list()
    mainplot.append(fig.add_subplot(2,2,1))
    mainplot.append(fig.add_subplot(2,2,3))
    subplot = list()
    subplot.append(fig.add_subplot(2,2,2))
    subplot.append(fig.add_subplot(2,2,4))
   
    mainplot[0].set_yscale('log')
    mainplot[1].set_yscale('log')
    readHISTandPlot(fnameAAT_ZBLuniv, mainplot[0], 1, 1, 'transmitted', '-', 'C0')
    readHISTandPlot(fnameAAT_ZBLspec, mainplot[0], 1, 1, '', ':', 'C1')
    readHISTandPlot(fnameAAT_NLH, mainplot[0], 1, 1, '', '-.', 'C2')
    readHISTandPlot(fnameAAT_NLHlin, mainplot[0], 1, 1,'', '--', 'C3')
    readHISTandPlot(fnameAAB_ZBLuniv, mainplot[0], 1, 1, 'backscattered', '-', 'C4')
    readHISTandPlot(fnameAAB_ZBLspec, mainplot[0], 1, 1, '', ':', 'C5')
    readHISTandPlot(fnameAAB_NLH, mainplot[0], 1, 1, '', '-.', 'C6')
    readHISTandPlot(fnameAAB_NLHlin, mainplot[0], 1, 1,'', '--', 'C7')
    readHISTandPlot(fnameET_ZBLuniv, mainplot[1], 1, 1, '', '-', 'C0')
    readHISTandPlot(fnameET_ZBLspec, mainplot[1], 1, 1, '', ':', 'C1')
    readHISTandPlot(fnameET_NLH, mainplot[1], 1, 1, '', '-.', 'C2')
    readHISTandPlot(fnameET_NLHlin, mainplot[1], 1, 1, '', '--', 'C3') 
    readHISTandPlot(fnameEB_ZBLuniv, mainplot[1], 1, 1, '', '-', 'C4')
    readHISTandPlot(fnameEB_ZBLspec, mainplot[1], 1, 1, '', ':', 'C5')
    readHISTandPlot(fnameEB_NLH, mainplot[1], 1, 1, '', ':', 'C6')
    readHISTandPlot(fnameEB_NLHlin, mainplot[1], 1, 1, '', '--', 'C7') 

    addtolegend('rrobin_coef_ZBLuniv', mainplot[0], '-')
    addtolegend('rrobin_coef_ZBLspec', mainplot[0], ':')
    addtolegend('rrobin_coef_NLH', mainplot[0], '-.')
    addtolegend('rrobin_coef_NLHlin', mainplot[0], '--')
    titleFromINP(root_ZBLuniv + '.inp', mainplot[0])
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
