import numpy as np
import matplotlib.pyplot as plt
import random as r
import os, sys
path = os.path.dirname(os.path.abspath(__file__))
print("path:", path)
basename = os.path.basename(path)
sys.path.append(os.path.join(path, '..', 'PLOT'))
path += '/'
from plotMODULE import addtolegend, readHISandPlot, Subtract2HISabs,  \
                       plotERRORfrom2HIS, titleFromINP


def main():
    datapath = 'asneti%d.his'
    pathletters = [1, 2]
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
    addtolegend('Current Data', mainplot[0], '-')
    addtolegend('Reference', mainplot[0], '--')
    for pathletter in pathletters:
        readHISandPlot(path+datapath%pathletter, mainplot[pathletter-1], 1, 'Ions', '-', 'b')
        readHISandPlot(path+datapath%pathletter, mainplot[pathletter-1], 2, 'Interstitials', '-', 'r')
        Subtract2HISabs(path+datapath%pathletter, path+datapath%pathletter, mainplot[pathletter-1],
                        2, 4, 'abs(Net I-V)', '-', 'g')
        if oldpath:
            readHISandPlot(oldpath+datapath%pathletter, mainplot[pathletter-1], 1, '', '--', 'b')
            readHISandPlot(oldpath+datapath%pathletter, mainplot[pathletter-1], 2, '', '--', 'r')
            Subtract2HISabs(oldpath+datapath%pathletter, oldpath+datapath%pathletter, 
                            mainplot[pathletter-1], 2, 4, '', '--', 'g')
            plotERRORfrom2HIS(path+datapath%pathletter,oldpath+datapath%pathletter, 
                              subplot[pathletter-1], 1, 'Ions', '--', 'b')
            plotERRORfrom2HIS(path+datapath%pathletter, oldpath+datapath%pathletter, 
                              subplot[pathletter-1], 2, 'Interstitials', '--', 'r')
            plotERRORfrom2HIS(path+datapath%pathletter, oldpath+datapath%pathletter, 
                              subplot[pathletter-1], 4, 'Vacancies', '--', 'r')
        titleFromINP('asneti%d.inp'%pathletter, mainplot[pathletter-1])
    
    mainplot[0].axis((0.0,1000.0,1e2,1e8))
    mainplot[1].axis((0.0,1000.0,1e17,1e23))
    mainplot[0].set_ylabel('Concentration [cm$^{-3}$]')
    mainplot[1].set_ylabel('Concentration [cm$^{-3}$]')
    mainplot[0].set_xlabel(r'Depth [$\AA$]')
    mainplot[1].set_xlabel(r'Depth [$\AA$]')
    subplot[0].axis((0.0, 1000.0, -10, 10))
    subplot[1].axis((0.0, 1000.0, -10, 10))
    subplot[0].set_xlabel(r'Depth [$\AA$]')
    #plt.savefig( sys.argv[0][:-3]+'.png', dpi=150)
    if warn:
        print(warn)
    plt.show()

main()
