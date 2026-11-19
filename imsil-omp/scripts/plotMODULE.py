#############################################################################
# python test script to take .his files and produce pngs
# setup options from pchannel.bak
# pchannel.bak simplified into english here:
# Title: IMSIL703 - bench 7: pchannel
# axes X=0..8000 titles Depth [A]
# axes Y=1e16..1e19 titles Concentration [cm^3]
# plot pchannel1.his a-Si
# plot pchannel2.his (100)-Si
# plot pchannel3.his (110)-Si
##############################################################################
import numpy as np
import matplotlib.pyplot as plt
import random as r
import sys, os

scripts_path = os.path.join( os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts')
sys.path.append(scripts_path)

filelist = []  #for color selection
commonTokens =[]


def colorcycler( anInteger ):
    '''takes an integer and returns a color'''
    colorlist = ['b', 'g', 'r', 'c', 'm', 'y', 'DarkSeaGreen', 'GreenYellow', 'MidnightBlue']
    return colorlist[anInteger%len(colorlist)]


def PlotOLDTXTData(filename, plottarget, columnToPlot=1, label='', linestyle='-', color='auto'):
    '''reads a .his file and plots the results'''
    shortfile = stripPath(filename)[1]
    print('plotting: '+filename)
    if shortfile not in filelist:
        filelist.append(shortfile)
    filenumber = filelist.index(shortfile) 
    if color == 'auto':
        color = colorcycler(int(filenumber))
    filereader = open(filename, 'r')
    columnToPlot = columnToPlot + 2    #old soitec data files use a different format
    Xvalue = []
    Yvalue = []
    filereader.readline()
    filereader.readline()
    filereader.readline()
    while (True):
        line = filereader.readline ()
        if not line:
            break
        if line.isspace():
            break
        line = line.replace('D', 'E')  #format change for old fortran code 
        data = line.split()
        Xvalue.append(float(data[2]))
        Yvalue.append(float(data[columnToPlot]))
    #print(Xvalue , Yvalue)
    filereader.close()
    Xvalue = np.array(Xvalue)
    Xvalue = Xvalue * 10.0
    l, = plottarget.plot(Xvalue,Yvalue, color=color, label= label, ls=linestyle)
    if l.get_linestyle() != '-':
        if l.get_linestyle() == '--':
            l.set_dashes([5,2,10,3])
        if l.get_linestyle() == ':':
            l.set_dashes([3,1,4,1])
        if l.get_linestyle() == '-.':
            l.set_dashes([11,2,4,2])
    if label != '':
            alegen = plottarget.legend(loc='best')
            alegen.get_frame().set_alpha(0.5)
            for t in alegen.get_texts():
                t.set_fontsize('x-small')

    return np.array(Xvalue), np.array(Yvalue)


def PlotOLDData(filename, plottarget, columnToPlot=1, label='', linestyle='-', color='auto'):
    '''reads a .his file and plots the results'''
    shortfile = stripPath(filename)[1]
    print('plotting: '+filename)
    if shortfile not in filelist:
        filelist.append(shortfile)
    filenumber = filelist.index(shortfile) 
    if color == 'auto':
        color = colorcycler(int(filenumber))
    filereader = open(filename, 'r')
    Xvalue = []
    Yvalue = []
    filereader.readline()
    filereader.readline()
    while (True):
        line = filereader.readline ()
        if not line:
            break
        line = line.replace('D', 'E')  #format change for old fortran code 
        data = line.split()
        Xvalue.append(float(data[0]))
        Yvalue.append(float(data[columnToPlot]))
    #print(Xvalue , Yvalue)
    filereader.close()
    Xvalue = np.array(Xvalue)
    Xvalue = Xvalue *10000.0
    l, = plottarget.plot(Xvalue,Yvalue, color=color, label= label, ls=linestyle)
    if l.get_linestyle() != '-':
        if l.get_linestyle() == '--':
            l.set_dashes([5,2,10,3])
        if l.get_linestyle() == ':':
            l.set_dashes([3,1,4,1])
        if l.get_linestyle() == '-.':
            l.set_dashes([11,2,4,2])
    if label != '':
            alegen = plottarget.legend(loc='best')
            alegen.get_frame().set_alpha(0.5)
            for t in alegen.get_texts():
                t.set_fontsize('x-small')

    return np.array(Xvalue), np.array(Yvalue)

def readHISTandPlot(filename, plottarget, columnToPlot=1, chunkToPlot=1, label='', linestyle='-', color='auto', scalefactorX=1, scalefactorY=1, polar=False):
    '''reads a .hist file and plots the results'''
    currentChunk=1
    shortfile = stripPath(filename)[1]
    print('plotting: '+filename)
    if shortfile not in filelist:
        filelist.append(shortfile)
    filenumber = filelist.index(shortfile)
    if color == 'auto':
        color = colorcycler(int(filenumber))
    filereader = open(filename, 'r')
    Xvalue = []
    Yvalue = []
    filereader.readline()
    filereader.readline()
    while (True):
        line = filereader.readline ()
        if not line:
            break
        if '.' not in line.split()[-1]:
            line = filereader.readline() 
            line = filereader.readline() 
            currentChunk = currentChunk + 1
        if currentChunk == chunkToPlot:
            line = line.replace('D', 'E')  #format change for old fortran code 
            data = line.split()
            Xvalue.append(float(data[0])*scalefactorX)
            Yvalue.append(float(data[columnToPlot])*scalefactorY)
    #print(Xvalue , Yvalue)
    filereader.close()
    Xvalue = np.array(Xvalue)
    Yvalue = np.array(Yvalue)
    if polar:
        plottarget.set_aspect('equal')
        Xvalue, Yvalue = Yvalue*360/np.pi * (np.sin(Xvalue*np.pi/180), np.cos(Xvalue*np.pi/180))
    print('label=', label)
    l, = plottarget.plot(Xvalue, Yvalue, color=color, label=label, ls=linestyle)
    if l.get_linestyle() != '-':
        if l.get_linestyle() == '--':
            l.set_dashes([5,2,10,3])
        if l.get_linestyle() == '-.':
            l.set_dashes([1,2,1,3])
    return (Xvalue, Yvalue, color)

def plotERROR( his_param1, his_param2, plottarget, label='', linestyle='-', color = 'auto'): 
    ''' reads two .his files and plots the error as determined by n1-n2/sqrt(n2)'''
    XvalueN1, YvalueN1, colorN1 = his_param1   
    XvalueN2, YvalueN2, colorN2 = his_param2   
    if color == 'auto':
        color = colorN1 
    XvalueN2 = np.array(XvalueN2)
    XvalueN1 = np.array(XvalueN1)
    YvalueN2 = np.array(YvalueN2)
    YvalueN1 = np.array(YvalueN1)
    XmaxN2 = np.nanmax(XvalueN2)
    XmaxN1 = np.nanmax(XvalueN1)
    XminN2 = np.nanmin(XvalueN2)
    XminN1 = np.nanmin(XvalueN1)
    if XmaxN1 > XmaxN2:
        bmax = XmaxN2
    else:
        bmax = XmaxN1
    if XminN1 > XminN2:
        bmin = XminN1
    else:
        bmin = XminN2
 
    xvals = np.linspace(bmin,bmax, 20000)
    yintN1 = np.interp(xvals, XvalueN1, YvalueN1)
    yintN2 = np.interp(xvals, XvalueN2, YvalueN2)
    #print(Xvalue , Yvalue)
    result = (yintN1- yintN2) / (np.sqrt(yintN1))
    C = np.sqrt(YvalueN1[-1]*10000)
    result = result/C
    l,= plottarget.plot( xvals, result,color=color, label = label)
    #plottarget.axhspan(-C,C, facecolor = (r.random(),r.random(), r.random()), alpha = 0.2)
    #plottarget.axhline(y= -C, ls=linestyle, color = 'k', label = 'normal range for '+label)
    #plottarget.axhline(y= C, ls=linestyle, color = 'k')
    plottarget.axhline(y= 0, ls=linestyle, color = 'k')
    colorUSED=l.get_color()
    plottarget.fill_between(xvals, 1, result, where=result>1, facecolor= colorUSED, alpha = 0.2)
    plottarget.fill_between(xvals, -1, result, where=result<-1, facecolor = colorUSED, alpha = 0.2)
    ylable = plottarget.set_ylabel(r' $(n_1 - n_2) / \sqrt{C \cdot n_1}$' )
    ylable.set_fontsize('xx-large')
    

def readHISandPlot(filename, plottarget, columnToPlot=1, label='', linestyle='-', color='auto', record=1):
    '''reads a .his file and plots the results'''
    shortfile = stripPath(filename)[1]
    if shortfile not in filelist:
        filelist.append(shortfile)
    filenumber = filelist.index(shortfile) 
    if color == 'auto':
        color = colorcycler(int(filenumber))
    filereader = open(filename, 'r')
    for _ in range(record):
        Xvalue = []
        Yvalue = []
        line = filereader.readline()
        try:
            ncol, npoints = line.split()
        except ValueError:
            print(filename+':', 'record', record, 'not available.')
            print('Skipping this file ...')
            return np.array([]), np.array([])
        filereader.readline()
        for ip in range(int(npoints)):
            line = filereader.readline ()
            if not line:
                break
            line = line.replace('D', 'E')  #format change for old fortran code 
            data = line.split()
            Xvalue.append(float(data[0]))
            Yvalue.append(float(data[columnToPlot]))
            #print(Xvalue[-1], Yvalue[-1])
    filereader.close()
    print('plotting:', filename, '('+color+')')
    l, = plottarget.plot(Xvalue, Yvalue, color=color, label= label, ls=linestyle)
    if l.get_linestyle() != '-':
        if l.get_linestyle() == '--':
            l.set_dashes([5,2,10,3])
        if l.get_linestyle() == '-.':
            l.set_dashes([1,2,1,3])
    if label != '':
        alegen = plottarget.legend(loc='best')
        alegen.get_frame().set_alpha(0.5)
        for t in alegen.get_texts():
            t.set_fontsize('small')
    return np.array(Xvalue), np.array(Yvalue)

def Subtract2HISabs(filename1, filename2, plottarget,column2ToPlot=2, column1ToPlot=1, 
                    label='', linestyle='-', color = 'auto'): 
    ''' reads two .his files and plots the error as determined by abs(n1-n2)'''
    shortfile = stripPath(filename1)[1]
    if shortfile not in filelist:
        filelist.append(shortfile)
    filenumber = filelist.index(shortfile) 
    if color == 'auto':
        color = colorcycler(int(filenumber))
    filereader = open(filename1, 'r')
    XvalueN1 = []
    YvalueN1 = []
    filereader.readline()
    filereader.readline()
    while (True):
        line = filereader.readline ()
        if not line:
            break
        line = line.replace('D', 'E')  #format change for old fortran code 
        data = line.split()
        XvalueN1.append(float(data[0]))
        YvalueN1.append(float(data[column1ToPlot]))
    #print(Xvalue , Yvalue)
    filereader.close()
    filereader = open(filename2, 'r')
    XvalueN2 = []
    YvalueN2 = []
    filereader.readline()
    filereader.readline()
    while (True):
        line = filereader.readline ()
        if not line:
            break
        line = line.replace('D', 'E')  #format change for old fortran code 
        data = line.split()
        XvalueN2.append(float(data[0]))
        YvalueN2.append(float(data[column2ToPlot]))
    XvalueN2 = np.array(XvalueN2)
    XvalueN1 = np.array(XvalueN1)
    YvalueN2 = np.array(YvalueN2)
    YvalueN1 = np.array(YvalueN1)
    XmaxN2 = np.nanmax(XvalueN2)
    XmaxN1 = np.nanmax(XvalueN1)
    XminN2 = np.nanmin(XvalueN2)
    XminN1 = np.nanmin(XvalueN1)
    if XmaxN1 > XmaxN2:
        bmax = XmaxN2
    else:
        bmax = XmaxN1
    if XminN1 > XminN2:
        bmin = XminN1
    else:
        bmin = XminN2
 
    xvals = np.linspace(bmin,bmax, 20000)
    yintN1 = np.interp(xvals, XvalueN1, YvalueN1)
    yintN2 = np.interp(xvals, XvalueN2, YvalueN2)
    #print(Xvalue , Yvalue)
    filereader.close()
    result = np.abs((yintN1- yintN2) )
    l,= plottarget.plot( xvals, result,color=color, label = label, ls=linestyle)
    if l.get_linestyle() != '-':
        if l.get_linestyle() == '--':
            l.set_dashes([5,2,10,3])
        if l.get_linestyle() == '-.':
            l.set_dashes([1,2,1,3])
    if label != '':
        alegen = plottarget.legend(loc='best')
        alegen.get_frame().set_alpha(0.5)
        for t in alegen.get_texts():
            t.set_fontsize('small')
    
   
def plotERRORfrom2HIS( filename1, filename2, plottarget, columnToPlot=1, label='', linestyle='-', color = 'auto'): 
    ''' reads two .his files and plots the error as determined by n1-n2/sqrt(n2)'''
    shortfile = stripPath(filename1)[1]
    if shortfile not in filelist:
        filelist.append(shortfile)
    filenumber = filelist.index(shortfile) 
    if color == 'auto':
        color = colorcycler(int(filenumber))
    filereader = open(filename1, 'r')
    XvalueN1 = []
    YvalueN1 = []
    filereader.readline()
    filereader.readline()
    while (True):
        line = filereader.readline ()
        if not line:
            break
        line = line.replace('D', 'E')  #format change for old fortran code 
        data = line.split()
        XvalueN1.append(float(data[0]))
        YvalueN1.append(float(data[columnToPlot]))
    #print(Xvalue , Yvalue)
    filereader.close()
    filereader = open(filename2, 'r')
    XvalueN2 = []
    YvalueN2 = []
    filereader.readline()
    filereader.readline()
    while (True):
        line = filereader.readline ()
        if not line:
            break
        line = line.replace('D', 'E')  #format change for old fortran code 
        data = line.split()
        XvalueN2.append(float(data[0]))
        YvalueN2.append(float(data[columnToPlot]))
    XvalueN2 = np.array(XvalueN2)
    XvalueN1 = np.array(XvalueN1)
    YvalueN2 = np.array(YvalueN2)
    YvalueN1 = np.array(YvalueN1)
    XmaxN2 = np.nanmax(XvalueN2)
    XmaxN1 = np.nanmax(XvalueN1)
    XminN2 = np.nanmin(XvalueN2)
    XminN1 = np.nanmin(XvalueN1)
    if XmaxN1 > XmaxN2:
        bmax = XmaxN2
    else:
        bmax = XmaxN1
    if XminN1 > XminN2:
        bmin = XminN1
    else:
        bmin = XminN2
 
    xvals = np.linspace(bmin,bmax, 20000)
    yintN1 = np.interp(xvals, XvalueN1, YvalueN1)
    yintN2 = np.interp(xvals, XvalueN2, YvalueN2)
    #print(Xvalue , Yvalue)
    filereader.close()
    result = (yintN1- yintN2) / (np.sqrt(yintN1))
    C = np.sqrt(YvalueN1[-1]*10000)
    result = result/C
    l,= plottarget.plot( xvals, result,color=color, label = label)
    #plottarget.axhspan(-C,C, facecolor = (r.random(),r.random(), r.random()), alpha = 0.2)
    #plottarget.axhline(y= -C, ls=linestyle, color = 'k', label = 'normal range for '+label)
    #plottarget.axhline(y= C, ls=linestyle, color = 'k')
    plottarget.axhline(y= 0, ls=linestyle, color = 'k')
    colorUSED=l.get_color()
    plottarget.fill_between(xvals, 1, result, where=result>1, facecolor= colorUSED, alpha = 0.2)
    plottarget.fill_between(xvals, -1, result, where=result<-1, facecolor = colorUSED, alpha = 0.2)
    ylable = plottarget.set_ylabel(r' $(n_1 - n_2) /\sqrt{C \cdot n_1}$' )
    ylable.set_fontsize('xx-large')
    
   
def addtolegend(label, plottarget, linestyle):
    '''adds an invisible point to the graph so as to add entries to the legend'''
    print('label=', label)
    l, = plottarget.plot((0.0,0.0),(0.0,0.1), color='k', label= label, ls=linestyle)
    if l.get_linestyle() != '-':
        l.set_dashes([5,2,10,3])
    alegen = plottarget.legend(loc='best')
    alegen.get_frame().set_alpha(0.5)
    for t in alegen.get_texts():
        t.set_fontsize('small')

def stripPath(filename):
    '''takes a path string and splits it by / and then returns a rebuilt string without the filename'''
    path = filename.split('/')
    thefile = path.pop(-1)
    stripped = ''
    while(len(path)>0):
       stripped =stripped + path.pop(0) + '/' 
        
    return stripped, thefile 

def plotFromTEM(filename,  plottarget, column=1, label='', linestyle='-', error='', ext='', record=1):
    ''' takes a TEM file and parses it and then plots the column of requested '''
    parfilename = filename[:-3] + 'par'
    path, shortname = stripPath(filename)
    title = shortname[:-3].upper() + 'tem: '
    filereader = open( parfilename, 'r')
    initShortenLabels(filename)
    parString = filereader.read()
    if parString.strip().startswith('FILE'):
        filereader.seek(0)
        inpFiles = []
        while True:
            line = filereader.readline ()
            if not line:
                break
            line = line.split()
            if line != []:
                if line[0] == 'FILE:':
                    inpFiles.append(line[1])
    else:
        from collections import OrderedDict
        param = OrderedDict()
        exec(parString)
        inpFiles = param.keys()
    for inpFile in inpFiles:
        plotFromINP(path+inpFile+'.inp', plottarget, column, title+label, 
                    linestyle, error, ext, record)

def plotYieldFromTEM(filename,  plottarget, column=1, label='', linestyle='-', color='auto', error='', ext=''):
    ''' takes a TEM file, parses it, extracts yields from the .out files, and then plots them.'''
    import read_output
    parfilename = filename[:-3] + 'par'
    path, shortname = stripPath(filename)
    if shortname not in filelist:
        filelist.append(shortname)
    filenumber = filelist.index(shortname) 
    if color == 'auto':
        color = colorcycler(int(filenumber))
    title = shortname[:-3].upper() + 'tem: '
    filereader = open( parfilename, 'r')
    initShortenLabels(filename)
    parString = filereader.read()
    if parString.strip().startswith('FILE'):
        filereader.seek(0)
        inpFiles = []
        while True:
            line = filereader.readline ()
            if not line:
                break
            line = line.split()
            if line != []:
                if line[0] == 'FILE:':
                    inpFiles.append(line[1])
    else:
        from collections import OrderedDict
        param = OrderedDict()
        exec(parString)
        inpFiles = param.keys()
    print('inpFiles:', inpFiles)
    Xvalues = []
    Yvalues = []
    Errors = []
    for inpFile in inpFiles:
        print('inpFile:', path+inpFile)
        outfilename = path+inpFile+'.out'
        print('outfilename',  outfilename)
        if linestyle == '-' and \
           os.path.getmtime(path+inpFile+'.inp') > os.path.getmtime(outfilename):
            print('')
            print(72*'*')
            print('*' + 8*' WARNING')
            print('* "' + path + inpFile + '.inp" is newer than "' + outfilename + '"')
            print(72*'*')
            input( '\nPress <CR> to continue: ' )
        if 'eclipse' in path:
            outfilereader = open(outfilename, 'r')
            while True:
                outline = outfilereader.readline()
                if not outline:
                    break
                outline = outline.replace('|', '').split()
                if 'Backscattered' in outline:
                    try:
                        Xvalues.append(float(inpFile.lstrip(shortname[:-4])))
                    except ValueError:
                        Xvalues.append(float(inpFile.lstrip(shortname[:-4]+'_')))
                    if Yvalues == []:
                        Yvalues = [[float(out)] for out in outline[1:]]
                    else:
                        for out, yvalues in zip(outline[1:], Yvalues):
                            yvalues.append(out)
                    break

                if linestyle == '-':
                    INPfilereader = open( outfilename[:-3]+'inp', 'r')
                    plottitle = INPfilereader.readline().strip() + '\n'
                    filelabel = INPfilereader.readline().strip()
                    filelabel, subtitle = shortenLabels(filelabel)
                    plottarget.set_title((plottitle+subtitle).strip('\n'))
        else:
            try:
                Xvalues.append(float(inpFile.lstrip(shortname[:-4])))
            except ValueError:
                try:
                    Xvalues.append(float(inpFile.lstrip(shortname[:-4]+'_')))
                except ValueError:
                    print('Could not extract x value from filename', inpFile)
                    print('Skipping this file ...')
                    continue
            Yvalue = read_output.read_out(outfilename, ('Yield moments', 'Backscattered', 'mean'), column)
            Yvalues.append(Yvalue)
#            print('Xvalues:', Xvalues)
#            print('Yvalues:', Yvalues)
    
            if linestyle == '-':
                Err = read_output.read_out(outfilename, ('Yield moments', 'Backscattered', 'mean', 'mean'), column)
                Errors.append(Err)
#                print('Errors:', Errors)

                INPfilereader = open( outfilename[:-3]+'inp', 'r')
                plottitle = INPfilereader.readline().strip() + '\n'
                filelabel = INPfilereader.readline().strip()
                filelabel, subtitle = shortenLabels(filelabel)
                plottarget.set_title((plottitle+subtitle).strip('\n'))

#    print(Xvalues, Yvalues[column-1])
    if linestyle != '-':
        label = ''
    if 'eclipse' in path:
        l, = plottarget.plot(Xvalues, Yvalues[column-1], color=color, label=label, ls=linestyle)
    else:
        if linestyle == '-':
            l, _, __ = plottarget.errorbar(Xvalues, Yvalues, Errors, color=color, label=label, ls=linestyle)
        else:
            l, = plottarget.plot(Xvalues, Yvalues, color=color, label=label, ls=linestyle)
    if l.get_linestyle() != '-':
        if l.get_linestyle() == '--':
            l.set_dashes([5,2,10,3])
        if l.get_linestyle() == '-.':
            l.set_dashes([1,2,1,3])
    if label != '':
        alegen = plottarget.legend(loc='best')
        alegen.get_frame().set_alpha(0.5)
        for t in alegen.get_texts():
            t.set_fontsize('small')

#    return np.array(Xvalues), np.array(Yvalues[column-1])
    return np.array(Xvalues), np.array(Yvalues)
        

def plotFromINP(filename, plottarget, column=1, title='', linestyle='-', error='', ext='', record=1):
    ''' takes an INP file and plots the .HIS and parses titles from INP'''
    if ext == '':
        HISfilename = filename[:-3] + 'his'
    else:
        HISfilename = filename[:-3] + ext
    if not os.path.exists(HISfilename):
        print('"' + HISfilename + '" does not exist')
        print('Skipping this file ...')
        return
    if os.path.getmtime(filename) > os.path.getmtime(HISfilename):
        print('')
        print(72*'*')
        print('*' + 8*' WARNING')
        print('* "' + filename + '" is newer than "' + HISfilename + '"')
        print(72*'*')
        input( '\nPress <CR> to continue: ' )
    INPfilereader = open( filename, 'r')

    if error == '':  #normal mode don't create error hist
        plottitle = INPfilereader.readline().strip() + '\n'
        filelabel = INPfilereader.readline().strip()
        filelabel, subtitle = shortenLabels(filelabel)
        if linestyle != '-':
            filelabel = ''
        readHISandPlot(HISfilename, plottarget, column, filelabel, linestyle, record=record)
        plottarget.set_title((title+plottitle+subtitle).strip('\n'))
  
    else:    # an alternate file name has been specified in error so we want to plot relative error
        plottitle = INPfilereader.readline()
        filelabel = shortenLabels(INPfilereader.readline().strip())[0]
        #shortlabel = filelabel.split(",")[0]
        shortlabel = title + filelabel
        otherpath = stripPath(error)[0]+stripPath(HISfilename)[1]
        plotERRORfrom2HIS(HISfilename, otherpath, plottarget, column, shortlabel, linestyle)

def initShortenLabels(filename):
    '''initialize the commonTokens in the shortenLabels function with a filename different 
       from the first file to be plotted.''' 
    temreader  = open(filename)
    temreader.readline().strip()
    shortenLabels(temreader.readline().strip())   #prime the shortener
            
def shortenLabels(labelstring):
    '''builds a list and checks a string vs the list to determine terms that never change'''
    global commonTokens
    setunique = set(labelstring.split(','))
    if commonTokens == []:    #first use
        commonTokens = set(labelstring.split(','))
    else:
        commonTokens.intersection_update(setunique) 
        setunique.difference_update(commonTokens)
    returnstring = ''
    for token in labelstring.split(','):
        if token in setunique:
            returnstring = returnstring + token + ', '
    returnstring = returnstring.strip()
    returnstring = returnstring.rstrip(',')
    commonstring = ''
    for token in labelstring.split(','):
        if token in commonTokens:
            commonstring = commonstring + token + ', '
    commonstring = commonstring.strip()
    commonstring = commonstring.rstrip(',')
    return returnstring, commonstring

def titleFromINP( filename, plottarget):
    '''reads an input file and parses the first two lines to produce a title'''
    INPfilereader = open( filename, 'r')
    plottitle = INPfilereader.readline().strip() + '\n'
    filelabel = INPfilereader.readline().strip()
    plottarget.set_title((plottitle+filelabel).strip('\n'))

