#! /usr/bin/env python3

"""
Plot data given in one of several file formats.

Usage: python3 plot.py [type] filename

where type is the file type ('cell2' or 'surf') and filename the name of the 
file containing the data
"""

import sys, os, shutil
from copy import deepcopy
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.transforms as mtransforms
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize


class FileReader(object):
    """
    Manage reading from a file (common code for specialized readers).
    """
    
    def __init__(self, filename):
        self.f = open(filename)

    def line(self):
        """
        Read and return next line.
        """
        return self.f.readline()
        
    def close(self):
        """
        Close the file.
        """
        self.f.close()

    def dose(self):
        """
        Read the dose from cell file.
        """
        # determine grid type
        line = self.line()
        if not line:
            raise EOFError
        c, d, dose = line.split()
        dose = float(dose)
        return dose


class Cell2FileReader(FileReader):
    """
    Manage reading from the surface file.
    """
    
    def grid(self):
        """
        Read grid from cell file.
        """
        # determine grid type
        line = self.line()
        if not line:
            raise EOFError
        c, g, grid_type = line.split()
        if g != 'Grid:':
            sys.exit('Cannot read grid from CELL2 file.')

        # depending on grid type, read or determine coordinates
        if grid_type == 'regular':
            line = self.line()
            nx, ny = [int(var) for var in line.split()]
            line = self.line()
            x0, y0 = [float(var) for var in line.split()]
            line = self.line()
            wx, wy = [float(var) for var in line.split()]
            x = [x0+i*wx for i in range(nx+1)]
            y = [y0+i*wy for i in range(ny+1)]
            x, y = np.meshgrid(x, y)
        elif grid_type == 'curvilinear':
            line = self.line()
            nx, ny = [int(var) for var in line.split()]
            x = []
            y = []
            for iy in range(ny+1):
                xx = []
                yy = []
                for ix in range(nx+1):
                    line = self.line()
                    xxx, yyy = [float(var) for var in line.split()]
                    xx.append(xxx)
                    yy.append(yyy)
                x.append(xx)
                y.append(yy)
            x = np.array(x)
            y = np.array(y)
        else:
            exit('Unknown grid type: ' + grid_type)

        return x, y
    
    def densities(self):
        """
        Read densities from cell file.
        """
        # skip first line
        line = self.line()
        if not line:
            raise EOFError
        if 'Values:' not in line:
            sys.exit('Cannot read density values from CELL2 file.')
        
        # parse array dimensions 
        line = self.line()
        if not line:
            raise EOFError
        ncol, nc = [int(var) for var in line.split()]

        # ignore next line
        line = self.line()
        
        # read points
        dens = []
        for i in range(nc):
            line = self.line()
            dens.append([float(var) for var in line.split()])
        dens = np.array(dens)
    
        # return columns as tuple of arrays with shape (nx,ny)
        return [dens[:,icol] for icol in range(ncol)]


class SurfaceFileReader(FileReader):
    """
    Manage reading from the surface file.
    """
    
    def fill_factors(self):
        """
        Read fill factors from the surface file.
        
        This is no longer used since fill factors are no longer stored on 
        SURF files.
        """
        # read array dimensions 
        line = self.line()
        nx, ny = [int(var) for var in line.split()]
    
        # read points
        points = []
        for i in range(nx*ny):
            line = self.line()
            points.append([float(var) for var in line.split()])
        points = np.array(points)
    
        # return columns as tuple of arrays with shape (nx,ny)
        return [points[:,i].reshape((nx,ny)) for i in range(points.shape[1])]
    
    def contours(self):
        """
        Read a contour from the surface file.
        """
        # read number of connected components
        line = self.line()
        ncc = int(line)
    
        for icc in range(ncc):
            # read number of points
            line = self.line()
            n = int(line)
    
            # read points
            points = []
            for i in range(n):
                line = self.line()
                points.append([float(var) for var in line.split()])
            points = np.array(points)
    
            # return columns as tuple of size np
            yield points[:,0], points[:,1]


class GridFileReader(FileReader):
    """
    Manage reading from the grid file.
    """

    def __init__(self, filename):
        FileReader.__init__(self, os.path.splitext(filename)[0] + '.grid')
    
    def grid_lines(self):
        """
        Read grid lines from the grid file.
        """
        # read header
        line = self.line()
        nxy = int(line.split()[1])
        line = self.line()
        
        # read points
        points = []
        for i in range(nxy):
            line = self.line()
            points.append([float(var) for var in line.split()])
        points = np.array(points)
        
        # return x- and y-values as tuples
        return points[:,0], points[:,1]


class Plotter(object):
    """
    Manage reading from file and displaying data in the plots.
    
    Common code for specialized plotters
    """
    
    def __init__(self, fig, ax, cax, options):
        self.fig = fig
        self.ax = ax
        self.cax = cax
        self.first = options.first
        self.last = options.last
        self.skip = options.skip
        self.fps = options.fps
        self.title = options.title
        self.scalebar = options.scalebar
        self.text = options.text
        self.frameon = options.frameon
        self.xmin = options.xmin
        self.xmax = options.xmax
        self.ymin = options.ymin
        self.ymax = options.ymax
        self.vmin = options.vmin
        self.vmax = options.vmax
        self.nm = options.nm
        self.record = 0
        self.plot_colormap = not options.nocolormap
        self.plot_isolines = options.isolines
        if options.gray:
            self.cmap = colors.LinearSegmentedColormap.from_list('inverted-gray', 
                                                                 ['white', 'black'])
        else:
            self.cmap = 'jet'
        self.plot_grid = not options.nogrid
        self.plot_oldgrid = not options.nooldgrid
        self.plot_surface = not options.nosurface
        self.plot_oldsurface = not options.nooldsurface
        
    def next_record(self, event):
        """
        Read and plot next data record.
        """
        # determine number of records to be skipped
        if event:
            if event.key in '0123456789':
                nskip = int(event.key)
            elif event.key == 'a' or event.key == 'A':
                nskip = 15
            elif event.key == 'b' or event.key == 'B':
                nskip = 31
            elif event.key == 'c' or event.key == 'C':
                nskip = 63
            elif event.key == 'd' or event.key == 'D':
                nskip = 127
            elif event.key == 'e' or event.key == 'E':
                nskip = 255
            elif event.key == 'f' or event.key == 'F':
                nskip = 511
            elif event.key == 'g' or event.key == 'G':
                nskip = 1023
            elif event.key == 'h' or event.key == 'H':
                nskip = 2047
            elif event.key == 'i' or event.key == 'I':
                nskip = 4095
            elif event.key == 'j' or event.key == 'J':
                nskip = 8191
            elif event.key == 'k' or event.key == 'K':
                nskip = 16383
            else:
                nskip = 0
        else:
            nskip = 0

        try:
            self.plot_record(nskip)
        except EOFError:
            print('No more data.')

        # update plot
        plt.draw()
        
    def make_movie(self, moviefile):
        """
        Make a movie from several records.
        """
        # check whether temporary directory already exists and create a new one
        if os.path.exists('_tmp'):
            prompt = '"_tmp" already exists. Delete it and all its contents? (y/n) '
            answer = input(prompt)
            if answer == 'y':
                shutil.rmtree('_tmp')
            else:
                return
        os.mkdir('_tmp')

        # first plot
        irec = self.first
        try:
            self.plot_record(irec-1)
        except EOFError:
            print("Could not read first record of file.")
            return
        
        # plot files and read other files
        while True:
        
            # save plot
            fname = os.path.join('_tmp', '%04d.png'%irec)
            self.fig.savefig(fname)
            
            # next plot
            try:
                self.plot_record(self.skip)
            except EOFError:
                break

            irec += self.skip + 1
            if irec >= self.last:
                break
        
        # make the movie out of the snapshot files
        print("Making movie " + moviefile + "...")
        os.system("mencoder 'mf://_tmp/*.png' -mf type=png:fps=%d " % self.fps +
              "-ovc lavc -lavcopts vcodec=wmv2 -oac copy -o " + moviefile)
        
        # clean up temporary files
        shutil.rmtree('_tmp')

    def make_sequence(self):
        """
        Plot a sequence of records.
        """
        # first plot
        irec = self.first
        try:
            self.plot_record(irec-1)
        except EOFError:
            print("Could not read first record of file.")
            return
        
        # plot other records
        while True:
        
            # next plot
            try:
                self.plot_record(self.skip, clear=False)
            except EOFError:
                break

            irec += self.skip + 1
            if irec >= self.last:
                break
        
    def close(self):
        self.file_reader.close()

    def _center_to_box(self, xc):
        """
        Return the box boundaries given the box centers.
        """
        # create array with size 1 larger in first dimension
        x1 = np.zeros((xc.shape[0]+1, xc.shape[1]))
        x1[0,:] = 1.5*xc[0,:] - 0.5*xc[1,:]
        x1[1:-1,:] = 0.5*(xc[:-1,:]+xc[1:,:])
        x1[-1,:] = 1.5*xc[-1,:] - 0.5*xc[-2,:]
    
        # create array with size 1 larger also in second dimension
        x = np.zeros((x1.shape[0], x1.shape[1]+1))
        x[:,0] = 1.5*x1[:,0] - 0.5*x1[:,1]
        x[:,1:-1] = 0.5*(x1[:,:-1]+x1[:,1:])
        x[:,-1] = 1.5*x1[:,-1] - 0.5*x1[:,-2]
    
        return x
    
    def _plot_text(self):
        """
        Plot text.
        """
        text, loc = self.text

        try:
            loc = int(loc)
        except:
            print('Invalid text location')
            loc = 3
            
        from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
        at = AnchoredText(text, loc, pad=0.5, borderpad=0.5, prop=dict(size=16), 
                        frameon=self.frameon)
    
        self.ax.add_artist(at)

    def _plot_scalebar(self):
        """
        Plot a scale bar.
        """
#        barx, bary, barlength = self.scalebar
#        ebarerr = 0.5*barlength
#        ebarx = barx + ebarerr
#        ebary = bary
#        if self.nm:
#            ebartext="%dnm"%barlength
#        else:
#            ebartext="%dA"%barlength
#
#        self.ax.errorbar(ebarx, ebary, xerr=ebarerr, ecolor="black", linewidth=3, mew=3, 
#                        capsize=6)
#        offset = mtransforms.offset_copy(self.ax.transData, fig=self.fig, x=0, y=10, 
#                                        units='points')
#        self.ax.text(ebarx, ebary, ebartext, horizontalalignment="center", 
#                    verticalalignment="bottom", transform=offset)
        
        from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar

        barlength, loc = self.scalebar
        if self.nm:
            bartext="%dnm"%barlength
        else:
            bartext="%dA"%barlength
        
        try:
            loc = int(loc)
        except:
            print('Invalid scalebar location')
            loc = 4

        asb =  AnchoredSizeBar(self.ax.transData, 
                            barlength, 
                            bartext,
                             loc=int(loc), 
                             pad=0.5, borderpad=0.5, sep=5, 
                             frameon=self.frameon)
        self.ax.add_artist(asb)


class Cell2Plotter(Plotter):
    """
    Manage reading from surface file and displaying data in the plots.
    """
    
    def __init__(self, filename, fig, ax, cax, options):
        Plotter.__init__(self, fig, ax, cax, options)
        self.ielem = options.variable
        if options.variable2:
            self.ielem2 = abs(options.variable2)
            self.ielem2_sign = math.copysign(1.0, options.variable2)
        else:
            self.ielem2 = None
            self.ielem2_sign = None
        self.cell_file_reader = Cell2FileReader(filename)
        surface_filename = os.path.splitext(filename)[0] + '.surf'
        self.surface_file_reader = SurfaceFileReader(surface_filename)
        self.old_cells = None
        self.cells = None
        self.old_surface = None
        self.surface = None
        
    def plot_record(self, nskip, clear=True):
        """
        Plot record with given index counted from current file position.
        """
        self.old_cells = deepcopy(self.cells)
        self.old_surface = deepcopy(self.surface)
        
        # read data
        for irec in range(nskip+1):
            
            self.record += 1

            # read grid and densities
            try:
                self.dose = self.cell_file_reader.dose()
                x, y = self.cell_file_reader.grid()
                dens = self.cell_file_reader.densities()
                if self.ielem2:
                    dens = dens[self.ielem] + self.ielem2_sign * dens[self.ielem2]
                else:
                    dens = dens[self.ielem]
                nx, ny = x.shape
                self.cells = (x, y, dens.reshape(nx-1,ny-1))
            except EOFError:
                if irec == 0:
                    raise EOFError
                else:
                    print('Tried to read beyond end of CELL2 file.')
                    break
            
            # read surface
            try:
                self.dose = self.surface_file_reader.dose()
                self.surface = []
                for x, y in self.surface_file_reader.contours():
                    self.surface.append((x, -y))
            except EOFError:
                print('No more surface.')

        # report current record and dose
        print('Record ', self.record, ', Dose ', self.dose)
        
        # clear figure
        if clear:
            self.ax.clear()
            if self.plot_colormap:
                self.cax.clear()
        
        # title
        if self.title:
            self.ax.set_title(self.title)
        else:
            #self.ax.set_title('record=%d' % self.record)
            #self.ax.set_title('%d ions/nm' % (self.dose*1e-14))
            self.ax.set_title('Dose=%.2e/cm2' % self.dose)

        # x- and y-values
        if self.nm:
            x = x / 10.
            y = y / 10.

        if self.scalebar:
            # indicate length scale by a bar
            self._plot_scalebar()
            self.ax.xaxis.set_visible(False)
            self.ax.yaxis.set_visible(False)
        else:
            # indicate length scale by axes
            if self.nm:
                self.ax.set_xlabel('lateral [nm]')
                self.ax.set_ylabel('vertical [nm]')
            else:
                self.ax.set_xlabel('lateral [A]')
                self.ax.set_ylabel('vertical [A]')
        
        if self.text:
            self._plot_text()

        # plot isolines
        if self.plot_isolines:
            if self.cells is not None:
                x, y, values = self.cells
                if self.nm:
                    x = x / 10.
                    y = y / 10.
                if self.vmin is not None and self.vmax is not None:
                    ISO = self.ax.contour(x, -y, values, 
                                          levels=np.arange(self.vmin, self.vmax, 
                                                           (self.vmax-self.vmin)/5), 
                                          colors='black')
                else:
                    ISO = self.ax.contour(x, -y, values, colors='black')
                self.ax.clabel(ISO, inline=1, fontsize=20)
            if self.old_cells is not None:
                x, y, values = self.old_cells
                if self.nm:
                    x = x / 10.
                    y = y / 10.
                if self.vmin is not None and self.vmax is not None:
                    self.ax.contour(x, -y, values, 
                                    levels=np.arange(self.vmin, self.vmax, 
                                                     (self.vmax-self.vmin)/5), 
                                    colors='black', linestyle='dotted')
                else:
                    self.ax.contour(x, -y, values, colors='black', linestyle='dotted')

        # plot densities
        if self.plot_colormap:
            if self.cells is not None:
                if self.plot_grid:
                    edgecolors = 'r'
                else:
                    edgecolors = 'none'
                x, y, values = self.cells
                if self.nm:
                    x = x / 10.
                    y = y / 10.
                c = self.ax.pcolor(x, -y, values, vmin=self.vmin, vmax=self.vmax, 
                                   cmap=self.cmap, edgecolors=edgecolors)#, linewidth=1)

            if self.old_cells is not None and self.plot_oldgrid:
                x, y, values = self.old_cells
                if self.nm:
                    x = x / 10.
                    y = y / 10.
                self.ax.pcolor(x, -y, values, vmin=self.vmin, vmax=self.vmax, 
                               cmap=self.cmap, facecolors='none', linestyle='dotted',
                               edgecolors='r')
                
            # draw color bar
            cbar = self.fig.colorbar(c, cax=self.cax)
            if self.plot_isolines:
                cbar.add_lines(ISO)

        # plot displacements
        if False and self.old_cells is not None:
            x, y, values = self.cells
            x_old, y_old, values = self.old_cells
            ux = x - x_old
            uy = y - y_old
            if self.nm:
                x = x / 10.
                y = y / 10.
                ux = ux / 10.
                uy = uy / 10.
#            self.ax.quiver(x[::2,::2], -y[::2,::2], ux[::2,::2], -uy[::2,::2], 
            self.ax.quiver(x, -y, ux, -uy, np.sqrt(ux**2+uy**2), pivot='middle')#,
                           #units='xy', scale=0.1)
            
        # plot surface
        if self.surface is not None and self.plot_surface:
            for x, y in self.surface:
                if self.nm:
                    x = x / 10.
                    y = y / 10.
                self.ax.plot(x, y, color='black', linewidth=2)#, label='visco')
        if self.old_surface is not None and self.plot_oldsurface:
            for x, y in self.old_surface:
                if self.nm:
                    x = x / 10.
                    y = y / 10.
                self.ax.plot(x, y, color='black', linestyle='dotted')
                #self.ax.plot(x, y, color='white', linestyle='dashed', linewidth=2, label='no visco')

        #legend = self.ax.legend(loc='right')
        #frame  = legend.get_frame()
        #frame.set_facecolor('0.80')
        
        # set axis limits (must do this last!)
        self.ax.set_xlim(self.xmin, self.xmax)
        self.ax.set_ylim(self.ymin, self.ymax)


class SurfacePlotter(Plotter):
    """
    Manage reading from surface file and displaying data in the plots.
    """
    
    def __init__(self, filename, fig, ax, cax, options):
        Plotter.__init__(self, fig, ax, cax, options)
        self.surface_reader = SurfaceFileReader(filename)
        self.old_contour = None
        self.new_contour = None
        self.grid_reader = GridFileReader(filename)
        
    def plot_record(self, nskip, clear=True):
        """
        Plot record with given index counted from current file position.
        """
        # read surface data
        for irec in range(nskip+1):
            
            # check if EOF has been reached; assume first line is a comment
            line = self.surface_reader.line()
            if not line:
                if irec == 0:
                    raise EOFError
                else:
                    print('Tried to read beyond end of surface file.')
                    break
            self.record += 1

            # read fill factors
#            ff = self.surface_reader.fill_factors()

            if self.record % 2 == 0:
                # save old contour
                self.old_contour = deepcopy(self.new_contour)
    
                # read contour
                self.new_contour = []
                for x, y in self.surface_reader.contours():
                    self.new_contour.append((x, -y))
            else:
                # read contour and ignore
                for x, y in self.surface_reader.contours():
                    pass

#                # read grid data
#                if self.grid:
#                    xgrid, ygrid = self.grid_reader.grid_lines()

        # report current record
        print('Record ', self.record)
        
        # clear figure
        if clear:
            self.ax.clear()
            if self.plot_colormap:
                self.cax.clear()
        
#        # plot fill factors
#        x, y, z = ff
#        x = self._center_to_box(x)
#        y = self._center_to_box(y)
#        if self.nm:
#            x = x / 10.
#            y = y / 10.
#            self.ax.set_xlabel('lateral [nm]')
#            self.ax.set_ylabel('vertical [nm]')
#        else:
#            self.ax.set_xlabel('lateral [A]')
#            self.ax.set_ylabel('vertical [A]')
#        
#        if self.plot_colormap:
#            c = self.ax.pcolor(x, -y, z, vmin=self.vmin, vmax=self.vmax, cmap=self.cmap)
#            # draw color bar
#            self.fig.colorbar(c, cax=self.cax)
#
#        # these quantities are needed to plot the grid...
#        x0 = x[0,0]
#        dx = x[1,0] - x[0,0]
#        y0 = y[0,0]
#        dy = y[0,1] - y[0,0]

        # plot old contour
        if self.old_contour is not None:
            for x, y in self.old_contour:
                if self.nm:
                    x = x / 10.
                    y = y / 10.
                self.ax.plot(x, y, linestyle='dotted', color='black')
    
        # plot new contour
        if self.new_contour is not None:
            for x, y in self.new_contour:
                if self.nm:
                    x = x / 10.
                    y = y / 10.
                self.ax.plot(x, y, color='black')
    
        # plot grid
#        if self.grid and self.record % 2 == 1:
##            print('x0=', x0)
##            print('dx=', dx)
##            print('xgrid=', xgrid)
#            x = x0 + dx*xgrid
#            y = y0 + dy*ygrid
#            self.ax.plot(x, -y, color='black')
            
    
        # set axis limits (must do this last!)
        self.ax.set_xlim(self.xmin, self.xmax)
        self.ax.set_ylim(self.ymin, self.ymax)


def parse_arguments(args):
    """
    Parse the arguments of the function or of the command line.
    """
    p = argparse.ArgumentParser( )
    p.add_argument('filename')
    p.add_argument('--movie', '-m', action='store_true')
    p.add_argument('--sequence', action='store_true')
    p.add_argument('--fps', type=float, default=10)
    p.add_argument('--first', type=int, default=1)
    p.add_argument('--last', type=int, default=10000)
    p.add_argument('--skip', type=int, default=0)
    p.add_argument('--isolines', action='store_true')
    p.add_argument('--nocolormap', action='store_true')
    p.add_argument('--title')
    p.add_argument('--scalebar', type=float, nargs=2)            # args: length, loc (1..10)
    p.add_argument('--text', nargs=2)                            # args: text, loc (1..10)
    p.add_argument('--frameon', action='store_true', default=False)
    p.add_argument('--xmin', type=float)
    p.add_argument('--xmax', type=float)
    p.add_argument('--ymin', type=float)
    p.add_argument('--ymax', type=float)
    p.add_argument('--vmin', type=float)
    p.add_argument('--vmax', type=float)
    p.add_argument('--nm', action='store_true')
    p.add_argument('--gray', action='store_true')
    p.add_argument('--variable', type=int, default=0)
    p.add_argument('--variable2', type=int, default=None)        # neg: subtract variable2
    p.add_argument('--nogrid', action='store_true')
    p.add_argument('--nooldgrid', action='store_true')
    p.add_argument('--nosurface', action='store_true')
    p.add_argument('--nooldsurface', action='store_true')
    
    print(args)
    options = p.parse_args(args)
    filename = options.filename
    return options, filename


def plot(data_type, args):
    """
    Plot the data given on file filename.
    """
    # check on valid data type
    if data_type not in ['surf', 'cell2']:
        sys.exit('Invalid data type: ' + data_type)

    # parse arguments
    options, filename = parse_arguments(args)
    print(options, data_type, filename)

    # open figure
    fig = plt.figure(figsize=(6,6))
    # open axes for main plot
    ax = fig.add_subplot(111, aspect='equal')

    # open axes for color bar
    if options.nocolormap:
        cax = None
    elif data_type == 'surf' or data_type == 'cell2':
#        fig.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#        cax = fig.add_axes([0.85, 0.1, 0.05, 0.8])
##        cax = fig.add_axes([0.75, 0.1, 0.05, 0.8])
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal(size="3%", pad="3%") 
        fig.add_axes(cax)       


    # create plotter
    if data_type == 'surf':
        plotter = SurfacePlotter(filename, fig, ax, cax, options)
    elif data_type == 'cell2':
        plotter = Cell2Plotter(filename, fig, ax, cax, options)
    
    if options.movie:
        # make a movie
        if options.variable2:
            plotter.make_movie(filename + '_%d'%options.variable + 
                               '%d'%options.variable2 + '.mpg')
        else:
            plotter.make_movie(filename + '_%d'%options.variable + '.mpg')
    elif options.sequence:
        # plot sequence of records
        plotter.make_sequence()
        plt.show()
    else:
        # plot first surface
        plotter.next_record(None)
    
        # connect key press event to plotter function and show plot
        fig.canvas.mpl_connect('key_press_event', plotter.next_record)
        plt.show()
    
    plotter.close()
    
if __name__ == '__main__':
#    plot('cell2', sys.argv[1:])
    plot(sys.argv[1], sys.argv[2:])
#    plot_surf('surf.dat')
else:
    print('plot.py loaded as a module')