import numpy as np
from mayavi import mlab
import random as r
import os, sys
path = os.path.dirname(__file__)
print("path:", path)
sys.path.append(os.path.join(path, '..', '..', 'scripts'))

from read_output import read_his3


def bspot3d_fun(fname, atom=1):

    inp_file = os.path.join(path, os.path.splitext(fname)[0]+'.inp')
    f = open(inp_file)
    text = f.readline()
    text += 'Iso-surface 0.5, 1.5, 2.5, 3.5, 4.5 orders of magnitude below maximum'
    
    his3_file = os.path.join(path, fname)
    x, y, z, hist, _ = read_his3(his3_file, atom=atom)    # returns hist(z,y,x)
    hist = hist.transpose((2,1,0))               # now hist(x,y,z)
    hist_max = np.max(hist)
    hist = np.log10(np.maximum(hist_max/1e5, hist))
    hist_max = np.log10(hist_max)
    hist_max_xy = np.max(np.max(hist, axis=1), axis=0)
    hist_max_xz = np.max(np.max(hist, axis=2), axis=0)
    hist_max_yz = np.max(np.max(hist, axis=2), axis=1)
    xmin = x[hist_max_yz > hist_max-5][0]
    xmax = x[hist_max_yz > hist_max-5][-1]
    ymin = y[hist_max_xz > hist_max-5][0]
    ymax = y[hist_max_xz > hist_max-5][-1]
    zmin = z[hist_max_xy > hist_max-5][0]
    zmax = z[hist_max_xy > hist_max-5][-1]
    x, y, z = np.meshgrid(x, y, z, indexing='ij')
    xshift = 0.1*(xmax-xmin) - xmin
    extent = (xmin+xshift, xmax+xshift, ymin, ymax, zmin, zmax)
    mlab.contour3d(x+xshift, y, z, hist, opacity=0.6,
                   contours=[hist_max-4.5, hist_max-3.5, hist_max-2.5, 
                             hist_max-1.5, hist_max-0.5],
                   extent=extent) 
    mlab.outline(extent=extent)
    mlab.axes(extent=extent)
    mlab.xlabel(f'X+{xshift:.2f}')
    mlab.title(text, size=1, height=0.9)

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

    if oldpath:
        his3_file = os.path.join(oldpath, fname)
        x, y, z, hist, _ = read_his3(his3_file, atom=atom)    # returns hist(z,y,x)
        hist = hist.transpose((2,1,0))               # now hist(x,y,z)
        hist_max = np.max(hist)
        hist = np.log10(np.maximum(hist_max/1e5, hist))
        hist_max = np.log10(hist_max)
        x, y, z = np.meshgrid(x, y, z, indexing='ij')
        xshift = - 0.1*(xmax-xmin) - xmax
        extent = (xmin+xshift, xmax+xshift, ymin, ymax, zmin, zmax)
        mlab.contour3d(x+xshift, y, z, hist, opacity=0.6,
                       contours=[hist_max-4.5, hist_max-3.5, hist_max-2.5, 
                                 hist_max-1.5, hist_max-0.5],
                       extent=extent) 
        mlab.outline(extent=extent)
        mlab.axes(extent=extent)
        mlab.xlabel(f'X+{xshift:.2f}')
    else:
        print(warn)
    
    mlab.show()
