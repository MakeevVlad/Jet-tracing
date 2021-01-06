import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDrawingArea
import numpy as np
import scipy.interpolate as spi
import math
from astropy.visualization import (astropy_mpl_style, SqrtStretch, ImageNormalize)
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
from utils import *


def hms_pixels_to_deg(pix_num, inc, center_pix, center_deg):
    h_c = int(center_deg.split(':')[0])
    m_c = int(center_deg.split(':')[1])
    s_c = float(center_deg.split(':')[2])

    h = h_c
    m = m_c
    s = (pix_num - center_pix)*inc*3600 + s_c
    while(s > 59):
        s-=60
        m+=1   
    while(m > 59):
        m-=60
        h+=1

    return '%(h)d:%(m)d:%(s)s' % {"h": h, "m" : m, "s": str(round(s, 3))}

def pixels_to_deg(pix_num, inc, center_pix, center_deg):

    h = 0
    m = 0
    s = ((pix_num - center_pix)*inc + center_deg)*3600
    while(s > 59):
        s-=60
        m+=1   
    while(m > 59):
        m-=60
        h+=1

    return '%(h)d:%(m)d:%(s)s' % {"h": h, "m" : m, "s": str(round(s, 3))}



def plot_simple(data, trace=[[]], extra=[]):

    #print(trace)
    #print([ arr[0] for arr in trace ])


    fig, ax = plt.subplots()#num=None, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    ax.imshow(data, norm=colors.LogNorm(), cmap=plt.cm.viridis)
    #ax.imshow(data, cmap=plt.cm.viridis)
    if(trace != [[]]):
        plt.scatter(x = [ arr[0] for arr in trace ], y = [ arr[1] for arr in trace ], c='red', marker = 'x')
    
    if(extra != []):
        for extr in extra:
            plt.scatter(x = [ arr[0] for arr in extr ], y = [ arr[1] for arr in extr ], c='green', marker = '.', s = 4)
    plt.show()

def plot_in_orb_coords(fits_file, trace=[[]], extra=[]):

    data = np.array(fits_file[0].data[0, 0]) 
    #data += abs(np.amin(data))
    
    object_name = fits_file[0].header['OBJECT']
    center_pix = [float(fits_file[0].header['CRPIX1']), float(fits_file[0].header['CRPIX2'])]
    #center_hms = [str(fits_file[0].header['RA']), str(fits_file[0].header['DEC'])] #[RA, dec]
    center_deg = [fits_file[0].header['CRVAL1'], (fits_file[0].header['CRVAL2'])] #[RA, dec]

    ra_inc = float(fits_file[0].header['CDELT1'])
    dec_inc = float(fits_file[0].header['CDELT2'])
    
    fig, ax = plt.subplots()#num=None, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    
    #ax.imshow(data, norm=colors.LogNorm(), cmap='viridis')
    x = np.arrange()
    ax.contour()

    if(trace != [[]]):
        plt.scatter(x = [ arr[0] for arr in trace ], y = [ arr[1] for arr in trace ], c='red', marker = 'x')
    
    if(extra != []):
        for extr in extra:
            plt.scatter(x = [ arr[0] for arr in extr ], y = [ arr[1] for arr in extr ], c='green', marker = '.', s = 4)
    
    ax.set_title(object_name)
    plt.show()

def plot(fits_file, trace=[[]], extra=[], data=[]):
    #clean_up_np = np.vectorize(clean_up, otypes=[float])
    if data == []:
        data = np.array(fits_file[0].data[0, 0]) 

    
    object_name = fits_file[0].header['OBJECT']
    center_pix = [float(fits_file[0].header['CRPIX1']), float(fits_file[0].header['CRPIX2'])]
    #center_hms = [str(fits_file[0].header['RA']), str(fits_file[0].header['DEC'])] #[RA, dec]
    center_deg = [fits_file[0].header['CRVAL1'], (fits_file[0].header['CRVAL2'])] #[RA, dec]
    
    ra_inc = float(fits_file[0].header['CDELT1'])*3600*1000
    dec_inc = float(fits_file[0].header['CDELT2'])*3600*1000

    bmaj = float(fits_file[0].header['BMAJ'])*3600*1000
    bmin = float(fits_file[0].header['BMIN'])*3600*1000

    bpa = float(fits_file[0].header['BPA'])
    
    fig, ax = plt.subplots()#num=None, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    extent = [-center_pix[0]*(ra_inc), (len(data) - center_pix[0])*(ra_inc), 
                -center_pix[1]*(dec_inc), (len(data) - center_pix[1])*(dec_inc)]

    #ax.imshow(data, norm=colors.LogNorm(), cmap='viridis', origin = 'lower', extent = extent)
    ax.imshow(data, norm=colors.SymLogNorm(linthresh = np.abs(np.amin(data))), cmap='jet', origin = 'lower', extent = extent)
    levels1 = np.logspace(0, np.amax(data), num = 8, base = 2)
    
    #ax.contour(data, -2, norm = colors.LogNorm())
    

    ax.set_xlabel('Relative R.A., [mas]')

    
    ax.set_ylabel('Relative Dec., [mas]')
    #ax.imshow(data, cmap=plt.cm.viridis)

    if(trace != [[]]):
        plt.scatter(x = [ arr[0] * (extent[1]-extent[0]) / len(data) + extent[0] for arr in trace ], 
                    y = [ - arr[1] * (extent[2]-extent[3]) / len(data) + extent[2] for arr in trace ], c='red', marker = 'x')
    
    if(extra != []):
        #for extr in extra:
        plt.scatter(x = [ arr[0] * (extent[1]-extent[0]) / len(data) + extent[0] for arr in extra ], 
                y = [ - arr[1] * (extent[2]-extent[3]) / len(data) + extent[2] for arr in extra ], c='green', marker = 'o')

    ax.add_patch(Ellipse((-len(data)*0.9*0.5*(ra_inc), -len(data)*0.9*0.5*(dec_inc)), 2*bmin, 2*bmaj, angle = -bpa,  fc = 'red', ec = 'black', lw = 1))

    ax.set_title(object_name)
    plt.show()


def ridgelines_plot(fits_files, ridgelines):
    data = []
    ra_inc = []
    dec_inc = []
    center_pix = []
    extent = []
    for fits_file in fits_files:
        data.append(fits_file[0].data[0, 0])

        ra_inc.append(float(fits_file[0].header['CDELT1'])*3600*1000)
        dec_inc.append(float(fits_file[0].header['CDELT2'])*3600*1000)

        center_pix.append([float(fits_file[0].header['CRPIX1']), float(fits_file[0].header['CRPIX2'])])

        extent.append([-center_pix[-1][0]*(ra_inc[-1]), (len(data[-1]) - center_pix[-1][0])*(ra_inc[-1]), 
                -center_pix[-1][1]*(dec_inc[-1]), (len(data[-1]) - center_pix[-1][1])*(dec_inc[-1])])

    fig, (ax, ax_pol) = plt.subplots(ncols = 2)
    ax.imshow(data[0],norm=colors.SymLogNorm(linthresh = np.abs(np.amin(data[0]))), cmap='jet', origin = 'lower', extent = extent[0])


    ext_data = []
    ext_data_sections =[]
    for i in range(len(ridgelines)):
        x = [ arr[0] * (extent[i][1]-extent[i][0]) / len(data[i]) + extent[i][0] for arr in ridgelines[i] ]
        y = [-arr[1] * (extent[i][2]-extent[i][3]) / len(data[i]) + extent[i][2] for arr in ridgelines[i] ]
        ax.scatter(x = x, y = y, marker = 'x')
        ext_data_sections.append([x, y])
        for j in range(len(ridgelines[i])):
            ext_data.append([x[j], y[j]])
    
    
    #spline_data = polar_spline_fit(ext_data, 3)
    
    #ax.plot([m[0] for m in spline_data], [m[1] for m in spline_data], c='blue')   

    ax.set_xlabel('Relative R.A., [mas]')
    object_name = fits_files[0][0].header['OBJECT']
    ax.set_title(object_name)
    ax.set_ylabel('Relative Dec., [mas]')

    
    polar_data = []
    for i in ext_data:
        polar_data.append(coversion_from_dec_to_polar_(i))



    polar_data_ = np.array(polar_data)

    polar_data_ = polar_data_[np.argsort([m[1] for m in polar_data])]
    #polar_data_ = polar_data_[np.array([m[1] for m in polar_data_]) < 1.5]
    

    #linsp = np.linspace(0, 2*np.pi, 100)
    #spline = spi.UnivariateSpline(x = [m[1] /6 for m in polar_data_] , y = [m[0] / 190 for m in polar_data_] , k = 2)
    
    #ax_pol.plot([m[1] for m in polar_spline_data], [m[0] for m in polar_spline_data], color = 'blue')
    ax_pol.scatter(x = [m[1] for m in polar_data_], y = [m[0] for m in polar_data_], color =  'red', marker ='x')
    
    '''
    for dat in ext_data_sections:
        r = []
        phi = []
        for i in range(len(dat)):
            r.append(coversion_from_dec_to_polar_([dat[1][i], dat[0][i]])[0])
            phi.append(coversion_from_dec_to_polar_([dat[1][i], dat[0][i]])[1]) 
        ax_pol.scatter(x=phi, y=r, marker='x')
    '''
    ax_pol.set_xlabel('Rotation angle, [rad]')
    ax_pol.set_ylabel('R, [mas]')


    plt.show()
    