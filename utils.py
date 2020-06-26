import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
import numpy as np
import math
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS



def plot(data, trace=[[]], extra=[]):

    #print(trace)
    #print([ arr[0] for arr in trace ])


    fig, ax = plt.subplots()#num=None, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    #ax.imshow(data, norm=colors.LogNorm(), cmap=plt.cm.viridis)
    ax.imshow(data, cmap=plt.cm.viridis)
    if(trace != [[]]):
        plt.scatter(x = [ arr[0] for arr in trace ], y = [ arr[1] for arr in trace ], c='red', marker = 'x')
    
    if(extra != []):
        for extr in extra:
            plt.scatter(x = [ arr[0] for arr in extr ], y = [ arr[1] for arr in extr ], c='green', marker = '.', s = 4)
    plt.show()

def plot_1(fits_file, trace=[[]], extra=[]):

    data = np.array(fits_file[0].data[0, 0])
    '''
    object_name = fits_file[0].header('OBJECT')
    center_pix = [float(fits_file[0].header('CRPIX1')), float(fits_file[0].header('CRPIX2'))]
    center_deg = [(fits_file[0].header('RA')), (fits_file[0].header('DEC'))] #[RA, dec]

    ra_inc = float(fits_file[0].header('CDELT1'))
    dec_inc = float(fits_file[0].header('CDELT2'))
    '''
    #fig, ax = plt.subplots()#num=None, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    #ax.imshow(data, norm=colors.LogNorm(), cmap=plt.cm.viridis)
    hdu = fits_file[0]
    
    wcs = WCS(hdu.header)

    plt.subplot(projection=wcs)
    plt.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
    plt.grid(color='white', ls='solid')
    plt.xlabel('Galactic Longitude')
    plt.ylabel('Galactic Latitude')


    #ax.imshow(data, cmap=plt.cm.viridis)
    if(trace != [[]]):
        plt.scatter(x = [ arr[0] for arr in trace ], y = [ arr[1] for arr in trace ], c='red', marker = 'x')
    
    if(extra != []):
        for extr in extra:
            plt.scatter(x = [ arr[0] for arr in extr ], y = [ arr[1] for arr in extr ], c='green', marker = '.', s = 4)
    
    #ax.set_title(object_name)
    plt.show()



#simple jet-tracing (only for maps with high quality)
def trace_1(data):
    x = len(data)/2 - 1 #firs position coordinates (center of the object)
    y = len(data)/2 - 1
    r1 = 40 #inner radius of the search area
    r2 = 50 #outter radius of the search area
    counter = 0
    mxms = [[x, y]]
    while(check(counter) == 1):
        counter += 1
        max = 0
        mx = 0
        my = 0
        for i in range(len(data)):
            for j in range(len(data)):
                if (search_area(i, j, x, y, r1, r2, mxms[counter-2], counter)):
                    if (max <= data[i][j]):
                        max = data[i][j]
                        mx = i
                        my = j
        mxms.append([mx, my])
        x = mx
        y = my
        print(x, y)
    return mxms


#returns 1 if pixel is in a search area
def search_area(i, j, x, y, r1, r2, v, c):
    if ((x-i)**2 + (y-j)**2 >= r1**2 and (x-i)**2 + (y-j)**2 <= r2**2  and c <= 1):
        return 1
    elif((x-i)**2 + (y-j)**2 >= r1**2 and (x-i)**2 + (y-j)**2 <= r2**2 and cos([i - x, j - y], [x - v[0], y - v[1]]) > 0):
        return 1
    else:
        return 0

#returns cos value between two vectors
def cos(a, b):
    return (a[0]*b[0] + a[1]*b[1])/(np.sqrt((a[0]**2+a[1]**2)*(b[0]**2+b[1]**2)))


#just a counter
def check(i):
        if i < 10:
            return 1
        else:
            return 0

#returns standard deviation of map from map angle-located rectangles (area of 1 rectangle = area of the map * ratio)
def map_noise(data, ratio=0.124):
    data = np.array(data)
    trig = int(np.floor(len(data)*ratio))
    trigm = len(data) - 1
    s1 = np.std(data[:trig, :trig].ravel())
    s2 = np.std(data[:trig, trigm-trig:trigm+1].ravel())
    s3 = np.std(data[trigm-trig:trigm+1, :trig].ravel())
    s4 = np.std(data[trigm-trig:trigm+1, trigm-trig:trigm+1].ravel())
    return np.median([s1, s2, s3, s4])

#returns coordinates of the map centre (the brightest point)
def map_center(data):
    return np.unravel_index(np.argmax(data, axis=None), data.shape)

def map_center_ref(firs_file):
    return [float(fits_file[0].HEADER('CRPIX1')), float(fits_file[0].HEADER('CRPIX2'))]

#cleans your map from signals less then noise*rougness
def clean_up(value, noise, roughness=3, mode='default'):

    if (value > roughness*noise):
        if(mode=='default'):
            return value
        if(mode=='full'):
            return 1
    else:
        return 0

#return point the next step will be to
def first_step(data, pos, step=1, accuracy=1000, magnitude=10):
    dphi = 2*np.pi/accuracy

    dirs = np.array([[n*dphi, 0] for n in range(accuracy)])
    max = [0, 0]
    for dir in dirs:
        value = data[int(pos[0] + np.around(dir[1]*np.cos(dir[0])))][int(pos[1] + np.around(dir[1]*np.sin(dir[0])))]
        while (value > 0):
            dir[1]+=step
            value = data[int(pos[0] + np.around(dir[1]*np.cos(dir[0])))][int(pos[1] + np.around(dir[1]*np.sin(dir[0])))]
        if dir[1] >= max[1]:
            max = dir.copy()
    return [int(pos[0] + np.around(magnitude*np.sin(max[0]))), int(pos[1] + np.around(magnitude*np.cos(max[0])))]

#returns vector's magnitude
def magnitude_flat(v):
    return np.sqrt(v[0]**2 + v[1]**2)

def trace_2(data):
    #data #= fits_file[0].data[0, 0]
    clean_up_np = np.vectorize(clean_up, otypes=[float])
    noise_lvl = map_noise(data)
    data_ = clean_up_np(data, noise_lvl, roughness=2, mode='full')

    trace = []
    check_dots = [[]] #For method vis.


    pos = np.array(map_center(data))
    pos_ = np.array([first_step(data_, pos=pos)[0], first_step(data_, pos=pos)[1]])

    trace.append(pos)
    trace.append(pos_)

    
    pos_tmp = np.array([])
    s = np.array([])
    p = np.array([])
    abs_s = 5 #step
    dt = 1 #sep of int.
    for step_num in range(20):
        ## print("Stem number" , step_num) ##
        s = ( (pos_ - pos) * abs_s / magnitude_flat(pos_ - pos) ).copy()
        p_0 = ( s + pos_ ).copy()
        a = np.array([s[1]/abs_s, s[0]/abs_s*-1])
        n = 0
        #print(s)
        ## print(p_0) ##
        #print(a)
        edge_top = np.array([])
        edge_bot = np.array([])


        for sign in [-1, 1]:
            p = p_0.copy()

            while(data_[ int(np.around(p[1])) , int(np.around(p[0])) ] > 0):
                ## print(p) ##
                ## print([ int(np.around(p[0])) , int(np.around(p[1])) ]) ##
                ## print(data_[ int(np.around(p[1])) , int(np.around(p[0])) ]) ##
                p = ( p_0 + a*sign*n*dt ).copy()
                check_dots[0].append(p) ##
                n+=1
            n = 0
            if sign > 0 :
                edge_top = p.copy()
                ## print("look up, top") ##
            else:
                edge_bot = p.copy() 
                ## print("look up, bottom") ##
        ## print(edge_top) ##
        ## print(edge_bot) ##
        
        pos_tmp = np.around(0.5*(edge_top + edge_bot))
        pos = pos_.copy()
        pos_ = pos_tmp.copy()

        trace.append(pos_)

    #print(pos)
    #print(pos_)

    plot(data_, trace, check_dots)
    #plot_1(fits_file, trace)

    return 0

#widt intergating
def trace_3(data):
    #data #= fits_file[0].data[0, 0]
    clean_up_np = np.vectorize(clean_up, otypes=[float])
    noise_lvl = map_noise(data)
    data_ = clean_up_np(data, noise_lvl, roughness=2, mode='full')

    trace = []
    check_dots = [[]] #For method vis.


    pos = np.array(map_center(data))
    pos_ = np.array([first_step(data_, pos=pos)[0], first_step(data_, pos=pos)[1]])

    trace.append(pos)
    trace.append(pos_)

    
    pos_tmp = np.array([])
    s = np.array([])
    p = np.array([])
    abs_s = 5 #step
    dt = 1 #sep of int.
    for step_num in range(50):
        ## print("Stem number" , step_num) ##
        s = ( (pos_ - pos) * abs_s / magnitude_flat(pos_ - pos) ).copy()
        p_0 = ( s + pos_ ).copy()
        a = np.array([s[1]/abs_s, s[0]/abs_s*-1])
        n = 0
        #print(s)
        ## print(p_0) ##
        #print(a)
        edge_top = np.array([])
        edge_bot = np.array([])

        integral = np.array( p_0.copy() * data[int(np.around(p_0[1])), int(np.around(p_0[0]))] )
        int_sum = data[int(np.around(p_0[1])), int(np.around(p_0[0]))]


        for sign in [-1, 1]:
            p = p_0.copy()

            while(data_[ int(np.around(p[1])) , int(np.around(p[0])) ] > 0):
                ## print(p) ##
                ## print([ int(np.around(p[0])) , int(np.around(p[1])) ]) ##
                ## print(data_[ int(np.around(p[1])) , int(np.around(p[0])) ]) ##
                p = ( p_0 + a*sign*n*dt ).copy()

                integral += p * data[ int(np.around(p[1])) , int(np.around(p[0])) ]
                int_sum += data[ int(np.around(p[1])) , int(np.around(p[0])) ]

                check_dots[0].append(p) ##
                n+=1
            n = 0
        
        pos_tmp = np.around(integral / int_sum)
        pos = pos_.copy()
        pos_ = pos_tmp.copy()

        trace.append(pos_)

    #print(pos)
    #print(pos_)

    plot(data_, trace, check_dots)
    #plot_1(fits_file, trace)

    return 0

#computes data mean, data - array whith data sets
def merge_1(fits_files):
    #Array with scales values for every map
    scale = np.array([ abs(float(dat[0].header['CDELT2'])) for dat in fits_files ])

    scale_factor = np.array([ np.amin(scale) / s for s in scale ])

    #Array with data for every map
    data =  np.array([ dat[0].data[0, 0] for dat in fits_files ])

    rsc_data = np.array([ rescale_1(data[i], scale_factor[i]) for i in range(len(data)) ])

    #Array with resolution for every map
    resol = np.array([len(dat) for dat in data])

    #Array with centre coordianates (in pixels)
    center = np.array([ map_center(dat) for dat in data ])

    for i in range(1, len(data)):
        rsc_data[0] += rsc_data[i]
    
    return rsc_data[0] / len(rsc_data)




#rescales image, doesnt change resolution, scale factor_between (0, 1)
def rescale_1(data, scale_factor):
    shape_data_0 = [ int(np.around(np.shape(data)[0] * scale_factor)), int(np.around(np.shape(data)[1] * scale_factor)) ]
    data_0 = np.full(shape_data_0, map_noise(data))
    data_  = np.full(np.shape(data), map_noise(data))
    center = map_center(data)

    for x in range(len(data_0)):
        for y in range(len(data_0)):
            data_0[x, y] = data[int(round(x/scale_factor)), int(round(y/scale_factor))]
    

    center_0 = map_center(data_0)
    for x in range(len(data_0)):
        for y in range(len(data_0)):
            data[center[0] - center_0[0] + x , center[1] - center_0[1] + y] = data_0[x , y]
    
    return(data)


    
