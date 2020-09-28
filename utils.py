import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
from matplotlib.patches import Ellipse
import numpy as np
import math
import scipy.optimize as so
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
import scipy.interpolate as spi




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

def map_center_ref(fits_file):
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

def cos_flat(v1, v2):
    return (v1[0]*v2[0] - v1[1]*v2[1]) / (magnitude_flat(v1) * magnitude_flat(v2))

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

    #plot(data_, trace, check_dots)
    #plot_1(fits_file, trace)

    return [data_, trace, check_dots]

#widt intergating
def trace_3(data):
    #data #= fits_file[0].data[0, 0]
    clean_up_np = np.vectorize(clean_up, otypes=[float])
    noise_lvl = map_noise(data)
    data_ = clean_up_np(data, noise_lvl, roughness=2, mode='default')

    trace = []
    check_dots = [[]] #For method vis.


    pos = np.array(map_center(data))
    pos_ = np.array([first_step(data_, pos=pos)[0], first_step(data_, pos=pos)[1]])

    trace.append(pos)
    trace.append(pos_)

    
    pos_tmp = np.array([])
    s = np.array([])
    p = np.array([])
    abs_s = 5 #step, default = 5
    dt = 1 #sep of int.
    for step_num in range(30):
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

    #plot_1(fits_file, trace)

    return [data_, trace, check_dots]

#computes data mean, data - array with data sets
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


#def merge_ridgline(fits_files):

    

#rescales image, doesnt change resolution, scale_factor between (0, 1)
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


def ridgeline(fits_file):
    data = fits_file[0].data[0, 0]

    clean_up_np = np.vectorize(clean_up, otypes=[float])
    noise_lvl = map_noise(data)
    data_ = clean_up_np(data, noise_lvl, roughness=2, mode='default')
    clean_data = np.copy(data_)

    

    mcenter = np.array(map_center(data))
    rline = [mcenter.copy()]
    rline_ = [mcenter.copy()]
    centers = mcenter.copy()

    
    check_dots = [mcenter.copy()]
  
    dr = 1 #radius step in px
    
    iterations = 80 #quantity of steps
    splits = 250 #circle fragmentation
    j_max = 8
    R0_koef = 1.5
    reflect_koef = 1 #now pointless
    R0 = R0_koef*circle_like_width(data_, mcenter)
    
    splits0 = splits

    j = 0
    phi = 0
    phi_prev = first_prev_phi(data_, rline_[0])
    for i in np.arange(0, iterations):
        R = j*dr + R0 #current radius
        j+=1
        
        
        
        accuracy_cur = 1
        I1 = 0
        I2 = 0


        if i == 0:
            center = [mcenter[0], mcenter[1]]

        try:
            for phi_tmp in np.linspace(0, 2*np.pi, splits, endpoint=False):
                I1 = 0
                I2 = 0
                ##print('Center')
                ##print(center)
                
                Alpha = np.linspace(phi_tmp - np.pi, phi_tmp, int(np.around(splits*0.5)), endpoint=False)
                
                inds =  np.array( [center[0] + np.around(R * np.sin(Alpha)),
                                                center[1] + np.around(R * np.cos(Alpha))], dtype = 'int16') 
                
                I1 = np.cumsum( data_[ inds[0], inds[1] ] )[-1]

                Alpha = np.linspace(phi_tmp, phi_tmp + np.pi, int(np.around(splits*0.5)), endpoint=False)

                inds =  np.array( [center[0] + np.around(R * np.sin(Alpha)),
                                                center[1] + np.around(R * np.cos(Alpha))], dtype = 'int16' ) 

                I2 =  np.cumsum( data_[ inds[0], inds[1] ] )[-1]
                '''
                for alpha in  np.linspace(phi_tmp - np.pi, phi_tmp, int(np.around(splits*0.5)), endpoint=False):
                    I1 += data_[center[0] + int(np.around(R * np.sin(alpha))), center[1] + int(np.around(R * np.cos(alpha)))]
                
                for alpha in  np.linspace(phi_tmp, phi_tmp + np.pi, int(np.around(splits*0.5)), endpoint=False):
                    I2 += data_[center[0] + int(np.around(R * np.sin(alpha))), center[1] + int(np.around(R * np.cos(alpha)))]
                '''
                if(I1 == 0 and I2 == 0):
                    break
                if(np.abs(I1-I2)/I2 <= accuracy_cur):
                    phi = phi_tmp
                    accuracy_cur = np.abs(I1-I2)/I2
        except IndexError:
            break

        print(i)
        print('phi_prev = ', phi_prev)
        print('phi = ', phi)
        
        if((I1 != 0 or I2 != 0) and I1 + I2 > 0.001):
    
            
            
            if np.cos(phi_prev - phi) < 0:
                phi += np.pi
                print('Reversed')

            rline_.append([center[0] + reflect_koef * int(np.around(R * np.sin(phi))), center[1] + reflect_koef * int(np.around(R * np.cos(phi)))])
            rline.append([center[1] + reflect_koef * int(np.around(R * np.cos(phi))), center[0] + reflect_koef * int(np.around(R * np.sin(phi)))])
            phi_prev = phi
            print(rline_[-1])
            print('I1 = ', I1)
            print('I2 = ', I2)

            if(j > j_max):
                print('New center')
                
                j = 0
                #splits = splits0
                data_ = circle_clean(data_, R-dr, center)
                

                center = rline_[-1]
                R0 = R0_koef*circle_like_width(data_, center)

                check_dots.append([center[1], center[0]])
                #centers.append(center)
            
            
        else:
            print('NNN')

        print('R = ', R)
        print('R0 = ', R0)
        print('====')

        if(stop_function(data, rline_[-1], noise_lvl)):
            print('R = ', R)
            print('R0 = ', R0)
            print('====')
            break

    return [data_, rline, check_dots, clean_data]


'''
def first_prev_phi(data, pos, dr=5, angles_num=180):
    alpha = np.arange(0, 2*np.pi, angles_num)
    r = 1
    prev_inds = np.array([], dtype = 'int16')
    while(1):
        inds = np.array([pos[0] + np.around(r * np.sin(alpha)), pos[1] + np.around(r * np.cos(alpha))], dtype = 'int16')
        if np.sum(data[inds[0], inds[1]] == 0) > len(inds[0]):
            return alpha[np.argmin( data[prev_inds[0], prev_inds[1]])]
        prev_inds = np.copy(inds)
        r+=dr
'''
def first_prev_phi(data, pos, dr=5, angles_num=180):
    max_alpha = 0
    max_len = 0
    for alpha in np.linspace(0, 2*np.pi, angles_num):
        r = 1
        while data[pos[0] + int(np.around(r * np.sin(alpha))), pos[1] + int(np.around(r * np.cos(alpha)))] != 0 :
            r+=dr
        if max_len < r:
            max_len = r
            max_alpha = alpha

    return max_alpha
    


def circle_clean(data, R, center, DR=5):
    '''
    for x in range(len(data)):
        for y in range(len(data)):
            if(magnitude_flat([center[0] - x, center[1] - y]) < R - 5):
                data[x, y] = 0
    '''
    for x in np.arange(center[0] - (R-DR), center[0] + (R-DR) + 1, dtype = 'int16'):
        y_ = int(np.floor(np.sqrt((R-DR)*(R-DR) - (x - center[0])*(x - center[0]))))
        y = np.arange( -y_ + center[1], y_ + center[1] + 1)
        data[x, y] = 0
    
    return data

def ocircle_clean(data, R, center, DR=5):
    for x in range(len(data)):
        for y in range(len(data)):
            if(magnitude_flat([center[0] - x, center[1] - y]) < R - DR):
                data[x, y] = 0
    
    return data


def circle_like_width(data, pos, dr=1, r0=1):
    r = r0
    while(1):
        alpha = np.linspace(0, 2*np.pi, 8*int(r))
        inds = np.array([pos[0] + np.around(r * np.sin(alpha)), pos[1] + np.around(r * np.cos(alpha))], dtype = 'int16')
        if np.sum(data[inds[0], inds[1]] == 0) > 0 :
            return 2*r
        r+=dr

def gauss_like_width(data, pos):
    return 0

def low_res_map_check(fits_file):
    return 1

def twoD_gauss(x, y, A, sigma_x, sigma_y, x_0, y_0, angle):
    x -= x_0
    y -= y_0
    angle*= np.pi/180
    angle -= 0.5*np.pi
    x_ = -(x)*np.cos(angle) + (y)*np.sin(angle)
    y_ = x*np.sin(angle) + y*np.cos(angle)
    return A * np.exp( -1 * ( ((x_) / sigma_x)**2 + ((y_) / sigma_y)**2 ) )

def oneD_gauss(x, A, sigma, x_0):
    return A * np.exp(-0.5 * ((x - x_0) / sigma)**2 )

def stop_function(data, pos, noise_lvl, roughness=2, size=0.05):
    x = np.arange( pos[0] - int(0.5*size*len(data)), pos[0] + int(0.5*size*len(data)), 1)
    y = np.arange( pos[1] - int(0.5*size*len(data)), pos[1] + int(0.5*size*len(data)), 1)
    xx, yy = np.meshgrid(x, y)
    
    
    integral = np.cumsum( data[xx, yy] )[-1]
    print('Integral = ', integral)
    if integral > noise_lvl*roughness*size*size:
        return 0
    
    else:
        return 1


def coversion_from_dec_to_polar_(pos, center=[0,0]):
    pos = [pos[0] - center[0], pos[1] - center[1]]

    return np.array([np.sqrt(pos[0]*pos[0] + pos[1]*pos[1]), arctan2(pos[0], pos[1])])

coversion_from_dec_to_polar = np.vectorize(coversion_from_dec_to_polar_)

def coversion_from_polar_to_dec_(pos, center=[0,0]):
    return np.array([pos[0] * np.sin(pos[1]), pos[0] * np.cos(pos[1])]) - np.array(center)
coversion_from_polar_to_dec = np.vectorize(coversion_from_polar_to_dec_)

def Arch_spiral_(x, a, b):
    return b*x + a
Arch_spiral = np.vectorize(Arch_spiral_)

def log_spiral_(x, a, b):
    return -1*a * (np.exp(b * x))
log_spiral = np.vectorize(log_spiral_)

def loglog_spiral_(x, a, b):
    return a*np.log(x) + b

def pol_spiral_(x, a, b, c, d, e, f):
    return np.dot(np.array([x**i for i in range(6)], np.array([a, b, c, d, e ,f])))
    
pol_spiral = np.vectorize(pol_spiral_)

def hyp_spiral_(x, a, b):
    return a/x + b

def arctan2(x, y):
    arc = np.arctan2(x, y)
    if arc < 0:
        arc += 2*np.pi
    return arc


def polar_fit(data, center=[0, 0], fit_function=Arch_spiral_):
    polar_data = []
    dec = data
    for i in data:
        polar_data.append(coversion_from_dec_to_polar_(i))
    
    #xdata = [a[1] for a in polar_data if a[0]>5]
    #ydata = [a[0] for a in polar_data if a[0]>5]
    xdata = [a[1] for a in polar_data]
    ydata = [a[0] for a in polar_data]
    popt, pcov = so.curve_fit(fit_function, xdata = xdata, ydata = ydata)

    phidata = np.linspace(0, 2*np.pi, 180)
    rdata = fit_function(phidata, *popt)
    new_data = []
    for i in range(len(phidata)):
        new_data.append(coversion_from_polar_to_dec_([rdata[i], phidata[[i]]]))

    '''
    xdata = np.array(xdata)
    z = np.polyfit(xdata,ydata,1)
    p = np.poly1d(z)
    #plt.scatter(x = xdata, y = p(xdata), color = 'blue', label = 'start')
    plt.scatter(x = xdata, y = ydata, color = 'yellow', label = 'train')

    plt.scatter(x = phidata, y = rdata, color = 'red', label = 'test')
    plt.legend()
    plt.show()
    '''
    
    return new_data

def polar_spline_fit(data, spline_deg=2):
    polar_data_ = []
    for i in data:
        polar_data_.append(coversion_from_dec_to_polar_(i))

    polar_data = np.array(polar_data_)
    polar_data = polar_data[np.argsort([m[1] for m in polar_data])]
    #polar_data = polar_data[np.array([m[1] for m in polar_data]) < 1.5]

    x_coef = 1 / np.amax([np.abs(m[1])  for m in polar_data])
    y_coef = 1 / np.amax([np.abs(m[0])  for m in polar_data])
    spline = spi.UnivariateSpline(x = [m[1] * x_coef  for m in polar_data] , y = [m[0] * y_coef for m in polar_data] , k = spline_deg)

    phi_data = np.linspace(np.amin([np.abs(m[1])  for m in polar_data]), (1 + 0.1)/x_coef)
    r_data = spline(phi_data * x_coef) / y_coef

    new_data = []
    for i in range(len(phi_data)):
        new_data.append(coversion_from_polar_to_dec_([r_data[i], phi_data[[i]]]))

    return new_data

def in_ellipse(x, y, center, a, b, p):
    p *= np.pi/180
    p = 0
    x -= center[0]
    y -= center[1]
    if (x * np.cos(p) + y*np.sin(p)) / a**2 + (y*np.cos(p) - x*np.sin(p)) / b**2 <= 1:
        return 1
    else:
        return 0


def ellipse_clean(data, center, a, b, angle, value=0):
    print(center)
    for x in np.arange(0, len(data), 1):
        for y in np.arange(0, len(data), 1):
            if in_ellipse(x, y, center, a, b, angle):
                data[x, y] = value
    return data


def gauss_clean(data, center, sx, sy, angle):
    Jm = np.amax(data)
    for x in np.arange(0, len(data), 1):
        for y in np.arange(0, len(data), 1):
            data[x, y] -= twoD_gauss(x, y, Jm, sx * 1.2, sy * 1.2, center[0], center[1], angle)
        print(x)
    return data

def is_blamba(fits_file):
    data = fits_file[0].data[0, 0]

    clean_up_np = np.vectorize(clean_up, otypes=[float])
    noise_lvl = map_noise(data)
    data_ = clean_up_np(data, noise_lvl, roughness=1, mode='default')
    clean_data = np.copy(data_)

    bmaj = fits_file[0].header['BMAJ']
    bmin = fits_file[0].header['BMIN']
    ra_inc = fits_file[0].header['CDELT1']
    dec_inc = fits_file[0].header['CDELT2']
    bpa = fits_file[0].header['BPA']
    center = map_center(data)

    fig, (ax1, ax2) = plt.subplots(1, 2)
    #
    #np.savetxt('dat_p', data)
    ax1.imshow(data, norm=colors.SymLogNorm(linthresh = np.abs(np.amin(data))), cmap='jet', origin = 'lower')
    data = gauss_clean(data, center, bmin/dec_inc, bmaj/ra_inc, bpa)
    np.savetxt('dat_a', data)
    #data = np.loadtxt('dat_a')
    ax2.imshow(data, norm=colors.SymLogNorm(linthresh = np.abs(np.amin(data))), cmap='jet', origin = 'lower')

    #x = np.linspace(0, 512, 100)
    #y = []
    #print(np.amax(data))
    #for i in x:
    #    y.append(twoD_gauss(i, center[1], np.amax(data), bmin/dec_inc * 1.2, bmaj/ra_inc * 1.2, center[0], center[1], bpa))
    #ax2.scatter(x, y)
    #ax2.scatter(np.arange(0, 512), data[np.arange(0, 512), center[1] ])
    #ax.add_patch(Ellipse((+len(data)*0.1, +len(data)*0.1),  
    #                4*bmin/dec_inc, 4*bmaj/ra_inc, fc = 'red', ec = 'black', lw = 1, angle = bpa))
    plt.show()

