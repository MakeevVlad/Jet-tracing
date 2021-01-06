from utils import *
import plotter as my_plt

plt.style.use(astropy_mpl_style)

def main():
    #opening fits file
    names = ['data/J0121/J0121+1149_S_2012_12_05_pet_map.fits', 'data/J0121/J0121+1149_X_2012_12_05_pet_map.fits', 'data/J0121/J0121+1149_L_2010_05_21_pus_map.fits'] #3
    names_2 = [ 'data/J0738/J0738+1742_L_2010_06_18_pus_map.fits', 'data/J0738/J0738+1742_L1_2017_08_11_pus_map.fits',
                'data/J0738/J0738+1742_L2_2017_08_11_pus_map.fits', 'data/J0738/J0738+1742_S_2017_08_11_pus_map.fits', 'data/J0738/J0738+1742_C_2017_08_11_pus_map.fits'] #5
    names_3 = ['data/J1130/J1130+3815_C_1996_08_22_tay_map.fits', 'data/J1130/J1130+3815_K_2002_12_26_fey_map.fits', 'data/J1130/J1130+3815_Q_2002_12_26_fey_map.fits',
                 'data/J1130/J1130+3815_S_2002_05_08_pus_map.fits', 
                'data/J1130/J1130+3815_U_1999_01_02_moj_map.fits', 'data/J1130/J1130+3815_X_2002_05_08_pus_map.fits'] #6
    
    blambas = [ 'J0249+2825_C_2015_10_30_pet_map', 'J0249+2825_X_2015_10_30_pet_map', 'J0259+0747_S_1997_01_10_fey_map',
                'J0259+0747_S_1997_01_10_pus_map', 'J0259+0747_S_2001_03_12_pus_map', 'J0259+0747_S_2003_07_09_nic_map',
                'J0259+0747_S_2003_07_09_pus_map', 'J0259+0747_S_2017_03_23_pet_map', 'J0259+0747_X_1997_01_10_pus_map',
                'J0351+7213_C_2015_11_09_pet_map', 'J0400+3247_C_2016_05_18_pet_map'] #11
    
    for i in range(len(blambas)):
        blambas[i] = 'data/blumbsamples/' + blambas[i]
        blambas[i] += '.fits'

    
    #image_fits = fits.open(names_3[5])
    

    #plt.imshow(image_fits[0].data[0, 0], norm=colors.LogNorm())

    #plt.show()
    
    

    #plot_ridgline_data = ridgeline(image_fits)

    '''
    ridgelines = []
    fits_files = []
    for name in names_2:
        fits_files.append(fits.open(name))
        ridgelines.append(ridgeline(fits_files[-1])[1])

    my_plt.ridgelines_plot(fits_files, ridgelines)
    '''
    bl_fits = []
    for name in blambas:
        bl_fits.append(fits.open(name))
    fig, axs = plt.subplots(3, 2)
    i = 0
    for arr in [['/mnt/d/Source/rellab/Data/alldata/images/J0011+0823/J0011+0823_S_1995_07_15_yyk_map.fits']]: #names, names_2, names_3, blambas
        for name in arr:
            #print(name)
            print(is_blamba_test(fits.open(name), name, axs, i))
            i+=1
    plt.show()

    #my_plt.plot(image_fits, plot_ridgline_data[1], plot_ridgline_data[2])# data=plot_ridgline_data[3])

from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

def main__():
    print(os.getcwd())
    pool = ThreadPool(8)

    catalog = '/mnt/d/Source/rellab/Data/alldata/images'
    #maps = find_maps(catalog)
    #np.savetxt('fitsmaps_samples', maps, fmt='%s')
    maps = np.loadtxt('fitsmaps_samples', dtype = 'str')


    
    maps = np.copy(maps[:500])
    fts = np.array([fits.open(map) for map in maps])
    
    result = np.array(pool.map(is_blamba, fts))

    not_blambas = np.array(maps[result > 5.1])
    blambas = np.array(maps[result < 5.1])
    '''
    for map in maps:
        print(i / size(maps))
        if (is_blamba(fits.open(map))):
            blambas = np.append(blambas, map)
        else:
            not_blambas = np.append(not_blambas, map)
    '''
    np.savetxt('asimmetry', result, fmt='%s')
    np.savetxt('blamb_names', blambas, fmt='%s')
    np.savetxt('notblamb_names', not_blambas, fmt='%s')
    
def show_pics():
    '''
    maps = np.loadtxt('fitsmaps_samples', dtype = 'str')
    ass = np.loadtxt('asimmetry')
    good_ass = np.array([])
    good_maps = np.array([])
    bad_ass = np.array([])
    bad_maps = np.array([])
    #datas = np.array([fits.open(map)[0].data[0, 0] for map in maps])
    j = 1
    for i in np.random.permutation(250):
        mape = maps[i]
        asi = ass[i]
        print('#' + str(j))
        j += 1
        print(asi)
        plt.imshow(fits.open(mape)[0].data[0, 0],norm=colors.SymLogNorm(linthresh = np.abs(np.amin(fits.open(mape)[0].data[0, 0]))), cmap='jet', origin = 'lower')
        plt.show()
        ans = input()
        if ans == '1':
            good_ass = np.append(good_ass, asi)
            good_maps = np.append(good_maps, mape)
        else:
            bad_ass = np.append(bad_ass, asi)
            bad_maps = np.append(bad_maps, mape)            

    
    np.savetxt('bad_asi', bad_ass, fmt='%s')
    np.savetxt('good_asi', good_ass, fmt='%s')
    np.savetxt('good_maps', good_maps, fmt='%s')
    np.savetxt('bad_maps', bad_maps, fmt='%s')
    '''
    asimmertry_ = np.vectorize(asimmertry)
    maps = np.loadtxt('good_maps', dtype = 'str')
    ass = asimmertry_(np.array([fits.open(mape)[0].data[0, 0] for mape in maps]))
    asi =  np.loadtxt('good_asi')
    asi1 = np.loadtxt('bad_asi')
    mfr = np.vectorize(map_flux_ratio)
    
    plt.hist(asi, bins = 50, label= 'Asymm. of extended sources')
    plt.hist(asi1, bins = 30, label= 'Asymm. of blambas', color = 'red')
    plt.legend()
    #plt.hist(flux_ratio, bins = 50)
    plt.show()

def flux_ratios_dist():
    bad = np.array([fits.open(mape) for mape in np.loadtxt('bad_maps', dtype = 'str')])[0:10]
    good = np.array([fits.open(mape) for mape in np.loadtxt('good_maps', dtype = 'str')])[0:10]
    print(bad)
    pool = ThreadPool(8)
    data_good_gcleaned = np.array(pool.map(gauss_clean_adv__, good))
    data_bad_gcleaned = np.array(pool.map(gauss_clean_adv__, bad))
   
        

    



#flux_ratios_dist()
    
def test():
    '''
    fig, (ax1, ax2) = plt.subplots(ncols = 2)
    arr1 = circle_clean(np.ones(shape = (50, 50)), 20, [15,29])
    arr2 = new_circle_clean(np.ones(shape = (50, 50)), 20, [15, 29])

    ax1.imshow(arr1, cmap = 'viridis')
    ax2.imshow(arr2 - arr1, cmap = 'viridis')

    ax2.set_title('New_cleaner')
    
    print(np.amin(arr2 - arr1))
    plt.show()
    '''
    bad = np.array([fits.open(mape) for mape in np.loadtxt('bad_maps', dtype = 'str')[0:10]])
    #good = np.array([fits.open(mape) for mape in np.loadtxt('good_maps', dtype = 'str')])[0:150]

    for m in bad:
        my_plt.plot(m)

    
test()

#main()
