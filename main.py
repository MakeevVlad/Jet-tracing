from utils import *
import plotter as my_plt
plt.style.use(astropy_mpl_style)

def main():
    #opening fits file
    names = ['data/J0121/J0121+1149_S_2012_12_05_pet_map.fits', 'data/J0121/J0121+1149_X_2012_12_05_pet_map.fits', 'data/J0121/J0121+1149_L_2010_05_21_pus_map.fits']
    names_2 = [ 'data/J0738/J0738+1742_L_2010_06_18_pus_map.fits', 'data/J0738/J0738+1742_L1_2017_08_11_pus_map.fits',
                'data/J0738/J0738+1742_L2_2017_08_11_pus_map.fits', 'data/J0738/J0738+1742_S_2017_08_11_pus_map.fits']
    names_3 = ['data/J1130/J1130+3815_C_1996_08_22_tay_map.fits', 'data/J1130/J1130+3815_K_2002_12_26_fey_map.fits', 'data/J1130/J1130+3815_Q_2002_12_26_fey_map.fits',
                 'data/J1130/J1130+3815_S_2002_05_08_pus_map.fits', 
                'data/J1130/J1130+3815_U_1999_01_02_moj_map.fits', 'data/J1130/J1130+3815_X_2002_05_08_pus_map.fits']
    
    blambas = [ 'J0249+2825_C_2015_10_30_pet_map', 'J0249+2825_X_2015_10_30_pet_map', 'J0259+0747_S_1997_01_10_fey_map',
                'J0259+0747_S_1997_01_10_pus_map', 'J0259+0747_S_2001_03_12_pus_map', 'J0259+0747_S_2003_07_09_nic_map',
                'J0259+0747_S_2003_07_09_pus_map', 'J0259+0747_S_2017_03_23_pet_map', 'J0259+0747_X_1997_01_10_pus_map',
                'J0351+7213_C_2015_11_09_pet_map', 'J0400+3247_C_2016_05_18_pet_map']
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
    
    is_blamba(bl_fits[0])

    #my_plt.plot(image_fits, plot_ridgline_data[1], plot_ridgline_data[2])# data=plot_ridgline_data[3])

    
def test():
    fig, (ax1, ax2) = plt.subplots(ncols = 2)
    arr1 = circle_clean(np.ones(shape = (50, 50)), 20, [15,29])
    arr2 = new_circle_clean(np.ones(shape = (50, 50)), 20, [15, 29])

    ax1.imshow(arr1, cmap = 'viridis')
    ax2.imshow(arr2 - arr1, cmap = 'viridis')

    ax2.set_title('New_cleaner')
    
    print(np.amin(arr2 - arr1))
    plt.show()
#test()

main()
