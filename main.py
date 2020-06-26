from utils import *
plt.style.use(astropy_mpl_style)

def main():
    #opening fits file
    names = ['J0121+1149_S_2012_12_05_pet_map.fits', 'J0121+1149_X_2012_12_05_pet_map.fits', 'J0121+1149_L_2010_05_21_pus_map.fits']
    image_file = get_pkg_data_filename(names[0])
    fits.info(image_file)
    image_fits = fits.open(names[0])

    image_data = image_fits[0].data[0, 0]
    noise_lvl = map_noise(image_data)
    
    

    #for name in names:
     #   plot(fits.open(name)[0].data[0, 0])
    #merged_data = merge_1([fits.open(filename_) for filename_ in names])
    #trace_2(image_data)

    

    #plot(fits.open(names[0]))
    trace_3(image_data)


main()
