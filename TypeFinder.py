def typing_kit(file_name) :

    '''
    Ellianna Schwab, Kelle Cruz

    typing_kit provides plots to qualitatively spectral type L-type brown dwarfs in NIR regime.
    Typing methods are taken from Cruz et al. 2017 and offer comparison to Kirkpatrick 2010 typing templates.

    To use input path to file_name in a string. A band-by-band grid of Cruz et al 2017 templates
    will be overplotted with the target spectrum. To compare to Kirkpatrick 2010, key in the spectral type
    number, '0' - '8'. Cruz et al. 2017 templates are shown band-by-band followed by Kirkpatrick 2010
    templates of the overall NIR spectrum. To save the plot, hit 's'.
    '''
    
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy.io import ascii

    hdulist = fits.open(file_name)
    spectrum = hdulist[0]
    wavelength = spectrum.data[0]
    flux = spectrum.data[1]
    
    fig1, axes1 = plt.subplots(9, 3, figsize=(6,11), sharey=True)

    #Create J, H and K bands for spectrum

    wavelength_J = []
    flux_J = []
    wavelength_H = []
    flux_H = []
    wavelength_K = []
    flux_K = []
    for jj in range(len(wavelength)):
        if wavelength[jj] >= 0.87 and wavelength[jj] <= 1.39:
            wavelength_J.append(wavelength[jj])
            flux_J.append(flux[jj])
        elif wavelength[jj] >= 1.41 and wavelength[jj] <= 1.89:
            wavelength_H.append(wavelength[jj])
            flux_H.append(flux[jj])
        elif wavelength[jj] >= 1.91 and wavelength[jj] <=2.39:
            wavelength_K.append(wavelength[jj])
            flux_K.append(flux[jj])


    #Normalize Each Band        

    wavelength_J = np.array(wavelength_J)
    flux_J = np.array(flux_J)
    wavelength_H = np.array(wavelength_H)
    flux_H = np.array(flux_H)
    wavelength_K = np.array(wavelength_K)
    flux_K = np.array(flux_K)

    flux_J = flux_J/np.mean(flux_J)
    flux_H = flux_H/np.mean(flux_H)
    flux_K = flux_K/np.mean(flux_K)


    #Create the plots

    for ii in range(9):
            ##J Band
            axes1[ii, 0].plot(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                              ascii.read("templates/L{}J_f.txt".format(ii))['col2'], c='red')
            axes1[ii, 0].fill_between(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                      ascii.read("templates/L{}J_f.txt".format(ii))['col4'], \
                                      ascii.read("templates/L{}J_f.txt".format(ii))['col5'], color='#c6c6c6')
            axes1[ii, 0].plot(wavelength_J, flux_J, c='k')
            axes1[ii, 0].annotate('L{}'.format(ii), xy=(0.1, 0.9), xycoords='axes fraction', color='k')
            axes1[ii, 0].axis('off')

            #H Band
            axes1[ii, 1].plot(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                              ascii.read("templates/L{}H_f.txt".format(ii))['col2'], c='red') 
            axes1[ii, 1].fill_between(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                      ascii.read("templates/L{}H_f.txt".format(ii))['col4'], \
                                      ascii.read("templates/L{}H_f.txt".format(ii))['col5'], color='#c6c6c6') 
            axes1[ii, 1].plot(wavelength_H, flux_H, c='k')
            axes1[ii, 1].axis('off')

            #K Band
            axes1[ii, 2].plot(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                              ascii.read("templates/L{}K_f.txt".format(ii))['col2'], c='red')
            axes1[ii, 2].fill_between(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                      ascii.read("templates/L{}K_f.txt".format(ii))['col4'], \
                                      ascii.read("templates/L{}K_f.txt".format(ii))['col5'], color='#c6c6c6')
            axes1[ii, 2].plot(wavelength_K, flux_K, c='k')
            axes1[ii, 2].axis('off')


    fig1.tight_layout()
    
    def on_key_press(event):
        hdulist = fits.open(file_name)
        spectrum = hdulist[0]
        wavelength = spectrum.data[0]
        flux = spectrum.data[1]


        #Create J, H and K bands for spectrum

        wavelength_J = []
        flux_J = []
        wavelength_H = []
        flux_H = []
        wavelength_K = []
        flux_K = []
        for jj in range(len(wavelength)):
            if wavelength[jj] >= 0.87 and wavelength[jj] <= 1.39:
                wavelength_J.append(wavelength[jj])
                flux_J.append(flux[jj])
            elif wavelength[jj] >= 1.41 and wavelength[jj] <= 1.89:
                wavelength_H.append(wavelength[jj])
                flux_H.append(flux[jj])
            elif wavelength[jj] >= 1.91 and wavelength[jj] <=2.39:
                wavelength_K.append(wavelength[jj])
                flux_K.append(flux[jj])


        #Normalize Each Band        

        wavelength_J = np.array(wavelength_J)
        flux_J = np.array(flux_J)
        wavelength_H = np.array(wavelength_H)
        flux_H = np.array(flux_H)
        wavelength_K = np.array(wavelength_K)
        flux_K = np.array(flux_K)

        flux_J = flux_J/np.mean(flux_J)
        flux_H = flux_H/np.mean(flux_H)
        flux_K = flux_K/np.mean(flux_K)

        
        #Normalize the Overall Spectra to Kirkpatrick Normalization
        norm_flux_array=[]
        for jj in range(len(wavelength)):
            if wavelength[jj] >= 1.28 and wavelength[jj] <= 1.39:
                norm_flux_array.append(flux[jj])
        norm_flux_array=np.array(norm_flux_array)

        wavelength=np.array(wavelength)
        flux=np.array(flux)
        norm_flux=flux/np.mean(norm_flux_array)


        #List the standards, so they can be called # "spectra/nir/prism_0835+1953_20050123_CHI06A.fits", this  one doesn't exist?
        NIR_standards = ["spectra/nir/U20165_0345+2540_new_w.fits", "spectra/nir/spex_prism_2130-0845_080713.fits", \
                         "spectra/nir/U10244_0408-1450.fits", "spectra/nir/u11291_1506+1321_050323.fits", \
                         "spectra/nir/U12101_2158-1550_davy.fits", "spectra/nir/spex_prism_2137+0808_U20909.fits", \
                         "spectra/nir/U10880.fits", "spectra/nir/u10721_050323.fits", "spectra/nir/2M1632.fits"]

        
        # If the key is in an easily bracketed part, or 1,2,3,4,5,6,7
        if event.key in '1234567':
            type_number = int(event.key)

            #Create the Plots
            fig2, axes2 = plt.subplots(
                nrows=3, ncols=4, sharex=False, sharey=False, 
                gridspec_kw={'width_ratios':[1,1,1,3]}
                )

            for jj, ii in zip([0,1,2], [type_number-1, type_number, type_number+1]):
                    ##J Band
                    axes2[jj, 0].plot(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                      ascii.read("templates/L{}J_f.txt".format(ii))['col2'], c='red')
                    axes2[jj, 0].fill_between(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                              ascii.read("templates/L{}J_f.txt".format(ii))['col4'], \
                                              ascii.read("templates/L{}J_f.txt".format(ii))['col5'], color='#c6c6c6')
                    axes2[jj, 0].plot(wavelength_J, flux_J, c='k')
                    axes2[jj, 0].annotate('L{}'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                    axes2[jj, 0].set_ylim([-0.5, 2])
                    axes2[jj, 0].axis('off')

                    #H Band
                    axes2[jj, 1].plot(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                      ascii.read("templates/L{}H_f.txt".format(ii))['col2'], c='red') 
                    axes2[jj, 1].fill_between(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                              ascii.read("templates/L{}H_f.txt".format(ii))['col4'], \
                                              ascii.read("templates/L{}H_f.txt".format(ii))['col5'], color='#c6c6c6') 
                    axes2[jj, 1].plot(wavelength_H, flux_H, c='k')
                    axes2[jj, 1].axis('off')


                    #K Band
                    axes2[jj, 2].plot(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                      ascii.read("templates/L{}K_f.txt".format(ii))['col2'], c='red')
                    axes2[jj, 2].fill_between(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                              ascii.read("templates/L{}K_f.txt".format(ii))['col4'], \
                                              ascii.read("templates/L{}K_f.txt".format(ii))['col5'], color='#c6c6c6')
                    axes2[jj, 2].plot(wavelength_K, flux_K, c='k')
                    axes2[jj, 2].axis('off')


                    #All Together Now! This is where the Kirkpatrick 10 one comes in.
                    temp_hdulist = fits.open(NIR_standards[ii])
                    temp_spectrum = temp_hdulist[0]
                    temp_wavelength = temp_spectrum.data[0]
                    temp_flux = temp_spectrum.data[1]

                    temp_norm_flux_array=[]
                    for kk in range(len(temp_wavelength)):
                        if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                            temp_norm_flux_array.append(temp_flux[kk])
                    temp_norm_flux_array=np.array(temp_norm_flux_array)

                    temp_wavelength=np.array(temp_wavelength)
                    temp_flux=np.array(temp_flux)
                    temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)


                    axes2[jj, 3].plot(temp_wavelength, temp_norm_flux, c='red')
                    axes2[jj, 3].plot(wavelength, norm_flux, c='k')
                    axes2[jj, 3].set_ylim([0, max(norm_flux) + 0.25])
                    axes2[jj, 3].set_xlim([min(wavelength), max(wavelength)])
                    axes2[jj, 3].axis('off')

        elif event.key == '0':
            type_number = int(event.key)
           
            #Create the Plots

            fig2, axes2 = plt.subplots(
                nrows=3, ncols=4, sharex=False, sharey=False, 
                gridspec_kw={'width_ratios':[1,1,1,3]}
                )
            
            #Begin with the M9 template and then iterate over the rest
            axes2[0, 0].annotate('M9', xy=(0.1, 0.8), xycoords='axes fraction', color='k')
            axes2[0, 0].axis('off')
            axes2[0, 1].axis('off')
            axes2[0, 2].axis('off')
            
            temp_hdulist = fits.open("spectra/M9V_LHS2924.fits")
            temp_spectrum = temp_hdulist[0]
            temp_wavelength = temp_spectrum.data[0]
            temp_flux = temp_spectrum.data[1]
            temp_norm_flux_array=[]
            for kk in range(len(temp_wavelength)):
                if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                    temp_norm_flux_array.append(temp_flux[kk])
            temp_norm_flux_array=np.array(temp_norm_flux_array)

            temp_wavelength=np.array(temp_wavelength)
            temp_flux=np.array(temp_flux)
            temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)
            axes2[0, 3].plot(temp_wavelength, temp_norm_flux, c='red')
            axes2[0, 3].plot(wavelength, norm_flux, c='k')
            axes2[0, 3].set_ylim([0, max(norm_flux) + 0.25])
            axes2[0, 3].set_xlim([min(wavelength), max(wavelength)])
            axes2[0, 3].axis('off')

            for jj, ii in zip([1,2], [type_number, type_number+1]):
                    ##J Band
                    axes2[jj, 0].plot(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                      ascii.read("templates/L{}J_f.txt".format(ii))['col2'], c='red')
                    axes2[jj, 0].fill_between(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                              ascii.read("templates/L{}J_f.txt".format(ii))['col4'], \
                                              ascii.read("templates/L{}J_f.txt".format(ii))['col5'], color='#c6c6c6')
                    axes2[jj, 0].plot(wavelength_J, flux_J, c='k')
                    axes2[jj, 0].annotate('L{}'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                    axes2[jj, 0].set_ylim([-0.5, 2])
                    axes2[jj, 0].axis('off')

                    #H Band
                    axes2[jj, 1].plot(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                      ascii.read("templates/L{}H_f.txt".format(ii))['col2'], c='red') 
                    axes2[jj, 1].fill_between(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                              ascii.read("templates/L{}H_f.txt".format(ii))['col4'], \
                                              ascii.read("templates/L{}H_f.txt".format(ii))['col5'], color='#c6c6c6') 
                    axes2[jj, 1].plot(wavelength_H, flux_H, c='k')
                    axes2[jj, 1].axis('off')


                    #K Band
                    axes2[jj, 2].plot(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                      ascii.read("templates/L{}K_f.txt".format(ii))['col2'], c='red')
                    axes2[jj, 2].fill_between(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                              ascii.read("templates/L{}K_f.txt".format(ii))['col4'], \
                                              ascii.read("templates/L{}K_f.txt".format(ii))['col5'], color='#c6c6c6')
                    axes2[jj, 2].plot(wavelength_K, flux_K, c='k')
                    axes2[jj, 2].axis('off')


                    #All Together Now! This is where the Kirkpatrick 10 one comes in.
                    temp_hdulist = fits.open(NIR_standards[ii])
                    temp_spectrum = temp_hdulist[0]
                    temp_wavelength = temp_spectrum.data[0]
                    temp_flux = temp_spectrum.data[1]

                    temp_norm_flux_array=[]
                    for kk in range(len(temp_wavelength)):
                        if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                            temp_norm_flux_array.append(temp_flux[kk])
                    temp_norm_flux_array=np.array(temp_norm_flux_array)

                    temp_wavelength=np.array(temp_wavelength)
                    temp_flux=np.array(temp_flux)
                    temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)


                    axes2[jj, 3].plot(temp_wavelength, temp_norm_flux, c='red')
                    axes2[jj, 3].plot(wavelength, norm_flux, c='k')
                    axes2[jj, 3].set_ylim([0, max(norm_flux) + 0.25])
                    axes2[jj, 3].set_xlim([min(wavelength), max(wavelength)])
                    axes2[jj, 3].axis('off')
                    
        elif event.key == '8':
            type_number = int(event.key)

            #List the standards, so they can be called # "spectra/nir/prism_0835+1953_20050123_CHI06A.fits", this  one doesn't exist?
            NIR_standards = ["spectra/nir/u10721_050323.fits", "spectra/nir/2M1632.fits"]

            #Create the Plots

            fig2, axes2 = plt.subplots(
                nrows=3, ncols=4, sharex=False, sharey=False, 
                gridspec_kw={'width_ratios':[1,1,1,3]}
                )

            for jj, ii in zip([0,1], [type_number-1, type_number]):
                    ##J Band
                    axes2[jj, 0].plot(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                      ascii.read("templates/L{}J_f.txt".format(ii))['col2'], c='red')
                    axes2[jj, 0].fill_between(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                              ascii.read("templates/L{}J_f.txt".format(ii))['col4'], \
                                              ascii.read("templates/L{}J_f.txt".format(ii))['col5'], color='#c6c6c6')
                    axes2[jj, 0].plot(wavelength_J, flux_J, c='k')
                    axes2[jj, 0].annotate('L{}'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                    axes2[jj, 0].set_ylim([-0.5, 2])
                    axes2[jj, 0].axis('off')

                    #H Band
                    axes2[jj, 1].plot(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                      ascii.read("templates/L{}H_f.txt".format(ii))['col2'], c='red') 
                    axes2[jj, 1].fill_between(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                              ascii.read("templates/L{}H_f.txt".format(ii))['col4'], \
                                              ascii.read("templates/L{}H_f.txt".format(ii))['col5'], color='#c6c6c6') 
                    axes2[jj, 1].plot(wavelength_H, flux_H, c='k')
                    axes2[jj, 1].axis('off')


                    #K Band
                    axes2[jj, 2].plot(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                      ascii.read("templates/L{}K_f.txt".format(ii))['col2'], c='red')
                    axes2[jj, 2].fill_between(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                              ascii.read("templates/L{}K_f.txt".format(ii))['col4'], \
                                              ascii.read("templates/L{}K_f.txt".format(ii))['col5'], color='#c6c6c6')
                    axes2[jj, 2].plot(wavelength_K, flux_K, c='k')
                    axes2[jj, 2].axis('off')


                    #All Together Now! This is where the Kirkpatrick 10 one comes in.
                    temp_hdulist = fits.open(NIR_standards[jj])
                    temp_spectrum = temp_hdulist[0]
                    temp_wavelength = temp_spectrum.data[0]
                    temp_flux = temp_spectrum.data[1]

                    temp_norm_flux_array=[]
                    for kk in range(len(temp_wavelength)):
                        if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                            temp_norm_flux_array.append(temp_flux[kk])
                    temp_norm_flux_array=np.array(temp_norm_flux_array)

                    temp_wavelength=np.array(temp_wavelength)
                    temp_flux=np.array(temp_flux)
                    temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)


                    axes2[jj, 3].plot(temp_wavelength, temp_norm_flux, c='red')
                    axes2[jj, 3].plot(wavelength, norm_flux, c='k')
                    axes2[jj, 3].set_ylim([0, max(norm_flux) + 0.25])
                    axes2[jj, 3].set_xlim([min(wavelength), max(wavelength)])
                    axes2[jj, 3].axis('off')
                    
            #Finish with T0 template
            axes2[2, 0].annotate('T0', xy=(0.1, 0.8), xycoords='axes fraction', color='k')
            axes2[2, 0].axis('off')
            axes2[2, 1].axis('off')
            axes2[2, 2].axis('off')
            
            temp_hdulist = fits.open("spectra/T0_2M0423_T0.fits")
            temp_spectrum = temp_hdulist[0]
            temp_wavelength = temp_spectrum.data[0]
            temp_flux = temp_spectrum.data[1]
            temp_norm_flux_array=[]
            for kk in range(len(temp_wavelength)):
                if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                    temp_norm_flux_array.append(temp_flux[kk])
            temp_norm_flux_array=np.array(temp_norm_flux_array) 
            temp_wavelength=np.array(temp_wavelength)
            temp_flux=np.array(temp_flux)
            temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)
            
            axes2[2, 3].plot(temp_wavelength, temp_norm_flux, c='red')
            axes2[2, 3].plot(wavelength, norm_flux, c='k')
            axes2[2, 3].set_ylim([0, max(norm_flux) + 0.25])
            axes2[2, 3].set_xlim([min(wavelength), max(wavelength)])
            axes2[2, 3].axis('off')
            
        return type_number
    fig1.canvas.draw()
    
    fig1.canvas.mpl_disconnect(fig1.canvas.manager.key_press_handler_id)
    fig1.canvas.mpl_connect('key_press_event', on_key_press)

def check_type(file_name, type_number) :
    '''
    Ellianna Schwab, Kelle Cruz

    check_type provides plots to qualitatively spectral type L-type brown dwarfs in NIR regime.
    Typing methods are taken from Cruz et al. 2017 and show comparison to Kirkpatrick 2010 typing templates.

    To use input path to file_name in a string followed by candidate spectral type as an integer. 
    Cruz et al. 2017 templates are shown band-by-band followed by Kirkpatrick 2010
    templates of the overall NIR spectrum. To save the plot, hit 's'.
    '''

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy.io import ascii

    hdulist = fits.open(file_name)
    spectrum = hdulist[0]
    wavelength = spectrum.data[0]
    flux = spectrum.data[1]


    #Create J, H and K bands for spectrum

    wavelength_J = []
    flux_J = []
    wavelength_H = []
    flux_H = []
    wavelength_K = []
    flux_K = []
    for jj in range(len(wavelength)):
        if wavelength[jj] >= 0.87 and wavelength[jj] <= 1.39:
            wavelength_J.append(wavelength[jj])
            flux_J.append(flux[jj])
        elif wavelength[jj] >= 1.41 and wavelength[jj] <= 1.89:
            wavelength_H.append(wavelength[jj])
            flux_H.append(flux[jj])
        elif wavelength[jj] >= 1.91 and wavelength[jj] <=2.39:
            wavelength_K.append(wavelength[jj])
            flux_K.append(flux[jj])


    #Normalize Each Band        

    wavelength_J = np.array(wavelength_J)
    flux_J = np.array(flux_J)
    wavelength_H = np.array(wavelength_H)
    flux_H = np.array(flux_H)
    wavelength_K = np.array(wavelength_K)
    flux_K = np.array(flux_K)

    flux_J = flux_J/np.mean(flux_J)
    flux_H = flux_H/np.mean(flux_H)
    flux_K = flux_K/np.mean(flux_K)

    
    #Normalize the Overall Spectra to Kirkpatrick Normalization
    norm_flux_array=[]
    for jj in range(len(wavelength)):
        if wavelength[jj] >= 1.28 and wavelength[jj] <= 1.39:
            norm_flux_array.append(flux[jj])
    norm_flux_array=np.array(norm_flux_array)

    wavelength=np.array(wavelength)
    flux=np.array(flux)
    norm_flux=flux/np.mean(norm_flux_array)


    #List the standards, so they can be called # "spectra/nir/prism_0835+1953_20050123_CHI06A.fits", this  one doesn't exist?
    NIR_standards = ["spectra/nir/U20165_0345+2540_new_w.fits", "spectra/nir/spex_prism_2130-0845_080713.fits", \
                     "spectra/nir/U10244_0408-1450.fits", "spectra/nir/u11291_1506+1321_050323.fits", \
                     "spectra/nir/U12101_2158-1550_davy.fits", "spectra/nir/spex_prism_2137+0808_U20909.fits", \
                     "spectra/nir/U10880.fits", "spectra/nir/u10721_050323.fits", "spectra/nir/2M1632.fits"]
    
    # If the number is in an easily bracketed part, or 1,2,3,4,5,6,7
    if type_number in range(1,8):

        #Create the Plots
        fig2, axes2 = plt.subplots(
            nrows=3, ncols=4, sharex=False, sharey=False, 
            gridspec_kw={'width_ratios':[1,1,1,3]}
            )

        for jj, ii in zip([0,1,2], [type_number-1, type_number, type_number+1]):
                ##J Band
                axes2[jj, 0].plot(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                  ascii.read("templates/L{}J_f.txt".format(ii))['col2'], c='red')
                axes2[jj, 0].fill_between(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                          ascii.read("templates/L{}J_f.txt".format(ii))['col4'], \
                                          ascii.read("templates/L{}J_f.txt".format(ii))['col5'], color='#c6c6c6')
                axes2[jj, 0].plot(wavelength_J, flux_J, c='k')
                axes2[jj, 0].annotate('L{}'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                axes2[jj, 0].set_ylim([-0.5, 2])
                axes2[jj, 0].axis('off')

                #H Band
                axes2[jj, 1].plot(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                  ascii.read("templates/L{}H_f.txt".format(ii))['col2'], c='red') 
                axes2[jj, 1].fill_between(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                          ascii.read("templates/L{}H_f.txt".format(ii))['col4'], \
                                          ascii.read("templates/L{}H_f.txt".format(ii))['col5'], color='#c6c6c6') 
                axes2[jj, 1].plot(wavelength_H, flux_H, c='k')
                axes2[jj, 1].axis('off')


                #K Band
                axes2[jj, 2].plot(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                  ascii.read("templates/L{}K_f.txt".format(ii))['col2'], c='red')
                axes2[jj, 2].fill_between(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                          ascii.read("templates/L{}K_f.txt".format(ii))['col4'], \
                                          ascii.read("templates/L{}K_f.txt".format(ii))['col5'], color='#c6c6c6')
                axes2[jj, 2].plot(wavelength_K, flux_K, c='k')
                axes2[jj, 2].axis('off')


                #All Together Now! This is where the Kirkpatrick 10 one comes in.
                temp_hdulist = fits.open(NIR_standards[ii])
                temp_spectrum = temp_hdulist[0]
                temp_wavelength = temp_spectrum.data[0]
                temp_flux = temp_spectrum.data[1]

                temp_norm_flux_array=[]
                for kk in range(len(temp_wavelength)):
                    if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                        temp_norm_flux_array.append(temp_flux[kk])
                temp_norm_flux_array=np.array(temp_norm_flux_array)

                temp_wavelength=np.array(temp_wavelength)
                temp_flux=np.array(temp_flux)
                temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)


                axes2[jj, 3].plot(temp_wavelength, temp_norm_flux, c='red')
                axes2[jj, 3].plot(wavelength, norm_flux, c='k')
                axes2[jj, 3].set_ylim([0, max(norm_flux) + 0.25])
                axes2[jj, 3].set_xlim([min(wavelength), max(wavelength)])
                axes2[jj, 3].axis('off')

    elif type_number == 0:
        #Create the Plots

        fig2, axes2 = plt.subplots(
            nrows=3, ncols=4, sharex=False, sharey=False, 
            gridspec_kw={'width_ratios':[1,1,1,3]}
            )
        
        #Begin with the M9 template and then iterate over the rest
        axes2[0, 0].annotate('M9', xy=(0.1, 0.8), xycoords='axes fraction', color='k')
        axes2[0, 0].axis('off')
        axes2[0, 1].axis('off')
        axes2[0, 2].axis('off')
        
        temp_hdulist = fits.open("spectra/M9V_LHS2924.fits")
        temp_spectrum = temp_hdulist[0]
        temp_wavelength = temp_spectrum.data[0]
        temp_flux = temp_spectrum.data[1]
        temp_norm_flux_array=[]
        for kk in range(len(temp_wavelength)):
            if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                temp_norm_flux_array.append(temp_flux[kk])
        temp_norm_flux_array=np.array(temp_norm_flux_array)

        temp_wavelength=np.array(temp_wavelength)
        temp_flux=np.array(temp_flux)
        temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)
        axes2[0, 3].plot(temp_wavelength, temp_norm_flux, c='red')
        axes2[0, 3].plot(wavelength, norm_flux, c='k')
        axes2[0, 3].set_ylim([0, max(norm_flux) + 0.25])
        axes2[0, 3].set_xlim([min(wavelength), max(wavelength)])
        axes2[0, 3].axis('off')

        for jj, ii in zip([1,2], [type_number, type_number+1]):
                ##J Band
                axes2[jj, 0].plot(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                  ascii.read("templates/L{}J_f.txt".format(ii))['col2'], c='red')
                axes2[jj, 0].fill_between(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                          ascii.read("templates/L{}J_f.txt".format(ii))['col4'], \
                                          ascii.read("templates/L{}J_f.txt".format(ii))['col5'], color='#c6c6c6')
                axes2[jj, 0].plot(wavelength_J, flux_J, c='k')
                axes2[jj, 0].annotate('L{}'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                axes2[jj, 0].set_ylim([-0.5, 2])
                axes2[jj, 0].axis('off')

                #H Band
                axes2[jj, 1].plot(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                  ascii.read("templates/L{}H_f.txt".format(ii))['col2'], c='red') 
                axes2[jj, 1].fill_between(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                          ascii.read("templates/L{}H_f.txt".format(ii))['col4'], \
                                          ascii.read("templates/L{}H_f.txt".format(ii))['col5'], color='#c6c6c6') 
                axes2[jj, 1].plot(wavelength_H, flux_H, c='k')
                axes2[jj, 1].axis('off')


                #K Band
                axes2[jj, 2].plot(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                  ascii.read("templates/L{}K_f.txt".format(ii))['col2'], c='red')
                axes2[jj, 2].fill_between(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                          ascii.read("templates/L{}K_f.txt".format(ii))['col4'], \
                                          ascii.read("templates/L{}K_f.txt".format(ii))['col5'], color='#c6c6c6')
                axes2[jj, 2].plot(wavelength_K, flux_K, c='k')
                axes2[jj, 2].axis('off')


                #All Together Now! This is where the Kirkpatrick 10 one comes in.
                temp_hdulist = fits.open(NIR_standards[ii])
                temp_spectrum = temp_hdulist[0]
                temp_wavelength = temp_spectrum.data[0]
                temp_flux = temp_spectrum.data[1]

                temp_norm_flux_array=[]
                for kk in range(len(temp_wavelength)):
                    if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                        temp_norm_flux_array.append(temp_flux[kk])
                temp_norm_flux_array=np.array(temp_norm_flux_array)

                temp_wavelength=np.array(temp_wavelength)
                temp_flux=np.array(temp_flux)
                temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)


                axes2[jj, 3].plot(temp_wavelength, temp_norm_flux, c='red')
                axes2[jj, 3].plot(wavelength, norm_flux, c='k')
                axes2[jj, 3].set_ylim([0, max(norm_flux) + 0.25])
                axes2[jj, 3].set_xlim([min(wavelength), max(wavelength)])
                axes2[jj, 3].axis('off')
                
    elif type_number == 8:
        #List the standards, so they can be called # "spectra/nir/prism_0835+1953_20050123_CHI06A.fits", this  one doesn't exist?
        NIR_standards = ["spectra/nir/u10721_050323.fits", "spectra/nir/2M1632.fits"]

        #Create the Plots

        fig2, axes2 = plt.subplots(
            nrows=3, ncols=4, sharex=False, sharey=False, 
            gridspec_kw={'width_ratios':[1,1,1,3]}
            )

        for jj, ii in zip([0,1], [type_number-1, type_number]):
                ##J Band
                axes2[jj, 0].plot(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                  ascii.read("templates/L{}J_f.txt".format(ii))['col2'], c='red')
                axes2[jj, 0].fill_between(ascii.read("templates/L{}J_f.txt".format(ii))['col1'], \
                                          ascii.read("templates/L{}J_f.txt".format(ii))['col4'], \
                                          ascii.read("templates/L{}J_f.txt".format(ii))['col5'], color='#c6c6c6')
                axes2[jj, 0].plot(wavelength_J, flux_J, c='k')
                axes2[jj, 0].annotate('L{}'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                axes2[jj, 0].set_ylim([-0.5, 2])
                axes2[jj, 0].axis('off')

                #H Band
                axes2[jj, 1].plot(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                  ascii.read("templates/L{}H_f.txt".format(ii))['col2'], c='red') 
                axes2[jj, 1].fill_between(ascii.read("templates/L{}H_f.txt".format(ii))['col1'], \
                                          ascii.read("templates/L{}H_f.txt".format(ii))['col4'], \
                                          ascii.read("templates/L{}H_f.txt".format(ii))['col5'], color='#c6c6c6') 
                axes2[jj, 1].plot(wavelength_H, flux_H, c='k')
                axes2[jj, 1].axis('off')


                #K Band
                axes2[jj, 2].plot(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                  ascii.read("templates/L{}K_f.txt".format(ii))['col2'], c='red')
                axes2[jj, 2].fill_between(ascii.read("templates/L{}K_f.txt".format(ii))['col1'], \
                                          ascii.read("templates/L{}K_f.txt".format(ii))['col4'], \
                                          ascii.read("templates/L{}K_f.txt".format(ii))['col5'], color='#c6c6c6')
                axes2[jj, 2].plot(wavelength_K, flux_K, c='k')
                axes2[jj, 2].axis('off')


                #All Together Now! This is where the Kirkpatrick 10 one comes in.
                temp_hdulist = fits.open(NIR_standards[jj])
                temp_spectrum = temp_hdulist[0]
                temp_wavelength = temp_spectrum.data[0]
                temp_flux = temp_spectrum.data[1]

                temp_norm_flux_array=[]
                for kk in range(len(temp_wavelength)):
                    if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                        temp_norm_flux_array.append(temp_flux[kk])
                temp_norm_flux_array=np.array(temp_norm_flux_array)

                temp_wavelength=np.array(temp_wavelength)
                temp_flux=np.array(temp_flux)
                temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)


                axes2[jj, 3].plot(temp_wavelength, temp_norm_flux, c='red')
                axes2[jj, 3].plot(wavelength, norm_flux, c='k')
                axes2[jj, 3].set_ylim([0, max(norm_flux) + 0.25])
                axes2[jj, 3].set_xlim([min(wavelength), max(wavelength)])
                axes2[jj, 3].axis('off')
                
        #Finish with T0 template
        axes2[2, 0].annotate('T0', xy=(0.1, 0.8), xycoords='axes fraction', color='k')
        axes2[2, 0].axis('off')
        axes2[2, 1].axis('off')
        axes2[2, 2].axis('off')
        
        temp_hdulist = fits.open("spectra/T0_2M0423_T0.fits")
        temp_spectrum = temp_hdulist[0]
        temp_wavelength = temp_spectrum.data[0]
        temp_flux = temp_spectrum.data[1]
        temp_norm_flux_array=[]
        for kk in range(len(temp_wavelength)):
            if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                temp_norm_flux_array.append(temp_flux[kk])
        temp_norm_flux_array=np.array(temp_norm_flux_array) 
        temp_wavelength=np.array(temp_wavelength)
        temp_flux=np.array(temp_flux)
        temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)
        
        axes2[2, 3].plot(temp_wavelength, temp_norm_flux, c='red')
        axes2[2, 3].plot(wavelength, norm_flux, c='k')
        axes2[2, 3].set_ylim([0, max(norm_flux) + 0.25])
        axes2[2, 3].set_xlim([min(wavelength), max(wavelength)])
        axes2[2, 3].axis('off')