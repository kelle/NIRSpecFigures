# to run this, at terminal, >import plot_optical

from astrodbkit import astrodb
import matplotlib.pyplot as plt
import numpy as np

db = astrodb.Database('/Users/kelle/Dropbox/BDNYCdb/BDNYCdev.db')

# path to where spectra should be written to
out_path = '/Users/kelle/Dropbox/Analysis/NIRtemplates/NIRSpecFigures/plots/opt/'

# Explictly decide which spectra to plot
spectra_ids = [1093, 1103, 470, 1125, 496, 578, 586, 1301, 400, 401, 1353, 403,
               874, 404, 405, 420, 408, 1521, 1503]

# print spectra_ids
names = []

# look up names, obsdates, and spectral types for the spectra
names2 = []
obsdates = []
spectral_types = []
for spectra_id in spectra_ids:
    source_id2 = db.query("SELECT source_id from spectra where id={} ".format(spectra_id))[0][0]
    name2 = db.query("SELECT shortname from sources where id={} ".format(source_id2))[0][0]
    obsdate = db.query("SELECT obs_date from spectra where id={} ".format(spectra_id))[0][0]
    spectral_type = db.query("SELECT spectral_type FROM spectral_types WHERE source_id={} AND adopted = 1".format(source_id2))
    names2.append(name2)
    obsdates.append(obsdate)
    print name2, source_id2, spectral_type, spectra_id, obsdate

for j, spectra_id in enumerate(spectra_ids):
        spectrum = db.query("SELECT spectrum from spectra where id='{}'".format(spectra_id), fetch='one')[0]
        wave = spectrum.data[0]
        flux = spectrum.data[1]

        plt.plot(wave, flux / np.mean(flux))
        # plt.show() #go to screen

        plt.xlim(6000, 10000)
        plt.ylim(0.0, 4.0)
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.ylabel('Normalized Flux')
        plt.text(6500, 3.5, '{}'.format(names2[j]))
        plt.text(6500, 3.2, 'obs date: ' '{}'.format(obsdates[j]), size='smaller')

        plt.savefig(out_path + '{}'.format(names2[j]) + '_{}'.format(spectra_id) + '.pdf')

        plt.close()
