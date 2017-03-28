# to run this, at terminal, >import plot_optical

from astrodbkit import astrodb
import matplotlib.pyplot as plt
import numpy as np

# setup matplotlib to use latex
from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

db = astrodb.Database('/Users/kelle/Dropbox/BDNYCdb/BDNYCdev.db')

# path to where spectra should be written to
out_path = '/Users/kelle/Dropbox/Analysis/NIRtemplates/NIRSpecFigures/plots/opt/'

# Explictly decide which spectra to plot
spectra_ids = [1093, 1103, 470, 1125, 496, 578, 586, 1301, 400, 401, 1353, 403,
               874, 404, 405, 420, 408, 1521, 1503]

# enter the strings by hand
spectral_type_strings = ["L0$\gamma$", "L2$\gamma$", "L1$\gamma$", "L0$\gamma$", "L0$\gamma$", "L1$\gamma$", "L2$\gamma$",
                         "L0", "L4$\gamma$", "L3.5$\gamma$", "L3.5$\gamma$", "L4$\gamma$", "L4$\gamma$", "L0$\gamma$",
                         "L3.5$\gamma$", "L1", "L5", "L4$\gamma$", "L0$\gamma$"]

ids_and_types = zip(spectra_ids, spectral_type_strings)

# look up names, obsdates, and spectral types for the spectra
names2 = []
obsdates = []
for spectra_id, spectral_type_string in ids_and_types:
    source_id2 = db.query("SELECT source_id from spectra where id={} ".format(spectra_id))[0][0]
    name2 = db.query("SELECT shortname from sources where id={} ".format(source_id2))[0][0]
    obsdate = db.query("SELECT obs_date from spectra where id={} ".format(spectra_id))[0][0]
    spectral_type = db.query(
        "SELECT spectral_type FROM spectral_types WHERE source_id={} AND adopted = 1".format(source_id2))
    sptype_gravity = db.query(
        "SELECT gravity FROM spectral_types WHERE source_id={} AND adopted = 1".format(source_id2))
    sptype_suffix = db.query(
        "SELECT suffix FROM spectral_types WHERE source_id={} AND adopted = 1".format(source_id2))
    names2.append(name2)
    obsdates.append(obsdate)
    print name2, source_id2, spectral_type, sptype_gravity, sptype_suffix, spectral_type_string, spectra_id, obsdate

for j, spectra_id in enumerate(spectra_ids):
        spectrum = db.query("SELECT spectrum from spectra where id='{}'".format(spectra_id), fetch='one')[0]
        wave = spectrum.data[0]
        flux = spectrum.data[1]

        plt.plot(wave, flux / np.mean(flux))
        # plt.show() #go to screen

        plt.xlim(6000, 10000)
        plt.ylim(0.0, 4.0)
        plt.xlabel(r'Wavelength ($\mbox{\AA}$)')
        plt.ylabel('Normalized Flux')
        plt.text(6500, 3.5, '{}'.format(names2[j]))
        plt.text(6500, 3.3, r'{}'.format(spectral_type_strings[j]))
        plt.text(6500, 3.0, 'obs date: ' '{}'.format(obsdates[j]), size='smaller')

        plt.savefig(out_path + '{}'.format(names2[j]) + '_{}'.format(spectra_id) + '.pdf')

        plt.close()
