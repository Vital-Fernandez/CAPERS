import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import lime

file_address = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/source/jwst_nirspec_prism_disp.fits'
print(fits.info(file_address))


redshift = 4.299
spec_address = '/home/vital/PycharmProjects/ceers-data/data/spectra/CEERs_DR0.9/nirspecDDT/prism/hlsp_ceers_jwst_nirspec_nirspecDDT-001586_prism_dr0.9_x1d.fits'
spec = lime.Spectrum.from_file(spec_address, instrument='nirspec', redshift=redshift)
# spec.plot.spectrum()

wave_arr = spec.wave.data

R_arr = np.empty(wave_arr.size)
deltalamb_arr = np.diff(wave_arr)
R_arr[:-1] = wave_arr[:-1] / deltalamb_arr
R_arr[-1] = R_arr[-2]


with fits.open(file_address) as hdul:
    data = hdul[1].data

'''
In the dispersion curve files from the NIRSpec instrument, the column labeled 'dlds' represents the dispersion value at
 each wavelength, expressed in micrometers per pixel (μm/pixel). This indicates how much the wavelength changes per pixel
  on the detector at a given point in the spectrum. These files typically contain three columns: wavelength (μm), dispersion
   (μm/pixel), and resolution (λ/Δλ, unitless).
'''

fig, ax = plt.subplots()
ax2 = ax.twinx()
disper_line = ax.plot(data['WAVELENGTH'], data['DLDS'], label='Dispersion', color='tab:green')
R_line = ax2.plot(data['WAVELENGTH'], data['R'], label='R')
R_line += ax2.plot(wave_arr, R_arr, label='LiMe R')
ax.update({'xlabel':r'wavelength $(\mu m)$', 'ylabel':r'Dispersion $(\mu m / pixel)$', 'title':'NIRSPEC Prism dispersion curve'})
ax2.update({'ylabel':r'$R = \frac{\lambda}{\Delta \lambda}$'})
lines = disper_line + R_line
labels = [item.get_label() for item in lines]
ax.legend(lines, labels, loc="upper left")
plt.tight_layout()
plt.show()

print(data)