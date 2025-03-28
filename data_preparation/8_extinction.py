import lime
import numpy as np
from pathlib import Path
from astropy.io import fits
from support.tools import nirspec_load_function
import specsy as sy
import pyneb as pn

# State data location
sample = 'CAPERS_UDS_V0.1'
tables_folder ='/home/vital/Dropbox/Astrophysics/Data/CAPERS/tables'
observations_folder = f'/home/vital/Astrodata/CAPERS'
flux_log_address = f'{tables_folder}/{sample}_flux_log.txt'

# Create observations sample variable
capers_sample = lime.Sample(flux_log_address, levels=["sample", "id", "file", "line"], norm_flux=1e-22,
                            load_function=nirspec_load_function, folder_obs=observations_folder)

# Index target object
idcs_object = capers_sample.frame.MPT_number == 22431
first_idx = capers_sample.frame.loc[idcs_object].index[0]
spec = capers_sample.get_spectrum(first_idx, z_column='z_bands')
lines_df = capers_sample.frame.loc[idcs_object].droplevel(['sample', 'id', 'file'])
spec.load_frame(lines_df)
spec.plot.spectrum(rest_frame=True)
# sy.extinction.cHbeta_from_log(lines_df, ref_line='H1_4862A', flux_entry='profile', show_plot=True)

cHbeta, cHbeta_err = sy.extinction.cHbeta_from_log(lines_df, ref_line='H1_4862A', flux_entry='profile',
                                            lines_ignore=['H1_3889A', 'H1_3971A'], show_plot=True)

lime.normalize_fluxes(lines_df, norm_list='H1_4862A')
flux_dict = sy.flux_distribution(lines_df, flux_type='line_flux')
# print(flux_dict)
print(lines_df.z_line)
pn.atomicData.setDataFile('o_iii_coll_AK99.dat')
O3 = pn.Atom('O', 3)
print(O3.lineList)
temp_dist = O3.getTemDen(flux_dict['O3_4364A']/(flux_dict['O3_4960A']+flux_dict['O3_5008A']), den=100, to_eval='L(4363)/(L(4959)+L(5007))')
TOIII, TOIII_err = np.nanmean(temp_dist), np.nanstd(temp_dist)
print(f'TOIII = {TOIII:0.1f} +/- {TOIII_err:0.1f}')

# temp_dist = O3.getTemDen(flux_dict['O3_1664A']/(flux_dict['O3_5008A']), den=100, to_eval='(L(1658)+L(1661)+L(1666))/(L(5007))')
# TOIII, TOIII_err = np.nanmean(temp_dist), np.nanstd(temp_dist)
# print(f'TOIII (1664A) = {TOIII:0.1f} +/- {TOIII_err:0.1f}')

O3_abund_dist = O3.getIonAbundance(flux_dict['O3_5008A'], tem=temp_dist, den=np.ones(flux_dict['O3_5008A'].size)*100, wave=5007, Hbeta=1)
O3_abund, O3_abund_err = np.nanmean(12 + np.log10(O3_abund_dist)),  np.nanstd(12 + np.log10(O3_abund_dist))
print(O3_abund, O3_abund_err)
# # Plot the data
# spec.plot.spectrum(rest_frame=True)
# spec.plot.grid()

