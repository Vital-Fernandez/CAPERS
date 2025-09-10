import lime
from support.tools import  nirspec_load_function
from pathlib import Path

from support.tools import create_backup, nirspec_load_function, save_detection_results, read_prism_r_curve, \
    catch_nan_spectra, recover_bands

lime.theme.set_style('dark')

# Load configuration
cfg_file = '../CAPERS_v2.toml'
capers_cfg = lime.load_cfg(cfg_file)

# Get sample information
version = capers_cfg['meta']['version']
sample = capers_cfg['meta']['sample']
preference_list = capers_cfg['file_structure']['ext_1d_fits']

observations_folder = Path(capers_cfg['file_structure']['observations_folder'])
data_folder = Path(capers_cfg['file_structure']['data_folder'])

tables_folder = data_folder/sample/'tables'
source_folder = data_folder/sample/'source'
bands_folder = data_folder/sample/'bands'
logs_folder = data_folder/sample/'line_logs'
comps_folder = data_folder/sample/'comps'

prism_bands_df = lime.load_frame(data_folder/capers_cfg['file_structure']['prism_bands'])
measure_bands_df = prism_bands_df
redshift_lines = capers_cfg['file_structure']['lines_redshift']

lime.theme.set_style('dark')
folder = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/22431_reductions'

fname = '/home/vital/Downloads/54833_alt_nod_reduction/CAPERS-EGS_combpoint_s000054833_x1d_optext.fits'

redshift = 0
norm_flux = None
log_list = []

spec = lime.Spectrum.from_file(fname, instrument='nirspec', redshift=redshift, norm_flux=norm_flux)
spec.unit_conversion(wave_units_out='AA', flux_units_out='FLAM')
# spec.plot.spectrum(rest_frame=True, include_err=True)
spec.check.redshift(reference_lines=redshift_lines)


# log_file = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/line_logs/CAPERS_UDS_P1_s000022431_x1d_optext_log.txt'
# spec.load_frame(log_file)


# file_addresses = ['CAPERS_UDS_P1_c1_s000022431_x1d_optext.fits',
#                   'CAPERS_UDS_P1_c2_s000022431_x1d_optext.fits',
#                   'CAPERS_UDS_P1_c3_s000022431_x1d_optext.fits']

#
# for i, fname in enumerate(file_addresses):
#     spec = lime.Spectrum.from_file(f'{folder}/{fname}', instrument='nirspec', redshift=redshift, norm_flux=norm_flux)
#     spec.unit_conversion(wave_units_out='AA', flux_units_out='FLAM')
#     # spec.plot.spectrum(rest_frame=True)
#     df_i = lime.load_frame(log_file)
#     df_i['file_path'] = file_addresses[i]
#     log_list.append(Path(f'{folder}/{Path(fname).stem}.txt'))
#     lime.save_frame(log_list[i], df_i)
#
# id_list = ['Visit 1', 'Visit 2', 'Visit 3']
# folder_obs = folder
# obs_list = file_addresses
# sample1 = lime.Sample.from_file(id_list, log_list, obs_list, folder_obs=folder_obs, load_function=nirspec_load_function,
#                                 norm_flux=norm_flux, redshift=redshift)
#
# sample1.plot.spectra(rest_frame=True, legend_handle='id', include_err=True)



# import lime
# from support.tools import  nirspec_load_function
# from pathlib import Path
#
# lime.theme.set_style('dark')
# folder = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/22431_reductions'
#
# fname = '/home/vital/Astrodata/CAPERS/CAPERS_UDS_V0.1/P1/CAPERS_UDS_P1_s000022431_x1d_optext.fits'
#
# file_addresses = ['CAPERS_UDS_P1_c1_s000022431_x1d_optext.fits',
#                   'CAPERS_UDS_P1_c2_s000022431_x1d_optext.fits',
#                   'CAPERS_UDS_P1_c3_s000022431_x1d_optext.fits']
#
# log_file = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/line_logs/CAPERS_UDS_P1_s000022431_x1d_optext_log.txt'
#
# redshift = 9.272523
# norm_flux = 1e-20
# log_list = []
#
# spec = lime.Spectrum.from_file(fname, instrument='nirspec', redshift=redshift, norm_flux=norm_flux)
# spec.unit_conversion(wave_units_out='AA', flux_units_out='FLAM')
# spec.load_frame(log_file)
# spec.plot.spectrum(rest_frame=True, include_err=True)
#
# for i, fname in enumerate(file_addresses):
#     spec = lime.Spectrum.from_file(f'{folder}/{fname}', instrument='nirspec', redshift=redshift, norm_flux=norm_flux)
#     spec.unit_conversion(wave_units_out='AA', flux_units_out='FLAM')
#     # spec.plot.spectrum(rest_frame=True)
#     df_i = lime.load_frame(log_file)
#     df_i['file_path'] = file_addresses[i]
#     log_list.append(Path(f'{folder}/{Path(fname).stem}.txt'))
#     lime.save_frame(log_list[i], df_i)
#
#
#
# id_list = ['Visit 1', 'Visit 2', 'Visit 3']
# folder_obs = folder
# obs_list = file_addresses
# sample1 = lime.Sample.from_file(id_list, log_list, obs_list, folder_obs=folder_obs, load_function=nirspec_load_function,
#                                 norm_flux=norm_flux, redshift=redshift)
#
# sample1.plot.spectra(rest_frame=True, legend_handle='id', include_err=True)