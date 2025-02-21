from pathlib import Path
import numpy as np
import pandas as pd
import lime
from support.tools import create_backup, nirspec_load_function, save_detection_results, read_prism_r_curve, catch_nan_spectra

# Load configuration
cfg_file = '../CAPERS.toml'
capers_cfg = lime.load_cfg('../CAPERS.toml')

# Get sample information
version = capers_cfg['meta']['version']
sample = capers_cfg['meta']['sample']
preference_list = capers_cfg['file_structure']['ext_1d_fits']

observations_folder =  Path(capers_cfg['file_structure']['observations_folder'])
data_folder = Path(capers_cfg['file_structure']['data_folder'])
tables_folder = data_folder/'tables'
comps_folder = data_folder/'components'
bands_folder = data_folder/'bands'
logs_folder = data_folder/'line_logs'
prism_bands_df = lime.load_frame(data_folder/capers_cfg['file_structure']['prism_bands'])
redshift_bands_df = lime.load_frame(data_folder/capers_cfg['file_structure']['redshift_bands'])
REF_LINEs = redshift_bands_df.index.to_numpy()
R_interpolator = read_prism_r_curve(data_folder/capers_cfg['file_structure']['prism_R_curve_file'])

# Line fitting data
norm_flux = capers_cfg['data']['norm_flux']
default_cfg_section = capers_cfg['data']['default_prism_cfg_section']
boundary_cfg_section = capers_cfg['data']['boundary_prism_cfg_section']

# Load the sample table (create daily back-up if necessary)
log_address = tables_folder/f'{sample}_files_log.txt'
create_backup(log_address)
files_sample = lime.Sample(log_address, levels=["sample", "id", "file"], load_function=nirspec_load_function,
                           norm_flux=norm_flux, folder_obs=observations_folder)

# List of files
log_list = list(logs_folder.glob(f'*_log.txt'))

# Containers for the data
id_list = [None] * len(log_list)
file_list = [None] * len(log_list)
# fits_path_list = [None] * len(log_list)

for i, log_address in enumerate(log_list):

    file_name = log_address.name
    fits_file = file_name[:file_name.rfind('_log')] + '.fits'
    name_comps = file_name.split('_')

    id_label = name_comps[3]
    MPT_number = int(id_label[1:])
    ext = name_comps[-2]

    id_list[i] = id_label
    file_list[i] = fits_file

# Save first version
flux_address = tables_folder/f'{sample}_flux_log.txt'
flux_sample = lime.Sample.from_file(id_list=id_list, log_list=log_list, file_list=file_list, load_function=nirspec_load_function)

# # Re-index with sample and sort
flux_sample.frame['sample'] = sample
flux_sample.frame['id'] = flux_sample.index.get_level_values('id')
flux_sample.frame['file'] = flux_sample.index.get_level_values('file')
flux_sample.frame['line'] = flux_sample.index.get_level_values('line')
flux_sample.frame.set_index(['sample', 'id', 'file', 'line'], inplace=True)
flux_sample.frame['MPT_number'] = np.nan
flux_sample.frame['z_bands'] = np.nan
flux_sample.frame['file_path'] = ''

# Add the file_path
idcs_entry = flux_sample.frame.index.droplevel('line').unique()
for idx in idcs_entry:
    flux_sample.frame.loc[idx, 'file_path'] = files_sample.frame.loc[idx, 'file_path']
    flux_sample.frame.loc[idx, 'z_bands'] = files_sample.frame.loc[idx, 'z_bands']
    flux_sample.frame.loc[idx, 'MPT_number'] = files_sample.frame.loc[idx, 'MPT_number']
flux_sample.save_frame(flux_address)
