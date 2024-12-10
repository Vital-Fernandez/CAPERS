import pandas as pd

import lime
import numpy as np
from support.tools import search_spectra_ceers, create_backup
from pathlib import Path

# Load configuration
cfg_file = '../CAPERS.toml'
capers_cfg = lime.load_cfg('../CAPERS.toml')

# Get sample information
version = capers_cfg['meta']['version']
sample = capers_cfg['meta']['sample']

observations_folder =  Path(capers_cfg['file_structure']['observations_folder'])
data_folder = Path(capers_cfg['file_structure']['data_folder'])
tables_folder = data_folder/'tables'
source_folder = data_folder/'source'

fits_ext_list = capers_cfg['file_structure']['fits_ext_list']
file_log_headers = capers_cfg['file_structure']['file_log_headers']
# ceers_redshift_name = capers_cfg['file_structure']['CEERs_redshift_file']
ceers_redshift_name = capers_cfg['file_structure']['CEERs_redshift_file']

# Create day back-ups
create_backup(tables_folder/f'{sample}_{version}_files_log.txt')

# Search the file names
raw_spectra_list = search_spectra_ceers(observations_folder/'CEERs_DR0.9', ext_list=fits_ext_list)

# Adjust the selection
# spectra_list = review_masked_files(raw_spectra_list)
spectra_list = raw_spectra_list

# Empty frame for the sample data
file_log = pd.DataFrame(columns=file_log_headers)
file_log.set_index(['sample', 'id', 'file'], inplace=True)

# Loop through the files and add the data
print(f'- Filling files log')
for idx_spec, fits_file in enumerate(spectra_list):

    # print(idx_spec, fits_file)
    dir_parts = fits_file.parts[1:-1]
    file_dir = Path(*dir_parts)
    file_name = fits_file.name

    file_path = Path('/'.join(fits_file.parts[7:-1]))/file_name
    pointing_dir = dir_parts[7]

    name_comps = file_name.split('_')
    disp = name_comps[5] #dir_parts[2]
    pointing_name, id_label = name_comps[4].split('-')

    # MSA = int(id_label)
    MSA_number = int(id_label)
    id_label = str(int(id_label)) if pointing_name != 'nirspecDDT' else 'D' + str(int(id_label))

    ext = name_comps[-1][0:-5]

    if pointing_name != pointing_dir:
        print(f'Missmatch folder: {pointing_dir} != {file_path}')

    file_log.loc[(sample, id_label, file_name), "MSA_number"] = MSA_number
    file_log.loc[(sample, id_label, file_name), "disp"] = disp
    file_log.loc[(sample, id_label, file_name), "pointing"] = pointing_name
    file_log.loc[(sample, id_label, file_name), "ext"] = ext
    file_log.loc[(sample, id_label, file_name), "file_path"] = file_path.as_posix()

# Adjust and organize the dataframe
idcs_ddt = file_log.pointing == 'nirspecDDT'
ceers_df = file_log.loc[~idcs_ddt].copy()
ddt_df = file_log.loc[idcs_ddt].copy()

ceers_df.sort_values(by=['MSA_number', 'disp', 'ext'], ascending=[True, True, True], inplace=True)
ddt_df.sort_values(by=['MSA_number', 'disp', 'ext'], ascending=[True, True, True], inplace=True)
file_log = pd.concat([ceers_df, ddt_df], axis=0)

# Add previous measurements
redshift_df = pd.read_csv('/home/vital/Dropbox/Astrophysics/Data/CEERs/RedShift_compilation/CEERS_NIRSpec_MSA_catalog_dr0.9_extended.csv', delimiter=';', index_col=0)
MSA_list = redshift_df.index.unique()
for MSA in MSA_list:
    idcs_obj = file_log.index.get_level_values('id') == MSA
    file_log.loc[idcs_obj, 'z_phot'] = redshift_df.loc[MSA, 'z_best']
    file_log.loc[idcs_obj, 'ra'] = redshift_df.loc[MSA, 'ra']
    file_log.loc[idcs_obj, 'dec'] = redshift_df.loc[MSA, 'dec']
lime.save_frame(tables_folder/f'{sample}_{version}_files_log.txt', file_log)


# redshift_df = pd.read_csv(source_folder/ceers_redshift_name, delimiter=';', index_col=0)
# MSA_list = redshift_df.index.unique()
# for MSA in MSA_list:
#
#     idcs_obj = file_log.index.get_level_values('id') == MSA
#     file_log.loc[idcs_obj, 'z_phot'] = redshift_df.loc[MSA, 'old_z_phot']
#     file_log.loc[idcs_obj, 'ra'] = redshift_df.loc[MSA, 'ra']
#     file_log.loc[idcs_obj, 'dec'] = redshift_df.loc[MSA, 'dec']
#
# # Save the log
# lime.save_frame(tables_folder/f'{sample}_{version}_files_log.txt', file_log)