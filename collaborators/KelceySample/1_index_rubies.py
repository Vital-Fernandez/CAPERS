import numpy as np
import pandas as pd
import lime

import lime
import pandas as pd
from support.tools import search_spectra_capers, create_backup
from pathlib import Path

# Load configuration
cfg_file = 'rubies_v3.toml'
capers_cfg = lime.load_cfg(cfg_file)

# Get sample information
version = capers_cfg['meta']['version']
sample = capers_cfg['meta']['sample']

observations_folder = Path(capers_cfg['file_structure']['observations_folder'])
data_folder = Path(capers_cfg['file_structure']['data_folder'])
source_folder = Path(capers_cfg['file_structure']['source_folder'])
tables_folder = data_folder/sample/'tables'

fits_ext_list = capers_cfg['file_structure']['fits_ext_list']

headers_id = capers_cfg['file_structure']['headers_id']
headers_sample = capers_cfg['file_structure']['headers_sample']
headers_analysis = capers_cfg['file_structure']['headers_analysis']

# Create day back-ups
create_backup(tables_folder/f'{sample}_files_log.txt')

# Load the sample data table
# redshift_df = pd.read_csv(source_folder/capers_cfg['file_structure'][f'{sample}_redshift_file'], delimiter=';', index_col=0)

# Search the file names
# raw_spectra_list = search_spectra_capers(observations_folder/sample, ext_list=fits_ext_list)
# fname = '/home/vital/Astrodata/KelceySample/rubies_prism_files.txt'
fname = '/home/vital/Astrodata/KelceySample/rubies_grat_files.txt'
raw_spectra_list = np.loadtxt(fname, dtype=str)

# Adjust the selection
# spectra_list = review_masked_files(raw_spectra_list)
spectra_list = raw_spectra_list

# Empty frame for the sample data
file_log = pd.DataFrame(columns=capers_cfg['file_structure']['headers_id']+
                                capers_cfg['file_structure']['headers_sample']+
                                capers_cfg['file_structure']['headers_analysis']+
                                capers_cfg['file_structure']['headers_photometry'])

file_log.set_index(['sample', 'id', 'pointing'], inplace=True)

# Default values
print(f'- Filling files log')
for idx_spec, fits_file in enumerate(spectra_list):

    if sample == 'rubies_prism':
        fits_file = Path(fits_file)
        dir_parts = fits_file.parts[:-1]
        file_dir = Path(*dir_parts)
        file_name = fits_file.name

        file_path = Path('/'.join(dir_parts))

        name_comps = file_name.split('_')
        disp = name_comps[3]
        pointing_name = name_comps[0]

        id_label = name_comps[4]
        MPT_number = int(id_label[1:])
        ext = name_comps[-1][0:-5]

        file_log.loc[(sample, id_label, pointing_name), "MPT_number"] = MPT_number
        file_log.loc[(sample, id_label, pointing_name), "disp"] = disp
        file_log.loc[(sample, id_label, pointing_name), "z_tier"] = 0
        file_log.loc[(sample, id_label, pointing_name), ext] = file_name
        file_log.loc[(sample, id_label, pointing_name), "file_path"] = file_path.as_posix()

    if sample == 'rubies_grat':
        fits_file = Path(fits_file)
        dir_parts = fits_file.parts[:-1]
        file_dir = Path(*dir_parts)
        file_name = fits_file.name

        file_path = Path('/'.join(dir_parts))

        name_comps = file_name.split('_')
        disp = name_comps[3]
        pointing_name = name_comps[0]

        id_label = name_comps[4]
        MPT_number = int(id_label[1:])
        ext = name_comps[-1][0:-5]

        file_log.loc[(sample, id_label, pointing_name), "MPT_number"] = MPT_number
        file_log.loc[(sample, id_label, pointing_name), "disp"] = disp
        file_log.loc[(sample, id_label, pointing_name), "z_tier"] = 0
        file_log.loc[(sample, id_label, pointing_name), ext] = file_name
        file_log.loc[(sample, id_label, pointing_name), "file_path"] = file_path.as_posix()



# # Adjust and organize the dataframe
# file_log.sort_values(by=['MPT_number'], ascending=[True], inplace=True)
#
# # Add previous measurements
# MPT_list = redshift_df.index.unique()
# redshift_headers  = redshift_df.columns.to_numpy()
# for MPT in MPT_list:
#     idcs_obj = file_log.MPT_number == MPT
#     for hdr in redshift_headers:
#         if hdr in file_log.columns:
#             formatted_entry = str(redshift_df.loc[MPT, hdr]).replace(' ','')
#             if 'flux_' not in hdr:
#                 file_log.loc[idcs_obj, hdr] = formatted_entry
#             else:
#                 file_log.loc[idcs_obj, hdr.lower()] = formatted_entry
#
# # Reassign previous redshifts
# if sample == 'CAPERS_EGS_V0.2.1':
#     previous_version_fname = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/sample/CAPERS_EGS_V0.2/tables/CAPERS_EGS_V0.2_files_log_DataSource.txt'
#     redshift_manual_df = lime.load_frame(previous_version_fname, levels=["sample", "id", "file"])
#     idcs_manual = redshift_manual_df.loc[pd.notnull(redshift_manual_df.z_manual)].index
#     for idx in idcs_manual:
#         pointing =  redshift_manual_df.loc[idx, 'pointing']
#         idx_new = (sample, idx[1], pointing)
#         if idx_new not in file_log.index:
#             print(f'WARNING {idx_new} NOTTTT THERE!!!!!')
#         else:
#             file_log.loc[idx_new, 'z_manual'] = redshift_manual_df.loc[idx, 'z_manual']
#             file_log.loc[idx_new, 'z_tier'] = 3
#
# if sample == 'CAPERS_UDS_V0.1':
#     previous_version_fname = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/sample/CAPERS_UDS_V0.1/tables/CAPERS_UDS_V0.1_files_log_data_source.txt'
#     redshift_manual_df = lime.load_frame(previous_version_fname, levels=["sample", "id", "file"])
#     idcs_manual = redshift_manual_df.loc[pd.notnull(redshift_manual_df.z_manual)].index
#     for idx in idcs_manual:
#         pointing =  redshift_manual_df.loc[idx, 'pointing']
#         idx_new = (idx[0], idx[1], pointing)
#         if idx_new not in file_log.index:
#             print(f'WARNING {idx_new} NOTTTT THERE!!!!!')
#         else:
#             file_log.loc[idx_new, 'z_manual'] = redshift_manual_df.loc[idx, 'z_manual']
#             file_log.loc[idx_new, 'z_tier'] = 3
#
# Save the log
lime.save_frame(tables_folder/f'{sample}_files_log.txt', file_log)

#
# fname = '/home/vital/Astrodata/KelceySample/rubies_prism_files.txt'
# file_arr = np.loadtxt(fname, dtype=str)
# print(type(file_arr))
