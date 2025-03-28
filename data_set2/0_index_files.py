import lime
import pandas as pd
from support.tools import search_spectra_capers, create_backup
from pathlib import Path

# Load configuration
cfg_file = '../CAPERS_v2.toml'
capers_cfg = lime.load_cfg(cfg_file)

# Get sample information
version = capers_cfg['meta']['version']
sample = capers_cfg['meta']['sample']

observations_folder = Path(capers_cfg['file_structure']['observations_folder'])
data_folder = Path(capers_cfg['file_structure']['data_folder'])
tables_folder = data_folder/sample/'tables'
source_folder = data_folder/sample/'source'

fits_ext_list = capers_cfg['file_structure']['fits_ext_list']

headers_id = capers_cfg['file_structure']['headers_id']
headers_sample = capers_cfg['file_structure']['headers_sample']
headers_analysis = capers_cfg['file_structure']['headers_analysis']

# Create day back-ups
create_backup(tables_folder/f'{sample}_files_log.txt')

# Load the sample data table
redshift_df = pd.read_csv(capers_cfg['file_structure']['CAPERs_redshift_file'], delimiter=',', index_col=0)

# Search the file names
raw_spectra_list = search_spectra_capers(observations_folder/sample, ext_list=fits_ext_list)

# Adjust the selection
# spectra_list = review_masked_files(raw_spectra_list)
spectra_list = raw_spectra_list

# Empty frame for the sample data
file_log = pd.DataFrame(columns=headers_id+headers_sample+headers_analysis)
file_log.set_index(['sample', 'id', 'file'], inplace=True)

# Default values

# Loop through the files and add the file name data
print(f'- Filling files log')
for idx_spec, fits_file in enumerate(spectra_list):

    # print(idx_spec, fits_file)
    dir_parts = fits_file.parts[1:-1]
    file_dir = Path(*dir_parts)
    file_name = fits_file.name

    file_path = Path('/'.join(fits_file.parts[5:-1]))/file_name
    pointing_dir = dir_parts[5]

    name_comps = file_name.split('_')
    disp = 'prism'
    pointing_name = name_comps[2]

    # MSA = int(id_label)
    id_label = name_comps[3]
    MPT_number = int(id_label[1:])
    ext = name_comps[-1][0:-5]

    if pointing_name != pointing_dir:
        print(f'Missmatch folder: {pointing_dir} != {file_path}')

    file_log.loc[(sample, id_label, file_name), "MPT_number"] = MPT_number
    file_log.loc[(sample, id_label, file_name), "disp"] = disp
    file_log.loc[(sample, id_label, file_name), "pointing"] = pointing_name
    file_log.loc[(sample, id_label, file_name), "ext"] = ext
    file_log.loc[(sample, id_label, file_name), "file_path"] = file_path.as_posix()
    file_log.loc[(sample, id_label, file_name), "z_tier"] = 0

# Adjust and organize the dataframe
file_log.sort_values(by=['z_UNICORN', 'ext'], ascending=[False, True], inplace=True)

# Add previous measurements
MPT_list = redshift_df.index.unique()
for MPT in MPT_list:
    idcs_obj = file_log.MPT_number == MPT
    file_log.loc[idcs_obj, 'z_med'] = redshift_df.loc[MPT, 'z_med']
    file_log.loc[idcs_obj, 'z_UNICORN'] = redshift_df.loc[MPT, 'z_UNICORN']
    file_log.loc[idcs_obj, 'ra'] = redshift_df.loc[MPT, 'ra']
    file_log.loc[idcs_obj, 'dec'] = redshift_df.loc[MPT, 'dec']

    for hdr in headers_sample:
        file_log.loc[idcs_obj, hdr] = redshift_df.loc[MPT, hdr]

# # Reasign previous redshifts
# previous_version_fname = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/tables/CAPERS_UDS_V0.1_files_log_ORIGINAL.txt'
# redshift_manual_df = lime.load_frame(previous_version_fname, levels=["sample", "id", "file"])
# idcs_manual = redshift_manual_df.loc[pd.notnull(redshift_manual_df.z_manual)].index
# for idx in idcs_manual:
#     if idx in file_log.index:
#         file_log.loc[idx, 'z_manual'] = redshift_manual_df.loc[idx, 'z_manual']

# redshift_df = pd.read_csv(source_folder/ceers_redshift_name, delimiter=';', index_col=0)
# MSA_list = redshift_df.index.unique()
# for MSA in MSA_list:
#
#     idcs_obj = file_log.index.get_level_values('id') == MSA
#     file_log.loc[idcs_obj, 'z_phot'] = redshift_df.loc[MSA, 'old_z_phot']
#     file_log.loc[idcs_obj, 'ra'] = redshift_df.loc[MSA, 'ra']
#     file_log.loc[idcs_obj, 'dec'] = redshift_df.loc[MSA, 'dec']

# # Save the log
# lime.save_frame(tables_folder/f'{sample}_{version}_files_log.txt', file_log)
lime.save_frame(tables_folder/f'{sample}_files_log.txt', file_log)
