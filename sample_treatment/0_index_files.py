import lime
import numpy as np
import pandas as pd
from support.tools import search_spectra_capers, create_backup
from pathlib import Path

# Load configuration
cfg_file = '../CAPERS_v3.toml'
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
if sample == 'CAPERS_UDS_V0.1':
    redshift_df = pd.read_csv(source_folder/capers_cfg['file_structure'][f'{sample}_redshift_file'], delimiter=',', index_col=0)
else:
    redshift_df = pd.read_csv(source_folder/capers_cfg['file_structure'][f'{sample}_redshift_file'], delimiter=';', index_col=0)

# Search the file names
raw_spectra_list = search_spectra_capers(observations_folder/sample, ext_list=fits_ext_list)

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

    if 'EGS' in sample:
        dir_parts = fits_file.parts[1:-1]
        file_dir = Path(*dir_parts)
        file_name = fits_file.name
        file_path = Path('/'.join(dir_parts[-2:]))

        name_comps = file_name.split('_')
        disp = 'prism'
        multipointing_check = True if 'MultiPointing' in fits_file.as_posix() else False
        pointing_name = 'MultiPointing' if multipointing_check else name_comps[2]

        id_label = name_comps[3] if multipointing_check is False else name_comps[2]
        MPT_number = int(id_label[1:])
        ext = name_comps[-1][0:-5]

        file_log.loc[(sample, id_label, pointing_name), "MPT_number"] = MPT_number
        file_log.loc[(sample, id_label, pointing_name), "disp"] = disp
        file_log.loc[(sample, id_label, pointing_name), "z_tier"] = 0
        file_log.loc[(sample, id_label, pointing_name), ext] = file_name
        file_log.loc[(sample, id_label, pointing_name), "file_path"] = file_path.as_posix()

    if 'COSMOS' in sample:
        dir_parts = fits_file.parts[1:-1]
        file_dir = Path(*dir_parts)
        file_name = fits_file.name
        file_path = Path('/'.join(dir_parts[-2:]))

        name_comps = file_name.split('_')
        disp = 'prism'
        multipointing_check = True if 'MultiPointing' in fits_file.as_posix() else False
        pointing_name = 'MultiPointing' if multipointing_check else name_comps[2]

        id_label = name_comps[3] if multipointing_check is False else name_comps[2]
        MPT_number = int(id_label[1:])
        ext = name_comps[-1][0:-5]

        # name_comps = file_name.split('_')
        # disp = 'prism'
        # pointing_name = name_comps[2]
        #
        # id_label = name_comps[3]
        # MPT_number = int(id_label[1:])
        # ext = name_comps[-1][0:-5]

        file_log.loc[(sample, id_label, pointing_name), "MPT_number"] = MPT_number
        file_log.loc[(sample, id_label, pointing_name), "disp"] = disp
        file_log.loc[(sample, id_label, pointing_name), "z_tier"] = 0
        file_log.loc[(sample, id_label, pointing_name), ext] = file_name
        file_log.loc[(sample, id_label, pointing_name), "file_path"] = file_path.as_posix()

    if 'UDS' in sample:
        dir_parts = fits_file.parts[1:-1]
        file_dir = Path(*dir_parts)
        file_name = fits_file.name
        file_path = Path('/'.join(dir_parts[-2:]))

        name_comps = file_name.split('_')
        disp = 'prism'
        pointing_name = name_comps[2]

        id_label = name_comps[3]
        MPT_number = int(id_label[1:])
        ext = name_comps[-1][0:-5]

        file_log.loc[(sample, id_label, pointing_name), "MPT_number"] = MPT_number
        file_log.loc[(sample, id_label, pointing_name), "disp"] = disp
        file_log.loc[(sample, id_label, pointing_name), "z_tier"] = 0
        file_log.loc[(sample, id_label, pointing_name), ext] = file_name
        file_log.loc[(sample, id_label, pointing_name), "file_path"] = file_path.as_posix()


# Adjust and organize the dataframe
file_log.sort_values(by=['MPT_number'], ascending=[True], inplace=True)

# Add CAPERS photometric predictions
MPT_list = redshift_df.index.unique()
redshift_headers  = redshift_df.columns.to_numpy()
missing_objects = []
for MPT in MPT_list:
    idcs_obj = file_log.MPT_number == MPT
    if np.any(idcs_obj):
        for hdr in redshift_headers:
            if hdr in file_log.columns:
                formatted_entry = str(redshift_df.loc[MPT, hdr]).replace(' ','')
                if 'flux_' not in hdr:
                    file_log.loc[idcs_obj, hdr] = formatted_entry
                else:
                    file_log.loc[idcs_obj, hdr.lower()] = formatted_entry
    else:
        missing_objects.append(MPT)

if len(missing_objects) > 0:
    print(f"These objects are not on the {capers_cfg['file_structure'][f'{sample}_redshift_file']} database:")
    print(f'{np.sort(missing_objects)}')

# Reassign previous redshifts
if sample == 'CAPERS_EGS_V0.2.2':
    previous_version_fname = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/sample/CAPERS_EGS_V0.2.2/tables/CAPERS_EGS_V0.2.1_files_log_the_source.txt'
    redshift_manual_df = lime.load_frame(previous_version_fname, levels=["sample", "id", "pointing"])

    for idx in redshift_manual_df.index:
        idx_new = ('CAPERS_EGS_V0.2.2', idx[1], idx[2])

        if idx_new in file_log.index:
            if redshift_manual_df.loc[idx, 'z_manual'] is not None:
                file_log.loc[idx_new, 'z_manual'] = redshift_manual_df.loc[idx, 'z_manual']
                file_log.loc[idx_new, 'z_tier'] = 3
            # else:
            #     file_log.loc[idx_new, 'z_manual'] = -1
            #     file_log.loc[idx_new, 'z_tier'] = -1
        else:
            print(f'WARNING {idx} NOTTTT THERE!!!!!')

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
# if sample == 'CAPERS_COSMOS_V0.2.1':
#     previous_sample = 'CAPERS_COSMOS_V0.2'
#     previous_version_fname = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/sample/CAPERS_COSMOS_V0.2/tables/CAPERS_COSMOS_V0.2_files_log.txt'
#     redshift_manual_df = lime.load_frame(previous_version_fname, levels=["sample", "id", "pointing"])
#     idcs_manual = redshift_manual_df.loc[pd.notnull(redshift_manual_df.z_manual)].index
#     for idx_old in idcs_manual:
#         pointing =  redshift_manual_df.loc[idx_old].name[2]# redshift_manual_df.loc[idx, 'pointing']
#         idx_new = (sample, idx_old[1], pointing)
#         if idx_new not in file_log.index:
#             print(f'WARNING {idx_new} NOTTTT THERE!!!!!')
#         else:
#             file_log.loc[idx_new, 'z_manual'] = redshift_manual_df.loc[idx_old, 'z_manual']
#             file_log.loc[idx_new, 'z_tier'] = 3

# Save the log
lime.save_frame(tables_folder/f'{sample}_files_log.txt', file_log)
