from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, rc_context

import lime
from support.tools import create_backup, capers_load_function

# Load configuration
cfg_file = '../CAPERS_v3.toml'
capers_cfg = lime.load_cfg(cfg_file)

# Get sample information
version = capers_cfg['meta']['version']
sample = capers_cfg['meta']['sample']
pref_ext_list = capers_cfg['file_structure']['ext_1d_fits']
sample_list = capers_cfg['meta']['core_sample_list']

observations_folder = Path(capers_cfg['file_structure']['observations_folder'])
data_folder = Path(capers_cfg['file_structure']['data_folder'])
source_folder = Path(capers_cfg['file_structure']['source_folder'])
norm_flux = capers_cfg['data']['norm_flux']

file_sample_fname = source_folder/capers_cfg['file_structure']['file_sample_address']
create_backup(file_sample_fname)
files_sample = lime.Sample(file_sample_fname, levels=["sample", "id", "pointing"],
                           load_function=capers_load_function, norm_flux=norm_flux, folder_obs=observations_folder)

flux_sample_fname =source_folder/capers_cfg['file_structure']['flux_sample_address']
create_backup(flux_sample_fname)
flux_sample = lime.Sample(flux_sample_fname, levels=["sample", "id", "file", "line"],
                        load_function=capers_load_function, norm_flux=norm_flux, folder_obs=observations_folder)

# Redshift table
z_df = lime.redshift_calculation(flux_sample, weight_parameter='profile_flux')

# Get the Objects which have issues
idcs_problems = ((z_df.z_mean < 0) | (z_df.z_mean > 15) | (z_df.z_std > 0.5))
conflicting_files = z_df.loc[idcs_problems].index.get_level_values('id').to_numpy()

print(f'Conflicting files: {conflicting_files.size} / {z_df.index.size}')
print(z_df.loc[idcs_problems, ['z_mean', 'z_std']])
print(f'- List: \n{np.sort(conflicting_files).tolist()}')

# Format and add more data:
# z_df.sort_values(['id', 'file'], inplace=True)
# z_df.rename(columns = {'z_mean': 'z_centroid', 'z_std': 'z_centroid_err'}, inplace=True)
# z_df.insert(2, 'z_tier', 0)
# z_df['z_manual'] = np.nan
# lime.save_frame(source_folder/'capers_redshift_frame_v0.txt', z_df)

# Adjust the files
out_df = files_sample.frame
columns_dead = ['z_med', 'disp', 'alias', 'z_best',
    "n_nods_vis1", "n_nods_vis2", "n_nods_vis3", "eff_exp_time", "shutter_centering", 'Notes',
    "z_tier", "n_emission", "n_doublet", 'MSA_weight',  'z_aspect_key',  'z_aspect_xor',  'z_aspect_points',  'z_bands',
    "n_cosmic-ray", "n_absorption", "n_dead-pixel", "n_lines", "conf_emission", "conf_doublet",
    "conf_cosmic-ray", "conf_absorption", "conf_dead-pixel", "file_path", "flux_f090w",
    "fluxerr_f090w", "flux_f105w", "fluxerr_f105w", "flux_f115w", "fluxerr_f115w",
    "flux_f125w", "fluxerr_f125w", "flux_f140w", "fluxerr_f140w", "flux_f150w",
    "fluxerr_f150w", "flux_f160w", "fluxerr_f160w", "flux_f200w", "fluxerr_f200w",
    "flux_f277w", "fluxerr_f277w", "flux_f356w", "fluxerr_f356w", "flux_f410m",
    "fluxerr_f410m", "flux_f435w", "fluxerr_f435w", "flux_f444w", "fluxerr_f444w",
    "flux_f606w", "fluxerr_f606w", "flux_f814w", "fluxerr_f814w", "confidence_emission",
    "confidence_cosmic-ray", "confidence_doublet", "confidence_absorption", "confidence_dead-pixel"]
out_df = out_df.drop(columns_dead, axis=1)
out_df.rename(columns = {'z_gaussian': 'z_centroid', 'z_gaussian_err': 'z_centroid_err'}, inplace=True)
out_df['measurements_fname'] = 'no_measurements'
out_df['line_list'] = 'none'

# Loop through the measurements and store them into the table
out_df['zq_quality'] = 0
out_df['zq_quality'] = out_df['zq_quality'].astype(int)
for idx in z_df.index:
    sample, obj_id, fname = idx
    idx_target = out_df[(out_df['x1d'] == fname) | (out_df['optext'] == fname)].index
    if len(idx_target) == 1:
        out_df.loc[idx_target, ['z_centroid', 'z_centroid_err']] = z_df.loc[idx, ['z_mean', 'z_std']].to_numpy()
        out_df.loc[idx_target, 'line_list'] = z_df.loc[idx, 'lines']
        out_df.loc[idx_target, 'measurements_fname'] = z_df.loc[idx].name[2]

        # Assign the quality
        n_lines = len(z_df.loc[idx, 'lines'].split(','))
        if n_lines == 1:
            z_quality = 9
        elif n_lines >= 3:
            z_quality = 4
        else:
            z_quality = 3
        out_df.loc[idx_target, 'zq_quality'] = z_quality

    else:
        raise KeyError(f'Multiple objects  for {obj_id}')


# Loop through the files with only manual measurements
mask = out_df['z_centroid'].isna() & out_df['z_manual'].notna()
idx_tier_2 = out_df[mask]
for idx_target in out_df[mask].index:
    out_df.loc[idx_target, 'zq_quality'] = 2

# Moving columns
col_data = out_df.pop('z_UNICORN')         # removes and returns the column
out_df.insert(6, 'z_UNICORN', col_data)    # inserts it at position i


# Histogram sample
title = (f'CAPERS (UDS_V0.1, COSMOS_V0.2.1, EGS_V0.2.2) line centroid redshifts\n'
         f'({out_df["z_centroid"].notnull().sum()}/{out_df.index.size} galaxies)')

with rc_context(lime.theme.fig_defaults()):
    fig, ax = plt.subplots()
    ax.hist(out_df['z_centroid'], edgecolor='black')  # Customize bins as needed
    ax.set_xlabel('Redshift')  # Adding x-axis label
    ax.set_ylabel('Frequency')
    ax.set_title(title)  # Adding title
    plt.savefig(source_folder/f'CAPERs_LineRedshifts_v0.png', bbox_inches='tight')
    # plt.show()


# Add object notes
out_df['Notes'] = 'none'
for sample_name in ['COSMOS', 'UDS', 'EGS']:
    fname = source_folder/f'CAPERS_{sample_name}_notes.txt'
    df_notes = pd.read_csv(fname, sep=',', header=0, index_col=0)
    for id_label in df_notes.index:
        mask = (out_df.index.get_level_values('id') == id_label) & (out_df.index.get_level_values('sample').str.contains(sample_name))
        matching_index = out_df.index[mask]
        out_df.loc[matching_index, 'Notes'] = df_notes.loc[id_label, 'Notes']
        # Or get just the matching indexes

# Save the file
lime.save_frame(source_folder/'CAPERS_fields_LineRedshifts_v0.txt', out_df)
lime.save_frame(source_folder/'CAPERS_fields_LineRedshifts_v0.csv', out_df)

