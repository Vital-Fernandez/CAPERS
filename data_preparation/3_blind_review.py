from pathlib import Path

from matplotlib import pyplot as plt, colors, rc_context
import numpy as np
import pandas as pd
import lime

lime.theme.set_style('dark')

from support.tools import create_backup, nirspec_load_function, save_detection_results, read_prism_r_curve

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
norm_flux = capers_cfg['data']['norm_flux']

# Load the sample table (create daily back-up if necessary)
log_address = tables_folder/f'{sample}_files_log.txt'
create_backup(log_address)
files_df = lime.load_frame(log_address, levels=["sample", "id", "file"])
files_df['Notes'] = ' '

# Table with the notes
fname = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/tables/table_spectra_notes.csv'
notes_df = pd.read_csv(fname, delimiter=';', header=0)
columns_to_clean = ['Spec ID', 'Notes']
notes_df[columns_to_clean] = notes_df[columns_to_clean].apply(lambda col: col.str.replace('\t', '', regex=False))
idcs_cont_fitting = notes_df['Notes'].str.contains('strong continuum|some continuum|continuum signal', case=False)
idcs_weird = ~idcs_cont_fitting
notes_df['Notes'] = notes_df['Notes'].str.replace(' ', '_')

lime.save_frame('/home/vital/Dropbox/Astrophysics/Data/CAPERS/tables/CAPERS_UDS_V01_objects_weak_emission_v2.txt', notes_df.loc[idcs_cont_fitting])

idcs_selection = pd.notnull(files_df["z_manual"]) & (files_df["ext"] != 's2d')
files_df = files_df.loc[idcs_selection]
files_df = files_df.drop_duplicates(subset='MPT_number', keep='first')

idcs_aspect = (files_df.ext.isin(['optext', 'x1d'])) & (pd.notnull(files_df["z_manual"]))
n_galaxies = files_df.loc[idcs_aspect].index.get_level_values('id').unique().size

# Creating the figure and axes
with rc_context(lime.theme.fig_defaults()):
    fig, ax = plt.subplots()
    ax.hist(files_df['z_aspect_brute'], edgecolor='black')  # Customize bins as needed
    ax.set_title(f'{sample} redshift histogram (v2): \n {n_galaxies} galaxies (ASPECT_v0 + 2nd review)')  # Adding title
    ax.set_xlabel('Redshift')  # Adding x-axis label
    ax.set_ylabel('Frequency')  # Adding y-axis label
    plt.show()

# Filter rows where A is not null
df_filtered = files_df[files_df['z_manual'].notnull()]
df_sorted = df_filtered.sort_values(by=['ext'])
result = df_sorted[~df_sorted.index.duplicated(keep='first')]

# Add the notes
for idx in notes_df.loc[idcs_weird, 'Spec ID'].index:
    mpt_id, comments = int(notes_df.loc[idx, 'Spec ID'][1:]), notes_df.loc[idx, 'Notes']
    idx_obj = result.MPT_number == mpt_id
    result.loc[idx_obj, 'Notes'] = comments

lime.save_frame('/home/vital/Dropbox/Astrophysics/Data/CAPERS/tables/CAPERS_UDS_V01_ASPECT_redshifts_v2.txt', result[['z_manual', 'Notes']])


# # Add the rows
# for idx in notes_df.loc[idcs_weird, 'Spec ID'].index:
#     mpt_id, comments = int(notes_df.loc[idx, 'Spec ID'][1:]), notes_df.loc[idx, 'Notes']
#     idx_obj = files_df.MPT_number == mpt_id
#     files_df.loc[idx_obj, 'Notes'] = comments
#
# files_df['z_manual'] = files_df['z_manual'].round(2)
# result = files_df[~files_df.index.duplicated(keep='first')]
# lime.save_frame('/home/vital/Dropbox/Astrophysics/Data/CAPERS/tables/CAPERS_UDS_V01_ASPECT_redshifts_v0.1.txt', result[['z_manual', 'Notes']])


# highest_value_index = files_df['z_manual'].idxmax()

# idx_obj = (files_df.MPT_number == 7348) & (files_df.ext == 'optext')
# idx_obj = files_df.loc[idx_obj].index[0]
# files_sample = lime.Sample(log_address, levels=["sample", "id", "file"], load_function=nirspec_load_function,
#                            norm_flux=norm_flux, folder_obs=observations_folder)
# spec = files_sample.get_spectrum(idx_obj, redshift=files_sample.loc[idx_obj, 'z_manual'])
# spec.plot.spectrum(rest_frame=True)
# print(files_df.columns)


# idcs = (files_df.disp == 'prism') & (files_df.ext == 'x1d')
# n_galaxies = files_df.loc[idcs].index.get_level_values('id').unique().size
# n_no_redshifts = files_df.loc[idcs & pd.isnull(files_df.z_phot)].index.get_level_values('id').unique().size
# print("Full CEERs number", n_galaxies)
# print("Negative redshift measurement", n_no_redshifts)
# print("Positive redshift measurement", n_galaxies - n_no_redshifts)
#
# z_target = "z_aspect_brute"
# idcs_aspect = (files_df.disp == 'prism') & (pd.notnull(files_df[z_target]))
# n_galaxies = files_df.loc[idcs_aspect].index.get_level_values('id').unique().size
# print("Full ASPECT measurement", n_galaxies)
#
#
# # Selection
# idcs_true_brute = pd.notnull(files_df.z_aspect_brute) & pd.notnull(files_df.z_phot) & (files_df.n_lines >= 1)
#
# x = files_df.loc[idcs_true_brute, 'z_aspect_brute_bound'].to_numpy()
# y = files_df.loc[idcs_true_brute, 'z_phot'].to_numpy()
# zq = files_df.loc[idcs_true_brute, 'zq_best'].to_numpy()
#
# unique_values = np.unique(zq)
# cmap = plt.cm.jet  # Choose a colormap
# norm = colors.Normalize(vmin=unique_values.min(), vmax=unique_values.max())
#
# x_guide = np.linspace(0, 12, 20)
# y_guide = x_guide
#
# fig, ax = plt.subplots()
# ax.plot(x_guide, y_guide, color='black')
# sc = ax.scatter(x, y, c=zq, cmap=cmap, norm=norm)
# cbar = plt.colorbar(sc, ticks=unique_values)
# cbar.ax.set_yticklabels([str(val) for val in unique_values])  # Ensure correct labels
# cbar.set_label('Quality label')
#
# ax.set_xlim(0, 12)
# ax.set_ylim(0, 12)
# ax.update({'title':f'CEERs and DDT redshift comparison ({idcs_true_brute.sum()} galaxies)', 'xlabel':'ASPECT redshift',
#            'ylabel': 'True redshift'})
# plt.tight_layout()
# plt.show()
#


# idcs_selection = np.isclose(files_df.z_aspect_brute, files_df.z_phot, atol=0.10) & (files_df.n_lines >= 1)
# # idcs_selection = ~idcs_selection & (pd.notnull(files_df.z_aspect_brute)) & (files_df.z_aspect_brute > files_df.z_phot)
#
# # idcs_selection = ~idcs_selection & (pd.notnull(files_df.z_aspect_brute)) & (files_df.z_aspect_brute > files_df.z_phot)
# # idcs_selection = ~idcs_selection & np.isclose(files_df.z_aspect_brute, files_df.z_phot, atol=1.00) LOW REDSHIFT MISSMATCH
#
# # idcs_selection = ~idcs_selection & (files_df.z_aspect_brute > 1.1 * files_df.z_phot) & (files_df.z_aspect_brute < 2.1 * files_df.z_phot)
#
# print(list(files_df.loc[idcs_selection].index.get_level_values('id').to_numpy()))
#
# x = files_df.loc[idcs_selection, 'z_aspect_brute']
# y = files_df.loc[idcs_selection, 'z_phot']
# x_guide = np.linspace(0, 12, 20)
# y_guide = x_guide
#
# fig, ax = plt.subplots()
# ax.plot(x_guide, y_guide, color='black')
# ax.scatter(x, y)
# ax.set_xlim(0, 12)
# ax.set_ylim(0, 12)
# ax.update({'title':f'CEERs and DDT redshift comparison ({idcs_selection.sum()} galaxies)', 'xlabel':'ASPECT redshift',
#            'ylabel': 'True redshift'})
# plt.tight_layout()
# plt.show()