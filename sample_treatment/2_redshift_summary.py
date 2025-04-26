from pathlib import Path
import numpy as np
import pandas as pd
import lime
from matplotlib import pyplot as plt, rc_context
from support.tools import create_backup, capers_load_function, read_prism_r_curve

# lime.theme.set_style('dark')

# Load configuration
cfg_file = '../CAPERS_v3.toml'
capers_cfg = lime.load_cfg(cfg_file)

# Get sample information
version = capers_cfg['meta']['version']
sample = capers_cfg['meta']['sample']
pref_ext_list = capers_cfg['file_structure']['ext_1d_fits']

observations_folder = Path(capers_cfg['file_structure']['observations_folder'])
data_folder = Path(capers_cfg['file_structure']['data_folder'])
source_folder = Path(capers_cfg['file_structure']['source_folder'])

# Output directories
tables_folder = data_folder/sample/'tables'
bands_folder = data_folder/sample/'bands'
logs_folder = data_folder/sample/'line_logs'
comps_folder = data_folder/sample/'comps'

# Spectra parameters
prism_bands_df = lime.load_frame(source_folder/'CAPERs_lines_database.txt')
redshift_lines = capers_cfg['lines_data']['lines_redshift']
lines_visualize = capers_cfg['lines_data']['lines_visualize']
R_interpolator = read_prism_r_curve(source_folder/capers_cfg['file_structure']['prism_R_curve_file'], units_factor=10000)
redshift_bands = prism_bands_df.loc[prism_bands_df.index.isin(redshift_lines)]

norm_flux = capers_cfg['data']['norm_flux']
default_cfg_section = capers_cfg['default_prism_line_fitting']

# Load the sample table (create daily back-up if necessary)
log_address = tables_folder/f'{sample}_files_log.txt'
create_backup(log_address)
files_sample = lime.Sample(log_address, levels=["sample", "id", "pointing"], load_function=capers_load_function,
                           norm_flux=norm_flux, folder_obs=observations_folder)

n_objs = files_sample.frame.MPT_number.unique().size
idcs_manual = files_sample.frame['z_manual'].notnull()
table_df = files_sample.frame.loc[idcs_manual]
title = f'{sample} \n ASPECT_v0.3 + manual inspection ({idcs_manual.sum()}/{n_objs} galaxies)'

# Histogram
# Creating the figure and axes
with rc_context(lime.theme.fig_defaults()):
    fig, ax = plt.subplots()
    ax.hist(table_df['z_manual'], edgecolor='black')  # Customize bins as needed
    ax.set_xlabel('Redshift')  # Adding x-axis label
    ax.set_ylabel('Frequency')
    ax.set_title(title)  # Adding title
    plt.savefig(source_folder/f'{sample}_redshifts_aspectv0.3_and_manual_inspection.png', bbox_inches='tight')
    # plt.show()


table_df = files_sample.frame.copy()
idcs_null_opt = pd.isnull(table_df.optext)
table_df.loc[idcs_null_opt, 'optext'] = table_df.loc[idcs_null_opt, 'x1d']
table_df.rename(columns={'optext': 'file_name'}, inplace=True)

table_df.loc[~idcs_manual, 'z_manual'] = -1
lime.save_frame(source_folder/f'{sample}_redshifts_aspectv0.3_and_manual_inspection.txt',
                table_df.loc[:, ['MPT_number', 'z_manual', 'file_name']])


# idcs_gaussian = files_sample.frame['z_gaussian'].notnull()
#
# x = files_sample.frame.loc[idcs_gaussian, 'z_manual'].to_numpy()
# y = files_sample.frame.loc[idcs_gaussian, 'z_gaussian'].to_numpy()
# y_err = files_sample.frame.loc[idcs_gaussian, 'z_gaussian_err'].to_numpy()
# fig, ax = plt.subplots()
# # ax.scatter(x, y)
# ax.errorbar(x, y, yerr=y_err, fmt='o', capsize=5)
# ax.update({'xlabel': 'z_manual', 'ylabel': 'z_gaussian'})
# plt.show()
#
# files_sample.loc[~idcs_manual, 'z_manual'] = -1
# files_sample.loc[~idcs_gaussian, 'z_gaussian'] = -1
#
# lime.save_frame('CAPERS_EGS_V0.2_lime_redshifts_v2.txt', files_sample.loc[:, ['MPT_number', 'z_manual', 'z_gaussian', 'optext']])


# files_sample.frame['id_str'] = files_sample.index.get_level_values('id')
# idcs = files_sample.frame['ext'] == 'optext'
# df_short = files_sample.frame.loc[idcs, ['id_str', 'pointing']]
# (df_short[['id_str', 'pointing']].duplicated()).any()
#
# # Plot life
# obj_list_total = files_sample.frame.index.get_level_values('id').unique()
# n_objs_total = obj_list_total.size
#
# obj_list = files_sample.frame.loc[files_sample.frame.ext == 'optext'].index.get_level_values('id').unique()
# n_objs = obj_list.size
#
# # # Create a helper column to prioritize preferred extension
# # files_sample.frame['priority'] = files_sample.frame['ext'].apply(lambda x: 0 if x == preference_list[0] else (1 if x == preference_list[1] else 2))
# # df_sorted = files_sample.frame.sort_values(['MPT_number', 'priority'])
# # df_sorted = df_sorted.drop_duplicates(subset='MPT_number', keep='first')
# # selected_indices = df_sorted.index
#
# idcs_3 = (files_sample.frame['z_tier'] == 3) & (files_sample.frame['ext'] == 'optext')
# table_df =  files_sample.frame.loc[idcs_3, ['z_tier', 'MPT_number', 'z_manual']].copy()
# table_df.sort_values(by=['MPT_number'], ascending=[True], inplace=True)
# lime.save_frame(tables_folder/f'{sample}_redshifts_aspectv0.3_inspected_optext_positive_detections.txt', table_df)
#
# # idcs_0 = (files_sample.frame['z_tier'] == 0) & (files_sample.frame['ext'] == 'optext')
# # table_df =  files_sample.frame.loc[idcs_0, ['z_tier', 'MPT_number', 'z_manual']].copy()
# # table_df.sort_values(by=['MPT_number'], ascending=[True], inplace=True)
# # lime.save_frame(tables_folder/f'{sample}_redshifts_aspectv0.3_inspected_optext_negative_detections.txt', table_df)
#
#
#
# idcs_ext = files_sample.frame['ext'] == 'optext'
# idcs_tier_dict = {i: files_sample.loc[idcs_ext, 'z_tier'] == i for i in range(5)}
# n_groups = 0
# for tier, idcs in idcs_tier_dict.items():
#     n_groups += idcs.sum()
#     print(f'z_tier = {tier}, n_obsj = {idcs.sum()}/{n_objs}')
# print(f'Total sum {n_groups} objects')
# print(f'Total {n_objs} objects')
# print(f'Total all extensions {n_objs_total} objects')
#
#
# # idcs_4 = files_sample.frame['z_tier'] == 4
# #
# # table = files_sample.frame.loc[idcs_4, ['MPT_number', 'pointing', 'z_aspect_key']]
# # table.sort_values(by=['MPT_number'], ascending=[True], inplace=True)
# # table.rename(columns={"z_aspect_key": "z_aspect"})
# # lime.save_frame(tables_folder/f'{sample}_redshifts_aspect_v0.3.txt', table)
# #
# # # # Show as list
# # #
# # Histogram
# # Creating the figure and axes
# with rc_context(lime.theme.fig_defaults()):
#     fig, ax = plt.subplots()
#     ax.hist(table_df['z_manual'], edgecolor='black')  # Customize bins as needed
#     ax.set_title(f'{sample} \n ASPECT_v0.3 + manual inspection ({idcs_tier_dict[3].sum()}/{n_objs} galaxies)')  # Adding title
#     ax.set_xlabel('Redshift')  # Adding x-axis label
#     ax.set_ylabel('Frequency')
#     plt.savefig(tables_folder/f'{sample}_redshifts_aspectv0.3_inspected_optext_positive_detections.png', bbox_inches='tight')
#     # # Adding y-axis label
#     plt.show()

