from pathlib import Path
import numpy as np
import pandas as pd
import lime

from support.tools import (create_backup, capers_load_function, save_detection_results, read_prism_r_curve,
                           z_selection)

lime.theme.set_style('dark')
lime.lineDB.set_database(vacuum_waves=True)

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

# Treatment selection
start_object = 0 #27 28 30
ASPECT_CHECK = False
ASPECT_REDSHIFT = False
REDSHIFT_CHECK = False
MEASURE_CHECK = True
BANDS_CHECK = False

# Determine the objects to treat
tier_df = files_sample.frame.sort_values(by='MPT_number', ascending=[False])
obj_list = tier_df.index.get_level_values('id').unique()

#EGS negative profile [7303, 22465, 33856, 35518, 36864, 92582, 92742]
idcs = tier_df['MPT_number'].isin([177138])
obj_list = tier_df.loc[idcs].index.get_level_values('id').unique()

# Loop through the objects
counter = 0 # DONT TOUCH
saving_iter_n = 20
failing_opening = {}
error_objects = {}
object_counter = 0
for i, obj in enumerate(obj_list):

    # Get list of its observations and their extensions
    idcs_obj_group = files_sample.loc[files_sample.index.get_level_values('id') == obj].index

    # # Exclude the pointings already done
    # mask = files_sample.loc[idcs_obj_group].index.get_level_values('pointing').isin(['P6', 'MultiPointing'])
    # idcs_obj_group = idcs_obj_group[mask]

    # Loop through the observations
    if  idcs_obj_group.size > 0:
        for j, idx_obj in enumerate(idcs_obj_group):
            if counter >= start_object:

                # Objects ID params
                MPT = idx_obj[1]
                pointing = idx_obj[2]
                disp = files_sample.loc[idx_obj, 'disp']
                if pd.notnull(files_sample.loc[idx_obj, "z_UNICORN"]): z_phot = files_sample.loc[idx_obj, "z_UNICORN"]
                else: z_phot = files_sample.loc[idx_obj, "z_med"]

                # Redshift value
                z_tier, z_obj = z_selection(files_sample, idx_obj)

                # Choose the file by aperture type
                fname_list = files_sample.loc[idx_obj, pref_ext_list].values
                fname = fname_list[0] if pd.notnull(fname_list[0]) else fname_list[1]
                file_path = Path(files_sample.loc[idx_obj, 'file_path'])/fname
                fits_stem = file_path.stem

                # Output file names
                obj_comps_file = comps_folder / f'{fits_stem}_components.txt'
                obj_bands_file = bands_folder / f'{fits_stem}_bands.txt'
                obj_logs_file = logs_folder/ f'{fits_stem}_log.txt'

                # Fitting configuration for these objects
                obj_cfg_section = f'{MPT}_{fits_stem}' if capers_cfg.get(f'{MPT}_{fits_stem}_line_fitting') \
                                                       else "boundary_prism"

                if z_obj is not None:
                    composite_key = 'line_composites_high_z' if z_obj > 5 else 'line_composites_medium_z'
                    composite_lines = capers_cfg['lines_data'][composite_key]
                else:
                    composite_lines = None

                # Load spectrum
                spec = files_sample.get_spectrum(idx_obj, redshift=z_obj, detection_file=obj_comps_file,
                                                 file_path=file_path)

                # Check for full nan spectra
                if spec is None: failing_opening[MPT] = file_path; continue

                # Prompt
                # print(f'\n{counter}: MPT {MPT}) Disp = {disp}, pointing = {pointing}, z_obj={z_obj} (z_tier={z_tier})')
                #
                # # Aspect detection
                # if ASPECT_CHECK:
                #     spec.infer.components(show_steps=False, exclude_continuum=True)
                #     save_detection_results(spec, idx_obj, files_sample, obj_comps_file)
                #
                # # Redshift from lines intervals (n_lines > 1)
                # if ASPECT_REDSHIFT:
                #     if files_sample.loc[idx_obj, 'n_lines'] >= 1:
                #
                #         # Compute the redshifts
                #         limits = (np.maximum(z_phot-2,0), np.minimum(z_phot+2, 16)) if z_phot < 3 else (None, None)
                #         z_key = spec.fit.redshift(redshift_bands,  sigma_factor=1, z_min=limits[0], z_max=limits[1],
                #                                   mode='key', plot_results=False,)
                #         z_xor = spec.fit.redshift(redshift_bands, sigma_factor=1, z_min=limits[0], z_max=limits[1],
                #                                   mode='xor', plot_results=False)
                #
                #         # Check redshift compatibility between techniques
                #         z_diff_check, z_tier_aspect = False, 1
                #         if z_key is not None and z_xor is not None:
                #             percent_diff = np.abs(z_key / z_xor - 1) * 100
                #             z_diff_check = percent_diff <= 5
                #             z_tier_aspect = 2 if z_diff_check else 1
                #
                #         # Save the measurements:
                #         files_sample.loc[idx_obj, 'z_aspect_key'] = z_key
                #         files_sample.loc[idx_obj, 'z_aspect_xor'] = z_xor
                #         if z_tier < z_tier_aspect:
                #             files_sample.loc[idx_obj, 'z_tier'] = z_tier_aspect
                #
                #         # Prompt
                #         print(f'- Aspect Check: z_phot = {z_phot}')
                #         print(f'-- Key: z_key = {z_key}')
                #         print(f'-- Key: z_xor = {z_xor}')
                #         print(f'--- Compatible: {z_diff_check}')
                #
                # # Redshift review
                # if REDSHIFT_CHECK:
                #
                #     z_tier, z_obj = z_selection(files_sample, idx_obj)
                #     print(f'- Manual check: z_manual = {z_obj}')
                #
                #     title = f'MPT {obj} ({disp}) z_photo={z_phot:0.3f}; '
                #     output_idcs = files_sample.index.get_level_values('id') == obj
                #     files_sample.check.redshift(idx_obj, lines_visualize,
                #                                 initial_z=z_obj, redshift_column='z_manual',
                #                                 title=title, legend_handle='optext', maximize=True,
                #                                 output_file_log=log_address, file_path=file_path)
                #
                #     # Establish the tier
                #     if pd.notnull(files_sample.loc[idx_obj, 'z_manual']):
                #         files_sample.loc[idx_obj, 'z_tier'] = 3


                if MEASURE_CHECK:

                    #Retrieve redshift
                    z_tier, z_obj = z_selection(files_sample, idx_obj)

                    if z_tier >= 3:

                        # Use the bands file available
                        if obj_bands_file.is_file():
                            spec.fit.frame(obj_bands_file, fit_cfg=capers_cfg, obj_cfg_prefix=obj_cfg_section,
                                           default_cfg_prefix='default_prism', cont_from_bands=False, err_from_bands=False)
                            spec.save_frame(obj_logs_file, skip_failed=True)

                        # Compute the object bands and measure them:
                        else:

                            line_list = ['O2_3726A', 'O2_3729A', 'He1_4471A', 'H1_4861A', 'O3_4959A', 'O3_5007A', 'He1_5876A',
                                         'H1_6563A', 'N2_6583A', 'N2_6548A', 'He1_10832A', 'Fe2_12570A']
                            object_bands = spec.retrieve.line_bands(line_list=line_list, fit_cfg=capers_cfg,
                                                                    **capers_cfg['bands_generation_parameters'])
                            line_list = ['H1-O3_4861A_b']
                            spec.fit.frame(object_bands, fit_cfg=capers_cfg, obj_cfg_prefix=obj_cfg_section,
                                           default_cfg_prefix='default_prism',
                                           # line_list=line_list,
                                           cont_from_bands=False, err_from_bands=False)
                            spec.save_frame(obj_logs_file, skip_failed=True)
                            spec.plot.spectrum(rest_frame=True, show_cont=False, show_err=True)

#                             # object_bands = spec.retrieve.line_bands(ref_bands=prism_bands_df, fit_cfg=capers_cfg,
#                             #                                         **capers_cfg['bands_generation_parameters'])
#                             #
#                             # spec.fit.continuum(degree_list=[3, 5, 6, 7], emis_threshold=[3, 2, 2, 1.5], plot_steps=False)
#                             # match_bands = spec.infer.peaks_troughs(object_bands, emission_type=True, sigma_threshold=3,
#                             #                                        plot_steps=False, log_scale=False, maximize=True)
#                             #
#                             # # Fit the lines and review the results if there are many
#                             # spec.fit.frame(match_bands, fit_cfg=capers_cfg, obj_cfg_prefix=obj_cfg_section,
#                             #                default_cfg_prefix='default_prism', cont_from_bands=False, err_from_bands=False)
#                             #
#                             # spec.plot.spectrum(rest_frame=True, show_cont=False, show_err=True)
#
#                             # # Review measurements
#                             # if spec.frame.index.size > 0:
#                             #
#                             #     repeat_measurements = False
#                             #     if ('H1_6563A_m' in spec.frame.index) and spec.frame.loc['H1_6563A_m', 'observations'] == 'No_errorbars':
#                             #         idcs_line = match_bands.index.isin(['H1-S2_6563A_b'])
#                             #         if np.any(idcs_line):
#                             #             match_bands.rename(index={'H1-S2_6563A_b': 'H1_6563A_m'}, inplace=True)
#                             #             match_bands.loc['H1_6563A_m', 'group_label'] = capers_cfg['default_prism_line_fitting']['H1_6563A_m']
#                             #             repeat_measurements = True
#                             #
#                             #     if ('O3_5007A' in spec.frame.index) and spec.frame.loc['O3_5007A', 'observations'] == 'No_errorbars':
#                             #         idcs_line = match_bands.index.isin(['H1_4861A_b'])
#                             #         if np.any(idcs_line):
#                             #             idx_line = match_bands.loc[idcs_line].index[0]
#                             #             match_bands.rename(index={idx_line: 'H1_4861A_m'}, inplace=True)
#                             #             match_bands.loc['H1_4861A_m', 'group_label'] = capers_cfg['default_prism_line_fitting']['H1_4861A_m']
#                             #             repeat_measurements = True
#                             #
#                             #         else:
#                             #             idcs_line = match_bands.index.isin(['O3_5007A_b'])
#                             #             if np.any(idcs_line):
#                             #                 idx_line = match_bands.loc[idcs_line].index[0]
#                             #                 match_bands.rename(index={idx_line: 'O3_5007A_m'}, inplace=True)
#                             #                 match_bands.loc['O3_5007A_m', 'group_label'] = capers_cfg['default_prism_line_fitting']['O3_5007A_m']
#                             #                 repeat_measurements = True
#                             #
#                             #     if ('O3_5007A_m' in spec.frame.index) and spec.frame.loc['O3_5007A_m', 'observations'] == 'No_errorbars':
#                             #         idcs_line = match_bands.index.isin(['H1-O3_4861A_b'])
#                             #         if np.any(idcs_line):
#                             #             idx_line = match_bands.loc[idcs_line].index[0]
#                             #             match_bands.rename(index={idx_line: 'H1_4861A_m'}, inplace=True)
#                             #             match_bands.loc['H1_4861A_m', 'group_label'] = capers_cfg['default_prism_line_fitting']['H1_4861A_m']
#                             #             repeat_measurements = True
#                             #
#                             #     if repeat_measurements:
#                             #         spec.frame.drop(spec.frame.index, inplace=True)
#                             #         spec.fit.frame(match_bands, fit_cfg=capers_cfg, obj_cfg_prefix=obj_cfg_section,
#                             #                        default_cfg_prefix='default_prism', cont_from_bands=False, err_from_bands=False)
#                             #
#                             #     spec.save_frame(obj_logs_file, skip_failed=True)
#
#                 # Manually adjust the
#                 if BANDS_CHECK:
#
#                     # Retrieve redshift
#                     z_tier, z_obj = z_selection(files_sample, idx_obj)
#                     if z_tier >= 3:
#
#                         object_counter += 1
#                         print(f'-- Bands treating: MPT {MPT}) Disp = {disp}, pointing = {pointing}')
#
#                         if obj_logs_file.is_file():
#                             spec.load_frame(obj_logs_file)
#                             spec.plot.spectrum(rest_frame=True, maximize=True)
#
#                         # Use the bands file available
#                         if obj_bands_file.is_file():
#                             input_bands = obj_bands_file
#
#                         else:
#                             object_bands = spec.retrieve.line_bands(ref_bands=prism_bands_df, fit_cfg=capers_cfg,
#                                                                     **capers_cfg['bands_generation_parameters'])
#
#                             spec.fit.continuum(degree_list=[3, 5, 6, 7], emis_threshold=[3, 2, 2, 1.5],
#                                                plot_steps=False)
#
#                             input_bands = spec.infer.peaks_troughs(object_bands, emission_type=True, sigma_threshold=3,
#                                                                    plot_steps=False, log_scale=False, maximize=True)
#
#                         # Manual review
#                         print(f'[{MPT}_{fits_stem}_line_fitting]', obj_cfg_section == f'{MPT}_{fits_stem}')
#                         spec.check.bands(obj_bands_file, bands_obj=input_bands, selected_by_default=True, ref_bands=prism_bands_df,
#                                          fit_cfg=capers_cfg, exclude_continua=True, maximize=True,
#                                          **capers_cfg['bands_generation_parameters'])
#
#                         # Second fit
#                         if obj_bands_file.is_file():
#                             spec.frame.drop(spec.frame.index, inplace=True)
#                             spec.fit.frame(obj_bands_file, fit_cfg=capers_cfg, obj_cfg_prefix=obj_cfg_section,
#                                            default_cfg_prefix='default_prism', cont_from_bands=False, err_from_bands=False)
#                             spec.save_frame(obj_logs_file, skip_failed=True)
#
#                         else:
#                             spec.load_frame(obj_logs_file)
#
#                         spec.plot.spectrum(rest_frame=True, maximize=True)
#
#
#                 # # Save the frame every 20 iterations
#                 # if ((i + 1) % saving_iter_n == 0) and (counter >= start_object):
#                 #     lime.save_frame(log_address, files_sample.frame)
#                 #     print(f'\n- FILES LOG SAVED ({i})\n')
#
#             counter += 1
#
# print(f'We are to treat: {object_counter}' )
# # Save the dataframe with the ASPECT redshifts
# lime.save_frame(log_address, files_sample.frame)

