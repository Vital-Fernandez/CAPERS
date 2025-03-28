from pathlib import Path
import numpy as np
import lime
import pandas as pd

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
lines_visualize = capers_cfg['file_structure']['lines_visualize']
R_interpolator = read_prism_r_curve(capers_cfg['file_structure']['prism_R_curve_file'], units_factor=10000)
redshift_bands = prism_bands_df.loc[prism_bands_df.index.isin(redshift_lines)]

# Line fitting data
norm_flux = capers_cfg['data']['norm_flux']
blended_cfg_section = capers_cfg['blended_prism_line_fitting']
default_cfg_section = capers_cfg['default_prism_line_fitting']
boundary_cfg_section = capers_cfg['boundary_prism_line_fitting']

# Load the sample table (create daily back-up if necessary)
log_address = tables_folder/f'{sample}_files_log.txt'
create_backup(log_address)
files_sample = lime.Sample(log_address, levels=["sample", "id", "file"], load_function=nirspec_load_function,
                           norm_flux=norm_flux, folder_obs=observations_folder)

# Treatment selection
start_object = 25
ASPECT_CHECK = True
ASPECT_REDSHIFT = True
REDSHIFT_CHECK = False
BANDS_CHECK = False
MEASURE_CHECK = True
RESET_redshift = False
SAVE_CHECK = False
AUTOMATIC_MODE = False

# Determine the objects to treat
obj_list = files_sample.index.get_level_values('id').unique()
error_objects = {}

# target_objs = [3772]
# obj_list = files_sample.frame.loc[files_sample.frame.MPT_number.isin(target_objs)].index.get_level_values('id').unique()

# Loop through the objects
counter = 0 # DONT TOUCH
for i, obj in enumerate(obj_list[:100]):

    # Get list of its observations and their extensions
    idcs_obj_group = files_sample.index.get_level_values('id') == obj
    ext_list = files_sample.frame.loc[idcs_obj_group].ext.unique()

    # Get the best reduction if available
    ext_target = next((ext for ext in preference_list if ext in ext_list), None)
    if ext_target is not None:

        idcs_obj = idcs_obj_group & (files_sample.frame.ext == ext_target)
        obj_sample = files_sample[idcs_obj]

        # Loop through the observations
        if  obj_sample.size > 0:
            for j, idx_obj in enumerate(obj_sample.index):
                if counter >= start_object:

                    # Objects ID params
                    MPT = idx_obj[1]
                    disp, pointing = obj_sample.loc[idx_obj, ['disp', 'pointing']]
                    z_phot = 0 if pd.isnull(obj_sample.loc[idx_obj, "z_UNICORN"]) else obj_sample.loc[idx_obj, "z_UNICORN"]
                    fits_stem = Path(obj_sample.loc[idx_obj, 'file_path']).stem
                    id_message = f'\n{counter}: MPT {MPT} ) Disp = {disp}, pointing = {pointing} '
                    print(id_message)

                    # output object files
                    obj_comps_file = comps_folder / f'{fits_stem}_components.txt'
                    obj_bands_file = bands_folder / f'{fits_stem}_bands.txt'
                    obj_logs_file = logs_folder/ f'{fits_stem}_log.txt'
                    spec = None

                    # Aspect detection
                    if ASPECT_CHECK:
                        z_obj = 0
                        spec, err = catch_nan_spectra(files_sample, idx_obj, z_obj)

                        if spec is None:
                            print(f'- Error opening: {err}')
                            error_objects[obj] = fits_stem
                        else:
                            # spec.plot.spectrum()
                            spec.infer.components(show_steps=False, exclude_continuum=True)
                            save_detection_results(spec, idx_obj, files_sample, obj_comps_file)

                    if ASPECT_REDSHIFT:
                        if files_sample.loc[idx_obj, 'n_lines'] >= 1:
                            res_power = R_interpolator(spec.wave.data)
                            spec = spec if spec is not None else files_sample.get_spectrum(idx_obj, redshift=z_obj)

                            print(f'- Aspect Check: z_phot = {z_phot}')
                            limits = (np.maximum(z_phot-1.5,0), np.minimum(z_phot+1.5, 16)) if z_phot < 3 else (None, None)
                            z_key = spec.fit.redshift(redshift_bands, res_power=res_power, plot_results=False, sigma_factor=1,
                                                        z_min=limits[0], z_max=limits[1])
                            z_xor = spec.fit.redshift(redshift_bands, res_power=res_power, plot_results=False, sigma_factor=1,
                                                        z_min=limits[0], z_max=limits[1], mode='xor')

                            print(f'-- Key: z_key = {z_key}')
                            print(f'-- Key: z_xor = {z_xor}')

                            # Check the quality
                            if not (z_key is None) and not (z_xor is None):
                                z_diff_check = np.abs(z_key / z_xor - 1) * 100 <= 5
                                print(f'--- Difference {(z_key / z_xor - 1)*100:0.02}% (check {z_diff_check})')
                                if z_diff_check:
                                    z_tier = 2
                                else:
                                    z_tier = 1
                            else:
                                z_tier = 1

                            # Save the values:
                            files_sample.loc[idcs_obj_group, 'z_aspect_key'] = z_key
                            files_sample.loc[idcs_obj_group, 'z_aspect_xor'] = z_xor
                            files_sample.loc[idcs_obj_group, 'z_tier'] = z_tier

                            # title =f'MPT {MPT}, z_phot = {z_phot:0.2f}, z_manual = {files_sample.loc[idx_obj, ["z_manual"]].values[0]}'
                            # spec.plot.spectrum(show_categories=True, include_err=True, ax_cfg={'title': title})

                    # Redshift review
                    if REDSHIFT_CHECK:
                        z_manual, z_aspect = files_sample.loc[idx_obj, ['z_manual','z_aspect_key']]
                        z_obj = z_manual if pd.notnull(z_manual) else z_aspect
                        print(f'- Manual check: z_manual = {z_obj}')

                        title = f'MPT {obj} ({disp}) z_photo={z_phot:0.3f}; '
                        output_idcs = files_sample.index.get_level_values('id') == obj
                        idcs__group = files_sample.index.get_level_values('id') == obj

                        files_sample.check.redshift(idcs_obj, lines_visualize, output_idcs=output_idcs,
                                                    initial_z=z_obj, redshift_column='z_manual',
                                                    title=title, legend_handle='file', maximize=False,
                                                    output_file_log=log_address)

                        # Establish the tier
                        if pd.notnull(z_manual):
                            files_sample.loc[idcs_obj_group, 'z_tier'] = 3

                    if BANDS_CHECK:
                        if files_sample.loc[idcs_obj_group, 'z_tier'] == 3:
                            z_obj = files_sample.loc[idx_obj, 'z_manual']
                        elif files_sample.loc[idcs_obj_group, 'z_tier'] == 2:
                            z_obj = files_sample.loc[idx_obj, 'z_aspect_key']
                        else:
                            z_obj = None

                        # If there is a redshift measurement
                        if pd.notnull(z_obj):

                            # Load the object
                            print(f'- Bands check: z_bands = {z_obj}')
                            spec = files_sample.get_spectrum(idx_obj, redshift=z_obj, detection_file=obj_comps_file)

                            # Review the bands manually
                            spec.check.bands(obj_bands_file, ref_bands=prism_bands_df, fit_conf=blended_cfg_section,
                                             band_vsigma=400, maximize=False, exclude_continua=True, vacuum_waves=True,
                                             components_detection=True)

                            # Save the z_bands to the dataframe
                            if pd.isnull(files_sample.loc[idx_obj, 'z_bands']) and obj_bands_file.is_file():
                                files_sample.loc[idx_obj, 'z_bands'] = z_obj

                    if MEASURE_CHECK:
                        z_tier = files_sample.loc[idx_obj, 'z_tier']
                        if z_tier == 3:
                            z_obj = files_sample.loc[idx_obj, 'z_manual']
                        elif z_tier == 2:
                            z_obj = files_sample.loc[idx_obj, 'z_aspect_key']
                        else:
                            z_obj = None

                        # If there are bands
                        if pd.notnull(z_obj):

                            # If first time, prepare some data
                            spec = files_sample.get_spectrum(idx_obj, redshift=z_obj, detection_file=obj_comps_file)

                            # Recover the line fitting configuration for these objects
                            obj_cfg_section = f'{MPT}_{fits_stem}'
                            obj_cfg = capers_cfg.get(f'{obj_cfg_section}_line_fitting')
                            obj_cfg_section = "boundary_prism" if obj_cfg is None else obj_cfg_section
                            print(f'-- {MPT}_{fits_stem}_line_fitting ({obj_cfg})')

                            if obj_bands_file.is_file():
                                input_bands = lime.load_bands(obj_bands_file)
                            else:
                                input_bands = spec.retrieve.line_bands(ref_bands=prism_bands_df, fit_conf=blended_cfg_section,
                                                                       obj_conf_prefix=obj_cfg_section, band_vsigma=400,
                                                                       vacuum_waves=True,  components_detection=True)

                            spec.plot.spectrum(bands=input_bands)

                            # Run the fitting
                            spec.fit.frame(input_bands, fit_conf=capers_cfg, obj_conf_prefix=obj_cfg_section,
                                           default_conf_prefix=default_cfg_section, progress_output=None)

                            # Show the results the lines
                            spec.plot.spectrum(show_categories=False, rest_frame=False)

                            # # Line redshift calculation
                            # z_df = lime.redshift_calculation(spec.frame, weight_parameter='profile_flux')
                            # z_gaussian, z_gaussian_err = z_df.loc['spec_0', ['z_mean', 'z_std']]
                            # files_sample.loc[idx_obj, ['z_gaussian', 'z_gaussian_err']] = z_gaussian, z_gaussian_err
                            # print(f'- line fitting: z_gaussian = {z_gaussian:0.3f}')
                            #
                            # # Save the results
                            # spec.save_frame(obj_logs_file)
                            # spec.plot.grid(output_address=logs_folder/f'{fits_stem}_grid.png')

                # Save the frame
                lime.save_frame(log_address, files_sample.frame)
                counter += 1

# Save the dataframe with the ASPECT redshifts
lime.save_frame(log_address, files_sample.frame)

# Show the objects with errors
for object, error_msg in error_objects.items():
    print(f'{object} : {error_msg}')

