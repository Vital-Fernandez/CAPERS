from pathlib import Path
import numpy as np
import pandas as pd
import lime

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

# Line fitting data
norm_flux = capers_cfg['data']['norm_flux']
default_cfg_section = capers_cfg['data']['default_prism_cfg_section']
boundary_cfg_section = capers_cfg['data']['boundary_prism_cfg_section']

# Load the sample table (create daily back-up if necessary)
log_address = tables_folder/f'{sample}_{version}_files_log.txt'
create_backup(log_address)
files_sample = lime.Sample(log_address, levels=["sample", "id", "file"], load_function=nirspec_load_function,
                           norm_flux=norm_flux, folder_obs=observations_folder)

# Crop the sample to certain objects and observations types
idcs_selection = files_sample.frame.disp == 'prism' # files_sample.index
files_sample = files_sample[idcs_selection]

# Treatment selection
start_object = 23#34
ASPECT_CHECK = True
REDSHIFT_CHECK = True
BANDS_CHECK = True
MEASURE_CHECK = True
RESET_redshift = False
SAVE_CHECK = False

# Loop through the objects
counter = 0
obj_list = files_sample.index.get_level_values('id').unique()
for i, obj in enumerate(obj_list):

    # Get list of its observations and their extensions
    idcs_obj = files_sample.index.get_level_values('id') == obj
    ext_list = files_sample.frame.loc[idcs_obj].ext.unique()

    # Get the best reduction if available
    ext_target = next((ext for ext in preference_list if ext in ext_list), None)
    if ext_target is not None:

        idcs_obj = idcs_obj & (files_sample.frame.ext == ext_target)
        obj_sample = files_sample[idcs_obj]

        # Loop through the observations
        if  obj_sample.size > 0:
            for j, idx_obj in enumerate(obj_sample.index):
                if counter >= start_object:

                    # Objects ID params
                    MSA = idx_obj[1]
                    disp, pointing = obj_sample.loc[idx_obj, ['disp', 'pointing']]
                    fits_stem = Path(obj_sample.loc[idx_obj, 'file_path']).stem
                    id_message = f'\n{counter}: MSA {MSA} ) Disp = {disp}, pointing = {pointing} '
                    print(id_message)

                    # output object files
                    obj_comps_file = comps_folder / f'{fits_stem}_components.txt'
                    obj_bands_file = bands_folder / f'{fits_stem}_bands.txt'
                    obj_logs_file = logs_folder/ f'{fits_stem}_log.txt'

                    # Aspect detection
                    if ASPECT_CHECK:
                        z_obj = 0
                        spec = files_sample.get_spectrum(idx_obj, redshift=z_obj)
                        spec.features.detection(show_steps=False, exclude_continuum=True)
                        save_detection_results(spec, idx_obj, files_sample, obj_comps_file)

                        # Get Aspect redshift
                        if files_sample.loc[idx_obj, 'n_lines'] > 1:
                            z_brute =  spec.fit.redshift(redshift_bands_df, plot_results=True)
                            files_sample.loc[idx_obj, 'z_aspect_brute'] = z_brute

                            print(f'- aspect check: z_brute = {z_brute}')
                            title = id_message + f', z_brute = {z_brute}'
                            spec.update_redshift(0 if z_brute is None else z_brute)
                            spec.plot.spectrum(show_categories=True, ax_cfg={'title': title}, bands=prism_bands_df)

                    # Focus on positive cases
                    PROCEED_CHECK = files_sample.loc[idx_obj, 'n_lines'] > 1
                    if PROCEED_CHECK:

                        # Redshift review
                        if REDSHIFT_CHECK:

                            # Get previous redshift
                            z_manual, z_aspect = files_sample.loc[idx_obj, ['z_manual','z_aspect_brute']]
                            z_obj = z_manual if pd.notnull(z_manual) else z_aspect
                            print(f'- redshift check: z_manual = {z_obj}')

                            title = f'MSA {obj} ({disp})'
                            output_idcs = files_sample.index.get_level_values('id') == obj
                            files_sample.check.redshift(idcs_obj, REF_LINEs, output_idcs=output_idcs,
                                                        initial_z=z_obj, redshift_column='z_manual',
                                                        title=title, legend_handle='file', maximize=False,
                                                        output_file_log=log_address)

                        if BANDS_CHECK:

                            # If first time, prepare some data
                            if pd.isnull(files_sample.loc[idx_obj, 'z_bands']):
                                files_sample.loc[idx_obj, 'z_bands'] = files_sample.loc[idx_obj, 'z_manual']
                            z_obj = files_sample.loc[idx_obj, 'z_bands']

                            # If there is a redshift measurement
                            if pd.notnull(z_obj):

                                # Load the object
                                print(f'- bands check: z_bands = {z_obj}')
                                spec = files_sample.get_spectrum(idx_obj, redshift=z_obj, detection_file=obj_comps_file)

                                # Prepare bands file if not available
                                if obj_bands_file.is_file():
                                    obj_bands = lime.load_frame(obj_bands_file)
                                else:
                                    obj_bands = spec.retrieve.line_bands(ref_bands=prism_bands_df, components_detection=True)
                                    lime.save_frame(obj_bands_file, obj_bands)

                                spec.check.bands(obj_bands_file, ref_bands=prism_bands_df, maximize=False, n_pixels=50,
                                                 ax_cfg={"title": f'MSA{MSA} ({disp})'})

                                # Save the z_bands to the dataframe
                                lime.save_frame(log_address, files_sample.frame)


                        if MEASURE_CHECK:

                            # If there are bands
                            if obj_bands_file.is_file():

                                # If first time, prepare some data
                                z_obj = files_sample.loc[idx_obj, 'z_bands']
                                spec = files_sample.get_spectrum(idx_obj, redshift=z_obj, detection_file=obj_comps_file)

                                # Recover the line fitting configuration for these objects
                                obj_cfg_section = f'{MSA}_{fits_stem}'
                                obj_cfg = capers_cfg.get(obj_cfg_section)
                                obj_cfg_section = "boundary_prism" if obj_cfg is None else obj_cfg_section
                                print(f'- {obj_cfg_section}_line_fitting')

                                spec.fit.frame(obj_bands_file, fit_conf=capers_cfg, id_conf_prefix=obj_cfg_section,
                                               default_conf_prefix=default_cfg_section, progress_output=None)
                                spec.plot.spectrum(show_categories=False)

                                # Line redshift calculation
                                z_df = lime.redshift_calculation(spec.frame, weight_parameter='profile_flux')
                                z_gaussian, z_gaussian_err = z_df.loc['spec_0', ['z_mean', 'z_std']]
                                files_sample.loc[idx_obj, ['z_gaussian', 'z_gaussian_err']] = z_gaussian, z_gaussian_err
                                print(f'- line fitting: z_gaussian = {z_gaussian:0.3f}')

                                # Save the results
                                lime.save_frame(log_address, files_sample.frame)
                                spec.save_frame(obj_logs_file)
                                spec.plot.grid(output_address=logs_folder/f'{fits_stem}_grid.png')

                counter += 1

