from pathlib import Path
import numpy as np
import pandas as pd
import lime

from support.tools import (create_backup, capers_load_function, save_detection_results, read_prism_r_curve,
                           z_selection)

lime.theme.set_style('dark')

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
map_bands_vsigma = capers_cfg['data']['map_bands_vsigma']
default_cfg_section = capers_cfg['default_prism_line_fitting']

# Load the sample table (create daily back-up if necessary)
log_address = tables_folder/f'{sample}_files_log.txt'
create_backup(log_address)
files_sample = lime.Sample(log_address, levels=["sample", "id", "pointing"], load_function=capers_load_function,
                           norm_flux=norm_flux, folder_obs=observations_folder)

# Treatment selection
start_object =  0#33# 38 36 33 32 7 3
ASPECT_CHECK = False
ASPECT_REDSHIFT = False
REDSHIFT_CHECK = False
BANDS_CHECK = True
MEASURE_CHECK = False
RESULTS_CHECK = False

# Determine the objects to treat
tier_df = files_sample.frame.sort_values(by='z_UNICORN', ascending=[False])
idcs = tier_df.MPT_number == 32172
obj_list = tier_df.loc[idcs].index.get_level_values('id').unique()

# idcs = (tier_df.z_manual > 4) & (tier_df.z_manual < 6)
# obj_list = tier_df.loc[idcs].index.get_level_values('id').unique()
# obj_list = tier_df.index.get_level_values('id').unique()

# Loop through the objects
counter = 0 # DONT TOUCH
saving_iter_n = 20
failing_opening = {}
error_objects = {}
for i, obj in enumerate(obj_list):

    # Get list of its observations and their extensions
    idcs_obj_group = files_sample.loc[files_sample.index.get_level_values('id') == obj].index

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
                obj_cfg_section = f'{MPT}_{fits_stem}'
                obj_cfg = capers_cfg.get(f'{obj_cfg_section}_line_fitting')
                obj_cfg_section = "boundary_prism" if obj_cfg is None else obj_cfg_section

                composite_key = 'line_composites_high_z' if z_obj > 5 else 'line_composites_medium_z'
                composite_lines = capers_cfg['lines_data'][composite_key]


                # Load spectrum
                spec = files_sample.get_spectrum(idx_obj, redshift=z_obj, detection_file=obj_comps_file,
                                                 file_path=file_path)

                # Check for full nan spectra
                if spec is None: failing_opening[MPT] = file_path; continue

                # Prompt
                print(f'\n{counter}: MPT {MPT}) Disp = {disp}, pointing = {pointing}, z_obj={z_obj} (z_tier={z_tier})')
                if RESULTS_CHECK:
                    print(f'- Results')

                    specVis = files_sample.get_spectrum(idx_obj, redshift=z_obj, detection_file=obj_comps_file,
                                                     file_path=file_path, lines_log_path=obj_logs_file)
                    specVis.plot.spectrum(maximize=True, rest_frame=True)

                # Aspect detection
                if ASPECT_CHECK:
                    spec.infer.components(show_steps=False, exclude_continuum=True)
                    save_detection_results(spec, idx_obj, files_sample, obj_comps_file)

                # Redshift from lines intervals (n_lines > 1)
                if ASPECT_REDSHIFT:
                    if files_sample.loc[idx_obj, 'n_lines'] >= 1:

                        # Compute the redshifts
                        limits = (np.maximum(z_phot-2,0), np.minimum(z_phot+2, 16)) if z_phot < 3 else (None, None)
                        z_key = spec.fit.redshift(redshift_bands,  sigma_factor=1, z_min=limits[0], z_max=limits[1],
                                                  mode='key', plot_results=False,)
                        z_xor = spec.fit.redshift(redshift_bands, sigma_factor=1, z_min=limits[0], z_max=limits[1],
                                                  mode='xor', plot_results=False)

                        # Check redshift compatibility between techniques
                        z_diff_check, z_tier_aspect = False, 1
                        if z_key is not None and z_xor is not None:
                            percent_diff = np.abs(z_key / z_xor - 1) * 100
                            z_diff_check = percent_diff <= 5
                            z_tier_aspect = 2 if z_diff_check else 1

                        # Save the measurements:
                        files_sample.loc[idx_obj, 'z_aspect_key'] = z_key
                        files_sample.loc[idx_obj, 'z_aspect_xor'] = z_xor
                        if z_tier < z_tier_aspect:
                            files_sample.loc[idx_obj, 'z_tier'] = z_tier_aspect

                        # Prompt
                        print(f'- Aspect Check: z_phot = {z_phot}')
                        print(f'-- Key: z_key = {z_key}')
                        print(f'-- Key: z_xor = {z_xor}')
                        print(f'--- Compatible: {z_diff_check}')

                # Redshift review
                if REDSHIFT_CHECK:

                    z_tier, z_obj = z_selection(files_sample, idx_obj)
                    print(f'- Manual check: z_manual = {z_obj}')

                    title = f'MPT {obj} ({disp}) z_photo={z_phot:0.3f}; '
                    output_idcs = files_sample.index.get_level_values('id') == obj
                    idcs__group = files_sample.index.get_level_values('id') == obj

                    files_sample.check.redshift(idx_obj, lines_visualize, output_idcs=output_idcs,
                                                initial_z=z_obj, redshift_column='z_manual',
                                                title=title, legend_handle='optext', maximize=True,
                                                output_file_log=log_address, file_path=file_path)

                    # Establish the tier
                    if pd.notnull(files_sample.loc[idx_obj, 'z_manual']):
                        files_sample.loc[idx_obj, 'z_tier'] = 3

                # Review the bands
                if BANDS_CHECK:

                    z_tier, z_obj = z_selection(files_sample, idx_obj)
                    if pd.notnull(z_obj):

                        # Load the object
                        print(f'- Bands check: z_bands = {z_obj}')
                        if spec.redshift != z_obj:
                            spec.update_redshift(z_obj)


                        bands = spec.retrieve.line_bands(ref_bands=prism_bands_df, fit_cfg=capers_cfg,
                                                         band_vsigma=200, n_sigma=4, instrumental_correction=True,
                                                         default_cfg_prefix='default_prism', automatic_grouping=True,
                                                         Rayleigh_threshold=3, map_band_vsigma=map_bands_vsigma)

                        spec.plot.spectrum(bands=bands, rest_frame=True)

                        # if bands.index.size > 0:
                        #     spec.plot.spectrum(bands=bands, rest_frame=True)

                        # Generate bands using tradiotional detection
                        # if not obj_bands_file:
                        #     input_bands = spec.retrieve.line_bands(ref_bands=prism_bands_df, fit_cfg=capers_cfg,
                        #                                            composite_lines=composite_lines,
                        #                                            default_cfg_prefix='default_prism',
                        #                                            obj_cfg_prefix=obj_cfg_section, band_vsigma=300,
                        #                                            vacuum_waves=True,  components_detection=False)
                        #     spec.fit.continuum(degree_list=[3, 5, 6, 7], emis_threshold=[3, 2, 2, 1.5], plot_steps=False)
                        #     input_bands = spec.infer.peaks_troughs(input_bands, emission_type=True, sigma_threshold=3,
                        #                                            plot_steps=False, log_scale=False, maximize=True)
                        #     lime.save_frame(obj_bands_file, input_bands)
                        # spec.plot.spectrum(rest_frame=True, bands=input_bands)

                        # # Check the bands
                        # print(obj_bands_file)
                        # # spec.plot.spectrum()
                        # spec.check.bands(obj_bands_file, ref_bands=prism_bands_df, fit_cfg=capers_cfg,
                        #                  band_vsigma=200, n_sigma=4, instrumental_correction=True,
                        #                  default_cfg_prefix='default_prism', automatic_grouping=True,
                        #                  Rayleigh_threshold=3, map_band_vsigma=map_bands_vsigma,
                        #                  vacuum_waves=True, exclude_continua=True, maximize=True)
                        #
                        #
                        #
                        # # spec.check.bands(obj_bands_file, ref_bands=prism_bands_df, band_vsigma=300,
                        # #                  fit_cfg=capers_cfg, default_cfg_prefix='default_prism',
                        # #                  obj_cfg_prefix=obj_cfg_section, vacuum_waves=True,
                        # #                  exclude_continua=True, composite_lines=composite_lines, maximize=True)
                        #
                        # # Save the z_bands to the dataframe
                        # if pd.isnull(files_sample.loc[idx_obj, 'z_bands']) and obj_bands_file.is_file():
                        #     files_sample.loc[idx_obj, 'z_bands'] = z_obj


            # Save the frame every 20 iterations
            if ((i + 1) % saving_iter_n == 0) and (counter >= start_object):
                lime.save_frame(log_address, files_sample.frame)
                print(f'\n- FILES LOG SAVED ({i})\n')

        counter += 1

# # Save the dataframe with the ASPECT redshifts
# lime.save_frame(log_address, files_sample.frame)
#
