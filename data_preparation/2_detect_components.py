from pathlib import Path
import numpy as np
import lime

from support.tools import create_backup, nirspec_load_function, save_detection_results, read_prism_r_curve, catch_nan_spectra

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
R_interpolator = read_prism_r_curve(data_folder/capers_cfg['file_structure']['prism_R_curve_file'], units_factor=10000)
redshift_bands_df.drop(index=["H1_1216A", "He2_1640A", "H1_4342A"], inplace=True)

# Line fitting data
norm_flux = capers_cfg['data']['norm_flux']
default_cfg_section = capers_cfg['data']['default_prism_cfg_section']
boundary_cfg_section = capers_cfg['data']['boundary_prism_cfg_section']

# Load the sample table (create daily back-up if necessary)
log_address = tables_folder/f'{sample}_files_log.txt'
create_backup(log_address)
files_sample = lime.Sample(log_address, levels=["sample", "id", "file"], load_function=nirspec_load_function,
                           norm_flux=norm_flux, folder_obs=observations_folder)

# Crop the sample to certain objects and observations types
idcs_selection = files_sample.frame.disp == 'prism' # files_sample.index
files_sample = files_sample[idcs_selection]

# Treatment selection
start_object = 0
ASPECT_CHECK = True
REDSHIFT_CHECK = True
BANDS_CHECK = True
MEASURE_CHECK = True
RESET_redshift = False
SAVE_CHECK = False

# Determine the objects to treat
obj_list = files_sample.index.get_level_values('id').unique()
error_objects = {}

# Loop through the objects
counter = 0 # DONT TOUCH
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
                    MPT = idx_obj[1]
                    disp, pointing = obj_sample.loc[idx_obj, ['disp', 'pointing']]
                    fits_stem = Path(obj_sample.loc[idx_obj, 'file_path']).stem
                    id_message = f'\n{counter}: MPT {MPT} ) Disp = {disp}, pointing = {pointing} '
                    print(id_message)

                    # output object files
                    obj_comps_file = comps_folder / f'{fits_stem}_components.txt'
                    obj_bands_file = bands_folder / f'{fits_stem}_bands.txt'
                    obj_logs_file = logs_folder/ f'{fits_stem}_log.txt'

                    # Aspect detection
                    if ASPECT_CHECK:
                        z_obj = 0
                        spec, err = catch_nan_spectra(files_sample, idx_obj, z_obj)

                        if spec is None:
                            print(f'- Error opening: {err}')
                            error_objects[obj] = fits_stem
                            # spec = files_sample.get_spectrum(idx_obj, redshift=z_obj)
                        else:
                            # spec.plot.spectrum()
                            spec.features.detection(show_steps=False, exclude_continuum=True)
                            save_detection_results(spec, idx_obj, files_sample, obj_comps_file)
                            # spec.plot.spectrum(show_categories=True)

                            # Get Aspect redshift
                            if files_sample.loc[idx_obj, 'n_lines'] >= 1:
                                z_phot = obj_sample.loc[idx_obj, "z_UNICORN"]
                                res_power = R_interpolator(spec.wave.data)
                                print(f'- Phot value: z_phot = {z_phot}')
                                z_brute = spec.fit.redshift(redshift_bands_df, res_power=res_power, plot_results=False, sigma_factor=1)
                                z_brute_bound =  spec.fit.redshift(redshift_bands_df, res_power=res_power, plot_results=False, sigma_factor=1,
                                                                   z_min=np.maximum(z_phot-1.5,0), z_max=np.minimum(z_phot+1.5, 12))
                                files_sample.loc[idx_obj, 'z_aspect_brute'] = z_brute
                                files_sample.loc[idx_obj, 'z_aspect_brute_bound'] = z_brute_bound

                                print(f'- aspect check: z_brute = {z_brute}')
                                title = id_message + f', z_brute = {z_brute}'
                                # spec.update_redshift(0 if z_brute is None else z_brute)
                                # spec.plot.spectrum(show_categories=True, ax_cfg={'title': title}, bands=prism_bands_df,
                                #                    # maximize=True)
                                #                    output_address=comps_folder/f'{fits_stem}_components.png')

                counter += 1

# Save the dataframe with the ASPECT redshifts
lime.save_frame(log_address, files_sample.frame)

# Show the objects with errors
for object, error_msg in error_objects.items():
    print(f'{object} : {error_msg}')

