from pathlib import Path
from datetime import datetime
from astropy.io import fits
import numpy as np
import lime
import aspect


def search_spectra_ceers(root_folder, ext_list):

    # Change type to path
    root_folder = Path(root_folder)

    # Loop through the files and store teh files you are interested
    file_list = []
    pointing_folders = list(root_folder.iterdir())
    for point_folder in pointing_folders:
        if point_folder.is_dir():
            disp_folder_list = list(point_folder.iterdir())
            for disp_folder in disp_folder_list:
                if disp_folder.is_dir():
                    for ext in ext_list:
                        file_list += list(disp_folder.glob(f'*{ext}'))

    return file_list


def review_masked_files(file_list):

    new_file_list = []
    for spec_address in file_list:

        file_str = spec_address.as_posix()
        if file_str.endswith('s2d.fits'):
            new_file_list.append(spec_address)

        elif file_str.endswith('x1d-masked.fits'):
            new_file_list.append(spec_address)

        else:
            masked_name = Path(file_str.replace('x1d.fits', 'x1d-masked.fits'))
            if masked_name not in file_list:
                new_file_list.append(spec_address)
                print(f'Missing object: {spec_address.name}')

    return new_file_list


def create_backup(file_path):

    """
    Create a backup of a file with the current date in the file name,
    only if a backup for the current day doesn't already exist.

    Parameters:
    - file_path (str or Path): The path to the original file.
    """

    file_path = Path(file_path)

    if not file_path.exists():
        print(f"Error: The file {file_path} does not exist.")
        return

    # Get today's date in DD-MM-YYYY format
    today = datetime.now().strftime("%Y-%m-%d")

    # Construct the backup file name
    backup_name = f"{file_path.stem}_backup_{today}{file_path.suffix}"
    backup_path = file_path.parent / backup_name

    # Check if today's backup already exists
    if backup_path.exists():
        return

    # Create the backup file
    else:
        with open(file_path, 'rb') as original_file:
            with open(backup_path, 'wb') as backup_file:
                backup_file.write(original_file.read())
        print(f"Backup created: {backup_path}")

    return


def nirspec_load_function(log_df, obs_idx, data_folder, **kwargs):

    # Use redshift provided
    if kwargs['redshift'] is not None:
        z_obj = kwargs['redshift']

    # Else use the log redshift from input column
    else:
        z_column = kwargs['z_column']

        z_obj = log_df.loc[obs_idx, z_column]
        z_obj = None if np.isnan(z_obj) else z_obj

    norm_flux = kwargs['norm_flux']

    # Interpolator for the R
    R_interpolator = kwargs.get('R_interpolator')

    # Get file address
    file_spec = Path(data_folder)/log_df.loc[obs_idx, 'file_path']

    # 2d files
    if "s2d" in file_spec.as_posix():
        with fits.open(file_spec) as hdu_list:
            header = (hdu_list[0].header, hdu_list[1].header)
            wave_array = np.linspace(header[1]['WAVSTART'], header[1]['WAVEND'], header[1]['NAXIS1'], endpoint=True) * 10000
            flux_array = hdu_list[1].data
            err_array = hdu_list[2].data

        objSpec = wave_array, flux_array, err_array, header

    # 1d files
    elif "x1d" in file_spec.as_posix():
        objSpec = lime.Spectrum.from_file(file_spec, instrument='nirspec', norm_flux=norm_flux, redshift=z_obj)
        objSpec.unit_conversion(wave_units_out='Angstrom', flux_units_out='FLAM', norm_flux=norm_flux)

        if R_interpolator is not None:
            objSpec.inst_FWHM = objSpec.wave.data/R_interpolator(objSpec.wave.data)

        with fits.open(file_spec) as hdu_list:
            header = hdu_list[1].header

        objSpec.header = header

        # Include detection arrays if available
        detection_file = kwargs.get('detection_file')
        if detection_file is not None:
            detection_file = Path(detection_file)
            if detection_file.is_file():
                pred_arr, conf_arr = np.loadtxt(detection_file, dtype=int, unpack=True)
                objSpec.features.pred_arr, objSpec.features.conf_arr = pred_arr, conf_arr

    # Unknown
    else:
        raise KeyError(f'Not recognizing the file type: {file_spec.as_posix()}')

    return objSpec


def save_detection_results(spectrum, idx, sample, obj_comps_file):

    # Extract prediction and confidence arrays
    pred_arr, conf_arr = spectrum.features.pred_arr, spectrum.features.conf_arr

    # Get indeces with predictions we are interested
    idcs_lines = (pred_arr == 3) | (pred_arr == 4)| (pred_arr == 7) | (pred_arr == 9) | (pred_arr == 10)

    # Identify where changes occur (edges of ones and zeros)
    if np.any(idcs_lines):
        edges = np.diff(np.concatenate(([0], idcs_lines, [0])))
        num_true_segm = np.sum(edges == 1)
        mean_conf = np.mean(conf_arr[idcs_lines])
        # start_indices = np.where(edges == 1)[0]
        # end_indices = np.where(edges == -1)[0] - 1

    else:
        num_true_segm, mean_conf = 0, 0

    # Save the results in table
    sample.loc[idx, 'n_lines'] = num_true_segm
    sample.loc[idx, 'emission_confidence'] = mean_conf

    # Save to a file
    np.savetxt(obj_comps_file, np.c_[(pred_arr.astype(int), conf_arr.astype(int))], fmt='%i')

    return


def read_prism_r_curve(fname):


    with fits.open(fname) as hdul:
        data = hdul[1].data

    wave_arr, R_arr, dlds_arr = data['WAVELENGTH'], data['R'], data['DLDS']

    import numpy as np
    from scipy.interpolate import interp1d
    import matplotlib.pyplot as plt

    # Convert to angstroms
    X = wave_arr * 10000
    Y = R_arr  # Dependent variable (e.g., Y = X^2)

    # Create the interpolation function
    interp_func = interp1d(X, Y, kind='linear')  # 'linear', 'quadratic', 'cubic', etc.

    # Interpolated values
    Y_new = interp_func(X)  # Interpolated Y values


    # '''
    # In the dispersion curve files from the NIRSpec instrument, the column labeled 'dlds' represents the dispersion value at
    #  each wavelength, expressed in micrometers per pixel (μm/pixel). This indicates how much the wavelength changes per pixel
    #   on the detector at a given point in the spectrum. These files typically contain three columns: wavelength (μm), dispersion
    #    (μm/pixel), and resolution (λ/Δλ, unitless).
    # '''

    # fig, ax = plt.subplots()
    # ax2 = ax.twinx()
    # disper_line = ax.plot(data['WAVELENGTH'], data['DLDS'], label='Dispersion', color='tab:green')
    # R_line = ax2.plot(data['WAVELENGTH'], data['R'], label='R')
    # # R_line += ax2.plot(wave_arr, R_arr, label='LiMe R')
    # R_line += ax2.plot(X/10000, Y_new, label='Interpolated', color='purple', linestyle=':')
    # ax.update({'xlabel': r'wavelength $(\mu m)$', 'ylabel': r'Dispersion $(\mu m / pixel)$',
    #            'title': 'NIRSPEC Prism dispersion curve'})
    # ax2.update({'ylabel': r'$R = \frac{\lambda}{\Delta \lambda}$'})
    # lines = disper_line + R_line
    # labels = [item.get_label() for item in lines]
    # ax.legend(lines, labels, loc="upper left")
    # plt.tight_layout()
    # plt.show()

    return interp_func