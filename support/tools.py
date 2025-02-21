from pathlib import Path
from datetime import datetime
from astropy.io import fits
import numpy as np
import lime
import aspect
from scipy.interpolate import interp1d
from astropy.visualization import ZScaleInterval
from matplotlib import pyplot as plt
from matplotlib.widgets import SpanSelector
from lime.io import LiMe_Error


Z_FUNC_CMAP = ZScaleInterval()

TARGET_COMPONENTS = np.array(['emission', 'cosmic-ray', 'doublet', 'absorption', 'dead-pixel'])


def catch_nan_spectra(log_df, id_file, redshift):

    try:
        spec = log_df.get_spectrum(id_file, redshift=redshift)
        err_message = None
    except LiMe_Error as e:
        spec = None
        err_message = e

    return spec, err_message

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


def search_spectra_capers(root_folder, ext_list):

    # Change type to path
    root_folder = Path(root_folder)

    # Loop through the files and store teh files you are interested
    file_list = []
    pointing_folders = list(root_folder.iterdir())
    for point_folder in pointing_folders:
        if point_folder.is_dir():
            for ext in ext_list:
                file_list += list(point_folder.glob(f'*{ext}'))

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


def unpack_nirspec_fits(fname):

    with fits.open(fname) as hdu_list:
        header = (hdu_list[0].header, hdu_list[1].header)
        wave_array = np.linspace(header[1]['WAVSTART'], header[1]['WAVEND'], header[1]['NAXIS1'], endpoint=True)
        flux_array = hdu_list[1].data
        err_array = hdu_list[2].data

    return wave_array, flux_array, err_array, header

def nirspec_load_function(log_df, obs_idx, data_folder, **kwargs):

    # Use redshift provided
    if kwargs.get('redshift') is not None:
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
        wave_array, flux_array, err_array, header = unpack_nirspec_fits(file_spec)
        wave_array = wave_array * 10000
        objSpec = wave_array, flux_array, err_array, header

        # with fits.open(file_spec) as hdu_list:
        #     header = (hdu_list[0].header, hdu_list[1].header)
        #     wave_array = np.linspace(header[1]['WAVSTART'], header[1]['WAVEND'], header[1]['NAXIS1'], endpoint=True) * 10000
        #     flux_array = hdu_list[1].data
        #     err_array = hdu_list[2].data

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

        # Add log if present
        print()

    # Unknown
    else:
        raise KeyError(f'Not recognizing the file type: {file_spec.as_posix()}')

    return objSpec


def save_detection_results(spectrum, idx, sample, obj_comps_file):

    # Extract prediction and confidence arrays
    pred_arr, conf_arr = spectrum.features.pred_arr, spectrum.features.conf_arr

    # Getting the indexes
    for shape in TARGET_COMPONENTS:
        idcs_lines = pred_arr == aspect.cfg["shape_number"][shape]

        # Count the number of this features:
        if np.any(idcs_lines):
            edges = np.diff(np.concatenate(([0], idcs_lines, [0])))
            num_true_segm = np.sum(edges == 1)
            mean_conf = np.mean(conf_arr[idcs_lines])
        else:
            num_true_segm, mean_conf = 0, 0

        # Store the results
        sample.loc[idx, f'n_{shape}'] = num_true_segm
        sample.loc[idx, f'confidence_{shape}'] = mean_conf

    # Get indeces with predictions we are interested
    # idcs_lines = (pred_arr == 3) | (pred_arr == 4) | (pred_arr == 7) | (pred_arr == 9) | (pred_arr == 10)

    # # Identify where changes occur (edges of ones and zeros)
    # if np.any(idcs_lines):
    #     edges = np.diff(np.concatenate(([0], idcs_lines, [0])))
    #     num_true_segm = np.sum(edges == 1)
    #     mean_conf = np.mean(conf_arr[idcs_lines])
    #     # start_indices = np.where(edges == 1)[0]
    #     # end_indices = np.where(edges == -1)[0] - 1
    #
    # else:
    #     num_true_segm, mean_conf = 0, 0

    # Save the results in table
    sample.loc[idx, 'n_lines'] = sample.loc[idx, f'n_emission'] + sample.loc[idx, f'n_doublet']
    # sample.loc[idx, 'n_lines'] = num_true_segm
    # sample.loc[idx, 'emission_confidence'] = mean_conf

    # Save to a file
    np.savetxt(obj_comps_file, np.c_[(pred_arr.astype(int), conf_arr.astype(int))], fmt='%i')

    return

def read_prism_r_curve(fname, units_factor=None):

    with fits.open(fname) as hdul:
        data = hdul[1].data

    # Create the interpolation function
    units_factor = units_factor if units_factor is not None else 1
    interp_func = interp1d(data['WAVELENGTH'] * units_factor, data['R'], kind='linear', fill_value="extrapolate")

    return interp_func

def round_off_rating(number):
    return round(number * 2) / 2

def maximize_center_fig(maximize_check=False, center_check=False):

    if maximize_check:

        # Windows maximize
        mng = plt.get_current_fig_manager()

        try:
            mng.window.showMaximized()
        except:
            try:
                mng.resize(*mng.window.maxsize())
            except:
                print()

    if center_check:

        try:
            mngr = plt.get_current_fig_manager()
            mngr.window.setGeometry(1100, 300, mngr.canvas.width(), mngr.canvas.height())
        except:
            print()

    return

class TraceSelection:

    def __init__(self, file_path, mpt, dispenser, database_addres=None):

        # Load the data
        self.wave_array, self.flux_array, err_array, header = unpack_nirspec_fits(file_path)
        self.wave_array = self.wave_array * 10000
        # self.wave_array, self.flux_array, err_array, header = load_nirspec_fits(file_path.as_posix())
        self.y_array = np.arange(self.flux_array.shape[0])
        self.idcs = [0, -1]
        self.mpt, self.disp = mpt, dispenser
        self.title = file_path.name
        self.df_address = database_addres
        self.df = lime.load_frame(database_addres)

        # Generate the figure
        self.fig, (self.axIm, self.axSelec) = plt.subplots(2, figsize=(8, 6))

        # Plot Total image and spec
        title = f'{self.title} {self.df.loc[self.title, "1D_extr_px1"]} {self.df.loc[self.title, "1D_extr_px2"]}'
        self.plot_2D(self.axIm, [0, -1], title=title)
        self.plot_1d(self.axSelec, self.idcs, crop=False)

        # Span selector function
        span = SpanSelector(self.axIm, self.on_select, "vertical", useblit=True, props=dict(alpha=0.5, facecolor="tab:blue"),
                            interactive=True, drag_from_anywhere=True, button=1)

        # Display the image
        maximize_center_fig(True)
        plt.show()

        return

    def plot_2D(self, ax,  idcs, title=None):

        z1, z2 = Z_FUNC_CMAP.get_limits(self.flux_array[:, idcs[0]:idcs[-1]])

        disp_low, disp_high = self.wave_array[idcs[0]], self.wave_array[idcs[-1]-1]
        spa_low, spa_high = 0, self.flux_array.shape[0]

        extend = np.array([disp_low, disp_high, spa_high, spa_low])

        im = ax.imshow(self.flux_array, cmap='gist_heat', vmin=z1, vmax=z2, aspect=2, origin='lower', interpolation='none')

        if title is not None:
            ax.set_title(title)

        return

    def plot_1d(self, ax, idcs, color='tab:blue', crop=True, title=None):

        flux_sum = np.sum(self.flux_array[idcs[0]:idcs[-1], :], axis=0)
        ax.step(self.wave_array, flux_sum, color=color)

        if crop:
            ax.set_xlim((self.wave_array[idcs[0]:idcs[-1]][0], self.wave_array[idcs[0]:idcs[-1]][-1]))

        if title is not None:
            ax.set_title(title)

        return

    def on_select(self, ymin, ymax):

        self.save_address(ymin, ymax)
        self.idcs = np.searchsorted(self.y_array, (ymin, ymax))

        if len(self.idcs) >= 2:

            # 1d selection spectrum
            self.axSelec.clear()
            self.plot_1d(self.axSelec, self.idcs, crop=False)

            # 2d selection
            title = f'{self.title} {self.df.loc[self.title, "1D_extr_px1"]} {self.df.loc[self.title, "1D_extr_px2"]}'
            self.axIm.clear()
            self.plot_2D(self.axIm, [0, -1], title)

            self.fig.canvas.draw_idle()

        return

    def save_address(self, ymin, ymax):

        if self.df_address is not None:
            min_round, max_round = round_off_rating(ymin), round_off_rating(ymax)
            self.df.loc[self.title, '1D_extr_px1':'1D_extr_px2'] = min_round, max_round
            lime.save_frame(self.df_address, self.df)

        return