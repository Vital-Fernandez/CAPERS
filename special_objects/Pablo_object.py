from pathlib import Path
import numpy as np
import pandas as pd
import lime

from support.tools import (create_backup, capers_load_function, save_detection_results, read_prism_r_curve,
                           z_selection)

# lime.theme.set_style('dark')
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


fits_address = '/home/vital/Downloads/CAPERS-COSMOS_p6_s000153965_x1d_optext_v0.3.fits'
bands_address = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/sample/CAPERS_COSMOS_V0.2.1/bands/CAPERS_COSMOS_P6_s000153965_x1d_optext_bands.txt'
z_value = 8.936

spec = lime.Spectrum.from_file(fits_address, 'nirspec', redshift=z_value)
spec.unit_conversion('AA', 'FLAM')

# spec.plot.spectrum(bands=bands_address)
#
spec.check.bands(bands_address, selected_by_default=True,
                 fit_cfg=capers_cfg, exclude_continua=True, maximize=True,
                 **capers_cfg['bands_generation_parameters'])


obj_cfg_section = 's000153965_CAPERS_COSMOS_P6_s000153965_x1d_optext'
spec.fit.frame(bands_address, fit_cfg=capers_cfg, obj_cfg_prefix=obj_cfg_section,
               default_cfg_prefix='default_prism', cont_from_bands=False, err_from_bands=False)
print()
print(spec.frame.z_line)
spec.save_frame('/home/vital/Downloads/s000153965_v0.3_line_measurements.txt')
spec.plot.spectrum(output_address='/home/vital/Downloads/s000153965_v0.3_spectrum.png')
spec.plot.grid(output_address='/home/vital/Downloads/s000153965_v0.3_grid.png')
