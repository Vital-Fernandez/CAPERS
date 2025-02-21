import lime
import numpy as np
from pathlib import Path
from astropy.io import fits
from support.tools import nirspec_load_function

# State data location
sample = 'CAPERS_UDS_V0.1'
tables_folder ='/home/vital/Dropbox/Astrophysics/Data/CAPERS/tables'
observations_folder = f'/home/vital/Astrodata/CAPERS'
flux_log_address = f'{tables_folder}/{sample}_flux_log.txt'

# Create observations sample variable
capers_sample = lime.Sample(flux_log_address, levels=["sample", "id", "file", "line"], norm_flux=1e-22,
                            load_function=nirspec_load_function, folder_obs=observations_folder)

# Index target object
idcs_object = capers_sample.frame.MPT_number == 22431
first_idx = capers_sample.frame.loc[idcs_object].index[0]
# Extract spectrum
spec = capers_sample.get_spectrum(first_idx, z_column='z_bands')

# Load its line measurements
lines_df = capers_sample.frame.loc[idcs_object].droplevel(['sample', 'id', 'file'])
spec.load_frame(lines_df)
print(spec.frame.z_line)

# Plot the data
spec.plot.spectrum(rest_frame=True)
spec.plot.grid()

