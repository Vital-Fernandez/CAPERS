from pathlib import Path
import numpy as np
import pandas as pd
import lime
from support.tools import TraceSelection, nirspec_load_function, create_backup

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

# Line fitting data
default_cfg_section = capers_cfg['data']['default_prism_cfg_section']
boundary_cfg_section = capers_cfg['data']['boundary_prism_cfg_section']

# Load the sample table (create daily back-up if necessary)
log_address = tables_folder/f'{sample}_{version}_files_log.txt'
create_backup(log_address)

files_sample = lime.Sample(log_address, levels=["sample", "id", "file"], load_function=nirspec_load_function,
                           norm_flux=capers_cfg['data']['norm_flux'], folder_obs=observations_folder)

# Crop the sample to certain objects and observations types
idcs_selection = (files_sample.frame.disp == 'prism') & (files_sample.frame.ext == 's2d') # files_sample.index
files_sample = files_sample[idcs_selection]
file_list = files_sample.index.get_level_values('file').to_numpy()

# Generate table with the files
trace_file = tables_folder/Path(capers_cfg['file_structure']['trace_file'])
if not trace_file.is_file():
    trace_df = pd.DataFrame(index=file_list, columns=('1D_extr_px1', '1D_extr_px2'))
    lime.save_frame(trace_file, trace_df)

# Start object
start_object = 0

# Loop throught the files
section_dict, spec_counter = {}, 0
for i, fits2d in enumerate(file_list):

    # Identifiers
    idx = files_sample.index.get_level_values('file') == fits2d
    MPT, disp, pointing, ext, fits_path = files_sample.loc[idx, ['MPT_number', 'disp', 'pointing', 'ext', 'file_path']].to_numpy()[0]
    print(f'\nSpec {i}) MPT {MPT}, file {fits2d}')

    if i >= start_object:

        # Interactive show
        # section_dict[MPT] = []
        aSelector = TraceSelection(observations_folder/fits_path, MPT, disp, trace_file)

        # # Run DS9 command to open the spectrum
        # ds9_command = f'start "ds9" "{ds9_path}" "{fits2d.as_posix()}" "-mode" "none" "-geometry" "1200x800" "-zscale" "-cmap" "Heat" "-zoom" "to" "fit"'
        # print(f'{file_name}\n{ds9_command}\n')
        # p = subprocess.Popen(ds9_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # out = p.communicate()[0]

    spec_counter += 1

