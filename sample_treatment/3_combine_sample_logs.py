from pathlib import Path
import numpy as np
import pandas as pd
import lime

from support.tools import create_backup, capers_load_function, save_detection_results, read_prism_r_curve, z_selection

lime.theme.set_style('dark')

# Load configuration
cfg_file = '../CAPERS_v3.toml'
capers_cfg = lime.load_cfg(cfg_file)

sample_list = ['CAPERS_EGS_V0.2.1', 'CAPERS_COSMOS_V0.2', 'CAPERS_UDS_V0.1']

df_list, sum_rows = [], 0
for sample in sample_list:
    data_folder = Path(capers_cfg['file_structure']['data_folder'])
    tables_folder = data_folder/sample/'tables'
    log_address = tables_folder/f'{sample}_files_log.txt'
    db = lime.load_frame(log_address, levels=["sample", "id", "pointing"])
    print(sample, db.shape)
    sum_rows += db.shape[0]
    df_list.append(db)

# Join the logs
df_combined = pd.concat(df_list, ignore_index=False)

print('Combined', df_combined.shape)
print('Sum total', sum_rows)
print('All unique', df_combined.index.is_unique)
df_combined.sort_values(by=['MPT_number'], inplace=True)

fname = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/source/CAPERS_sample_file_log.txt'
lime.save_frame(fname, df_combined)

fname = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/source/CAPERS_sample_file_log.csv'
lime.save_frame(fname, df_combined)
