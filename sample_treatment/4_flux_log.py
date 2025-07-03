from pathlib import Path
import numpy as np
import pandas as pd
import lime

from support.tools import create_backup, capers_load_function, z_selection

# Load configuration
cfg_file = '../CAPERS_v3.toml'
capers_cfg = lime.load_cfg(cfg_file)

# Get sample information
version = capers_cfg['meta']['version']
sample = capers_cfg['meta']['sample']
pref_ext_list = capers_cfg['file_structure']['ext_1d_fits']
sample_list = capers_cfg['meta']['core_sample_list']

observations_folder = Path(capers_cfg['file_structure']['observations_folder'])
data_folder = Path(capers_cfg['file_structure']['data_folder'])
source_folder = Path(capers_cfg['file_structure']['source_folder'])

df_dict = {}
for sample in sample_list:
    print(sample)

    logs_folder = data_folder / sample / 'line_logs'
    log_file_list = list(logs_folder.glob('*.txt'))
    fits_file_list = []
    id_list = []

    for i, log_file in enumerate(log_file_list):
        fits_name = log_file.with_name(log_file.name.replace('_log.txt', '.fits')).name

        parts_split = fits_name.split('_')
        idx_id = [i for i, s in enumerate(parts_split) if s.startswith('s0')]
        assert len(idx_id) == 1, f'File {log_file} does not match a expected Id location'
        id_number = int(parts_split[idx_id[0]][1:])

        fits_file_list.append(fits_name)
        id_list.append(id_number)

    flux_sample = lime.Sample.from_file(id_list=id_list, log_list=log_file_list,
                                        file_list=fits_file_list, load_function=None)

    # Re-index with sample and sort
    flux_sample.frame['sample'] = sample
    flux_sample.frame['id'] = flux_sample.index.get_level_values('id')
    flux_sample.frame['file'] = flux_sample.index.get_level_values('file')
    flux_sample.frame['line'] = flux_sample.index.get_level_values('line')
    flux_sample.frame.set_index(['sample', 'id', 'file', 'line'], inplace=True)
    flux_sample.frame['MPT_number'] = np.nan
    flux_sample.frame['z_bands'] = np.nan
    flux_sample.frame['file_path'] = ''

    # Store the df
    df_dict[sample] = flux_sample.frame

# Combined the dataframes
df_combined = pd.concat(list(df_dict.values()), ignore_index=False)

# Sort the dataframe
# df_combined.sort_values(by=[('sample', 'id'), 'wavelength'])

# Save the combined frame
flux_log_address = source_folder/f'CAPERS_sample_flux_log.csv'
lime.save_frame(flux_log_address, df_combined)

# # Load the sample table (create daily back-up if necessary)
# norm_flux = capers_cfg['data']['norm_flux']
# tables_folder = data_folder/sample/'tables'
# log_address = tables_folder/f'{sample}_files_log.txt'
# files_sample = lime.Sample(log_address, levels=["sample", "id", "pointing"], load_function=capers_load_function,
#                            norm_flux=norm_flux, folder_obs=observations_folder)
#
# # Output directories
# bands_folder = data_folder/sample/'bands'
# logs_folder = data_folder/sample/'line_logs'
# bands_file_arr = np.array(list(bands_folder.glob('*')), dtype=str)  # all files and folders in folder_path)
# logs_file_arr = np.array(list(logs_folder.glob('*')), dtype=str)
#
# tier_df = files_sample.frame.sort_values(by='MPT_number', ascending=[False])
# idcs = tier_df['MPT_number'].isin([2648, 3412, 3907, 4241, 4434, 4979, 6754, 6818, 6911, 6916, 7082, 7468, 7720, 8032, 8322, 8790, 8873, 9262, 9575, 9801, 9828, 9914, 10294, 10961, 11148, 12016, 14731, 16182, 16190, 16632, 17082, 17336, 17405, 17740, 17896, 18458, 18479, 19021, 19540, 21007, 21466, 22876, 24135, 24147, 24252, 25053, 26571, 26863, 27175, 27242, 27439, 27640, 27943, 27954, 27996, 28353, 29015, 29358, 31028, 31574, 31886, 32200, 32872, 33201, 33352, 33525, 33644, 33647, 34842, 35566, 36108, 36186, 36842, 36865, 37166, 37315, 37510, 37545, 38098, 38169, 38174, 38431, 38596, 38904, 38918, 38954, 38985, 39574, 39600, 40106, 40108, 40487, 40968, 41852, 42068, 42393, 42681, 42709, 42725, 42759, 42912, 42939, 43095, 43654, 44402, 44706, 44863, 44985, 45307, 45418, 45422, 45457, 45994, 46549, 46621, 46843, 47127, 48083, 49540, 51844, 53872, 68595, 70650, 71321, 72554, 73110, 73866, 74783, 78183, 79909, 83919, 84886, 86508, 86965, 87352, 88196, 88286, 88305, 89329, 89661, 89900, 89926, 90105, 90194, 90786, 91244, 91591, 92493, 92922, 93202, 93655, 93916, 94628, 94786, 95044, 95148, 95308, 95492, 95664, 96161, 96433, 96865, 97724, 97991, 98539, 100091, 100828, 100871, 101272, 101293, 101557, 102478, 102896, 102897, 103248, 103595, 104109, 104781, 105133, 106538, 106975, 107308, 108374, 108834, 108864, 109364, 109483, 109564, 109571, 109798, 110830, 110865, 110889, 111170, 111746, 111999, 112629, 112629, 113091, 113180, 113205, 114183, 114635, 115388, 116303, 117925, 118011, 119331, 124783, 129819, 131478, 131579, 132460, 132792, 134215, 134374, 136567, 138241, 138989, 139386, 139396, 140111, 140857, 140909, 141282, 143397, 144890, 146833, 147449, 147467, 147470, 147518, 147745, 148561, 149010, 149850, 150920, 152047, 152358, 152880, 153767, 153951, 154491, 155468, 155993, 156509, 156537, 156656, 156690, 156968, 159413, 160511, 161875, 162734, 164583, 166377, 166560, 167089, 167645, 167916, 168273, 168643, 168649, 169799, 170062, 170108, 170442, 170984, 171098, 171407, 171663, 171764, 172440, 173276, 173350, 175517, 175527, 175729, 176186, 176879, 176974, 177200, 177212, 177884, 177967, 179148, 179470, 179847, 179898, 181398, 181713, 181896, 182082, 189552])
# z_tiers = tier_df.loc[idcs, 'z_tier'].unique()
#
# # matches = np.char.find(file_paths, keyword) >= 0
#
# for i, idx in enumerate(tier_df.loc[idcs].index):
#     if tier_df.loc[idx, 'z_tier'] < 3:
#         obj_name = tier_df.loc[idx].name[1]
#         match_bands_files = np.char.find(bands_file_arr, obj_name) >= 0
#         match_log_files = np.char.find(logs_file_arr, obj_name) >= 0
#         if (match_bands_files.any() == False) and (match_log_files.any() == True):
#             for idx_file in match_log_files:
#                 file_path = logs_file_arr[match_log_files]




