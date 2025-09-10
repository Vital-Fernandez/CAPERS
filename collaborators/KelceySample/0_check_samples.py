import pandas as pd
import numpy as np

rubies_sample_file = '/home/vital/Astrodata/KelceySample/rubies_catalogue.csv'
prism_ids_file = '/home/vital/Astrodata/KelceySample/rubies_prism_id.txt'
gratings_ids_file = '/home/vital/Astrodata/KelceySample/rubies_gratings_id.txt'

sample_df = pd.read_csv(rubies_sample_file, sep=',', header=0)
sample_df = sample_df.iloc[:, 2:]
sample_df.set_index(['MPT ID'], inplace=True)
sample_df.sort_values(by=['MPT ID'], inplace=True)

prism_ids = np.loadtxt(prism_ids_file).astype(int)
grat_ids = np.loadtxt(gratings_ids_file).astype(int)

sample_dict = {'prism': {'obj_arr': prism_ids}}

for sample_cfg in ['prism']:

    cfg = sample_dict[sample_cfg]

    for i, obj in enumerate(cfg['obj_arr']):
        fpath = sample_df.loc[obj]
        print(f'{i}) {obj}')