import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lime

s_fname = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/source/group1_egs_final.csv'
s_df = pd.read_csv(s_fname)
s_df.rename(columns={'Galaxy': 'MPT_number'}, inplace=True)

v_fname = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/source/CAPERS_fields_LineRedshifts_v0.csv'
v_df = pd.read_csv(v_fname, index_col=0)

target_sample = 'CAPERS_EGS_V0.2.2'
v_df = v_df.loc[v_df['sample'] == target_sample]

columns_s = ['MPT_number', 'Assigned_Redshift', 'Assigned_Flag', 'photo_z']
merged = pd.merge(v_df, s_df[columns_s], on="MPT_number")

idcs = pd.notnull(merged['Assigned_Redshift']) & pd.notnull(merged['z_centroid'])
x_arr, y_arr = merged.loc[idcs, 'z_centroid'], merged.loc[idcs, 'Assigned_Redshift']
flag_arr = merged.loc[idcs, 'Assigned_Flag']  # use as color scale

# Scatter + 1:1 line
lo = 0
hi = 11#max(x_arr.max(), y_arr.max())

fig, ax = plt.subplots(figsize=(7, 7), dpi=300)

scatter = ax.scatter(x_arr, y_arr, c=flag_arr, cmap="tab10")

ax.plot([lo, hi], [lo, hi], linestyle="--")  # 1-to-1 line

cbar = fig.colorbar(scatter, ax=ax)
cbar.set_label("Assigned_Flag")

ax.set_xlabel("LiMe redshift", fontsize=14)
ax.set_ylabel("Sara's assigned redshift", fontsize=14)
ax.set_title(f"CAPERs EGS redshift comparison ({idcs.sum()}/{len(s_df.index)} galaxies)", fontsize=16)
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
plt.tight_layout()
# plt.show()

# Save the results
plt.savefig(f'CAPERS_{target_sample}_redshift_comparison.png')
lime.save_frame(f'CAPERS_{target_sample}_redshift_comparison.csv', merged)
lime.save_frame(f'CAPERS_{target_sample}_redshift_comparison.txt', merged)