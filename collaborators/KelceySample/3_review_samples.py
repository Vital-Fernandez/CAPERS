import pandas as pd
import numpy as np
import lime
import matplotlib.pyplot as plt

def error_division(x, sx, y, sy):

    """
    Propagate uncertainty for z = x / y.

    Parameters
    ----------
    x : float
        Value of numerator
    sx : float
        Uncertainty in numerator
    y : float
        Value of denominator
    sy : float
        Uncertainty in denominator

    Returns
    -------
    z : float
        Result of division
    sz : float
        Propagated uncertainty
    """
    z = x / y
    sz = np.abs(z) * np.sqrt((sx / x) ** 2 + (sy / y) ** 2)
    return z, sz

k_fname = '/home/vital/Downloads/master_data.csv'
v_fname = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/source/CAPERS_FIELDS_flux_log.csv'

k_df = pd.read_csv(k_fname, index_col=0)
k_df = k_df.drop('level_0', axis=1)
v_df = lime.load_frame(v_fname, levels=['sample', 'id', 'file', 'line'])
np.sort(v_df.index.get_level_values('line').unique())

np.sum(v_df.index.get_level_values('line') == 'H1_6565A')

line_k, line_v = 'Haflx', 'H1_6563A_m'
# line_k, line_v = 'Hbflx', 'H1_4861A'
idcs_k = (k_df.Type == 'CAPERS PRISM') & pd.notnull(k_df[line_k])

idcs_v = v_df.index.get_level_values('id').isin(k_df.loc[idcs_k, 'MPTID'].unique().astype(int))
idcs_v = idcs_v & (v_df.index.get_level_values('line') == line_v)
obj_list = np.sort(v_df.loc[idcs_v].index.get_level_values('id'))

# Store the data:
out_df = pd.DataFrame(columns=['Hbeta_k', 'Hbeta_k_err', 'Hbeta_v', 'Hbeta_v_err',
                               'Halpha_k', 'Halpha_k_err', 'Halpha_v', 'Halpha_v_err'])

fix, ax = plt.subplots(figsize=(7, 7), dpi=300)
for obj_id in obj_list:

    idx_k = k_df.MPTID == str(int(obj_id))
    k_flux, k_flux_err = k_df.loc[idx_k, [line_k, f'{line_k}_err']].values[0]

    idx_v = (v_df.index.get_level_values('id') == obj_id) & (v_df.index.get_level_values('line') == line_v)
    v_flux, v_err = v_df.loc[idx_v, ['profile_flux', 'profile_flux_err']].values[0]

    # ax.scatter(v_flux, k_flux, color='red')
    ax.errorbar(v_flux, k_flux, xerr=v_err, yerr=k_flux_err,  fmt='o', color='red')

    # Save the data
    out_df.loc[obj_id, 'Halpha_k'] = k_flux
    out_df.loc[obj_id, 'Halpha_k_err'] = k_flux_err
    out_df.loc[obj_id, 'Halpha_v'] = v_flux
    out_df.loc[obj_id, 'Halpha_v_err'] = v_err

    out_df.loc[obj_id, ['Hbeta_k', 'Hbeta_k_err']] = k_df.loc[idx_k, ['Hbflx', f'Hbflx_err']].values[0]

    idx_beta = (v_df.index.get_level_values('id') == obj_id) & (v_df.index.get_level_values('line') == 'H1_4861A')
    values = v_df.loc[idx_beta, ['profile_flux', 'profile_flux_err']].values
    print(values)
    if values.size > 0:
        out_df.loc[obj_id, ['Hbeta_v', 'Hbeta_v_err']] = values

# Get current limits
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
low = min(xmin, ymin)
high = max(xmax, ymax)
ax.plot([low, high], [low, high], 'k--', lw=1)  # black dashed line

ax.update({'xlabel': 'Vital', 'ylabel': 'Kelcey', 'title':r'$H\alpha$ comparison CAPERS'})
plt.savefig('Halpha_comparison_CAPERS.png')
lime.save_frame(f'comparison_flux_table.csv', out_df)
# ax.update({'xlabel': 'Vital', 'ylabel': 'Kelcey', 'title':r'$H\beta$ comparison CAPERS'})
# plt.savefig('Hbetaa_comparison_CAPERS.png')

plt.show()

# coso = '0,0,2,,,,,8.14907439325191,7.49041e-14,,CAPERS PRISM,nircam1-10370,1136.8455362706038,2.274346e-13,87842,8.38205844661602e-06,,8.170000000000018'
# str_lit = coso.split(',')
#
# k_df = pd.read_csv(k_fname, index_col=0)
# k_df = k_df.drop('level_0', axis=1)
# print(k_df.columns)
# print(k_df.Type.unique())
#
# idcs_capers = k_df.Type == 'CAPERS PRISM'
# idcs_ceers = k_df.Type == 'CEERS PRISM'
# IDs_capers = k_df.loc[idcs_capers, 'MPTID']
# IDs_ceers = k_df.loc[idcs_ceers, 'MPTID']
#
# x, x_err = k_df.loc[idcs_capers, 'zphot'].to_numpy(), None
# y, y_err = error_division(k_df.loc[idcs_capers, 'Haflx'].to_numpy(), k_df.loc[idcs_capers, 'Haflx_err'].to_numpy(),
#                           k_df.loc[idcs_capers, 'Hbflx'].to_numpy(), k_df.loc[idcs_capers, 'Hbflx_err'].to_numpy())
#
# print(y)
# print(y_err)
