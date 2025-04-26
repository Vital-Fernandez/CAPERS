import lime
from pathlib import Path
from support.tools import read_prism_r_curve

# Load configuration
cfg_file = '../CAPERS_v2.toml'
capers_cfg = lime.load_cfg(cfg_file)

# Get sample information
version = capers_cfg['meta']['version']
sample = capers_cfg['meta']['sample']
pref_ext_list = capers_cfg['file_structure']['ext_1d_fits']

observations_folder = Path(capers_cfg['file_structure']['observations_folder'])
data_folder = Path(capers_cfg['file_structure']['data_folder'])

tables_folder = data_folder/sample/'tables'
source_folder = data_folder/sample/'source'
bands_folder = data_folder/sample/'bands'
logs_folder = data_folder/sample/'line_logs'
comps_folder = data_folder/sample/'comps'
R_interpolator = read_prism_r_curve(capers_cfg['file_structure']['prism_R_curve_file'], units_factor=10000)

prism_bands_df = lime.load_frame(data_folder/capers_cfg['file_structure']['prism_bands'])

fname = '/home/vital/Downloads/CAPERS-COSMOS_p4_s000119334_x1d_optext.fits'
z = 9.287

spec = lime.Spectrum.from_file(fname, instrument='nirspec', redshift=z)
spec.unit_conversion('AA', 'FLAM')
spec.res_power = R_interpolator(spec.wave.data)

# spec.check.bands('119334_bands.txt', ref_bands=prism_bands_df, fit_cfg=cfg_file, obj_cfg_prefix='119334', exclude_continua=False)

spec.fit.frame('119334_bands.txt', fit_cfg=cfg_file, default_cfg_prefix='119334', cont_from_bands=True, err_from_bands=True)
spec.save_frame('119334_deblend2.txt')
spec.fit.report()
print(spec.frame.amp)           
print(spec.fit.line.amp)        
                                # spec.plot.spectrum(bands='119334_bands.txt', rest_frame=True)
spec.plot.bands('O3_5008A', rest_frame=True, show_err=True)
spec.plot.bands('H1_4342A', rest_frame=True, show_err=True)
# spec.fit.log

spec.plot.spectrum()
