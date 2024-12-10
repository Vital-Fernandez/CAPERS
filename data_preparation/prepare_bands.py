import lime
import numpy as np

fname = '/home/vital/Dropbox/Astrophysics/Data/CAPERS/PRISM_ref_bands.txt'

bands_capers = lime.load_frame(fname)

bands_ref_air = lime.line_bands(vacuum_waves=False, update_labels=False)
bands_ref_vacuum = lime.line_bands(vacuum_waves=True, update_labels=True)

# bands_df.insert(0, 'wavelength', np.nan)
#
for label in bands_capers.index:
    line = lime.Line(label)
    print(line, line.transition_comp)
    bands_capers.loc[label, 'transition'] = line.transition_comp[0]

# print(bands_ref_air)

lime.save_frame(fname, bands_capers)