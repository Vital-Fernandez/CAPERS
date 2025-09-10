from pathlib import Path
import lime

lime.theme.set_style('dark')

# Declare input files
file_list = ['CAPERS-EGS_p5_s000025297_x1d_optext_corrected.fits',
             "CAPERS-COSMOS_p1_s000052185_corrected_x1d_optext.fits",
             "CAPERS-COSMOS_p1_s000300001_corrected_x1d_optext.fits"]
redshift_list = [9.918, 9.808, 9.812]
cfg_fname = "Callum_sample.toml"

for i, file in enumerate(file_list):

    fname =  Path(file)
    stem = fname.stem
    redshift = redshift_list[i]

    # Load observation
    spec = lime.Spectrum.from_file(fname, instrument='nirspec', redshift=redshift)
    spec.unit_conversion('AA', "FLAM")

    # Generate bands
    # bands = spec.retrieve.line_bands(fit_cfg=cfg_fname, vacuum_waves=True)
    # spec.plot.spectrum(rest_frame=True, bands=bands)

    # Manual Review bands
    # spec.check.bands(f"{stem}_bands.txt", fit_cfg=cfg_fname, exclude_continua=True, vacuum_waves=True)

    # Automatic line review
    # spec.fit.continuum(degree_list=[3, 4, 4, 5], emis_threshold=[4, 3, 2, 1.5], plot_steps=True)
    # bands_peaks = spec.infer.peaks_troughs(bands_fname, plot_steps=True)

    # Fit lines
    spec.fit.frame(f"{stem}_bands.txt", cfg_fname, cont_from_bands=False)
    spec.plot.spectrum(log_scale=False, rest_frame=True, output_address=f"{stem}_plot_spectrum.png")
    spec.plot.grid(rest_frame=True, output_address=f"{stem}_plot_grid.png")
    spec.save_frame(f"{stem}_log.txt")