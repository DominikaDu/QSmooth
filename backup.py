	# Smooths an SDSS spec-PLATE-MJD-FIBER.fits file according to procedure outlined in Appendix B
    # of Durovcikova et al. 2019 (https://arxiv.org/abs/1912.01050). In this process, a mask can be 
    # used to reject some data from the spec file.

	hdu_raw = pf.open(str(path)+str(filename))
    if len(mask)>0:
	    hdu = hdu_raw[1].data[mask]
    else:
        hdu = hdu_raw[1].data
    # Either calibrate to rest wavelengths based on redshift from a particular line, e.g. Mg-II
    # https://data.sdss.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/PLATE4/RUN1D/spZline.html
	# x = hdu['loglam'] - np.log10(1+hdu_raw[3].data['LINEZ'][5])
    # Or calibrate to rest wavelengths based on the redshift from the SDSS pipelines.
    x = hdu['loglam'] - np.log10(1+hdu_raw[2].data['Z'])
	y = hdu['flux']
	err = hdu['ivar']