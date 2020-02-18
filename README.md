# QSmooth

QSmooth is a smoothing algorithm for quasar spectra developed by Dominika Durovcikova and implemented in Durovcikova et al. 2019 (https://arxiv.org/abs/1912.01050).

## Requirements:

Please find the list of all required Python packages in the "requirements.txt" file. In addition, a working Latex installation is required for the plotting.
This code has been developed using Python 2.7.

## Description:

QSmooth first computes a running median with a bin size of 50 data points to capture the main continuum and emission features in the spectrum. It then performs a peak-finding procedure using the SciPy Python library (Jones et al. 2001) above the aforementioned running median border and interpolates the peaks to construct an upper envelope of the spectrum. This envelope is then subtracted from the spectrum. QSmooth then applies the RANSAC regressor algorithm (Fischler & Bolles 1981) from the Scikit-Learn Python package (Pedregosa et al. 2012) on the residuals, thus rejecting most absorption features in the
spectrum. The data points that are flagged as inliers by RANSAC are interpolated and smoothed by computing a running median with a bin size of 20 (and an reduced bin size around the Lya peak), thus creating the final smooth flux fit of the spectrum.

For SDSS quasars from a spec-PLATE-MJD-FIBER.fits file, we also provide a routine to mask out all sky lines listed in Table 30 of Stoughton et al. (2002) as well as all pixels that were flagged as highly uncertain by the SDSS pipelines.

For a detailed account of the smoothing procedure, please refer to Appendix B in Durovcikova et al. 2019 (https://arxiv.org/abs/1912.01050).

## Example usage

Please refer to "example.py" for an example usage of the algorithm.

To smooth and plot an example spectrum from an SDSS spec-PLATE-MJD-FIBER.fits file, run:

python example.py

## Contact

Please contact Dominika Durovcikova at dominika.durovcikova@gmail.com in case of questions/issues.
