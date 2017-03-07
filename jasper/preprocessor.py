
import astropy.units as u
import numpy as np
from spectral_cube import SpectralCube
import astropy.io.fits as fits


def preprocessor(cube, min_intensity=0.0 * u.K,
                 channels_per_sigma=15, sigma_extent=4, norm_intensity=True,
                 norm_percentile=95):
    vaxis = cube.spectral_axis

    # Mean velocity
    meanvel = np.nanmean(cube.moment1().quantity)
    # Velocity dispersion
    dispersion = np.nanmean(cube.linewidth_sigma().quantity)

    normalized_axis = (vaxis - meanvel) / dispersion
    target_normalized_axis = \
        np.linspace(-sigma_extent, sigma_extent,
                    2 * sigma_extent * channels_per_sigma + 1)
    sortidx = np.argsort(normalized_axis)
    target_axis = np.interp(target_normalized_axis,
                            normalized_axis[sortidx],
                            vaxis[sortidx].value)[::-1] * vaxis.unit

    outcube = \
        cube.spectral_interpolate(target_axis, suppress_smooth_warning=False)

    # What we really want is the interpolated cube, but pretend it still has
    # the original spectral resolution.
    new_header = cube.header.copy()
    new_header["NAXIS3"] = outcube.shape[0]
    new_header["CRPIX3"] = (outcube.shape[0] / 2) - \
        int(meanvel.value / cube.header["CDELT3"])

    outcube = SpectralCube.read(fits.PrimaryHDU(outcube.filled_data[:].value,
                                                new_header))

    # Optionally normalize by the 90th percentile
    if norm_intensity:
        # We're going to keep the brightness unit after the normalization
        # The moment making class is expecting a reasonable unit to be attached
        outcube = outcube / outcube.percentile(norm_percentile).value

    return outcube
