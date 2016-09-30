from spectral_cube import SpectralCube
import astropy.units as u
import numpy as np
def preprocessor(cube, min_intensity = 0.0*u.K,
                 channels_per_sigma=15, sigma_extent = 4, meanscale = True):
    vaxis = cube.spectral_axis
    projection = (cube.with_mask(cube>min_intensity)).moment(axis=(1))
    meanspec = np.nansum(projection,axis=1)
    totalspec = np.nansum(meanspec)
    meanvel = np.nansum(meanspec*vaxis)/totalspec
    dispersion = (np.nansum(meanspec*(vaxis-meanvel)**2)/totalspec)**0.5
    normalized_axis = (vaxis-meanvel)/dispersion
    target_normalized_axis = np.linspace(-sigma_extent, sigma_extent,
                                         2*sigma_extent*channels_per_sigma+1)
    sortidx = np.argsort(normalized_axis)
    target_axis = np.interp(target_normalized_axis,normalized_axis[sortidx], vaxis[sortidx].value)*vaxis.unit
    outcube = cube.spectral_interpolate(target_axis,suppress_smooth_warning=True)
    return(outcube)

