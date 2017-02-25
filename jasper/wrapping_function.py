
'''
Wrapper to run all of the statistics between two data sets.
For large-scale comparisons, this is the function that should be called.
'''

import numpy as np
from copy import copy


from turbustat.statistics import (Wavelet_Distance, MVC_Distance,
                                  PSpec_Distance,
                                  BiSpectrum_Distance, GenusDistance,
                                  DeltaVariance_Distance, VCA_Distance,
                                  VCS_Distance, Tsallis_Distance,
                                  StatMoments_Distance, PCA_Distance,
                                  SCF_Distance, Cramer_Distance,
                                  DendroDistance, PDF_Distance,
                                  statistics_list)


def stats_wrapper(dataset1, dataset2, fiducial_models=None,
                  statistics=None, multicore=False, vca_break=None,
                  vcs_break=None, vcs_regrid=[False, False],
                  dendro_params=None,
                  periodic_bounds=[True, True],
                  noise_value=[-np.inf, -np.inf],
                  dendro_saves=[None, None],
                  inertial_range=[[None] * 2, [None] * 2],
                  spatial_range=[[None] * 2, [None] * 2]):
    '''
    Function to run all of the statistics on two datasets.
    Each statistic is run with set inputs. This function needs to be altered
    to change the inputs.

    Parameters
    ----------
    dataset1 : dict
        Contains the cube and all of its property arrays.
    dataset2 : dict
        See dataset1
    fiducial_models : dict, optional
        Models for dataset1. Avoids recomputing when comparing
        many sets to dataset1.
    statistics : list, optional
        List of all of the statistics to use. If None, all are run.
    multicore : bool, optional
        If the wrapper is being used in parallel, this disables
        returning model values for dataset1.
    vcs_break : float, optional
        Pass an initial guess for the location of the VCS break.
    vcs_regrid : list of bools, optional
        The simulated cubes lack information on the smallest spectral scales.
        When enabled, the cube is downsampled by a factor of 5 spectrally
        before running VCS.
    dendro_params : dict or list, optional
        Provides parameters to use when computing the initial dendrogram.
        If different parameters are required for each dataset, the
        the input should be a list containing the two dictionaries.
    periodic_bounds : list of bools
        Set whether the boundaries should be handled as 'continuous' (True) or
        not ('cut' or 'fill'; False).
    cleanup : bool, optional
        Delete distance classes after running.
    '''

    if statistics is None:  # Run them all
        statistics = statistics_list

    distances = {}

    # Calculate the fiducial case and return it for later use
    if fiducial_models is None:
        fiducial_models = {}
        for statistic in statistics:
            if "PDF" in statistic:
                statistic = "PDF"
            elif statistic == "Skewness" or statistic == "Kurtosis":
                statistic = "stat_moments"
            elif "Dendrogram" in statistic:
                statistic = "Dendrogram"
            elif "DeltaVariance_Centroid" in statistic:
                statistic = "DeltaVariance_Centroid"
            elif "DeltaVariance" in statistic and "Centroid" not in statistic:
                statistic = "DeltaVariance"
            fiducial_models[statistic] = None

    if any("Wavelet" in s for s in statistics):
        wavelet_distance = \
            Wavelet_Distance(dataset1["moment0"],
                             dataset2["moment0"],
                             fiducial_model=fiducial_models["Wavelet"],
                             xlow=spatial_range[0],
                             xhigh=spatial_range[1])
        wavelet_distance.distance_metric()
        distances["Wavelet"] = wavelet_distance.distance
        if not multicore:
            fiducial_models["Wavelet"] = copy(wavelet_distance.wt1)

        del wavelet_distance

    if any("MVC" in s for s in statistics):
        mvc_distance = \
            MVC_Distance(dataset1, dataset2,
                         fiducial_model=fiducial_models["MVC"],
                         low_cut=inertial_range[0],
                         high_cut=inertial_range[1])
        mvc_distance.distance_metric()
        distances["MVC"] = mvc_distance.distance
        if not multicore:
            fiducial_models["MVC"] = copy(mvc_distance.mvc1)

        del mvc_distance

    if any("PSpec" in s for s in statistics):
        pspec_distance = \
            PSpec_Distance(dataset1["moment0"],
                           dataset2["moment0"],
                           fiducial_model=fiducial_models['PSpec'],
                           low_cut=inertial_range[0],
                           high_cut=inertial_range[1])
        pspec_distance.distance_metric()
        distances["PSpec"] = pspec_distance.distance
        if not multicore:
            fiducial_models["PSpec"] = copy(pspec_distance.pspec1)

        del pspec_distance

    if any("Bispectrum" in s for s in statistics):
        bispec_distance = \
            BiSpectrum_Distance(dataset1["moment0"],
                                dataset2["moment0"],
                                fiducial_model=fiducial_models['Bispectrum'])
        bispec_distance.distance_metric()
        distances["Bispectrum"] = bispec_distance.distance
        if not multicore:
            fiducial_models["Bispectrum"] = copy(bispec_distance.bispec1)

        del bispec_distance

    if any("DeltaVariance_Slope" in s for s in statistics) or \
       any("DeltaVariance_Curve" in s for s in statistics):

        # Check for how boundaries should be handled.
        boundary1 = 'wrap' if periodic_bounds[0] else 'cut'
        boundary2 = 'wrap' if periodic_bounds[1] else 'cut'

        delvar_distance = \
            DeltaVariance_Distance(dataset1["moment0"],
                                   dataset2["moment0"],
                                   weights1=dataset1["moment0_error"][0],
                                   weights2=dataset2["moment0_error"][0],
                                   fiducial_model=fiducial_models["DeltaVariance"],
                                   xlow=spatial_range[0],
                                   xhigh=spatial_range[1],
                                   boundary=[boundary1, boundary2])

        delvar_distance.distance_metric()
        distances["DeltaVariance_Curve"] = delvar_distance.curve_distance
        distances["DeltaVariance_Slope"] = delvar_distance.slope_distance
        if not multicore:
            fiducial_models["DeltaVariance"] = copy(delvar_distance.delvar1)

        del delvar_distance

    if any("DeltaVariance_Centroid_Slope" in s for s in statistics) or \
       any("DeltaVariance_Centroid_Curve" in s for s in statistics):

        # Check for how boundaries should be handled.
        boundary1 = 'wrap' if periodic_bounds[0] else 'cut'
        boundary2 = 'wrap' if periodic_bounds[1] else 'cut'

        delvar_distance = \
            DeltaVariance_Distance(dataset1["centroid"],
                                   dataset2["centroid"],
                                   weights1=dataset1["centroid_error"][0],
                                   weights2=dataset2["centroid_error"][0],
                                   fiducial_model=fiducial_models["DeltaVariance_Centroid"],
                                   xlow=spatial_range[0],
                                   xhigh=spatial_range[1],
                                   boundary=[boundary1, boundary2])
        delvar_distance.distance_metric()
        distances["DeltaVariance_Centroid_Curve"] = \
            delvar_distance.curve_distance
        distances["DeltaVariance_Centroid_Slope"] = \
            delvar_distance.slope_distance
        if not multicore:
            fiducial_models["DeltaVariance_Centroid"] = \
                copy(delvar_distance.delvar1)

        del delvar_distance

    if any("Genus" in s for s in statistics):
        genus_distance = \
            GenusDistance(dataset1["moment0"],
                          dataset2["moment0"],
                          fiducial_model=fiducial_models['Genus'])
        genus_distance.distance_metric()
        distances["Genus"] = genus_distance.distance
        if not multicore:
            fiducial_models["Genus"] = copy(genus_distance.genus1)

        del genus_distance

    if any("VCS" in s for s in statistics):

        # Regrid the cube to lower spectral resolution
        if any(vcs_regrid):
            from spectral_cube import SpectralCube
            import astropy.io.fits as fits

            from jasper.analysis_funcs import spectral_regrid_cube

            if vcs_regrid[0]:
                cube1_hdu = fits.PrimaryHDU(dataset1["cube"][0],
                                            header=dataset1["cube"][1])
                cube1 = SpectralCube.read(cube1_hdu)

                cube1 = spectral_regrid_cube(cube1, 100)
            else:
                cube1 = dataset1["cube"]

            if vcs_regrid[1]:
                cube2_hdu = fits.PrimaryHDU(dataset2["cube"][0],
                                            header=dataset2["cube"][1])
                cube2 = SpectralCube.read(cube2_hdu)

                cube2 = spectral_regrid_cube(cube2, 100)
            else:
                cube2 = dataset2["cube"]

        else:
            cube1 = dataset1["cube"]
            cube2 = dataset2["cube"]

        vcs_distance = VCS_Distance(cube1,
                                    cube2,
                                    breaks=vcs_break,
                                    fiducial_model=fiducial_models['VCS'])
        vcs_distance.distance_metric()
        distances["VCS"] = vcs_distance.distance
        distances["VCS_Small_Scale"] = vcs_distance.small_scale_distance
        distances["VCS_Large_Scale"] = vcs_distance.large_scale_distance
        distances["VCS_Break"] = vcs_distance.break_distance
        if not multicore:
            fiducial_models["VCS"] = copy(vcs_distance.vcs1)

        del vcs_distance

    if any("VCA" in s for s in statistics):
        vca_distance = VCA_Distance(dataset1["cube"],
                                    dataset2["cube"],
                                    breaks=vca_break,
                                    fiducial_model=fiducial_models['VCA'],
                                    low_cut=inertial_range[0],
                                    high_cut=inertial_range[1])
        vca_distance.distance_metric()
        distances["VCA"] = vca_distance.distance
        if not multicore:
            fiducial_models["VCA"] = copy(vca_distance.vca1)

        del vca_distance

    if any("Tsallis" in s for s in statistics):
        tsallis_distance = \
            Tsallis_Distance(dataset1["moment0"],
                             dataset2["moment0"],
                             fiducial_model=fiducial_models['Tsallis'])
        tsallis_distance.distance_metric()
        distances["Tsallis"] = tsallis_distance.distance
        if not multicore:
            fiducial_models["Tsallis"] = copy(tsallis_distance.tsallis1)

        del tsallis_distance

    if any("Skewness" in s for s in statistics) or\
       any("Kurtosis" in s for s in statistics):
        moment_distance = \
            StatMoments_Distance(dataset1["moment0"],
                                 dataset2["moment0"], 5,
                                 fiducial_model=fiducial_models['stat_moments'])
        moment_distance.distance_metric()
        distances["Skewness"] = moment_distance.skewness_distance
        distances["Kurtosis"] = moment_distance.kurtosis_distance
        if not multicore:
            fiducial_models["stat_moments"] = \
                copy(moment_distance.moments1)

        del moment_distance

    if any("PCA" in s for s in statistics):
        pca_distance = \
            PCA_Distance(dataset1["cube"],
                         dataset2["cube"],
                         fiducial_model=fiducial_models['PCA'])
        pca_distance.distance_metric()
        distances["PCA"] = pca_distance.distance
        if not multicore:
            fiducial_models["PCA"] = pca_distance.pca1

        del pca_distance

    if any("SCF" in s for s in statistics):

        boundary1 = "continuous" if periodic_bounds[0] else 'cut'
        boundary2 = "continuous" if periodic_bounds[1] else 'cut'

        scf_distance = \
            SCF_Distance(dataset1["cube"],
                         dataset2["cube"],
                         boundary=[boundary1, boundary2],
                         fiducial_model=fiducial_models['SCF'])
        scf_distance.distance_metric()
        distances["SCF"] = scf_distance.distance
        if not multicore:
            fiducial_models["SCF"] = copy(scf_distance.scf1)

        del scf_distance

    if any("Cramer" in s for s in statistics):
        cramer_distance = \
            Cramer_Distance(dataset1["cube"],
                            dataset2["cube"],
                            noise_value1=noise_value[0],
                            noise_value2=noise_value[1]).distance_metric()
        distances["Cramer"] = cramer_distance.distance

        del cramer_distance

    if any("Dendrogram_Hist" in s for s in statistics) or \
       any("Dendrogram_Num" in s for s in statistics):

        if dendro_saves[0] is None:
            input1 = dataset1["cube"]

        elif isinstance(dendro_saves[0], str):
            input1 = dendro_saves[0]
        else:
            raise UserWarning("dendro_saves must be the filename of the"
                              " saved file.")

        if dendro_saves[1] is None:
            input2 = dataset2["cube"]
        elif isinstance(dendro_saves[1], str):
            input2 = dendro_saves[1]
        else:
            raise UserWarning("dendro_saves must be the filename of the"
                              " saved file.")

        dendro_distance = \
            DendroDistance(input1, input2,
                           dendro_params=dendro_params,
                           fiducial_model=fiducial_models['Dendrogram'],
                           periodic_bounds=periodic_bounds)
        dendro_distance.distance_metric()

        distances["Dendrogram_Hist"] = dendro_distance.histogram_distance
        distances["Dendrogram_Num"] = dendro_distance.num_distance
        if not multicore:
            fiducial_models["Dendrogram"] = copy(dendro_distance.dendro1)

        del dendro_distance

    if any("PDF_Hellinger" in s for s in statistics) or \
       any("PDF_KS" in s for s in statistics) or \
       any("PDF_Lognormal" in s for s in statistics):  # or \
       # any("PDF_AD" in s for s in statistics):
        pdf_distance = \
            PDF_Distance(dataset1["moment0"],
                         dataset2["moment0"],
                         min_val1=2 * noise_value[0],
                         min_val2=2 * noise_value[1])

        pdf_distance.distance_metric()

        distances["PDF_Hellinger"] = pdf_distance.hellinger_distance
        distances["PDF_KS"] = pdf_distance.ks_distance
        distances["PDF_Lognormal"] = pdf_distance.lognormal_distance
        # distances["PDF_AD"] = pdf_distance.ad_distance
        if not multicore:
                fiducial_models["PDF"] = copy(pdf_distance.PDF1)

        del pdf_distance

    if multicore:
        return distances
    else:
        return distances, fiducial_models
