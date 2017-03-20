
'''
Use common observational values as the "distance" (like total intensity, peak
intensity, average linewidth).
'''

from spectral_cube import SpectralCube
import astropy.units as u
import numpy as np
import pandas as pd
import statsmodels.formula.api as sm
import matplotlib.pyplot as p

from analysis import make_coefplots


def return_vals(filename, min_intensity=0.0):
    '''
    Load a cube, return some basic properties.
    '''

    cube = SpectralCube.read(filename)
    cube = cube.with_mask(cube > min_intensity * cube.unit)
    vaxis = cube.spectral_axis

    total_flux = cube.sum()

    peak_flux = cube.max()

    proj = cube.moment0(axis=1)
    totalspec = np.nansum(proj, axis=1)
    meanvel = np.nansum(totalspec * vaxis) / np.nansum(totalspec)
    linewidth_sigma = np.sqrt(np.nansum(totalspec * (vaxis - meanvel)**2) /
                              np.nansum(totalspec)) / cube.header["CDELT3"]

    return total_flux, peak_flux, linewidth_sigma


def append_design(design_file, df):
    '''
    Add the design column onto a DataFrame.
    '''

    des = pd.read_csv(design_file, index_col=0)
    param_names = []

    # Switch to the short hand names used in the rest of the analysis pipeline.
    factors = ['Solenoidal Fraction', 'Virial Parameter', 'k', 'Mach Number',
               'Plasma Beta']
    des = des[factors]

    param_names = ['sf', 'vp', 'k', 'm', 'pb']

    des.columns = pd.Index(param_names)

    for col in des.columns:
        df[col] = des[col]

    return df, param_names


def distance(value1, value2):
    return np.abs(value1 - value2)


if __name__ == "__main__":

    import sys
    import os

    from jasper.analysis_funcs import files_sorter

    try:
        path_to_data = sys.argv[1]
    except IndexError:
        path_to_data = raw_input("Give path to folder with simulation set: ")

    set_name = os.path.basename(os.path.normpath(path_to_data))
    print(set_name)

    try:
        design_file = sys.argv[2]
    except IndexError:
        design_file = raw_input("Give path and name of design file: ")

    try:
        output_path = sys.argv[3]
    except IndexError:
        output_path = raw_input("Give output path to save to: ")

    # Check if the products are already saved. If so, just make the summary
    # plots.
    dist_files = [os.path.join(output_path,
                               "{0}_face{1}_params.csv".format(set_name, face))
                  for face in [0, 1, 2]]
    if all([os.path.exists(file) for file in dist_files]):
        print("Results are already saved in the output directory. Skipping"
              " re-run. If a re-run is needed, delete the output files first.")

    else:

        fiducials, designs = files_sorter(path_to_data, timesteps='max',
                                          append_prefix=True)[:2]

        # Compute the time-averaged properties.
        design_params = dict.fromkeys(designs)
        fiducial_params = dict.fromkeys(fiducials)
        distances = dict.fromkeys(designs)

        # Fiducial properties
        for face in fiducials.keys():
            fiducial_params[face] = dict.fromkeys(fiducials[face])
            for sim in fiducials[face]:
                tfluxes = np.empty((len(fiducials[face][sim])))
                pfluxes = np.empty((len(fiducials[face][sim])))
                sigmas = np.empty((len(fiducials[face][sim])))
                for i, cube in enumerate(fiducials[face][sim]):
                    tflux, pflux, sigma = return_vals(cube)
                    tfluxes[i] = tflux.value
                    pfluxes[i] = pflux.value
                    sigmas[i] = sigma.value

                fiducial_params[face][sim] = {"tflux": tfluxes.mean(),
                                              "pflux": pfluxes.mean(),
                                              "sigma": sigmas.mean()}

        # Design properties
        for face in designs.keys():
            design_params[face] = dict.fromkeys(designs[face])
            for sim in designs[face]:
                tfluxes = np.empty((len(designs[face][sim])))
                pfluxes = np.empty((len(designs[face][sim])))
                sigmas = np.empty((len(designs[face][sim])))
                for i, cube in enumerate(designs[face][sim]):
                    tflux, pflux, sigma = return_vals(cube)
                    tfluxes[i] = tflux.value
                    pfluxes[i] = pflux.value
                    sigmas[i] = sigma.value

                design_params[face][sim] = {"tflux": tfluxes.mean(),
                                            "pflux": pfluxes.mean(),
                                            "sigma": sigmas.mean()}

                # design_params[face][sim] = {"tflux": tfluxes,
                #                             "pflux": pfluxes,
                #                             "sigma": sigmas}

        # Compute measures of distance
        for face in designs.keys():
            distances[face] = dict.fromkeys(fiducial_params[face])
            individ_dfs = []
            for fid in fiducial_params[face]:
                distances[face][fid] = dict.fromkeys(designs[face])
                for sim in designs[face]:
                    name = "{0}_{1}".format(sim, fid)
                    distances[face][fid][sim] = {}
                    distances[face][fid][sim]["dist_tflux"] = \
                        distance(design_params[face][sim]["tflux"],
                                 fiducial_params[face][fid]["tflux"])
                    distances[face][fid][sim]["dist_pflux"] = \
                        distance(design_params[face][sim]["pflux"],
                                 fiducial_params[face][fid]["pflux"])
                    distances[face][fid][sim]["dist_sigma"] = \
                        distance(design_params[face][sim]["sigma"],
                                 fiducial_params[face][fid]["sigma"])
                    # Track which fiducial as this is the random effect.
                    distances[face][fid][sim]["Cube"] = fid

                # Convert to dataframes
                df = pd.DataFrame(distances[face][fid]).T
                df, param_names = append_design(design_file, df)
                individ_dfs.append(df)

            df_params = pd.concat(individ_dfs, ignore_index=True)
            save_name = \
                os.path.join(output_path,
                             "{0}_face{1}_params.csv".format(set_name, face))
            df_params.to_csv(save_name)

            model = "*".join(param_names)

            result_tflux = sm.mixedlm(formula="dist_tflux ~ {}".format(model),
                                      data=df_params,
                                      groups=df_params["Cube"]).fit()
            save_name = \
                os.path.join(output_path,
                             "{0}_face{1}_tflux_fit.pkl".format(set_name, face))
            result_tflux.save(save_name)

            result_pflux = sm.mixedlm(formula="dist_pflux ~ {}".format(model),
                                      data=df_params,
                                      groups=df_params["Cube"]).fit()
            save_name = \
                os.path.join(output_path,
                             "{0}_face{1}_pflux_fit.pkl".format(set_name, face))
            result_pflux.save(save_name)

            result_sigma = sm.mixedlm(formula="dist_sigma ~ {}".format(model),
                                      data=df_params,
                                      groups=df_params["Cube"]).fit()
            save_name = \
                os.path.join(output_path,
                             "{0}_face{1}_sigma_fit.pkl".format(set_name, face))
            result_sigma.save(save_name)

            # Make some model plots
            df_tvals = pd.DataFrame({"dist_sigma": result_sigma.tvalues,
                                     "dist_pflux": result_pflux.tvalues,
                                     "dist_tflux": result_tflux.tvalues})

            # Order the terms
            ind = list(df_tvals.index)
            ind.sort(key=lambda x: x.count(":"))
            ind = pd.Index(ind)
            df_tvals = df_tvals.ix[ind]
            # Remove intercept
            df_tvals = df_tvals.drop("Intercept")
            df_tvals = df_tvals.drop("Intercept RE")

            # Save the t-values
            save_name = \
                os.path.join(output_path,
                             "{0}_face{1}_tvalues.csv".format(set_name, face))
            df_tvals.to_csv(save_name)

    # Load in face 0 and face 2
    df_params_0 = pd.read_csv(dist_files[0], index_col=0)
    # df_params_2 = pd.read_csv(dist_files[2], index_col=0)

    # Add a face column
    # df_params_0["fc"] = np.ones((df_params_0.shape[0]), dtype=np.int) * 0
    # df_params_2["fc"] = np.ones((df_params_0.shape[0]), dtype=np.int) * 1

    df_params = df_params_0  # pd.concat([df_params_0, df_params_2])

    statistics = ["dist_sigma", "dist_pflux", "dist_tflux"]

    basename = os.path.basename(os.path.realpath(path_to_data))

    # Getting cases with many significant terms, and no lower order
    # insignificant terms. Just make coefplots for now.
    # effect_plots(df_params, df_tvals, statistics=statistics,
    #              save=True, out_path=output_path,
    #              output_name=basename, params=['sf', 'vp', 'k', 'm', 'pb'])
    # p.close()


    import seaborn as sb
    # sb.set_context("poster", font_scale=1.2)
    sb.set_context("paper")
    sb.set(font='Times New Roman', style='ticks', font_scale=1.2)

    figsize = (4.2, 6.0)

    make_coefplots(df_params, save=True, out_path=output_path,
                   output_name=basename,
                   min_tvalue=3.46, endog_formula='m*k*pb*vp*sf',  # *fc',
                   statistics=statistics, figsize=figsize)
