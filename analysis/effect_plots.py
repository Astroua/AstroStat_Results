
'''
Read in the effects from the fit and plot the important effects.
'''
from __future__ import division

import matplotlib as mpl
import matplotlib.pyplot as p
import matplotlib.colors as cols
import matplotlib.cm as cm
import numpy as np
from pandas import read_csv
import warnings
import os

from turbustat.statistics import statistics_list

p.rcParams.update({'font.size': 14})


def effect_plots(distance_file, effects_file, min_tscore=2.0, statistics=None,
                 params=["fc", "pb", "m", "k", "sf", "vp"], save=False,
                 out_path=None, output_name=None):
    '''
    Creates a series of plots for the important effects in the model.
    '''

    if isinstance(distance_file, str):
        distances = read_csv(distance_file)
    else:
        distances = distance_file

    if isinstance(effects_file, str):
        effects = read_csv(effects_file)
    else:
        effects = effects_file

    # Extract the design matrix
    design = distances.T.ix[params].T

    # Now drop the design matrix and the Fiducial number
    distances = distances.T.drop(params).T
    try:
        distances = distances.T.drop(["Cube"]).T
    except ValueError:
        pass

    # Get the model effects from the index
    model_effects = effects.index

    # Replace the shorthand names
    rep_name = {'fc': "F", "pb": r'$\beta$', 'm': r'$\mathcal{M}$',
                'k': r'$k$', 'sf': r'$\zeta$',
                'vp': r'$\alpha$'}

    if statistics is None:
        statistics = distances.columns

    # Now loop through the statistics in each file
    for stat in statistics:
        try:
            response = effects[stat].abs()
        except KeyError:
            warnings.warn("%s effects not found. Skipping" % (stat))
            continue

        if response.max() == 0.0:
            warnings.warn("Largest t-value is 0. Assuming failed fit for"
                          " {}.".format(stat))
            continue

        # Ignore higher than 2nd order effects
        cut_to_2nd = int(len(params) + 0.5 * (len(params) * (len(params) - 1)))
        response = response[:cut_to_2nd]

        # Find the most important effects 2nd order effects
        # Values are ~ z-scores. Significant effects have an absolute value
        # greater than 2.
        imp_inters = response[len(params):].order(ascending=False)[:3]

        # Create colours based on the absolute values of the responses
        milagro = \
            colormap_milagro(np.log10(response.min()),
                             np.log10(response.max()),
                             np.log10(min_tscore))
        cNorm = cols.Normalize(vmin=np.log10(response.min()),
                               vmax=np.log10(response.max()))
        scalMap = cm.ScalarMappable(norm=cNorm, cmap=milagro)

        # Create plots for the main effects, regardless of importance
        fig, axes = p.subplots(3, 3, sharex=True)
        for i, (param, ax) in enumerate(zip(params + list(imp_inters.index),
                                            axes.flatten())):
            if i < len(params):
                low_data = distances[stat][design[param] == -1]
                high_data = distances[stat][design[param] == 1]

                ax.plot([-1, 1], [low_data.mean(), high_data.mean()],
                        marker="D",
                        color=scalMap.to_rgba(np.log10(response[param])),
                        lw=0)

                # Plot the slope
                ax.plot([-1, 1], [low_data.mean(), high_data.mean()],
                        color=scalMap.to_rgba(np.log10(response[param])),
                        label=rep_name[param])

                ax.set_title(rep_name[param])

            else:
                param1 = param.split(":")[0]
                param2 = param.split(":")[-1]

                low_low = distances[stat][np.logical_and(design[param1] == -1,
                                                         design[param2] == -1)]
                low_high = distances[stat][np.logical_and(design[param1] == 1,
                                                          design[param2] == -1)]
                high_low = distances[stat][np.logical_and(design[param1] == -1,
                                                          design[param2] == 1)]
                high_high = distances[stat][np.logical_and(design[param1] == 1,
                                                          design[param2] == 1)]

                ax.plot([-1, 1], [low_low.mean(), low_high.mean()], marker="D",
                        color=scalMap.to_rgba(np.log10(imp_inters[param])))
                ax.plot([-1, 1],
                        [high_low.mean(), high_high.mean()], marker="s",
                        color=scalMap.to_rgba(np.log10(imp_inters[param])))

                ax.plot([-1, 1], [low_low.mean(), low_high.mean()],
                        color=scalMap.to_rgba(np.log10(imp_inters[param])))
                ax.plot([-1, 1], [high_low.mean(), high_high.mean()],
                        color=scalMap.to_rgba(np.log10(imp_inters[param])))

                ax.set_title(rep_name[param1]+" : "+rep_name[param2])

            ax.set_xlim([-2, 2])
            ax.set_ylim([distances[stat].min(), distances[stat].max()])

            # Y-ticks at row beginnings
            if i == 0 or i == 3 or i == 6:
                pass
            else:
                ax.set_yticklabels([])
            # X-ticks on bottom row
            if i < 3:
                ax.set_xticks([])
            else:
                ax.set_xticks([-1, 1])

        fig.subplots_adjust(right=0.82)
        cax = fig.add_axes([0.85, 0.1, 0.02, 0.8])
        cb = mpl.colorbar.ColorbarBase(cax, cmap=milagro, norm=cNorm)
        cb.set_ticks(np.log10(response))
        # Avoid white lines in the pdf rendering
        cb.solids.set_edgecolor("face")

        max_resp = np.max(np.log10(response))
        min_resp = np.min(np.log10(response))
        text_posns = np.linspace(0.1, 0.9, len(response))[::-1]

        # model_effects = model_effects[:21][np.asarray(response).argsort()]
        # model_effects = model_effects[::-1]

        response = response.order(ascending=False)
        model_effects = response.index

        for effect, value, tpos in zip(model_effects, np.log10(response), text_posns):
            if ":" in effect:
                splitted = effect.split(":")
                param1 = splitted[0]
                param2 = splitted[-1]
                label = rep_name[param1] + " : " + rep_name[param2]
            else:
                label = rep_name[effect]
            cax.annotate(label,
                         xy=(0.88, (value - min_resp)/(max_resp-min_resp)),
                         xytext=(0.9, tpos), textcoords='figure fraction',
                         arrowprops=dict(facecolor='k',
                                         width=0.05, alpha=1.0, headwidth=0.1),
                         horizontalalignment='left',
                         verticalalignment='top')
        cb.set_ticklabels([])
        cb.ax.tick_params(labelsize=10, colors='white', length=10)

        cax.annotate(np.round(10**min_resp, 1), xy=(0.83, 0.07),
                     xytext=(0.83, 0.07), textcoords='figure fraction')

        cax.annotate(np.round(10**max_resp, 1), xy=(0.83, 0.91),
                     xytext=(0.83, 0.91), textcoords='figure fraction')

        if save:
            out_name = "full_factorial_{0}_modeleffects.pdf".format(stat)
            if output_name is not None:
                out_name = "{0}_{1}".format(output_name, out_name)

            if out_path is not None:
                out_name = os.path.join(out_path, out_name)
            fig.savefig(out_name)
            p.close()
        else:
            fig.canvas.set_window_title("Model results for: "+stat)
            fig.show()
            p.show()


def map_all_results(effects_file, min_tscore=2.0, max_tscore=10.0,
                    save_name=None,
                    max_order=2, statistics=statistics_list,
                    normed=False, out_path=None,
                    params={"fc": "F", "pb": r"$\beta$",
                            "m": r"$\mathcal{M}$", "k": r"$k$",
                            "sf": r"$\zeta$", "vp": r"$\alpha$"}):

    if isinstance(effects_file, str) or isinstance(effects_file, unicode):
        effects = read_csv(effects_file)
    else:
        effects = effects_file

    # We only car about the absolute value
    effects = effects.abs()

    # Get the model effects from the index
    model_effects = effects.index

    if statistics is None:
        statistics = list(effects.columns)

    # Alter non-latex friendly strings
    stat_labels = []
    for stat in statistics:
        stat_labels.append(stat.replace("_", " "))

    # Find the cutoff if a maximum order is given
    if max_order is not None:
        for posn, effect in enumerate(model_effects):
            splitted = effect.split(":")

            if len(splitted) > max_order:
                cutoff_posn = posn
                break

        model_effects = model_effects[:cutoff_posn]

    values = np.empty((len(statistics), len(model_effects)))

    for i, stat in enumerate(statistics):
        if normed:
            for j, effect in enumerate(model_effects):
                if ":" in effect:
                    splitted = effect.split(":")[::2]
                    norm_factor = 1
                    for param in splitted:
                        norm_factor *= effects[stat][param]
                    value = np.power(effects[stat][effect] / norm_factor,
                                     1 / float(len(splitted)))
                else:
                    value = effects[stat][effect]
                values[i, j] = value
        else:
            values[i, :] = effects[stat][:len(model_effects)]

    # Change parameter names to those provided in params
    if params is not None:

        for param in params:
            model_effects = [eff.replace(param, params[param])
                             for eff in model_effects]

    # Use the maximum from the whole set.
    if max_tscore is None:
        max_tscore = values.max()
        print(max_tscore)

    milagro = \
        colormap_milagro(0, max_tscore, min_tscore)

    w, h = mpl.rcParams["figure.figsize"]
    hsize = lambda n: n * (h / 2)
    wsize = lambda n: n * (w / (4 * (n / 5)))

    n_terms = values.shape[1]

    fig, ax = p.subplots(1, 1, figsize=(hsize(1.5), wsize(n_terms)))

    # Flip dimension to put lowest terms at the bottom
    p.imshow(values.T[::-1], vmin=0, vmax=max_tscore, cmap=milagro,
             interpolation="nearest")
    p.yticks(np.arange(len(model_effects)), model_effects[::-1], rotation=0)
    p.xticks(np.arange(len(statistics)), stat_labels, rotation=90)
    cbar = p.colorbar(fraction=0.05, shrink=0.9)
    cbar.ax.set_ylabel(r'$t$-value')
    # cbar.ax.tick_params(labelsize=fontsize)
    # Avoid white lines in the pdf rendering
    cbar.solids.set_edgecolor("face")

    p.tight_layout()

    # Save if save_name has been given
    if save_name is not None:

        if out_path is None:
            out_path = ""

        save_name = os.path.join(out_path, save_name)

        p.savefig(save_name)
        p.close()
    else:
        p.show()


def make_coefplots(data, endog_formula="fc*m*k*pb*vp*sf",
                   statistics=None, min_tvalue=2, save=False,
                   out_path=None, output_name=None):
    '''
    Create coefplots for all statistics in a given data file.
    '''

    if isinstance(data, str):
        data = read_csv(data)

    if statistics is None:
        statistics = statistics_list

    _reset_interactive = False
    if p.isinteractive():
        _reset_interactive = True
        p.ioff()

    for stat in statistics:

        if stat not in data.columns:
            warnings.warn("{} is not in the given data.".format(stat))
            continue

        coefplot("{0} ~ {1}".format(stat, endog_formula), data,
                 min_tvalue=min_tvalue)

        if save:
            out_name = "full_factorial_{0}_coefplot.pdf".format(stat)
            if output_name is not None:
                out_name = "{0}_{1}".format(output_name, out_name)

            if out_path is not None:
                out_name = os.path.join(out_path, out_name)
            p.savefig(out_name)
            p.close()
        else:
            p.show()

    if _reset_interactive:
        p.ion()


def coefplot(formula, data, intercept=False, ci=95, min_tvalue=2,
             mixed_effect="Cube", sig_color='r', nonsig_color='k'):
    """Plot the coefficients from a linear model.

    Parameters
    ----------
    formula : string
        patsy formula for ols model
    data : dataframe
        data for the plot; formula terms must appear in columns
    groupby : grouping object, optional
        object to group data with to fit conditional models
    intercept : bool, optional
        if False, strips the intercept term before plotting
    ci : float, optional
        size of confidence intervals
    palette : {seaborn color palette, color}, optional
        palette for the horizonal plots or a single color

    """

    try:
        import statsmodels.formula.api as sf
        _has_statsmodels = True
    except ImportError:
        _has_statsmodels = False

    if not _has_statsmodels:
        raise ImportError("The `coefplot` function requires statsmodels")

    import pandas as pd

    alpha = 1 - ci / 100
    if mixed_effect is None:
        model = sf.ols(formula, data).fit()
    else:
        model = sf.mixedlm(formula, data,
                           groups=data[mixed_effect]).fit(reml=False)
    coefs = model.params

    # Order by term order. This should be a default in statsmodels IMO.
    ind = list(coefs.index)
    ind.sort(key=lambda x: x.count(":"))
    ind = pd.Index(ind)
    coefs = coefs[ind]

    cis = model.conf_int(alpha).T[ind].T
    tvals = np.abs(model.tvalues)[ind]
    model_effects = ind

    # Possibly ignore the intercept
    if not intercept:
        coefs = coefs.ix[model_effects != "Intercept"]
        cis = cis.ix[model_effects != "Intercept"]
        tvals = tvals.ix[model_effects != "Intercept"]
        model_effects = model_effects[model_effects != "Intercept"]

        if mixed_effect is not None:
            coefs = coefs.ix[model_effects != "Intercept RE"]
            cis = cis.ix[model_effects != "Intercept RE"]
            tvals = tvals.ix[model_effects != "Intercept RE"]
            model_effects = model_effects[model_effects != "Intercept RE"]

    n_terms = len(coefs)

    rep_name = {'fc': "F", "pb": r'$\beta$', 'm': r'$\mathcal{M}$',
                'k': r'$k$', 'sf': r'$\zeta$',
                'vp': r'$\alpha$'}

    w, h = mpl.rcParams["figure.figsize"]
    hsize = lambda n: n * (h / 2)
    wsize = lambda n: n * (w / (4 * (n / 5)))

    fig, ax = p.subplots(1, 1, figsize=(hsize(1.5), wsize(n_terms)))
    for i, term in enumerate(coefs.index):
        if tvals[term] < min_tvalue:
            color = nonsig_color
            symbol = '^'
        else:
            color = sig_color
            symbol = 'o'

        low, high = cis.ix[term]
        ax.plot([low, high], [i, i], c=color,
                solid_capstyle="round", lw=2.5)
        ax.plot(coefs.ix[term], i, symbol, c=color, ms=8)
        # ax.plot([i, i], [low, high], c=color,
        #         solid_capstyle="round", lw=2.5)
        # ax.plot(i, coefs.ix[term], symbol, c=color, ms=8)
    ax.set_ylim(-.5, n_terms - .5)
    ax.set_xlabel("Coefficent Values")
    ax.axvline(0, ls="--", c="dimgray")
    ax.grid(True)

    # Change to the nice parameter labels
    altered_labels = []
    for mod in model_effects:
        if ":" in mod:
            new_mod = []
            for param in mod.split(":"):
                new_mod.append(rep_name[param])
        else:
            new_mod = [rep_name[mod]]
        altered_labels.append(":".join(new_mod))

    p.yticks(np.arange(len(model_effects)), altered_labels, rotation=0)

    # Add in a legend with the symbols wrt min_tval
    sig_artist = p.Line2D((0, 1), (0, 0), color=sig_color, marker='o',
                          linestyle='-')
    nonsig_artist = p.Line2D((0, 1), (0, 0), color=nonsig_color, marker='^',
                             linestyle='-')

    ax.legend([sig_artist, nonsig_artist],
              [r"$t$-value > {0:.2f}".format(min_tvalue),
               r"$t$-value < {0:.2f}".format(min_tvalue)], frameon=True,
              loc='best')

    p.tight_layout()


def colormap_milagro(vmin, vmax, vtransition, width=0.0001, huestart=0.6):

    """
    Gammapy License:

    Copyright (c) 2014, Gammapy Developers All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

        Redistributions of source code must retain the above copyright notice,
        this list of conditions and the following disclaimer.

        Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

        Neither the name of the Astropy Team nor the names of its contributors
        may be used to endorse or promote products derived from this software
        without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    --------------------------------------------------------------------------

    Colormap often used in Milagro collaboration publications.

    This colormap is gray below ``vtransition`` and similar to the jet
    colormap above.

    A sharp gray -> color transition is often used for significance images
    with a transition value of ``vtransition ~ 5`` or ``vtransition ~ 7``,
    so that the following effect is achieved:

    - gray: non-significant features are not well visible
    - color: significant features at the detection threshold ``vmid``

    Note that this colormap is often critizised for over-exaggerating small
    differences in significance below and above the gray - color transition
    threshold.

    Parameters
    ----------
    vmin : float
        Minimum value (color: black)
    vmax : float
        Maximum value
    vtransition : float
        Transition value (below: gray, above: color).
    width : float
        Width of the transition
    huestart : float
        Hue of the color at ``vtransition``

    Returns
    -------
    colormap : `matplotlib.colors.LinearSegmentedColormap`
        Colormap

    Examples
    --------
    >>> from gammapy.image import colormap_milagro
    >>> vmin, vmax, vtransition = -5, 15, 5
    >>> cmap = colormap_milagro(vmin=vmin, vmax=vmax, vtransition=vtransition)

    """

    from colorsys import hls_to_rgb

    from matplotlib.colors import LinearSegmentedColormap

    if not vtransition > vmin:
        raise ValueError("vtransition must be greater than vmin.")

    # Compute normalised red, blue, yellow values
    transition = float(vtransition - vmin) / (vmax - vmin)

    # Create custom colormap
    # List entries: (value, (H, L, S))
    colors = [(0, (0.5, 0, 0)),
              (transition - width, (0.5, 0, 0)),
              (transition, (huestart, 0.4, 0.5)),
              (transition + width, (huestart, 0.4, 1)),
              # (0.99, (0, 0.6, 1)),
              (1, (0, 0.5, 0.5)),
              ]
    # Convert HLS values to RGB values
    rgb_colors = [(val, hls_to_rgb(*hls)) for (val, hls) in colors]
    cmap = LinearSegmentedColormap.from_list(name='milagro', colors=rgb_colors)

    return cmap
