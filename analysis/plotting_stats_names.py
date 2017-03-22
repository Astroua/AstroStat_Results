
'''
A list of the statistics names to be used in the plots (namely map_all_results)
'''

from copy import copy

from turbustat.statistics import statistics_list

statistics = copy(statistics_list)
statistics.append("DeltaVariance_Centroid_Slope")
statistics.append("DeltaVariance_Centroid_Curve")

statistics_names = dict.fromkeys(statistics)

for stat in statistics_names:

    stat_spaced = stat.replace("_", " ")

    if "Dendrogram" in stat:
        statistics_names[stat] = \
            stat_spaced.replace("Dendrogram", "Dendro.") + "."
    elif "DeltaVariance" in stat:
        statistics_names[stat] = \
            stat_spaced.replace("DeltaVariance", "Del. Var.")
    elif stat == "PSpec":
        statistics_names[stat] = "Spatial Power Spec."
    else:
        statistics_names[stat] = stat_spaced
