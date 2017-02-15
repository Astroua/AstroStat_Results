
'''
Constrain the effects of viewing angle.

Requires comparisons between all sets of fiducial faces.
'''

import sys
import os
from pandas import DataFrame
import numpy as np
from multiprocessing import Manager
from datetime import datetime

from wrapping_function import stats_wrapper
from turbustat.statistics import statistics_list

from analysis_funcs import load_and_reduce, files_sorter

# Enable multiprocessing
multiprocess = True
if multiprocess:
    manager = Manager()

path_to_data = sys.argv[1]
face1 = int(sys.argv[2])
face2 = int(sys.argv[3])
results_dir = sys.argv[4]

if not os.path.exists(path_to_data):
    raise IOError("path_to_data ({}) does not exist.".format(path_to_data))

if not os.path.exists(results_dir):
    raise IOError("results_dir ({}) does not exist.".format(results_dir))

# Face 1 and Face 2 are the comparison to perform. There are 3 unique
# combinations
comp = "{0}_{1}".format(face1, face2)

stats = statistics_list
stats.remove("Tsallis")

if comp not in ["0_1", "0_2", "1_2"]:
    raise Warning("The unique combinations are 0_1, 0_2, 1_2. For the analysis"
                  " pipeline, only these choices are allowed.")

# Find and organize filenames
fiducials, designs = files_sorter(path_to_data, timesteps='max',
                                  append_prefix=True)[:2]


def runner(args):
    name1, name2, comp, stats_dict = args

    # print("On {0} {1}".format(name1, name2))

    dataset1 = load_and_reduce(name1)
    dataset2 = load_and_reduce(name2)

    output = stats_wrapper(dataset1, dataset2, statistics=stats)[0]

    for stat in output:
        stats_dict[stat].append([output[stat], os.path.basename(name1),
                                 os.path.basename(name2)])


print("Starting pool at {}".format(datetime.now()))

if multiprocess:

    use_mpi = False
    if use_mpi:
        from mpipool import MPIPool
        pool = MPIPool(loadbalance=False)

        if not pool.is_master():
            # Wait for instructions from the master process.
            pool.wait()
            sys.exit(0)
    else:
        from multiprocessing import cpu_count, Pool
        psize = cpu_count()
        print("Found {} CPUs to run on.".format(psize))
        pool = Pool(processes=psize)

print("Pool created at {}".format(datetime.now()))

# Run distance between comparisons

stats_dict = dict.fromkeys(stats)
for stat in stats_dict:
    if multiprocess:
        stats_dict[stat] = manager.list([])
    else:
        stats_dict[stat] = []

for fid in fiducials[int(comp[0])]:

    print("On fiducial {0} of {1}".format(fid, len(fiducials[int(comp[0])])))
    print(str(datetime.now()))

    if multiprocess:
        iterat = ((fiducials[int(comp[0])][fid][i],
                   fiducials[int(comp[-1])][fid][i], comp, stats_dict)
                  for i in xrange(len(fiducials[int(comp[0])][fid])))

        pool.map(runner, iterat)

    else:
        for i in range(len(fiducials[int(comp[0])][fid])):
            name1 = fiducials[int(comp[0])][fid][i]
            name2 = fiducials[int(comp[-1])][fid][i]
            print("On {0} {1}".format(name1, name2))

            dataset1 = load_and_reduce(name1)
            dataset2 = load_and_reduce(name2)

            output = stats_wrapper(dataset1, dataset2, statistics=stats)[0]

            for stat in output:
                stats_dict[stat].append(output[stat], os.path.basename(name1),
                                        os.path.basename(name2))

print("Finished fiducial at {}".format(datetime.now()))

for des in designs[int(comp[0])]:

    print("On design {0} of {1}".format(des, len(designs[int(comp[0])])))
    print(str(datetime.now()))

    if multiprocess:
        iterat = ((designs[int(comp[0])][des][i],
                   designs[int(comp[-1])][des][i], comp, stats_dict)
                  for i in xrange(len(designs[int(comp[0])][des])))

        pool.map(runner, iterat)

    else:
        for i in range(len(designs[int(comp[0])][des])):
            name1 = designs[int(comp[0])][des][i]
            name2 = designs[int(comp[-1])][des][i]
            print("On {0} {1}".format(name1, name2))

            dataset1 = load_and_reduce(name1)
            dataset2 = load_and_reduce(name2)

            output = stats_wrapper(dataset1, dataset2, statistics=stats)[0]

            for stat in output:
                stats_dict[stat].append(output[stat], os.path.basename(name1),
                                        os.path.basename(name2))

print("Finished designs at {}".format(datetime.now()))

# Convert to normal lists
if multiprocess:
    for stat in stats_dict:
        stats_dict[stat] = list(stats_dict[stat])

# Add the filename columns. The list order should be preserved since we're not
# using asynchronous mapping
for i, stat in enumerate(stats_dict.keys()):
    results_arr = np.array(stats_dict[stat]).T
    if i == 0:
        stats_dict["Data 1"] = results_arr[1]
        stats_dict["Data 2"] = results_arr[2]
    stats_dict[stat] = results_arr[0].astype(np.float)

# Now save the results
df = DataFrame(stats_dict)
df.to_csv(os.path.join(results_dir,
                       "view_angle_comparison_{}.csv".format(comp)))

if multiprocess:
    pool.close()

print("Done at {}".format(datetime.now()))
