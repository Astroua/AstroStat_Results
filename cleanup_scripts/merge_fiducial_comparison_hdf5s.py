
import sys
from glob import glob
from pandas import HDFStore
import os

folder_path = sys.argv[1]

faces = ["_0_0_", "_0_2_", "_2_0_", "_2_2_"]

for face in faces:
    old_comp = glob(os.path.join(folder_path,
                                 "*_comparisons_*" + face + "*.h5"))
    new_comp = glob(os.path.join(folder_path,
                                 "*8_fiducialfid_comp" + face + "*.h5"))

    print(old_comp)
    print(new_comp)
    assert len(old_comp) == 1
    assert len(new_comp) == 1

    old_result = HDFStore(old_comp[0])
    new_result = HDFStore(new_comp[0])

    for key in old_result.keys():
        if key in new_result.keys():
            continue
        new_result[key] = old_result[key].copy()

    print("New file keys: " + str(new_result.keys()))

    old_result.close()
    new_result.close()
