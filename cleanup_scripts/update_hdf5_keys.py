import sys
from glob import glob
from pandas import HDFStore
import os

folder_path = sys.argv[1]

# All HDF5 files in the path
all_files = glob(os.path.join(folder_path, "*.h5"))

# Remove the PDF only ones (relevant results are in the PDF_KS and
# PDF_Hellinger keywords)
remove_keys = ["PDF"]

# Rename keys
rename_keys = {"VCS_Density": "VCS_Small_Scale",
               "VCS_Velocity": "VCS_Large_Scale"}

for f in all_files:
    store = HDFStore(f)

    # Removals
    for key in remove_keys:
        if "/"+key in store.keys():
            del store[key]

    # Rename
    for old_key in rename_keys:
        if "/"+old_key in store.keys():
            store[rename_keys[old_key]] = store[old_key].copy()
            del store[old_key]

    print("Final keys: " + str(store.keys()))

    store.close()
