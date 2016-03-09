
'''
Remove a key from an HDF5 file.
'''

from pandas import HDFStore
import os
import sys
import glob


path = sys.argv[1]
key = sys.argv[2]

# Get HDF5 files in the given path
hdf5_files = glob.glob(os.path.join(path, "*.h5"))

for hfile in hdf5_files:
    print(hfile)
    store = HDFStore(hfile)

    if key in store:
        del store[key]

    store.close()
