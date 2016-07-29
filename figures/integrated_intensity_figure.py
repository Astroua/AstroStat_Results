
'''
Load up some moment 0 arrays and make a nice plot
'''

import glob
import os
from astropy.io import fits
from astropy.visualization import scale_image
import matplotlib.pyplot as p
import numpy as np


moment_direc = "/media/eric/Data_3/SimSuite8/moments/"
face = 0
face_select = "_0"+str(face)+"_"


timestep = 30

fiducials = []
for i in range(5):
    # Get all moment0 for that Fiducial

    fid_select = "Fiducial"+str(i)

    ntime = timestep

    while True:

        time_select = "00"+str(timestep)
        select = "*"+fid_select+"*"+time_select+"*"+face_select+"*moment0*"

        mom0s = glob.glob(os.path.join(moment_direc, select))

        if len(mom0s) != 0:
            break

        ntime -= ntime

    fiducials.append(fits.getdata(mom0s[0]))


designs = []
for i in range(32):
    # Get all moment0 for that Fiducial

    des_select = "Design"+str(i)+"_"

    ntime = timestep

    while True:

        time_select = "00"+str(ntime)
        select = "*"+des_select+"*"+time_select+"*"+face_select+"*moment0*"

        mom0s = glob.glob(os.path.join(moment_direc, select))

        if len(mom0s) != 0:
            break

        ntime -= 1

        if ntime < 21:
            raise Warning("Not working for Design "+str(i))

    designs.append(fits.getdata(mom0s[0]))

all_laststeps = fiducials + designs

# for mom0 in all_laststeps:

#     mom0_arr = fits.getdata(mom0)

#     p.imshow(mom0_arr, origin='lower')
#     p.title(mom0.split("/")[-1])
#     raw_input("Continue?")

vmin = 100.0
vmax = 0.0

viewed = [fiducials[1], fiducials[2], designs[2], designs[3], designs[6],
          designs[8], designs[22], designs[24]]

for arr in viewed:
    arr_max = np.nanmax(arr)
    arr_min = np.nanmin(arr)

    print arr_max, arr_min

    if arr_max > vmax:
        vmax = arr_max
    if arr_min < vmin:
        vmin = arr_min

for i in range(len(viewed)):
    viewed[i][np.isnan(viewed[i])] = 0.0


param_str = r"$\zeta \mathcal{M} \alpha \beta k$"

fig, axes = p.subplots(2, 4)

axes[0, 0].imshow(np.log10(viewed[0][:-1, :-1]), vmin=np.log10(vmin),
                  vmax=np.log10(vmax), origin='lower')
axes[0, 0].set_title("Fiducial")
axes[0, 0].set_xticks([])
axes[0, 0].set_yticks([])

axes[1, 0].imshow(np.log10(viewed[1][:-1, :-1]), vmin=np.log10(vmin),
                  vmax=np.log10(vmax), origin='lower')
axes[1, 0].set_title("Fiducial")
axes[1, 0].set_xticks([])
axes[1, 0].set_yticks([])

axes[0, 1].imshow(np.log10(viewed[2][:-1, :-1]), vmin=np.log10(vmin),
                  vmax=np.log10(vmax), origin='lower')
axes[0, 1].set_title(param_str + "\n01000")
axes[0, 1].set_xticks([])
axes[0, 1].set_yticks([])

axes[0, 2].imshow(np.log10(viewed[3][:-1, :-1]), vmin=np.log10(vmin),
                  vmax=np.log10(vmax), origin='lower')
axes[0, 2].set_title(param_str + "\n01001")
axes[0, 2].set_xticks([])
axes[0, 2].set_yticks([])

axes[0, 3].imshow(np.log10(viewed[4][:-1, :-1]), vmin=np.log10(vmin),
                  vmax=np.log10(vmax), origin='lower')
axes[0, 3].set_title(param_str + "\n01010")
axes[0, 3].set_xticks([])
axes[0, 3].set_yticks([])

axes[1, 1].imshow(np.log10(viewed[5][:-1, :-1]), vmin=np.log10(vmin),
                  vmax=np.log10(vmax), origin='lower')
axes[1, 1].set_title(param_str + "\n10000")
axes[1, 1].set_xticks([])
axes[1, 1].set_yticks([])

axes[1, 2].imshow(np.log10(viewed[6][:-1, :-1]), vmin=np.log10(vmin),
                  vmax=np.log10(vmax), origin='lower')
axes[1, 2].set_title(param_str + "\n01111")
axes[1, 2].set_xticks([])
axes[1, 2].set_yticks([])

axes[1, 3].imshow(np.log10(viewed[7][:-1, :-1]), vmin=np.log10(vmin),
                  vmax=np.log10(vmax), origin='lower')
axes[1, 3].set_title(param_str + "\n10100")
axes[1, 3].set_xticks([])
axes[1, 3].set_yticks([])
fig.subplots_adjust(hspace=0.025, wspace=0.025)
fig.show()
