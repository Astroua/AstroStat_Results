
import numpy as np
import matplotlib.pyplot as plt
import yt
import matplotlib as mpl

'''
Create the turbulent KE power spectrum to investigate the inertial range.

Code is from: http://yt-project.org/doc/cookbook/power_spectrum_example.py
with modifications by @low-sky

'''

plt.ion()

mpl.rcParams['font.family']='serif'
mpl.rcParams['font.serif'] = 'FreeSerif'


def makefig():

    fig = plt.figure(figsize=(4.5,4))
    ds = yt.load('Fiducial0/DD0020/DD0020')
    doit(ds,label=r'Fiducial 1 $128^3$')
    ds3 = yt.load('Fiducial1/DD0020/DD0020')
    doit(ds3, lw=3,alpha=0.5, label='Fiducial 2 $128^3$')
    ds2 = yt.load('Fiducial_256/DD0020/DD0020')
    doit(ds2,linestyle='dashed',label=r'Fiducial $256^3$')

    plt.ylim(0.06,2)
    plt.xlim(1,200)

    plt.fill_between([2,8],0.06,2,alpha=0.5,facecolor='gray')
    plt.axvline(1.0)
    plt.legend(loc=3)
    plt.tight_layout()
    # plt.savefig('powerspec.pdf')
    # plt.clf()

def doit(ds, **kwargs):

    # a FFT operates on uniformly gridded data.  We'll use the yt
    # covering grid for this.

    max_level = ds.index.max_level

    ref = 1 #int(np.product(ds.ref_factors[0:max_level]))

    low = ds.domain_left_edge
    dims = ds.domain_dimensions*ref

    print("dims: {}".format(dims))

    nx, ny, nz = dims

    nindex_rho = 1./3.

    Kk = np.zeros( (nx/2+1, ny/2+1, nz/2+1))

    for vel in [("gas", "velocity_x"), ("gas", "velocity_y"),
                ("gas", "velocity_z")]:

        Kk += 0.5*fft_comp(ds, ("gas", "density"), vel,
                           nindex_rho, max_level, low, dims)

    # wavenumbers
    L = (ds.domain_right_edge - ds.domain_left_edge).d

    print("L: {}".format(L))

    kx = np.fft.rfftfreq(nx)*nx/L[0]
    ky = np.fft.rfftfreq(ny)*ny/L[1]
    kz = np.fft.rfftfreq(nz)*nz/L[2]

    # physical limits to the wavenumbers
    kmin = np.min(1.0/L)
    kmax = np.min(0.5*dims/L)

    print("kmin: {}".format(kmin))
    print("kmax: {}".format(kmax))

    kbins = np.arange(kmin, kmax, kmin)
    N = len(kbins)

    # bin the Fourier KE into radial kbins
    kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
    k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

    whichbin = np.digitize(k.flat, kbins)
    ncount = np.bincount(whichbin)

    E_spectrum = np.zeros(len(ncount)-1)

    for n in range(1,len(ncount)):
        E_spectrum[n-1] = np.sum(Kk.flat[whichbin==n])


    k = 0.5*(kbins[0:N-1] + kbins[1:N])
    E_spectrum = E_spectrum[1:N]
    E_spectrum = E_spectrum[1:-1]
    k=k[1:-1]  # Drop last bin.  It's wiggy.
    index = np.argmax(E_spectrum)
    kmax = k[index]
    Emax = E_spectrum[index]

    plt.loglog(k, E_spectrum/Emax*(k/kmax)**(5/3.),**kwargs)
    # plt.loglog(k, Emax*(k/kmax)**(-5./3.), ls=":", color="0.5")
    # plt.axhline(Emax,alpha=0.5,lw=4)

    plt.xlabel(r"$k$")
    plt.ylabel(r"$E(k) k^{5/3}\, dk$")

#    plt.savefig("spectrum.png")
#    plt.clf()

def fft_comp(ds, irho, iu, nindex_rho, level, low, delta ):

    cube = ds.covering_grid(level, left_edge=low,
                            dims=delta,
                            fields=[irho, iu])

    rho = cube[irho].d
    u = cube[iu].d

    nx, ny, nz = rho.shape

    # do the FFTs -- note that since our data is real, there will be
    # too much information here.  fftn puts the positive freq terms in
    # the first half of the axes -- that's what we keep.  Our
    # normalization has an '8' to account for this clipping to one
    # octant.
    ru = np.fft.fftn(rho**nindex_rho * u)[0:nx/2+1,0:ny/2+1,0:nz/2+1]
    ru = 8.0*ru/(nx*ny*nz)

    return np.abs(ru)**2
