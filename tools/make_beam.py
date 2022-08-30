import numpy as np
import healpy as hp

lmax = 2000
fwhm = 5 #arcmin
nside = 512 

bl = hp.gauss_beam(np.deg2rad(fwhm/60.),lmax,pol=True)
wpix = hp.pixwin(nside,pol=True,lmax=lmax)

beam = np.empty((lmax+1,3))

beam[:,0] = bl[:,0]*wpix[0]
beam[:,1] = bl[:,1]*wpix[1]
beam[:,2] = bl[:,2]*wpix[1]

hp.write_cl('../inputs/beam.fits',beam.T)
