import numpy as np
import healpy as hp

rootglobal='/marconi_work/INF22_lspe/lpagano0/test_alpha_lm/alpha_lm/'

rootindir = rootglobal+'inputs/'
rootoutdir = rootglobal+'outputs/'
rootparamdir = rootglobal+'params/'
rootslurmdir = rootglobal+'slurms/'
codepath = rootglobal+'EB_estimator'


namerun = 'wmaplike_Giorgia' #'plancklike_Giorgia' #'cmbpollike_Glu' #'plancklike_Glu'

nside = 512
lmax_gen = 3*nside
fwhm = 21.0 #5.0 #7.1 #arcmin
noiseT = 40.0 #61.4
noiseP = 894.6 #589.0#40.0*np.sqrt(2.) #86.8  

nsims=1

lmin = 2
lmax = 1024#1500

Lmin = 0
Lmax = 500

#files
beamfile = rootindir+'beam_'+namerun+'.fits'
fiducialfile = rootindir+'ffp10_lensedCls.dat'
paramfile = rootparamdir+'params_'+namerun+'.ini'
slurmfile = rootslurmdir+'slurm_'+namerun+'.sl' 

if nsims < 2:
    mapfile = rootindir+'map_'+namerun+'.fits'
    almfile = rootoutdir+'alphalm_'+namerun+'.fits'
    clfile = rootoutdir+'alphacl_'+namerun+'.txt'
    biasfile = rootoutdir+'alphabias_'+namerun+'.txt'
else:
    mapfile = rootindir+'map_'+namerun+'_'
    almfile = rootoutdir+'alphalm_'+namerun+'_'
    clfile = rootoutdir+'alphacl_'+namerun+'_'
    biasfile = rootoutdir+'alphabias_'+namerun+'_'

sigmafile = rootoutdir+'sigma_'+namerun+'.txt'

#beam
bl = hp.gauss_beam(np.deg2rad(fwhm/60.),lmax_gen,pol=True)
wpix = hp.pixwin(nside,pol=True,lmax=lmax_gen)

beam = np.empty((lmax_gen+1,3))

beam[:,0] = bl[:,0]*wpix[0]
beam[:,1] = bl[:,1]*wpix[1]
beam[:,2] = bl[:,2]*wpix[1]

hp.write_cl(beamfile,beam.T,overwrite=True)

#fiducial
TT,EE,BB,TE=np.loadtxt(fiducialfile,unpack=True,usecols=(1,2,3,4))

cl=np.zeros((4,lmax_gen+1))

l=np.arange(lmax_gen-1)+2
ll=l*(l+1)/2/np.pi

NlTT = (noiseT*np.pi/180./60.)**2
NlEE = (noiseP*np.pi/180./60.)**2
NlBB = (noiseP*np.pi/180./60.)**2

cl[0,2:lmax_gen+1]=TT[0:lmax_gen-1]/ll[0:lmax_gen-1]*beam[2:lmax_gen+1,0]*beam[2:lmax_gen+1,0]+NlTT
cl[1,2:lmax_gen+1]=EE[0:lmax_gen-1]/ll[0:lmax_gen-1]*beam[2:lmax_gen+1,1]*beam[2:lmax_gen+1,1]+NlEE
cl[2,2:lmax_gen+1]=BB[0:lmax_gen-1]/ll[0:lmax_gen-1]*beam[2:lmax_gen+1,2]*beam[2:lmax_gen+1,2]+NlBB
cl[3,2:lmax_gen+1]=TE[0:lmax_gen-1]/ll[0:lmax_gen-1]*beam[2:lmax_gen+1,0]*beam[2:lmax_gen+1,1]

if nsims < 2:
    m=hp.synfast(cl,nside,lmax=lmax_gen,pol=True,pixwin=False,new=True,verbose=False)
    hp.write_map(mapfile,m,overwrite=True)
else:
    for isim in np.arange(nsims):
        print(isim)
        m=hp.synfast(cl,nside,lmax=lmax_gen,pol=True,pixwin=False,new=True,verbose=False)
        hp.write_map(mapfile+str(isim).zfill(4)+'.fits',m,overwrite=True)

params = """
feedback = 4

compute_alpha_lm = F
compute_alpha_cl = F
compute_alpha_bias = F
subtract_bias = F

ellmin = {lmin}
ellmax = {lmax}

Lmin = {Lmin} 
Lmax = {Lmax}

cl_file = {fiducialfile}
beam_file = {beamfile}
noise_file = '' 
noise_E = {noiseP}

map_file = {mapfile}

n_sims = {nsims}
first_sim = 1

output_sigma = {sigmafile}
output_alm = {almfile}
output_cl = {clfile}
output_bias = {biasfile}
"""
params = params.format(**locals())
f = open(paramfile, "wt")
f.write(params)
f.close()

slurm = """#!/bin/bash -l
#SBATCH -p skl_usr_dbg
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH -t 0:30:00
#SBATCH -J run_job
#SBATCH -o run_job.log
#SBATCH -A INF22_lspe
#SBATCH --export=ALL
#SBATCH --mem=86000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pgnlcu@unife.it

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1
srun {codepath} {paramfile} 

"""
slurm = slurm.format(**locals())
f = open(slurmfile, "wt")
f.write(slurm)
f.close()
