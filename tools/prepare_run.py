import numpy as np
import healpy as hp

rootcode='/marconi_work/INF22_lspe/lpagano0/test_alpha_lm/alpha_lm/'
rootglobal='/marconi_work/INF22_lspe/lpagano0/test_alpha_lm/runs/'

rootindir = rootglobal+'inputs/'
rootoutdir = rootglobal+'outputs/'
rootparamdir = rootglobal+'params/'
rootslurmdir = rootglobal+'slurms/'
codepath = rootcode+'EB_estimator'


namerun = 'litebirdlikerot'#'wmaplike_Giorgia' #'plancklike_Giorgia' #'cmbpollike_Glu' #'plancklike_Glu'

nside = 256
fwhm = 30.0 #5.0 #7.1 #arcmin
noiseT = 2.0 #61.4
noiseP = 2.0*np.sqrt(2.0) #589.0#40.0*np.sqrt(2.) #86.8  

nsims=100

rotate_maps = True
sd = 0.006*np.deg2rad(1)**2
save_alpha = True

lmin = 2
lmax = 512#1500

Lmin = 0
Lmax = 500

lmax_gen = lmax

lmax_beam = lmax+Lmax

#files
beamfile = rootindir+'beam_'+namerun+'.fits'
fiducialfile = rootindir+'ffp10_lensedCls.dat'
paramfile = rootparamdir+'params_'+namerun+'.ini'
slurmfile = rootslurmdir+'slurm_'+namerun+'.sl' 

if nsims < 2:
    mapfile = rootindir+'map_'+namerun+'.fits'
    almfile = '!'+rootoutdir+'alphalm_'+namerun+'.fits'
    clfile = rootoutdir+'alphacl_'+namerun+'.txt'
    biasfile = rootoutdir+'alphabias_'+namerun+'.txt'
    if(save_alpha):
        alphafile = rootindir+'alphamap_'+namerun+'.fits'
else:
    mapfile = rootindir+'map_'+namerun+'_'
    almfile = '!'+rootoutdir+'alphalm_'+namerun+'_'
    clfile = rootoutdir+'alphacl_'+namerun+'_'
    biasfile = rootoutdir+'alphabias_'+namerun+'_'
    if(save_alpha):
        alphafile = rootindir+'alphamap_'+namerun+'_'

sigmafile = rootoutdir+'sigma_'+namerun+'.txt'

#beam
bl = hp.gauss_beam(np.deg2rad(fwhm/60.),lmax_beam,pol=True)
wpix = hp.pixwin(nside,pol=True,lmax=lmax_beam)

beam = np.empty((lmax_beam+1,3))

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

cl[0,2:lmax_gen+1]=TT[0:lmax_gen-1]/ll[0:lmax_gen-1]*beam[2:lmax_gen+1,0]*beam[2:lmax_gen+1,0]
cl[1,2:lmax_gen+1]=EE[0:lmax_gen-1]/ll[0:lmax_gen-1]*beam[2:lmax_gen+1,1]*beam[2:lmax_gen+1,1]
cl[2,2:lmax_gen+1]=BB[0:lmax_gen-1]/ll[0:lmax_gen-1]*beam[2:lmax_gen+1,2]*beam[2:lmax_gen+1,2]
cl[3,2:lmax_gen+1]=TE[0:lmax_gen-1]/ll[0:lmax_gen-1]*beam[2:lmax_gen+1,0]*beam[2:lmax_gen+1,1]

nl=np.zeros((4,lmax_gen+1))
nl[0,:]=NlTT
nl[1,:]=NlEE
nl[2,:]=NlBB

if rotate_maps:
    lalpha = np.arange(lmax_gen+1)
    llalpha = lalpha*(lalpha+1)/2/np.pi
    clalpha = sd/llalpha
    clalpha[0] = 0.0

if nsims < 2:
    m=hp.synfast(cl,nside,lmax=lmax_gen,pol=True,pixwin=False,new=True,verbose=False)
    noise=hp.synfast(nl,nside,lmax=lmax_gen,pol=True,pixwin=False,new=True,verbose=False)
    if rotate_maps:
        alpha=hp.synfast(clalpha,nside,lmax=lmax_gen,pol=False,pixwin=False,verbose=False)
        if save_alpha:
            hp.write_map(alphafile,alpha,overwrite=True)
        cos_a = np.cos(2*alpha)
        sin_a = np.sin(2*alpha)      
        m_R = np.zeros_like(m)
        m_R[0] = m[0]
        m_R[1] =  m[1]*cos_a+m[2]*sin_a
        m_R[2] = -m[1]*sin_a+m[2]*cos_a
        hp.write_map(mapfile,m_R+noise,overwrite=True)
    else:
        hp.write_map(mapfile,m+noise,overwrite=True) 
else:
    for isim in np.arange(nsims):
        print(isim)
        m=hp.synfast(cl,nside,lmax=lmax_gen,pol=True,pixwin=False,new=True,verbose=False)
        noise=hp.synfast(nl,nside,lmax=lmax_gen,pol=True,pixwin=False,new=True,verbose=False)
        if rotate_maps:
            alpha=hp.synfast(clalpha,nside,lmax=lmax_gen,pol=False,pixwin=False,verbose=False)
            if save_alpha:
                hp.write_map(alphafile+str(isim).zfill(4)+'.fits',alpha,overwrite=True)
            cos_a = np.cos(2*alpha)
            sin_a = np.sin(2*alpha)
            m_R = np.zeros_like(m)
            m_R[0] = m[0]
            m_R[1] =  m[1]*cos_a+m[2]*sin_a
            m_R[2] = -m[1]*sin_a+m[2]*cos_a
            hp.write_map(mapfile+str(isim).zfill(4)+'.fits',m_R+noise,overwrite=True)
        else:
            hp.write_map(mapfile+str(isim).zfill(4)+'.fits',m+noise,overwrite=True)



params = """
feedback = 4

do_cross = F

compute_alpha_lm = T
compute_alpha_cl = T
compute_alpha_bias = T
subtract_bias = F

number_of_iterations = 3

ellmin = {lmin}
ellmax = {lmax}

Lmin = {Lmin} 
Lmax = {Lmax}

cl_file = {fiducialfile}
beam_file1 = {beamfile}
noise_file1 = '' 
noise_E1 = {noiseP}

input_map1 = {mapfile}

n_sims = {nsims}
first_sim = 0

output_sigma1 = {sigmafile}
output_alm1 = {almfile}
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
