import numpy as np
import healpy as hp

rootcode='/marconi_work/INF22_lspe/lpagano0/test_alpha_lm/alpha_lm/'
rootglobal='/marconi_work/INF22_lspe/lpagano0/test_alpha_lm/runs/'

rootindir = rootglobal+'inputs/'
rootinmaps = '/marconi_scratch/userexternal/lpagano0/commander/' 
rootoutdir = rootglobal+'outputs/'
rootparamdir = rootglobal+'params/'
rootslurmdir = rootglobal+'slurms/'
codepath = rootcode+'EB_estimator'

read_alms = True

nside = 2048
fwhm = 5.0

nsimsnoise = 300

nsims = 301
start_sim = 0

lmin = 50
lmax = 750

Lmin = 0
Lmax = 500

namerun = 'commander_ns2048_Lmin'+str(Lmin).zfill(4)+'_Lmax'+str(Lmax).zfill(4)+'_lmin'+str(lmin).zfill(4)+'_lmax'+str(lmax).zfill(4)+'_maskstd'

lmax_beam = lmax

#files
beamfile = rootindir+'beam_'+namerun+'.fits'
fiducialfile = rootindir+'ffp10_lensedCls.dat'
paramfile = rootparamdir+'params_'+namerun+'.ini'
slurmfile = rootslurmdir+'slurm_'+namerun+'.sl' 

mask_file = '/marconi_work/INF22_indark/mbortola/masks/COM_Mask_CMB-common-HM-Misspix-Mask-Pol_2048_R3.00.fits'

if read_alms:   
    mapfile1 = rootinmaps+'dx12_v3_commander_cmb_noise_hm1_mc_alm_lmax'+str(lmax).zfill(4)+'_maskstd_'
    mapfile2 = rootinmaps+'dx12_v3_commander_cmb_noise_hm2_mc_alm_lmax'+str(lmax).zfill(4)+'_maskstd_'
    par_read_alms = 'T'
else:
    mapfile1 = rootinmaps+'dx12_v3_commander_cmb_noise_hm1_mc_'
    mapfile2 = rootinmaps+'dx12_v3_commander_cmb_noise_hm2_mc_'
    par_read_alms = 'F'

suffix_map1 = '.fits'
suffix_map2 = '.fits'
 
almfile1 = '!'+rootoutdir+'alphalm_'+namerun+'_hm1_'
almfile2 = '!'+rootoutdir+'alphalm_'+namerun+'_hm2_'

clfile = rootoutdir+'alphacl_'+namerun+'_'
biasfile = rootoutdir+'alphabias_'+namerun+'_'

sigmafile1 = rootoutdir+'sigma_'+namerun+'_hm1.txt'
sigmafile2 = rootoutdir+'sigma_'+namerun+'_hm2.txt'

noisefile1 = rootindir+'noise_'+namerun+'_hm1.txt'
noisefile2 = rootindir+'noise_'+namerun+'_hm2.txt'


#beam
bl = hp.gauss_beam(np.deg2rad(fwhm/60.),lmax_beam,pol=True)
wpix = hp.pixwin(nside,pol=True,lmax=lmax_beam)

beam = np.empty((lmax_beam+1,3))

beam[:,0] = bl[:,0]*wpix[0]
beam[:,1] = bl[:,1]*wpix[1]
beam[:,2] = bl[:,2]*wpix[1]

hp.write_cl(beamfile,beam.T,overwrite=True)

nl = np.zeros((lmax+1,2))
for i in range(nsimsnoise):
    d = np.load('/marconi_work/INF22_indark/mbortola/spectra/spectra_anafast/PR3_mc_commander_onlyhm1noise_autospectra/spectra_nside2048_lmax1500_'+str(i).zfill(5)+'.npy')
    nl[:,0] += d[0:lmax+1,2]
    nl[:,1] += d[0:lmax+1,3]

nl /= nsimsnoise

#nl[:,0]*=beam[0:lmax+1,1]*beam[0:lmax+1,1] 
#nl[:,1]*=beam[0:lmax+1,2]*beam[0:lmax+1,2]

np.savetxt(noisefile1,np.column_stack([np.arange(lmax+1),nl]))

nl = np.zeros((lmax+1,2))
for i in range(nsimsnoise):
    d = np.load('/marconi_work/INF22_indark/mbortola/spectra/spectra_anafast/PR3_mc_commander_onlyhm2noise_autospectra/spectra_nside2048_lmax1500_'+str(i).zfill(5)+'.npy')
    nl[:,0] += d[0:lmax+1,2]
    nl[:,1] += d[0:lmax+1,3]

nl /= nsimsnoise

#nl[:,0]*=beam[0:lmax+1,1]*beam[0:lmax+1,1]
#nl[:,1]*=beam[0:lmax+1,2]*beam[0:lmax+1,2]

np.savetxt(noisefile2,np.column_stack([np.arange(lmax+1),nl]))


params = """
feedback = 4

do_cross = T

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
beam_file2 = {beamfile}
noise_file1 = {noisefile1}
noise_file2 = {noisefile2}

input_mask1 = {mask_file}
input_mask2 = {mask_file}

read_precomputed_alms = {par_read_alms}

input_map1 = {mapfile1}
input_map2 = {mapfile2}
suffix_map1 = {suffix_map1}
suffix_map2 = {suffix_map2}

n_sims = {nsims}
first_sim = {start_sim}

zero_fill = 5


output_sigma1 = {sigmafile1}
output_sigma2 = {sigmafile2}
output_alm1 = {almfile1}
output_alm2 = {almfile2}
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
#SBATCH --mem=182000
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
