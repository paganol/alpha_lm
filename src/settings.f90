module settings
  
  use healpix_types
  use rngmod
  
  !initialize random handler
  type(planck_rng) :: rng_handle

  integer,parameter,public :: myTT=1,myEE=2,myBB=3,myTE=4,myTB=5,myEB=6
  integer,parameter,public :: myT=1,myE=2,myB=3
  !parameters
  Type Params
     integer :: feedback
     character(len=FILENAMELEN) :: inbeamfile1,endnamemap1,endnamealm1,endnamecl,innoisefile1
     character(len=FILENAMELEN) :: inbeamfile2,endnamemap2,endnamealm2,innoisefile2
     character(len=FILENAMELEN) :: inmaskfile1,inmaskfile2
     integer :: Lmax,ellmin,ellmax,zerofill
     integer :: nsims,ssim,niter,nside
     integer :: nsims_mask
     character(len=FILENAMELEN) :: inmapfile1, inclfile, outalmfile1,outsigmafile1,outsigmafile2,outclfile
     character(len=FILENAMELEN) :: inmapfile2, outalmfile2
     character(len=FILENAMELEN) :: outbiasfile, endnamebias,outfskyfile
     logical :: compute_alphalm, compute_alphacl, compute_biasalpha, subtract_bias, do_cross, read_precomputed_alms
     logical :: compute_fskyl
     real(dp) :: noiseE1,noiseB1,noiseE2,noiseB2,ampsignal
  end Type Params

contains
  
end module settings
