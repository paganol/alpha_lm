module settings
  
  use healpix_types
  use rngmod
  
  !initialize random handler
  type(planck_rng) :: rng_handle

  integer,parameter,public :: myTT=1,myEE=2,myBB=3,myTE=4,myTB=5,myEB=6

  !parameters
  Type Params
     integer :: feedback
     character(len=FILENAMELEN) :: inbeamfile, endnamemap,endnamealm,endnamecl,innoisefile
     integer :: Lmin,Lmax,ellmin,ellmax,zerofill
     integer :: elloffset
     integer :: nsims,ssim
     character(len=FILENAMELEN) :: inmapfile, inclfile, outalmfile,outsigmafile,outclfile
     character(len=FILENAMELEN) :: outbiasfile, endnamebias
     logical :: compute_alphalm, compute_alphacl, compute_biasalpha
     real(dp) :: noiseE,noiseB
  end Type Params

contains
  
end module settings
