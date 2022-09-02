module driver

  use paramfile_io 
  use settings
  
  
contains
  
  subroutine read_parameter_file(P)
    Type(Params) :: P
    integer :: argc
    type(paramfile_handle) :: handle
    character(len=FILENAMELEN) :: parfile
    
    
    argc = command_argument_count()     
    if (argc /= 1) then
       stop 'Run ./code <param.ini>'
    else
       CALL get_command_argument(1,parfile)
    endif
    
    handle = parse_init(trim(parfile),.false.)
    
    P%feedback=parse_int(handle,'feedback',1)
    if (P%feedback .gt. 3 ) write(*,*) 'Read input file: STARTING'

    P%compute_alphalm=parse_lgt(handle,'compute_alpha_lm',.true.)

    P%ellmin=parse_int(handle,'ellmin',2)
    P%ellmax=parse_int(handle,'ellmax',100)
    P%Lmin=parse_int(handle,'Lmin',1)
    P%Lmax=parse_int(handle,'Lmax',100)

    P%elloffset=parse_int(handle,'ell_offset',P%Lmax+1)

    P%inclfile=parse_string(handle,'cl_file','inputs/cls.dat') !order TT EE BB TE

    P%inbeamfile=parse_string(handle,'beam_file','inputs/beam.fits')

    P%outsigmafile=parse_string(handle,'output_sigma','outputs/sigma.dat')

    P%innoisefile=parse_string(handle,'noise_file','inputs/noise.dat')
    if (len(P%innoisefile) .gt. 0) then
       P%noiseE=parse_double(handle,'noise_E',2.0_dp)
       P%noiseB=parse_double(handle,'noise_B',P%noiseE)
    endif

    if (P%compute_alphalm) then
       P%inmapfile=parse_string(handle,'map_file','inputs/map.fits')
       P%outalmfile=parse_string(handle,'output_alm','outputs/almhat.fits') 
       P%nsims=parse_int(handle,'n_sims',1)
       P%ssim=parse_int(handle,'first_sim',1)
       if (P%nsims .gt. 1) then
          P%zerofill=parse_int(handle,'zero_fill',4)
          P%endnamemap=parse_string(handle,'suffix_map','.fits')
          P%endnamealm=parse_string(handle,'suffix_alm','.fits')
       endif
       P%compute_alphacl=parse_lgt(handle,'compute_alpha_cl',.true.)
       if (P%compute_alphacl) then
           P%outclfile=parse_string(handle,'output_cl','outputs/clhat.fits')
           P%endnamecl=parse_string(handle,'suffix_cl','.txt')
       endif
       P%compute_biasalpha=parse_lgt(handle,'compute_alpha_bias',.true.)
       if (P%compute_biasalpha) then
          P%subtract_bias=parse_lgt(handle,'subtract_bias',.true.)
          P%outbiasfile=parse_string(handle,'output_bias','outputs/bias.txt')
          P%endnamebias=parse_string(handle,'suffix_bias','.txt')
       endif
    endif
      
    call parse_finish (handle)
    
    if (P%feedback .gt. 1 ) write(*,*) 'Read input file: DONE'
    

  end subroutine read_parameter_file
  
  
end module driver
