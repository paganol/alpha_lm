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

    P%do_cross=parse_lgt(handle,'do_cross',.false.)

    P%inclfile=parse_string(handle,'cl_file','inputs/cls.dat') !order TT EE BB TE

    P%inbeamfile1=parse_string(handle,'beam_file1','inputs/beam1.fits')
    if (P%do_cross) P%inbeamfile2=parse_string(handle,'beam_file1','inputs/beam1.fits')

    P%outsigmafile1=parse_string(handle,'output_sigma1','outputs/sigma1.dat')
    if (P%do_cross) P%outsigmafile2=parse_string(handle,'output_sigma2','outputs/sigma2.dat')

    P%innoisefile1=parse_string(handle,'noise_file1','inputs/noise1.dat')
    if (len(P%innoisefile1) .gt. 0) then
       P%noiseE1=parse_double(handle,'noise_E1',2.0_dp)
       P%noiseB1=parse_double(handle,'noise_B1',P%noiseE1)
    endif

    if (P%do_cross) then 
       P%innoisefile2=parse_string(handle,'noise_file2',P%innoisefile1)
       if (len(P%innoisefile2) .gt. 0) then
          P%noiseE2=parse_double(handle,'noise_E2',P%noiseE1)
          P%noiseB2=parse_double(handle,'noise_B2',P%noiseE2)
       endif
    endif

    P%compute_biasalpha=parse_lgt(handle,'compute_alpha_bias',.true.)
    if (P%compute_biasalpha) then
       P%subtract_bias=parse_lgt(handle,'subtract_bias',.true.)
       P%outbiasfile=parse_string(handle,'output_bias','outputs/bias.txt')
       P%endnamebias=parse_string(handle,'suffix_bias','.txt')
    endif

    if (P%compute_alphalm) then
       P%inmapfile1=parse_string(handle,'input_map1','inputs/map1.fits')
       P%outalmfile1=parse_string(handle,'output_alm1','outputs/almhat1.fits') 
       P%inmaskfile1=parse_string(handle,'input_mask1','')
       if (P%do_cross) then
          P%inmapfile2=parse_string(handle,'input_map2','inputs/map2.fits')
          P%outalmfile2=parse_string(handle,'output_alm2','outputs/almhat2.fits')
          P%inmaskfile2=parse_string(handle,'input_mask2',P%inmaskfile1)
       endif
       P%nsims=parse_int(handle,'n_sims',1)
       P%ssim=parse_int(handle,'first_sim',1)
       if (P%nsims .gt. 1) then
          P%endnamemap1=parse_string(handle,'suffix_map1','.fits')
          P%endnamealm1=parse_string(handle,'suffix_alm1','.fits')
          if (P%do_cross) then
             P%endnamemap2=parse_string(handle,'suffix_map2','.fits')
             P%endnamealm2=parse_string(handle,'suffix_alm2','.fits')
          endif
       endif
       P%compute_alphacl=parse_lgt(handle,'compute_alpha_cl',.true.)
       if (P%compute_alphacl) then
           P%outclfile=parse_string(handle,'output_cl','outputs/clhat.fits')
           P%endnamecl=parse_string(handle,'suffix_cl','.txt')
       endif
    endif

    if (P%compute_biasalpha .or. P%compute_alphalm) then
       P%zerofill=parse_int(handle,'zero_fill',4)
       P%niter=parse_int(handle,'number_of_iterations',3)
    endif
     
    call parse_finish (handle)
    
    if (P%feedback .gt. 1 ) write(*,*) 'Read input file: DONE'
    

  end subroutine read_parameter_file
  
  
end module driver
