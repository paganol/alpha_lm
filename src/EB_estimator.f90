program EB_estimator
  
  use mpi
  use alm_tools
  use pix_tools
  use fitstools
  use head_fits
  use settings
  use driver
  use utils
  use wigner 
!  use wigner1
  
  implicit none
  
  !variables
  integer :: myid,nproc,mpierr
  integer :: ind,indmax,indLMmax,usedellpmax,lm(2)
  integer :: iellp,ellpmin,ellpmax,iell
  integer :: isim,iL,iM,maxsimperproc,simproc
  integer(i8b) :: npix,ipix
  logical :: apply_mask1=.false.,apply_mask2=.false.
  real(dp), allocatable, dimension(:) :: one_o_var1,one_o_var2,wig2,fsky_l,fsky_l_den,alm_amp_fid_fsky_l
  real(dp), allocatable, dimension(:) :: F_EB1,F_BE1,clEEfid,clEEobs1,clBBobs1
  real(dp), allocatable, dimension(:) :: F_EB2,F_BE2,clEEobs2,clBBobs2
  real(dp), allocatable, dimension(:,:) :: clEEmap1,clBBmap1,clEEmap12,clBBmap12,biasalpha
  complex(dpc), allocatable, dimension(:,:) :: almE1,almB1,almE2,almB2
  complex(dpc), allocatable, dimension(:,:,:) :: alm1,alm2
  complex(spc), allocatable, dimension(:,:,:) :: alm
  complex(dpc), allocatable, dimension(:) :: almalpha1,almalpha2
  real(dp), allocatable,dimension(:,:) :: clfid,wl1,bl1,nl1,wl2,bl2,nl2
  real(sp), allocatable,dimension(:) :: maskt,map
  real(dp), allocatable,dimension(:,:) :: mask1,mask2
  real(sp), allocatable,dimension(:,:) :: numerator,denominator
  real(dp), allocatable,dimension(:,:) :: wpix
  real(dp) :: factor,Gl,fsky1,fsky2,amp,hsqrt2
  integer :: t0,t1,t2,t3,t4,clock_rate,clock_max,myunit,ct
  character(len=FILENAMELEN) :: mapname, filename
  character(len=16) :: simstr
  character(len=1) :: strzerofill
!  integer :: IER
  
  !input parameters
  Type(Params) :: Par
  
  !  !! ================================ Init ================================ !!
  call mpi_init(mpierr)
  call mpi_comm_size(mpi_comm_world, nproc, mpierr)
  call mpi_comm_rank(mpi_comm_world, myid, mpierr)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !! Master reads parameters and does some stuff: reads maps and computes alms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  if (myid .eq. 0) then 
     call read_parameter_file(Par)
     if (Par%feedback .gt. 0) write(0,*) 'Reading parameters'
     if (Par%feedback .gt. 2) call system_clock(t0,clock_rate,clock_max)

     allocate(bl1(0:Par%ellmax,3),wl1(0:Par%ellmax,6))
     call read_beam(Par%inbeamfile1,bl1,wl1)
     if (Par%do_cross) then
        allocate(bl2(0:Par%ellmax,3),wl2(0:Par%ellmax,6))
        call read_beam(Par%inbeamfile2,bl2,wl2)
     endif
     
     allocate(clfid(0:Par%ellmax,6))
     call read_cl(Par%inclfile,clfid)
     allocate(clEEfid(0:Par%ellmax))
     clEEfid = clfid(:,myEE)
     
     allocate(nl1(0:Par%ellmax,2))
     call make_noise(Par%innoisefile1,nl1,Par%noiseE1,Par%noiseB1,Par%feedback) 
     if (Par%do_cross) then     
        allocate(nl2(0:Par%ellmax,2))
        call make_noise(Par%innoisefile2,nl2,Par%noiseE2,Par%noiseB2,Par%feedback)
     endif
     
     allocate(clEEobs1(0:Par%ellmax))
     allocate(clBBobs1(0:Par%ellmax))
     clEEobs1 = clfid(:,myEE) * wl1(:,myEE) + nl1(:,1)
     clBBobs1 = clfid(:,myBB) * wl1(:,myBB) + nl1(:,2)
     deallocate(nl1,wl1)
 
     if (Par%do_cross) then
        allocate(clEEobs2(0:Par%ellmax))
        allocate(clBBobs2(0:Par%ellmax))
        clEEobs2 = clfid(:,myEE) * wl2(:,myEE) + nl2(:,1)
        clBBobs2 = clfid(:,myBB) * wl2(:,myBB) + nl2(:,2)
        deallocate(nl2,wl2)
     endif
   
     deallocate(clfid)
     
  endif
  
  call mpi_barrier(mpi_comm_world, mpierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !! Broadcast relevant parameters and variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  if ((myid .eq. 0) .and. (Par%feedback .gt. 1)) write(*,*) 'Broadcast parameters and variables'
  
  call mpi_bcast(Par%feedback,1,mpi_integer,0, mpi_comm_world, mpierr)
  call mpi_bcast(Par%Lmax,1,mpi_integer,0,mpi_comm_world, mpierr)
  call mpi_bcast(Par%ellmin,1,mpi_integer,0,mpi_comm_world, mpierr)
  call mpi_bcast(Par%ellmax,1,mpi_integer,0,mpi_comm_world, mpierr)
  call mpi_bcast(Par%compute_alphalm,1,mpi_logical,0,mpi_comm_world,mpierr)
  call mpi_bcast(Par%compute_biasalpha,1,mpi_logical,0,mpi_comm_world,mpierr)
  if (Par%compute_biasalpha) call mpi_bcast(Par%subtract_bias,1,mpi_logical,0,mpi_comm_world,mpierr)
  call mpi_bcast(Par%do_cross,1,mpi_logical,0,mpi_comm_world,mpierr)
  
  if (myid .ne. 0) then
     allocate(clEEfid(0:Par%ellmax),clEEobs1(0:Par%ellmax),clBBobs1(0:Par%ellmax))
     allocate(bl1(0:Par%ellmax,3))
  endif
  call mpi_bcast(clEEfid,Par%ellmax+1,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(clEEobs1,Par%ellmax+1,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(clBBobs1,Par%ellmax+1,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(bl1,3*(Par%ellmax+1),mpi_real8,0,mpi_comm_world, mpierr)

  if (Par%do_cross) then
     if (myid .ne. 0) then
        allocate(clEEobs2(0:Par%ellmax),clBBobs2(0:Par%ellmax))
        allocate(bl2(0:Par%ellmax,3))
     endif
     call mpi_bcast(clEEobs2,Par%ellmax+1,mpi_real8,0,mpi_comm_world, mpierr)
     call mpi_bcast(clBBobs2,Par%ellmax+1,mpi_real8,0,mpi_comm_world, mpierr)
     call mpi_bcast(bl2,3*(Par%ellmax+1),mpi_real8,0,mpi_comm_world, mpierr)
  endif
  
  call mpi_barrier(mpi_comm_world, mpierr)
  if (Par%compute_alphalm .or. Par%compute_biasalpha)  then
     call mpi_bcast(Par%nsims,1,mpi_integer,0,mpi_comm_world,mpierr)
     call mpi_bcast(Par%ssim,1,mpi_integer,0,mpi_comm_world,mpierr)
     call mpi_bcast(Par%inmapfile1,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
     call mpi_bcast(Par%endnamemap1,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
     call mpi_bcast(Par%zerofill,1,mpi_integer,0,mpi_comm_world,mpierr)
     call mpi_bcast(Par%niter,1,mpi_integer,0,mpi_comm_world,mpierr)
     call mpi_bcast(Par%compute_fskyl,1,mpi_logical,0,mpi_comm_world,mpierr)
     if (Par%compute_alphalm) then
        call mpi_bcast(Par%nside,1,mpi_integer,0,mpi_comm_world,mpierr)
        call mpi_bcast(Par%outalmfile1,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
        call mpi_bcast(Par%endnamealm1,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
        call mpi_bcast(Par%outclfile,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
        call mpi_bcast(Par%endnamecl,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
        if (Par%do_cross) then
           call mpi_bcast(Par%inmapfile2,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
           call mpi_bcast(Par%endnamemap2,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
           call mpi_bcast(Par%outalmfile2,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
           call mpi_bcast(Par%endnamealm2,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
        endif  
     endif

     if (Par%compute_fskyl) then
        allocate(wpix(4*Par%nside+1,1:1))
        if (myid .eq. 0) call pixel_window(wpix,Par%nside)
     endif

     if (Par%compute_fskyl) then
        call mpi_bcast(Par%ampsignal,1,mpi_real8,0,mpi_comm_world,mpierr)
        call mpi_bcast(Par%amp_in_dl,1,mpi_logical,0,mpi_comm_world,mpierr)
        call mpi_bcast(Par%nsims_mask,1,mpi_integer,0,mpi_comm_world,mpierr)
        call mpi_bcast(Par%outfskyfile,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
        call mpi_bcast(wpix,4*Par%nside+1,mpi_real8,0,mpi_comm_world,mpierr)
     endif

     call mpi_bcast(Par%read_precomputed_alms,1,mpi_logical,0,mpi_comm_world,mpierr)

     !mask
     if (myid .eq. 0) then
        if (Par%inmaskfile1 .ne. '') then
           if (Par%feedback .gt. 1) write(*,*) 'Read mask 1'
           apply_mask1 = .true.
           npix = getsize_fits(trim(Par%inmaskfile1))
           allocate(mask1(0:npix-1,1:3)) 
           call read_mask_and_compute_fsky(Par%inmaskfile1,mask1,fsky1)
           if ((myid .eq. 0) .and. (Par%feedback .gt. 1)) write(*,*) 'fsky mask 1: ',fsky1
        endif
        if (Par%do_cross .and. (Par%inmaskfile2 .ne. '')) then
           if (Par%feedback .gt. 1) write(*,*) 'Read mask 2'
           apply_mask2 = .true.
           npix = getsize_fits(trim(Par%inmaskfile2))
           allocate(mask2(0:npix-1,1:3))     
           call read_mask_and_compute_fsky(Par%inmaskfile2,mask2,fsky2)
           if ((myid .eq. 0) .and. (Par%feedback .gt. 1)) write(*,*) 'fsky mask2: ',fsky2
        endif
     endif

     call mpi_barrier(mpi_comm_world, mpierr)  
     call mpi_bcast(apply_mask1,1,mpi_logical,0,mpi_comm_world,mpierr)
     if (apply_mask1) then
        if ((.not. Par%read_precomputed_alms) .or. (Par%compute_fskyl))then
           call mpi_bcast(npix,1,mpi_integer8,0,mpi_comm_world,mpierr)
           if (myid .ne. 0) allocate(mask1(0:npix-1,1:3))
           call mpi_bcast(mask1,npix*3,mpi_real8,0,mpi_comm_world,mpierr)
        else
           if (myid .eq. 0) deallocate(mask1)
        endif
        call mpi_bcast(fsky1,1,mpi_real8,0,mpi_comm_world,mpierr)
     endif
     if (Par%do_cross) then
        call mpi_bcast(apply_mask2,1,mpi_logical,0,mpi_comm_world,mpierr)
        if (apply_mask2) then
           if ((.not. Par%read_precomputed_alms) .or. (Par%compute_fskyl)) then
              call mpi_bcast(npix,1,mpi_integer8,0,mpi_comm_world,mpierr)
              if (myid .ne. 0) allocate(mask2(0:npix-1,1:3))
              call mpi_bcast(mask2,npix*3,mpi_real8,0,mpi_comm_world,mpierr)
           else
              if (myid .eq. 0) deallocate(mask2)
           endif
           call mpi_bcast(fsky2,1,mpi_real8,0,mpi_comm_world,mpierr)
        endif
     endif
     call mpi_barrier(mpi_comm_world, mpierr)

     if ((myid .eq. 0) .and. (Par%compute_fskyl) .and. (Par%feedback .gt. 1)) write(*,*) 'Compute fsky_l'
     if ((Par%compute_fskyl) .and. (apply_mask1 .or. apply_mask2)) then
        allocate(fsky_l(0:Par%Lmax),fsky_l_den(0:Par%Lmax))
        fsky_l = 0.d0
        fsky_l_den = 0.d0
        allocate(numerator(0:Par%Lmax,1:1),denominator(0:Par%Lmax,1:1))
        allocate(map(0:npix-1))
        allocate(alm(1:1,0:Par%Lmax,0:Par%Lmax))
        allocate(alm_amp_fid_fsky_l(0:Par%Lmax))
        
        if (Par%amp_in_dl) then
           do iL=1,Par%Lmax
              alm_amp_fid_fsky_l(iL) = sqrt(Par%ampsignal/iL/(iL+1)*TWOPI)
           enddo
           alm_amp_fid_fsky_l(0) = alm_amp_fid_fsky_l(1)
        else
           alm_amp_fid_fsky_l = sqrt(Par%ampsignal)
        endif

        allocate(maskt(0:npix-1))
        maskt = 1.d0
        if (apply_mask1) maskt = maskt*mask1(:,2)
        if (apply_mask2) maskt = maskt*mask2(:,2)
        call rand_init(rng_handle, 12345+myid, 6789012+myid)

        hsqrt2 = SQRT2 / 2.0_dp

        do isim=1,Par%nsims_mask
            if (myid .eq. mod(isim-1,nproc)) then
               if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' computing fsky_l of sim: ', isim
               do iL=0,Par%Lmax
                  alm(1,iL,0) = cmplx(rand_gauss(rng_handle),0.d0) * alm_amp_fid_fsky_l(iL) * wpix(iL,1)
                  do iM=1,iL 
                     alm(1,iL,iM) = cmplx(rand_gauss(rng_handle) * hsqrt2,rand_gauss(rng_handle) * hsqrt2) * alm_amp_fid_fsky_l(iL) * wpix(iL,1)
                  enddo
               enddo 
               call alm2map(Par%nside, Par%Lmax, Par%Lmax, alm, map)           

!               call map2alm(Par%nside, Par%Lmax, Par%Lmax, map, alm) 

               call alm2cl(Par%Lmax, Par%Lmax, alm, denominator)
               fsky_l_den = fsky_l_den + denominator(:,1)

               call map2alm(Par%nside, Par%Lmax, Par%Lmax, map*maskt, alm)
               call alm2cl(Par%Lmax, Par%Lmax, alm, numerator)
               fsky_l = fsky_l + numerator(:,1)
            endif
        enddo
        deallocate(alm_amp_fid_fsky_l,maskt,map,alm)
        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_allreduce(mpi_in_place,fsky_l_den,Par%Lmax+1,mpi_real8,mpi_sum,mpi_comm_world,mpierr)
        call mpi_allreduce(mpi_in_place,fsky_l,Par%Lmax+1,mpi_real8,mpi_sum,mpi_comm_world,mpierr)
        call mpi_barrier(mpi_comm_world, mpierr)
        if (myid .eq. 0) then
           fsky_l = fsky_l/fsky_l_den
           open(newunit=myunit,file=trim(Par%outfskyfile),status='replace',form='formatted')
           do iL=0,Par%Lmax
              write(myunit,'(I4,*(E15.7))') iL,fsky_l(iL),fsky_l_den(iL)
           enddo
           close(myunit)
        endif
        deallocate(fsky_l,fsky_l_den,numerator,denominator)
        if (Par%read_precomputed_alms) then
           if (apply_mask1) deallocate(mask1)
           if (Par%do_cross .and. apply_mask2) deallocate(mask2) 
        endif
     endif

     !allocation alms
     indmax = lm2index(Par%ellmax,Par%ellmax,Par%ellmax)
     if (Par%nsims .eq. 1) then
        if (myid .eq. 0) then
           allocate(almE1(1,0:indmax))
           allocate(almB1(1,0:indmax))
           almE1=0
           almB1=0
        endif
     else
        maxsimperproc = ceiling(Par%nsims/real(nproc))
        allocate(almE1(maxsimperproc,0:indmax))
        allocate(almB1(maxsimperproc,0:indmax))
        almE1=0
        almB1=0
     endif

     if (Par%do_cross) then
        if (Par%nsims .eq. 1) then
           if (myid .eq. 0) then
              allocate(almE2(1,0:indmax))
              allocate(almB2(1,0:indmax))
              almE2=0
              almB2=0
           endif
        else
           allocate(almE2(maxsimperproc,0:indmax))
           allocate(almB2(maxsimperproc,0:indmax))
           almE2=0
           almB2=0
        endif
     endif

     if ((myid .eq. 0) .and. (Par%feedback .gt. 1)) write(*,*) 'Read maps'
     if (Par%nsims .eq. 1) then
        if (myid .eq. 0) then
           if (.not. Par%read_precomputed_alms) then
              if (apply_mask1) then
                 call read_map_and_compute_alms(Par%inmapfile1,Par%niter,almE1(1,:),almB1(1,:),Par%ellmax,mask=mask1)
              else
                 call read_map_and_compute_alms(Par%inmapfile1,Par%niter,almE1(1,:),almB1(1,:),Par%ellmax)
              endif
           else
              call read_precomputed_alms(Par%inmapfile1,almE1(1,:),almB1(1,:))
           endif
        endif
        if (Par%do_cross) then
           if (myid .eq. 0) then
              if (.not. Par%read_precomputed_alms) then
                 if (apply_mask2) then
                    call read_map_and_compute_alms(Par%inmapfile2,Par%niter,almE2(1,:),almB2(1,:),Par%ellmax,mask=mask2)
                 else
                    call read_map_and_compute_alms(Par%inmapfile2,Par%niter,almE2(1,:),almB2(1,:),Par%ellmax)
                 endif
              else
                 call read_precomputed_alms(Par%inmapfile2,almE2(1,:),almB2(1,:))
              endif
           endif
        endif
     else
        ct=1
        write(strzerofill,fmt='(i1)') Par%zerofill
        do isim=Par%ssim,Par%ssim+Par%nsims-1
           if (myid .eq. mod(ct-1,nproc)) then
              if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' reading map ', isim
              write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
              mapname=trim(Par%inmapfile1)//trim(simstr)//trim(Par%endnamemap1)
              simproc = ceiling(ct/real(nproc))
              if (.not. Par%read_precomputed_alms) then
                 if (apply_mask1) then
                    call read_map_and_compute_alms(mapname,Par%niter,almE1(simproc,:),almB1(simproc,:),Par%ellmax,mask=mask1)
                 else
                    call read_map_and_compute_alms(mapname,Par%niter,almE1(simproc,:),almB1(simproc,:),Par%ellmax)
                 endif
              else
                 call read_precomputed_alms(mapname,almE1(simproc,:),almB1(simproc,:))
              endif
           endif
           ct=ct+1
        enddo
        
        if (Par%do_cross) then
           ct=1
           do isim=Par%ssim,Par%ssim+Par%nsims-1
              if (myid .eq. mod(ct-1,nproc)) then
                 if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' reading map ', isim
                 write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
                 mapname=trim(Par%inmapfile2)//trim(simstr)//trim(Par%endnamemap2)
                 simproc = ceiling(ct/real(nproc))
                 if (.not. Par%read_precomputed_alms) then
                    if (apply_mask2) then
                       call read_map_and_compute_alms(mapname,Par%niter,almE2(simproc,:),almB2(simproc,:),Par%ellmax,mask=mask2)
                    else
                       call read_map_and_compute_alms(mapname,Par%niter,almE2(simproc,:),almB2(simproc,:),Par%ellmax)
                    endif
                 else
                    call read_precomputed_alms(mapname,almE2(simproc,:),almB2(simproc,:))
                 endif
              endif
              ct=ct+1
           enddo
        endif
     endif
     if (apply_mask1 .and. (.not. Par%read_precomputed_alms)) deallocate(mask1)
     if (apply_mask2 .and. (.not. Par%read_precomputed_alms)) deallocate(mask2)
  endif

  if (Par%compute_biasalpha) then
     if (Par%do_cross) then
        allocate(clEEmap12(Par%nsims,0:Par%ellmax))
        allocate(clBBmap12(Par%nsims,0:Par%ellmax))
        clEEmap12 = 0
        clBBmap12 = 0
        if (Par%nsims .eq. 1) then 
           if (myid .eq. 0) call compute_cls_from_alms(almE1(1,:),almB1(1,:),Par%ellmax,clEEmap12(1,:),clBBmap12(1,:),almE2(1,:),almB2(1,:))
        else
           ct=1
           do isim=Par%ssim,Par%ssim+Par%nsims-1
              if (myid .eq. mod(ct-1,nproc)) then
                 simproc = ceiling(ct/real(nproc))
                 call compute_cls_from_alms(almE1(simproc,:),almB1(simproc,:),Par%ellmax,clEEmap12(ct,:),clBBmap12(ct,:),almE2(simproc,:),almB2(simproc,:))
              endif
              ct=ct+1
           enddo
        endif
        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_allreduce(mpi_in_place,clEEmap12,Par%nsims*(Par%ellmax+1),mpi_real8,mpi_sum,mpi_comm_world,mpierr)
        call mpi_allreduce(mpi_in_place,clBBmap12,Par%nsims*(Par%ellmax+1),mpi_real8,mpi_sum,mpi_comm_world,mpierr)
        if (apply_mask1) then
           clEEmap12 = clEEmap12 / sqrt(fsky1)
           clBBmap12 = clBBmap12 / sqrt(fsky1)
        endif
        if (apply_mask2) then
           clEEmap12 = clEEmap12 / sqrt(fsky2)
           clBBmap12 = clBBmap12 / sqrt(fsky2)
        endif
     else
        allocate(clEEmap1(Par%nsims,0:Par%ellmax))
        allocate(clBBmap1(Par%nsims,0:Par%ellmax))
        clEEmap1 = 0
        clEEmap1 = 0 
        if (Par%nsims .eq. 1) then
           if (myid .eq. 0) call compute_cls_from_alms(almE1(1,:),almB1(1,:),Par%ellmax,clEEmap1(1,:),clBBmap1(1,:))
        else
           ct=1
           do isim=Par%ssim,Par%ssim+Par%nsims-1
              if (myid .eq. mod(ct-1,nproc)) then
                 simproc = ceiling(ct/real(nproc))
                 call compute_cls_from_alms(almE1(simproc,:),almB1(simproc,:),Par%ellmax,clEEmap1(ct,:),clBBmap1(ct,:))
              endif
              ct=ct+1
           enddo
        endif
        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_allreduce(mpi_in_place,clEEmap1,Par%nsims*(Par%ellmax+1),mpi_real8,mpi_sum,mpi_comm_world,mpierr)
        call mpi_allreduce(mpi_in_place,clBBmap1,Par%nsims*(Par%ellmax+1),mpi_real8,mpi_sum,mpi_comm_world,mpierr)
        if (apply_mask1) then
           clEEmap1 = clEEmap1 / fsky1
           clBBmap1 = clBBmap1 / fsky1
        endif
     endif
  endif

  if (Par%compute_alphalm) indLMmax = lm2index(Par%Lmax,Par%Lmax,Par%Lmax)

  if (Par%compute_biasalpha .and. (.not. Par%compute_alphalm)) then
    if (Par%nsims .eq. 1) then
       if (myid .eq. 0) then
          deallocate(almE1,almB1)
          if (Par%do_cross) deallocate(almE2,almB2)
       endif
    else
       deallocate(almE1,almB1)
       if (Par%do_cross) deallocate(almE2,almB2)
    endif
  endif
  
  call mpi_barrier(mpi_comm_world, mpierr)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !! Start sigma and bias computations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'Starting Computation'
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock(t1, clock_rate, clock_max)
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write (0,*) 'Elapsed real time for initialization = ', real( t1 - t0 ) / real ( clock_rate )
  !! Compute sigma
  
  allocate(one_o_var1(0:Par%Lmax))
  one_o_var1 = 0.0

  if (Par%do_cross) then
    allocate(one_o_var2(0:Par%Lmax))
    one_o_var2 = 0.0
  endif  

  if (Par%compute_biasalpha) then
     allocate(biasalpha(1:Par%nsims,0:Par%Lmax))
     biasalpha = 0.0
  endif

  if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'Computing sigma and bias'

  do iL=0,Par%Lmax
     if (myid .eq. mod(iL,nproc)) then
        if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' doing multipole ', iL
        do iell=Par%ellmin,Par%ellmax
           allocate(wig2(iL+iell+1))
           call Wigner3j(wig2, ellpmin, ellpmax, iell, iL, -2, 2, 0)
           usedellpmax = min(ellpmax,Par%ellmax)
           allocate(F_EB1(ellpmin:usedellpmax),F_BE1(ellpmin:usedellpmax))
           F_EB1 = 2 * wig2(1:usedellpmax-ellpmin+1) * clEEfid(iell)*bl1(iell,myE)*bl1(ellpmin:usedellpmax,myE)
           F_BE1 = 2 * wig2(1:usedellpmax-ellpmin+1) * clEEfid(ellpmin:usedellpmax)*bl1(ellpmin:usedellpmax,myE)*bl1(iell,myE)
           if (Par%do_cross) then
              allocate(F_EB2(ellpmin:usedellpmax),F_BE2(ellpmin:usedellpmax))
              F_EB2 = 2 * wig2(1:usedellpmax-ellpmin+1) * clEEfid(iell)*bl2(iell,myE)*bl2(ellpmin:usedellpmax,myE)
              F_BE2 = 2 * wig2(1:usedellpmax-ellpmin+1) * clEEfid(ellpmin:usedellpmax)*bl2(ellpmin:usedellpmax,myE)*bl2(iell,myE)
           endif
           Gl = (2*iell + 1)/FOURPI
           do iellp = max(ellpmin,Par%ellmin,iell),usedellpmax              
              if ((jmod(iL+iell+iellp,2) .eq. 0) .and. (wig2(iellp-ellpmin+1) .ne. 0.0d0)) then
!              if ((iellp .ge. iell) .and. (jmod(iL+iell+iellp,2) .eq. 0) .and. (iL .ge. abs(iellp-iell)) .and. (iL .le. iellp+iell)) then
                 if (iellp .eq. iell) then
                    factor = 0.5 * Gl * (2.0*iellp + 1.0)
                 else 
                    factor = Gl * (2.0*iellp + 1.0)
                 endif
                 one_o_var1(iL) = one_o_var1(iL) + factor * &
                       (F_EB1(iellp)**2/clEEobs1(iell)/clBBobs1(iellp) + &
                        F_BE1(iellp)**2/clBBobs1(iell)/clEEobs1(iellp))                  
                 if (Par%do_cross) then
                    one_o_var2(iL) = one_o_var2(iL) + factor * &
                         (F_EB2(iellp)**2/clEEobs2(iell)/clBBobs2(iellp) + &
                          F_BE2(iellp)**2/clBBobs2(iell)/clEEobs2(iellp))
                 endif
                 !bias computation
                 if (Par%compute_biasalpha) then
                    if (iellp .eq. iell) then
                       factor = 0.25 * Gl * (2.0*iellp + 1.0)
                    else
                       factor = Gl * (2.0*iellp + 1.0)
                    endif
                    if (Par%do_cross) then 
                        biasalpha(:,iL) = biasalpha(:,iL) + factor * &
                          (F_EB1(iellp)*F_EB2(iellp) * clEEmap12(:,iell) * clBBmap12(:,iellp) / &
                             (clBBobs1(iellp) * clBBobs2(iellp) * clEEobs1(iell) * clEEobs2(iell)) + &
                              F_BE1(iellp)*F_BE2(iellp) * clEEmap12(:,iellp) * clBBmap12(:,iell) / &
                             (clBBobs1(iell) * clBBobs1(iell) * clEEobs2(iellp) * clEEobs2(iellp)))
                    else
                        biasalpha(:,iL) = biasalpha(:,iL) + factor * &
                          (F_EB1(iellp)**2 * clEEmap1(:,iell) * clBBmap1(:,iellp) / &
                             (clBBobs1(iellp) * clBBobs1(iellp) * clEEobs1(iell) * clEEobs1(iell)) + &
                              F_BE1(iellp)**2 * clEEmap1(:,iellp) * clBBmap1(:,iell) / &
                             (clBBobs1(iell) * clBBobs1(iell) * clEEobs1(iellp) * clEEobs1(iellp)))
                    endif 
                 endif
              endif
           enddo
           deallocate(wig2,F_EB1,F_BE1)
           if (Par%do_cross) deallocate(F_EB2,F_BE2)
        enddo
     endif
  enddo

  if (Par%compute_biasalpha) then
     if (Par%do_cross) then
        deallocate(clEEmap12,clBBmap12) 
     else
        deallocate(clEEmap1,clBBmap1)
     endif
  endif 

  call mpi_barrier(mpi_comm_world, mpierr)
  if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'Computation of sigma done'
  !compute timing
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock( t2, clock_rate, clock_max )
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write (0,*) 'Elapsed real time for sigma computation = ', real ( t2 - t1 ) / real ( clock_rate )    

  call mpi_barrier(mpi_comm_world, mpierr)
  call mpi_allreduce(mpi_in_place,one_o_var1,Par%Lmax+1,mpi_real8,mpi_sum,mpi_comm_world,mpierr)

  if (Par%do_cross) then 
     call mpi_barrier(mpi_comm_world, mpierr)
     call mpi_allreduce(mpi_in_place,one_o_var2,Par%Lmax+1,mpi_real8,mpi_sum,mpi_comm_world,mpierr)
  endif
  
  if (myid .eq. 0) then
     if (Par%feedback .gt. 1) write(0,*) 'Write out sigma'
     open(newunit=myunit,file=trim(Par%outsigmafile1),status='replace',form='formatted')
     do iL=0,Par%Lmax
        write(myunit,'(I4,*(E15.7))') iL,sqrt(1.0/one_o_var1(iL))
     enddo
     close(myunit)
     if (Par%do_cross) then
        open(newunit=myunit,file=trim(Par%outsigmafile2),status='replace',form='formatted')
        do iL=0,Par%Lmax
           write(myunit,'(I4,*(E15.7))') iL,sqrt(1.0/one_o_var2(iL))
        enddo
        close(myunit)
     endif
  endif

  if (Par%compute_biasalpha) then
     call mpi_barrier(mpi_comm_world, mpierr)
     call mpi_allreduce(mpi_in_place,biasalpha,(Par%Lmax+1)*Par%nsims,mpi_real8,mpi_sum,mpi_comm_world,mpierr)
  endif

  call mpi_barrier(mpi_comm_world, mpierr)

  if (Par%do_cross) then
     do iL=0,Par%Lmax
        biasalpha(:,iL) = biasalpha(:,iL) / one_o_var1(iL) / one_o_var2(iL)
     enddo
  else
     do iL=0,Par%Lmax
        biasalpha(:,iL) = biasalpha(:,iL) / one_o_var1(iL)**2
     enddo
  endif

  call mpi_barrier(mpi_comm_world, mpierr)

  if (Par%compute_biasalpha .and. (myid .eq. 0)) then
     if (Par%feedback .gt. 1) write(0,*) 'Write out bias'
     call write_out_cls(Par%outbiasfile,Par%ssim,Par%zerofill,Par%endnamebias,biasalpha)
  endif

  call mpi_barrier(mpi_comm_world, mpierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !! Start alphalm computation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  if (Par%compute_alphalm) then

     if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'Computing alpha_LM'
      
     if (Par%nsims .eq. 1) then
        if (myid .eq. 0) then
           do ind=0,indmax
              lm=index2lm(Par%ellmax,ind)
              almE1(1,ind) = conjg(almE1(1,ind)) * clEEfid(lm(1)) * bl1(lm(1),myE) / clEEobs1(lm(1))
              almB1(1,ind) = conjg(almB1(1,ind)) * bl1(lm(1),myB) / clBBobs1(lm(1))
           enddo
           if (Par%do_cross) then
              do ind=0,indmax
                 lm=index2lm(Par%ellmax,ind)
                 almE2(1,ind) = conjg(almE2(1,ind)) * clEEfid(lm(1)) * bl2(lm(1),myE) / clEEobs2(lm(1))
                 almB2(1,ind) = conjg(almB2(1,ind)) * bl2(lm(1),myB) / clBBobs2(lm(1))
              enddo
           endif
        endif
     else
        do ind=0,indmax
           lm=index2lm(Par%ellmax,ind)
           almE1(:,ind) = conjg(almE1(:,ind)) * clEEfid(lm(1)) * bl1(lm(1),myE) / clEEobs1(lm(1))
           almB1(:,ind) = conjg(almB1(:,ind)) * bl1(lm(1),myB) / clBBobs1(lm(1))
        enddo
        if (Par%do_cross) then
           do ind=0,indmax
              lm=index2lm(Par%ellmax,ind)
              almE2(:,ind) = conjg(almE2(:,ind)) * clEEfid(lm(1)) * bl2(lm(1),myE) / clEEobs2(lm(1))
              almB2(:,ind) = conjg(almB2(:,ind)) * bl2(lm(1),myB) / clBBobs2(lm(1))
           enddo
        endif
     endif

     deallocate(clEEfid,bl1,clEEobs1,clBBobs1)
     if (Par%do_cross) deallocate(bl2,clEEobs2,clBBobs2)

     allocate(almalpha1(0:indLMmax))
     almalpha1=0
     if (Par%do_cross) then
        allocate(almalpha2(0:indLMmax))
        almalpha2=0
     endif

     call mpi_barrier(mpi_comm_world, mpierr)
     allocate(alm1(1:1,0:Par%Lmax,0:Par%Lmax))
     if (Par%do_cross) allocate(alm2(1:1,0:Par%Lmax,0:Par%Lmax))
     
     if (Par%nsims .eq. 1) then
        if (myid .eq. 0) then
           call compute_alphalm(almE1(1,:),almB1(1,:),Par%Lmax,Par%ellmin,Par%ellmax,Par%nside,almalpha1)
           call reorder_and_normalize_alms(almalpha1,one_o_var1,alm1(1,:,:))
           call write_out_alm(Par%outalmfile1,alm1(1,:,:)) 
           if (Par%do_cross) then
              call compute_alphalm(almE2(1,:),almB2(1,:),Par%Lmax,Par%ellmin,Par%ellmax,Par%nside,almalpha2)
              call reorder_and_normalize_alms(almalpha2,one_o_var2,alm2(1,:,:))
              call write_out_alm(Par%outalmfile2,alm2(1,:,:))
           endif           
           if (Par%do_cross) then 
              if (Par%compute_biasalpha .and. Par%subtract_bias) then
                 call compute_and_write_cl(Par%outclfile,alm1,bias=biasalpha(1,:),alm2=alm2)
              else
                 call compute_and_write_cl(Par%outclfile,alm1,alm2=alm2)                 
              endif
           else
              if (Par%compute_biasalpha .and. Par%subtract_bias) then
                 call compute_and_write_cl(Par%outclfile,alm1,bias=biasalpha(1,:))
              else
                 call compute_and_write_cl(Par%outclfile,alm1)
              endif
           endif
        endif
     else
        write(strzerofill,fmt='(i1)') Par%zerofill
        ct=1
        do isim=Par%ssim,Par%ssim+Par%nsims-1
           if (myid .eq. mod(ct-1,nproc)) then
              if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' processing sim ', isim
              simproc = ceiling(ct/real(nproc))
              call compute_alphalm(almE1(simproc,:),almB1(simproc,:),Par%Lmax,Par%ellmin,Par%ellmax,Par%nside,almalpha1)
              call reorder_and_normalize_alms(almalpha1,one_o_var1,alm1(1,:,:))
              write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
              filename=trim(Par%outalmfile1)//trim(simstr)//trim(Par%endnamealm1)
              call write_out_alm(filename,alm1(1,:,:))
              if (Par%do_cross) then
                 call compute_alphalm(almE2(simproc,:),almB2(simproc,:),Par%Lmax,Par%ellmin,Par%ellmax,Par%nside,almalpha2)
                 call reorder_and_normalize_alms(almalpha2,one_o_var2,alm2(1,:,:))
                 write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
                 filename=trim(Par%outalmfile2)//trim(simstr)//trim(Par%endnamealm2)
                 call write_out_alm(filename,alm2(1,:,:))
              endif

              write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
              filename=trim(Par%outclfile)//trim(simstr)//trim(Par%endnamecl)

              if (Par%do_cross) then
                 if (Par%compute_biasalpha .and. Par%subtract_bias) then
                    call compute_and_write_cl(filename,alm1,bias=biasalpha(ct,:),alm2=alm2)
                 else
                    call compute_and_write_cl(filename,alm1,alm2=alm2)
                 endif
              else
                 if (Par%compute_biasalpha .and. Par%subtract_bias) then
                    call compute_and_write_cl(filename,alm1,bias=biasalpha(ct,:))
                 else
                    call compute_and_write_cl(filename,alm1)
                 endif
              endif
           endif
           ct=ct+1
        enddo
     endif
  
     call mpi_barrier(mpi_comm_world, mpierr)
     
     if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'Computation of alphalm done'
     if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock( t3, clock_rate, clock_max )
     if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write(0,*) 'Elapsed real time for alphalm computation = ', real ( t3 - t2 ) / real ( clock_rate )

     call mpi_barrier(mpi_comm_world, mpierr)
     
     deallocate(almalpha1,alm1)
     if (Par%do_cross) then
        deallocate(almalpha2,alm2)
     endif
     
  else
     deallocate(clEEfid,clEEobs1,clBBobs1,bl1)
     if (Par%do_cross) deallocate(clEEobs2,clBBobs2,bl2)     
  endif

  deallocate(one_o_var1)
  if (Par%do_cross) deallocate(one_o_var2)

  if (Par%compute_biasalpha) deallocate(biasalpha)

  call mpi_barrier(mpi_comm_world, mpierr)
  
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock(t4, clock_rate, clock_max)
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write(0,*) 'Elapsed real time for total computation = ', real(t4 - t1) / real(clock_rate)
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write(0,*) 'Elapsed real time for computation and initialization = ', real(t4 - t0)/real(clock_rate)
  if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'End Program'  
  call mpi_finalize(mpierr)
  
end program EB_estimator
