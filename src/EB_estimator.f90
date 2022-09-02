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
  
  implicit none
  
  !variables
  integer :: myid,nproc,mpierr
  integer :: j,jmin,jmax, iell
  integer :: jminall,jmaxall
  integer :: isim,Lcount,iL,nele,elementcount
  integer :: iemm, iM, iemmp
  logical :: emmppos
  real(dp), allocatable, dimension(:) :: one_o_var,red_one_o_var,wig2,wigall,csi
  real(dp), allocatable, dimension(:) :: F_EB,F_BE,clEEfid,clEEobs,clBBobs
  real(dp), allocatable, dimension(:,:) :: clEEmap,clBBmap,biasalpha,red_biasalpha
  complex(dpc), allocatable, dimension(:,:,:) :: almE,almB,almalpha,red_almalpha
  complex(dpc), allocatable, dimension(:) :: curralmE,curralmB,curralmEstar,curralmBstar
  complex(dpc) :: EB_csi, BE_csi
  real(dp), allocatable,dimension(:,:) :: clfid,wl,bl,nl
  real(dp) :: factor,exp_m,Gl,norm
  integer :: t1=0,t2,t3,t4,clock_rate,clock_max,myunit
  character(len=16) :: simstr
  character(len=1) :: strzerofill
  
  !input parameters
  Type(Params) :: Par
  
  !  !! ================================ Init ================================ !!
  call mpi_init(mpierr)
  call mpi_comm_size(mpi_comm_world, nproc, mpierr)
  call mpi_comm_rank(mpi_comm_world, myid, mpierr)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !! Master reads parameters and does some stuff: reads maks and covmats, computes inverted masks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  if (myid .eq. 0) then 
     if (Par%feedback .gt. 0) write(*,*) 'Reading parameters'
     
     call read_parameter_file(Par)
     
     allocate(bl(0:Par%ellmax+Par%elloffset,3),wl(0:Par%ellmax+Par%elloffset,6))
     call read_beam(bl,wl,Par%inbeamfile)
     
     allocate(clfid(0:Par%ellmax+Par%elloffset,6))
     call read_cl(Par%inclfile,clfid)

     allocate(nl(0:Par%ellmax+Par%elloffset,2))
     call make_noise(Par%innoisefile,nl,Par%noiseE,Par%noiseB) 
     
     allocate(clEEfid(0:Par%ellmax+Par%elloffset))
     allocate(clEEobs(0:Par%ellmax+Par%elloffset))
     allocate(clBBobs(0:Par%ellmax+Par%elloffset))
     clEEfid = clfid(:,myEE)
     clEEobs = clfid(:,myEE) * wl(:,myEE) + nl(:,1)
     clBBobs = clfid(:,myBB) * wl(:,myBB) + nl(:,2)
     deallocate(clfid,nl)

     if (Par%compute_alphalm) then
        if (Par%feedback .gt. 1) write(*,*) 'Read maps'
        allocate(almE(Par%nsims,0:Par%ellmax+Par%elloffset,0:Par%ellmax+Par%elloffset)) 
        allocate(almB(Par%nsims,0:Par%ellmax+Par%elloffset,0:Par%ellmax+Par%elloffset))
        
        call read_maps_and_compute_alms(Par%inmapfile,Par%ssim,Par%zerofill,Par%endnamemap,almE,almB)
        if (Par%compute_biasalpha) then
           allocate(clEEmap(Par%nsims,0:Par%ellmax+Par%elloffset))
           allocate(clBBmap(Par%nsims,0:Par%ellmax+Par%elloffset))
           call compute_cls_from_alms(almE,almB,clEEmap,clBBmap)
        endif
     endif
  endif

  call mpi_barrier(mpi_comm_world, mpierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !! Broadcast relevant parameters and variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  if ((myid .eq. 0) .and. (Par%feedback .gt. 1)) write(*,*) 'Broadcast parameters and variables'
  
  call mpi_bcast(Par%feedback,1,mpi_integer,0, mpi_comm_world, mpierr)
  call mpi_bcast(Par%Lmin,1,mpi_integer,0,mpi_comm_world, mpierr)
  call mpi_bcast(Par%Lmax,1,mpi_integer,0,mpi_comm_world, mpierr)
  call mpi_bcast(Par%ellmin,1,mpi_integer,0,mpi_comm_world, mpierr)
  call mpi_bcast(Par%ellmax,1,mpi_integer,0,mpi_comm_world, mpierr)
  call mpi_bcast(Par%elloffset,1,mpi_integer,0,mpi_comm_world, mpierr)
  call mpi_bcast(Par%compute_alphalm,1,mpi_logical,0,mpi_comm_world,mpierr)

  if (Par%compute_alphalm) then
     call mpi_bcast(Par%nsims,1,mpi_integer,0,mpi_comm_world,mpierr)
     call mpi_bcast(Par%ssim,1,mpi_integer,0,mpi_comm_world,mpierr)
     call mpi_bcast(Par%compute_biasalpha,1,mpi_logical,0,mpi_comm_world,mpierr)
  endif

  nele = Par%ellmax+Par%elloffset+1
  if (myid .ne. 0) then
     if (Par%compute_alphalm) then
        allocate(almE(Par%nsims,0:Par%ellmax+Par%elloffset,0:Par%ellmax+Par%elloffset))
        allocate(almB(Par%nsims,0:Par%ellmax+Par%elloffset,0:Par%ellmax+Par%elloffset))
        if (Par%compute_biasalpha) then
           allocate(clEEmap(Par%nsims,0:Par%ellmax+Par%elloffset))
           allocate(clBBmap(Par%nsims,0:Par%ellmax+Par%elloffset))
        endif
     endif
     allocate(clEEfid(0:Par%ellmax+Par%elloffset),clEEobs(0:Par%ellmax+Par%elloffset),clBBobs(0:Par%ellmax+Par%elloffset))
     allocate(bl(0:Par%ellmax+Par%elloffset,3),wl(0:Par%ellmax+Par%elloffset,6))
  endif
  
  if (Par%compute_alphalm) then
     call mpi_bcast(almE,Par%nsims*nele*nele,mpi_double_complex,0,mpi_comm_world, mpierr)
     call mpi_bcast(almB,Par%nsims*nele*nele,mpi_double_complex,0,mpi_comm_world, mpierr)
     if (Par%compute_biasalpha) then
        call mpi_bcast(clEEmap,Par%nsims*nele,mpi_real8,0,mpi_comm_world, mpierr)
        call mpi_bcast(clBBmap,Par%nsims*nele,mpi_real8,0,mpi_comm_world, mpierr)
     endif
  endif
  
  call mpi_bcast(clEEfid,nele,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(clEEobs,nele,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(clBBobs,nele,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(bl,3*nele,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(wl,6*nele,mpi_real8,0,mpi_comm_world, mpierr)
  
  call mpi_barrier(mpi_comm_world, mpierr)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !! Start E operator computation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(*,*) 'Starting Computation'
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock ( t1, clock_rate, clock_max)
  
  !! Compute sigma
  
  allocate(one_o_var(Par%Lmin:Par%Lmax))
  one_o_var = 0.0
  
  Lcount = 0
  do iL=Par%Lmin,Par%Lmax
     if (myid .eq. mod(Lcount,nproc)) then
        if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' doing multipole ', iL
        do iell=Par%ellmin,Par%ellmax
           allocate(wig2(iL+iell+1))
           call Wigner3j(wig2, jmin, jmax, iell, iL, -2, 2, 0)
           allocate(F_EB(jmin:jmax),F_BE(jmin:jmax))
           F_EB = (2 * wig2(1:jmax-jmin+1) * clEEfid(iell)*bl(iell,1)*bl(jmin:jmax,1))**2
           F_BE = (2 * wig2(1:jmax-jmin+1) * clEEfid(jmin:jmax)*bl(jmin:jmax,1)*bl(iell,1))**2
           Gl = (2*iell + 1)/FOURPI
           do j = jmin,jmax
              if (j .eq. iell) then
                 factor = 0.5 
              else 
                 factor = 1
              endif
              one_o_var(iL) = one_o_var(iL) + factor * &
                   Gl * (2.0*j + 1.0) * &
                   (F_EB(j)/clEEobs(iell)/clBBobs(j) + &
                    F_BE(j)/clBBobs(iell)/clEEobs(j))
           enddo
           deallocate(wig2,F_EB,F_BE)
        enddo
     endif
     Lcount = Lcount + 1
  enddo

  call mpi_barrier(mpi_comm_world, mpierr)
  if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'Computation of sigma done'
  !compute timing
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock ( t2, clock_rate, clock_max )
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write (*,*) 'Elapsed real time for sigma computation = ', real ( t2 - t1 ) / real ( clock_rate )    
  
  if (myid .eq. 0) then
     allocate(red_one_o_var(Par%Lmin:Par%Lmax))
     red_one_o_var=0.0
  endif
  
  call mpi_barrier(mpi_comm_world, mpierr) 
  call mpi_reduce(one_o_var,red_one_o_var,Par%Lmax-Par%Lmin+1,mpi_real8,mpi_sum,0,mpi_comm_world, mpierr)
  deallocate(one_o_var)
  call mpi_barrier(mpi_comm_world, mpierr)
  
  if (myid .eq. 0) then
     if (Par%feedback .gt. 1) write(*,*) 'Write out sigma'
     open(newunit=myunit,file=trim(Par%outsigmafile),status='replace',form='formatted')
     do iL=Par%Lmin,Par%Lmax
        write(myunit,'(I4,*(E15.7))') iL,sqrt(1.0/red_one_o_var(iL))
     enddo
     close(myunit)
  endif
  
  !! Compute alphalm
  if (Par%compute_alphalm) then
     
     call mpi_barrier(mpi_comm_world, mpierr)
     allocate(almalpha(1:Par%nsims,0:Par%Lmax,0:Par%Lmax))
     allocate(curralmE(Par%nsims),curralmB(Par%nsims),curralmEstar(Par%nsims),curralmBstar(Par%nsims))
     almalpha=0
        
     if (Par%compute_biasalpha) then 
        allocate(biasalpha(1:Par%nsims,Par%Lmin:Par%Lmax))
        biasalpha = 0
     endif

     elementcount = 0
     do iL=Par%Lmin,Par%Lmax
        do iM=0,iL
           if (myid .eq. mod(elementcount,nproc)) then
              if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' doing L=', iL, ' and M=',iM
              !loop ell
              do iell=Par%ellmin,Par%ellmax
                 allocate(wig2(iL+iell+1),wigall(iL+iell+1))
                 call Wigner3j(wig2, jmin, jmax, iell, iL, -2, 2, 0)
                 allocate(F_EB(jmin:jmax),F_BE(jmin:jmax))
                 F_EB = 2 * wig2(1:jmax-jmin+1) * clEEfid(iell)*bl(iell,1)*bl(jmin:jmax,1)
                 F_BE = 2 * wig2(1:jmax-jmin+1) * clEEfid(jmin:jmax)*bl(jmin:jmax,1)*bl(iell,1)
                 !compute bias if requested
                 if (Par%compute_biasalpha) then
                    Gl = (2.0*iell + 1)/FOURPI
                    do j = jmin,jmax
                       if (j .eq. iell) then
                          factor = 0.25
                       else
                          factor = 1
                       endif
                       do isim=1,Par%nsims
                          biasalpha(isim,iL) = biasalpha(isim,iL) + factor * &
                             Gl * (2.0*j + 1.0) * & 
                             (F_EB(j)**2 * clEEmap(isim,iell) * clBBmap(isim,j) / &
                             (clBBobs(j) * clBBobs(j) * clEEobs(iell) * clEEobs(iell)) + &
                              F_BE(j)**2 * clEEmap(isim,j) * clBBmap(isim,iell) / &
                             (clBBobs(iell) * clBBobs(iell) * clEEobs(j) * clEEobs(j)))
                       enddo
                    enddo
                 endif

                 norm = sqrt((2.0*iell + 1)*(2.0*iL + 1)/FOURPI)
                 !loop emm
                 do iemm=-iell,iell
                    iemmp = iemm-iM
                    call Wigner3j(wigall, jminall, jmaxall, iell, iL, iemmp , -iemm, iM)

                    if (iemm .ge. 0) then
                       curralmE = almE(:,iell,iemm)
                       curralmB = almB(:,iell,iemm)
                    else 
                       curralmE = conjg(almE(:,iell,-iemm))
                       curralmB = conjg(almB(:,iell,-iemm))
                    endif

                    emmppos = .true.
                    if (iemmp .lt. 0) emmppos=.false.

                    exp_m = -1.0**iemm
                    allocate(csi(jminall:jmaxall))
                    do j=jminall,jmaxall
                       csi(j) = sqrt(2.0*j + 1) * wigall(j-jminall+1)
                    enddo
                    csi = csi * exp_m * norm

                    !loop ell'
                    do j = max(jminall,jmin),jmaxall
                       if (j .eq. iell) then
                          factor = 0.5
                       else
                          factor = 1
                       endif

                       if (emmppos) then
                          curralmEstar = conjg(almE(:,j,iemmp))
                          curralmBstar = conjg(almB(:,j,iemmp))
                       else
                          curralmEstar = almE(:,j,-iemmp)
                          curralmBstar = almB(:,j,-iemmp)
                       endif
                        
                       do isim=1,Par%nsims
                          EB_csi = csi(j) * curralmE(isim) * curralmBstar(isim)
                          BE_csi = csi(j) * curralmB(isim) * curralmEstar(isim)
                          almalpha(isim,iL,iM) = almalpha(isim,iL,iM) + factor * &
                               (F_EB(j) * EB_csi / clBBobs(j) / clEEobs(iell) + &
                                F_BE(j) * BE_csi / clBBobs(iell) / clEEobs(j))
                       enddo
                    enddo
                    deallocate(csi)
                 enddo
                 deallocate(wig2,wigall,F_EB,F_BE)
              enddo
           endif
           elementcount = elementcount + 1
        enddo
     enddo
     deallocate(curralmE,curralmB,curralmEstar,curralmBstar)

     call mpi_barrier(mpi_comm_world, mpierr)
     if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'Computation of alphalm done'
     if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock ( t3, clock_rate, clock_max )
     if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write (*,*) 'Elapsed real time for alphalm computation = ', real ( t3 - t2 ) / real ( clock_rate )

     call mpi_barrier(mpi_comm_world, mpierr)
     
     if (myid .eq. 0) then
        allocate(red_almalpha(1:Par%nsims,0:Par%Lmax,0:Par%Lmax))
        red_almalpha=0.0
     endif

     call mpi_barrier(mpi_comm_world, mpierr)
     nele = Par%Lmax+1
     call mpi_reduce(almalpha,red_almalpha,nele*nele*Par%nsims,mpi_double_complex,mpi_sum,0,mpi_comm_world,mpierr)
     deallocate(almalpha)
     call mpi_barrier(mpi_comm_world, mpierr)

     if (Par%compute_biasalpha) then
        if (myid .eq. 0) then
           allocate(red_biasalpha(1:Par%nsims,Par%Lmin:Par%Lmax))
           red_biasalpha=0.0
        endif
        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_reduce(biasalpha,red_biasalpha,nele*Par%nsims,mpi_real8,mpi_sum,0,mpi_comm_world,mpierr)
        deallocate(biasalpha)
     endif

     call mpi_barrier(mpi_comm_world, mpierr)
     
     if (myid .eq. 0) then
        do iL=Par%Lmin,Par%Lmax
           red_almalpha(:,iL,:) = red_almalpha(:,iL,:) / red_one_o_var(iL)
        enddo
        if (Par%compute_biasalpha) then
           do iL=Par%Lmin,Par%Lmax
              red_biasalpha(:,iL) = red_biasalpha(:,iL) / red_one_o_var(iL)**2
           enddo
           if (Par%feedback .gt. 1) write(*,*) 'Write out bias'
           call write_out_cls(Par%outbiasfile,Par%ssim,Par%zerofill,Par%endnamebias,red_biasalpha,Par%Lmin)
        endif
        if (Par%feedback .gt. 1) write(*,*) 'Write out alphalm'
        call write_out_alms(Par%outalmfile,Par%ssim,Par%zerofill,Par%endnamealm,red_almalpha)
        if (Par%compute_alphacl) then
           if (Par%compute_biasalpha) then
              call compute_and_write_cl(Par%outclfile,Par%ssim,Par%zerofill,Par%endnamecl,red_almalpha,Par%Lmin,red_biasalpha)
           else
              call compute_and_write_cl(Par%outclfile,Par%ssim,Par%zerofill,Par%endnamecl,red_almalpha,Par%Lmin)
           endif
        endif
        deallocate(red_almalpha,red_one_o_var)
        if (Par%compute_biasalpha) deallocate(red_biasalpha)
     endif
  else
     if (myid .eq. 0) deallocate(red_one_o_var)
  endif
  
  call mpi_barrier(mpi_comm_world, mpierr)
  
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock ( t4, clock_rate, clock_max )
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write (*,*) 'Elapsed real time for total computation = ', real ( t4 - t1 ) / real ( clock_rate )
  if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(*,*) 'End Program'  
  call mpi_finalize(mpierr)
  
end program EB_estimator
