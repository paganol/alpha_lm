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
  integer(i4b) :: nmaps
  integer(i8b) :: npix
  logical :: emmppos,apply_mask1=.false.,apply_mask2=.false.
  real(dp), allocatable, dimension(:) :: one_o_var1,red_one_o_var1,one_o_var2,red_one_o_var2,wig2,wigall,csi
  real(dp), allocatable, dimension(:) :: F_EB1,F_BE1,clEEfid,clEEobs1,clBBobs1
  real(dp), allocatable, dimension(:) :: F_EB2,F_BE2,clEEobs2,clBBobs2
  real(dp), allocatable, dimension(:,:) :: clEEmap1,clBBmap1,clEEmap12,clBBmap12,biasalpha,red_biasalpha
  complex(spc), allocatable, dimension(:,:,:) :: almE1,almB1,almE2,almB2
  complex(spc), allocatable, dimension(:,:,:) :: almalpha1,red_almalpha1,almalpha2,red_almalpha2
  complex(dpc), allocatable, dimension(:) :: curralmE1,curralmB1,curralmE1star,curralmB1star
  complex(dpc), allocatable, dimension(:) :: curralmE2,curralmB2,curralmE2star,curralmB2star
  real(dp), allocatable,dimension(:,:) :: clfid,wl1,bl1,nl1,wl2,bl2,nl2,mask1,mask2
  real(dp) :: factor,Gl,norm,fsky1,fsky2  
  integer :: t0,t1,t2,t3,t4,clock_rate,clock_max,myunit,ct,zerofill
  character(len=FILENAMELEN) :: mapname
  character(len=16) :: simstr
  character(len=1) :: strzerofill
  
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

     allocate(bl1(0:Par%ellmax+Par%Lmax,3),wl1(0:Par%ellmax+Par%Lmax,6))
     call read_beam(Par%inbeamfile1,bl1,wl1)
     if (Par%do_cross) then
        allocate(bl2(0:Par%ellmax+Par%Lmax,3),wl2(0:Par%ellmax+Par%Lmax,6))
        call read_beam(Par%inbeamfile2,bl2,wl2)
     endif

     allocate(clfid(0:Par%ellmax+Par%Lmax,6))
     call read_cl(Par%inclfile,clfid)
     allocate(clEEfid(0:Par%ellmax+Par%Lmax))
     clEEfid = clfid(:,myEE)

     allocate(nl1(0:Par%ellmax,2))
     call make_noise(Par%innoisefile1,nl1,Par%noiseE1,Par%noiseB1) 
     if (Par%do_cross) then     
        allocate(nl2(0:Par%ellmax,2))
        call make_noise(Par%innoisefile2,nl2,Par%noiseE2,Par%noiseB2)
     endif   

     allocate(clEEobs1(0:Par%ellmax))
     allocate(clBBobs1(0:Par%ellmax))
     clEEobs1 = clfid(0:Par%ellmax,myEE) * wl1(0:Par%ellmax,myEE) + nl1(:,1)
     clBBobs1 = clfid(0:Par%ellmax,myBB) * wl1(0:Par%ellmax,myBB) + nl1(:,2)
     deallocate(nl1)
 
     if (Par%do_cross) then
        allocate(clEEobs2(0:Par%ellmax))
        allocate(clBBobs2(0:Par%ellmax))
        clEEobs2 = clfid(0:Par%ellmax,myEE) * wl2(0:Par%ellmax,myEE) + nl2(:,1)
        clBBobs2 = clfid(0:Par%ellmax,myBB) * wl2(0:Par%ellmax,myBB) + nl2(:,2)
        deallocate(nl2)
     endif
   
     deallocate(clfid)

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
  call mpi_bcast(Par%compute_alphalm,1,mpi_logical,0,mpi_comm_world,mpierr)
  call mpi_bcast(Par%compute_biasalpha,1,mpi_logical,0,mpi_comm_world,mpierr)
  call mpi_bcast(Par%do_cross,1,mpi_logical,0,mpi_comm_world,mpierr)

  nele = Par%ellmax+Par%Lmax+1

  if (myid .ne. 0) then
     allocate(clEEfid(0:Par%ellmax+Par%Lmax),clEEobs1(0:Par%ellmax),clBBobs1(0:Par%ellmax))
     allocate(bl1(0:Par%ellmax+Par%Lmax,3),wl1(0:Par%ellmax+Par%Lmax,6))
  endif
  call mpi_bcast(clEEfid,nele,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(clEEobs1,Par%ellmax+1,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(clBBobs1,Par%ellmax+1,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(bl1,3*nele,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(wl1,6*nele,mpi_real8,0,mpi_comm_world, mpierr)

  if (Par%do_cross) then
     if (myid .ne. 0) then
        allocate(clEEobs2(0:Par%ellmax),clBBobs2(0:Par%ellmax))
        allocate(bl2(0:Par%ellmax+Par%Lmax,3),wl2(0:Par%ellmax+Par%Lmax,6))
     endif
  call mpi_bcast(clEEobs2,Par%ellmax+1,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(clBBobs2,Par%ellmax+1,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(bl2,3*nele,mpi_real8,0,mpi_comm_world, mpierr)
  call mpi_bcast(wl2,6*nele,mpi_real8,0,mpi_comm_world, mpierr)
  endif

  call mpi_barrier(mpi_comm_world, mpierr)
  if (Par%compute_alphalm .or. Par%compute_biasalpha)  then
     call mpi_bcast(Par%nsims,1,mpi_integer,0,mpi_comm_world,mpierr)
     call mpi_bcast(Par%ssim,1,mpi_integer,0,mpi_comm_world,mpierr)
     call mpi_bcast(Par%inmapfile1,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
     call mpi_bcast(Par%endnamemap1,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
     call mpi_bcast(Par%zerofill,1,mpi_integer,0,mpi_comm_world,mpierr)
     if (Par%do_cross) then
        call mpi_bcast(Par%inmapfile2,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
        call mpi_bcast(Par%endnamemap2,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
     endif

     if ((myid .eq. 0) .and. (Par%feedback .gt. 1)) write(*,*) 'Read maps'
     nele = Par%ellmax+1
     allocate(almE1(Par%nsims,0:Par%ellmax,0:Par%ellmax))
     allocate(almB1(Par%nsims,0:Par%ellmax,0:Par%ellmax))     
     if (Par%do_cross) then
        allocate(almE2(Par%nsims,0:Par%ellmax,0:Par%ellmax))
        allocate(almB2(Par%nsims,0:Par%ellmax,0:Par%ellmax))
     endif

     !mask
     if (myid .eq. 0) then
        if (Par%inmaskfile1 .ne. '') then
           apply_mask1 = .true.
           npix = getsize_fits(trim(Par%inmaskfile1))
           allocate(mask1(0:npix-1,1:3)) 
           call read_mask_and_compute_fsky(Par%inmaskfile1,mask1,fsky1)
        endif
        if (Par%do_cross .and. (Par%inmaskfile2 .ne. '')) then
           apply_mask2 = .true.
           npix = getsize_fits(trim(Par%inmaskfile2))
           allocate(mask2(0:npix-1,1:3))     
           call read_mask_and_compute_fsky(Par%inmaskfile2,mask2,fsky2)
        endif
     endif

     call mpi_bcast(apply_mask1,1,mpi_logical,0,mpi_comm_world,mpierr)
     if (apply_mask1) call mpi_bcast(fsky1,1,mpi_real8,0,mpi_comm_world, mpierr)
     if (Par%do_cross) call mpi_bcast(apply_mask2,1,mpi_logical,0,mpi_comm_world,mpierr)
     if (Par%do_cross .and. apply_mask2) call mpi_bcast(fsky2,1,mpi_real8,0,mpi_comm_world, mpierr)

     if (Par%nsims .eq. 1) then
        if (myid .eq. 0) then
           if (apply_mask1) then
              call read_map_and_compute_alms(Par%inmapfile1,Par%niter,almE1,almB1,1,mask=mask1)
           else
              call read_map_and_compute_alms(Par%inmapfile1,Par%niter,almE1,almB1,1)
           endif
        endif
        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_bcast(almE1,Par%nsims*nele*nele,mpi_complex,0,mpi_comm_world, mpierr)
        call mpi_bcast(almB1,Par%nsims*nele*nele,mpi_complex,0,mpi_comm_world, mpierr)       
        if (Par%do_cross) then
           if (myid .eq. 0) then
              if (apply_mask2) then
                 call read_map_and_compute_alms(Par%inmapfile2,Par%niter,almE2,almB2,1,mask=mask2)
              else
                 call read_map_and_compute_alms(Par%inmapfile2,Par%niter,almE2,almB2,1)
              endif
           endif
           call mpi_barrier(mpi_comm_world, mpierr)
           call mpi_bcast(almE2,Par%nsims*nele*nele,mpi_complex,0,mpi_comm_world,mpierr)
           call mpi_bcast(almB2,Par%nsims*nele*nele,mpi_complex,0,mpi_comm_world,mpierr)
        endif
     else
       almE1=0
       almB1=0
       ct=1
       write(strzerofill,fmt='(i1)') Par%zerofill
       do isim=Par%ssim,Par%ssim+Par%nsims-1
          if (myid .eq. mod(ct,nproc)) then
             if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' reading map ', isim
             write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
             mapname=trim(Par%inmapfile1)//trim(simstr)//trim(Par%endnamemap1)
             if (apply_mask1) then
                call read_map_and_compute_alms(mapname,Par%niter,almE1,almB1,ct,mask=mask1)
             else
                call read_map_and_compute_alms(mapname,Par%niter,almE1,almB1,ct)
             endif
          endif
          ct=ct+1
       enddo
       call mpi_barrier(mpi_comm_world, mpierr)
       call mpi_allreduce(mpi_in_place,almE1,Par%nsims*nele*nele,mpi_complex,mpi_sum,mpi_comm_world,mpierr)
       call mpi_allreduce(mpi_in_place,almB1,Par%nsims*nele*nele,mpi_complex,mpi_sum,mpi_comm_world,mpierr)
       if (Par%do_cross) then
          almE2=0
          almB2=0
          ct=1
          do isim=Par%ssim,Par%ssim+Par%nsims-1
             if (myid .eq. mod(ct,nproc)) then
                if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' reading map ', isim
                write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
                mapname=trim(Par%inmapfile2)//trim(simstr)//trim(Par%endnamemap2)
                if (apply_mask2) then
                   call read_map_and_compute_alms(mapname,Par%niter,almE2,almB2,ct,mask=mask2)
                else
                   call read_map_and_compute_alms(mapname,Par%niter,almE2,almB2,ct)
                endif
             endif
             ct=ct+1
          enddo
          call mpi_barrier(mpi_comm_world, mpierr)
          call mpi_allreduce(mpi_in_place,almE2,Par%nsims*nele*nele,mpi_complex,mpi_sum,mpi_comm_world,mpierr)
          call mpi_allreduce(mpi_in_place,almB2,Par%nsims*nele*nele,mpi_complex,mpi_sum,mpi_comm_world,mpierr)
       endif
    endif
  endif

  if (Par%compute_biasalpha) then
     if (Par%do_cross) then
        allocate(clEEmap12(Par%nsims,0:Par%ellmax))
        allocate(clBBmap12(Par%nsims,0:Par%ellmax))
        if (myid .eq. 0) call compute_cls_from_alms(almE1,almB1,clEEmap12,clBBmap12,almE2,almB2)
        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_bcast(clEEmap12,Par%nsims*nele,mpi_real8,0,mpi_comm_world,mpierr)
        call mpi_bcast(clBBmap12,Par%nsims*nele,mpi_real8,0,mpi_comm_world,mpierr)
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
        if  (myid .eq. 0) call compute_cls_from_alms(almE1,almB1,clEEmap1,clBBmap1)
        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_bcast(clEEmap1,Par%nsims*nele,mpi_real8,0,mpi_comm_world,mpierr)
        call mpi_bcast(clBBmap1,Par%nsims*nele,mpi_real8,0,mpi_comm_world,mpierr)
        if (apply_mask1) then
           clEEmap1 = clEEmap1 / fsky1
           clBBmap1 = clBBmap1 / fsky1
        endif
     endif
  endif

  if (Par%compute_biasalpha .and.(.not. Par%compute_alphalm)) then
    deallocate(almE1,almB1)
    if (Par%do_cross) deallocate(almE2,almB2)
  endif
  
  call mpi_barrier(mpi_comm_world, mpierr)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !! Start E operator computation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'Starting Computation'
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock(t1, clock_rate, clock_max)
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write (0,*) 'Elapsed real time for initialization = ', real( t1 - t0 ) / real ( clock_rate )
  !! Compute sigma
  
  allocate(one_o_var1(Par%Lmin:Par%Lmax))
  one_o_var1 = 0.0

  if (Par%do_cross) then
    allocate(one_o_var2(Par%Lmin:Par%Lmax))
    one_o_var2 = 0.0
  endif  

  if (Par%compute_biasalpha) then
     allocate(biasalpha(1:Par%nsims,Par%Lmin:Par%Lmax))
     biasalpha = 0
  endif

  if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'Computing sigma and bias'
  Lcount = 0
  do iL=Par%Lmin,Par%Lmax
     if (myid .eq. mod(Lcount,nproc)) then
        if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' doing multipole ', iL
        do iell=Par%ellmin,Par%ellmax
           allocate(wig2(iL+iell+1))
           call Wigner3j(wig2, jmin, jmax, iell, iL, -2, 2, 0)
           allocate(F_EB1(jmin:jmax),F_BE1(jmin:jmax))
           F_EB1 = 2 * wig2(1:jmax-jmin+1) * clEEfid(iell)*bl1(iell,1)*bl1(jmin:jmax,1)
           F_BE1 = 2 * wig2(1:jmax-jmin+1) * clEEfid(jmin:jmax)*bl1(jmin:jmax,1)*bl1(iell,1)
           if (Par%do_cross) then
              allocate(F_EB2(jmin:jmax),F_BE2(jmin:jmax))
              F_EB2 = 2 * wig2(1:jmax-jmin+1) * clEEfid(iell)*bl2(iell,1)*bl2(jmin:jmax,1)
              F_BE2 = 2 * wig2(1:jmax-jmin+1) * clEEfid(jmin:jmax)*bl2(jmin:jmax,1)*bl2(iell,1)
           endif
           Gl = (2*iell + 1)/FOURPI
           do j = jmin,min(jmax,Par%ellmax)
              if ((j .ge. iell) .and. (mod(iL+iell+j,2).eq.0)) then
                 if (j .eq. iell) then
                    factor = 0.5 * Gl * (2.0*j + 1.0)
                 else 
                    factor = Gl * (2.0*j + 1.0)
                 endif
                 one_o_var1(iL) = one_o_var1(iL) + factor * &
                       (F_EB1(j)**2/clEEobs1(iell)/clBBobs1(j) + &
                        F_BE1(j)**2/clBBobs1(iell)/clEEobs1(j))                  
                 if (Par%do_cross) then
                    one_o_var2(iL) = one_o_var2(iL) + factor * &
                         (F_EB2(j)**2/clEEobs2(iell)/clBBobs2(j) + &
                          F_BE2(j)**2/clBBobs2(iell)/clEEobs2(j))
                 endif
                 !bias computation
                 if (Par%compute_biasalpha) then
                    if (j .eq. iell) then
                       factor = 0.25 * Gl * (2.0*j + 1.0)
                    else
                       factor = Gl * (2.0*j + 1.0)
                    endif
                    if (Par%do_cross) then 
                        biasalpha(:,iL) = biasalpha(:,iL) + factor * &
                          (F_EB1(j)*F_EB2(j) * clEEmap12(:,iell) * clBBmap12(:,j) / &
                             (clBBobs1(j) * clBBobs2(j) * clEEobs1(iell) * clEEobs2(iell)) + &
                              F_BE1(j)*F_BE2(j) * clEEmap12(:,j) * clBBmap12(:,iell) / &
                             (clBBobs1(iell) * clBBobs1(iell) * clEEobs2(j) * clEEobs2(j)))
                    else
                        biasalpha(:,iL) = biasalpha(:,iL) + factor * &
                          (F_EB1(j)**2 * clEEmap1(:,iell) * clBBmap1(:,j) / &
                             (clBBobs1(j) * clBBobs1(j) * clEEobs1(iell) * clEEobs1(iell)) + &
                              F_BE1(j)**2 * clEEmap1(:,j) * clBBmap1(:,iell) / &
                             (clBBobs1(iell) * clBBobs1(iell) * clEEobs1(j) * clEEobs1(j)))
                    endif 
                 endif
              endif
           enddo
           deallocate(wig2,F_EB1,F_BE1)
           if (Par%do_cross) deallocate(F_EB2,F_BE2)
        enddo
     endif
     Lcount = Lcount + 1
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

  if (myid .eq. 0) then
     allocate(red_one_o_var1(Par%Lmin:Par%Lmax))
     red_one_o_var1=0.0
  endif
  call mpi_barrier(mpi_comm_world, mpierr)
  call mpi_reduce(one_o_var1,red_one_o_var1,Par%Lmax-Par%Lmin+1,mpi_real8,mpi_sum,0,mpi_comm_world, mpierr)
  deallocate(one_o_var1)
  call mpi_barrier(mpi_comm_world, mpierr)

  if (Par%do_cross) then 
     if (myid .eq. 0) then
        allocate(red_one_o_var2(Par%Lmin:Par%Lmax))
        red_one_o_var2=0.0
     endif
     call mpi_barrier(mpi_comm_world, mpierr)
     call mpi_reduce(one_o_var2,red_one_o_var2,Par%Lmax-Par%Lmin+1,mpi_real8,mpi_sum,0,mpi_comm_world,mpierr)
     deallocate(one_o_var2)
     call mpi_barrier(mpi_comm_world, mpierr)  
  endif
  
  if (myid .eq. 0) then
     if (Par%feedback .gt. 1) write(0,*) 'Write out sigma'
     open(newunit=myunit,file=trim(Par%outsigmafile1),status='replace',form='formatted')
     do iL=Par%Lmin,Par%Lmax
        write(myunit,'(I4,*(E15.7))') iL,sqrt(1.0/red_one_o_var1(iL))
     enddo
     close(myunit)
     if (Par%do_cross) then
        open(newunit=myunit,file=trim(Par%outsigmafile2),status='replace',form='formatted')
        do iL=Par%Lmin,Par%Lmax
           write(myunit,'(I4,*(E15.7))') iL,sqrt(1.0/red_one_o_var2(iL))
        enddo
        close(myunit)
     endif
  endif

  call mpi_barrier(mpi_comm_world, mpierr) 

  if (Par%compute_biasalpha) then
     if (myid .eq. 0) then
        allocate(red_biasalpha(1:Par%nsims,Par%Lmin:Par%Lmax))
        red_biasalpha=0.0
     endif
     call mpi_barrier(mpi_comm_world, mpierr)
     call mpi_reduce(biasalpha,red_biasalpha,(Par%Lmax-Par%Lmin+1)*Par%nsims,mpi_real8,mpi_sum,0,mpi_comm_world,mpierr)
     deallocate(biasalpha)
  endif

  call mpi_barrier(mpi_comm_world, mpierr)

  if (Par%compute_biasalpha .and. (myid .eq. 0)) then
     if (Par%do_cross) then
        do iL=Par%Lmin,Par%Lmax
           red_biasalpha(:,iL) = red_biasalpha(:,iL) / red_one_o_var1(iL) / red_one_o_var2(iL)
        enddo
     else
        do iL=Par%Lmin,Par%Lmax
           red_biasalpha(:,iL) = red_biasalpha(:,iL) / red_one_o_var1(iL)**2
        enddo
     endif
     if (Par%feedback .gt. 1) write(0,*) 'Write out bias'
     call write_out_cls(Par%outbiasfile,Par%ssim,Par%zerofill,Par%endnamebias,red_biasalpha,Par%Lmin)
  endif

  !! Compute alphalm
  if (Par%compute_alphalm) then

     if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'Computing alpha_LM'
      
     call mpi_barrier(mpi_comm_world, mpierr)
     allocate(almalpha1(1:Par%nsims,0:Par%Lmax,0:Par%Lmax))
     almalpha1=0
     allocate(curralmE1(Par%nsims),curralmB1(Par%nsims),curralmE1star(Par%nsims),curralmB1star(Par%nsims))

     if (Par%do_cross) then
        allocate(almalpha2(1:Par%nsims,0:Par%Lmax,0:Par%Lmax))
        almalpha2=0
        allocate(curralmE2(Par%nsims),curralmB2(Par%nsims),curralmE2star(Par%nsims),curralmB2star(Par%nsims))
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
                 allocate(F_EB1(jmin:jmax),F_BE1(jmin:jmax))
                 F_EB1 = 2 * wig2(1:jmax-jmin+1) * clEEfid(iell)*bl1(iell,1)*bl1(jmin:jmax,1)
                 F_BE1 = 2 * wig2(1:jmax-jmin+1) * clEEfid(jmin:jmax)*bl1(jmin:jmax,1)*bl1(iell,1)
                 if (Par%do_cross) then
                    allocate(F_EB2(jmin:jmax),F_BE2(jmin:jmax))
                    F_EB2 = 2 * wig2(1:jmax-jmin+1) * clEEfid(iell)*bl2(iell,1)*bl2(jmin:jmax,1)
                    F_BE2 = 2 * wig2(1:jmax-jmin+1) * clEEfid(jmin:jmax)*bl2(jmin:jmax,1)*bl2(iell,1)
                 endif
                 norm = sqrt((2.0*iell + 1.0)*(2.0*iL + 1.0)/FOURPI)
                 !loop emm
                 do iemm=-iell,iell
                    iemmp = iemm-iM
                    call Wigner3j(wigall, jminall, jmaxall, iell, iL, iemmp , -iemm, iM)

                    if (iemm .ge. 0) then
                       curralmE1 = almE1(:,iell,iemm)
                       curralmB1 = almB1(:,iell,iemm)
                    else 
                       curralmE1 = conjg(almE1(:,iell,-iemm))
                       curralmB1 = conjg(almB1(:,iell,-iemm))
                    endif

                    if (Par%do_cross) then
                       if (iemm .ge. 0) then
                          curralmE2 = almE2(:,iell,iemm)
                          curralmB2 = almB2(:,iell,iemm)
                       else
                          curralmE2 = conjg(almE2(:,iell,-iemm))
                          curralmB2 = conjg(almB2(:,iell,-iemm))
                       endif
                    endif

                    emmppos = .true.
                    if (iemmp .lt. 0) emmppos=.false.

                    allocate(csi(jminall:jmaxall))
                    do j=jminall,jmaxall
                       csi(j) = sqrt(2.0*j + 1.0) * wigall(j-jminall+1)
                    enddo
                    csi = csi * norm * -1.0**iemm

                    !loop ell'
                    do j = max(jminall,jmin,Par%ellmin),min(jmaxall,Par%ellmax)
                       if ((j .ge. iell) .and. (mod(iL+iell+j,2) .eq. 0)) then
                          if (j .eq. iell) then
                             factor = 0.5
                          else
                             factor = 1.0
                          endif

                          if (emmppos) then
                             curralmE1star = conjg(almE1(:,j,iemmp))
                             curralmB1star = conjg(almB1(:,j,iemmp))
                          else
                             curralmE1star = almE1(:,j,-iemmp)
                             curralmB1star = almB1(:,j,-iemmp)
                          endif
                          almalpha1(:,iL,iM) = almalpha1(:,iL,iM) + factor * &
                             (F_EB1(j) * csi(j) * curralmE1 * curralmB1star / clBBobs1(j) / clEEobs1(iell) + &
                              F_BE1(j) * csi(j) * curralmB1 * curralmE1star / clBBobs1(iell) / clEEobs1(j))

                          if (Par%do_cross) then
                             if (emmppos) then
                                curralmE2star = conjg(almE2(:,j,iemmp))
                                curralmB2star = conjg(almB2(:,j,iemmp))
                             else
                                curralmE2star = almE2(:,j,-iemmp)
                                curralmB2star = almB2(:,j,-iemmp)
                             endif                    
                             almalpha2(:,iL,iM) = almalpha2(:,iL,iM) + factor * &
                                (F_EB2(j) * csi(j) * curralmE2 * curralmB2star / clBBobs2(j) / clEEobs2(iell) + &
                                 F_BE2(j) * csi(j) * curralmB2 * curralmE2star / clBBobs2(iell) / clEEobs2(j)) 
                          endif
                       endif
                    enddo
                    deallocate(csi)
                 enddo
                 deallocate(wig2,wigall,F_EB1,F_BE1)
                 if (Par%do_cross) deallocate(F_EB2,F_BE2)
              enddo
           endif
           elementcount = elementcount + 1
        enddo
     enddo
     deallocate(curralmE1,curralmB1,curralmE1star,curralmB1star)
     if (Par%do_cross) deallocate(curralmE2,curralmB2,curralmE2star,curralmB2star)
     deallocate(almE1,almB1)
     if (Par%do_cross) deallocate(almE2,almB2) 
  
     call mpi_barrier(mpi_comm_world, mpierr)
     if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'Computation of alphalm done'
     if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock( t3, clock_rate, clock_max )
     if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write(0,*) 'Elapsed real time for alphalm computation = ', real ( t3 - t2 ) / real ( clock_rate )

     call mpi_barrier(mpi_comm_world, mpierr)
     
     if (myid .eq. 0) then
        allocate(red_almalpha1(1:Par%nsims,0:Par%Lmax,0:Par%Lmax))
        red_almalpha1=0.0
     endif

     call mpi_barrier(mpi_comm_world, mpierr)
     nele = Par%Lmax+1
     call mpi_reduce(almalpha1,red_almalpha1,nele*nele*Par%nsims,mpi_complex,mpi_sum,0,mpi_comm_world,mpierr)
     deallocate(almalpha1)
     call mpi_barrier(mpi_comm_world, mpierr)

     if (Par%do_cross) then
        if (myid .eq. 0) then
           allocate(red_almalpha2(1:Par%nsims,0:Par%Lmax,0:Par%Lmax))
           red_almalpha2=0.0
        endif

        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_reduce(almalpha2,red_almalpha2,nele*nele*Par%nsims,mpi_complex,mpi_sum,0,mpi_comm_world,mpierr)
        deallocate(almalpha2)
        call mpi_barrier(mpi_comm_world, mpierr)
     endif
     
     if (myid .eq. 0) then
        do iL=Par%Lmin,Par%Lmax
           red_almalpha1(:,iL,:) = red_almalpha1(:,iL,:) / red_one_o_var1(iL)
        enddo
        if (Par%do_cross) then
           do iL=Par%Lmin,Par%Lmax
              red_almalpha2(:,iL,:) = red_almalpha2(:,iL,:) / red_one_o_var2(iL)
           enddo
        endif
        if (Par%feedback .gt. 1) write(0,*) 'Write out alphalm'
        call write_out_alms(Par%outalmfile1,Par%ssim,Par%zerofill,Par%endnamealm1,red_almalpha1)
        if (Par%do_cross) call write_out_alms(Par%outalmfile2,Par%ssim,Par%zerofill,Par%endnamealm2,red_almalpha2)
        if (Par%compute_alphacl) then
           if (Par%do_cross) then
              if (Par%compute_biasalpha .and. Par%subtract_bias) then
                 call compute_and_write_cl(Par%outclfile,Par%ssim,Par%zerofill,Par%endnamecl,red_almalpha1,Par%Lmin,bias=red_biasalpha,alms2=red_almalpha2)
              else
                 call compute_and_write_cl(Par%outclfile,Par%ssim,Par%zerofill,Par%endnamecl,red_almalpha1,Par%Lmin,alms2=red_almalpha2)
              endif
           else
              if (Par%compute_biasalpha .and. Par%subtract_bias) then
                 call compute_and_write_cl(Par%outclfile,Par%ssim,Par%zerofill,Par%endnamecl,red_almalpha1,Par%Lmin,bias=red_biasalpha)
              else
                 call compute_and_write_cl(Par%outclfile,Par%ssim,Par%zerofill,Par%endnamecl,red_almalpha1,Par%Lmin)
              endif

           endif
        endif
        deallocate(red_almalpha1,red_one_o_var1)
        if (Par%do_cross) deallocate(red_almalpha2,red_one_o_var2)
        if (Par%compute_biasalpha) deallocate(red_biasalpha)
     endif
  else
     if (myid .eq. 0) deallocate(red_one_o_var1)
     if ((myid .eq. 0) .and. Par%do_cross) deallocate(red_one_o_var2)
  endif
  
  call mpi_barrier(mpi_comm_world, mpierr)
  
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock( t4, clock_rate, clock_max )
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write (0,*) 'Elapsed real time for total computation = ', real ( t4 - t1 ) / real ( clock_rate )
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write (0,*) 'Elapsed real time for computation and initialization = ', real ( t4 - t0 ) / real ( clock_rate )
  if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'End Program'  
  call mpi_finalize(mpierr)
  
end program EB_estimator
