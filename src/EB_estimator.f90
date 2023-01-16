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
  integer :: indmax,indLMmax,ind,indLM,LM(2)
  integer :: iellp,ellpmin,ellpmax, iell
  integer :: ellpminall,ellpmaxall,usedellpmax
  integer :: isim,Lcount,iL,nele,elementcount,procelementcount
  integer :: iemm, iM, iemmp
  integer :: selellpmin,selellpmax 
  integer(i8b) :: npix
  logical :: emmppos,apply_mask1=.false.,apply_mask2=.false.
  real(dp), allocatable, dimension(:) :: one_o_var1,red_one_o_var1,one_o_var2,red_one_o_var2,wig2,wigall
  real(dp), allocatable, dimension(:) :: F_EB1,F_BE1,clEEfid,clEEobs1,clBBobs1
  real(dp), allocatable, dimension(:) :: F_EB2,F_BE2,clEEobs2,clBBobs2
  real(dp), allocatable, dimension(:,:) :: clEEmap1,clBBmap1,clEEmap12,clBBmap12,biasalpha,red_biasalpha
  complex(spc), allocatable, dimension(:,:) :: almE1,almB1,almE2,almB2
  complex(spc), allocatable, dimension(:,:,:) :: alm1,alm2
  complex(spc), allocatable, dimension(:,:) :: almalpha1,red_almalpha1,almalpha2,red_almalpha2
  complex(spc), allocatable, dimension(:,:) :: procalmalpha1,procalmalpha2
  complex(dpc), allocatable, dimension(:) :: curralmE1,curralmB1,curralmE1star,curralmB1star
  complex(dpc), allocatable, dimension(:) :: curralmE2,curralmB2,curralmE2star,curralmB2star
  real(dp), allocatable,dimension(:,:) :: clfid,wl1,bl1,nl1,wl2,bl2,nl2
  real(sp), allocatable,dimension(:,:) :: mask1,mask2
  real(dp) :: factor,Gl,norm,fsky1,fsky2,m1_to_emm,m1_to_memm,m1_to_memmp,abswigmax 
  integer :: t0,t1,t2,t3,t4,clock_rate,clock_max,myunit,ct
  character(len=FILENAMELEN) :: mapname
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
  call mpi_bcast(Par%Lmin,1,mpi_integer,0,mpi_comm_world, mpierr)
  call mpi_bcast(Par%Lmax,1,mpi_integer,0,mpi_comm_world, mpierr)
  call mpi_bcast(Par%ellmin,1,mpi_integer,0,mpi_comm_world, mpierr)
  call mpi_bcast(Par%ellmax,1,mpi_integer,0,mpi_comm_world, mpierr)
  call mpi_bcast(Par%compute_alphalm,1,mpi_logical,0,mpi_comm_world,mpierr)
  call mpi_bcast(Par%compute_biasalpha,1,mpi_logical,0,mpi_comm_world,mpierr)
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
     if (Par%do_cross) then
        call mpi_bcast(Par%inmapfile2,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
        call mpi_bcast(Par%endnamemap2,FILENAMELEN,mpi_character,0,mpi_comm_world,mpierr)
     endif
     call mpi_bcast(Par%read_precomputed_alms,1,mpi_logical,0,mpi_comm_world,mpierr)

     indmax = lm2index(Par%ellmax,Par%ellmax,Par%ellmax) 
     nele = indmax+1
     allocate(almE1(Par%nsims,0:indmax))
     allocate(almB1(Par%nsims,0:indmax))     
     almE1=0
     almB1=0
     
     if (Par%do_cross) then
        allocate(almE2(Par%nsims,0:indmax))
        allocate(almB2(Par%nsims,0:indmax))
        almE2=0
        almB2=0
     endif

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
        if (.not. Par%read_precomputed_alms) then
           call mpi_bcast(npix,1,mpi_integer8,0,mpi_comm_world,mpierr)
           if (myid .ne. 0) allocate(mask1(0:npix-1,1:3))
           call mpi_bcast(mask1,npix*3,mpi_real4,0,mpi_comm_world,mpierr)
        else
           if (myid .eq. 0) deallocate(mask1)
        endif
        call mpi_bcast(fsky1,1,mpi_real8,0,mpi_comm_world,mpierr)
     endif
     if (Par%do_cross) then
        call mpi_bcast(apply_mask2,1,mpi_logical,0,mpi_comm_world,mpierr)
        if (apply_mask2) then
           if (.not. Par%read_precomputed_alms) then
              call mpi_bcast(npix,1,mpi_integer8,0,mpi_comm_world,mpierr)
              if (myid .ne. 0) allocate(mask2(0:npix-1,1:3))
              call mpi_bcast(mask2,npix*3,mpi_real4,0,mpi_comm_world,mpierr)
           else
              if (myid .eq. 0) deallocate(mask2)
           endif
           call mpi_bcast(fsky2,1,mpi_real8,0,mpi_comm_world,mpierr)
        endif
     endif
     call mpi_barrier(mpi_comm_world, mpierr)

     if ((myid .eq. 0) .and. (Par%feedback .gt. 1)) write(*,*) 'Read maps'
     if (Par%nsims .eq. 1) then
        if (myid .eq. 0) then
           if (.not. Par%read_precomputed_alms) then
              if (apply_mask1) then
                 call read_map_and_compute_alms(Par%inmapfile1,Par%niter,almE1,almB1,Par%ellmax,1,mask=mask1)
              else
                 call read_map_and_compute_alms(Par%inmapfile1,Par%niter,almE1,almB1,Par%ellmax,1)
              endif
           else
              call read_precomputed_alms(Par%inmapfile1,almE1,almB1,1)
           endif
        endif
        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_bcast(almE1,nele,mpi_complex,0,mpi_comm_world, mpierr)
        call mpi_bcast(almB1,nele,mpi_complex,0,mpi_comm_world, mpierr)       
        if (Par%do_cross) then
           if (myid .eq. 0) then
              if (.not. Par%read_precomputed_alms) then
                 if (apply_mask2) then
                    call read_map_and_compute_alms(Par%inmapfile2,Par%niter,almE2,almB2,Par%ellmax,1,mask=mask2)
                 else
                    call read_map_and_compute_alms(Par%inmapfile2,Par%niter,almE2,almB2,Par%ellmax,1)
                 endif
              else
                 call read_precomputed_alms(Par%inmapfile2,almE2,almB2,1)
              endif
           endif
           call mpi_barrier(mpi_comm_world, mpierr)
           call mpi_bcast(almE2,nele,mpi_complex,0,mpi_comm_world,mpierr)
           call mpi_bcast(almB2,nele,mpi_complex,0,mpi_comm_world,mpierr)
        endif
     else
        ct=1
        write(strzerofill,fmt='(i1)') Par%zerofill
        do isim=Par%ssim,Par%ssim+Par%nsims-1
           if (myid .eq. mod(ct-1,nproc)) then
              if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' reading map ', isim
              write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
              mapname=trim(Par%inmapfile1)//trim(simstr)//trim(Par%endnamemap1)
              if (.not. Par%read_precomputed_alms) then
                 if (apply_mask1) then
                    call read_map_and_compute_alms(mapname,Par%niter,almE1,almB1,Par%ellmax,ct,mask=mask1)
                 else
                    call read_map_and_compute_alms(mapname,Par%niter,almE1,almB1,Par%ellmax,ct)
                 endif
              else
                 call read_precomputed_alms(mapname,almE1,almB1,ct)
              endif
           endif
           ct=ct+1
        enddo
        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_allreduce(mpi_in_place,almE1,Par%nsims*nele,mpi_complex,mpi_sum,mpi_comm_world,mpierr)
        call mpi_allreduce(mpi_in_place,almB1,Par%nsims*nele,mpi_complex,mpi_sum,mpi_comm_world,mpierr)

        call mpi_barrier(mpi_comm_world, mpierr)     
        
        if (Par%do_cross) then
           ct=1
           do isim=Par%ssim,Par%ssim+Par%nsims-1
              if (myid .eq. mod(ct-1,nproc)) then
                 if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' reading map ', isim
                 write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
                 mapname=trim(Par%inmapfile2)//trim(simstr)//trim(Par%endnamemap2)
                 if (.not. Par%read_precomputed_alms) then
                    if (apply_mask2) then
                       call read_map_and_compute_alms(mapname,Par%niter,almE2,almB2,Par%ellmax,ct,mask=mask2)
                    else
                       call read_map_and_compute_alms(mapname,Par%niter,almE2,almB2,Par%ellmax,ct)
                    endif
                 else
                    call read_precomputed_alms(mapname,almE2,almB2,ct)
                 endif
              endif
              ct=ct+1
           enddo
           call mpi_barrier(mpi_comm_world, mpierr)
           call mpi_allreduce(mpi_in_place,almE2,Par%nsims*nele,mpi_complex,mpi_sum,mpi_comm_world,mpierr)
           call mpi_allreduce(mpi_in_place,almB2,Par%nsims*nele,mpi_complex,mpi_sum,mpi_comm_world,mpierr)
        endif
     endif
     if (apply_mask1 .and. (.not. Par%read_precomputed_alms)) deallocate(mask1)
     if (apply_mask2 .and. (.not. Par%read_precomputed_alms)) deallocate(mask2)
  endif

  if (Par%compute_biasalpha) then
     if (Par%do_cross) then
        allocate(clEEmap12(Par%nsims,0:Par%ellmax))
        allocate(clBBmap12(Par%nsims,0:Par%ellmax))
        if (myid .eq. 0) call compute_cls_from_alms(almE1,almB1,Par%ellmax,clEEmap12,clBBmap12,almE2,almB2)
        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_bcast(clEEmap12,Par%nsims*(Par%ellmax+1),mpi_real8,0,mpi_comm_world,mpierr)
        call mpi_bcast(clBBmap12,Par%nsims*(Par%ellmax+1),mpi_real8,0,mpi_comm_world,mpierr)
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
        if  (myid .eq. 0) call compute_cls_from_alms(almE1,almB1,Par%ellmax,clEEmap1,clBBmap1)
        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_bcast(clEEmap1,Par%nsims*(Par%ellmax+1),mpi_real8,0,mpi_comm_world,mpierr)
        call mpi_bcast(clBBmap1,Par%nsims*(Par%ellmax+1),mpi_real8,0,mpi_comm_world,mpierr)
        if (apply_mask1) then
           clEEmap1 = clEEmap1 / fsky1
           clBBmap1 = clBBmap1 / fsky1
        endif
     endif
  endif

  if (Par%compute_alphalm) indLMmax = lm2index(Par%Lmax,Par%Lmax,Par%Lmax)

  if (Par%compute_biasalpha .and. (.not. Par%compute_alphalm)) then
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

  call mpi_barrier(mpi_comm_world, mpierr)
  allocate(red_one_o_var1(Par%Lmin:Par%Lmax))
  red_one_o_var1=0.0
  call mpi_barrier(mpi_comm_world, mpierr)
  call mpi_reduce(one_o_var1,red_one_o_var1,Par%Lmax-Par%Lmin+1,mpi_real8,mpi_sum,0,mpi_comm_world, mpierr)
  deallocate(one_o_var1)
  call mpi_barrier(mpi_comm_world, mpierr)

  if (Par%do_cross) then 
     allocate(red_one_o_var2(Par%Lmin:Par%Lmax))
     red_one_o_var2=0.0
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
     allocate(red_biasalpha(1:Par%nsims,Par%Lmin:Par%Lmax))
     red_biasalpha=0.0
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

     !! count the alms alpha and allocate local variables 
     elementcount = 0
     procelementcount = 0 
     do iL=Par%Lmin,Par%Lmax
        do iM=0,iL
           if (myid .eq. mod(elementcount,nproc)) then
              procelementcount = procelementcount + 1
           endif
           elementcount = elementcount + 1
        enddo
     enddo
     allocate(procalmalpha1(1:Par%nsims,procelementcount))
     procalmalpha1 = 0

     if (Par%do_cross) then
        allocate(procalmalpha2(1:Par%nsims,procelementcount))
        procalmalpha2=0
     endif

     allocate(curralmE1(Par%nsims),curralmB1(Par%nsims),curralmE1star(Par%nsims),curralmB1star(Par%nsims))
     if (Par%do_cross) allocate(curralmE2(Par%nsims),curralmB2(Par%nsims),curralmE2star(Par%nsims),curralmB2star(Par%nsims))

     call mpi_barrier(mpi_comm_world, mpierr)
   
     elementcount = 0
     procelementcount = 0
     do iL=Par%Lmin,Par%Lmax
        do iM=0,iL
           if (myid .eq. mod(elementcount,nproc)) then
              if (Par%feedback .gt. 3) write(0,*) 'Proc ', myid,' doing L=', iL, ' and M=',iM
              procelementcount = procelementcount + 1
              indLM = lm2index(Par%Lmax,iL,iM)
              !loop ell
              do iell=Par%ellmin,Par%ellmax
                 allocate(wig2(iL+iell+1),wigall(iL+iell+1))
                 call Wigner3j(wig2, ellpmin, ellpmax, iell, iL, -2, 2, 0)
                 usedellpmax = min(ellpmax,Par%ellmax)
                 allocate(F_EB1(ellpmin:usedellpmax),F_BE1(ellpmin:usedellpmax))
                 ! F_XY divided here by clEEobs(l) or clBBobs(l). It saves some computation time 
                 F_EB1 = 2 * wig2(1:usedellpmax-ellpmin+1) * clEEfid(iell)*bl1(iell,myE)*bl1(ellpmin:usedellpmax,myE) / clEEobs1(iell) 
                 F_BE1 = 2 * wig2(1:usedellpmax-ellpmin+1) * clEEfid(ellpmin:usedellpmax)*bl1(ellpmin:usedellpmax,myE)*bl1(iell,myE) / clBBobs1(iell)
                 if (Par%do_cross) then
                    allocate(F_EB2(ellpmin:usedellpmax),F_BE2(ellpmin:usedellpmax))
                    F_EB2 = 2 * wig2(1:usedellpmax-ellpmin+1) * clEEfid(iell)*bl2(iell,myE)*bl2(ellpmin:usedellpmax,myE) / clEEobs2(iell) 
                    F_BE2 = 2 * wig2(1:usedellpmax-ellpmin+1) * clEEfid(ellpmin:usedellpmax)*bl2(ellpmin:usedellpmax,myE)*bl2(iell,myE) / clBBobs2(iell)
                 endif
                 deallocate(wig2)
                 norm = sqrt((2.0*iell + 1.0)*(2.0*iL + 1.0)/FOURPI)
                 !loop emm
                 do iemm=-iell,iell
                    iemmp = iemm-iM
                    m1_to_memm=(-1.0)**(-iemm)
                    m1_to_emm=(-1.0)**iemm

                    call Wigner3j(wigall, ellpminall, ellpmaxall, iell, iL, iemmp , -iemm, iM)
!                    call DRC3JJ(real(iell,kind=dp), real(iL,kind=dp), real(-iemm,kind=dp), real(iM,kind=dp), &
!                           ellpminall, ellpmaxall, wigall, iL+iell+1,IER)
                    if (iemm .ge. 0) then
                       ind = lm2index(Par%ellmax,iell,iemm)
                       curralmE1 = almE1(:,ind)
                       curralmB1 = almB1(:,ind)
                    else
                       ind = lm2index(Par%ellmax,iell,-iemm)
                       curralmE1 = conjg(almE1(:,ind))*m1_to_memm
                       curralmB1 = conjg(almB1(:,ind))*m1_to_memm
                    endif

                    if (Par%do_cross) then
                       if (iemm .ge. 0) then
                          !ind is the same of map1
                          curralmE2 = almE2(:,ind)
                          curralmB2 = almB2(:,ind)
                       else
                          !ind is the same of map1
                          curralmE2 = conjg(almE2(:,ind))*m1_to_memm
                          curralmB2 = conjg(almB2(:,ind))*m1_to_memm
                       endif
                    endif

                    emmppos = .true.
                    if (iemmp .lt. 0) emmppos=.false.
                    
                    m1_to_memmp=(-1.0)**(-iemmp) 

                    selellpmin = max(ellpminall,ellpmin,Par%ellmin,iell)
                    selellpmax = min(ellpmaxall,usedellpmax)
                    abswigmax = maxval(abs(wigall(selellpmin-ellpminall+1:selellpmax-ellpminall+1)))
                    
                    do iellp = selellpmin,selellpmax
                       if ((jmod(iL+iell+iellp,2) .eq. 0) .and. (abs(wigall(iellp-ellpminall+1)) .gt. WIGRATIO*abswigmax)) then
                          if (iellp .eq. iell) then
                             factor = 0.5 * norm * m1_to_emm * sqrt(2.0*iellp + 1.0) * wigall(iellp-ellpminall+1)
                          else
                             factor = norm * m1_to_emm * sqrt(2.0*iellp + 1.0) * wigall(iellp-ellpminall+1)
                          endif

                          if (emmppos) then
                             ind = lm2index(Par%ellmax,iellp,iemmp)
                             curralmE1star = conjg(almE1(:,ind))
                             curralmB1star = conjg(almB1(:,ind))
                          else
                             ind = lm2index(Par%ellmax,iellp,-iemmp)
                             curralmE1star = almE1(:,ind)*m1_to_memmp
                             curralmB1star = almB1(:,ind)*m1_to_memmp
                          endif
                          ! F_XY already divided by clEEobs(l) or clBBobs(l)
                          procalmalpha1(:,procelementcount) = procalmalpha1(:,procelementcount) + factor * &
                             (F_EB1(iellp) * curralmE1 * curralmB1star / clBBobs1(iellp) + &
                              F_BE1(iellp) * curralmB1 * curralmE1star / clEEobs1(iellp))

                          if (Par%do_cross) then
                             if (emmppos) then
                                !ind is the same of map1
                                curralmE2star = conjg(almE2(:,ind))
                                curralmB2star = conjg(almB2(:,ind))
                             else
                                !ind is the same of map1
                                curralmE2star = almE2(:,ind)*m1_to_memmp
                                curralmB2star = almB2(:,ind)*m1_to_memmp
                             endif
                             ! F_XY already divided by clEEobs(l) or clBBobs(l)
                             procalmalpha2(:,procelementcount) = procalmalpha2(:,procelementcount) + factor * &
                                (F_EB2(iellp) * curralmE2 * curralmB2star / clBBobs2(iellp) + &
                                 F_BE2(iellp) * curralmB2 * curralmE2star / clEEobs2(iellp)) 
                          endif
                       endif
                    enddo
                 enddo
                 deallocate(wigall,F_EB1,F_BE1)
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

     !! reorder the alms alpha
     allocate(almalpha1(1:Par%nsims,0:indLMmax))
     if (Par%do_cross) allocate(almalpha2(1:Par%nsims,0:indLMmax))
     elementcount = 0
     procelementcount = 0
     do iL=Par%Lmin,Par%Lmax
        do iM=0,iL
           if (myid .eq. mod(elementcount,nproc)) then
              procelementcount = procelementcount + 1
              indLM = lm2index(Par%Lmax,iL,iM)
              almalpha1(:,indLM) = procalmalpha1(:,procelementcount)
              if (Par%do_cross) almalpha2(:,indLM) = procalmalpha2(:,procelementcount) 
           endif
           elementcount = elementcount + 1
        enddo
     enddo
     
     deallocate(procalmalpha1)
     if (Par%do_cross) deallocate(procalmalpha2)

     call mpi_barrier(mpi_comm_world, mpierr)
     if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'Computation of alphalm done'
     if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock( t3, clock_rate, clock_max )
     if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write(0,*) 'Elapsed real time for alphalm computation = ', real ( t3 - t2 ) / real ( clock_rate )

     call mpi_barrier(mpi_comm_world, mpierr)
     
     allocate(red_almalpha1(1:Par%nsims,0:indLMmax))
     red_almalpha1=0.0

     call mpi_barrier(mpi_comm_world, mpierr)
     nele = indLMmax+1
     call mpi_reduce(almalpha1,red_almalpha1,Par%nsims*nele,mpi_complex,mpi_sum,0,mpi_comm_world,mpierr)
     deallocate(almalpha1)
     call mpi_barrier(mpi_comm_world, mpierr)

     if (Par%do_cross) then
        allocate(red_almalpha2(1:Par%nsims,0:indLMmax))
        red_almalpha2=0.0

        call mpi_barrier(mpi_comm_world, mpierr)
        call mpi_reduce(almalpha2,red_almalpha2,Par%nsims*nele,mpi_complex,mpi_sum,0,mpi_comm_world,mpierr)
        deallocate(almalpha2)
        call mpi_barrier(mpi_comm_world, mpierr)
     endif
     
     if (myid .eq. 0) then
        allocate(alm1(1:Par%nsims,0:Par%Lmax,0:Par%Lmax))
        alm1 = 0
        do indLM=0,indLMmax
           LM = index2lm(Par%Lmax,indLM)
           alm1(:,lm(1),lm(2)) = red_almalpha1(:,indLM) 
        enddo
        deallocate(red_almalpha1)
        do iL=Par%Lmin,Par%Lmax
           alm1(:,iL,:) = alm1(:,iL,:) / red_one_o_var1(iL)
        enddo
        if (Par%do_cross) then
           allocate(alm2(1:Par%nsims,0:Par%Lmax,0:Par%Lmax))
           alm2 = 0
           do indLM=0,indLMmax
              LM = index2lm(Par%Lmax,indLM)
              alm2(:,lm(1),lm(2)) = red_almalpha2(:,indLM)
           enddo
           deallocate(red_almalpha2)
           do iL=Par%Lmin,Par%Lmax
              alm2(:,iL,:) = alm2(:,iL,:) / red_one_o_var2(iL)
           enddo
        endif
        if (Par%feedback .gt. 1) write(0,*) 'Write out alphalm'
        call write_out_alms(Par%outalmfile1,Par%ssim,Par%zerofill,Par%endnamealm1,alm1)
        if (Par%do_cross) call write_out_alms(Par%outalmfile2,Par%ssim,Par%zerofill,Par%endnamealm2,alm2)
        if (Par%compute_alphacl) then
           if (Par%do_cross) then
              if (Par%compute_biasalpha .and. Par%subtract_bias) then
                 call compute_and_write_cl(Par%outclfile,Par%ssim,Par%zerofill,Par%endnamecl,alm1,Par%Lmin,bias=red_biasalpha,alms2=alm2)
              else
                 call compute_and_write_cl(Par%outclfile,Par%ssim,Par%zerofill,Par%endnamecl,alm1,Par%Lmin,alms2=alm2)
              endif
           else
              if (Par%compute_biasalpha .and. Par%subtract_bias) then
                 call compute_and_write_cl(Par%outclfile,Par%ssim,Par%zerofill,Par%endnamecl,alm1,Par%Lmin,bias=red_biasalpha)
              else
                 call compute_and_write_cl(Par%outclfile,Par%ssim,Par%zerofill,Par%endnamecl,alm1,Par%Lmin)
              endif

           endif
        endif
        deallocate(alm1,red_one_o_var1)
        if (Par%do_cross) deallocate(alm2,red_one_o_var2)
        if (Par%compute_biasalpha) deallocate(red_biasalpha)
     endif
  else
     if (myid .eq. 0) deallocate(red_one_o_var1)
     if ((myid .eq. 0) .and. Par%do_cross) deallocate(red_one_o_var2)
  endif

  deallocate(clEEfid,clEEobs1,clBBobs1)
  if (Par%do_cross) deallocate(clEEobs2,clBBobs2)  

  call mpi_barrier(mpi_comm_world, mpierr)
  
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) call system_clock(t4, clock_rate, clock_max)
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write(0,*) 'Elapsed real time for total computation = ', real(t4 - t1) / real(clock_rate)
  if ((myid .eq. 0) .and. (Par%feedback .gt. 2)) write(0,*) 'Elapsed real time for computation and initialization = ', real(t4 - t0)/real(clock_rate)
  if ((myid .eq. 0) .and. (Par%feedback .gt. 0)) write(0,*) 'End Program'  
  call mpi_finalize(mpierr)
  
end program EB_estimator
