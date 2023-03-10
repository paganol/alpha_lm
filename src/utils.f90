module utils
  use healpix_types 
  use settings
  use alm_tools
  use pix_tools
  use fitstools
  use head_fits 
  
  implicit none
  
contains
  
  subroutine make_noise(filename,nl,noiseE,noiseB,feedback)
    character(len=FILENAMELEN) :: filename
    real(dp),dimension(0:,1:) , intent(inout) :: nl
    real(dp) :: noiseE,noiseB
    integer :: lmax,myunit,ll,l,lstart,feedback
    character(len=2048) :: fakestring
    
    if (len(trim(filename)) .gt. 0) then
       if (feedback .gt. 1) write(0,*) 'Reading noise spectra'       
       lmax=size(nl,dim=1)-1
       nl=0.0d0
       
       open(newunit=myunit,file=trim(filename),status='old',form='formatted')
       read(myunit,'(a)') fakestring
       if (scan(fakestring,'#') .eq. 0) rewind(myunit)

       read(myunit,*) lstart,noiseE,noiseB
       nl(lstart,1) = noiseE
       nl(lstart,2) = noiseB
       do l=lstart+1,lmax
          read(myunit,*) ll,nl(l,1),nl(l,2)
       enddo
       close(myunit)
    else
       if (feedback .gt. 1) write(0,*) 'Building noise spectra'
       nl(:,1) = (noiseE*DEG2RAD/60.0)**2
       nl(:,2) = (noiseB*DEG2RAD/60.0)**2
    endif

  end subroutine make_noise

  subroutine read_mask_and_compute_fsky(filename,mask,fsky)
    character(len=FILENAMELEN) :: filename
    integer(i4b) :: nmaps
    integer(i8b) :: npix
    real(sp),intent(out),dimension(0:,1:) :: mask
    real(dp),allocatable,dimension(:,:) :: maps
    real(dp),intent(out) :: fsky

    npix = getsize_fits(trim(filename),nmaps=nmaps)
    allocate(maps(0:npix-1,1:nmaps))
    call input_map(trim(filename),maps, npix, nmaps)
    if (nmaps .eq. 1) then
       mask(:,1) = maps(:,1)
       mask(:,2) = maps(:,1)
       mask(:,3) = maps(:,1)
    endif
    if (nmaps .eq. 2) then
       mask(:,1) = maps(:,1)
       mask(:,2) = maps(:,2)
       mask(:,3) = maps(:,2)
    endif
    if (nmaps .eq. 3) mask = maps

    fsky = sum(mask(:,2))/real(npix)

    deallocate(maps)

  end subroutine read_mask_and_compute_fsky

  subroutine read_map_and_compute_alms(filename,iter,almE,almB,lmax,sim,mask)
    character(len=FILENAMELEN) :: filename
    integer(i4b) :: lmax,nside,iter,sim
    integer(i8b) :: npix
    integer :: imax,ind,lm(2)
    real(sp),allocatable,dimension(:,:) :: maps
    real(sp),optional,dimension(0:,1:) :: mask
    complex(spc),allocatable,dimension(:,:,:) :: alms
    complex(spc),dimension(1:,0:) :: almE,almB

    imax=size(almE,dim=2)-1
    allocate(alms(1:3,0:lmax,0:lmax))

    npix = getsize_fits(trim(filename),nside=nside)
    allocate(maps(0:npix-1,1:3))
    call input_map(trim(filename),maps, npix, 3)
 
    if (present(mask)) then
       call map2alm_iterative(nside, lmax, lmax, iter, maps, alms,mask=mask)
    else
       call map2alm_iterative(nside, lmax, lmax, iter, maps, alms)
    endif

    do ind=0,imax 
       lm = index2lm(lmax,ind)
       almE(sim,ind)=alms(2,lm(1),lm(2))
       almB(sim,ind)=alms(3,lm(1),lm(2))
    enddo

    deallocate(maps,alms)

  end subroutine read_map_and_compute_alms

  subroutine read_precomputed_alms(filename,almE,almB,sim)
    character(len=FILENAMELEN) :: filename
    integer(i4b) :: sim
    integer(i4b) :: ind,imax,nalm
    real(dp),allocatable,dimension(:,:,:) :: alms
    complex(spc),dimension(1:,0:) :: almE,almB
    character(len=80), dimension(80,2) :: header

    nalm = getsize_fits(trim(filename))
    allocate(alms(1:nalm,1:4,1:2))
    
    call fits2alms(trim(filename),nalm,alms,3,header,80,2)

    imax=size(almE,dim=2)-1

    do ind=0,imax
       almE(sim,ind) = cmplx(alms(ind+1,3,1),alms(ind+1,4,1),spc)
       almB(sim,ind) = cmplx(alms(ind+1,3,2),alms(ind+1,4,2),spc) 
    enddo

    deallocate(alms)

  end subroutine read_precomputed_alms
   
  subroutine read_cl(filename,cl)
    real(dp),dimension(0:,1:) :: cl
    character(len=FILENAMELEN),intent(in) ::  filename
    integer :: lmax,myunit,ll,l
    character(len=2048) :: fakestring
    
    lmax=size(cl,dim=1)-1
    cl=0.0d0
 
    open(newunit=myunit,file=trim(filename),status='old',form='formatted')
    read(myunit,'(a)') fakestring
    if (scan(fakestring,'#') .eq. 0) rewind(myunit)
    do l=2,lmax
       read(myunit,*) ll,cl(l,myTT),cl(l,myEE),cl(l,myBB),cl(l,myTE)
       cl(l,:) = cl(l,:)/ll/(ll+1)*TWOPI
    enddo
    close(myunit)

  end subroutine read_cl


  subroutine read_beam(beamfile,beam,window)
    real(dp),dimension(0:,1:) , intent(out) :: window, beam
    integer(i4b) :: lmax
    character(len=80), DIMENSION(74) :: headercl
    character(len=FILENAMELEN) , intent(in) :: beamfile 
   
    lmax=size(beam,dim=1)-1

    call fits2cl(trim(beamfile),beam,lmax,3,headercl)
    
    window(:,myTT)=beam(:,1)*beam(:,1)
    window(:,myTE)=beam(:,1)*beam(:,2)
    window(:,myTB)=beam(:,1)*beam(:,3)
    window(:,myEE)=beam(:,2)*beam(:,2)
    window(:,myEB)=beam(:,2)*beam(:,3)
    window(:,myBB)=beam(:,3)*beam(:,3)
    
  end subroutine read_beam
  
  subroutine write_out_cls(filename,ssim,zerofill,endname,cls,lmin)
    real(dp), dimension(:,:) :: cls
    character(len=FILENAMELEN),intent(in) ::  filename,endname
    character(len=FILENAMELEN) :: clname
    integer :: myunit,lmin,nell,ct
    character(len=16) :: simstr
    character(len=1) :: strzerofill
    integer(i4b) :: zerofill,nsims,ssim,isim,il

    nsims=Size(cls,dim=1)
    nell=Size(cls,dim=2)

    if (nsims .eq. 1) then
       clname=filename
       open(newunit=myunit,file=trim(clname),status='replace',form='formatted')
       do il=1,nell
          write(myunit,'(I4,*(E15.7))') lmin+il-1,cls(1,il)
       enddo
       close(myunit)
    else
       write(strzerofill,fmt='(i1)') zerofill
       ct=1
       do isim=ssim,ssim+nsims-1
          write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
          clname=trim(filename)//trim(simstr)//trim(endname)
          open(newunit=myunit,file=trim(clname),status='replace',form='formatted')
          do il=1,nell
             write(myunit,'(I4,*(E15.7))') lmin+il-1,cls(ct,il)
          enddo
          close(myunit)
          ct=ct+1
       enddo
    endif
    
  end subroutine write_out_cls
  
  subroutine write_out_alms(filename,ssim,zerofill,endname,alms)
    complex(spc), dimension(1:,0:,0:) :: alms
    character(len=FILENAMELEN) :: filename,almname,endname
    integer(i4b) :: lmax
    character(len=16) :: simstr
    character(len=1) :: strzerofill
    integer(i4b) :: ct,zerofill,nsims,ssim,isim
    character(len=80), dimension(1:60) :: header 
    
    
    nsims=Size(alms,DIM=1)
    lmax=Size(alms,DIM=2)-1
    
    call write_minimal_header(header, 'ALM', nlmax=lmax,nmmax=lmax,polar=.False.) 
    
    if (nsims .eq. 1) then
       almname=filename
       call dump_alms(trim(almname),alms(1,:,:),lmax,header,60,0)
    else
       write(strzerofill,fmt='(i1)') zerofill
       ct=1
       do isim=ssim,ssim+nsims-1
          write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
          almname=trim(filename)//trim(simstr)//trim(endname)
          call dump_alms(trim(almname),alms(ct,:,:),lmax,header,60,0)          
          ct=ct+1
       enddo
    endif
    
  end subroutine write_out_alms

  subroutine compute_and_write_cl(filename,ssim,zerofill,endname,alms1,lmin,bias,alms2)
    complex(spc), dimension(1:,0:,0:) :: alms1
    complex(spc), optional, dimension(1:,0:,0:) :: alms2
    real(dp), optional, dimension(:,:) :: bias
    real(sp),allocatable, dimension(:,:) :: cl
    character(len=FILENAMELEN) :: filename,clname,endname
    integer(i4b) :: lmax,lmin
    character(len=16) :: simstr
    character(len=1) :: strzerofill
    integer(i4b) :: ct,zerofill,nsims,ssim,isim,il,myunit

    nsims=Size(alms1,DIM=1)
    lmax=Size(alms1,DIM=2)-1

    allocate(cl(0:lmax,1:1))

    if (nsims .eq. 1) then
       clname=filename
       if (present(alms2)) then
          call alm2cl(lmax,lmax,alms1(1:1,:,:),alms2(1:1,:,:),cl)
       else
          call alm2cl(lmax,lmax,alms1(1:1,:,:),cl)
       endif
       if (present(bias)) cl(lmin:lmax,1) = cl(lmin:lmax,1) - bias(1,:) 
       open(newunit=myunit,file=trim(clname),status='replace',form='formatted')
       do il=lmin,lmax
          write(myunit,'(I4,*(E15.7))') il,cl(il,1)
       enddo
       close(myunit)
    else
       write(strzerofill,fmt='(i1)') zerofill
       ct=1
       do isim=ssim,ssim+nsims-1
          write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
          clname=trim(filename)//trim(simstr)//trim(endname)
          if (present(alms2)) then
             call alm2cl(lmax,lmax,alms1(ct:ct,:,:),alms2(ct:ct,:,:),cl)
          else
             call alm2cl(lmax,lmax,alms1(ct:ct,:,:),cl)
          endif
          if (present(bias)) cl(lmin:lmax,1) = cl(lmin:lmax,1) - bias(ct,:)
          open(newunit=myunit,file=trim(clname),status='replace',form='formatted')
          do il=lmin,lmax
             write(myunit,'(I4,*(E15.7))') il,cl(il,1)
          enddo
          close(myunit)
          ct=ct+1
       enddo
    endif

    deallocate(cl)

  end subroutine compute_and_write_cl

  subroutine compute_cls_from_alms(almE1,almB1,lmax,clEE,clBB,almE2,almB2)
    complex(spc), dimension(1:,0:) :: almE1,almB1
    complex(spc), optional, dimension(1:,0:) :: almE2,almB2
    complex(spc),allocatable, dimension(:,:,:) :: alm1,alm2
    real(dp), dimension(1:,0:) :: clEE, clBB
    real(sp),allocatable, dimension(:,:) :: cl
    integer(i4b) :: nsims, lmax, isim
    integer :: ind, imax, lm(2)

    nsims=Size(almE1,dim=1)
    imax=Size(almE1,dim=2)-1    

    allocate(cl(0:lmax,1:1))

    if (present(almE2) .and. present(almB2)) then
       allocate(alm1(1:1,0:lmax,0:lmax))
       allocate(alm2(1:1,0:lmax,0:lmax))    
       do isim=1,nsims
          do ind=0,imax
             lm=index2lm(lmax,ind)
             alm1(1,lm(1),lm(2)) = almE1(isim,ind) 
             alm2(1,lm(1),lm(2)) = almE2(isim,ind)
          enddo
          call alm2cl(lmax,lmax,alm1,alm2,cl)
          clEE(isim,:) = cl(:,1)
          do ind=0,imax
             lm=index2lm(lmax,ind)
             alm1(1,lm(1),lm(2)) = almB1(isim,ind)           
             alm2(1,lm(1),lm(2)) = almB2(isim,ind)
          enddo
          call alm2cl(lmax,lmax,alm1,alm2,cl)
          clBB(isim,:) = cl(:,1)
       enddo
       deallocate(alm1,alm2)
    else
       allocate(alm1(1:1,0:lmax,0:lmax))
       do isim=1,nsims
          do ind=0,imax
             lm=index2lm(lmax,ind)
             alm1(1,lm(1),lm(2)) = almE1(isim,ind)           
          enddo
          call alm2cl(lmax,lmax,alm1,cl)
          clEE(isim,:) = cl(:,1)
          do ind=0,imax
             lm=index2lm(lmax,ind)
             alm1(1,lm(1),lm(2)) = almB1(isim,ind)
          enddo
          call alm2cl(lmax,lmax,alm1,cl)
          clBB(isim,:) = cl(:,1)
       enddo
       deallocate(alm1)
    endif 
    deallocate(cl)
  end subroutine compute_cls_from_alms

  function lm2index(lmax,l,m) result(ind)
    integer, intent(in) :: lmax,l,m
    integer :: ind

    ind = floor(m*(2.d0*lmax+1-m)/2.d0)+l
  end function lm2index

  function index2lm(lmax,ind) result(lm)
    integer, intent(in) :: lmax,ind
    integer, dimension(2) :: lm
    real(dp) :: twolmaxp1
   
    twolmaxp1 = 2.d0 * lmax +1.d0
    lm(2) = ceiling((twolmaxp1-sqrt(twolmaxp1*twolmaxp1-8*(ind-lmax)))/2.d0)
    lm(1) = ind - floor(lm(2)*(twolmaxp1-lm(2))/2.d0)
    
  end function index2lm

end module utils
