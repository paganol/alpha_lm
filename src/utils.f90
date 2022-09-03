module utils
  use healpix_types 
  use settings
  use alm_tools
  use pix_tools
  use fitstools
  use head_fits 
  
  implicit none
  
contains
  
  subroutine make_noise(filename,nl,noiseE,noiseB)
    character(len=FILENAMELEN) :: filename
    real(dp),dimension(0:,1:) , intent(inout) :: nl
    real(dp) :: noiseE,noiseB
    integer :: lmax,myunit,ll,l
    character(len=2048) :: fakestring
    
    if (len(trim(filename)) .gt. 0) then
       
       lmax=size(nl,dim=1)-1
       nl=0.0d0
       
       open(newunit=myunit,file=trim(filename),status='old',form='formatted')
       read(myunit,'(a)') fakestring
       if (scan(fakestring,'#') .eq. 0) rewind(myunit)
       do l=2,lmax
          read(myunit,*) ll,nl(l,1),nl(l,2)
       enddo
       close(myunit)
       
    else
       nl(:,1) = (noiseE*DEG2RAD/60.0)**2
       nl(:,2) = (noiseB*DEG2RAD/60.0)**2
    endif
  end subroutine make_noise
  
  subroutine read_maps_and_compute_alms(filename,ssim,zerofill,endname,iter,almE,almB)
    character(len=FILENAMELEN) :: filename,mapname,endname
    integer(i4b) :: lmax,nside,iter
    integer(i8b) :: npix
    character(len=16) :: simstr
    character(len=1) :: strzerofill
    integer(i4b) :: ct,zerofill,nsims,ssim,isim
    real(dp),allocatable,dimension(:,:) :: maps
    complex(dpc),allocatable,dimension(:,:,:) :: alms   
    complex(dpc),dimension(1:,0:,0:) :: almE,almB
    
    nsims=size(almE,dim=1)
    lmax=size(almE,dim=2)-1    

    allocate(alms(1:3,0:lmax,0:lmax))
    
    if (nsims .eq. 1) then
       mapname=filename
       npix = getsize_fits(trim(mapname),nside=nside)
       allocate(maps(0:npix-1,1:3))
       call input_map(trim(mapname),maps, npix, 3)   
       call map2alm_iterative(nside, lmax, lmax, iter, maps, alms)
       almE(1,:,:)=alms(2,:,:)
       almB(1,:,:)=alms(3,:,:)
    else
       ct=1
       write(strzerofill,fmt='(i1)') zerofill
       write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') ssim
       mapname=trim(filename)//trim(simstr)//trim(endname)
       npix = getsize_fits(trim(mapname),nside=nside)
       allocate(maps(0:npix-1,1:3))
       do isim=ssim,ssim+nsims-1
          write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
          mapname=trim(filename)//trim(simstr)//trim(endname)
          call input_map(trim(mapname),maps, npix, 3)   
          call map2alm_iterative(nside, lmax, lmax, iter, maps, alms)
          almE(ct,:,:)=alms(2,:,:)
          almB(ct,:,:)=alms(3,:,:)
          ct=ct+1
       enddo
    endif
    deallocate(maps,alms)
    
  end subroutine read_maps_and_compute_alms

  subroutine read_map_and_compute_alms(filename,iter,almE,almB,sim)
    character(len=FILENAMELEN) :: filename
    integer(i4b) :: lmax,nside,iter,sim
    integer(i8b) :: npix
    real(dp),allocatable,dimension(:,:) :: maps
    complex(dpc),allocatable,dimension(:,:,:) :: alms
    complex(dpc),dimension(:,0:,0:) :: almE,almB

    lmax=size(almE,dim=2)-1

    allocate(alms(1:3,0:lmax,0:lmax))

    npix = getsize_fits(trim(filename),nside=nside)
    allocate(maps(0:npix-1,1:3))
    call input_map(trim(filename),maps, npix, 3)
    call map2alm_iterative(nside, lmax, lmax, iter, maps, alms)
    almE(sim,:,:)=alms(2,:,:)
    almB(sim,:,:)=alms(3,:,:)
    deallocate(maps,alms)

  end subroutine read_map_and_compute_alms
   
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
    complex(dpc), dimension(1:,0:,0:) :: alms
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

  subroutine compute_and_write_cl(filename,ssim,zerofill,endname,alms,lmin,bias)
    complex(dpc), dimension(1:,0:,0:) :: alms
    real(dp), dimension(:,:),optional :: bias
    real(dp),allocatable, dimension(:,:) :: cl
    character(len=FILENAMELEN) :: filename,clname,endname
    integer(i4b) :: lmax,lmin
    character(len=16) :: simstr
    character(len=1) :: strzerofill
    integer(i4b) :: ct,zerofill,nsims,ssim,isim,il,myunit

    nsims=Size(alms,DIM=1)
    lmax=Size(alms,DIM=2)-1

    allocate(cl(0:lmax,1:1))

    if (nsims .eq. 1) then
       clname=filename
       call alm2cl(lmax,lmax,alms(1:1,:,:),cl)
       if (present(bias)) cl(:,1) = cl(:,1) - bias(1,:) 
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
          call alm2cl(lmax,lmax,alms(ct:ct,:,:),cl)
          if (present(bias)) cl(:,1) = cl(:,1) - bias(ct,:)
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

  subroutine compute_cls_from_alms(almE,almB,clEE,clBB)
    complex(dpc), dimension(1:,0:,0:) :: almE,almB
    real(dp), dimension(1:,0:) :: clEE, clBB
    real(dp),allocatable, dimension(:,:) :: cl
    integer(i4b) :: nsims, lmax, isim

    nsims=Size(almE,dim=1)
    lmax=Size(almE,dim=2)-1    

    allocate(cl(0:lmax,1:1))

    do isim=1,nsims
       call alm2cl(lmax,lmax,almE(isim:isim,:,:),cl)
       clEE(isim,:) = cl(:,1)
       call alm2cl(lmax,lmax,almB(isim:isim,:,:),cl)
       clBB(isim,:) = cl(:,1)
    enddo
 
    deallocate(cl)
  end subroutine compute_cls_from_alms

end module utils
