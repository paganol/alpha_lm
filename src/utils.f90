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
  
  subroutine read_maps_and_compute_alms(filename,ssim,zerofill,endname,almE,almB)
    character(len=FILENAMELEN) :: filename,mapname,endname
    integer(i4b) :: lmax,nside,iter
    integer(i8b) :: npix
    character(len=16) :: simstr
    character(len=1) :: strzerofill
    integer(i4b) :: ct,zerofill,nsims,ssim,isim
    real(dp),allocatable,dimension(:,:) :: maps
    complex(dpc),allocatable,dimension(:,:,:) :: alms   
    complex(dpc),dimension(1:,0:,0:) :: almE,almB
    
    iter = 1
    
    nsims=size(almE,dim=1)
    lmax=size(almE,dim=2)-1    

    allocate(alms(1:3,0:lmax,0:lmax))
    
    if (nsims .eq. 1) then
       mapname=filename
       npix = getsize_fits(trim(mapname),nside=nside)
       allocate(maps(0:npix-1,1:3))
       call input_map(trim(mapname),maps, npix, 3)   
       call map2alm(nside, lmax, lmax, maps, alms)
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
          call map2alm(nside, lmax, lmax, maps, alms)
          almE(ct,:,:)=alms(2,:,:)
          almB(ct,:,:)=alms(3,:,:)
          ct=ct+1
       enddo
    endif
    deallocate(maps,alms)
    
  end subroutine read_maps_and_compute_alms
   
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


  subroutine read_beam(beam,window,beamfile)
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
  
  subroutine write_out_cls(outfile,cls,lmin)
    real(dp), dimension(:) :: cls
    character(len=FILENAMELEN),intent(in) ::  outfile
    integer :: myunit,n,i,j,lmin

    n=Size(cls,DIM=1)
    
    open(newunit=myunit,file=trim(outfile),status='replace',form='formatted') 
    do i=1,n
       write(myunit,'(I4,*(E15.7))') lmin+i-1,cls(i)  
    enddo
    
    close(myunit)
    
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

  subroutine compute_and_write_cl(filename,ssim,zerofill,endname,alms,lmin)
    complex(dpc), dimension(1:,0:,0:) :: alms
    real(dp),allocatable, dimension(:,:) :: cl
    character(len=FILENAMELEN) :: filename,clname,endname
    integer(i4b) :: lmax,lmin
    character(len=16) :: simstr
    character(len=1) :: strzerofill
    integer(i4b) :: ct,zerofill,nsims,ssim,isim

    nsims=Size(alms,DIM=1)
    lmax=Size(alms,DIM=2)-1

    allocate(cl(0:lmax,1:1))

    if (nsims .eq. 1) then
       clname=filename
       call alm2cl(lmax,lmax,alms(1:1,:,:),cl)
       call write_out_cls(clname,cl(lmin:,1),lmin) 
    else
       write(strzerofill,fmt='(i1)') zerofill
       ct=1
       do isim=ssim,ssim+nsims-1
          write (simstr,fmt='(i'//trim(strzerofill)//'.'//trim(strzerofill)//')') isim
          clname=trim(filename)//trim(simstr)//trim(endname)
          call alm2cl(lmax,lmax,alms(ct:ct,:,:),cl)
          call write_out_cls(clname,cl(lmin:,1),lmin)
          ct=ct+1
       enddo
    endif

    deallocate(cl)


  end subroutine compute_and_write_cl

end module utils
