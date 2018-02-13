module mod_loadflds

  use decomp_2d
  use decomp_2d_io
  use mod_param 
  use mod_common_mpi
  use mod_common

  implicit none

  private
  public loadflds, init_from_fld

contains
  
  subroutine loadflds(in,nr)
    !
    ! Note that the dimensions are all (1:imax,1:jmax,1:kmax)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    integer :: in,nr, fh, l
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    character(len=7) :: istepchar
    real, dimension(3) :: fldinfo
    real, dimension(imax,jmax,kmax) :: temp
    
    !-- read data --!
    
    if (in.eq.0) then
       write(istepchar,'(i7.7)') nr
       call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fld'//istepchar, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
       call MPI_FILE_CLOSE(fh,error)
       call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fld'//istepchar, &
            MPI_MODE_RDONLY, MPI_INFO_NULL,fh, error)
       disp = 0_MPI_OFFSET_KIND
     
       call decomp_2d_read_var(fh,disp,3,temp)
       unew(1:imax,1:jmax,1:kmax) = temp
   
       call decomp_2d_read_var(fh,disp,3,temp)
       vnew(1:imax,1:jmax,1:kmax) = temp
   
       call decomp_2d_read_var(fh,disp,3,temp)
       wnew(1:imax,1:jmax,1:kmax) = temp
   
       call decomp_2d_read_var(fh,disp,3,temp)
       pnew(1:imax,1:jmax,1:kmax) = temp
   
     !  call decomp_2d_read_var(fh,disp,3,temp)
     !  ynew(1:imax,1:jmax,1:kmax) = temp
   
       do l = 1,lmax
         call decomp_2d_read_var(fh,disp,3,temp)
         lvset(1:imax,1:jmax,1:kmax,l) = temp
       enddo
   
       call decomp_2d_read_var(fh,disp,3,temp)
       jiemian(1:imax,1:jmax,1:kmax) = temp
     
       call decomp_2d_read_scalar(fh,disp,3,fldinfo)
       time = fldinfo(1)
       nr = int(fldinfo(2))
       dt = fldinfo(3)
       call MPI_FILE_CLOSE(fh,error)
    endif
  
    !-- write data --!
    
    if (in.eq.1) then
       write(istepchar,'(i7.7)') nr
       fldinfo = (/time,1.*nr,dt/)
       call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fld'//istepchar, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
       filesize = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
       disp = 0_MPI_OFFSET_KIND
     
       temp = unew(1:imax,1:jmax,1:kmax)
       call decomp_2d_write_var(fh,disp,3,temp)
   
       temp = vnew(1:imax,1:jmax,1:kmax)
       call decomp_2d_write_var(fh,disp,3,temp)
   
       temp = wnew(1:imax,1:jmax,1:kmax)
       call decomp_2d_write_var(fh,disp,3,temp)
   
       temp = pnew(1:imax,1:jmax,1:kmax)
       call decomp_2d_write_var(fh,disp,3,temp)
   
     !  temp = ynew(1:imax,1:jmax,1:kmax)
     !  call decomp_2d_write_var(fh,disp,3,temp)
   
       do l = 1,lmax
         temp = lvset(1:imax,1:jmax,1:kmax,l)
         call decomp_2d_write_var(fh,disp,3,temp)
       enddo
   
       temp = jiemian(1:imax,1:jmax,1:kmax)  
       call decomp_2d_write_var(fh,disp,3,temp)
   
       call decomp_2d_write_scalar(fh,disp,3,fldinfo)
       call MPI_FILE_CLOSE(fh,error)
    endif

  
    return
  end subroutine loadflds


  !------------------------------------------------


  subroutine init_from_fld
    !
    ! Initialize unew,vnew,wnew,pnew from init_fld (some previous state).
    ! Note that the dimensions are all (1:imax,1:jmax,1:kmax).
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    real, dimension(imax,jmax,kmax) :: temp
    
    !-- read fld_init --!
    
    call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fld_init', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    call MPI_FILE_CLOSE(fh,error)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fld_init', &
         MPI_MODE_RDONLY, MPI_INFO_NULL,fh, error)
    disp = 0_MPI_OFFSET_KIND
    
    call decomp_2d_read_var(fh,disp,3,temp)
    unew(1:imax,1:jmax,1:kmax) = temp
   
    call decomp_2d_read_var(fh,disp,3,temp)
    vnew(1:imax,1:jmax,1:kmax) = temp
   
    call decomp_2d_read_var(fh,disp,3,temp)
    wnew(1:imax,1:jmax,1:kmax) = temp
   
    call decomp_2d_read_var(fh,disp,3,temp)
    pnew(1:imax,1:jmax,1:kmax) = temp
    
    call MPI_FILE_CLOSE(fh,error)

  
    return
  end subroutine init_from_fld


end module mod_loadflds
