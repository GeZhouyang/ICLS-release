module mod_write_vtk

 implicit none

 private
 public write_vtk, disp_cg, disp_speed
 contains

! #1 Main
  subroutine write_vtk(istep)

    use mod_param
    use mod_common
    use decomp_2d
    use mpi
    use mod_common_mpi

    implicit none

    integer, intent(in) :: istep
    real, dimension(1:itot, 1:jtot, 1:ktot) :: u,v,w,p,jam
    !# jam stands for jiemian (interface)


    !---- Generate global coordinate variables ----!

    call globe4vtk(u,v,w,p,jam) ! all cell-centered

    !---- Estimate average speed of the dispersed entity ----!

    call disp_speed

    !---- Write vtk files with these variables ----!

    if (myid .eq. 0) call vtk(u,v,w,p,jam,istep) 


    return
  end subroutine write_vtk




! #2 
  subroutine globe4vtk(u_tot,v_tot,w_tot, p_tot,jam_tot)

    use mod_param
    use mod_common
    use decomp_2d
    use mpi
    use mod_common_mpi

    implicit none

    !----- I/O -----!

    real, dimension(1:itot, 1:jtot, 1:ktot), intent(out) :: u_tot,v_tot,w_tot, p_tot,jam_tot

    !---- Buffers ----!

!    real, dimension(1:itot, 1:jtot, 1:ktot) :: u_tot,v_tot,w_tot
    real, dimension(1:imax, 1:jmax, 1:kmax) :: u,v,w
    real, dimension(0:i1, 0:j1, 0:k1) :: p,jam

    !----- Miscellaneous -----!

    integer :: i,j,k
    integer :: npoints1,npoints2, sub_rank
    integer, dimension(2) :: l_coords  ! starts from 0
    integer :: l_x, l_y, i_g, j_g  ! global indices


    u_tot = 0.
    v_tot = 0.
    w_tot = 0.
    p_tot = 0.
    jam_tot = 0.

    npoints1 = imax*jmax*kmax ! total # of (internal, cell-center) points in each processor
    npoints2 = (i1+1)*(j1+1)*(k1+1) ! total # of points in each processor

    !---- Loop through all sub-domains ----!

    do l_x = 0,dims(1)-1
      l_coords(1) = l_x

      do l_y = 0,dims(2)-1   
        l_coords(2) = l_y   

        if (l_coords(1) .eq. 0 .and. l_coords(2) .eq. 0) then
          ! do nothing for rank 0 at the moment (fill it later)
        else

          !---- Determine the rank of current sub-domain ----!

          call MPI_CART_RANK(comm_cart,l_coords,sub_rank,error)

          !---- Pass data (stored in common.f90) from sub_rank to rank 0 ----!
    
          if (myid .eq. sub_rank) then
    
            call MPI_SSEND(uccnew,  npoints1,MPI_REAL8, 0, 61,comm_cart,error)
            call MPI_SSEND(vccnew,  npoints1,MPI_REAL8, 0, 62,comm_cart,error)
            call MPI_SSEND(wccnew,  npoints1,MPI_REAL8, 0, 63,comm_cart,error)
            call MPI_SSEND(pnew,    npoints2,MPI_REAL8, 0, 64,comm_cart,error)
            call MPI_SSEND(jiemian, npoints2,MPI_REAL8, 0, 65,comm_cart,error)
    
          endif
    
          if (myid .eq. 0) then

            !---- Processor 0 recieves data ----!
    
            call MPI_RECV(u,   npoints1,MPI_REAL8, sub_rank, 61,comm_cart,status,error)
            call MPI_RECV(v,   npoints1,MPI_REAL8, sub_rank, 62,comm_cart,status,error)
            call MPI_RECV(w,   npoints1,MPI_REAL8, sub_rank, 63,comm_cart,status,error)
            call MPI_RECV(p,   npoints2,MPI_REAL8, sub_rank, 64,comm_cart,status,error)
            call MPI_RECV(jam, npoints2,MPI_REAL8, sub_rank, 65,comm_cart,status,error)
    
            !---- Processor 0 distributes global values ----!

            do k = 1,kmax
              do j = 1,jmax
                j_g = j  + l_coords(2)*jmax
                do i = 1,imax
                  i_g = i + l_coords(1)*imax
    
                  u_tot(i_g,j_g,k)   = u(i,j,k)
                  v_tot(i_g,j_g,k)   = v(i,j,k)
                  w_tot(i_g,j_g,k)   = w(i,j,k)
                  p_tot(i_g,j_g,k)   = p(i,j,k)
                  jam_tot(i_g,j_g,k) = jam(i,j,k)
    
                enddo
              enddo
            enddo
    
          endif

        endif
  
      enddo
    enddo

    !---- Finally processor 0 updates itself ----!

    if (myid .eq. 0) then
      do k = 1,kmax
        do j = 1,jmax
          do i = 1,imax
  
            u_tot(i,j,k)   = uccnew(i,j,k)
            v_tot(i,j,k)   = vccnew(i,j,k)
            w_tot(i,j,k)   = wccnew(i,j,k)
            p_tot(i,j,k)   = pnew(i,j,k)
            jam_tot(i,j,k) = jiemian(i,j,k)
  
          enddo
        enddo
      enddo
    endif

!    !---- Obtain cell-center velocity ----!
!
!    do k = 1,ktot
!      do j = 1,jtot
!
!        i = 1
!        u_tot_c(i,j,k) = ( u_tot(i,j,k) + u_tot(itot,j,k) )/2.  ! periodicity
!
!        do i = 2,itot
!          u_tot_c(i,j,k) = ( u_tot(i,j,k) + u_tot(i-1,j,k) )/2.
!        enddo
!
!      enddo
!    enddo
!
!    do k = 1,ktot
!
!      j = 1
!      do i = 1,itot
!        v_tot_c(i,j,k) = ( v_tot(i,j,k) + v_tot(i,jtot,k) )/2.  ! periodicity
!      enddo
!
!      do j = 2,jtot
!        do i = 1,itot
!          v_tot_c(i,j,k) = ( v_tot(i,j,k) + v_tot(i,j-1,k) )/2.
!        enddo
!      enddo
!
!    enddo
!
!    k = 1
!    do j = 1,jtot
!      do i = 1,itot
!        w_tot_c(i,j,k) = ( w_tot(i,j,k) + w_tot(i,j,ktot) )/2.  ! periodicity
!      enddo
!    enddo
!
!    do k = 2,ktot
!      do j = 1,jtot
!        do i = 1,itot
!          w_tot_c(i,j,k) = ( w_tot(i,j,k) + w_tot(i,j,k-1) )/2.
!        enddo
!      enddo
!    enddo


    return
  end subroutine globe4vtk


! #3
  subroutine vtk(vx,vy,vz,pr,jam_tot,istep)
  !-- Adapted from JC's version (tpls_io.f90)

    use mod_param
    use mod_common
    use decomp_2d
    use mpi
    use mod_common_mpi

    implicit none

    !----- Inputs (global variables) -----!

    integer, intent(in) :: istep
    real, dimension(1:itot, 1:jtot, 1:ktot), intent(in) :: vx, vy, vz, pr, jam_tot

    !----- Miscellaneous -----!

    integer :: i, j, k, ifich
    integer :: iostatus
    character :: q
    character(len=1), parameter :: newline = char(10)
    character(len=100) :: s_buffer
    real :: coorz
    real, dimension(1:itot, 1:jtot, 1:ktot) :: v_undisturbed

    integer :: output_unit, input_unit
    parameter ( input_unit = 8, output_unit = 9 )

    character(len=7) :: istepchar
    character(len=20) :: filename

    !----- Create file -----!

    ifich = 10
    q = char(34)

    write(istepchar,'(i7.7)') istep
    filename = datadir//'para'//istepchar//'.vtk'

    open( unit = ifich , file = filename , form = 'unformatted' , access = 'stream' , &
          action = 'write' , convert = 'BIG_ENDIAN' , iostat = iostatus )

    write(unit = ifich, iostat = iostatus) '# vtk DataFile Version 3.0' // newline
    write(unit = ifich, iostat = iostatus) 'test file' // newline
    write(unit = ifich, iostat = iostatus) 'BINARY' // newline
    write(unit = ifich, iostat = iostatus) newline
    write(unit = ifich, iostat = iostatus) 'DATASET ' // 'STRUCTURED_POINTS' // newline

    !---- Define the domain ----!

    write(s_buffer, FMT = '(A,3I12)', iostat = iostatus) 'DIMENSIONS', itot, jtot, ktot
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(s_buffer, FMT = '(A,3E14.6E2)', iostat = iostatus) 'ORIGIN', dx/2., dy/2., dz/2.
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(s_buffer, FMT = '(A,3E14.6E2)', iostat = iostatus) 'SPACING', dx, dy, dz
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline

    !---- Write data ----!

    write(unit = ifich, iostat = iostatus) newline
    write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'POINT_DATA', itot*jtot*ktot
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline

    !---- 1. U-velocity ----!

    write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS U double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

    do k = 1,ktot
       do j = 1,jtot
          do i = 1,itot

             write(unit = ifich, iostat = iostatus) vx(i,j,k)

          enddo
       enddo
    enddo

    !---- 2. V-velocity ----!

    write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS V double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

    do k = 1,ktot
       do j = 1,jtot
          do i = 1,itot

             write(unit = ifich, iostat = iostatus) vy(i,j,k)

          enddo
       enddo
    enddo

    !---- 3. W-velocity ----!

    write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS W double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

    do k = 1,ktot
       do j = 1,jtot
          do i = 1,itot

             write(unit = ifich, iostat = iostatus) vz(i,j,k)

          enddo
       enddo
    enddo

    !---- 4. Velocity magnitude ----!

    write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS abs(Vel) double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

    do k = 1,ktot
       do j = 1,jtot
          do i = 1,itot
             write(unit = ifich, iostat = iostatus) sqrt( vx(i,j,k)**2 &
                                                        + vy(i,j,k)**2 &
                                                        + vz(i,j,k)**2 )
          enddo
       enddo
    enddo

!    !---- 5. Relative u-velocity ----!

!    write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS r_U double', 1
!    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
!    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

!    do k = 1,ktot
!       do j = 1,jtot
!          do i = 1,itot

!             write(unit = ifich, iostat = iostatus) vx(i,j,k) -u_disp(1)

!          enddo
!       enddo
!    enddo

    !---- 6. Relative v-velocity ----!

    write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS V_r double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

    do k = 1,ktot

      coorz = (k-0.5)*dz/lz ! normalised channel height [0,1]
      v_undisturbed(:,:,k) = 6.*v_bulk_desired*coorz*(1.-coorz)  ! mean value = 1

       do j = 1,jtot
          do i = 1,itot

             write(unit = ifich, iostat = iostatus) vy(i,j,k) -v_undisturbed(i,j,k)

          enddo
       enddo
    enddo

!    !---- 7. Relative w-velocity ----!

!    write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS r_W double', 1
!    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
!    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

!    do k = 1,ktot
!       do j = 1,jtot
!          do i = 1,itot

!             write(unit = ifich, iostat = iostatus) vz(i,j,k) -w_disp(1)

!          enddo
!       enddo
!    enddo

    !---- 8. Pressure ----!

    write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS Pressure double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

    do k = 1,ktot
       do j = 1,jtot
          do i = 1,itot

             write(unit = ifich, iostat = iostatus) pr(i,j,k)

          enddo
       enddo
    enddo

    !---- 9. Interface ----!

    write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS Interface double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

    do k = 1,ktot
       do j = 1,jtot
          do i = 1,itot

            write(unit = ifich, iostat = iostatus) jam_tot(i,j,k)

          enddo
       enddo
    enddo


    write(unit = ifich, iostat = iostatus) newline
    close(ifich)
   
    
    return
  end subroutine vtk


!-------------------------------------------------------!


  subroutine disp_cg

  !-- Calculate the center of gravity of the dispersed bubbles/droplets.

    use mod_param
    use mod_common
    use mod_common
    use mod_common_mpi

    implicit none

    integer :: i,j,k,l
    real :: coorx, coory, coorz
    real :: lowerbound,leftbound
    real :: counter, counter_all, xd_all, yd_all, zd_all

    !---- Bounds of Each Core ----!

    lowerbound =  coords(2) *ly/(dims(2)*1.)
    leftbound  =  coords(1) *lx/(dims(1)*1.)

    do l = 1,lmax

      counter = 0.
      x_disp(l) = 0.
      y_disp(l) = 0.
      z_disp(l) = 0.

      do k = 1,kmax
        coorz = (k-0.5)/dzi  ! [dz/2, lz-dz/2]
        do j = 1,jmax
          coory = lowerbound + (j-0.5)/dyi  ! [dy/2, ly-dy/2]
          do i = 1,imax
            coorx = leftbound + (i-0.5)/dxi  ! [dx/2, lx-dx/2]

            if (lvset(i,j,k,l) .le. 1.5*dz) then

              counter = counter + 1.

              x_disp(l) = x_disp(l) + coorx
              y_disp(l) = y_disp(l) + coory
              z_disp(l) = z_disp(l) + coorz

            endif

          enddo
        enddo
      enddo

      call mpi_allreduce(counter,counter_all,1,mpi_real8,mpi_sum,comm_cart,error)
      call mpi_allreduce(x_disp(l),xd_all,1,mpi_real8,mpi_sum,comm_cart,error)
      call mpi_allreduce(y_disp(l),yd_all,1,mpi_real8,mpi_sum,comm_cart,error)
      call mpi_allreduce(z_disp(l),zd_all,1,mpi_real8,mpi_sum,comm_cart,error)

      x_disp(l) = xd_all/counter_all
      y_disp(l) = yd_all/counter_all
      z_disp(l) = zd_all/counter_all

    enddo 

    x_disp_avr = sum(x_disp)/lmax
    y_disp_avr = sum(y_disp)/lmax
    z_disp_avr = sum(z_disp)/lmax


    return
  end subroutine disp_cg




  subroutine disp_speed

  !-- Calculate (estimate) the overal speed of the dispersed bubbles/droplets.

    use mod_param
    use mod_common
    use mod_common_mpi

    implicit none

    integer :: i,j,k,l
    real :: u, v, w, ud_all, vd_all, wd_all

    do l = 1,lmax

      u_disp(l) = 0.
      v_disp(l) = 0.
      w_disp(l) = 0.

      do k = 1,kmax
        do j = 1,jmax
          do i = 1,imax
  
            !---- Velocity integral of dispersed phase ----!

            if (lvset(i,j,k,l) .le. 1.5*dz) then

              u = uccnew(i,j,k)
              v = vccnew(i,j,k)
              w = wccnew(i,j,k)

              u_disp(l) = u_disp(l) + u*dx*dy*dz*(1.-hnew_mls(i,j,k,l))
              v_disp(l) = v_disp(l) + v*dx*dy*dz*(1.-hnew_mls(i,j,k,l))
              w_disp(l) = w_disp(l) + w*dx*dy*dz*(1.-hnew_mls(i,j,k,l))

            endif
  
          enddo
        enddo
      enddo

      call mpi_allreduce(u_disp(l),ud_all,1,mpi_real8,mpi_sum,comm_cart,error)
      call mpi_allreduce(v_disp(l),vd_all,1,mpi_real8,mpi_sum,comm_cart,error)
      call mpi_allreduce(w_disp(l),wd_all,1,mpi_real8,mpi_sum,comm_cart,error)

      u_disp(l) = ud_all/(lx*ly*lz*init_volume(l))
      v_disp(l) = vd_all/(lx*ly*lz*init_volume(l))
      w_disp(l) = wd_all/(lx*ly*lz*init_volume(l))

    enddo

    u_disp_avr = sum(u_disp)/lmax
    v_disp_avr = sum(v_disp)/lmax
    w_disp_avr = sum(w_disp)/lmax


    return
  end subroutine disp_speed


end module mod_write_vtk
