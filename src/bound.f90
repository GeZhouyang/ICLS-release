module mod_bound

  use mod_param

  implicit none

  private
  public bounduvw,boundp,boundc, boundfu,boundfv,boundfw, boundsls,boundj,boundmls,boundvof, &
         updthalos, updthalos_1d, updthalos_ls, updthalos_mls

contains

  !-------------------------
  !
  subroutine bounduvw(u,v,w)

    use mpi
    use mod_common
    use mod_common_mpi

    integer :: i,j,k
    real :: frontbound,backbound, coory,coorz
    real :: v_lower,v_upper
    real, dimension(0:,0:,0:) :: u,v,w

    select case(BC_in_z)
    case('Dirichl')

       if (ns_np) then
          v_lower = 0.0
          v_upper = 0.0    
       else
          v_lower =-1.5  ! lower plate velocity
          v_upper = 1.5  ! upper plate velocity
       endif
       
       do j=0,j1
          do i=0,i1
             u(i,j,0)    = -u(i,j,1)                ! no-slip
             u(i,j,k1)   = -u(i,j,kmax)             ! no-slip
             v(i,j,0)    = 2.*v_lower - v(i,j,1)    ! prescribed
             v(i,j,k1)   = 2.*v_upper - v(i,j,kmax) ! prescribed
             w(i,j,0)    = 0.                       ! no-penetration
             w(i,j,kmax) = 0.                       ! no-penetration
             w(i,j,k1)   = w(i,j,kmax-1)            ! not used
          enddo
       enddo

    case('shr-drv')
       
       v_lower = 0.  ! lower plate velocity
       v_upper = 1.  ! upper plate *shear rate* for v

       do j=0,j1
          do i=0,i1
             u(i,j,0)    = -u(i,j,1)                ! no-slip
             u(i,j,k1)   = -u(i,j,kmax)             ! no-slip    
             v(i,j,0)    = 2.*v_lower - v(i,j,1)    ! prescribed
             v(i,j,k1)   = v(i,j,kmax) + dz*v_upper ! prescribed normal derivative
             w(i,j,0)    = 0.                       ! no-penetration
             w(i,j,kmax) = 0.                       ! no-penetration
             w(i,j,k1)   = w(i,j,kmax-1)            ! not used
          enddo
       enddo

    end select

    ! pass x halo data
    
    call updthalos(u,1)
    call updthalos(v,1)
    call updthalos(w,1)

    ! pass y halo data

    call updthalos(u,2)
    call updthalos(v,2)
    call updthalos(w,2)


    select case(BC_in_y)
    case('Walls')

       frontbound  = (coords(2)  )*ly/(1.*dims(2)) ! front boundary of process myid
       backbound   = (coords(2)+1)*ly/(1.*dims(2)) ! back  boundary of process myid

       if (frontbound .eq. 0.) then
          do k=1,kmax
             coorz = (k - 0.5)*dz/lz
             do i=1,imax

                u(i,0,k) = -u(i,1,k)
                w(i,0,k) = -w(i,1,k)                   
                v(i,0,k) = 0.

             enddo
          enddo
       endif

       if (backbound .eq. ly) then
          do k=1,kmax
             do i=1,imax

                u(i,j1,k) = -u(i,jmax,k)
                v(i,jmax,k) = 0.
                v(i,j1,k) = v(i,jmax,k)
                w(i,j1,k) = -w(i,jmax,k)

             enddo
          enddo
       endif

    case('InOut')

       frontbound  = (coords(2)  )*ly/(1.*dims(2)) ! front boundary of process myid
       backbound   = (coords(2)+1)*ly/(1.*dims(2)) ! back  boundary of process myid

       if (frontbound .eq. 0.) then
          do k=1,kmax
             coorz = (k - 0.5)*dz/lz
             do i=1,imax

                u(i,0,k) = -u(i,1,k)
                w(i,0,k) = -w(i,1,k)   
                v(i,0,k) = v_inflow(i,k)                          ! specified inflow profile
                !v(i,0,k) = 6.*v_bulk_desired*coorz*(1.-coorz)     ! Poiseuille flow

             enddo
          enddo
       endif

       if (backbound .eq. ly) then
          do k=1,kmax
             do i=1,imax

                u(i,j1,k) = u(i,jmax,k)
                v(i,j1,k) = 2.*v(i,jmax,k)-v(i,jmax-1,k)  ! compatible div-free
                w(i,j1,k) = w(i,jmax,k)

             enddo
          enddo
       endif

    case('OutOut')
    
       frontbound  = (coords(2)  )*ly/(1.*dims(2)) ! front boundary of process myid
       backbound   = (coords(2)+1)*ly/(1.*dims(2)) ! back  boundary of process myid
    
       if (frontbound .eq. 0.) then
          do k=1,kmax
             do i=1,imax
    
                u(i,0,k) = u(i,1,k)
                w(i,0,k) = w(i,1,k)
                v(i,0,k) = v(i,1,k) + dy/dz*(w(i,1,k)-w(i,1,k-1))  ! div-free zero-Neumann
                
             enddo
          enddo
       endif
    
       if (backbound .eq. ly) then
          do k=1,kmax
             do i=1,imax
    
                u(i,j1,k) = u(i,jmax,k)
                v(i,j1,k) = 2.*v(i,jmax,k)-v(i,jmax-1,k)  ! compatible div-free
                w(i,j1,k) = w(i,jmax,k)
    
             enddo
          enddo
       endif

    end select


    return
  end subroutine bounduvw
  

  !-------------------
  !
  subroutine boundp(p)

    use mpi
    use mod_common_mpi

    integer :: i,j,k
    real, dimension(0:,0:,0:) :: p
    real :: frontbound,backbound

    !$omp parallel default(none) shared(p) private(i,j)
    !$omp do
    do j=0,j1
       do i=0,i1
          p(i,j,0) = p(i,j,1)     ! Newmann (consistent with no/free-slip)
          p(i,j,k1) = p(i,j,kmax) ! Newmann (consistent with no/free-slip)
       enddo
    enddo
    !$omp end parallel

    ! pass x and y halo data

    call updthalos(p,1)
    call updthalos(p,2)


    select case(BC_in_y)
    case('Walls')

       frontbound  = (coords(2)  )*ly/(1.*dims(2)) ! front boundary of process myid
       backbound   = (coords(2)+1)*ly/(1.*dims(2)) ! back  boundary of process myid

       if (frontbound .eq. 0.) then
          do k=1,kmax
             do i=1,imax

                p(i,0,k) = p(i,1,k)  ! zero Neumann

             enddo
          enddo
       endif

       if (backbound .eq. ly) then
          do k=1,kmax
             do i=1,imax

                p(i,j1,k) = p(i,jmax,k)  ! zero Neumann

             enddo
          enddo
       endif

    case('InOut')

       frontbound  = (coords(2)  )*ly/(1.*dims(2)) ! front boundary of process myid
       backbound   = (coords(2)+1)*ly/(1.*dims(2)) ! back  boundary of process myid

       if (frontbound .eq. 0.) then
          do k=1,kmax
             do i=1,imax

                p(i,0,k) = p(i,1,k)  ! zero Neumann

             enddo
          enddo
       endif

       if (backbound .eq. ly) then
          do k=1,kmax
             do i=1,imax

                p(i,j1,k) = -p(i,jmax,k)  ! zero Dirichlet

             enddo
          enddo
       endif

    case('OutOut')

       frontbound  = (coords(2)  )*ly/(1.*dims(2)) ! front boundary of process myid
       backbound   = (coords(2)+1)*ly/(1.*dims(2)) ! back  boundary of process myid

       if (frontbound .eq. 0.) then
          do k=1,kmax
             do i=1,imax

                p(i,0,k) = -p(i,1,k)  ! zero Dirichlet

             enddo
          enddo
       endif

       if (backbound .eq. ly) then
          do k=1,kmax
             do i=1,imax

                p(i,j1,k) = -p(i,jmax,k)  ! zero Dirichlet

             enddo
          enddo
       endif

    end select


    return
  end subroutine boundp


  !----------------------
  !
  subroutine boundc(curv)

    real, dimension(0:,0:,0:) :: curv

    ! do nothing in the wall normal direction

    ! pass x and y halo data

    call updthalos(curv,1)
    call updthalos(curv,2)


    return
  end subroutine boundc


  !------------------------
  !
  subroutine boundfu(flux)

    use mpi
    use mod_common_mpi

    integer :: i,j,k
    real, dimension(-2:,-2:,-2:) :: flux

    write(*,*) '---------------------'
    write(*,*) 'Bummer! WENO5 is not fully implemented yet (ghost points).'
    write(*,*) 'It may produce erroneous results when used with Inflow/Outflow B.C. or at high Re.'
    write(*,*) 'This should be fixed later.'
    write(*,*) 'Program aborted.'
    stop
    
    !if(isfreeslip) then
    if (.not. ns_np) then
       do j=0,j1
          do i=0,i1
             flux(i,j,0)    = flux(i,j,1)   
             flux(i,j,k1)   = flux(i,j,kmax)
          enddo
       enddo
    else
       do j=0,j1
          do i=0,i1
             flux(i,j,0)    = -flux(i,j,1)   
             flux(i,j,k1)   = -flux(i,j,kmax)
          enddo
       enddo
    endif

    !-- Periodic BC in x and y --!

    call updthalos_ls(flux,1)
    call updthalos_ls(flux,2)


    return
  end subroutine boundfu

  subroutine boundfv(flux)

    use mpi
    use mod_common_mpi

    integer :: i,j,k
    real, dimension(-2:,-2:,-2:) :: flux


    return
  end subroutine boundfv


  !-----------------------
  !
  subroutine boundfw(flux)

    use mpi
    use mod_common_mpi

    integer :: i,j,k
    real, dimension(-2:,-2:,-2:) :: flux

    !-- Periodic BC in x and y --!

    call updthalos_ls(flux,1)
    call updthalos_ls(flux,2)


    return
  end subroutine boundfw


  !-------------------------
  !
  subroutine boundsls(y_sls)

    use mpi
    use mod_common_mpi
    use mod_common

    integer :: i,j,k
    real :: frontbound,backbound, phi,dphidz,d2phidz2
    real, dimension(10) :: phi_v
    real, dimension(-2:,-2:,-2:) :: y_sls
    character(len=4), parameter :: extrap = 'line'

    !-- extrapolation of level set at walls --!

    select case(extrap)

    case('line')

       do j=-2,j1+2
          do i=-2,i1+2

             y_sls(i,j,0)  = 2.*y_sls(i,j,1)- y_sls(i,j,2)     
             y_sls(i,j,k1) = 2.*y_sls(i,j,kmax)- y_sls(i,j,kmax-1)  

          enddo
       enddo

    case('quad')

       do j=-2,j1+2
          do i=-2,i1+2

             ! lower side

             phi_v = matmul(pseudo_inv_Ab_lower, reshape(y_sls(i-1:i+1,j-1:j+1,1:3), (/27/)))

             phi      = phi_v(1)
             dphidz   = phi_v(4)
             d2phidz2 = phi_v(7)

             y_sls(i,j,0)  = phi - dphidz*dz + 0.5*d2phidz2*dz**2

             ! upper side

             phi_v = matmul(pseudo_inv_Ab_upper, reshape(y_sls(i-1:i+1,j-1:j+1,kmax-2:kmax), (/27/)))

             phi      = phi_v(1)
             dphidz   = phi_v(4)
             d2phidz2 = phi_v(7)

             y_sls(i,j,k1)  = phi + dphidz*dz + 0.5*d2phidz2*dz**2

          enddo
       enddo

    end select

    !-- Periodic BC in x and y --!

    call updthalos_ls(y_sls,1)
    call updthalos_ls(y_sls,2)

    !-- If not periodic in y  --!

    if (open_end_in_y) then
       
       frontbound  = (coords(2)  )*ly/(1.*dims(2)) ! front boundary of process myid
       backbound   = (coords(2)+1)*ly/(1.*dims(2)) ! back  boundary of process myid

       if (frontbound .eq. 0.) then
          do k=1,kmax
             do i=1,imax

                y_sls(i,0,k) = 2.*y_sls(i,1,k)- y_sls(i,2,k)  ! linear extrapolation
                y_sls(i,-1,k) = 2.*y_sls(i,0,k)- y_sls(i,1,k)
                y_sls(i,-2,k) = 2.*y_sls(i,-1,k)- y_sls(i,0,k)

             enddo
          enddo
       endif

       if (backbound .eq. ly) then
          do k=1,kmax
             do i=1,imax

                y_sls(i,j1,k) = 2.*y_sls(i,jmax,k)- y_sls(i,jmax-1,k)  ! linear extrapolation
                y_sls(i,j1+1,k) = 2.*y_sls(i,j1,k)- y_sls(i,jmax,k)
                y_sls(i,j1+2,k) = 2.*y_sls(i,j1+1,k)- y_sls(i,j1,k)

             enddo
          enddo
       endif

    endif


    return
  end subroutine boundsls


  !----
  subroutine boundj(jump)

    use mpi
    use mod_common_mpi
    use mod_newtypes

    type(real_3by3_matrix), dimension(-2:,-2:,-2:) :: jump
    real ,dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: temp
    integer :: i,j

    do i = 1,3
       do j = 1,3

          temp(:,:,:) = jump(:,:,:)%grad(i,j)

          !-- Periodic BC in x and y --!

          call updthalos_ls(temp,1)
          call updthalos_ls(temp,2)

          jump(:,:,:)%grad(i,j) = temp(:,:,:)

       enddo
    enddo

    return
  end subroutine boundj


  !-------------------------
  !
  subroutine boundmls(y_mls)

    use mpi
    use mod_common_mpi
    use mod_common

    integer :: i,j,k,l
    real :: frontbound,backbound, phi,dphidz,d2phidz2
    real, dimension(10) :: phi_v
    real, dimension(-2:,-2:,-2:,1:) :: y_mls
    real ,dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: y_temp
    character(len=4), parameter :: extrap = 'line'

    !-- extrapolation of level set at walls --!

    do l=1,lmax

       select case(extrap)
          
       case('line')
          
          do j=-2,j1+2
             do i=-2,i1+2

                y_mls(i,j, 0,l)  = 2.*y_mls(i,j, 1,l)- y_mls(i,j,2,l)
                y_mls(i,j,-1,l)  = 2.*y_mls(i,j, 0,l)- y_mls(i,j,1,l)
                y_mls(i,j,-2,l)  = 2.*y_mls(i,j,-1,l)- y_mls(i,j,0,l)
                
                y_mls(i,j,k1,  l) = 2.*y_mls(i,j,kmax,l)- y_mls(i,j,kmax-1,l)
                y_mls(i,j,k1+1,l) = 2.*y_mls(i,j,k1  ,l)- y_mls(i,j,kmax  ,l)
                y_mls(i,j,k1+2,l) = 2.*y_mls(i,j,k1+1,l)- y_mls(i,j,k1    ,l)

             enddo
          enddo

       case('quad')
          
          y_temp(:,:,:) = y_mls(:,:,:,l)

          do j=-2,j1+2
             do i=-2,i1+2

                ! lower side
                
                phi_v = matmul(pseudo_inv_Ab_lower, reshape(y_temp(i-1:i+1,j-1:j+1,1:3), (/27/)))

                phi      = phi_v(1)
                dphidz   = phi_v(4)
                d2phidz2 = phi_v(7)

                y_mls(i,j,0,l)  = phi - dphidz*dz + 0.5*d2phidz2*dz**2

                ! upper side

                phi_v = matmul(pseudo_inv_Ab_upper, reshape(y_temp(i-1:i+1,j-1:j+1,kmax-2:kmax), (/27/)))

                phi      = phi_v(1)
                dphidz   = phi_v(4)
                d2phidz2 = phi_v(7)

                y_mls(i,j,k1,l)  = phi + dphidz*dz + 0.5*d2phidz2*dz**2
                
             enddo
          enddo

       end select
       
    enddo

    !-- Periodic BC in x and y --!

    do l=1,lmax
       y_temp(:,:,:) = y_mls(:,:,:,l)
       call updthalos_ls(y_temp,1)
       call updthalos_ls(y_temp,2)
       y_mls(:,:,:,l) = y_temp(:,:,:)
    enddo

    !-- If not periodic in y  --!

    if (open_end_in_y) then
       
       frontbound  = (coords(2)  )*ly/(1.*dims(2)) ! front boundary of process myid
       backbound   = (coords(2)+1)*ly/(1.*dims(2)) ! back  boundary of process myid

       if (frontbound .eq. 0.) then
          do l=1,lmax
             do k=1,kmax
                do i=1,imax

                   y_mls(i,0,k,l) = 2.*y_mls(i,1,k,l)- y_mls(i,2,k,l)  ! linear extrapolation

                enddo
             enddo
          enddo
       endif
       
       if (backbound .eq. ly) then
          do l=1,lmax
             do k=1,kmax
                do i=1,imax

                   y_mls(i,j1,k,l) = 2.*y_mls(i,jmax,k,l)- y_mls(i,jmax-1,k,l)  ! linear extrapolation

                enddo
             enddo
          enddo
       endif

    endif
    

    return
  end subroutine boundmls

  
  !--------------------------
  !
  subroutine boundvof(y_mvof)

    use mpi
    use mod_common_mpi

    integer :: l
    real, dimension(0:,0:,0:,1:) :: y_mvof
    real ,dimension(0:i1,0:j1,0:k1) :: vof_temp


    !-- Periodic BC in x and y --!

    do l=1,lmax
       vof_temp(:,:,:) = y_mvof(:,:,:,l)
       call updthalos(vof_temp,1)
       call updthalos(vof_temp,2)
       y_mvof(:,:,:,l) = vof_temp(:,:,:)
    enddo


    return
  end subroutine boundvof



  
  !~~~~~~~~~~~~~~~~! 
  !                !
  !  Halo passing  !
  !                !
  !~~~~~~~~~~~~~~~~!
  
  subroutine updthalos(var,dir)

    use mpi
    use mod_common_mpi

    real, dimension(0:,0:,0:), intent(inout) :: var
    integer, intent(in) :: dir
    integer :: requests(4), statuses(MPI_STATUS_SIZE,4)

    !  This subroutine updates the halos that store info
    !  from the neighboring computational sub-domain

    select case(dir)
    case(1) ! x direction
       call MPI_SENDRECV(var(1,0,0),1,xhalo,left,0,   &
                         var(i1,0,0),1,xhalo,right,0, &
                         comm_cart,status,error)
       call MPI_SENDRECV(var(imax,0,0),1,xhalo,right,0, &
                         var(0,0,0),1,xhalo,left,0,     &
                         comm_cart,status,error)
       !call MPI_IRECV(var(0,0,0),1,xhalo,left,1, &
       !               comm_cart,requests(2),error)
       !call MPI_IRECV(var(i1,0,0),1,xhalo,right,0, &
       !               comm_cart,requests(1),error)
       !call MPI_ISSEND(var(imax,0,0),1,xhalo,right,1, &
       !               comm_cart,requests(4),error)
       !call MPI_ISSEND(var(1,0,0),1,xhalo,left,0, &
       !               comm_cart,requests(3),error)
       !call MPI_WAITALL(4, requests, statuses, error)
    case(2) ! y direction
       call MPI_SENDRECV(var(0,1,0),1,yhalo,front,0, &
                         var(0,j1,0),1,yhalo,back,0, &
                         comm_cart,status,error)
       call MPI_SENDRECV(var(0,jmax,0),1,yhalo,back,0, &
                         var(0,0,0),1,yhalo,front,0,   &
                         comm_cart,status,error)
       !call MPI_IRECV(var(0,j1,0),1,yhalo,back,0, &
       !               comm_cart,requests(1),error)
       !call MPI_IRECV(var(0,0,0),1,yhalo,front,1, &
       !               comm_cart,requests(2),error)
       !call MPI_ISSEND(var(0,1,0),1,yhalo,front,0, &
       !               comm_cart,requests(3),error)
       !call MPI_ISSEND(var(0,jmax,0),1,yhalo,back,1, &
       !               comm_cart,requests(4),error)
       !call MPI_WAITALL(4, requests, statuses, error)
    end select

    return
  end subroutine updthalos



  subroutine updthalos_1d(var,dir)

    use mpi
    use mod_common_mpi

    real, dimension(0:), intent(inout) :: var
    integer, intent(in) :: dir
    integer :: requests(4), statuses(MPI_STATUS_SIZE,4)

    !  This subroutine updates the halos of a 1D array
    ! 
    !  It uses ?halo_1d defined in initmpi.f90

    select case(dir)
    case(1) ! x direction

       call MPI_SENDRECV(var(1),1,xhalo_1d,left, 0, &
            var(i1),1,xhalo_1d,right,0, &
            comm_cart,status,error)

       call MPI_SENDRECV(var(imax),1,xhalo_1d,right,0, &
            var(0),1,xhalo_1d,left,0, &
            comm_cart,status,error)

    case(2) ! y direction
       call MPI_SENDRECV(var(1),1,yhalo_1d,front,0, &
            var(j1),1,yhalo_1d,back,0, &
            comm_cart,status,error)

       call MPI_SENDRECV(var(jmax),1,yhalo_1d,back,0, &
            var(0),1,yhalo_1d,front,0, &
            comm_cart,status,error)

    end select

    return
  end subroutine updthalos_1d



  subroutine updthalos_ls(var,dir)

    use mpi
    use mod_common_mpi

    real, dimension(-2:,-2:,-2:), intent(inout) :: var
    integer, intent(in) :: dir
    integer :: requests(4), statuses(MPI_STATUS_SIZE,4)

    !  This subroutine updates the halos that store level set
    !  from the neighboring computational sub-domain
    ! 
    !  It uses ?halo_ls defined in initmpi.f90

    select case(dir)
    case(1) ! x direction

       call MPI_SENDRECV(var(1,-2,-2), 1,xhalo_ls,left, 0, &
            var(i1,-2,-2),1,xhalo_ls,right,0, &
            comm_cart,status,error)

       call MPI_SENDRECV(var(imax-2,-2,-2),1,xhalo_ls,right,0, &
            var(-2,-2,-2),    1,xhalo_ls,left, 0, &
            comm_cart,status,error)

    case(2) ! y direction
       call MPI_SENDRECV(var(-2,1,-2), 1,yhalo_ls,front,0, &
            var(-2,j1,-2),1,yhalo_ls,back, 0, &
            comm_cart,status,error)

       call MPI_SENDRECV(var(-2,jmax-2,-2),1,yhalo_ls,back, 0, &
            var(-2,-2,-2),    1,yhalo_ls,front,0, &
            comm_cart,status,error)

    end select

    return
  end subroutine updthalos_ls


  subroutine updthalos_mls(var,dir,l)

    use mpi
    use mod_common_mpi

    real, dimension(-2:,-2:,-2:,1:), intent(inout) :: var
    integer, intent(in) :: dir, l
    integer :: requests(4), statuses(MPI_STATUS_SIZE,4)

    !  This subroutine updates the halos that store multiple level set
    !  from the neighboring computational sub-domain
    ! 
    !  It uses ?halo_ls defined in initmpi.f90

    select case(dir)
    case(1) ! x direction

       call MPI_SENDRECV(var(1,-2,-2,l), 1,xhalo_ls,left, 0, &
            var(i1,-2,-2,l),1,xhalo_ls,right,0, &
            comm_cart,status,error)

       call MPI_SENDRECV(var(imax-2,-2,-2,l),1,xhalo_ls,right,0, &
            var(-2,-2,-2,l),    1,xhalo_ls,left, 0, &
            comm_cart,status,error)

    case(2) ! y direction
       call MPI_SENDRECV(var(-2,1,-2,l), 1,yhalo_ls,front,0, &
            var(-2,j1,-2,l),1,yhalo_ls,back, 0, &
            comm_cart,status,error)

       call MPI_SENDRECV(var(-2,jmax-2,-2,l),1,yhalo_ls,back, 0, &
            var(-2,-2,-2,l),    1,yhalo_ls,front,0, &
            comm_cart,status,error)

    end select

    return
  end subroutine updthalos_mls


end module mod_bound
