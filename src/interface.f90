module mod_interface

  use mod_common
  use mod_common_mpi
  use mod_param
  use mod_misc

  implicit none
  
  private
  public get_curvature, re_curvature, icurvature, exact_curvature, ipot, Least_Squares
  
contains


  !-------------
  !-------------
  
  
  subroutine get_curvature(phi,curvature)  !# fast version

    !---- Compute curvature at time level n+1 (1st order on the ø-grid, using Marchandise JCP2007): 
    ! 
    !                grad ø
    !     k = - div ------- 
    !                 |ø|
    ! 
    !     Negative curvature for a circle

    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in)  :: phi
    real, dimension(0:i1,0:j1,0:k1),          intent(out) :: curvature

    integer :: i, j, k, n
    integer :: dist_x, dist_y, dist_z, counter

    real :: max_kappa, norm_val, term1, term2, term3, dirac_phi
    real :: dphidx, dphidy, dphidz
    real :: d2phidx, d2phidy, d2phidz
    real :: d2phidxdy, d2phidxdz, d2phidydz

    real, dimension(10) :: x ! smoothed level set and its derivatives

    !#
    normal%m = 1.
    normal%x = 0.
    normal%y = 0.
    normal%z = 0.
    !#
    curvature(:,:,:) = 0.

    max_kappa = 1./dz

    do k=1,kmax
       do j=1,jmax
          do i=1,imax

             if (abs(phi(i,j,k)) .le. NarrowBand_2 -1.*dz) then ! wide support for mass-correction
 
                x = matmul(pseudo_inv_A, reshape(phi(i-1:i+1,j-1:j+1,k-1:k+1), (/27/)))

                dphidx = x(2)
                dphidy = x(3)
                dphidz = x(4)

                d2phidx = x(5)
                d2phidy = x(6)
                d2phidz = x(7)

                d2phidxdy = x(8)
                d2phidxdz = x(9)
                d2phidydz = x(10)

                norm_val = sqrt( dphidx**2 &
                               + dphidy**2 &
                               + dphidz**2 )

                !#
                normal(i,j,k)%m = norm_val
                normal(i,j,k)%x = dphidx/norm_val
                normal(i,j,k)%y = dphidy/norm_val
                normal(i,j,k)%z = dphidz/norm_val
                !#

         
                !         if (norm_val**2 .lt. 1.e-10) then
                if (norm_val .lt. 0.9) then  !# the first filter for bad cells
                   curvature(i,j,k) = 0.
                else
                   term1 = (d2phidx + d2phidy + d2phidz)/norm_val

                   term2 = dphidx*dphidx*d2phidx  + dphidy*dphidy*d2phidy  + dphidz*dphidz*d2phidz
                   term2 = term2/norm_val**3

                   term3 = dphidx*dphidy*d2phidxdy + dphidx*dphidz*d2phidxdz + dphidy*dphidz*d2phidydz
                   term3 = 2.*term3/norm_val**3

                   curvature(i,j,k) = -(term1-term2-term3)
                endif
 
                ! limiter due to resolution
                
                if ( curvature(i,j,k) .lt. -max_kappa ) then
                   curvature(i,j,k) = -max_kappa
                elseif ( curvature(i,j,k) .gt. max_kappa ) then
                   curvature(i,j,k) = max_kappa

                elseif ( abs(curvature(i,j,k)) .lt. 1e-9 ) then  ! for flat interface
                   curvature(i,j,k) = 0.
                endif
      
             endif

          enddo
       enddo
    enddo


    return
  end subroutine get_curvature


  !-------------
  !-------------
  
  
  subroutine icurvature(nb, ph,cv, icurv)

    ! Estimation of interface curvature (assuming constant curv nearby).
    ! Must ensure the same calculation for two adjacent cells,
    ! hence consistent in the projection/correction steps.

    integer, intent(in) :: nb
    real, dimension(-1:1), intent(in) :: ph, cv

    real, intent(out) :: icurv
    real :: sp, eps_0,eps_nb,eps_op
    real :: r_0,r_nb,r_op, r_i, a1,a2, b1,b2,b3

    if (cv(0) .lt. 0.) sp = -1.  ! concave
    if (cv(0) .gt. 0.) sp =  1.  ! convex

    eps_0  = sp*ph(0)
    eps_nb = sp*ph(nb)
    eps_op = sp*ph(-nb)

    if (TwoD) then

       if (cv(0)*cv(nb) .gt. 0.) then ! "good" cells

          r_0  = -1./cv(0)  + eps_0
          r_nb = -1./cv(nb) + eps_nb

          r_i = (r_0 + r_nb)/2.
          icurv = -1./r_i

       else !#

          icurv = 0.

       endif

    else  ! 3D

       if (cv(0)*cv(nb) .gt. 0.) then ! "good" cells

          a1 = eps_nb*cv(0)
          a2 = eps_0*cv(nb)
          icurv = (a1 - a2)/(eps_nb - eps_0)

       else !#

          icurv = 0.

       endif

    endif

    return
  end subroutine icurvature


  !-------------
  !-------------
  
  
  subroutine re_curvature(phi,curv,re_curv)

    !-- Reconstruct curvature from the least squares of 27 nodal values
    !-- so that it's less sensitive to small oscillations

    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in)  :: phi
    real, dimension(0:i1,0:j1,0:k1),          intent(in)  :: curv
    real, dimension(0:i1,0:j1,0:k1),          intent(out) :: re_curv

    integer :: i, j, k
    real, dimension(10) :: x ! array of reconstructed values


    re_curv(:,:,:) = 0.

    do k=1,kmax
       do j=1,jmax
          do i=1,imax

             if (abs(phi(i,j,k)) .lt. 2.*dz) then ! one more layer than pressure jump

                x = matmul(pseudo_inv_A, reshape(curv(i-1:i+1,j-1:j+1,k-1:k+1), (/27/)))

                re_curv(i,j,k) = x(1)

             endif
          enddo
       enddo
    enddo


    return
  end subroutine re_curvature


  !-------------
  !-------------
  
  
  subroutine ipot

    !-- Find the minimal distance, dmin, between two interfaces, m and n.

    integer :: i,j,k,l,m,n
    real :: cutoff, d1

    real, dimension(lmax,lmax) :: dtemp
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phim,phin,super_phi

    dtemp = 0.
    dmin = 0.

    cutoff = NarrowBand_2  ! max distance considered

    do m = 1,lmax
       phim = lvset(:,:,:,m)  ! buffer

       do n = 1,lmax
          phin = lvset(:,:,:,n)  ! buffer

          if (n .ne. m) then

             super_phi = 1e3  ! initialize superimposed level set w/ a big number

             where( phim .gt.-1.*dz .and. &
                  phim .lt. 1.*dz .and. &
                  phin .lt. cutoff)

                super_phi = phim + phin

             end where

             d1 = minval(super_phi)
             call mpi_allreduce(mpi_in_place,d1,1,mpi_real8,mpi_min,comm_cart,error)
             dtemp(m,n) = d1

          endif

       enddo
    enddo

    do m = 1,lmax-1
       do n = m+1,lmax

          dmin(m,n) = 0.5*( dtemp(m,n)+dtemp(n,m) )

       enddo
    enddo


    return
  end subroutine ipot


  !-------------
  !-------------
  
  
  subroutine Least_Squares(A,x,b,nr,nc)

    !---- Construct level set function from 27 neighboring cells (including itself)
    !     to get 2nd order accurate curvature. (Copied from JC.)
    !     (Need Lapack library)
    ! 

    integer :: nr, nc
    real, dimension(nr,nc) :: A
    real, dimension(nr)    :: b
    real, dimension(nc)    :: x

    !----- Variables for LAPACK DGELS -----!

    character*1 :: TRANS = 'N'
    integer :: m
    integer :: n
    integer :: nrhs = 1
    real, dimension(nr,nc) :: A_
    integer :: LDA
    integer :: LDB
    integer :: LWORK
    real, dimension(:), allocatable :: WORK
    integer INFO

    !----- Solving the least square problem -----!
    
    m = nr
    n = nc

    A_ = A
    LDA = max(1,m)
    LDB = max(1,m,n)
    LWORK = max( 1, m*n + max(m*n,nrhs) )

    allocate(WORK(max(1,LWORK)))

    call DGELS( TRANS, M, N, NRHS, A_ , LDA, B, LDB, WORK, LWORK, INFO )

    deallocate(WORK)

    x = 0.
    x = b(1:nc)

  end subroutine Least_Squares


  !-------------
  !-------------
  
  
  subroutine exact_curvature(curvature)

    !---- Compute analytically the exact mean curvature of some simple shapes
    !     Circle  k = -1/R
    !     Sphere  k = -1/R1 - 1/R2

    real, dimension(0:i1,0:j1,0:k1), intent(out) :: curvature
    real :: radius

    radius = .2
    curvature = -2./radius  ! sphere

    return
  end subroutine exact_curvature


end module mod_interface
