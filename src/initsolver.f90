! Initialization of the pressure solver using FFT/cosine transform.
!
! Author: Anthony Ge (KTH)
! Date:   25 January 2017
!------------------------------------------------------------------

module mod_initsolver

  use mod_solver_common
  use mod_param
  use decomp_2d
  
  implicit none
  
  private
  public initsolver
  
contains

  subroutine initsolver
    
    integer :: i,j,k,iv,jv
    real :: xrt(itot), yrt(jtot), zrt(kmax)
    real :: bb

    select case(BC_in_y)

    case('Periodic')

       ! Generate tridiagonal matrix

       allocate (a(kmax))
       allocate (b(kmax))
       allocate (c(kmax))
       
       do k=1,kmax
          a(k) = dzi*dzi
          c(k) = dzi*dzi
          b(k) = -(a(k) + c(k))
       enddo

       ! Neumann boundary condition for correction pressure in z-direction;
       ! consistent with no-penetration condition for prediction velocity 
       ! at solid walls (Channel)

       b(1) = b(1) + a(1)
       b(kmax) = b(kmax) + c(kmax)
       a(1) = 0.
       c(kmax) = 0.

       ! set lookup tables.

       call vrffti(itot,wi)
       call vrffti(jtot,wj)

       ! generate eigenvalues ( xrt and yrt ).

       ! x direction

       do i=3,itot,2
          xrt(i-1) = -4.*dxi*dxi*(sin(float((i-1))*pi/(2.*itot)))**2
          xrt(i) = xrt(i-1)
       enddo
       xrt(1   ) = 0.
       xrt(itot) = -4.*dxi*dxi

       ! y direction

       do j=3,jtot,2
          yrt(j-1) = -4.*dyi*dyi*(sin(float((j-1))*pi/(2.*jtot)))**2
          yrt(j  ) = yrt(j-1)
       enddo
       yrt(1   ) = 0.
       yrt(jtot) = -4.*dyi*dyi

       do j=1,jmax
          jv = j + zstart(2) - 1
          do i=1,imax
             iv = i + zstart(1) - 1
             xyrt(i,j) = xrt(iv)+yrt(jv)
          enddo
       enddo

    case('Walls','InOut','OutOut')

       ! Generate tridiagonal matrix

       allocate (a(jtot))
       allocate (b(jtot))
       allocate (c(jtot))
       
       do j=1,jtot
          a(j) = dyi*dyi
          c(j) = dyi*dyi
          b(j) = -(a(j) + c(j))
       enddo

       ! Neumann boundary condition for correction pressure in z-direction;
       ! consistent with no-penetration condition for prediction velocity 
       ! at solid walls (Channel)

       select case(BC_in_y)
       case('Walls')

          b(1) = b(1) + a(1)
          b(jtot) = b(jtot) + c(jtot)
          a(1) = 0.
          c(jtot) = 0.

       case('InOut')

          b(1) = b(1) + a(1)
          b(jtot) = b(jtot) - c(jtot)
          a(1) = 0.
          c(jtot) = 0.

       case('OutOut') !####

          b(1) = b(1) - a(1)
          b(jtot) = b(jtot) - c(jtot)
          a(1) = 0.
          c(jtot) = 0.

       end select

       ! set lookup tables.

       call vrffti(itot,wi)
       call vcosqi(kmax,wk)

       ! generate eigenvalues ( xrt and yrt ).

       ! x direction

       do i=3,itot,2
          xrt(i-1) = -4.*dxi*dxi*(sin(float((i-1))*pi/(2.*itot)))**2
          xrt(i) = xrt(i-1)
       enddo
       xrt(1   ) = 0.
       xrt(itot) = -4.*dxi*dxi

       ! z direction

       do k=1,kmax
          zrt(k) = (dzi**2)*( -2. + 2.*cos(real(k-1)*pi/kmax) )
       enddo

       do j=1,jmax
          do k=1,kmax
             do i=1,imax
                iv = i + zstart(1) - 1
                xzrt(i,j,k) = xrt(iv)+zrt(k)
             enddo
          enddo
       enddo

    end select
    
    return
  end subroutine initsolver

end module mod_initsolver
