module mod_fillps

  use mod_param
  use mod_common

  implicit none

  private
  public fillps

contains

  subroutine fillps(p)
    !                                                                       
    ! Fill the RHS of the Poisson equation for the correction pressure.     
    !                                                                       
    ! The discrete divergence is:                                           
    !                                                                       
    ! w(i,j,k)-w(i,j,k-1)   v(i,j,k)-v(i,j-1,k)   u(i,j,k)-u(i-1,j,k)       
    ! ------------------- + ------------------- + -------------------  = div
    !         dz                    dy                    dx                
    !                                                                       
    ! Note: in this subroutine p is not the correction pressure, but        
    ! the rhs of the Poisson equation, i.e.                                 
    !                                                                       
    ! p = RHS = rho0 / dt * div                                             
    !-----------------------------------------------------------------------

    implicit none

    integer :: i,j,k,im,jm,km
    real :: coef, vel_x,vel_y,vel_z
    real, intent(out), dimension(0:,0:,0:) :: p


    coef = rho0/dt

    do k=1,kmax
       km = k-1
       do j=1,jmax
          jm = j-1
          do i=1,imax
             im = i-1

             vel_x = (dudt(i,j,k)-dudt(im,j,k))*dxi
             vel_y = (dvdt(i,j,k)-dvdt(i,jm,k))*dyi
             vel_z = (dwdt(i,j,k)-dwdt(i,j,km))*dzi
             
             p(i,j,k) = coef*(vel_x + vel_y + vel_z)

          enddo
       enddo
    enddo


    return
  end subroutine fillps


end module mod_fillps
