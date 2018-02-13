module mod_correc

  use mod_param
  use mod_common

  implicit none
  
  private
  public correc
  
contains
  !
  ! Corrects the velocity so that it is divergence free
  !
  subroutine correc(p)

    real, intent(in), dimension(0:,0:,0:) :: p
    real, dimension(0:i1,0:j1,0:k1) :: dpdx, dpdy, dpdz
    integer :: i,j,k

    !$omp parallel default(none) &
    !$omp&shared(factori,factorj,factork) &
    !$omp&shared(p,unew,vnew,wnew,dudt,dvdt,dwdt) &
    !$omp&private(i,j,k)
    !$omp do

    do k=1,kmax
       do j=1,jmax
          do i=1,imax

             !---- Remove pressure jump (GFM) ----!

             unew(i,j,k) = dudt(i,j,k) - dt*dxi * ( p(i+1,j,k)-p(i,j,k) -p_x(i,j,k) )/rho0 
             vnew(i,j,k) = dvdt(i,j,k) - dt*dyi * ( p(i,j+1,k)-p(i,j,k) -p_y(i,j,k) )/rho0
             wnew(i,j,k) = dwdt(i,j,k) - dt*dzi * ( p(i,j,k+1)-p(i,j,k) -p_z(i,j,k) )/rho0

          enddo
       enddo
    enddo
    !$omp end parallel


    return
  end subroutine correc


end module mod_correc
