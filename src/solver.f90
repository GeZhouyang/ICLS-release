!********************  DIRECT POISSON SOLVER **********************
!*****                                                        *****
!***               p_xx + p_yy + p_zz  = f(x,y,z)               ***
!****                                                         *****
!******************************************************************
! Fourier Transforms of the Poisson equation in the x and y
! directions yield:
!
! FFT_j[FTT_i[ p_xx + p_yy + p_zz ]] = FFT_j[FTT_i[ f(x,y,z) ]],
! a^2 p' + b^2 p' + p'_zz = f'(kx,ky,z),
!
! where a and b are the known eigenvalues, kx and ky are the
! wavenumbers in the x and y direction, and p'_zz is given by:
!
! p'_zz =[p'_{kx,ky,k+1} -2*p'_{kx,ky,k} +p'_{kx,ky,k-1}]/(dz*dz).
!
! The equation above results in a tridiagonal system in k:
!
! a^2 p' + b^2 p' + p'_zz =
! [p'_{kx,ky,k+1}-(2+a^2+b^2)*p'_{kx,ky,k}+p'_{kx,ky,k-1}]/(dz*dz)
! =f'(kx,ky,z).
!
! This can be solved with Gaussian elimination. The pressure p in
! physical space is obtained from 2 inverse Fourier Transforms
! according to: p = iFFT_j[ iFFT_i[ p' ] ].
!
!******************************************************************
!****   Programmers: Bendiks Jan Boersma                     ******
!****                Wim-Paul Breugem (1D parallellisation)  ******
!****                Pedro Costa (2D parallelisation)        ******
!****   email      : w.p.breugem@tudelft.nl                  ******
!****   Modified by: Anthony Ge (KTH) on 25 Jan 2017         ******
!****   USES       : VFFTPACK   (netlib)                     ******
!****                2DECOMP&FFT (www.2decomp.org)           ******
!******************************************************************

module mod_solver_common

  use mod_param

  implicit none

  real(mytype) :: wi(itot+15), wj(jtot+15), wk(2*kmax+15)
  real, dimension(imax,jmax) :: xyrt
  real, dimension(imax,jmax,kmax) :: xzrt
  
  real, dimension(:), allocatable :: a,b,c
  

end module mod_solver_common


!------------------------------------------------------------


module mod_solver
  
  use mod_param
  use mod_solver_common
  
  implicit none
  
  private
  public solver2d, solver2d_cos, solver1d
  
contains

  subroutine solver2d(pz)
    
    use decomp_2d

    real, intent(inout), dimension(1:,1:,1:) :: pz
    real, dimension(itot/dims(1),jtot,ktot/dims(2)) :: py
    real, dimension(itot,jtot/dims(1),ktot/dims(2)) :: px
    real :: bb
    real :: z,d(imax,jmax,kmax)
    real :: di(itot),dj(jtot)
    integer :: i,j,k

    !-- forward FFT in the x direction --!
    
    call transpose_z_to_y(pz,py)
    call transpose_y_to_x(py,px)

    !$omp parallel default(none) &
    !$omp& shared(px,xsize) private(i,j,k,di) &
    !$omp&firstprivate(wi)
    !$omp do
    do k=1,xsize(3)
       do j=1,xsize(2)
          call vrfftf(1,itot,px(1:itot,j,k),di,1,wi)
       enddo
    enddo
    !$omp end parallel

    !-- forward FFT in the y direction --!
    
    call transpose_x_to_y(px,py)

    !$omp parallel default(none) &
    !$omp& shared(py,ysize) private(i,j,k,dj) &
    !$omp&firstprivate(wj)
    !$omp do
    do k=1,ysize(3)
       do i=1,ysize(1)
          call vrfftf(1,jtot,py(i,1:jtot,k),dj,1,wj)
       enddo
    enddo
    !$omp end parallel

    !-- Gauss elimination in the z direction --!
    
    call transpose_y_to_z(py,pz)

    !$omp parallel default(none) &
    !$omp&shared(pz,d,a,b,c,xyrt) &
    !$omp&private(i,j,k,z,bb)
    !$omp do
    do j=1,jmax
       do i=1,imax
          z        = 1./(b(1)+xyrt(i,j))
          d(i,j,1) = c(1)*z
          pz(i,j,1) = pz(i,j,1)*z
       enddo
    enddo
    !$omp barrier
    do k=2,kmax-1
       !$omp do
       do j=1,jmax
          do i=1,imax
             bb       = b(k)+xyrt(i,j)
             z        = 1./(bb-a(k)*d(i,j,k-1))
             d(i,j,k) = c(k)*z
             pz(i,j,k) = (pz(i,j,k)-a(k)*pz(i,j,k-1))*z
          enddo
       enddo
       !$omp barrier
    enddo
    !$omp do
    do j=1,jmax
       do i=1,imax
          bb       = b(kmax)+xyrt(i,j)
          z        = bb-a(kmax)*d(i,j,kmax-1)
          if(z.ne.0.) then
             pz(i,j,kmax) = (pz(i,j,kmax)-a(kmax)*pz(i,j,kmax-1))/z
          else
             pz(i,j,kmax) =0.
          endif
       enddo
    enddo

    do k=kmax-1,1,-1
       !$omp do
       do j=1,jmax
          do i=1,imax
             pz(i,j,k) = pz(i,j,k)-d(i,j,k)*pz(i,j,k+1)
          enddo
       enddo
       !$omp barrier
    enddo
    !$omp end parallel

    !-- backward FFT in the y direction --!
    
    call transpose_z_to_y(pz,py)
    !$omp parallel default(none) &
    !$omp& shared(py,ysize) private(i,j,k,dj) &
    !$omp&firstprivate(wj)
    !$omp do
    do k=1,ysize(3)
       do i=1,ysize(1)
          call vrfftb(1,jtot,py(i,1:jtot,k),dj,1,wj)
       enddo
    enddo
    !$omp end parallel

    !-- backward FFT in the x direction --!
    
    call transpose_y_to_x(py,px)
    
    !$omp parallel default(none) &
    !$omp& shared(px,xsize) private(i,j,k,di) &
    !$omp&firstprivate(wi)
    !$omp do
    do k=1,xsize(3)
       do j=1,xsize(2)
          call vrfftb(1,itot,px(1:itot,j,k),di,1,wi)
       enddo
    enddo
    !$omp end parallel
    
    call transpose_x_to_y(px,py)
    call transpose_y_to_z(py,pz)
    
    return
  end subroutine solver2d
  

  !------------------------------------------------------------


  subroutine solver2d_cos(pz)
    
    use decomp_2d
    
    real, intent(inout), dimension(1:,1:,1:) :: pz
    real, dimension(itot/dims(1),jtot,ktot/dims(2)) :: py
    real, dimension(itot/dims(1),jtot,ktot/dims(2)) :: xzrt_y
    real, dimension(itot,jtot/dims(1),ktot/dims(2)) :: px
    real :: bb, z
    real, dimension(itot/dims(1),jtot,ktot/dims(2)) :: d
    real :: di(itot),dk(kmax)
    integer :: i,j,k
    
    !-- forward FFT in the x direction --!
    
    call transpose_z_to_y(pz,py)
    call transpose_y_to_x(py,px)

    !$omp parallel default(none) &
    !$omp& shared(px,xsize) private(i,j,k,di) &
    !$omp&firstprivate(wi)
    !$omp do
    do k=1,ktot/dims(2)
       do j=1,jtot/dims(1)
          call vrfftf(1,itot,px(1:itot,j,k),di,1,wi)
       enddo
    enddo
    !$omp end parallel

    !-- forward cosine transform in the z direction --!
    
    call transpose_x_to_y(px,py)
    call transpose_y_to_z(py,pz)
    
    !$omp parallel default(none) &
    !$omp& shared(py,ysize) private(i,j,k,dj) &
    !$omp&firstprivate(wj)
    !$omp do
    do j=1,jmax
       do i=1,imax
          call vcosqb(1,kmax,pz(i,j,1:kmax),dk,1,wk)
       enddo
    enddo
    !$omp end parallel

    !-- Gauss elimination in the y direction --!

    call transpose_z_to_y(pz,py)
    call transpose_z_to_y(xzrt,xzrt_y)

    do k=1,ktot/dims(2)
       do i=1,imax
          z        = 1./(b(1)+xzrt_y(i,1,k))
          d(i,1,k) = c(1)*z
          py(i,1,k) = py(i,1,k)*z
       enddo
    enddo

    do j=2,jtot-1     
       do k=1,ktot/dims(2)
          do i=1,imax
             bb       = b(j)+xzrt_y(i,1,k)
             z        = 1./(bb-a(j)*d(i,j-1,k))
             d(i,j,k) = c(j)*z
             py(i,j,k) = (py(i,j,k)-a(j)*py(i,j-1,k))*z
          enddo
       enddo
    enddo

    do k=1,ktot/dims(2)
       do i=1,imax
          bb       = b(jtot)+xzrt_y(i,1,k)
          z        = bb-a(jtot)*d(i,jtot-1,k)
          if(z.ne.0.) then
             py(i,jtot,k) = (py(i,jtot,k)-a(jtot)*py(i,jtot-1,k))/z
          else
             py(i,jtot,k) =0.
          endif
       enddo
    enddo

    do j=jtot-1,1,-1
       do k=1,ktot/dims(2)
          do i=1,imax
             py(i,j,k) = py(i,j,k)-d(i,j,k)*py(i,j+1,k)
          enddo
       enddo
    enddo

    !-- backward cosine transform in the z direction --!

    call transpose_y_to_z(py,pz)

    do j=1,jmax
       do i=1,imax
          call vcosqf(1,kmax,pz(i,j,1:kmax),dk,1,wk)
       enddo
    enddo

    !-- backward FFT in the x direction --!

    call transpose_z_to_y(pz,py)
    call transpose_y_to_x(py,px)

    do k=1,ktot/dims(2)
       do j=1,jtot/dims(1)
          call vrfftb(1,itot,px(1:itot,j,k),di,1,wi)
       enddo
    enddo

    call transpose_x_to_y(px,py)
    call transpose_y_to_z(py,pz)


    return
  end subroutine solver2d_cos


  !------------------------------------------------------------
  

subroutine solver1d(p)
use mod_zredistribute
implicit none
integer, parameter :: nprocs = dims(1)*dims(2)
integer, parameter :: ksol = kmax/nprocs
real, intent(inout) :: p(0:,0:,0:)
real,dimension(itot,jtot,ksol) :: p2
real :: bb
real :: z,d(imax,jmax,kmax)
real :: di(itot),dj(jtot)
integer i,j,k
!
! Redistribute the pressure so that it is only distributed
! in the z-direction (starting from 0 to nprocs)
!
call zredistribute(p,p2,0)
!$omp parallel default(none) &
!$omp&shared(p2) &
!$omp&private(i,j,k,di,dj) &
!$omp&firstprivate(wi,wj)
!$omp do
do k=1,ksol
  !  FFT  ---> I direction
  do j=1,jtot
    call vrfftf(1,itot,p2(1:itot,j,k),di,1,wi)
  enddo
  !  FFT  ---> J direction
  do i=1,itot
    call vrfftf(1,jtot,p2(i,1:jtot,k),dj,1,wj)
  enddo
enddo
!$omp end parallel
!
! Redistribute the pressure so that it is distributed
! on the cartesian grid
!
call zredistribute(p,p2,1)
!$omp parallel default(none) &
!$omp&shared(a,p,b,c,d,xyrt) &
!$omp&private(i,j,k,z,bb)
k=1
!$omp do
do j=1,jmax
  do i=1,imax
    z        = 1./(b(1)+xyrt(i,j))
    d(i,j,1) = c(1)*z
    p(i,j,1) = p(i,j,1)*z
  enddo
enddo
!$omp barrier
do k=2,kmax-1
!$omp do
  do j=1,jmax
    do i=1,imax
      bb       = b(k)+xyrt(i,j)
      z        = 1./(bb-a(k)*d(i,j,k-1))
      d(i,j,k) = c(k)*z
      p(i,j,k) = (p(i,j,k)-a(k)*p(i,j,k-1))*z
    enddo
  enddo
!$omp barrier
enddo
k = kmax
!$omp do
do j=1,jmax
  do i=1,imax
    bb       = b(kmax)+xyrt(i,j)
    z        = bb-a(kmax)*d(i,j,kmax-1)
    if(z.ne.0.) then
      p(i,j,kmax) = (p(i,j,kmax)-a(kmax)*p(i,j,kmax-1))/z
    else
      p(i,j,kmax) =0.
    endif
  enddo
enddo
do k=kmax-1,1,-1
!$omp do
  do j=1,jmax
    do i=1,imax
      p(i,j,k) = p(i,j,k)-d(i,j,k)*p(i,j,k+1)
    enddo
  enddo
!$omp barrier
enddo
!$omp end parallel
!
call zredistribute(p,p2,0)
!
!$omp parallel default(none) &
!$omp&shared(p2) &
!$omp&private(i,j,k,di,dj) &
!$omp&firstprivate(wi,wj)
!$omp do
do k=1,ksol
  ! BACKWARD FFT ---> J direction
  do i=1,itot
    call vrfftb(1,jtot,p2(i,1:jtot,k),dj,1,wj)
  enddo
  ! BACKWARD FFT ---> I direction
  do j=1,jtot
    call vrfftb(1,itot,p2(1:itot,j,k),di,1,wi)
  enddo
enddo
!$omp end parallel
call zredistribute(p,p2,1)
!
return
end subroutine solver1d


end module mod_solver
