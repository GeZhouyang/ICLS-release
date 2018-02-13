module mod_post ! BCs missing
!
! Probably not bug-free
!
use mpi
use mod_param
use mod_common
use mod_common_mpi
use decomp_2d
implicit none
private
public vorticity, enstrophy, strain_rate, q_criterion, compute_r, swirl
contains
!
subroutine vorticity(ux,uy,uz,vox,voy,voz)
implicit none
real, intent(in), dimension(1:,1:,1:) :: ux,uy,uz
real, intent(out), dimension(imax,jmax,kmax) :: vox,voy,voz
integer :: i,j,k
!
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      !
      ! x component of the vorticity at cell center
      !
      vox(i,j,k) = 0.25*( &
                         (uz(i,j+1,k  )-uz(i,j  ,k  ))*dyi - (uy(i,j  ,k+1)-uy(i,j  ,k  ))*dzi + &
                         (uz(i,j+1,k-1)-uz(i,j  ,k-1))*dyi - (uy(i,j  ,k  )-uy(i,j  ,k-1))*dzi + &
                         (uz(i,j  ,k  )-uz(i,j-1,k  ))*dyi - (uy(i,j-1,k+1)-uy(i,j-1,k  ))*dzi + &
                         (uz(i,j  ,k-1)-uz(i,j-1,k-1))*dyi - (uy(i,j-1,k  )-uy(i,j-1,k-1))*dzi &
                        )
      !
      ! y component of the vorticity at cell center
      !
      voy(i,j,k) = 0.25*( &
                         (ux(i  ,j,k+1)-ux(i  ,j,k  ))*dzi - (uz(i+1,j,k  )-uz(i  ,j,k  ))*dxi + &
                         (ux(i  ,j,k  )-ux(i  ,j,k-1))*dzi - (uz(i+1,j,k-1)-uz(i  ,j,k-1))*dxi + &
                         (ux(i-1,j,k+1)-ux(i-1,j,k  ))*dzi - (uz(i  ,j,k  )-uz(i-1,j,k  ))*dxi + &
                         (ux(i-1,j,k  )-ux(i-1,j,k-1))*dzi - (uz(i  ,j,k-1)-uz(i-1,j,k-1))*dxi &
                        )
      !
      ! z component of the vorticity at cell center
      !
      voz(i,j,k) = 0.25*( &
                         ((uy(i+1,j  ,k)-uy(i  ,j  ,k))*dxi - ux(i  ,j+1,k)-ux(i  ,j  ,k))*dyi + &
                         ((uy(i+1,j-1,k)-uy(i  ,j-1,k))*dxi - ux(i  ,j  ,k)-ux(i ,j-1 ,k))*dyi + &
                         ((uy(i  ,j  ,k)-uy(i-1,j  ,k))*dxi - ux(i-1,j+1,k)-ux(i-1,j  ,k))*dyi + &
                         ((uy(i  ,j-1,k)-uy(i-1,j-1,k))*dxi - ux(i-1,j  ,k)-ux(i-1,j-1,k))*dyi &
                        )
    enddo
  enddo
enddo
!
return
end subroutine vorticity
!
subroutine strain_rate(ux,uy,uz,str)
implicit none
real, intent(in), dimension(1:,1:,1:) :: ux,uy,uz
real, intent(out), dimension(imax,jmax,kmax) :: str
real :: s11,s22,s33,s12,s13,s23
integer :: i,j,k
!
! compute sijsij, where sij = (1/2)(du_i/dx_j + du_j/dx_i)
!
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      s11 = ((ux(i,j,k)-ux(i-1,j,k))*dxi)**2.
      s22 = ((uy(i,j,k)-uy(i,j-1,k))*dyi)**2.
      s33 = ((uz(i,j,k)-uz(i,j,k-1))*dzi)**2.
      s12 = .25*( &
                 ((ux(i  ,j+1,k)-ux(i  ,j  ,k))*dyi + (uy(i+1,j  ,k)-uy(i  ,j  ,k))*dxi)**2. + &
                 ((ux(i  ,j  ,k)-ux(i ,j-1 ,k))*dyi + (uy(i+1,j-1,k)-uy(i  ,j-1,k))*dxi)**2. + &
                 ((ux(i-1,j+1,k)-ux(i-1,j  ,k))*dyi + (uy(i  ,j  ,k)-uy(i-1,j  ,k))*dxi)**2. + &
                 ((ux(i-1,j  ,k)-ux(i-1,j-1,k))*dyi + (uy(i  ,j-1,k)-uy(i-1,j-1,k))*dxi)**2. &
                )*.25
      s13 = .25*( &
                 ((ux(i  ,j,k+1)-ux(i  ,j,k  ))*dzi + (uz(i+1,j,k  )-uz(i  ,j,k  ))*dxi)**2. + &
                 ((ux(i  ,j,k  )-ux(i  ,j,k-1))*dzi + (uz(i+1,j,k-1)-uz(i  ,j,k-1))*dxi)**2. + &
                 ((ux(i-1,j,k+1)-ux(i-1,j,k  ))*dzi + (uz(i  ,j,k  )-uz(i-1,j,k  ))*dxi)**2. + &
                 ((ux(i-1,j,k  )-ux(i-1,j,k-1))*dzi + (uz(i  ,j,k-1)-uz(i-1,j,k-1))*dxi)**2. &
                )*.25
      s23 = .25*( &
                 ((uy(i,j  ,k+1)-uy(i,j  ,k  ))*dzi + (uz(i,j+1,k  )-uz(i,j  ,k  ))*dyi)**2. + &
                 ((uy(i,j  ,k  )-uy(i,j  ,k-1))*dzi + (uz(i,j+1,k-1)-uz(i,j  ,k-1))*dyi)**2. + &
                 ((uy(i,j-1,k+1)-uy(i,j-1,k  ))*dzi + (uz(i,j  ,k  )-uz(i,j-1,k  ))*dyi)**2. + &
                 ((uy(i,j-1,k  )-uy(i,j-1,k-1))*dzi + (uz(i,j  ,k-1)-uz(i,j-1,k-1))*dyi)**2. &
                )*.25
      str(i,j,k) = s11+s22+s33 + 2*(s12+s13+s23)
    enddo
  enddo
enddo
!
return
end subroutine strain_rate
!

subroutine enstrophy(ux,uy,uz,ens)
implicit none
real, intent(in), dimension(1:,1:,1:) :: ux,uy,uz
real, intent(out), dimension(imax,jmax,kmax) :: ens
real :: e11,e22,e33,e12,e13,e23
integer :: i,j,k
!
! compute wijwij, where wij = (1/2)(du_i/dx_j - du_j/dx_i)
!
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      e12 = .25*( &
                 ((ux(i  ,j+1,k)-ux(i  ,j  ,k))*dyi - (uy(i+1,j  ,k)-uy(i  ,j  ,k))*dxi)**2. + &
                 ((ux(i  ,j  ,k)-ux(i ,j-1 ,k))*dyi - (uy(i+1,j-1,k)-uy(i  ,j-1,k))*dxi)**2. + &
                 ((ux(i-1,j+1,k)-ux(i-1,j  ,k))*dyi - (uy(i  ,j  ,k)-uy(i-1,j  ,k))*dxi)**2. + &
                 ((ux(i-1,j  ,k)-ux(i-1,j-1,k))*dyi - (uy(i  ,j-1,k)-uy(i-1,j-1,k))*dxi)**2. &
                )*.25
      e13 = .25*( &
                 ((ux(i  ,j,k+1)-ux(i  ,j,k  ))*dzi - (uz(i+1,j,k  )-uz(i  ,j,k  ))*dxi)**2. + &
                 ((ux(i  ,j,k  )-ux(i  ,j,k-1))*dzi - (uz(i+1,j,k-1)-uz(i  ,j,k-1))*dxi)**2. + &
                 ((ux(i-1,j,k+1)-ux(i-1,j,k  ))*dzi - (uz(i  ,j,k  )-uz(i-1,j,k  ))*dxi)**2. + &
                 ((ux(i-1,j,k  )-ux(i-1,j,k-1))*dzi - (uz(i  ,j,k-1)-uz(i-1,j,k-1))*dxi)**2. &
                )*.25
      e23 = .25*( &
                 ((uy(i,j  ,k+1)-uy(i,j  ,k  ))*dzi - (uz(i,j+1,k  )-uz(i,j  ,k  ))*dyi)**2. + &
                 ((uy(i,j  ,k  )-uy(i,j  ,k-1))*dzi - (uz(i,j+1,k-1)-uz(i,j  ,k-1))*dyi)**2. + &
                 ((uy(i,j-1,k+1)-uy(i,j-1,k  ))*dzi - (uz(i,j  ,k  )-uz(i,j-1,k  ))*dyi)**2. + &
                 ((uy(i,j-1,k  )-uy(i,j-1,k-1))*dzi - (uz(i,j  ,k-1)-uz(i,j-1,k-1))*dyi)**2. &
                )*.25
      ens(i,j,k) =  2*(e12+e13+e23)
    enddo
  enddo
enddo
!
return
end subroutine enstrophy

subroutine q_criterion(ens,str,qq)
implicit none
real, intent(in), dimension(1:,1:,1:) :: ens, str
real, intent(out), dimension(imax,jmax,kmax) :: qq
real :: s11,s22,s33,s12,s13,s23
integer :: i,j,k
!
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      qq(i,j,k) = .5*(ens(i,j,k)-str(i,j,k))
    enddo
  enddo
enddo
!
return
end subroutine q_criterion
subroutine compute_r(ux,uy,uz,rr)
implicit none
real, intent(in), dimension(1:,1:,1:) :: ux,uy,uz
real, intent(out), dimension(imax,jmax,kmax) :: rr
real :: s11,s22,s33,s12,s13,s23,vox,voy,voz
real :: s21,s31,s32
integer :: i,j,k
!
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      s11 = ((ux(i,j,k)-ux(i-1,j,k))*dxi)
      s22 = ((uy(i,j,k)-uy(i,j-1,k))*dyi)
      s33 = ((uz(i,j,k)-uz(i,j,k-1))*dzi)
      s12 = .25*( &
                 ((ux(i  ,j+1,k)-ux(i  ,j  ,k))*dyi + (uy(i+1,j  ,k)-uy(i  ,j  ,k))*dxi) + &
                 ((ux(i  ,j  ,k)-ux(i ,j-1 ,k))*dyi + (uy(i+1,j-1,k)-uy(i  ,j-1,k))*dxi) + &
                 ((ux(i-1,j+1,k)-ux(i-1,j  ,k))*dyi + (uy(i  ,j  ,k)-uy(i-1,j  ,k))*dxi) + &
                 ((ux(i-1,j  ,k)-ux(i-1,j-1,k))*dyi + (uy(i  ,j-1,k)-uy(i-1,j-1,k))*dxi) &
                )*.25
      s13 = .25*( &
                 ((ux(i  ,j,k+1)-ux(i  ,j,k  ))*dzi + (uz(i+1,j,k  )-uz(i  ,j,k  ))*dxi) + &
                 ((ux(i  ,j,k  )-ux(i  ,j,k-1))*dzi + (uz(i+1,j,k-1)-uz(i  ,j,k-1))*dxi) + &
                 ((ux(i-1,j,k+1)-ux(i-1,j,k  ))*dzi + (uz(i  ,j,k  )-uz(i-1,j,k  ))*dxi) + &
                 ((ux(i-1,j,k  )-ux(i-1,j,k-1))*dzi + (uz(i  ,j,k-1)-uz(i-1,j,k-1))*dxi) &
                )*.25
      s23 = .25*( &
                 ((uy(i,j  ,k+1)-uy(i,j  ,k  ))*dzi + (uz(i,j+1,k  )-uz(i,j  ,k  ))*dyi) + &
                 ((uy(i,j  ,k  )-uy(i,j  ,k-1))*dzi + (uz(i,j+1,k-1)-uz(i,j  ,k-1))*dyi) + &
                 ((uy(i,j-1,k+1)-uy(i,j-1,k  ))*dzi + (uz(i,j  ,k  )-uz(i,j-1,k  ))*dyi) + &
                 ((uy(i,j-1,k  )-uy(i,j-1,k-1))*dzi + (uz(i,j  ,k-1)-uz(i,j-1,k-1))*dyi) &
                )*.25
      s21 = s12
      s31 = s13
      s32 = s23
      !
      ! x component of the vorticity at cell center
      !
      vox = 0.25*( &
                         (uz(i,j+1,k  )-uz(i,j  ,k  ))*dyi - (uy(i,j  ,k+1)-uy(i,j  ,k  ))*dzi + &
                         (uz(i,j+1,k-1)-uz(i,j  ,k-1))*dyi - (uy(i,j  ,k  )-uy(i,j  ,k-1))*dzi + &
                         (uz(i,j  ,k  )-uz(i,j-1,k  ))*dyi - (uy(i,j-1,k+1)-uy(i,j-1,k  ))*dzi + &
                         (uz(i,j  ,k-1)-uz(i,j-1,k-1))*dyi - (uy(i,j-1,k  )-uy(i,j-1,k-1))*dzi &
                        )
      !
      ! y component of the vorticity at cell center
      !
      voy = 0.25*( &
                         (ux(i  ,j,k+1)-ux(i  ,j,k  ))*dzi - (uz(i+1,j,k  )-uz(i  ,j,k  ))*dxi + &
                         (ux(i  ,j,k  )-ux(i  ,j,k-1))*dzi - (uz(i+1,j,k-1)-uz(i  ,j,k-1))*dxi + &
                         (ux(i-1,j,k+1)-ux(i-1,j,k  ))*dzi - (uz(i  ,j,k  )-uz(i-1,j,k  ))*dxi + &
                         (ux(i-1,j,k  )-ux(i-1,j,k-1))*dzi - (uz(i  ,j,k-1)-uz(i-1,j,k-1))*dxi &
                        )
      !
      ! z component of the vorticity at cell center
      !
      voz = 0.25*( &
                         ((uy(i+1,j  ,k)-uy(i  ,j  ,k))*dxi - ux(i  ,j+1,k)-ux(i  ,j  ,k))*dyi + &
                         ((uy(i+1,j-1,k)-uy(i  ,j-1,k))*dxi - ux(i  ,j  ,k)-ux(i ,j-1 ,k))*dyi + &
                         ((uy(i  ,j  ,k)-uy(i-1,j  ,k))*dxi - ux(i-1,j+1,k)-ux(i-1,j  ,k))*dyi + &
                         ((uy(i  ,j-1,k)-uy(i-1,j-1,k))*dxi - ux(i-1,j  ,k)-ux(i-1,j-1,k))*dyi &
                        )
      rr(i,j,k) = (1./3.)*(s11**3.+s22**3.+s33**3.) +       &
                  s33*s32*s23 + s22*s23*s32 + s11*s12*s21 + &
                  s11*s13*s31 + s12*s22*s21 + s13*s33*s31 + &
                  s13*s32*s21 + s12*s23*s31 +               &
                  (3./4.)*( vox**2.*s11 + voy**2.*s22 + voz**2.*s33 + &
                            vox*voy*s12 + vox*voz*s13 + voy*vox*s21 + voy*voz*s23 + voz*vox*s31 + voz*voy*s32 &
                          )
    enddo
  enddo
enddo
!
return
end subroutine compute_r
!
subroutine swirl(qq,rr,ss)
implicit none
real, intent(in), dimension(1:,1:,1:) :: qq,rr
real, intent(out), dimension(imax,jmax,kmax) :: ss
integer :: i,j,k
real :: discr,s,t
!
! compute , lambda_ci 
!
ss = 0.
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      discr = qq(i,j,k)**3.+rr(i,j,k)**2.
      if(discr.gt.0) then
        s = (rr(i,j,k)+sqrt(discr))**(1./3.)
        t = (rr(i,j,k)-sqrt(discr))**(1./3.)
        ss(i,j,k) = .5*sqrt(3.)*(s-t)
      endif
    enddo
  enddo
enddo
!
return
end subroutine swirl
!
end module mod_post
