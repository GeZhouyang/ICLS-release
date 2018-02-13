module mod_chkdt

implicit none

 private
 public chkdt

 contains


   subroutine chkdt(dtmax)

    use mpi
    use mod_param
    use mod_common
    use mod_common_mpi

    implicit none


      real, intent(inout) :: dtmax
      real :: speed_conv, speed_grav, speed_visc, speed_surf
      integer :: i,j,k
      real, dimension(1:i1,1:j1,1:k1) :: u,v,w
      real :: u_max,v_max,w_max, u_max_all,v_max_all,w_max_all
      real :: g1,g2,g3,g4,g5,g6, g_max,h_min
      real :: v1,v2
      real :: curv_max
      real :: s1,speed


      !-- Convective dt restriction --!

      u = 0.
      v = 0.
      w = 0.

      do k = 1,k1
        do j = 1,j1
          do i = 1,i1

            u(i,j,k) = abs( (unew(i,j,k) + unew(i-1,j,k))/2. )
            v(i,j,k) = abs( (vnew(i,j,k) + vnew(i,j-1,k))/2. )
            w(i,j,k) = abs( (wnew(i,j,k) + wnew(i,j,k-1))/2. )

          enddo
        enddo
      enddo

      u_max = maxval(u)
      v_max = maxval(v)
      w_max = maxval(w)

      call mpi_allreduce(u_max,u_max_all,1,mpi_real8,mpi_max,comm_cart,error)
      call mpi_allreduce(v_max,v_max_all,1,mpi_real8,mpi_max,comm_cart,error)
      call mpi_allreduce(w_max,w_max_all,1,mpi_real8,mpi_max,comm_cart,error)

      u_max = u_max_all
      v_max = v_max_all
      w_max = w_max_all

      speed_conv = u_max/dx + v_max/dy + w_max/dz


      !-- Gravitational dt restriction --!

      g1 = abs( iFr*(1.-rho_avr/rho1)*g_unit(1) )
      g2 = abs( iFr*(1.-rho_avr/rho2)*g_unit(1) )
      g3 = abs( iFr*(1.-rho_avr/rho1)*g_unit(2) )
      g4 = abs( iFr*(1.-rho_avr/rho2)*g_unit(2) )
      g5 = abs( iFr*(1.-rho_avr/rho1)*g_unit(3) )
      g6 = abs( iFr*(1.-rho_avr/rho2)*g_unit(3) )

      g_max = max(g1,g2,g3,g4,g5,g6)
      h_min = min(dx,dy,dz)

      speed_grav = sqrt(g_max/h_min)


      !-- Viscous dt restriction --!

      v1 = max( miu1/rho1,miu2/rho2 ) *visc
      v2 = 2./dx**2 + 2./dy**2 + 2./dz**2

      speed_visc = v1*v2


      !-- Capillary dt restriction --!

      if (TwoD) then
        curv_max = 1./h_min
      else
        curv_max = 2./h_min
      endif

      speed_surf = sqrt( curv_max/(We*rho0*h_min**2) )


      !-- Overall constraint --!

      s1 = speed_conv + speed_visc
      speed = ( s1 + sqrt(s1**2 + 4.*speed_grav**2 + 4.*speed_surf**2) )/2.

      dtmax = 0.5/speed



!!     Calculates the timestep dt based on stability conditions
!!     for AB2 scheme. (See: P. Wesseling, 'Principles of computational
!!     fluid dynamics'. Initially from Mehdi.)


!      real, intent(inout) :: dtmax
!      real ::   velo,dxi2
!      real ::   t1,t2,t3,n1,n2,n3,u2,v2,w2
!      real ::   dt1,dt1a,dt1b,dt2,dt2a,dt3,viscpm
!      real ::   dtmaxold,dtmax_all(2,1)
!      real ::   dtmaxu(2,1)
!      real ::   dtmaxv(2,1)
!      real ::   dtmaxw(2,1)
!      integer  ::  counter1,counter2,counter3
!      integer  ::  counter1_all,counter2_all,counter3_all
!      integer  ::  counter1u,counter2u,counter3u
!      integer  ::  counter1v,counter2v,counter3v
!      integer  ::  counter1w,counter2w,counter3w
!      integer  ::  maxi,maxj,maxk
!      integer  ::  maxiu,maxju,maxku
!      integer  ::  maxiv,maxjv,maxkv
!      integer  ::  maxiw,maxjw,maxkw
!      integer  ::  flag,flagmax
!      real     ::  flagmaxu,flagmaxv,flagmaxw
!      real     ::  dummy
!      integer  ::  i,j,k
!      real     ::  dtsurf


!!     spanwise velocity

!      dtmaxu(1,1) = 899999.
!      dtmaxu(2,1) = myid*1.
!      dtmaxold    = dtmaxu(1,1)
!      counter1    = 0
!      counter2    = 0
!      counter3    = 0
!      maxi        = 1000
!      maxj        = 1000
!      maxk        = 1000
!      flag        = 0
!      flagmax     = 0
!      do k=1,kmax
!        do j=1,jmax
!          do i=1,imax

!!           diffusion

!            dxi2  = dxi**2 + dyi**2 + dzi**2
!            dt1   = 1./( 4.*visc*dxi2 )
!            dt1a  = 1./( 4.*visc2*dxi2 ) !# 2nd phase
!            dt1b  = Reb*dx*dx/6.         !# Dodd & Ferrante
!            dt1   = min(dt1,dt1a,dt1b)

!!           convection

!            u2 = ( unew(i,j,k) )**2.
!            v2 = ( 0.25*(vnew(i,j,k)+vnew(i,j-1,k)+vnew(i+1,j,k)+vnew(i+1,j-1,k)) )**2.
!            w2 = ( 0.25*(wnew(i,j,k)+wnew(i,j,k-1)+wnew(i+1,j,k)+wnew(i+1,j,k-1)) )**2.
!            velo = sqrt(u2+v2+w2)
!            if (abs(1.*velo) .gt. 1.e-12) then
!              dt3 = 1./( sqrt(3.)*velo*sqrt(dxi2) )
!            else
!              dt3 = 99999.
!            endif

!!           Courant-number weighted with cell-Peclet number
!!           sum (((Pe_i)**(1./3.))*u_i*dt/dx_i) < 2**(1./3.) = 1.26
!!           With Pe_i < 2 (?) --> sum (u_i*dt/dx_i) < 1
!!           This condition is therefore expected to be less
!!           restrictive than condition 2.

!            t1 = ( abs( unew(i,j,k) ) )**(4./3.)
!            t2 = ( abs( 0.25*(vnew(i,j,k)+vnew(i,j-1,k)+vnew(i+1,j,k)+vnew(i+1,j-1,k))) )**(4./3.)
!            t3 = ( abs( 0.25*(wnew(i,j,k)+wnew(i,j,k-1)+wnew(i+1,j,k)+wnew(i+1,j,k-1))) )**(4./3.)
!            if (abs(1.*velo) .gt. 1.e-12) then
!              dt2 = ((2.*visc)**(1./3.))/ &
!                    (t1*(dxi**(2./3.)) + t2*(dyi**(2./3.)) + t3*(dyi**(2./3.)))
!              dt2a= ((2.*visc2)**(1./3.))/ &
!                    (t1*(dxi**(2./3.)) + t2*(dyi**(2./3.)) + t3*(dyi**(2./3.))) !# 2nd phase
!	      dt2   = min(dt2,dt2a)
!            else
!              dt2 = 99999.
!            endif
!            if (dt1 .lt. dt2) then
!              if (dt1 .lt. dt3) then
!                flag = 1
!                counter1 = counter1 + 1
!              else
!                flag = 3
!                counter3 = counter3 + 1
!              endif
!            else
!              if (dt2 .lt. dt3) then
!                flag = 2
!                counter2 = counter2 + 1
!              else
!                flag = 3
!                counter3 = counter3 + 1
!              endif
!            endif
!            dummy       = dtmaxu(1,1)
!            dtmaxu(1,1) = min(dummy,dt1,dt2,dt3)
!            if (dtmaxu(1,1) .lt. dtmaxold) then
!              maxi     = i
!              maxj     = j
!              maxk     = k
!              flagmax  = flag
!              dtmaxold = dtmaxu(1,1)
!            endif
!          enddo
!        enddo
!      enddo
!      call mpi_allreduce(dtmaxu,dtmax_all,1,mpi_2double_precision,mpi_minloc,comm_cart,error)
!      call mpi_allreduce(counter1,counter1_all,1,mpi_integer,mpi_sum,comm_cart,error)
!      call mpi_allreduce(counter2,counter2_all,1,mpi_integer,mpi_sum,comm_cart,error)
!      call mpi_allreduce(counter3,counter3_all,1,mpi_integer,mpi_sum,comm_cart,error)
!      dtmaxu(1,1) = dtmax_all(1,1)
!      dtmaxu(2,1) = dtmax_all(2,1)
!      maxiu       = maxi
!      maxju       = maxj
!      maxku       = maxk
!      counter1u   = counter1_all
!      counter2u   = counter2_all
!      counter3u   = counter3_all
!      flagmaxu    = flagmax
!!
!!     streamwise velocity
!!
!      dtmaxv(1,1) = 899999.
!      dtmaxv(2,1) = myid*1.
!      dtmaxold    = dtmaxv(1,1)
!      counter1    = 0
!      counter2    = 0
!      counter3    = 0
!      maxi        = 1000
!      maxj        = 1000
!      maxk        = 1000
!      flag        = 0
!      flagmax     = 0
!      do k=1,kmax
!        do j=1,jmax
!          do i=1,imax

!!           diffusion

!            dxi2  = dxi**2 + dyi**2 + dzi**2
!            dt1   = 1./( 4.*visc*dxi2 )
!            dt1a  = 1./( 4.*visc2*dxi2 ) !# 2nd phase
!	    dt1   = min(dt1,dt1a)
! 
!!           convection

!            u2 = ( 0.25*(unew(i,j,k)+unew(i,j+1,k)+unew(i-1,j+1,k)+unew(i-1,j,k)) )**2.
!            v2 = ( vnew(i,j,k) )**2.
!            w2 = ( 0.25*(wnew(i,j,k)+wnew(i,j+1,k)+wnew(i,j+1,k-1)+wnew(i,j,k-1)) )**2.
!            velo = sqrt(u2+v2+w2)
!            if (abs(1.*velo) .gt. 1.e-12) then
!              dt3 = 1./( sqrt(3.)*velo*sqrt(dxi2) )
!            else
!              dt3 = 99999.
!            endif

!!           Courant-number weighted with cell-Peclet number
!!           sum (((Pe_i)**(1./3.))*u_i*dt/dx_i) < 2**(1./3.) = 1.26
!!           With Pe_i < 2 (?) --> sum (u_i*dt/dx_i) < 1
!!           This condition is therefore expected to be less
!!           restrictive than condition 2.

!            t1 = ( abs( 0.25*(unew(i,j,k)+unew(i,j+1,k)+unew(i-1,j+1,k)+unew(i-1,j,k))) )**(4./3.)
!            t2 = ( abs(vnew(i,j,k)) )**(4./3.)
!            t3 = ( abs( 0.25*(wnew(i,j,k)+wnew(i,j+1,k)+wnew(i,j+1,k-1)+wnew(i,j,k-1))) )**(4./3.)
!            if (abs(1.*velo) .gt. 1.e-12) then
!              dt2   = ((2.*visc)**(1./3.))/ &
!                      (t1*(dxi**(2./3.)) + t2*(dyi**(2./3.)) + t3*(dyi**(2./3.)))
!              dt2a  = ((2.*visc2)**(1./3.))/ &
!                      (t1*(dxi**(2./3.)) + t2*(dyi**(2./3.)) + t3*(dyi**(2./3.))) !# 2nd phase
!	      dt2   = min(dt2,dt2a)
!            else
!              dt2 = 99999.
!            endif
!            if (dt1 .lt. dt2) then
!              if (dt1 .lt. dt3) then
!                flag = 1
!                counter1 = counter1 + 1
!              else
!                flag = 3
!                counter3 = counter3 + 1
!              endif
!            else
!              if (dt2 .lt. dt3) then
!                flag = 2
!                counter2 = counter2 + 1
!              else
!                flag = 3
!                counter3 = counter3 + 1
!              endif
!            endif
!            dummy       = dtmaxv(1,1)
!            dtmaxv(1,1) = min(dummy,dt1,dt2,dt3)
!            if (dtmaxv(1,1) .lt. dtmaxold) then
!              maxi     = i
!              maxj     = j
!              maxk     = k
!              flagmax  = flag
!              dtmaxold = dtmaxv(1,1)
!            endif
!          enddo
!        enddo
!      enddo
!      call mpi_allreduce(dtmaxv,dtmax_all,1,mpi_2double_precision,mpi_minloc,comm_cart,error)
!      call mpi_allreduce(counter1,counter1_all,1,mpi_integer,mpi_sum,comm_cart,error)
!      call mpi_allreduce(counter2,counter2_all,1,mpi_integer,mpi_sum,comm_cart,error)
!      call mpi_allreduce(counter3,counter3_all,1,mpi_integer,mpi_sum,comm_cart,error)
!      dtmaxv(1,1) = dtmax_all(1,1)
!      dtmaxv(2,1) = dtmax_all(2,1)
!      maxiv       = maxi
!      maxjv       = maxj
!      maxkv       = maxk
!      counter1v   = counter1_all
!      counter2v   = counter2_all
!      counter3v   = counter3_all
!      flagmaxv    = flagmax
!!
!!     wall-normal velocity
!!
!      dtmaxw(1,1) = 899999.
!      dtmaxw(2,1) = myid*1.
!      dtmaxold    = dtmaxw(1,1)
!      counter1    = 0
!      counter2    = 0
!      counter3    = 0
!      maxi        = 1000
!      maxj        = 1000
!      maxk        = 1000
!      flag        = 0
!      flagmax     = 0
!      do k=1,kmax
!        do j=1,jmax
!          do i=1,imax

!!           diffusion

!            dxi2  = dxi**2 + dyi**2 + dzi**2
!            dt1   = 1./( 4.*visc*dxi2 )
!            dt1a  = 1./( 4.*visc2*dxi2 ) !# 2nd phase
!	    dt1   = min(dt1,dt1a)

!!           convection

!            u2 = ( 0.25*(unew(i,j,k)+unew(i-1,j,k)+unew(i-1,j,k+1)+unew(i,j,k+1)) )**2.
!            v2 = ( 0.25*(vnew(i,j,k)+vnew(i,j-1,k)+vnew(i,j-1,k+1)+vnew(i,j,k+1)) )**2.
!            w2 = ( wnew(i,j,k) )**2.
!            velo = sqrt(u2+v2+w2)
!            if (abs(1.*velo) .gt. 1.e-12) then
!              dt3 = 1./( sqrt(3.)*velo*sqrt(dxi2) )
!            else
!              dt3 = 99999.
!            endif

!!           Courant-number weighted with cell-Peclet number
!!           sum (((Pe_i)**(1./3.))*u_i*dt/dx_i) < 2**(1./3.) = 1.26
!!           With Pe_i < 2 (?) --> sum (u_i*dt/dx_i) < 1
!!           This condition is therefore expected to be less
!!           restrictive than condition 2.

!            t1 = ( abs( 0.25*(unew(i,j,k)+unew(i-1,j,k)+unew(i-1,j,k+1)+unew(i,j,k+1))) )**(4./3.)
!            t2 = ( abs( 0.25*(vnew(i,j,k)+vnew(i,j-1,k)+vnew(i,j-1,k+1)+vnew(i,j,k+1))) )**(4./3.)
!            t3 = ( abs(wnew(i,j,k)) )**(4./3.)
!            if (abs(1.*velo) .gt. 1.e-12) then
!              dt2 = ((2.*visc)**(1./3.))/ &
!                    (t1*(dxi**(2./3.)) + t2*(dyi**(2./3.)) + t3*(dyi**(2./3.)))
!              dt2a= ((2.*visc2)**(1./3.))/ &
!                    (t1*(dxi**(2./3.)) + t2*(dyi**(2./3.)) + t3*(dyi**(2./3.))) !# 2nd phase
!	      dt2   = min(dt2,dt2a)
!            else
!              dt2 = 99999.
!            endif
!            if (dt1 .lt. dt2) then
!              if (dt1 .lt. dt3) then
!                flag = 1
!                counter1 = counter1 + 1
!              else
!                flag = 3
!                counter3 = counter3 + 1
!              endif
!            else
!              if (dt2 .lt. dt3) then
!                flag = 2
!                counter2 = counter2 + 1
!              else
!                flag = 3
!                counter3 = counter3 + 1
!              endif
!            endif
!            dummy       = dtmaxw(1,1)
!            dtmaxw(1,1) = min(dummy,dt1,dt2,dt3)
!            if (dtmaxw(1,1) .lt. dtmaxold) then
!              maxi     = i
!              maxj     = j
!              maxk     = k
!              flagmax  = flag
!              dtmaxold = dtmaxw(1,1)
!            endif
!          enddo
!        enddo
!      enddo
!      call mpi_allreduce(dtmaxw,dtmax_all,1,mpi_2double_precision,mpi_minloc,comm_cart,error)
!      call mpi_allreduce(counter1,counter1_all,1,mpi_integer,mpi_sum,comm_cart,error)
!      call mpi_allreduce(counter2,counter2_all,1,mpi_integer,mpi_sum,comm_cart,error)
!      call mpi_allreduce(counter3,counter3_all,1,mpi_integer,mpi_sum,comm_cart,error)
!      dtmaxw(1,1) = dtmax_all(1,1)
!      dtmaxw(2,1) = dtmax_all(2,1)
!      maxiw       = maxi
!      maxjw       = maxj
!      maxkw       = maxk
!      counter1w   = counter1_all
!      counter2w   = counter2_all
!      counter3w   = counter3_all
!      flagmaxw    = flagmax

!      if (dtmaxu(1,1) .lt. dtmaxv(1,1)) then
!        if (dtmaxu(1,1) .lt. dtmaxw(1,1)) then
!          dtmax = dtmaxu(1,1)
!!          if (myid .eq. nint(dtmaxu(2,1))) then
!!            write(6,*) 'STABILITY FOR spanwise MOMENTUM EQUATION'
!!            write(6,*) 'Maximum allowed delta t, i, j, k = ', dtmax,maxiu+coords(1)*itot/dims(1), &
!!                        maxju+coords(2)*jtot/dims(2),maxku
!!            write(6,*) '%1,%2,%3,flag = ',    &
!!       100.*real(counter1u)/(itot*jtot*kmax), &
!!       100.*real(counter2u)/(itot*jtot*kmax), &
!!       100.*real(counter3u)/(itot*jtot*kmax), &
!!       flagmaxu
!!          endif
!        else
!          dtmax = dtmaxw(1,1)
!!          if (myid .eq. nint(dtmaxw(2,1))) then
!!            write(6,*) 'STABILITY FOR WALL-NORMAL MOMENTUM EQUATION'
!!            write(6,*) 'Maximum allowed delta t, i, j, k = ', dtmax,maxiw+coords(1)*itot/dims(1), &
!!                        maxjw+coords(2)*jtot/dims(2),maxkw
!!            write(6,*) '%1,%2,%3,flag = ',    &
!!       100.*real(counter1w)/(itot*jtot*kmax), &
!!       100.*real(counter2w)/(itot*jtot*kmax), &
!!       100.*real(counter3w)/(itot*jtot*kmax), &
!!       flagmaxw
!!          endif
!        endif
!      else
!        dtmax = dtmaxv(1,1)
!!        if (myid .eq. nint(dtmaxv(2,1))) then
!!          write(6,*) 'STABILITY FOR streamwise MOMENTUM EQUATION'
!!          write(6,*) 'Maximum allowed delta t, i, j, k = ',dtmax,maxiv+coords(1)*itot/dims(1), &
!!                      maxjv+coords(2)*jtot/dims(2),maxkv
!!          write(6,*) '%1,%2,%3,flag = ',      &
!!       100.*real(counter1v)/(itot*jtot*kmax), &
!!       100.*real(counter2v)/(itot*jtot*kmax), &
!!       100.*real(counter3v)/(itot*jtot*kmax), &
!!       flagmaxv
!!        endif
!      endif

!!# Commented above is info of locations that requires the smallest dt.

!      !# time step restriction due to surface tension (Dodd & Ferrante JCP2014)

!      dtsurf = sqrt(We*(rho1+rho2)*dx**3 /4.*pi)
!      if (dtsurf .lt. dtmax) dtmax = dtsurf






!! calculates the timestep dt based on stability conditions
!! for RK3 scheme (see: P. Wesseling, 'Principles of computational
!! fluid dynamics').
!!
!real, intent(inout) :: dtmax
!real :: velo,dxi2
!real :: u2,v2,w2
!real :: dt1,dt2
!real :: dtmaxold,dtmax_all(2,1)
!real :: dtmaxu(2,1)
!real :: dtmaxv(2,1)
!real :: dtmaxw(2,1)
!real  :: flagmaxu,flagmaxv,flagmaxw
!real  :: dummy
!integer :: counter1,counter2
!integer :: counter1_all,counter2_all
!integer :: counter1u,counter2u
!integer :: counter1v,counter2v
!integer :: counter1w,counter2w
!integer :: maxi,maxj,maxk
!integer :: maxiu,maxju,maxku
!integer :: maxiv,maxjv,maxkv
!integer :: maxiw,maxjw,maxkw
!integer :: flag,flagmax
!integer :: i,j,k
!!
!! streamwise velocity
!!
!dtmaxu(1,1) = 899999.
!dtmaxu(2,1) = myid*1.
!dtmaxold    = dtmaxu(1,1)
!counter1    = 0
!counter2    = 0
!maxi        = 1000
!maxj        = 1000
!maxk        = 1000
!flag        = 0
!flagmax     = 0
!do k=1,kmax
!  do j=1,jmax
!    do i=1,imax
!      ! diffusion
!      dxi2  = dxi**2 + dyi**2 + dzi**2
!      dt1   = 1.65/( 4.*visc*dxi2 )  ! 2.5127/( 4.*visc*dxi2 )
!      ! convection
!      u2 = ( unew(i,j,k) )**2.
!      v2 = ( 0.25*(vnew(i,j,k)+vnew(i,j-1,k)+vnew(i+1,j,k)+vnew(i+1,j-1,k)) )**2.
!      w2 = ( 0.25*(wnew(i,j,k)+wnew(i,j,k-1)+wnew(i+1,j,k)+wnew(i+1,j,k-1)) )**2.
!      velo = sqrt(u2+v2+w2)
!      if (abs(1.*velo) .gt. 1.e-12) then
!         dt2   = ( sqrt(1.5) )/( velo*sqrt(dxi2) )
!        dt2   = sqrt(3.)/( sqrt(u2)*dxi + sqrt(v2)*dyi + sqrt(w2)*dzi )
!      else
!        dt2 = 99999.
!      endif
!      if (dt1 .lt. dt2) then
!        flag = 1
!        counter1 = counter1 + 1
!      else
!        flag = 2
!        counter2 = counter2 + 1
!      endif
!      dummy       = dtmaxu(1,1)
!      dtmaxu(1,1) = min(dummy,dt1,dt2)
!      if (dtmaxu(1,1) .lt. dtmaxold) then
!        maxi     = i
!        maxj     = j
!        maxk     = k
!        flagmax  = flag
!        dtmaxold = dtmaxu(1,1)
!      endif
!    enddo
!  enddo
!enddo
!call mpi_allreduce(dtmaxu,dtmax_all,1,mpi_2double_precision,mpi_minloc,comm_cart,error)
!call mpi_allreduce(counter1,counter1_all,1,mpi_integer,mpi_sum,comm_cart,error)
!call mpi_allreduce(counter2,counter2_all,1,mpi_integer,mpi_sum,comm_cart,error)
!dtmaxu(1,1) = dtmax_all(1,1)
!dtmaxu(2,1) = dtmax_all(2,1)
!maxiu       = maxi
!maxju       = maxj
!maxku       = maxk
!counter1u   = counter1_all
!counter2u   = counter2_all
!flagmaxu    = flagmax
!!
!! spanwise velocity
!!
!dtmaxv(1,1) = 899999.
!dtmaxv(2,1) = myid*1.
!dtmaxold    = dtmaxv(1,1)
!counter1    = 0
!counter2    = 0
!maxi        = 1000
!maxj        = 1000
!maxk        = 1000
!flag        = 0
!flagmax     = 0
!do k=1,kmax
!  do j=1,jmax
!    do i=1,imax
!      ! diffusion
!      dxi2  = dxi**2 + dyi**2 + dzi**2
!      dt1   = 1.65/( 4.*visc*dxi2 ) ! 2.5127/( 4.*visc*dxi2 )
!      ! convection
!      u2 = ( 0.25*(unew(i,j,k)+unew(i,j+1,k)+unew(i-1,j+1,k)+unew(i-1,j,k)) )**2.
!      v2 = ( vnew(i,j,k) )**2.
!      w2 = ( 0.25*(wnew(i,j,k)+wnew(i,j+1,k)+wnew(i,j+1,k-1)+wnew(i,j,k-1)) )**2.
!      velo = sqrt(u2+v2+w2)
!      if (abs(1.*velo) .gt. 1.e-12) then
!         dt2   = ( sqrt(1.5) )/( velo*sqrt(dxi2) )
!        dt2   = sqrt(3.)/( sqrt(u2)*dxi + sqrt(v2)*dyi + sqrt(w2)*dzi )
!      else
!        dt2 = 99999.
!      endif
!      if (dt1 .lt. dt2) then
!        flag = 1
!        counter1 = counter1 + 1
!      else
!        flag = 2
!        counter2 = counter2 + 1
!      endif
!      dummy       = dtmaxv(1,1)
!      dtmaxv(1,1) = min(dummy,dt1,dt2)
!      if (dtmaxv(1,1) .lt. dtmaxold) then
!        maxi     = i
!        maxj     = j
!        maxk     = k
!        flagmax  = flag
!        dtmaxold = dtmaxv(1,1)
!      endif
!    enddo
!  enddo
!enddo
!call mpi_allreduce(dtmaxv,dtmax_all,1,mpi_2double_precision,mpi_minloc,comm_cart,error)
!call mpi_allreduce(counter1,counter1_all,1,mpi_integer,mpi_sum,comm_cart,error)
!call mpi_allreduce(counter2,counter2_all,1,mpi_integer,mpi_sum,comm_cart,error)
!dtmaxv(1,1) = dtmax_all(1,1)
!dtmaxv(2,1) = dtmax_all(2,1)
!maxiv       = maxi
!maxjv       = maxj
!maxkv       = maxk
!counter1v   = counter1_all
!counter2v   = counter2_all
!flagmaxv    = flagmax
!!
!! wall-normal velocity
!!
!dtmaxw(1,1) = 899999.
!dtmaxw(2,1) = myid*1.
!dtmaxold    = dtmaxw(1,1)
!counter1    = 0
!counter2    = 0
!maxi        = 1000
!maxj        = 1000
!maxk        = 1000
!flag        = 0
!flagmax     = 0
!do k=1,kmax
!  do j=1,jmax
!    do i=1,imax
!      ! diffusion
!      dxi2  = dxi**2 + dyi**2 + dzi**2
!      dt1   = 1.65/( 4.*visc*dxi2 ) ! 2.5127/( 4.*visc*dxi2 )
!      ! convection
!      u2 = ( 0.25*(unew(i,j,k)+unew(i-1,j,k)+unew(i-1,j,k+1)+unew(i,j,k+1)) )**2.
!      v2 = ( 0.25*(vnew(i,j,k)+vnew(i,j-1,k)+vnew(i,j-1,k+1)+vnew(i,j,k+1)) )**2.
!      w2 = ( wnew(i,j,k) )**2.
!      velo = sqrt(u2+v2+w2)
!      if (abs(1.*velo) .gt. 1.e-12) then
!         dt2   = ( sqrt(1.5) )/( velo*sqrt(dxi2) )
!        dt2   = sqrt(3.)/( sqrt(u2)*dxi + sqrt(v2)*dyi + sqrt(w2)*dzi )
!      else
!        dt2 = 99999.
!      endif
!      if (dt1 .lt. dt2) then
!        flag = 1
!        counter1 = counter1 + 1
!      else
!        flag = 2
!        counter2 = counter2 + 1
!      endif
!      dummy       = dtmaxw(1,1)
!      dtmaxw(1,1) = min(dummy,dt1,dt2)
!      if (dtmaxw(1,1) .lt. dtmaxold) then
!        maxi     = i
!        maxj     = j
!        maxk     = k
!        flagmax  = flag
!        dtmaxold = dtmaxw(1,1)
!      endif
!    enddo
!  enddo
!enddo
!call mpi_allreduce(dtmaxw,dtmax_all,1,mpi_2double_precision,mpi_minloc,comm_cart,error)
!call mpi_allreduce(counter1,counter1_all,1,mpi_integer,mpi_sum,comm_cart,error)
!call mpi_allreduce(counter2,counter2_all,1,mpi_integer,mpi_sum,comm_cart,error)
!dtmaxw(1,1) = dtmax_all(1,1)
!dtmaxw(2,1) = dtmax_all(2,1)
!maxiw       = maxi
!maxjw       = maxj
!maxkw       = maxk
!counter1w   = counter1_all
!counter2w   = counter2_all
!flagmaxw    = flagmax
!!
!! determine maximum time step
!!
!dtmax = min(dtmaxu(1,1), dtmaxv(1,1), dtmaxw(1,1),Reb*dx*dx/6.)  !# I add the last term


!!
!! x-momentum
!!
!if ( dtmaxu(1,1) .eq. dtmax ) then
!  if (myid .eq. nint(dtmaxu(2,1))) then
!    write(6,*) 'STABILITY FOR STREAMWISE MOMENTUM EQUATION'
!    write(6,*) 'Maximum allowed delta t, i, j, k = ', dtmaxu(1,1), maxiu+coords(1)*itot/dims(1), &
!                maxju+coords(2)*jtot/dims(2), maxku 
!    write(6,*) '%1, %2, flag = ', &
!                100*real(counter1u)/(itot*jtot*kmax), &
!                100*real(counter2u)/(itot*jtot*kmax),flagmaxu
!  endif
!endif
!!
!! y-momentum
!!
!if ( dtmaxv(1,1) .eq. dtmax ) then
!  if (myid .eq. nint(dtmaxv(2,1))) then
!    write(6,*) 'STABILITY FOR SPANWISE MOMENTUM EQUATION'
!    write(6,*) 'Maximum allowed delta t, i, j, k = ', dtmaxv(1,1), maxiv+coords(1)*itot/dims(1), &
!             maxjv+coords(2)*jtot/dims(2), maxkv
!    write(6,*) '%1, %2, flag = ', &
!              100*real(counter1v)/(itot*jtot*kmax), &
!              100*real(counter2v)/(itot*jtot*kmax),flagmaxv
!  endif
!endif
!!
!! z-momentum
!!
!if ( dtmaxw(1,1) .eq. dtmax ) then
!  if (myid .eq. nint(dtmaxw(2,1))) then
!    write(6,*) 'STABILITY FOR WALL-NORMAL MOMENTUM EQUATION'
!    write(6,*) 'Maximum allowed delta t, i, j, k = ', dtmaxw(1,1), maxiw+coords(1)*itot/dims(1), &
!               maxjw+coords(2)*jtot/dims(2), maxkw
!    write(6,*) '%1, %2, flag = ', &
!                100*real(counter1w)/(itot*jtot*kmax), &
!                100*real(counter2w)/(itot*jtot*kmax),flagmaxw
!  endif
!endif
!!


    return
  end subroutine chkdt


end module mod_chkdt
