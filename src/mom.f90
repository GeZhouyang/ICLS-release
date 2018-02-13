module mod_mom

 use mod_param

 implicit none
 
 private
 public momxad,momyad,momzad,momsad, momxp,momyp,momzp, mom_main
 
 contains
   
   ! #0   
   subroutine mom_main(RHSu,RHSv,RHSw)

     use mod_common
     use mod_gfm

     real, dimension(1:imax,1:jmax,1:kmax) :: RCu,RCv,RCw
     real, dimension(1:imax,1:jmax,1:kmax) :: RDu,RDv,RDw
     real, dimension(1:imax,1:jmax,1:kmax) :: RGu,RGv,RGw
     
     real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: RHSu,RHSv,RHSw

     
     !-- Convection --!

     select case (NS_conv)
     case('CNTR2')
        call CNTR2_conv(RCu,RCv,RCw, unew,vnew,wnew)
        
     case('WENO5')
        call WENO5_conv('x-dir',RCu)
        call WENO5_conv('y-dir',RCv)
        call WENO5_conv('z-dir',RCw)
     end select
     
     
     !-- Diffusion --!
     
     select case (NS_diff_discretization)
     case('del')
        call delta_diff(RDu,RDv,RDw, unew,vnew,wnew, miu, rho_u,rho_v,rho_w)

     case('GFM')
        call gfm4diff(RDu,RDv,RDw, unew,vnew,wnew, rho_u,rho_v,rho_w)
     end select


     !-- Gravity --!

     call NS_grav(RGu,RGv,RGw)

     
     !-- Output --!
     
     RHSu = RCu + RDu*visc + RGu
     RHSv = RCv + RDv*visc + RGv
     RHSw = RCw + RDw*visc + RGw


     return
   end subroutine mom_main


   
   ! #1
   subroutine CNTR2_conv(RCu,RCv,RCw, u,v,w)

     integer :: i,j,k, ip,jp,kp, im,jm,km

     real :: uuip,uuim, uvjp,uvjm, uwkp,uwkm
     real :: uvip,uvim, vvjp,vvjm, wvkp,wvkm
     real :: uwip,uwim, vwjp,vwjm, wwkp,wwkm

     real, dimension(0:,0:,0:), intent(in) :: u,v,w
     real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: RCu,RCv,RCw

     
     do k = 1,kmax
        do j = 1,jmax
           do i = 1,imax

              ip = i + 1
              jp = j + 1
              kp = k + 1
              im = i - 1
              jm = j - 1
              km = k - 1

              !-- x direction

              uuip  = 0.25 *( u(ip,j ,k ) +u(i,j,k) )*( u(ip,j ,k ) +u(i,j ,k ) )
              uuim  = 0.25 *( u(im,j ,k ) +u(i,j,k) )*( u(im,j ,k ) +u(i,j ,k ) )
              uvjp  = 0.25 *( u(i ,jp,k ) +u(i,j,k) )*( v(ip,j ,k ) +v(i,j ,k ) )
              uvjm  = 0.25 *( u(i ,jm,k ) +u(i,j,k) )*( v(ip,jm,k ) +v(i,jm,k ) )
              uwkp  = 0.25 *( u(i ,j ,kp) +u(i,j,k) )*( w(ip,j ,k ) +w(i,j ,k ) )
              uwkm  = 0.25 *( u(i ,j ,km) +u(i,j,k) )*( w(ip,j ,km) +w(i,j ,km) )

              RCu(i,j,k) = (-uuip +uuim)/dx + (-uvjp +uvjm)/dy + (-uwkp +uwkm)/dz

              !-- y direction
              
              uvip  = 0.25 *(u(i ,j,k ) +u(i ,jp,k ) ) *(v(i,j,k ) +v(ip,j ,k) )
              uvim  = 0.25 *(u(im,j,k ) +u(im,jp,k ) ) *(v(i,j,k ) +v(im,j ,k) )
              vvjp  = 0.25 *(v(i ,j,k ) +v(i ,jp,k ) ) *(v(i,j,k ) +v(i ,jp,k) )
              vvjm  = 0.25 *(v(i ,j,k ) +v(i ,jm,k ) ) *(v(i,j,k ) +v(i ,jm,k) )
              wvkp  = 0.25 *(w(i ,j,k ) +w(i ,jp,k ) ) *(v(i,j,kp) +v(i ,j ,k) )
              wvkm  = 0.25 *(w(i ,j,km) +w(i ,jp,km) ) *(v(i,j,km) +v(i ,j ,k) )

              RCv(i,j,k) = (-uvip +uvim)/dx + (-vvjp +vvjm)/dy + (-wvkp +wvkm)/dz

              !-- z direction

              uwip  = 0.25 *(w(i,j,k) +w(ip,j ,k ) ) *(u(i ,j ,k) +u(i ,j ,kp) )
              uwim  = 0.25 *(w(i,j,k) +w(im,j ,k ) ) *(u(im,j ,k) +u(im,j ,kp) )
              vwjp  = 0.25 *(w(i,j,k) +w(i ,jp,k ) ) *(v(i ,j ,k) +v(i ,j ,kp) )
              vwjm  = 0.25 *(w(i,j,k) +w(i ,jm,k ) ) *(v(i ,jm,k) +v(i ,jm,kp) )
              wwkp  = 0.25 *(w(i,j,k) +w(i ,j ,kp) ) *(w(i ,j ,k) +w(i ,j ,kp ))
              wwkm  = 0.25 *(w(i,j,k) +w(i ,j ,km) ) *(w(i ,j ,k) +w(i ,j ,km ))

              RCw(i,j,k) = (-uwip +uwim)/dx + (-vwjp +vwjm)/dy + (-wwkp +wwkm)/dz

           enddo
        enddo
     enddo

     if (TwoD) RCu = 0.

     return
   end subroutine CNTR2_conv
   
   
   ! #2
   subroutine delta_diff(RDu,RDv,RDw, u,v,w, miu, rho_u,rho_v,rho_w)

     real, dimension(0:,0:,0:), intent(in) :: u,v,w, miu
     real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: rho_u,rho_v,rho_w
     real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: RDu,RDv,RDw
     
     integer :: i,j,k, im,ip, jm,jp, km,kp
     
     real :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
     real :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
     real :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
     real :: ma1,ma2,ma3,ma4,ma5,ma6,ma7,ma8,ma9
     real :: rhou,rhov,rhow
     real :: RD1,RD2,RD3,RD4,RD5
     

     do k = 1,kmax
        do j = 1,jmax
           do i = 1,imax

              ip = i + 1
              jp = j + 1
              kp = k + 1
              im = i - 1
              jm = j - 1
              km = k - 1

              !-- u-component
              !
              ! first-derivatives at u-grid

              dudxp = (u(ip,j,k)-u(i ,j,k))/dx
              dudxm = (u(i ,j,k)-u(im,j,k))/dx

              dudyp = (u(i,jp,k)-u(i,j ,k))/dy
              dudym = (u(i,j ,k)-u(i,jm,k))/dy

              dvdxp = (v(ip,j ,k)-v(i,j ,k))/dx     
              dvdxm = (v(ip,jm,k)-v(i,jm,k))/dx  

              dudzp = (u(i,j,kp)-u(i,j,k ))/dz
              dudzm = (u(i,j,k )-u(i,j,km))/dz

              dwdxp = (w(ip,j,k )-w(i,j,k ))/dx   
              dwdxm = (w(ip,j,km)-w(i,j,km))/dx

              ! averaging at face centers

              ma1  = 0.25 *( miu(i,j,k)+miu(i,jm,k )+miu(ip,jm,k )+miu(ip,j,k) )   ! miu(i+0.5,j-0.5,k)
              ma2  = 0.25 *( miu(i,j,k)+miu(i,jp,k )+miu(ip,jp,k )+miu(ip,j,k) )   ! miu(i+0.5,j+0.5,k)
              ma3  = 0.25 *( miu(i,j,k)+miu(i,j ,km)+miu(ip,j ,km)+miu(ip,j,k) )   ! miu(i+0.5,j,k-0.5)
              ma4  = 0.25 *( miu(i,j,k)+miu(i,j ,kp)+miu(ip, j,kp)+miu(ip,j,k) )   ! miu(i+0.5,j,k+0.5)

              RD1 = 2.*(dudxp*miu(ip,j,k)-dudxm*miu(i,j,k))/dx  ! d/dx(2*miu*du/dx)
              RD2 = (dudyp*ma2-dudym*ma1)/dy                    ! d/dy(miu*du/dy)
              RD3 = (dvdxp*ma2-dvdxm*ma1)/dy                    ! d/dy(miu*dv/dx)
              RD4 = (dudzp*ma4-dudzm*ma3)/dz                    ! d/dz(miu*du/dz)
              RD5 = (dwdxp*ma4-dwdxm*ma3)/dz                    ! d/dz(miu*dw/dx)

              RDu(i,j,k) = (RD1+RD2+RD3+RD4+RD5)/rho_u(i,j,k)
              

              !-- v-component
              !
              ! first-derivatives at v-grid

              dvdxp = (v(ip,j,k)-v(i ,j,k))/dx
              dvdxm = (v(i ,j,k)-v(im,j,k))/dx

              dudyp = (u(i ,jp,k)-u(i ,j,k))/dy     
              dudym = (u(im,jp,k)-u(im,j,k))/dy  

              dvdyp = (v(i,jp,k)-v(i,j ,k))/dy
              dvdym = (v(i,j ,k)-v(i,jm,k))/dy

              dvdzp = (v(i,j,kp)-v(i,j,k ))/dz
              dvdzm = (v(i,j,k )-v(i,j,km))/dz

              dwdyp = (w(i,jp,k )-w(i,j,k ))/dy    
              dwdym = (w(i,jp,km)-w(i,j,km))/dy  

              ! averaging at face centers

              ma2  = 0.25 *( miu(i,j,k)+miu(i,jp,k )+miu(ip,jp,k )+miu(ip,j ,k) )   ! miu(i+0.5,j+0.5,k)
              ma5  = 0.25 *( miu(i,j,k)+miu(i,jp,k )+miu(im,jp,k )+miu(im,j ,k) )   ! miu(i-0.5,j+0.5,k)
              ma6  = 0.25 *( miu(i,j,k)+miu(i,j ,km)+miu(i ,jp,km)+miu(i ,jp,k) )   ! miu(i,j+0.5,k-0.5)
              ma7  = 0.25 *( miu(i,j,k)+miu(i,j ,kp)+miu(i ,jp,kp)+miu(i ,jp,k) )   ! miu(i,j+0.5,k+0.5)      

              RD1 = (dvdxp*ma2-dvdxm*ma5)/dx                       ! d/dx(miu*dv/dx)
              RD2 = (dudyp*ma2-dudym*ma5)/dx                       ! d/dx(miu*du/dy)
              RD3 = 2.*(dvdyp*miu(i,jp,k)-dvdym*miu(i,j,k))/dy     ! d/dy(2*miu*dv/dy)
              RD4 = (dvdzp*ma7-dvdzm*ma6)/dz                       ! d/dz(miu*dv/dz)
              RD5 = (dwdyp*ma7-dwdym*ma6)/dz                       ! d/dz(miu*dw/dy)

              RDv(i,j,k) = (RD1+RD2+RD3+RD4+RD5)/rho_v(i,j,k)

              
              !-- w-component
              !
              ! first-derivatives at w-grid

              dwdxp = (w(ip,j,k)-w(i ,j,k))/dx
              dwdxm = (w(i ,j,k)-w(im,j,k))/dx

              dudzp = (u(i ,j,kp)-u(i ,j,k))/dz   
              dudzm = (u(im,j,kp)-u(im,j,k))/dz  

              dwdyp = (w(i,jp,k)-w(i,j ,k))/dy
              dwdym = (w(i,j ,k)-w(i,jm,k))/dy

              dvdzp = (v(i,j ,kp)-v(i,j ,k))/dz     
              dvdzm = (v(i,jm,kp)-v(i,jm,k))/dz  

              dwdzp = (w(i,j,kp)-w(i,j,k ))/dz
              dwdzm = (w(i,j,k )-w(i,j,km))/dz

              ! averaging at cell faces

              ma4  = 0.25 *(miu(i,j,k)+miu(i,j,kp)+miu(ip,j ,kp)+miu(ip,j ,k) )   ! miu(i+0.5,j,k+0.5)
              ma8  = 0.25 *(miu(i,j,k)+miu(i,j,kp)+miu(im,j ,kp)+miu(im,j ,k) )   ! miu(i-0.5,j,k+0.5)
              ma7  = 0.25 *(miu(i,j,k)+miu(i,j,kp)+miu(i ,jp,kp)+miu(i ,jp,k) )   ! miu(i,j+0.5,k+0.5)
              ma9  = 0.25 *(miu(i,j,k)+miu(i,j,kp)+miu(i ,jm,kp)+miu(i ,jm,k) )   ! miu(i,j-0.5,k+0.5)

              RD1 = (dwdxp*ma4-dwdxm*ma8)/dx                    ! d/dx(miu*dw/dx)
              RD2 = (dudzp*ma4-dudzm*ma8)/dx                    ! d/dx(miu*du/dx)
              RD3 = (dwdyp*ma7-dwdym*ma9)/dy                    ! d/dy(miu*dw/dy)
              RD4 = (dvdzp*ma7-dvdzm*ma9)/dy                    ! d/dy(miu*dv/dy)
              RD5 = 2.*(dwdzp*miu(i,j,kp)-dwdzm*miu(i,j,k))/dz  ! d/dz(2*miu*dw/dz)

              RDw(i,j,k) = (RD1+RD2+RD3+RD4+RD5)/rho_w(i,j,k)

           enddo
        enddo
     enddo

     if (TwoD) RDu = 0.
     

     return
   end subroutine delta_diff


   ! #3
   subroutine NS_grav(RGu,RGv,RGw)

     use mod_common

     integer :: i,j,k

     real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: RGu,RGv,RGw

     do k = 1,kmax
        do j = 1,jmax
           do i = 1,imax     

              RGu(i,j,k) = iFr*(1.-rho_avr/rho_u(i,j,k))*g_unit(1)
              RGv(i,j,k) = iFr*(1.-rho_avr/rho_v(i,j,k))*g_unit(2)
              RGw(i,j,k) = iFr*(1.-rho_avr/rho_w(i,j,k))*g_unit(3)

           enddo
        enddo
     enddo

     if (TwoD) RGu = 0.

     
     return
   end subroutine NS_grav

   

   ! #1+
   subroutine WENO5_conv(dir,Lconv)
   !                                              !
   ! Compute the NS convective terms using WENO5. !      
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   
   use mod_common
   use mod_bound
   
   character(len=5), intent(in) :: dir
   
   real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: Lconv
   real, dimension(1:imax,1:jmax,1:kmax) :: wind_u, wind_v, wind_w
   real, dimension(1:imax,1:jmax,1:kmax) :: convx, convy, convz
   
   !-- Buffer flux matrix
   
   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: fx,fy,fz, fdir
   
   !-- Parameters
   
   real :: eps_weno
   parameter ( eps_weno = 1.e-6 )
   real :: sigma_1, sigma_2, sigma_3
   parameter ( sigma_1 = 0.1 )
   parameter ( sigma_2 = 0.6 )
   parameter ( sigma_3 = 0.3 )
   integer :: i, j, k
   real :: v1, v2, v3, v4, v5
   real :: S1, S2, S3
   real :: a1, a2, a3
   real :: w1, w2, w3
   real :: vel_vertice, f_x,f_y,f_z
   
   
    !-- Initialize fluxes --!
   
    fx = 0.
    fy = 0.
    fz = 0.
   
    fx(0:i1,0:j1,0:k1) = unew
    fy(0:i1,0:j1,0:k1) = vnew
    fz(0:i1,0:j1,0:k1) = wnew
   
    call boundfuv(fx)
    call boundfuv(fy)
    call boundfw(fz)
   
   
    select case (dir)
   
      case('x-dir')
   
        if (TwoD) then
           Lconv = 0.
           return
        endif
        
        !-- Direction flux
   
        fdir = fx
   
        !-- Face-velocity
   
        do k = 1,kmax
          do j = 1,jmax
            do i = 1,imax
   
              !-- U at unew(i,j,k)
   
              wind_u(i,j,k) = fx(i,j,k)
              
              !-- V at unew(i,j,k)
   
              vel_vertice = (fy(i,j  ,k) + fy(i+1,j  ,k) + &
                             fy(i,j-1,k) + fy(i+1,j-1,k) )/4.
   
              wind_v(i,j,k) = vel_vertice 
   
              !-- W at unew(i,j,k)
   
              if ( (k .ge. 3) .and. (k .le. kmax-2) ) then
   
                vel_vertice = (fz(i,j,k  ) + fz(i+1,j,k  ) + &
                               fz(i,j,k-1) + fz(i+1,j,k-1) )/4.
   
                wind_w(i,j,k) = vel_vertice
   
              else  ! 2nd-order central difference
   
                wind_w(i,j,k) = 0.  ! not used
   
                convz(i,j,k) = ( (fx(i,j,k)+fx(i,j,k+1))/2.*(fz(i,j,k  )+fz(i+1,j,k  ))/2. - &
                                 (fx(i,j,k)+fx(i,j,k-1))/2.*(fz(i,j,k-1)+fz(i+1,j,k-1))/2. )/dz
   
              endif
   
            enddo
          enddo
        enddo
   
      case('y-dir')
   
        !-- Direction flux
   
        fdir = fy
   
        !-- Face-velocity
   
        do k = 1,kmax
          do j = 1,jmax
            do i = 1,imax
   
              !-- U at vnew(i,j,k)
   
              vel_vertice = (fx(i  ,j,k) + fx(i  ,j+1,k) + &
                             fx(i-1,j,k) + fx(i-1,j+1,k) )/4.
   
              wind_u(i,j,k) = vel_vertice
   
              !-- V at vnew(i,j,k)
   
              wind_v(i,j,k) = fy(i,j,k)
   
              !-- W at vnew(i,j,k)
   
              if ( (k .ge. 3) .and. (k .le. kmax-2) ) then
   
                vel_vertice = (fz(i,j,k  ) + fz(i,j+1,k  ) + &
                               fz(i,j,k-1) + fz(i,j+1,k-1) )/4.
   
                wind_w(i,j,k) = vel_vertice
   
              else  ! 2nd-order central difference
   
                wind_w(i,j,k) = 0.  ! not used
   
                convz(i,j,k) = ( (fy(i,j,k)+fy(i,j,k+1))/2.*(fz(i,j,k  )+fz(i,j+1,k  ))/2. - &
                                 (fy(i,j,k)+fy(i,j,k-1))/2.*(fz(i,j,k-1)+fz(i,j+1,k-1))/2. )/dz
              
              endif
   
            enddo
          enddo
        enddo
   
   
      case('z-dir')
   
        !-- Direction flux
   
        fdir = fz
   
        !-- Face-velocity
   
        do k = 1,kmax
          do j = 1,jmax
            do i = 1,imax
   
              !-- U at wnew(i,j,k)
   
              vel_vertice = (fx(i  ,j,k) + fx(i  ,j,k+1) + &
                             fx(i-1,j,k) + fx(i-1,j,k+1) )/4.
   
              wind_u(i,j,k) = vel_vertice
   
              !-- V at wnew(i,j,k)
   
              vel_vertice = (fy(i,j  ,k) + fy(i,j  ,k+1) + &
                             fy(i,j-1,k) + fy(i,j-1,k+1) )/4.
   
              wind_v(i,j,k) = vel_vertice
   
              !-- W at wnew(i,j,k)
   
              if ( (k .ge. 3) .and. (k .le. kmax-2) ) then
   
                wind_w(i,j,k) = fz(i,j,k)
   
              else  ! 2nd-order central difference
   
                wind_w(i,j,k) = 0.  ! not used
   
                convz(i,j,k) = fz(i,j,k) * (fz(i,j,k+1) - fz(i,j,k-1))/2. /dz
   
              endif
   
            enddo
          enddo
        enddo
   
    end select
   
   
    do k = 1,kmax
      do j = 1,jmax
        do i = 1,imax
     
          !-- X-derivative
   
          if (wind_u(i,j,k) .gt. 0.) then
        
            v1 = (fdir(i-2,j,k) - fdir(i-3,j,k))/dx 
            v2 = (fdir(i-1,j,k) - fdir(i-2,j,k))/dx
            v3 = (fdir(i  ,j,k) - fdir(i-1,j,k))/dx
            v4 = (fdir(i+1,j,k) - fdir(i  ,j,k))/dx
            v5 = (fdir(i+2,j,k) - fdir(i+1,j,k))/dx
       
            !-- Smoothness indicators
   
            S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
                 + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
        
            S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                 + (1./4.) * (v2 - v4)**2
        
            S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
                 + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
   
            !-- Nonlinear weights
        
            a1 = sigma_1 / (eps_weno + S1)**2
            a2 = sigma_2 / (eps_weno + S2)**2
            a3 = sigma_3 / (eps_weno + S3)**2
        
            w1 = a1 / (a1 + a2 + a3)
            w2 = a2 / (a1 + a2 + a3)
            w3 = a3 / (a1 + a2 + a3)
       
            !-- Minus flux
   
            f_x = w1 * ( (1./3.)  * v1   &
                       - (7./6.)  * v2   &
                       + (11./6.) * v3 ) &
                + w2 * ( (-1./6.) * v2   &
                       + (5./6.)  * v3   &
                       + (1./3.)  * v4 ) &
                + w3 * ( (1./3.)  * v3   &
                       + (5./6.)  * v4   &
                       - (1./6.)  * v5 )
          else
        
            v1 = (fdir(i+3,j,k) - fdir(i+2,j,k))/dx 
            v2 = (fdir(i+2,j,k) - fdir(i+1,j,k))/dx
            v3 = (fdir(i+1,j,k) - fdir(i  ,j,k))/dx
            v4 = (fdir(i  ,j,k) - fdir(i-1,j,k))/dx
            v5 = (fdir(i-1,j,k) - fdir(i-2,j,k))/dx
               
            !++ Smoothness indicators
   
            S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
                 + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
        
            S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                         + (1./4.) * (v2 - v4)**2
        
            S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
                 + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
   
            !++ Nonlinear weights
   
            a1 = sigma_1 / (eps_weno + S1)**2
            a2 = sigma_2 / (eps_weno + S2)**2
            a3 = sigma_3 / (eps_weno + S3)**2
        
            w1 = a1 / (a1 + a2 + a3)
            w2 = a2 / (a1 + a2 + a3)
            w3 = a3 / (a1 + a2 + a3)
   
            !++ Plus flux
        
            f_x = w1 * ( (1./3.)  * v1   &
                       - (7./6.)  * v2   &
                       + (11./6.) * v3 ) &
                + w2 * ( (-1./6.) * v2   &
                       + (5./6.)  * v3   &
                       + (1./3.)  * v4 ) &
                + w3 * ( (1./3.)  * v3   &
                       + (5./6.)  * v4   &
                       - (1./6.)  * v5 )
   
          endif
   
          ! Time derivative (x-component)
   
          convx(i,j,k) = wind_u(i,j,k)*f_x
   
   
          !-- Y-derivative
   
          if (wind_v(i,j,k) .gt. 0.) then
   
            v1 = (fdir(i,j-2,k) - fdir(i,j-3,k))/dy 
            v2 = (fdir(i,j-1,k) - fdir(i,j-2,k))/dy
            v3 = (fdir(i,j  ,k) - fdir(i,j-1,k))/dy
            v4 = (fdir(i,j+1,k) - fdir(i,j  ,k))/dy
            v5 = (fdir(i,j+2,k) - fdir(i,j+1,k))/dy
   
            !-- Smoothness indicators
   
            S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
                 + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
                      
            S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                           + (1./4.) * (v2 - v4)**2
                      
            S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
                 + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
   
            !-- Nonlinear weights
                
            a1 = sigma_1 / (eps_weno + S1)**2
            a2 = sigma_2 / (eps_weno + S2)**2
            a3 = sigma_3 / (eps_weno + S3)**2
                      
            w1 = a1 / (a1 + a2 + a3)
            w2 = a2 / (a1 + a2 + a3)
            w3 = a3 / (a1 + a2 + a3)
   
            !-- Minus flux
                
            f_y = w1 * ( (1./3.)  * v1   &
                       - (7./6.)  * v2   &
                       + (11./6.) * v3 ) &
                + w2 * ( (-1./6.) * v2   &
                       + (5./6.)  * v3   &
                       + (1./3.)  * v4 ) &
                + w3 * ( (1./3.)  * v3   &
                       + (5./6.)  * v4   &
                       - (1./6.)  * v5 )
                      
          else
              
            v1 = (fdir(i,j+3,k) - fdir(i,j+2,k))/dy 
            v2 = (fdir(i,j+2,k) - fdir(i,j+1,k))/dy 
            v3 = (fdir(i,j+1,k) - fdir(i,j  ,k))/dy
            v4 = (fdir(i,j  ,k) - fdir(i,j-1,k))/dy 
            v5 = (fdir(i,j-1,k) - fdir(i,j-2,k))/dy 
                        
            !++ Smoothness indicators
   
            S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
                 + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
                    
            S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                           + (1./4.) * (v2 - v4)**2
                      
            S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
                 + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
   
            !++ Nonlinear weights
                
            a1 = sigma_1 / (eps_weno + S1)**2
            a2 = sigma_2 / (eps_weno + S2)**2
            a3 = sigma_3 / (eps_weno + S3)**2
                      
            w1 = a1 / (a1 + a2 + a3)
            w2 = a2 / (a1 + a2 + a3)
            w3 = a3 / (a1 + a2 + a3)
   
            !++ Plus flux
    
            f_y = w1 * ( (1./3.)  * v1   &
                       - (7./6.)  * v2   &
                       + (11./6.) * v3 ) &
                + w2 * ( (-1./6.) * v2   &
                       + (5./6.)  * v3   &
                       + (1./3.)  * v4 ) &
                + w3 * ( (1./3.)  * v3   &
                       + (5./6.)  * v4   &
                       - (1./6.)  * v5 )
   
          endif
   
          ! Time derivative (y-component)
           
          convy(i,j,k) = wind_v(i,j,k)*f_y
   
   
          !-- Z-derivative
   
          if ( (k .ge. 3) .and. (k .le. kmax-2) ) then
   
            if (wind_w(i,j,k) .gt. 0.) then
   
              v1 = (fdir(i,j,k-2) - fdir(i,j,k-3))/dz  
              v2 = (fdir(i,j,k-1) - fdir(i,j,k-2))/dz
              v3 = (fdir(i,j,k  ) - fdir(i,j,k-1))/dz
              v4 = (fdir(i,j,k+1) - fdir(i,j,k  ))/dz
              v5 = (fdir(i,j,k+2) - fdir(i,j,k+1))/dz
   
              !-- Smoothness indicators
   
              S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
                   + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
          
              S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                   + (1./4.) * (v2 - v4)**2
          
              S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
                   + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
   
              !-- Nonlinear weights
          
              a1 = sigma_1 / (eps_weno + S1)**2
              a2 = sigma_2 / (eps_weno + S2)**2
              a3 = sigma_3 / (eps_weno + S3)**2
          
              w1 = a1 / (a1 + a2 + a3)
              w2 = a2 / (a1 + a2 + a3)
              w3 = a3 / (a1 + a2 + a3)
   
              !-- Minus flux
          
              f_z = w1 * ( (1./3.)  * v1   &
                         - (7./6.)  * v2   &
                         + (11./6.) * v3 ) &
                  + w2 * ( (-1./6.) * v2   &
                         + (5./6.)  * v3   &
                         + (1./3.)  * v4 ) &
                  + w3 * ( (1./3.)  * v3   &
                         + (5./6.)  * v4   &
                         - (1./6.)  * v5 )
          
            else !++
   
              v1 = (fdir(i,j,k+3) - fdir(i,j,k+2))/dz 
              v2 = (fdir(i,j,k+2) - fdir(i,j,k+1))/dz
              v3 = (fdir(i,j,k+1) - fdir(i,j,k  ))/dz
              v4 = (fdir(i,j,k  ) - fdir(i,j,k-1))/dz
              v5 = (fdir(i,j,k-1) - fdir(i,j,k-2))/dz
   
              !++ Smoothness indicators
   
              S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
                   + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
          
              S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                   + (1./4.) * (v2 - v4)**2
          
              S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
                   + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
   
              !++ Nonlinear weights
          
              a1 = sigma_1 / (eps_weno + S1)**2
              a2 = sigma_2 / (eps_weno + S2)**2
              a3 = sigma_3 / (eps_weno + S3)**2
          
              w1 = a1 / (a1 + a2 + a3)
              w2 = a2 / (a1 + a2 + a3)
              w3 = a3 / (a1 + a2 + a3)
   
              !++ Plus flux
          
              f_z = w1 * ( (1./3.)  * v1   &
                         - (7./6.)  * v2   &
                         + (11./6.) * v3 ) &
                  + w2 * ( (-1./6.) * v2   &
                         + (5./6.)  * v3   &
                         + (1./3.)  * v4 ) &
                  + w3 * ( (1./3.)  * v3   &
                         + (5./6.)  * v4   &
                         - (1./6.)  * v5 )
            endif
   
            ! Time derivative (z-component)
           
            convz(i,j,k) = wind_w(i,j,k)*f_z
   
          endif
    
        enddo
      enddo
    enddo
   
   
    !-- Output convective fluxes --!
   
    Lconv = - convx - convy - convz
   
      
    return
   end subroutine WENO5_conv


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

   
  




  !----------------
  ! 
  !     Subroutines mom?p compute pressure corection terms split from the pressure gradient in ? direction
  !     as used in Dodd & Ferrante (the rest computed with FFT)
  !     
  !     INPUT  : Predicted pressure and density at the next time level
  !     OUTPUT : "Acceleration" at an intermediate time level (AB2)
  ! 
  !     Key internal variables
  !     ----------------------
  !     rho?         : Averaged density at face centers
  !     rho0         : Min density of phase 1 and 2
  ! 
  !-----
  
  subroutine momxp(dudt,p,p_j,rho_u)
  
    real, dimension(0:,0:,0:), intent(in) :: p,p_j
    real, dimension(1:,1:,1:), intent(in) :: rho_u
    real, dimension(1:,1:,1:), intent(out) :: dudt
    integer :: i,j,k
    integer :: ip
  
  
    do k=1,kmax
       do j=1,jmax
          do i=1,imax

             ip = i + 1
             dudt(i,j,k) = - (1./rho_u(i,j,k)-1./rho0)*dxi*( p(ip,j,k)-p_j(i,j,k) -p(i,j,k) )
  
          enddo
       enddo
    enddo
  
  
    return
  end subroutine momxp
  
  
  subroutine momyp(dvdt,p,p_j,rho_v)
  
    real, dimension(0:,0:,0:), intent(in) :: p,p_j
    real, dimension(1:,1:,1:), intent(in) :: rho_v
    real, dimension(1:,1:,1:), intent(out) :: dvdt
    integer :: i,j,k
    integer :: jp
  
  
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             
             jp = j + 1
             dvdt(i,j,k) = - (1./rho_v(i,j,k)-1./rho0)*dyi*( p(i,jp,k)-p_j(i,j,k) -p(i,j,k) )
  
          enddo
       enddo
    enddo
  
  
    return
  end subroutine momyp
  
  
  subroutine momzp(dwdt,p,p_j,rho_w)
  
    real, dimension(0:,0:,0:), intent(in) :: p,p_j
    real, dimension(1:,1:,1:), intent(in) :: rho_w
    real, dimension(1:,1:,1:), intent(out) :: dwdt
    integer :: kp
    integer :: i,j,k
  
  
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             
             kp = k + 1
             dwdt(i,j,k) = - (1./rho_w(i,j,k)-1./rho0)*dzi*( p(i,j,kp)-p_j(i,j,k) -p(i,j,k) )
  
          enddo
       enddo
    enddo
  
  
    return
  end subroutine momzp






!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Not used below

 ! 
 ! Subroutines mom?ad compute convection, diffusion, and gravity terms on the RHS of NS in ? direction
 ! using a second-order accurate central difference scheme.
 ! 
 ! INPUT  : Velocity, viscosity, and density at a known time level
 ! OUTPUT : "Acceleration" at an intermediate time level (AB2) or at the previous time level (restarting)
 ! 
 ! Key internal variables
 ! ----------------------
 ! ma?          : Averaged viscosity at face centers
 ! rho?         : Averaged density at face centers
 ! RC           : Negative convection terms
 ! sum of RD?   : Diffusion terms
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ! #1

   subroutine momxad(dudt, u,v,w,miu,rho,g_unit)   

     real, intent(in) :: g_unit(1)
     real, dimension(0:,0:,0:), intent(in) :: u,v,w,miu,rho  
     real, dimension(0:,0:,0:), intent(out) :: dudt

     integer :: im,ip,jm,jp,km,kp,i,j,k
     
     real :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
     real :: dudxp,dudxm,dudyp,dudym,dvdxp,dvdxm,dudzp,dudzm,dwdxp,dwdxm
     real :: ma1,ma2,ma3,ma4
     real :: rhou  !#
     real :: RD1,RD2,RD3,RD4,RD5
     
     real, dimension(1:imax,1:jmax,1:kmax) :: RC,RD,RG

     !$omp parallel default(shared) &
     !$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
     !$omp&private(uuip,uuim,uvjp,uvjm,uwkp,uwkm,dudxp,dudxm,dudyp,dudym,dudzp,dudzm)
     !$omp do

     if (TwoD) then
        dudt = 0.
        return
     endif

     !-- Convection --!

     select case (NS_conv)

     case('CNTR2')

        do k = 1,kmax
           do j = 1,jmax
              do i = 1,imax
                 
                 ip = i + 1
                 jp = j + 1
                 kp = k + 1
                 im = i - 1
                 jm = j - 1
                 km = k - 1

                 uuip  = 0.25 *( u(ip,j ,k ) +u(i,j,k) )*( u(ip,j ,k ) +u(i,j ,k ) )
                 uuim  = 0.25 *( u(im,j ,k ) +u(i,j,k) )*( u(im,j ,k ) +u(i,j ,k ) )
                 uvjp  = 0.25 *( u(i ,jp,k ) +u(i,j,k) )*( v(ip,j ,k ) +v(i,j ,k ) )
                 uvjm  = 0.25 *( u(i ,jm,k ) +u(i,j,k) )*( v(ip,jm,k ) +v(i,jm,k ) )
                 uwkp  = 0.25 *( u(i ,j ,kp) +u(i,j,k) )*( w(ip,j ,k ) +w(i,j ,k ) )
                 uwkm  = 0.25 *( u(i ,j ,km) +u(i,j,k) )*( w(ip,j ,km) +w(i,j ,km) )

                 RC(i,j,k) = (-uuip +uuim)/dx + (-uvjp +uvjm)/dy + (-uwkp +uwkm)/dz
                 
              enddo
           enddo
        enddo


     case('WENO5')

        call WENO5_conv('x-dir',RC)

     end select

     !-- Diffusion --!

     do k = 1,kmax
        do j = 1,jmax
           do i = 1,imax

              ip = i + 1
              jp = j + 1
              kp = k + 1
              im = i - 1
              jm = j - 1
              km = k - 1

              ! velocity gradients

              dudxp = (u(ip,j,k)-u(i ,j,k))/dx
              dudxm = (u(i ,j,k)-u(im,j,k))/dx

              dudyp = (u(i,jp,k)-u(i,j ,k))/dy
              dudym = (u(i,j ,k)-u(i,jm,k))/dy

              dvdxp = (v(ip,j ,k)-v(i,j ,k))/dx     
              dvdxm = (v(ip,jm,k)-v(i,jm,k))/dx  

              dudzp = (u(i,j,kp)-u(i,j,k ))/dz
              dudzm = (u(i,j,k )-u(i,j,km))/dz

              dwdxp = (w(ip,j,k )-w(i,j,k ))/dx   
              dwdxm = (w(ip,j,km)-w(i,j,km))/dx

              ! averaging at face centers

              ma1  = 0.25 *( miu(i,j,k)+miu(i,jm,k )+miu(ip,jm,k )+miu(ip,j,k) )   ! miu(i+0.5,j-0.5,k)
              ma2  = 0.25 *( miu(i,j,k)+miu(i,jp,k )+miu(ip,jp,k )+miu(ip,j,k) )   ! miu(i+0.5,j+0.5,k)
              ma3  = 0.25 *( miu(i,j,k)+miu(i,j ,km)+miu(ip,j ,km)+miu(ip,j,k) )   ! miu(i+0.5,j,k-0.5)
              ma4  = 0.25 *( miu(i,j,k)+miu(i,j ,kp)+miu(ip, j,kp)+miu(ip,j,k) )   ! miu(i+0.5,j,k+0.5)

              RD1 = 2.*(dudxp*miu(ip,j,k)-dudxm*miu(i,j,k))/dx  ! d/dx(2*miu*du/dx)
              RD2 = (dudyp*ma2-dudym*ma1)/dy                    ! d/dy(miu*du/dy)
              RD3 = (dvdxp*ma2-dvdxm*ma1)/dy                    ! d/dy(miu*dv/dx)
              RD4 = (dudzp*ma4-dudzm*ma3)/dz                    ! d/dz(miu*du/dz)
              RD5 = (dwdxp*ma4-dwdxm*ma3)/dz                    ! d/dz(miu*dw/dx)

              rhou = (rho(i,j,k)+rho(ip,j,k))/2.

              RD(i,j,k) = (RD1+RD2+RD3+RD4+RD5)/rhou

           enddo
        enddo
     enddo

     !-- Gravity --!

     do k = 1,kmax
        do j = 1,jmax
           do i = 1,imax

              rhou = (rho(i,j,k)+rho(i+1,j,k))/2.
              
              RG(i,j,k) = iFr*(1.-rho_avr/rhou)*g_unit(1)

           enddo
        enddo
     enddo


     !-- Output  --!

     dudt(1:imax,1:jmax,1:kmax) = RC + RD*visc + RG
     
     !$omp end parallel

     return
   end subroutine momxad


   ! #2
   
   subroutine momyad(dvdt, u,v,w,miu,rho,g_unit)   

     real, intent(in) :: g_unit(2)
     real, dimension(0:,0:,0:), intent(in) :: u,v,w,miu,rho   
     real, dimension(0:,0:,0:), intent(out) :: dvdt

     integer :: im,ip,jm,jp,km,kp,i,j,k

     real :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
     real :: dvdxp,dvdxm,dudyp,dudym,dvdyp,dvdym,dvdzp,dvdzm,dwdyp,dwdym
     real :: ma2,ma5,ma6,ma7 
     real :: rhov  !#
     real :: RD1,RD2,RD3,RD4,RD5
     
     real, dimension(1:imax,1:jmax,1:kmax) :: RC,RD,RG

     !$omp parallel default(shared) &
     !$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
     !$omp&private(uvip,uvim,vvjp,vvjm,wvkp,wvkm,dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm)
     !$omp do


     !-- Convection --!

     select case (NS_conv)

     case('CNTR2')

        do k = 1,kmax
           do j = 1,jmax
              do i = 1,imax
                 
                 ip = i + 1
                 jp = j + 1
                 kp = k + 1
                 im = i - 1
                 jm = j - 1
                 km = k - 1

                 uvip  = 0.25 *(u(i ,j,k ) +u(i ,jp,k ) ) *(v(i,j,k ) +v(ip,j ,k) )
                 uvim  = 0.25 *(u(im,j,k ) +u(im,jp,k ) ) *(v(i,j,k ) +v(im,j ,k) )
                 vvjp  = 0.25 *(v(i ,j,k ) +v(i ,jp,k ) ) *(v(i,j,k ) +v(i ,jp,k) )
                 vvjm  = 0.25 *(v(i ,j,k ) +v(i ,jm,k ) ) *(v(i,j,k ) +v(i ,jm,k) )
                 wvkp  = 0.25 *(w(i ,j,k ) +w(i ,jp,k ) ) *(v(i,j,kp) +v(i ,j ,k) )
                 wvkm  = 0.25 *(w(i ,j,km) +w(i ,jp,km) ) *(v(i,j,km) +v(i ,j ,k) )

                 RC(i,j,k) = (-uvip +uvim)/dx + (-vvjp +vvjm)/dy + (-wvkp +wvkm)/dz
                 
              enddo
           enddo
        enddo


     case('WENO5')

        call WENO5_conv('y-dir',RC)

     end select

     !-- Diffusion --!

     do k = 1,kmax
        do j = 1,jmax
           do i = 1,imax
              
              ip = i + 1
              jp = j + 1
              kp = k + 1
              im = i - 1
              jm = j - 1
              km = k - 1

              ! velocity gradients

              dvdxp = (v(ip,j,k)-v(i ,j,k))/dx
              dvdxm = (v(i ,j,k)-v(im,j,k))/dx

              dudyp = (u(i ,jp,k)-u(i ,j,k))/dy     
              dudym = (u(im,jp,k)-u(im,j,k))/dy  

              dvdyp = (v(i,jp,k)-v(i,j ,k))/dy
              dvdym = (v(i,j ,k)-v(i,jm,k))/dy

              dvdzp = (v(i,j,kp)-v(i,j,k ))/dz
              dvdzm = (v(i,j,k )-v(i,j,km))/dz

              dwdyp = (w(i,jp,k )-w(i,j,k ))/dy    
              dwdym = (w(i,jp,km)-w(i,j,km))/dy  

              ! averaging at face centers

              ma2  = 0.25 *( miu(i,j,k)+miu(i,jp,k )+miu(ip,jp,k )+miu(ip,j ,k) )   ! miu(i+0.5,j+0.5,k)
              ma5  = 0.25 *( miu(i,j,k)+miu(i,jp,k )+miu(im,jp,k )+miu(im,j ,k) )   ! miu(i-0.5,j+0.5,k)
              ma6  = 0.25 *( miu(i,j,k)+miu(i,j ,km)+miu(i ,jp,km)+miu(i ,jp,k) )   ! miu(i,j+0.5,k-0.5)
              ma7  = 0.25 *( miu(i,j,k)+miu(i,j ,kp)+miu(i ,jp,kp)+miu(i ,jp,k) )   ! miu(i,j+0.5,k+0.5)      

              RD1 = (dvdxp*ma2-dvdxm*ma5)/dx                       ! d/dx(miu*dv/dx)
              RD2 = (dudyp*ma2-dudym*ma5)/dx                       ! d/dx(miu*du/dx)
              RD3 = 2.*(dvdyp*miu(i,jp,k)-dvdym*miu(i,j,k))/dy     ! d/dy(2*miu*dv/dy)
              RD4 = (dvdzp*ma7-dvdzm*ma6)/dz                       ! d/dz(miu*dv/dz)
              RD5 = (dwdyp*ma7-dwdym*ma6)/dz                       ! d/dz(miu*dw/dy)

              rhov = (rho(i,j,k)+rho(i,jp,k))/2.

              RD(i,j,k) = (RD1+RD2+RD3+RD4+RD5)/rhov

           enddo
        enddo
     enddo

     !-- Gravity --!

     do k = 1,kmax
        do j = 1,jmax
           do i = 1,imax

              rhov = (rho(i,j,k)+rho(i,j+1,k))/2.

              RG(i,j,k) = iFr*(1.-rho_avr/rhov)*g_unit(2)

           enddo
        enddo
     enddo
     !$omp end parallel


     !-- Output --!

     dvdt(1:imax,1:jmax,1:kmax) = RC + RD*visc + RG


     return
   end subroutine momyad


   ! #3
   
   subroutine momzad(dwdt, u,v,w,miu,rho,g_unit)  

     real, intent(in) :: g_unit(3)
     real, dimension(0:,0:,0:), intent(in) :: u,v,w,miu,rho   
     real, dimension(0:,0:,0:), intent(out) :: dwdt

     integer :: im,ip,jm,jp,km,kp,i,j,k
     
     real :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
     real :: dwdxp,dwdxm,dudzp,dudzm,dwdyp,dwdym,dvdzp,dvdzm,dwdzp,dwdzm
     real :: ma4,ma8,ma7,ma9 
     real :: rhow  !#
     real :: RD1,RD2,RD3,RD4,RD5
     
     real, dimension(1:imax,1:jmax,1:kmax) :: RC,RD,RG

     !$omp parallel default(shared) &
     !$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
     !$omp&private(uwip,uwim,vwjp,vwjm,wwkp,wwkm,dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm)
     !$omp do


     !-- Convection --!

     select case (NS_conv)

     case('CNTR2')

        do k = 1,kmax
           do j = 1,jmax
              do i = 1,imax
                 
                 ip = i + 1
                 jp = j + 1
                 kp = k + 1
                 im = i - 1
                 jm = j - 1
                 km = k - 1

                 uwip  = 0.25 *(w(i,j,k) +w(ip,j ,k ) ) *(u(i ,j ,k) +u(i ,j ,kp) )
                 uwim  = 0.25 *(w(i,j,k) +w(im,j ,k ) ) *(u(im,j ,k) +u(im,j ,kp) )
                 vwjp  = 0.25 *(w(i,j,k) +w(i ,jp,k ) ) *(v(i ,j ,k) +v(i ,j ,kp) )
                 vwjm  = 0.25 *(w(i,j,k) +w(i ,jm,k ) ) *(v(i ,jm,k) +v(i ,jm,kp) )
                 wwkp  = 0.25 *(w(i,j,k) +w(i ,j ,kp) ) *(w(i ,j ,k) +w(i ,j ,kp ))
                 wwkm  = 0.25 *(w(i,j,k) +w(i ,j ,km) ) *(w(i ,j ,k) +w(i ,j ,km ))

                 RC(i,j,k) = (-uwip +uwim)/dx + (-vwjp +vwjm)/dy + (-wwkp +wwkm)/dz
                 
              enddo
           enddo
        enddo


     case('WENO5')

        call WENO5_conv('z-dir',RC)

     end select

     
     !-- Diffusion --!

     do k = 1,kmax
        do j = 1,jmax
           do i = 1,imax
              
              ip = i + 1
              jp = j + 1
              kp = k + 1
              im = i - 1
              jm = j - 1
              km = k - 1

              ! velocity gradients

              dwdxp = (w(ip,j,k)-w(i ,j,k))/dx
              dwdxm = (w(i ,j,k)-w(im,j,k))/dx

              dudzp = (u(i ,j,kp)-u(i ,j,k))/dz   
              dudzm = (u(im,j,kp)-u(im,j,k))/dz  

              dwdyp = (w(i,jp,k)-w(i,j ,k))/dy
              dwdym = (w(i,j ,k)-w(i,jm,k))/dy

              dvdzp = (v(i,j ,kp)-v(i,j ,k))/dz     
              dvdzm = (v(i,jm,kp)-v(i,jm,k))/dz  

              dwdzp = (w(i,j,kp)-w(i,j,k ))/dz
              dwdzm = (w(i,j,k )-w(i,j,km))/dz

              ! averaging at cell faces

              ma4  = 0.25 *(miu(i,j,k)+miu(i,j,kp)+miu(ip,j ,kp)+miu(ip,j ,k) )   ! miu(i+0.5,j,k+0.5)
              ma8  = 0.25 *(miu(i,j,k)+miu(i,j,kp)+miu(im,j ,kp)+miu(im,j ,k) )   ! miu(i-0.5,j,k+0.5)
              ma7  = 0.25 *(miu(i,j,k)+miu(i,j,kp)+miu(i ,jp,kp)+miu(i ,jp,k) )   ! miu(i,j+0.5,k+0.5)
              ma9  = 0.25 *(miu(i,j,k)+miu(i,j,kp)+miu(i ,jm,kp)+miu(i ,jm,k) )   ! miu(i,j-0.5,k+0.5)

              RD1 = (dwdxp*ma4-dwdxm*ma8)/dx                    ! d/dx(miu*dw/dx)
              RD2 = (dudzp*ma4-dudzm*ma8)/dx                    ! d/dx(miu*du/dx)
              RD3 = (dwdyp*ma7-dwdym*ma9)/dy                    ! d/dy(miu*dw/dy)
              RD4 = (dvdzp*ma7-dvdzm*ma9)/dy                    ! d/dy(miu*dv/dy)
              RD5 = 2.*(dwdzp*miu(i,j,kp)-dwdzm*miu(i,j,k))/dz  ! d/dz(2*miu*dw/dz)

              rhow = (rho(i,j,k)+rho(i,j,kp))/2.

              RD(i,j,k) = (RD1+RD2+RD3+RD4+RD5)/rhow

           enddo
        enddo
     enddo              

     !-- Gravity --!
     
     do k = 1,kmax
        do j = 1,jmax
           do i = 1,imax     

              rhow = (rho(i,j,k)+rho(i,j,k+1))/2.
              
              RG(i,j,k) = iFr*(1.-rho_avr/rhow)*g_unit(3)

           enddo
        enddo
     enddo
     !$omp end parallel


     !-- Output --!

     dwdt(1:imax,1:jmax,1:kmax) = RC + RD*visc + RG


     return
   end subroutine momzad






!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!# convection of a scalar
subroutine momsad(dydt,u,v,w,y)
implicit none

real, dimension(0:,0:,0:), intent(in) :: u,v,w
real, dimension(-2:,-2:,-2:), intent(in) :: y
real, dimension(-2:,-2:,-2:), intent(out) :: dydt
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: yuip,yuim,yvjp,yvjm,ywkp,ywkm
real :: dydxp,dydxm,dydyp,dydym,dydzp,dydzm

!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(uwip,uwim,vwjp,vwjm,wwkp,wwkm,dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm)
!$omp do
 do k=1,kmax
   do j=1,jmax
     do i=1,imax
       ip = i + 1
       jp = j + 1
       kp = k + 1
       im = i - 1
       jm = j - 1
       km = k - 1
       yuip  = 0.5 * (Y(i,j,k)+Y(ip,j,k))*U(i,j,k)
       yuim  = 0.5 * (Y(i,j,k)+Y(im,j,k))*U(im,j,k)
       yvjp  = 0.5 * (Y(i,j,k)+Y(i,jp,k))*V(i,j,k)
       yvjm  = 0.5 * (Y(i,j,k)+Y(i,jm,k))*V(i,jm,k)
       ywkp  = 0.5 * (Y(i,j,k)+Y(i,j,kp))*W(i,j,k)
       ywkm  = 0.5 * (Y(i,j,k)+Y(i,j,km))*W(i,j,km)

       ! Momentum balance

       dydt(i,j,k) = dxi*( -yuip + yuim ) + &
                     dyi*( -yvjp + yvjm ) + &
                     dzi*( -ywkp + ywkm ) 
     enddo
   enddo
 enddo
!$omp end parallel

return
end subroutine momsad


end module mod_mom
