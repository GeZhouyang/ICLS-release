module mod_hyperPDE
  !
  ! This module is dedicated to those super long schemes for hyperbolic PDEs.
  !
  use mod_common
  use mod_common_mpi
  use mod_param

  implicit none

  private
  public UPWD2, semiLagrangian, HOUC5_nb, WENO5_nb, mWENO5c_nb, mWENO5nc_nb, &
         RussoSmereka_init,RussoSmereka_iterate_1, RussoSmereka_iterate_2, &
         RussoSmereka_init_4, RussoSmereka_iterate_4, &
         WENO5_reinit, WENO5z_reinit, WENO5_deri
  

contains


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  !  2nd-order Upwinding (unvalidated)
  !  Note that there is no order reduction near proc boundarier

  subroutine UPWD2(vel_u,vel_v,vel_w, phi, dfdt)
    
    !-- Advection velocities (cell-center)
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in) :: vel_u,vel_v,vel_w

    !-- Scalar to be transported
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in) :: phi

    !-- Temporal derivative
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(out) :: dfdt

    integer :: i,j,k
    real :: udfdx,vdfdy,wdfdz, dfx,dfy,dfz

    dfdt = 0.

    do k = -2,k1+2
       do j = -2,j1+2
          do i = -2,i1+2

             !-- Flux in x --!

             if (TwoD .or. vel_u(i,j,k) .eq. 0. ) then
                udfdx = 0.
             else
                if ( vel_u(i,j,k) .gt. 0. ) then                  
                   dfx = 0.5*phi(i-2,j,k) - 2.*phi(i-1,j,k) + 1.5*phi(i,j,k)                   
                else                  
                   dfx = -1.5*phi(i,j,k) + 2.*phi(i+1,j,k) - 0.5*phi(i+2,j,k)                  
                endif
                udfdx = -vel_u(i,j,k)*dfx/dx
             endif

             !-- Flux in y --!

             if ( vel_v(i,j,k) .eq. 0. ) then
                vdfdy = 0.
             else
                if ( vel_v(i,j,k) .gt. 0. ) then                  
                   dfy = 0.5*phi(i,j-2,k) - 2.*phi(i,j-1,k) + 1.5*phi(i,j,k)                 
                else                 
                   dfy = -1.5*phi(i,j,k) + 2.*phi(i,j+1,k) - 0.5*phi(i,j+2,k)               
                endif
                vdfdy = -vel_v(i,j,k)*dfy/dy
             endif
             
             !-- Flux in z --!

             if ( vel_w(i,j,k) .eq. 0. ) then
                wdfdz = 0.
             else
                if ( vel_w(i,j,k) .gt. 0. ) then                
                   dfz = 0.5*phi(i,j,k-2) - 2.*phi(i,j,k-1) + 1.5*phi(i,j,k)                  
                else                   
                   dfz = -1.5*phi(i,j,k) + 2.*phi(i,j,k+1) - 0.5*phi(i,j,k+2)                  
                endif
                wdfdz = -vel_w(i,j,k)*dfz/dz
             endif

             !---- Rate of Change ----!
   
             dfdt(i,j,k) = udfdx + vdfdy + wdfdz
             
          enddo
       enddo
    enddo
    
    
    return
  end subroutine UPWD2


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  !  Semi-Lagrangian scheme (implemented in a dimension-by-dimension fashion)
  !
  !  2D only
  
  subroutine semiLagrangian(l,phi, phi_new)

    integer, intent(in) :: l
    
    !-- Old level set
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in) :: phi

    !-- New level set
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(out) :: phi_new

    integer :: i,j,k, cnt
    integer :: im,ip,jm,jp,km,kp

    real :: vel_u,vel_v,vel_w
    real :: vel_ui,vel_vi,vel_wi
    real :: delta_x,delta_y,delta_z, xs,ys,zs
    real :: phi_i, phiyy_min,phizz_min, dyy,dzz

    real, dimension(0:1,0:1) :: coef
    real, dimension(0:i1,0:j1,0:k1) :: phixx,phiyy,phizz


    !---- Initialization ----!

    phi_new = phi
    
    !---- Spacetime integration ----!

    do cnt = 1,nb_cnt

       !-- narrow band indices
       
       i = nb_i(cnt)
       j = nb_j(cnt)
       k = nb_k(cnt)

       !-- current characteristics

       vel_u = uccnew(i,j,k)
       vel_v = vccnew(i,j,k)
       vel_w = wccnew(i,j,k)

       !-- move backwards, first RK step

       delta_x = - vel_u*timestep/2.
       delta_y = - vel_v*timestep/2.
       delta_z = - vel_w*timestep/2.

       !-- search for bounding cell (mid-position)
       
       call bound_search(i,delta_x,dx, im,ip,xs)
       call bound_search(j,delta_y,dy, jm,jp,ys)
       call bound_search(k,delta_z,dz, km,kp,zs)

       call bi_linear_coef(ys,zs, coef)

       !-- interpolate characteristics within cell

       call bi_linear_interp(i,jm,jp,km,kp,l, coef, 'u', vel_ui)
       call bi_linear_interp(i,jm,jp,km,kp,l, coef, 'v', vel_vi)
       call bi_linear_interp(i,jm,jp,km,kp,l, coef, 'w', vel_wi)

       !-- move backwards, second RK step
       
       delta_x = - vel_ui*timestep
       delta_y = - vel_vi*timestep
       delta_z = - vel_wi*timestep

       !-- search for bounding cell (final position)
       
       call bound_search(i,delta_x,dx, im,ip,xs)
       call bound_search(j,delta_y,dy, jm,jp,ys)
       call bound_search(k,delta_z,dz, km,kp,zs)

       call bi_linear_coef(ys,zs, coef)
       
       !-- interpolate level set within cell

       call bi_linear_interp(i,jm,jp,km,kp,l, coef, 'p', phi_i)

       !-- min second derivatives (2D on y-z plane, level set 1)

       call min_2nd_deriv(1, jm,jp,km,kp, phiyy_min,phizz_min)

       !-- stable quadratic interpolation (Min&Gibou JCP2007)

       dyy = - phiyy_min*ys*(1.-ys)/2.
       dzz = - phizz_min*zs*(1.-zs)/2.

       phi_i = phi_i + dyy +dzz
       
       !-- update
       
       phi_new(i,j,k) = phi_i
       
    enddo

    !---- Update old velocities ----!

    if (l .eq. lmax) then
       uccold = uccnew
       vccold = vccnew
       wccold = wccnew
    endif

    return
  end subroutine semiLagrangian

  !-- Coefficients of the bi-linear interpolation

  subroutine bi_linear_coef(ys,zs, coef)

    real, intent(in) :: ys,zs
    real, dimension(0:1,0:1), intent(out) :: coef

    coef(0,0) = (1.-ys)*(1.-zs)
    coef(0,1) = (1.-ys)*zs
    coef(1,0) = ys*(1.-zs)
    coef(1,1) = ys*zs

    return
  end subroutine bi_linear_coef

  !-- Direction-independent bi-linear interpolation

  subroutine bi_linear_interp(i,jm,jp,km,kp,l, coef, flag, qi)

    integer, intent(in) :: i,jm,jp,km,kp,l
    real, dimension(0:1,0:1), intent(in) :: coef
    character(len=1), intent(in) :: flag

    real, intent(out) :: qi

    real :: q00,q01,q10,q11

    select case (flag)
    case('u')
       
       q00 = 1.5*uccnew(i,jm,km) - 0.5*uccold(i,jm,km)
       q01 = 1.5*uccnew(i,jm,kp) - 0.5*uccold(i,jm,kp)
       q10 = 1.5*uccnew(i,jp,km) - 0.5*uccold(i,jp,km)
       q11 = 1.5*uccnew(i,jp,kp) - 0.5*uccold(i,jp,kp)
       
    case('v')

       q00 = 1.5*vccnew(i,jm,km) - 0.5*vccold(i,jm,km)
       q01 = 1.5*vccnew(i,jm,kp) - 0.5*vccold(i,jm,kp)
       q10 = 1.5*vccnew(i,jp,km) - 0.5*vccold(i,jp,km)
       q11 = 1.5*vccnew(i,jp,kp) - 0.5*vccold(i,jp,kp)

    case('w')

       q00 = 1.5*wccnew(i,jm,km) - 0.5*wccold(i,jm,km)
       q01 = 1.5*wccnew(i,jm,kp) - 0.5*wccold(i,jm,kp)
       q10 = 1.5*wccnew(i,jp,km) - 0.5*wccold(i,jp,km)
       q11 = 1.5*wccnew(i,jp,kp) - 0.5*wccold(i,jp,kp)
       
    case('p') ! for level set l

       q00 = lvset(i,jm,km,l)
       q01 = lvset(i,jm,kp,l)
       q10 = lvset(i,jp,km,l)
       q11 = lvset(i,jp,kp,l)    

    end select

    qi = coef(0,0)*q00 + coef(0,1)*q01 + coef(1,0)*q10 + coef(1,1)*q11
    
    return
  end subroutine bi_linear_interp

  !-- Direction-independent bounding cell search
  
  subroutine bound_search(ae,delta_ae,dae, aem,aep,aes)

    integer, intent(in) :: ae
    real, intent(in) :: delta_ae, dae
    
    integer, intent(out) :: aem,aep
    real, intent(out) :: aes

    aem = 0
    aep = 0
    aes = 0.

    if (abs(delta_ae) .le. 1.*dae) then ! within +/-1 cell
          
          if (delta_ae .gt. 0.) then
             aem = ae
             aep = ae+1
             aes = delta_ae
          else
             aem = ae-1
             aep = ae
             aes = dae + delta_ae
          endif

       elseif (abs(delta_ae) .le. 2.*dae) then ! within +/-2 cell

          if (delta_ae .gt. 0.) then
             aem = ae+1
             aep = ae+2
             aes = delta_ae - dae
          else
             aem = ae-2
             aep = ae-1
             aes = 2.*dae + delta_ae
          endif

       elseif (abs(delta_ae) .le. 3.*dae) then ! within +/-3 cell

          if (delta_ae .gt. 0.) then
             aem = ae+2
             aep = ae+3
             aes = delta_ae - 2.*dae
          else
             aem = ae-3
             aep = ae-2
             aes = 3.*dae + delta_ae
          endif

       else

          if (myid .eq. 0) then
             write(6,*) 'Failed to locate the bounding cell (semi-Lagrangian).'
             write(6,*) 'Perhaps the time step is too big.'
             write(6,*) 'Program aborted.'
          endif
          call decomp_2d_finalize
          call MPI_FINALIZE(error)
          stop

       endif

       aes = aes/dae
    
    return
  end subroutine bound_search

  !-- Evaluate the min second derivatives modulus of the bounding cell

  subroutine min_2nd_deriv(l, jm,jp,km,kp, yy_min,zz_min)

    integer, intent(in) :: l, jm,jp,km,kp

    real, intent(out) :: yy_min,zz_min

    integer :: i,j,k

    real, dimension(imax/2-1:imax/2+1,jm-1:jp+1,km-1:kp+1) :: phi
    real, dimension(imax/2,jm:jp,km:kp) :: phiyy,phizz
    real, dimension(10) :: x ! smoothed level set and its derivatives


    !- Initialization
    
    phi = lvset(imax/2-1:imax/2+1,jm-1:jp+1,km-1:kp+1,l)

    !phixx = 0.
    phiyy = 0.
    phizz = 0.

    !- Compute the absolute diagonal 2nd derivatives
    !- using least squares from +/-1 neighbours

    i = imax/2  ! 2D at the moment
    
    do k = km,kp
       do j = jm,jp

          x = matmul(pseudo_inv_A, reshape(phi(i-1:i+1,j-1:j+1,k-1:k+1), (/27/)))

          !phixx(i,j,k) = abs( x(5) )
          phiyy(i,j,k) = abs( x(6) )
          phizz(i,j,k) = abs( x(7) )

       enddo
    enddo

    !- Determine the minimal absolute derivatives

    yy_min = minval(phiyy)
    zz_min = minval(phizz)

    
    return
  end subroutine min_2nd_deriv

  
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! 
 !  High-Order Upstream Central, 5th-order

 subroutine HOUC5_nb(vel_u,vel_v,vel_w, phi, dfdt)

  !-- Advection velocities (cell-center)
  real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: vel_u,vel_v,vel_w

  !-- Scalar to be transported
  real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in) :: phi

  !-- Temporal derivative
  real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: dfdt

  integer :: i,j,k, cnt
  real :: udfdx,vdfdy,wdfdz, dfx,dfy,dfz

  dfdt = 0.
  
  do cnt = 1,nb_cnt

    i = nb_i(cnt)
    j = nb_j(cnt)
    k = nb_k(cnt)  
   
    !---- Flux in x ----!
   
    if (TwoD) then
   
      udfdx = 0.
  
    else
   
      if ( vel_u(i,j,k) .GE. 0. ) then              
        dfx = 1./60.*( -2.*phi(i-3,j,k)+15.*phi(i-2,j,k)-60.*phi(i-1,j,k)+ &
                       20.*phi(i,j,k)  +30.*phi(i+1,j,k)-3. *phi(i+2,j,k)  )         
      else             
        dfx = 1./60.*(  2.*phi(i+3,j,k)-15.*phi(i+2,j,k)+60.*phi(i+1,j,k)- &
                       20.*phi(i,j,k)  -30.*phi(i-1,j,k)+3. *phi(i-2,j,k)  )         
      endif
       
      udfdx = -vel_u(i,j,k)*dfx/dx
   
    endif
  
    !---- Flux in y ----!
  
    if ( vel_v(i,j,k) .GE. 0. ) then
      dfy = 1./60.*( -2.*phi(i,j-3,k)+15.*phi(i,j-2,k)-60.*phi(i,j-1,k)+ &
                     20.*phi(i,j,k)  +30.*phi(i,j+1,k)-3. *phi(i,j+2,k)  )        
    else
      dfy = 1./60.*(  2.*phi(i,j+3,k)-15.*phi(i,j+2,k)+60.*phi(i,j+1,k)- &
                     20.*phi(i,j,k)  -30.*phi(i,j-1,k)+3. *phi(i,j-2,k)  )          
    endif
   
    vdfdy = -vel_v(i,j,k)*dfy/dy
   
    !---- Flux in z ----!
   
    if (k .gt. 2 .and. k .lt. kmax-1) then   
   
      if ( vel_w(i,j,k) .GE. 0. ) then
        dfz = 1./60.*( -2.*phi(i,j,k-3)+15.*phi(i,j,k-2)-60.*phi(i,j,k-1)+ &
                       20.*phi(i,j,k)  +30.*phi(i,j,k+1)-3. *phi(i,j,k+2)  )          
      else
        dfz = 1./60.*(  2.*phi(i,j,k+3)-15.*phi(i,j,k+2)+60.*phi(i,j,k+1)- &
                       20.*phi(i,j,k)  -30.*phi(i,j,k-1)+3. *phi(i,j,k-2)  )                          
      endif
   
      wdfdz = -vel_w(i,j,k)*dfz/dz
   
    elseif (k .eq. 2 .or. k .eq. kmax-1) then  ! HOUC3 near boundaries
   
      if ( vel_w(i,j,k) .GE. 0. ) then
        dfz = (1./6.)*phi(i,j,k-2)-phi(i,j,k-1)+(1./2.)*phi(i,j,k)+(1./3.)*phi(i,j,k+1)  
      else
        dfz =-(1./6.)*phi(i,j,k+2)+phi(i,j,k+1)-(1./2.)*phi(i,j,k)-(1./3.)*phi(i,j,k-1)                          
      endif
   
      wdfdz = -vel_w(i,j,k)*dfz/dz
   
    else  ! first order upwinding
   
      if ( vel_w(i,j,k) .GE. 0. ) then
        dfz = phi(i,j,k) - phi(i,j,k-1)
      else
        dfz = phi(i,j,k+1) - phi(i,j,k)                    
      endif
   
      wdfdz = -vel_w(i,j,k)*dfz/dz
   
    endif

   
    !---- Rate of Change ----!
   
    dfdt(i,j,k) = udfdx + vdfdy + wdfdz

  enddo

  return
 end subroutine HOUC5_nb


 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! 
 !  Weighted Essentially-Non-Oscillation, 5th-order
 !  (3rd-order near discontinuities)

 subroutine WENO5_nb(vel_u,vel_v,vel_w, phi, dfdt)

  !-- Advection velocities (cell-center)
  real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: vel_u,vel_v,vel_w

  !-- Scalar to be transported
  real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in) :: phi

  !-- Temporal derivative
  real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: dfdt

  integer :: i,j,k, cnt
  real :: eps_weno
  parameter ( eps_weno = 1.e-6 )
  real :: sigma_1, sigma_2, sigma_3
  parameter ( sigma_1 = 0.1 )
  parameter ( sigma_2 = 0.6 )
  parameter ( sigma_3 = 0.3 )
  real :: v1, v2, v3, v4, v5
  real :: S1, S2, S3
  real :: a1, a2, a3
  real :: w1, w2, w3
  real :: udfdx,vdfdy,wdfdz, dfx,dfy,dfz

  dfdt = 0.

  do cnt = 1,nb_cnt

    i = nb_i(cnt)
    j = nb_j(cnt)
    k = nb_k(cnt)  
   
    !---- Flux in x ----!
   
    if (TwoD) then
   
      udfdx = 0.
  
    else
   
      if ( vel_u(i,j,k) .GE. 0. ) then    
          
        v1 = ( phi(i-2,j,k) - phi(i-3,j,k) ) / dx
        v2 = ( phi(i-1,j,k) - phi(i-2,j,k) ) / dx
        v3 = ( phi(i,j,k)   - phi(i-1,j,k) ) / dx
        v4 = ( phi(i+1,j,k) - phi(i,j,k)   ) / dx
        v5 = ( phi(i+2,j,k) - phi(i+1,j,k) ) / dx
       
        S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
             + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
       
        S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
             + (1./4.) * (v2 - v4)**2
       
        S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
             + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
       
        a1 = sigma_1 / (eps_weno + S1)**2
        a2 = sigma_2 / (eps_weno + S2)**2
        a3 = sigma_3 / (eps_weno + S3)**2
       
        w1 = a1 / (a1 + a2 + a3)
        w2 = a2 / (a1 + a2 + a3)
        w3 = a3 / (a1 + a2 + a3)
       
        dfx = w1 * ( (1./3.)  * v1   &
                   - (7./6.)  * v2   &
                   + (11./6.) * v3 ) &
            + w2 * ( (-1./6.) * v2   &
                   + (5./6.)  * v3   &
                   + (1./3.)  * v4 ) &
            + w3 * ( (1./3.)  * v3   &
                   + (5./6.)  * v4   &
                   - (1./6.)  * v5 )
      else             

        v1 = ( phi(i+3,j,k) - phi(i+2,j,k) ) / dx
        v2 = ( phi(i+2,j,k) - phi(i+1,j,k) ) / dx
        v3 = ( phi(i+1,j,k) - phi(i,j,k)   ) / dx
        v4 = ( phi(i,j,k)   - phi(i-1,j,k) ) / dx
        v5 = ( phi(i-1,j,k) - phi(i-2,j,k) ) / dx
       
       
        S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
             + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
       
        S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                     + (1./4.) * (v2 - v4)**2
       
        S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
             + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
       
        a1 = sigma_1 / (eps_weno + S1)**2
        a2 = sigma_2 / (eps_weno + S2)**2
        a3 = sigma_3 / (eps_weno + S3)**2
       
        w1 = a1 / (a1 + a2 + a3)
        w2 = a2 / (a1 + a2 + a3)
        w3 = a3 / (a1 + a2 + a3)
       
        dfx = w1 * ( (1./3.)  * v1   &
                   - (7./6.)  * v2   &
                   + (11./6.) * v3 ) &
            + w2 * ( (-1./6.) * v2   &
                   + (5./6.)  * v3   &
                   + (1./3.)  * v4 ) &
            + w3 * ( (1./3.)  * v3   &
                   + (5./6.)  * v4   &
                   - (1./6.)  * v5 )
      endif
       
      udfdx = -vel_u(i,j,k)*dfx
   
    endif
  
    !---- Flux in y ----!
  
    if ( vel_v(i,j,k) .GE. 0. ) then

      v1 = ( phi(i,j-2,k) - phi(i,j-3,k) ) / dx
      v2 = ( phi(i,j-1,k) - phi(i,j-2,k) ) / dx
      v3 = ( phi(i,j,k)   - phi(i,j-1,k) ) / dx
      v4 = ( phi(i,j+1,k) - phi(i,j,k)   ) / dx
      v5 = ( phi(i,j+2,k) - phi(i,j+1,k) ) / dx
      
      S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
           + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
                
      S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                     + (1./4.) * (v2 - v4)**2
                
      S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
           + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
                
      a1 = sigma_1 / (eps_weno + S1)**2
      a2 = sigma_2 / (eps_weno + S2)**2
      a3 = sigma_3 / (eps_weno + S3)**2
                
      w1 = a1 / (a1 + a2 + a3)
      w2 = a2 / (a1 + a2 + a3)
      w3 = a3 / (a1 + a2 + a3)
                
      dfy = w1 * ( (1./3.)  * v1   &
                 - (7./6.)  * v2   &
                 + (11./6.) * v3 ) &
          + w2 * ( (-1./6.) * v2   &
                 + (5./6.)  * v3   &
                 + (1./3.)  * v4 ) &
          + w3 * ( (1./3.)  * v3   &
                 + (5./6.)  * v4   &
                 - (1./6.)  * v5 )
    else

      v1 = ( phi(i,j+3,k) - phi(i,j+2,k) ) / dx
      v2 = ( phi(i,j+2,k) - phi(i,j+1,k) ) / dx
      v3 = ( phi(i,j+1,k) - phi(i,j,k)   ) / dx
      v4 = ( phi(i,j,k)   - phi(i,j-1,k) ) / dx
      v5 = ( phi(i,j-1,k) - phi(i,j-2,k) ) / dx
                
      S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
           + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
              
      S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                     + (1./4.) * (v2 - v4)**2
                
      S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
           + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
                
      a1 = sigma_1 / (eps_weno + S1)**2
      a2 = sigma_2 / (eps_weno + S2)**2
      a3 = sigma_3 / (eps_weno + S3)**2
                
      w1 = a1 / (a1 + a2 + a3)
      w2 = a2 / (a1 + a2 + a3)
      w3 = a3 / (a1 + a2 + a3)
      
      dfy = w1 * ( (1./3.)  * v1   &
                 - (7./6.)  * v2   &
                 + (11./6.) * v3 ) &
          + w2 * ( (-1./6.) * v2   &
                 + (5./6.)  * v3   &
                 + (1./3.)  * v4 ) &
          + w3 * ( (1./3.)  * v3   &
                 + (5./6.)  * v4   &
                 - (1./6.)  * v5 )
    endif
   
    vdfdy = -vel_v(i,j,k)*dfy
   
    !---- Flux in z ----!
   
    if (k .gt. 2 .and. k .lt. kmax-1) then   
   
      if ( vel_w(i,j,k) .GE. 0. ) then

        v1 = ( phi(i,j,k-2) - phi(i,j,k-3) ) / dz
        v2 = ( phi(i,j,k-1) - phi(i,j,k-2) ) / dz
        v3 = ( phi(i,j,k)   - phi(i,j,k-1) ) / dz
        v4 = ( phi(i,j,k+1) - phi(i,j,k)   ) / dz
        v5 = ( phi(i,j,k+2) - phi(i,j,k+1) ) / dz
        
        S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
             + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
        
        S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
             + (1./4.) * (v2 - v4)**2
        
        S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
             + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
        
        a1 = sigma_1 / (eps_weno + S1)**2
        a2 = sigma_2 / (eps_weno + S2)**2
        a3 = sigma_3 / (eps_weno + S3)**2
        
        w1 = a1 / (a1 + a2 + a3)
        w2 = a2 / (a1 + a2 + a3)
        w3 = a3 / (a1 + a2 + a3)
        
        dfz = w1 * ( (1./3.)  * v1   &
                   - (7./6.)  * v2   &
                   + (11./6.) * v3 ) &
            + w2 * ( (-1./6.) * v2   &
                   + (5./6.)  * v3   &
                   + (1./3.)  * v4 ) &
            + w3 * ( (1./3.)  * v3   &
                   + (5./6.)  * v4   &
                   - (1./6.)  * v5 )
      else

        v1 = ( phi(i,j,k+3) - phi(i,j,k+2) ) / dz
        v2 = ( phi(i,j,k+2) - phi(i,j,k+1) ) / dz
        v3 = ( phi(i,j,k+1) - phi(i,j,k)   ) / dz
        v4 = ( phi(i,j,k)   - phi(i,j,k-1) ) / dz
        v5 = ( phi(i,j,k-1) - phi(i,j,k-2) ) / dz
        
        S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
             + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
        
        S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
             + (1./4.) * (v2 - v4)**2
        
        S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
             + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
        
        a1 = sigma_1 / (eps_weno + S1)**2
        a2 = sigma_2 / (eps_weno + S2)**2
        a3 = sigma_3 / (eps_weno + S3)**2
        
        w1 = a1 / (a1 + a2 + a3)
        w2 = a2 / (a1 + a2 + a3)
        w3 = a3 / (a1 + a2 + a3)
        
        dfz = w1 * ( (1./3.)  * v1   &
                   - (7./6.)  * v2   &
                   + (11./6.) * v3 ) &
            + w2 * ( (-1./6.) * v2   &
                   + (5./6.)  * v3   &
                   + (1./3.)  * v4 ) &
            + w3 * ( (1./3.)  * v3   &
                   + (5./6.)  * v4   &
                   - (1./6.)  * v5 )
      endif
   
      wdfdz = -vel_w(i,j,k)*dfz
   
    elseif (k .eq. 2 .or. k .eq. kmax-1) then  ! HOUC3 near boundaries
   
      if ( vel_w(i,j,k) .GE. 0. ) then

        dfz = (1./6.)*phi(i,j,k-2)-phi(i,j,k-1)+(1./2.)*phi(i,j,k)+(1./3.)*phi(i,j,k+1)
      else

        dfz =-(1./6.)*phi(i,j,k+2)+phi(i,j,k+1)-(1./2.)*phi(i,j,k)-(1./3.)*phi(i,j,k-1)
      endif
   
      wdfdz = -vel_w(i,j,k)*dfz/dz
   
    else  ! first order upwinding
   
      if ( vel_w(i,j,k) .GE. 0. ) then
        dfz = phi(i,j,k) - phi(i,j,k-1)
      else
        dfz = phi(i,j,k+1) - phi(i,j,k)                    
      endif
   
      wdfdz = -vel_w(i,j,k)*dfz/dz
   
    endif

   
    !---- Rate of Change ----!
   
    dfdt(i,j,k) = udfdx + vdfdy + wdfdz

  enddo

  return
 end subroutine WENO5_nb


 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! 
 !  Modified WENO5, in the conservative form
 !  (Kateryna PhD thesis, 2002). 
 
 subroutine mWENO5c_nb(u,v,w, phi, dfdt)
 
  real, dimension(0:i1,0:j1,0:k1), intent(in) :: u,v,w
  real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in) :: phi
 
  !---- Time derivative, div(u*phi)
  
  real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: dfdt
  
  !----- Parameters
  
  integer :: i, j, k, cnt, m
  real :: eps_weno
  parameter ( eps_weno = 1.e-6 )
  real :: A, B, dfx, dfy, dfz, vel
  real, dimension(-1:0) :: wind
  real, dimension(5) :: va, vb
  real, dimension(3) :: Sm, Rm, Pw
  real, dimension(3) :: sigma


  sigma(1) = 1.
  sigma(2) = 6.
  sigma(3) = 3.

  dfdt = 0.

  do cnt = 1,nb_cnt

    i = nb_i(cnt)
    j = nb_j(cnt)
    k = nb_k(cnt)  
 
  
!== x ==!

    if (TwoD) then

      dfx = 0.

    else

      wind(-1:0) = u(i-1:i,j,k)  ! velocity buffer

      !-- i+1/2 --!

      if (wind(0) .ge. 0.) then  ! for u > 0

        do m = 1,5
          va(m) = phi(i-3+m,j,k)  ! level set buffer
        enddo

        !-- Smoothness indicator
    
        Sm(1) = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
               + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2) = 13.*( va(2) -2.*va(3) + va(4) )**2 &
               + 3.*( va(2) -   va(4) )**2
    
        Sm(3) = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
               + 3.*(3.*va(3) -4.*va(4) + va(5) )**2
    
        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        A = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

      else  ! for u < 0

        do m = 1,5
          va(m) = phi(i+4-m,j,k)
        enddo

        !-- Smoothness indicator

        Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
                + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
                + 3.*( va(2) -   va(4) )**2
    
        Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
                + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        A = sum( Pw(:)*Rm(:) )/sum( Pw(:) ) 

      endif

      !-- i-1/2 --!

      if (wind(-1) .ge. 0.) then  ! for u > 0
      
        do m = 1,5
          va(m) = phi(i-4+m,j,k)  ! level set buffer
        enddo

        !-- Smoothness indicator

        Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
                + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
                + 3.*( va(2) -   va(4) )**2
    
        Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
                + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        B = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

      else  ! u < 0

        do m = 1,5
          va(m) = phi(i+3-m,j,k)  ! level set buffer
        enddo

        !-- Smoothness indicator

        Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
                + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
                + 3.*( va(2) -   va(4) )**2
    
        Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
                + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        B = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

      endif

      dfx = (wind(0)*A - wind(-1)*B)/dx

    endif


!== y ==!
         
    wind(-1:0) = v(i,j-1:j,k)  ! velocity buffer

    !-- j+1/2 --!

    if (wind(0) .ge. 0.) then  ! for v > 0

      do m = 1,5
        va(m) = phi(i,j-3+m,k)  ! level set buffer
      enddo

      !-- Smoothness indicator
    
      Sm(1) = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
             + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
      Sm(2) = 13.*( va(2) -2.*va(3) + va(4) )**2 &
             + 3.*( va(2) -   va(4) )**2
    
      Sm(3) = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
             + 3.*(3.*va(3) -4.*va(4) + va(5) )**2
    
      Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
      Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
      Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

      !-- Weights

      Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

      !-- Undivided difference

      A = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

    else  ! for v < 0

      do m = 1,5
        va(m) = phi(i,j+4-m,k)
      enddo

      !-- Smoothness indicator

      Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
              + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
      Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
              + 3.*( va(2) -   va(4) )**2
    
      Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
              + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

      Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
      Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
      Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

      !-- Weights

      Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

      !-- Undivided difference

      A = sum( Pw(:)*Rm(:) )/sum( Pw(:) ) 

    endif

    !-- j-1/2 --!

    if (wind(-1) .ge. 0.) then  ! for v > 0
    
      do m = 1,5
        va(m) = phi(i,j-4+m,k)  ! level set buffer
      enddo

      !-- Smoothness indicator

      Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
              + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
      Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
              + 3.*( va(2) -   va(4) )**2
    
      Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
              + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

      Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
      Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
      Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

      !-- Weights

      Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

      !-- Undivided difference

      B = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

    else  ! v < 0

      do m = 1,5
        va(m) = phi(i,j+3-m,k)  ! level set buffer
      enddo

      !-- Smoothness indicator

      Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
              + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
      Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
              + 3.*( va(2) -   va(4) )**2
    
      Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
              + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

      Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
      Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
      Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

      !-- Weights

      Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

      !-- Undivided difference

      B = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

    endif

    dfy = (wind(0)*A - wind(-1)*B)/dy

!== z ==!

    wind(-1:0) = w(i,j,k-1:k)  ! velocity buffer

    if ( ( k .GE. 3 ) .AND. ( k .LE. kmax-2 ) ) then

      !-- k+1/2 --!

      if (wind(0) .ge. 0.) then  ! w > 0

        do m = 1,5
          va(m) = phi(i,j,k-3+m)  ! level set buffer
        enddo

        !-- Smoothness indicator
    
        Sm(1) = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
               + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2) = 13.*( va(2) -2.*va(3) + va(4) )**2 &
               + 3.*( va(2) -   va(4) )**2
    
        Sm(3) = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
               + 3.*(3.*va(3) -4.*va(4) + va(5) )**2
    
        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        A = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

      else  ! w < 0

        do m = 1,5
          va(m) = phi(i,j,k+4-m)
        enddo

        !-- Smoothness indicator

        Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
                + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
                + 3.*( va(2) -   va(4) )**2
    
        Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
                + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        A = sum( Pw(:)*Rm(:) )/sum( Pw(:) ) 

      endif

      !-- k-1/2 --!

      if (wind(-1) .ge. 0.) then  ! w > 0
      
        do m = 1,5
          va(m) = phi(i,j,k-4+m)  ! level set buffer
        enddo

        !-- Smoothness indicator

        Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
                + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
                + 3.*( va(2) -   va(4) )**2
    
        Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
                + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        B = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

      else  ! w < 0

        do m = 1,5
          va(m) = phi(i,j,k+3-m)  ! level set buffer
        enddo

        !-- Smoothness indicator

        Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
                + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
                + 3.*( va(2) -   va(4) )**2
    
        Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
                + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        B = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

      endif

      dfz = (wind(0)*A - wind(-1)*B)/dz


    else  !# 1st order upwinding

      vel = sum(wind)/2.

      if (vel .ge. 0.) then
        dfz = vel*( phi(i,j,k) - phi(i,j,k-1) )/dz
      else
        dfz = vel*( phi(i,j,k+1) - phi(i,j,k) )/dz
      endif

    endif


    !-- Time derivative --!

    dfdt(i,j,k) = - dfx - dfy -dfz

  enddo
     
  return
 end subroutine mWENO5c_nb



 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! 
 !  Modified WENO5, in the non-conservative form
 
 subroutine mWENO5nc_nb(u,v,w, phi, dfdt)
 
  real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: u,v,w
  real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in) :: phi
 
  !---- Time derivative, div(u*phi)
  
  real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: dfdt
  
  !----- Parameters
  
  integer :: i, j, k, cnt, m
  real :: eps_weno
  parameter ( eps_weno = 1.e-6 )
  real :: A, B, dfx, dfy, dfz
  real, dimension(5) :: va, vb
  real, dimension(3) :: Sm, Rm, Pw
  real, dimension(3) :: sigma


  sigma(1) = 1.
  sigma(2) = 6.
  sigma(3) = 3.

  dfdt = 0.

  do cnt = 1,nb_cnt

    i = nb_i(cnt)
    j = nb_j(cnt)
    k = nb_k(cnt)  
 
  
!== x ==!

    if (TwoD) then

      dfx = 0.

    else

      if (u(i,j,k) .ge. 0.) then  ! for u > 0

        !-- i+1/2 --!

        do m = 1,5
          va(m) = phi(i-3+m,j,k)  ! level set buffer
        enddo

        !-- Smoothness indicator
    
        Sm(1) = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
               + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2) = 13.*( va(2) -2.*va(3) + va(4) )**2 &
               + 3.*( va(2) -   va(4) )**2
    
        Sm(3) = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
               + 3.*(3.*va(3) -4.*va(4) + va(5) )**2
    
        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        A = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

        !-- i-1/2 --!

        do m = 1,5
          va(m) = phi(i-4+m,j,k)  ! level set buffer
        enddo

        !-- Smoothness indicator

        Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
                + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
                + 3.*( va(2) -   va(4) )**2
    
        Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
                + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        B = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

      else  ! for u < 0

        !-- i+1/2 --!

        do m = 1,5
          va(m) = phi(i+4-m,j,k)
        enddo

        !-- Smoothness indicator

        Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
                + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
                + 3.*( va(2) -   va(4) )**2
    
        Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
                + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        A = sum( Pw(:)*Rm(:) )/sum( Pw(:) ) 

        !-- i-1/2 --!

        do m = 1,5
          va(m) = phi(i+3-m,j,k)  ! level set buffer
        enddo

        !-- Smoothness indicator

        Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
                + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
                + 3.*( va(2) -   va(4) )**2
    
        Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
                + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        B = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

      endif

      dfx = u(i,j,k)*(A - B)/dx

    endif


!== y ==!         

    if (v(i,j,k) .ge. 0.) then  ! for v > 0

      !-- j+1/2 --!

      do m = 1,5
        va(m) = phi(i,j-3+m,k)  ! level set buffer
      enddo

      !-- Smoothness indicator
    
      Sm(1) = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
             + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
      Sm(2) = 13.*( va(2) -2.*va(3) + va(4) )**2 &
             + 3.*( va(2) -   va(4) )**2
    
      Sm(3) = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
             + 3.*(3.*va(3) -4.*va(4) + va(5) )**2
    
      Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
      Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
      Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

      !-- Weights

      Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

      !-- Undivided difference

      A = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

      !-- j-1/2 --!

      do m = 1,5
        va(m) = phi(i,j-4+m,k)  ! level set buffer
      enddo

      !-- Smoothness indicator

      Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
              + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
      Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
              + 3.*( va(2) -   va(4) )**2
    
      Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
              + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

      Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
      Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
      Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

      !-- Weights

      Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

      !-- Undivided difference

      B = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

    else  ! for v < 0

      !-- j+1/2 --!

      do m = 1,5
        va(m) = phi(i,j+4-m,k)
      enddo

      !-- Smoothness indicator

      Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
              + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
      Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
              + 3.*( va(2) -   va(4) )**2
    
      Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
              + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

      Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
      Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
      Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

      !-- Weights

      Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

      !-- Undivided difference

      A = sum( Pw(:)*Rm(:) )/sum( Pw(:) ) 

      !-- j-1/2 --!

      do m = 1,5
        va(m) = phi(i,j+3-m,k)  ! level set buffer
      enddo

      !-- Smoothness indicator

      Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
              + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
      Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
              + 3.*( va(2) -   va(4) )**2
    
      Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
              + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

      Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
      Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
      Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

      !-- Weights

      Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

      !-- Undivided difference

      B = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

    endif

    dfy = v(i,j,k)*(A - B)/dy


!== z ==!

    if ( ( k .GE. 3 ) .AND. ( k .LE. kmax-2 ) ) then

      if (w(i,j,k) .ge. 0.) then  ! w > 0

        !-- k+1/2 --!

        do m = 1,5
          va(m) = phi(i,j,k-3+m)  ! level set buffer
        enddo

        !-- Smoothness indicator
    
        Sm(1) = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
               + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2) = 13.*( va(2) -2.*va(3) + va(4) )**2 &
               + 3.*( va(2) -   va(4) )**2
    
        Sm(3) = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
               + 3.*(3.*va(3) -4.*va(4) + va(5) )**2
    
        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        A = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

        !-- k-1/2 --!

        do m = 1,5
          va(m) = phi(i,j,k-4+m)  ! level set buffer
        enddo

        !-- Smoothness indicator

        Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
                + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
                + 3.*( va(2) -   va(4) )**2
    
        Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
                + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        B = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

      else  ! w < 0

        !-- k+1/2 --!

        do m = 1,5
          va(m) = phi(i,j,k+4-m)
        enddo

        !-- Smoothness indicator

        Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
                + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
                + 3.*( va(2) -   va(4) )**2
    
        Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
                + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        A = sum( Pw(:)*Rm(:) )/sum( Pw(:) ) 

        !-- k-1/2 --!

        do m = 1,5
          va(m) = phi(i,j,k+3-m)  ! level set buffer
        enddo

        !-- Smoothness indicator

        Sm(1)  = 13.*( va(1) -2.*va(2) +   va(3) )**2 &
                + 3.*( va(1) -4.*va(2) +3.*va(3) )**2
    
        Sm(2)  = 13.*( va(2) -2.*va(3) + va(4) )**2 &
                + 3.*( va(2) -   va(4) )**2
    
        Sm(3)  = 13.*(   va(3) -2.*va(4) + va(5) )**2 &
                + 3.*(3.*va(3) -4.*va(4) + va(5) )**2

        Rm(1) =  1./3.*va(1) -7./6.*va(2) +11./6.*va(3)
        Rm(2) = -1./6.*va(2) +5./6.*va(3) +1./3. *va(4)
        Rm(3) =  1./3.*va(3) +5./6.*va(4) -1./6. *va(5)

        !-- Weights

        Pw(:) = sigma(:)/( eps_weno + Sm(:) )**2

        !-- Undivided difference

        B = sum( Pw(:)*Rm(:) )/sum( Pw(:) )

      endif

      dfz = w(i,j,k)*(A - B)/dz


    else  !# 1st order upwinding

      if (w(i,j,k) .ge. 0.) then
        dfz = w(i,j,k)*( phi(i,j,k) - phi(i,j,k-1) )/dz
      else
        dfz = w(i,j,k)*( phi(i,j,k+1) - phi(i,j,k) )/dz
      endif

    endif


    !-- Time derivative --!

    dfdt(i,j,k) = - dfx - dfy -dfz

  enddo
     
  return
 end subroutine mWENO5nc_nb


 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !
 !   Redistance discretization considering the interface location
 !   i.e. sub-cell fix (Russo & Smereka JCP2000)

 subroutine RussoSmereka_init(phi_0, sgn,region,Dist)

   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in)  :: phi_0
   real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: sgn, region, Dist
   
   integer :: i,j,k

   real :: pijk,p_im,p_ip,p_jm,p_jp,p_km,p_kp
   real :: s
   real :: d1,d2,d3,d4,d5,d6,d7,d8,d9,d0
   real, parameter :: eps = 1e-6
   
   
   
   do k = 1,kmax
      do j = 1,jmax
         do i = 1,imax

            !-- cell p0(i,j,k) and its neighbors

            pijk = phi_0(i  ,j,k)
            p_im = phi_0(i-1,j,k)
            p_ip = phi_0(i+1,j,k)
            p_jm = phi_0(i,j-1,k)
            p_jp = phi_0(i,j+1,k)
            p_km = phi_0(i,j,k-1)
            p_kp = phi_0(i,j,k+1)

            !-- sign function of p0

            if (pijk .gt. 0.) then
               sgn(i,j,k) = 1.
            elseif(pijk .lt. 0.) then
               sgn(i,j,k) = -1.
            else
               sgn(i,j,k) = 0.
            endif

            !-- determine if p0 is close to the interface

            if ( pijk*p_im .lt. 0. .or. pijk*p_ip .lt. 0. .or. &
                 pijk*p_jm .lt. 0. .or. pijk*p_jp .lt. 0. .or. &
                 pijk*p_km .lt. 0. .or. pijk*p_kp .lt. 0.) then

               region(i,j,k) = 1.

               ! compute a signed distance

               d1 = abs(p_ip - p_im)/2.
               d2 = abs(p_jp - p_jm)/2.
               d3 = abs(p_kp - p_km)/2.
               d4 = abs(p_ip - pijk)
               d5 = abs(p_jp - pijk)
               d6 = abs(p_kp - pijk)
               d7 = abs(pijk - p_im)
               d8 = abs(pijk - p_jm)
               d9 = abs(pijk - p_km)

               !d0 = max(d1,d2,d3,d4,d5,d6,d7,d8,d9,eps)  ! 1D version (surprisingly works for 2D as well)
               d0 = max( sqrt(d1**2+d2**2+d3**2),sqrt(d4**2+d5**2+d6**2),sqrt(d7**2+d8**2+d9**2),eps )

               Dist(i,j,k) = dz*pijk/d0

            else

               region(i,j,k) = 0.
               Dist(i,j,k) = 0.  ! dummy

            endif

         enddo
      enddo
   enddo


   return
 end subroutine RussoSmereka_init


 subroutine RussoSmereka_iterate_1(sgn,region,Dist,phi, RHS)
   
   real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: sgn, region, Dist
   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in)  :: phi

   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(out) :: RHS

   integer :: i,j,k

   real :: s,pijk
   real :: a,b,c,d,e,f
   real :: ap,bp,cp,dp,ep,fp
   real :: am,bm,cm,dm,em,fm
   real :: apbm, cpdm, epfm
   real :: ambp, cmdp, emfp
   real :: G

   RHS = 0.

   do k = 1,kmax
      do j = 1,jmax
         do i = 1,imax

            s = sgn(i,j,k)
            pijk = phi(i,j,k)
           
            if (region(i,j,k) .eq. 1.) then  ! sub-cell fix

               RHS(i,j,k) = -1./dx*(s*abs(pijk)-Dist(i,j,k))

            else  ! Godunov scheme

               a = (phi(i,j,k)   - phi(i-1,j,k))/dx
               b = (phi(i+1,j,k) - phi(i,j,k)  )/dx
               c = (phi(i,j,k)   - phi(i,j-1,k))/dy
               d = (phi(i,j+1,k) - phi(i,j,k)  )/dy
               e = (phi(i,j,k)   - phi(i,j,k-1))/dz
               f = (phi(i,j,k+1) - phi(i,j,k)  )/dz

               ap = max(a,0.)
               bp = max(b,0.)
               cp = max(c,0.)
               dp = max(d,0.)
               ep = max(e,0.)
               fp = max(f,0.)

               am = min(a,0.)
               bm = min(b,0.)
               cm = min(c,0.)
               dm = min(d,0.)
               em = min(e,0.)
               fm = min(f,0.)

               apbm = max(ap**2,bm**2)
               cpdm = max(cp**2,dm**2)
               epfm = max(ep**2,fm**2)

               ambp = max(am**2,bp**2)
               cmdp = max(cm**2,dp**2)
               emfp = max(em**2,fp**2)

               if (s .gt. 0.) then
                  G = sqrt(apbm + cpdm + epfm) -1.
               else
                  G = sqrt(ambp + cmdp + emfp) -1.
               endif

               RHS(i,j,k) = -s*G

            endif

         enddo
      enddo
   enddo
            

   return
 end subroutine RussoSmereka_iterate_1


 subroutine RussoSmereka_iterate_2(sgn,region,Dist,phi, RHS)
   
   real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: sgn, region, Dist
   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in)  :: phi

   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(out) :: RHS

   integer :: i,j,k,m

   real :: s,pijk
   real :: a,b,c,d,e,f
   real :: ap,bp,cp,dp,ep,fp
   real :: am,bm,cm,dm,em,fm
   real :: apbm, cpdm, epfm
   real :: ambp, cmdp, emfp
   real :: G
   real, dimension(-2:2) :: x_p,f_p

   RHS = 0.

   do k = 1,kmax
      do j = 1,jmax
         do i = 1,imax

            s = sgn(i,j,k)
            pijk = phi(i,j,k)

            if (region(i,j,k) .eq. 1.) then  ! sub-cell fix

               RHS(i,j,k) = -1./dx*(s*abs(pijk)-Dist(i,j,k))

            else  ! second order ENO scheme
               
               do m = -2,2
                  x_p(m) = m*dx
                  f_p(m) = phi(i+m,j,k)
               enddo
               call ENO2_deri(x_p,f_p, a,b)  ! x-dir
               
               do m = -2,2
                  x_p(m) = m*dy
                  f_p(m) = phi(i,j+m,k)
               enddo
               call ENO2_deri(x_p,f_p, c,d)  ! y-dir
               
               do m = -2,2
                  x_p(m) = m*dz
                  f_p(m) = phi(i,j,k+m)
               enddo
               call ENO2_deri(x_p,f_p, e,f)  ! z-dir

               ap = max(a,0.)
               bp = max(b,0.)
               cp = max(c,0.)
               dp = max(d,0.)
               ep = max(e,0.)
               fp = max(f,0.)
               
               am = min(a,0.)
               bm = min(b,0.)
               cm = min(c,0.)
               dm = min(d,0.)
               em = min(e,0.)
               fm = min(f,0.)

               apbm = max(ap**2,bm**2)
               cpdm = max(cp**2,dm**2)
               epfm = max(ep**2,fm**2)
               
               ambp = max(am**2,bp**2)
               cmdp = max(cm**2,dp**2)
               emfp = max(em**2,fp**2)
               
               if (s .gt. 0.) then
                  G = sqrt(apbm + cpdm + epfm) -1.
               else
                  G = sqrt(ambp + cmdp + emfp) -1.
               endif

               RHS(i,j,k) = -s*G

            endif
            
         enddo
      enddo
   enddo
            

   return
 end subroutine RussoSmereka_iterate_2


 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !
 !  4th order extension of Russo & Smereka (du Chene etal JSC2000)

 subroutine RussoSmereka_init_4(phi_0, sgn,region,dist_a, flag_xo,flag_yo,flag_zo, dist_xo,dist_yo,dist_zo)

   use mod_bound
 
   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in)  :: phi_0
 
   real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: sgn,region,dist_a
   real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: flag_xo,flag_yo,flag_zo
   real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: dist_xo,dist_yo,dist_zo
   
   integer :: i,j,k, m

   logical :: notzero
   real :: pijk,p_im,p_ip,p_jm,p_jp,p_km,p_kp, pabs,pg
   real :: s, dist, eps
   real :: d1,d2,d3,d4,d5,d6,d7,d8,d9,d0
   real, dimension(0:3) :: ph
   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: flag_x,flag_y,flag_z, shift
   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: dist_x,dist_y,dist_z
 
 
   dist_x = 0.
   dist_y = 0.
   dist_z = 0.
   flag_x = 0.
   flag_y = 0.
   flag_z = 0.
   shift  = 1.
   dist_a = 0.

   eps = 1e-6
   
   do k = 1,kmax
      do j = 1,jmax
         do i = 1,imax
            
            pijk = phi_0(i  ,j,k)
            p_im = phi_0(i-1,j,k)
            p_ip = phi_0(i+1,j,k)
            p_jm = phi_0(i,j-1,k)
            p_jp = phi_0(i,j+1,k)
            p_km = phi_0(i,j,k-1)
            p_kp = phi_0(i,j,k+1)
 
            !-- sign function
 
            if (pijk .gt. 0.) then
               sgn(i,j,k) = 1.
            elseif(pijk .lt. 0.) then
               sgn(i,j,k) = -1.
            else
               sgn(i,j,k) = 0.
            endif
 
            !-- segment length to the interpolated interface in each direction

            pabs = abs(pijk)
            pg = 0.
 
            if (pabs .lt. 6.*dx) then ! near the interface
 
               !-- x-dir

               flag_x(i,j,k) = 0.
               dist_x(i,j,k) = 0.
               
               !-- y-dir
               
               if (pijk*p_jp .lt. 0.) then

                  !-- interpolate distance between j and interface
                  
                  ph(0:3) = phi_0(i,j-1:j+2,k)
                  call cubic_int(ph,dist)
                  dist = dist - dy

                  !-- update distance for neighbors of j
                  
                  do m = -3,3  ! -3 is dummy
                     dist_y(i,j+m,k) = dist - m*dy
                  enddo

                  !-- tag point to replace

                  if (pabs .lt. dy**2) then  ! interface too close to j
                     pg = 1. !#
                     do m = -2,3
                        flag_y(i,j+m,k) = -1.*m
                        shift(i,j+m,k) = 0.
                     enddo
                  elseif (pabs .gt. dy-dy**2) then  ! interface too close to j+1
                     do m = -2,3
                        flag_y(i,j+m,k) = -1.*m + 1.
                        shift(i,j+m,k) = 0.
                     enddo
                  else  ! interface in between
                     do m = -2,0
                        flag_y(i,j+m,k) = -1.*m + 1.
                     enddo
                     do m = 1,3
                        flag_y(i,j+m,k) = -1.*m
                     enddo
                  endif
                  
               endif

               if (pijk*p_jm .lt. 0. .and. pabs .lt. dy**2) pg = 1. !#

               
               
               !-- z-dir

               if (pijk*p_kp .lt. 0.) then

                  !-- interpolate distance between k and interface
                  
                  ph(0:3) = phi_0(i,j,k-1:k+2)
                  call cubic_int(ph,dist)
                  dist = dist - dz

                  !-- update distance for neighbors of k
                  
                  do m = -3,3  ! -3 is dummy
                     dist_z(i,j,k+m) = dist - m*dz
                  enddo

                  !-- tag point to replace

                  if (pabs .lt. dz**2) then  ! interface too close to k
                     pg = 1. !#
                     do m = -2,3
                        flag_z(i,j,k+m) = -1.*m
                        shift(i,j,k+m) = 0.
                     enddo
                  elseif (pabs .gt. dz-dz**2) then  ! interface too close to k+1
                     do m = -2,3
                        flag_z(i,j,k+m) = -1.*m + 1.
                        shift(i,j,k+m) = 0.
                     enddo
                  else  ! interface in between
                     do m = -2,0
                        flag_z(i,j,k+m) = -1.*m + 1.
                     enddo
                     do m = 1,3
                        flag_z(i,j,k+m) = -1.*m
                     enddo
                  endif
                  
               endif

               if (pijk*p_km .lt. 0. .and. pabs .lt. dy**2) pg = 1. !#
               

               !-- distance for most adjacent node
               
               if ( pg .eq. 1.) then
               
                  d1 = abs(p_ip - p_im)/2.
                  d2 = abs(p_jp - p_jm)/2.
                  d3 = abs(p_kp - p_km)/2.
                  d4 = abs(p_ip - pijk)
                  d5 = abs(p_jp - pijk)
                  d6 = abs(p_kp - pijk)
                  d7 = abs(pijk - p_im)
                  d8 = abs(pijk - p_jm)
                  d9 = abs(pijk - p_km)
               
                  d0 = max( sqrt(d1**2+d2**2+d3**2), sqrt(d4**2+d5**2+d6**2),sqrt(d7**2+d8**2+d9**2),eps )
                  !d0 = max(d1,d2,d3,d4,d5,d6,d7,d8,d9,eps)
                  
                  dist_a(i,j,k) = dz*pijk/d0
               
               endif
 
            endif
 
         enddo
      enddo
   enddo

   call boundsls(dist_x)
   call boundsls(dist_y)
   call boundsls(dist_z)
   call boundsls(flag_x)
   call boundsls(flag_y)
   call boundsls(flag_z)
   call boundsls(shift)

   flag_xo(:,:,:) = flag_x(1:imax,1:jmax,1:kmax)
   flag_yo(:,:,:) = flag_y(1:imax,1:jmax,1:kmax)
   flag_zo(:,:,:) = flag_z(1:imax,1:jmax,1:kmax)
   dist_xo(:,:,:) = dist_x(1:imax,1:jmax,1:kmax)
   dist_yo(:,:,:) = dist_y(1:imax,1:jmax,1:kmax)
   dist_zo(:,:,:) = dist_z(1:imax,1:jmax,1:kmax)
   region(:,:,:) = shift(1:imax,1:jmax,1:kmax) !#
   
 
   return
 end subroutine RussoSmereka_init_4


 subroutine RussoSmereka_iterate_4(sgn,region,dist_a, flag_x,flag_y,flag_z,dist_x,dist_y,dist_z,phi, RHS)

   real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: sgn,region,dist_a
   real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: flag_x,flag_y,flag_z
   real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: dist_x,dist_y,dist_z
   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in)  :: phi

   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(out) :: RHS

   integer :: i,j,k,m

   real :: s, shif,flag,dist, G
   real :: a,b,c,d,e,f
   real :: ap,bp,cp,dp,ep,fp
   real :: am,bm,cm,dm,em,fm
   real :: apbm, cpdm, epfm
   real :: ambp, cmdp, emfp
   real, dimension(0:3) :: ph
   real, dimension(-3:3) :: stcl,func

   RHS = 0.
   
   do k = 1,kmax
      do j = 1,jmax
         do i = 1,imax
            
            s     = sgn(i,j,k)
            ph(0) = phi(i,j,k)

            if (dist_a(i,j,k) .ne. 0.) then  ! (i,j,k) too close to interface
            
               RHS(i,j,k) = -1./dx*(s*abs(ph(0))-dist_a(i,j,k))
            
            else  ! fourth order ENO w/ modified stencil
            
               
               shif = region(i,j,k) !#

               !-- x-dir
               
               a = 0.
               b = 0.

               !-- y-dir

               do m = -3,3
                  stcl(m) = m*dy
                  func(m) = phi(i,j+m,k)
               enddo

               flag = flag_y(i,j,k)
               dist = dist_y(i,j,k)

               if (flag .ne. 0.) call update_stencil(shif,flag,dist,stcl,func)
               !#if (flag .eq. 2. .and. i .eq. 1 .and. k .eq. 1) write(*,*) j,phi(i,j,k),stcl(3),func(3) !!##

               call ENO4_deri(stcl,func, c,d)
            
               !-- z-dir
               
               do m = -3,3
                  stcl(m) = m*dz
                  func(m) = phi(i,j,k+m)
               enddo

               flag = flag_z(i,j,k)
               dist = dist_z(i,j,k)

               if (flag .ne. 0.) call update_stencil(shif,flag,dist,stcl,func)

               call ENO4_deri(stcl,func, e,f)

               !-- Godunov scheme
               
               ap = max(a,0.)
               bp = max(b,0.)
               cp = max(c,0.)
               dp = max(d,0.)
               ep = max(e,0.)
               fp = max(f,0.)

               am = min(a,0.)
               bm = min(b,0.)
               cm = min(c,0.)
               dm = min(d,0.)
               em = min(e,0.)
               fm = min(f,0.)

               apbm = max(ap**2,bm**2)
               cpdm = max(cp**2,dm**2)
               epfm = max(ep**2,fm**2)

               ambp = max(am**2,bp**2)
               cmdp = max(cm**2,dp**2)
               emfp = max(em**2,fp**2)

               if (s .gt. 0.) then
                  G = sqrt(apbm + cpdm + epfm) -1.
               else
                  G = sqrt(ambp + cmdp + emfp) -1.
               endif

               RHS(i,j,k) = -s*G

               !if (shif .eq. 0. .and. abs(ph(0)) .lt. dy**2) then
               !   RHS(i,j,k) = 0.
               !endif

               !if (abs(ph(0)) .lt. NarrowBand_2 .and. abs(G)/dx .gt. 2. .and. i .eq. 1) then
               !   !write(*,*) 'ahhhhh',j,k,ph(0)/dx,c,d,e,f,G/dx  !!!!###
               !   write(*,*) 'bahhhh',j,k,ph(0)/dx,ep,fm,G/dx  !!!!###
               !endif

            endif
               
         enddo
      enddo
   enddo
            

   return
 end subroutine RussoSmereka_iterate_4


 function zero_check(a,b) result(c)

   real, intent(in) :: a,b
   logical :: c

   if (abs(a)/dy .gt. 1e-6 .and. abs(b)/dy .gt. 1e-6) then
   !if (abs(a)/dy .gt. 5.*dy .and. abs(b)/dy .gt. 5.*dy) then
      c = .true.
   else
      c = .false.
   endif

   return
 end function zero_check

 
 subroutine update_stencil(shift,flag,d,s,g)

   real, intent(in) :: shift,flag,d
   real, dimension(-3:3), intent(out) :: s,g


   if (shift .eq. 1.) then
      
      if (flag .eq. 3.) then ! (2,3)
         s(3) = d
         g(3) = 0.
      elseif (flag .eq. 2.) then ! (1,2)
         s(3) = s(2)
         g(3) = g(2)
         s(2) = d
         g(2) = 0.
      elseif (flag .eq. 1.) then ! (0,1)
         s(3) = s(2)
         g(3) = g(2)
         s(2) = s(1)
         g(2) = g(1)
         s(1) = d
         g(1) = 0.
      elseif (flag .eq. -1.) then ! (-1,0)
         s(-3) = s(-2)
         g(-3) = g(-2)
         s(-2) = s(-1)
         g(-2) = g(-1)
         s(-1) = d
         g(-1) = 0.
      elseif (flag .eq. 2.) then ! (-2,-1)
         s(-3) = s(-2)
         g(-3) = g(-2)
         s(-2) = d
         g(-2) = 0.
      elseif (flag .eq. -3.) then ! (-3,-2)
         s(-3) = d
         g(-3) = 0.
      endif

   else

      if (flag .eq. 3.) then ! (2,3)
         s(3) = d
         g(3) = 0.
      elseif (flag .eq. 2.) then ! (1,2)
         s(2) = d
         g(2) = 0.
      elseif (flag .eq. 1.) then ! (0,1)
         s(1) = d
         g(1) = 0.
      elseif (flag .eq. -1.) then ! (-1,0)
         s(-1) = d
         g(-1) = 0.
      elseif (flag .eq. 2.) then ! (-2,-1)
         s(-2) = d
         g(-2) = 0.
      elseif (flag .eq. -3.) then ! (-3,-2)
         s(-3) = d
         g(-3) = 0.
      endif

   endif
   
   return
 end subroutine update_stencil











 

 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 function interpolate_dist(f,h,m)

   real, dimension(0:3), intent(in) :: f
   real, intent(in) :: h,m
   integer :: i
   real :: fi, dist,flag
   real, dimension(2) :: interpolate_dist

   ! remove zero nodes (otherwise get NaN)
   do i = 0,3
      fi = abs(f(i))
      if (fi .lt. 1e-9) then
         interpolate_dist = [0.,0.]
         return
      endif
   enddo

   if (f(0)*f(3) .lt. 0.) then  ! between 0 and 3m
      if (f(0)*f(2) .gt. 0.) then  ! between 2m and 3m
         dist = linear_int(f(2),f(3),h)
         dist = (2.*h + dist)*m
         !call cubic_int(f,dist)
         flag = 1.
      else
         if (f(0)*f(1) .gt. 0.) then  ! between 1m and 2m
            dist = linear_int(f(1),f(2),h)
            dist = (1.*h + dist)*m
            !call cubic_int(f,dist)
            flag = 1.
         else  ! between 0 and 1m
            dist = linear_int(f(0),f(1),h)
            dist = dist*m
            !call cubic_int(f,dist)
            flag = 1.
         endif
      endif
   else  ! no interface between 0 and 3m
      dist = 0.
      flag = 0.
   endif

   interpolate_dist = [dist,flag]
                  
   return
 end function interpolate_dist

 function linear_int(a,b,delta) result(dist)

   real, intent(in) :: a,b,delta
   real :: dist

   dist = delta/(1.+abs(b)/abs(a))
   
   return
 end function linear_int
 

 subroutine cubic_int(a,x0)

   use mod_misc

   real, dimension(0:3), intent(in) :: a
   real, intent(out) :: x0

   integer :: i
   real, dimension(0:3) :: x,b
   real, dimension(0:3,0:3) :: a_mat, a_mat_inv

   do i = 0,3
      a_mat(i,0) = 1. 
      a_mat(i,1) = a(i)
      a_mat(i,2) = a(i)*a_mat(i,1)
      a_mat(i,3) = a(i)*a_mat(i,2)
      b(i) = i*dy
   enddo

   a_mat_inv = inv(a_mat)
   x = matmul(a_mat_inv,b)
   x0 = x(0)
   
   return
 end subroutine cubic_int
 
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 



 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! 
 !  ENO scheme for spatial derivatives

 subroutine ENO2_deri(x_p,f_p, a,b)

   real, dimension(-2:2), intent(in) :: x_p,f_p

   real, intent(out) :: a,b

   integer :: m
   real :: cm,cp

   real, dimension(-2:1) :: Ph1
   real, dimension(-2:0) :: Ph2

   
   !-- divided difference, degree 1

   do m = -2,1
      Ph1(m) = (f_p(m+1)-f_p(m))/(x_p(m+1)-x_p(m))
   enddo

   !-- divided difference, degree 2

   do m = -2,0
      Ph2(m) = (Ph1(m+1)-Ph1(m))/(x_p(m+2)-x_p(m))
   enddo

   cm = MinMod( Ph2(-2), Ph2(-1) )
   cp = MinMod( Ph2(-1), Ph2(0)  )

   a = Ph1(-1) + cm*(x_p(0)-x_p(-1))
   b = Ph1(0)  + cp*(x_p(0)-x_p(1) )
   

   return
 end subroutine ENO2_deri


 subroutine ENO4_deri(x_p,f_p, a,b)

   real, dimension(-3:3), intent(in) :: x_p,f_p

   real, intent(out) :: a,b

   integer :: m
   real :: a1,a2,a3, b1,b2,b3
   real :: ca,cb, d1,d2,d3

   real, dimension(-5:5) :: Ph1
   real, dimension(-2:2) :: Ph2
   real, dimension(-3:3) :: Ph3

   
   !-- divided difference, degree 1

   Ph1 = 0.
   do m = -3,2
      Ph1(2*m+1) = (f_p(m+1)-f_p(m))/(x_p(m+1)-x_p(m))
   enddo

   !-- divided difference, degree 2

   do m = -2,2
      Ph2(m) = (Ph1(2*m+1)-Ph1(2*m-1))/(x_p(m+1)-x_p(m-1))
   enddo

   !-- divided difference, degree 3

   Ph3 = 0.
   do m = -2,1
      Ph3(2*m+1) = (Ph2(m+1)-Ph2(m))/(x_p(m+2)-x_p(m-1))
   enddo

   !-- term 1
   
   a1 = Ph1(-1)
   b1 = Ph1( 1)

   !-- term 2
   
   ca = MinMod( Ph2(-1), Ph2(0) )
   cb = MinMod( Ph2( 0), Ph2(1) )

   a2 = ca*(x_p(0)-x_p(-1))
   b2 = cb*(x_p(0)-x_p( 1))
   
   !-- term 3
   
   d1 = MinAbs( Ph3(-3), Ph3(-1) )
   d2 = MinAbs( Ph3(-1), Ph3( 1) )
   d3 = MinAbs( Ph3( 1), Ph3( 3) )

   if ( abs(Ph2(-1)) .lt. abs(Ph2(0)) ) then
      a3 = (x_p(0)-x_p(-1))*(x_p(0)-x_p(-2))*d1
   else
      a3 = (x_p(0)-x_p(-1))*(x_p(0)-x_p(1) )*d2
   endif

   if ( abs(Ph2(0)) .lt. abs(Ph2(1)) ) then
      b3 = (x_p(0)-x_p(-1))*(x_p(0)-x_p(1) )*d2
   else
      b3 = (x_p(0)-x_p(1) )*(x_p(0)-x_p(2) )*d3
   endif

   a = a1 + a2 + a3
   b = b1 + b2 + b3
   

   return
 end subroutine ENO4_deri


 function MinMod(a,b) result(mm)

   real, intent(in) :: a,b
   real :: a1,b1,ab, mm

   a1 = abs(a)
   b1 = abs(b)
   ab = a*b
   
   if (ab .gt. 0.) then
      if (a1 .le. b1) then
         mm = a
      else
         mm = b
      endif
   else
      mm = 0.
   endif
   
   return
 end function MinMod

 
 function MinAbs(a,b) result(mm)

   real, intent(in) :: a,b
   real :: a1,b1, mm

   a1 = abs(a)
   b1 = abs(b)
   
   if (a1 .lt. b1) then
      mm = a
   else
      mm = b
   endif
   
   return
 end function MinAbs
 
 
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! 
 !  WENO5 stencil selection for the reinitialization
 
 subroutine WENO5_reinit(sgn,phi,RHS)
   
   ! Compute the reinit R.H.S using WENO5, adapted from JC.
   ! Originally he quoted Sussman & Fatemi.            
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   integer :: i, j, k

   real :: dydx,dydy,dydz

   real, dimension(1:imax,1:jmax,1:kmax) :: A, B, C, D, E, F
   real, dimension(1:imax,1:jmax,1:kmax),    intent(in)  :: sgn
   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in)  :: phi
   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(out) :: RHS


   RHS = 0.

   call WENO5_deri(phi, A,B,C,D,E,F)

   do k = 1,kmax
      do j = 1,jmax
         do i = 1,imax

            if (abs(phi(i,j,k)) .lt. NarrowBand_2 + dz) then ! slightly beyond to reshape higher gradients

               !+++++ In the x-direction +++++!

               if (TwoD) then

                  dydx = 0.

               else

                  if ( ( B(i,j,k)*sgn(i,j,k) .LT. 0. ) .AND. &
                       ((A(i,j,k) + B(i,j,k))*sgn(i,j,k) .LT. 0.) ) then

                     dydx = B(i,j,k)

                  elseif (( A(i,j,k)*sgn(i,j,k) .GT. 0. ) .AND. &
                       ((A(i,j,k) + B(i,j,k))*sgn(i,j,k) .GT. 0.) ) then

                     dydx = A(i,j,k)

                  else

                     dydx = 0.5*(A(i,j,k)+B(i,j,k))

                  endif

               endif

               !+++++ In the y-direction +++++!

               if (( D(i,j,k)*sgn(i,j,k) .LT. 0. ) .AND. &
                    ((C(i,j,k) + D(i,j,k))*sgn(i,j,k) .LT. 0.) ) then

                  dydy = D(i,j,k)

               elseif ((C(i,j,k)*sgn(i,j,k) .GT. 0. ) .AND. &
                    ((C(i,j,k) + D(i,j,k))*sgn(i,j,k) .GT. 0.) ) then

                  dydy = C(i,j,k)

               else

                  dydy = 0.5*(C(i,j,k)+D(i,j,k))

               endif

               !+++++ In the z-direction +++++!

               if ( ( F(i,j,k)*sgn(i,j,k) .LT. 0. ) .AND. &
                    ((E(i,j,k) + F(i,j,k))*sgn(i,j,k) .LT. 0.) ) then

                  dydz = F(i,j,k)

               elseif ( (E(i,j,k)*sgn(i,j,k) .GT. 0. ) .AND. &
                    ( (E(i,j,k) + F(i,j,k))*sgn(i,j,k) .GT. 0.) ) then

                  dydz = E(i,j,k)

               else

                  dydz = 0.5*(E(i,j,k)+F(i,j,k))

               endif

               !+++++ Compute right hand side +++++!

               RHS(i,j,k) = sgn(i,j,k) * ( 1.-sqrt(dydx**2 + dydy**2 + dydz**2) )

            endif ! end narrow band

         enddo
      enddo
   enddo


   return
 end subroutine WENO5_reinit


 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! 
 !  WENO5z stencil selection for the reinitialization
 
 subroutine WENO5z_reinit(sgn,phi,RHS)
                                                
   ! Compute the reinit R.H.S using WENO5z (Tanguy, private communication)        
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   integer :: i,j,k, cnt

   real :: fm,fp,sign,smoy
   real :: dydx,dydy,dydz, dfdx,dfdy,dfdz

   real, dimension(1:imax,1:jmax,1:kmax) :: A, B, C, D, E, F
   real, dimension(1:imax,1:jmax,1:kmax),    intent(in)  :: sgn
   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in)  :: phi
   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(out) :: RHS


   call WENO5z_deri(phi, A,B,C,D,E,F)

   do k = 1,kmax
      do j = 1,jmax
         do i = 1,imax
   
            if (abs(phi(i,j,k)) .lt. NarrowBand_2 +1.*dz) then ! slightly beyond to reshape higher gradients
      
               sign = sgn(i,j,k)

               !---- x-direction ----!

               if (TwoD) then

                  dfdx = 0.

               else

                  fm = A(i,j,k)
                  fp = B(i,j,k)

                  if (fm*sign .lt. 0. .and. fp*sign .lt. 0.) then

                     dfdx = fp

                  else

                     if (fm*sign .gt. 0.  .and. fp*sign .gt. 0.) then

                        dfdx = fm

                     else

                        if (fm*sign .lt. 0. .and. fp*sign .gt. 0.) then

                           dfdx = 0.

                        else

                           if (fm*sign .gt. 0. .and. fp*sign .lt. 0.) then

                              smoy = sign*(abs(fp)-abs(fm))/(fp-fm)

                              if (smoy .gt. 0.) then

                                 dfdx = fm

                              else

                                 dfdx = fp

                              endif
                           endif
                        endif
                     endif
                  endif
               endif

               !---- y-direction ----!

               fm = C(i,j,k)
               fp = D(i,j,k)

               if (fm*sign .lt. 0. .and. fp*sign .lt. 0.) then

                  dfdy = fp

               else

                  if (fm*sign .gt. 0.  .and. fp*sign .gt. 0.) then

                     dfdy = fm

                  else

                     if (fm*sign .lt. 0. .and. fp*sign .gt. 0.) then

                        dfdy = 0.

                     else

                        if (fm*sign .gt. 0. .and. fp*sign .lt. 0.) then

                           smoy = sign*(abs(fp)-abs(fm))/(fp-fm)

                           if (smoy .gt. 0.) then

                              dfdy = fm

                           else

                              dfdy = fp

                           endif
                        endif
                     endif
                  endif
               endif

               !---- z-direction ----!

               fm = E(i,j,k)
               fp = F(i,j,k)

               if (fm*sign .lt. 0. .and. fp*sign .lt. 0.) then

                  dfdz = fp

               else

                  if (fm*sign .gt. 0.  .and. fp*sign .gt. 0.) then

                     dfdz = fm

                  else

                     if (fm*sign .lt. 0. .and. fp*sign .gt. 0.) then

                        dfdz = 0.

                     else

                        if (fm*sign .gt. 0. .and. fp*sign .lt. 0.) then

                           smoy = sign*(abs(fp)-abs(fm))/(fp-fm)

                           if (smoy .gt. 0.) then

                              dfdz = fm

                           else

                              dfdz = fp

                           endif
                        endif
                     endif
                  endif
               endif

               !---- Compute the right hand side ----!

               RHS(i,j,k) = sign * ( 1.-sqrt(dfdx**2 + dfdy**2 + dfdz**2) )

            else

               RHS(i,j,k) = 0.
            
            endif ! end narrow band
   
         enddo
      enddo
   enddo


   return
 end subroutine WENO5z_reinit


 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 !  WENO5 for reinitialization
 !  (both derivatives are required)
 
 subroutine WENO5_deri(phi, A,B,C,D,E,F)
 
  real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in) :: phi
  
  !---- 5th order left and right derivatives of phi
  
  real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: A, B, C, D, E, F
  
  !----- Parameters
  
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
 
 
   A = 0.
   B = 0.
   C = 0.
   D = 0.
   E = 0.
   F = 0.
 
   do k = 1,kmax
     do j = 1,jmax
       do i = 1,imax
   
         if (abs(phi(i,j,k)) .lt. NarrowBand_2 +dz) then  ! slightly beyond to reshape higher gradients
   
           !----- In x-direction -----!

           if (.not. TwoD) then
               
           v1 = ( phi(i-2,j,k) - phi(i-3,j,k) ) / dx
           v2 = ( phi(i-1,j,k) - phi(i-2,j,k) ) / dx
           v3 = ( phi(i,j,k)   - phi(i-1,j,k) ) / dx
           v4 = ( phi(i+1,j,k) - phi(i,j,k)   ) / dx
           v5 = ( phi(i+2,j,k) - phi(i+1,j,k) ) / dx
       
           S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
                + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
       
           S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                + (1./4.) * (v2 - v4)**2
       
           S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
                + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
       
           a1 = sigma_1 / (eps_weno + S1)**2
           a2 = sigma_2 / (eps_weno + S2)**2
           a3 = sigma_3 / (eps_weno + S3)**2
       
           w1 = a1 / (a1 + a2 + a3)
           w2 = a2 / (a1 + a2 + a3)
           w3 = a3 / (a1 + a2 + a3)
       
           A(i,j,k) = w1 * ( (1./3.)  * v1   &
                           - (7./6.)  * v2   &
                           + (11./6.) * v3 ) &
                    + w2 * ( (-1./6.) * v2   &
                           + (5./6.)  * v3   &
                           + (1./3.)  * v4 ) &
                    + w3 * ( (1./3.)  * v3   &
                           + (5./6.)  * v4   &
                           - (1./6.)  * v5 )
       
       
           v1 = ( phi(i+3,j,k) - phi(i+2,j,k) ) / dx
           v2 = ( phi(i+2,j,k) - phi(i+1,j,k) ) / dx
           v3 = ( phi(i+1,j,k) - phi(i,j,k)   ) / dx
           v4 = ( phi(i,j,k)   - phi(i-1,j,k) ) / dx
           v5 = ( phi(i-1,j,k) - phi(i-2,j,k) ) / dx
       
       
           S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
                + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
       
           S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                        + (1./4.) * (v2 - v4)**2
       
           S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
                + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
       
           a1 = sigma_1 / (eps_weno + S1)**2
           a2 = sigma_2 / (eps_weno + S2)**2
           a3 = sigma_3 / (eps_weno + S3)**2
       
           w1 = a1 / (a1 + a2 + a3)
           w2 = a2 / (a1 + a2 + a3)
           w3 = a3 / (a1 + a2 + a3)
       
           B(i,j,k) = w1 * ( (1./3.)  * v1   &
                           - (7./6.)  * v2   &
                           + (11./6.) * v3 ) &
                    + w2 * ( (-1./6.) * v2   &
                           + (5./6.)  * v3   &
                           + (1./3.)  * v4 ) &
                    + w3 * ( (1./3.)  * v3   &
                           + (5./6.)  * v4   &
                           - (1./6.)  * v5 )

           endif
        
           !----- In y-direction -----!
       
           v1 = ( phi(i,j-2,k) - phi(i,j-3,k) ) / dx
           v2 = ( phi(i,j-1,k) - phi(i,j-2,k) ) / dx
           v3 = ( phi(i,j,k)   - phi(i,j-1,k) ) / dx
           v4 = ( phi(i,j+1,k) - phi(i,j,k)   ) / dx
           v5 = ( phi(i,j+2,k) - phi(i,j+1,k) ) / dx
       
           S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
                + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
                     
           S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                          + (1./4.) * (v2 - v4)**2
                     
           S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
                + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
                     
           a1 = sigma_1 / (eps_weno + S1)**2
           a2 = sigma_2 / (eps_weno + S2)**2
           a3 = sigma_3 / (eps_weno + S3)**2
                     
           w1 = a1 / (a1 + a2 + a3)
           w2 = a2 / (a1 + a2 + a3)
           w3 = a3 / (a1 + a2 + a3)
                     
           C(i,j,k) = w1 * ( (1./3.)  * v1   &
                           - (7./6.)  * v2   &
                           + (11./6.) * v3 ) &
                    + w2 * ( (-1./6.) * v2   &
                           + (5./6.)  * v3   &
                           + (1./3.)  * v4 ) &
                    + w3 * ( (1./3.)  * v3   &
                           + (5./6.)  * v4   &
                           - (1./6.)  * v5 )
                     
                     
           v1 = ( phi(i,j+3,k) - phi(i,j+2,k) ) / dx
           v2 = ( phi(i,j+2,k) - phi(i,j+1,k) ) / dx
           v3 = ( phi(i,j+1,k) - phi(i,j,k)   ) / dx
           v4 = ( phi(i,j,k)   - phi(i,j-1,k) ) / dx
           v5 = ( phi(i,j-1,k) - phi(i,j-2,k) ) / dx
                     
           S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
                + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
                   
           S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                          + (1./4.) * (v2 - v4)**2
                     
           S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
                + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
                     
           a1 = sigma_1 / (eps_weno + S1)**2
           a2 = sigma_2 / (eps_weno + S2)**2
           a3 = sigma_3 / (eps_weno + S3)**2
                     
           w1 = a1 / (a1 + a2 + a3)
           w2 = a2 / (a1 + a2 + a3)
           w3 = a3 / (a1 + a2 + a3)
           
           D(i,j,k) = w1 * ( (1./3.)  * v1   &
                           - (7./6.)  * v2   &
                           + (11./6.) * v3 ) &
                    + w2 * ( (-1./6.) * v2   &
                           + (5./6.)  * v3   &
                           + (1./3.)  * v4 ) &
                    + w3 * ( (1./3.)  * v3   &
                           + (5./6.)  * v4   &
                           - (1./6.)  * v5 )
           
           !----- In z-direction -----!
           
           if ( ( k .GE. 3 ) .AND. ( k .LE. kmax-2 ) ) then
           
             v1 = ( phi(i,j,k-2) - phi(i,j,k-3) ) / dz
             v2 = ( phi(i,j,k-1) - phi(i,j,k-2) ) / dz
             v3 = ( phi(i,j,k)   - phi(i,j,k-1) ) / dz
             v4 = ( phi(i,j,k+1) - phi(i,j,k)   ) / dz
             v5 = ( phi(i,j,k+2) - phi(i,j,k+1) ) / dz
           
             S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
                  + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
           
             S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                  + (1./4.) * (v2 - v4)**2
           
             S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
                  + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
           
             a1 = sigma_1 / (eps_weno + S1)**2
             a2 = sigma_2 / (eps_weno + S2)**2
             a3 = sigma_3 / (eps_weno + S3)**2
           
             w1 = a1 / (a1 + a2 + a3)
             w2 = a2 / (a1 + a2 + a3)
             w3 = a3 / (a1 + a2 + a3)
           
             E(i,j,k) = w1 * ( (1./3.)  * v1   &
                             - (7./6.)  * v2   &
                             + (11./6.) * v3 ) &
                      + w2 * ( (-1./6.) * v2   &
                             + (5./6.)  * v3   &
                             + (1./3.)  * v4 ) &
                      + w3 * ( (1./3.)  * v3   &
                             + (5./6.)  * v4   &
                             - (1./6.)  * v5 )
           
           
             v1 = ( phi(i,j,k+3) - phi(i,j,k+2) ) / dz
             v2 = ( phi(i,j,k+2) - phi(i,j,k+1) ) / dz
             v3 = ( phi(i,j,k+1) - phi(i,j,k)   ) / dz
             v4 = ( phi(i,j,k)   - phi(i,j,k-1) ) / dz
             v5 = ( phi(i,j,k-1) - phi(i,j,k-2) ) / dz
             
             S1 = (13./12.) * ( v1 - 2.*v2 + v3 )**2 &
                  + (1./4.) * ( v1 - 4.*v2 + 3.*v3)**2
           
             S2 = (13./12.) * (v2 - 2.*v3 + v4)**2 &
                  + (1./4.) * (v2 - v4)**2
           
             S3 = (13./12.) * (v3 - 2.*v4 + v5)**2 &
                  + (1./4.) * (3.*v3 - 4.*v4 + v5)**2
           
             a1 = sigma_1 / (eps_weno + S1)**2
             a2 = sigma_2 / (eps_weno + S2)**2
             a3 = sigma_3 / (eps_weno + S3)**2
           
             w1 = a1 / (a1 + a2 + a3)
             w2 = a2 / (a1 + a2 + a3)
             w3 = a3 / (a1 + a2 + a3)
           
             F(i,j,k) = w1 * ( (1./3.)  * v1   &
                             - (7./6.)  * v2   &
                             + (11./6.) * v3 ) &
                      + w2 * ( (-1./6.) * v2   &
                             + (5./6.)  * v3   &
                             + (1./3.)  * v4 ) &
                      + w3 * ( (1./3.)  * v3   &
                             + (5./6.)  * v4   &
                             - (1./6.)  * v5 )
           
           else  !# 1st order upwinding difference
           
              E(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / dz
              F(i,j,k) = ( phi(i,j,k+1) - phi(i,j,k) ) / dz
              
           endif
   
         endif ! end narrow band
   
       enddo
     enddo
   enddo
 
     
  return
 end subroutine WENO5_deri



 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 !  WENO5z for reinitialization
 !  (both derivatives are required)
 
 subroutine WENO5z_deri(phi, A,B,C,D,E,F)

   real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in) :: phi

   !---- 5th order left and right derivatives of phi

   real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: A, B, C, D, E, F

   !----- Parameters

   real :: eps_weno
   parameter ( eps_weno = 1.e-40 )
   real :: sigma_1, sigma_2, sigma_3
   parameter ( sigma_1 = 0.1 )
   parameter ( sigma_2 = 0.6 )
   parameter ( sigma_3 = 0.3 )
   integer :: i, j, k
   real :: d13p12,d1p4,d1p3,d7p6,d11p6,d1p6,d5p6
   real :: vm1,vm2,vm3,vm4,vm5, vp1,vp2,vp3,vp4,vp5
   real :: S1, S2, S3, tau5
   real :: a1, a2, a3, asum
   real :: w1, w2, w3


   !-- Initialization --!
   
   A = 0.
   B = 0.
   C = 0.
   D = 0.
   E = 0.
   F = 0.

   !-- Pre-computate some real divisions --!

   d13p12 = 13./12.
   d1p4   =  1./4.
   d1p3   =  1./3.
   d7p6   =  7./6.
   d11p6  = 11./6.
   d1p6   =  1./6.
   d5p6   =  5./6.

   !-- WENO5-Z --!
 
   do k = 1,kmax
      do j = 1,jmax
         do i = 1,imax

            if (abs(phi(i,j,k)) .lt. NarrowBand_2 +1.*dz) then  ! slightly beyond to reshape higher gradients

               !----- In x-direction -----!

               if (.not. TwoD) then
       
                  vm1 = ( phi(i-2,j,k) - phi(i-3,j,k) ) / dx
                  vm2 = ( phi(i-1,j,k) - phi(i-2,j,k) ) / dx
                  vm3 = ( phi(i,j,k)   - phi(i-1,j,k) ) / dx
                  vm4 = ( phi(i+1,j,k) - phi(i,j,k)   ) / dx
                  vm5 = ( phi(i+2,j,k) - phi(i+1,j,k) ) / dx

                  S1 = d13p12 * (   vm1 - 2.*vm2 +    vm3)**2 &
                       + d1p4 * (   vm1 - 4.*vm2 + 3.*vm3)**2
                                       
                  S2 = d13p12 * (   vm2 - 2.*vm3 +    vm4)**2 &
                       + d1p4 * (   vm2 -    vm4         )**2
                                                         
                  S3 = d13p12 * (   vm3 - 2.*vm4 +    vm5)**2 &
                       + d1p4 * (3.*vm3 - 4.*vm4 +    vm5)**2

                  tau5 = abs(S1-S3)
                  
                  a1 = sigma_1 *(1. + tau5/(eps_weno + S1))
                  a2 = sigma_2 *(1. + tau5/(eps_weno + S2))
                  a3 = sigma_3 *(1. + tau5/(eps_weno + S3))
                  asum = a1 + a2 + a3

                  w1 = a1 / asum
                  w2 = a2 / asum
                  w3 = a3 / asum

                  A(i,j,k) = w1 * ( d1p3  * vm1   &
                                  - d7p6  * vm2   &
                                  + d11p6 * vm3 ) &
                           + w2 * ( -d1p6 * vm2   &
                                  + d5p6  * vm3   &
                                  + d1p3  * vm4 ) &
                           + w3 * ( d1p3  * vm3   &
                                  + d5p6  * vm4   &
                                  - d1p6  * vm5 )


                  vp1 = ( phi(i+3,j,k) - phi(i+2,j,k) ) / dx
                  vp2 = vm5
                  vp3 = vm4
                  vp4 = vm3
                  vp5 = vm2

                  S1 = d13p12 * (   vp1 - 2.*vp2 +    vp3)**2 &
                       + d1p4 * (   vp1 - 4.*vp2 + 3.*vp3)**2
                                                          
                  S2 = d13p12 * (   vp2 - 2.*vp3 +    vp4)**2 &
                       + d1p4 * (   vp2 -             vp4)**2
                                                          
                  S3 = d13p12 * (   vp3 - 2.*vp4 +    vp5)**2 &
                       + d1p4 * (3.*vp3 - 4.*vp4 +    vp5)**2

                  tau5 = abs(S1-S3)

                  a1 = sigma_1 *(1. + tau5/(eps_weno + S1))
                  a2 = sigma_2 *(1. + tau5/(eps_weno + S2))
                  a3 = sigma_3 *(1. + tau5/(eps_weno + S3))
                  asum = a1 + a2 + a3

                  w1 = a1 / asum
                  w2 = a2 / asum
                  w3 = a3 / asum

                  B(i,j,k) = w1 * ( d1p3  * vp1   &
                                  - d7p6  * vp2   &
                                  + d11p6 * vp3 ) &
                           + w2 * ( -d1p6 * vp2   &
                                  + d5p6  * vp3   &
                                  + d1p3  * vp4 ) &
                           + w3 * ( d1p3  * vp3   &
                                  + d5p6  * vp4   &
                                  - d1p6  * vp5 )

               endif
       
               !----- In y-direction -----!
               
               vm1 = ( phi(i,j-2,k) - phi(i,j-3,k) ) / dy
               vm2 = ( phi(i,j-1,k) - phi(i,j-2,k) ) / dy
               vm3 = ( phi(i,j,  k) - phi(i,j-1,k) ) / dy
               vm4 = ( phi(i,j+1,k) - phi(i,j,  k) ) / dy
               vm5 = ( phi(i,j+2,k) - phi(i,j+1,k) ) / dy

               S1 = d13p12 * (   vm1 - 2.*vm2 +    vm3)**2 &
                    + d1p4 * (   vm1 - 4.*vm2 + 3.*vm3)**2
                                    
               S2 = d13p12 * (   vm2 - 2.*vm3 +    vm4)**2 &
                    + d1p4 * (   vm2 -    vm4         )**2
                                                      
               S3 = d13p12 * (   vm3 - 2.*vm4 +    vm5)**2 &
                    + d1p4 * (3.*vm3 - 4.*vm4 +    vm5)**2

               tau5 = abs(S1-S3)
               
               a1 = sigma_1 *(1. + tau5/(eps_weno + S1))
               a2 = sigma_2 *(1. + tau5/(eps_weno + S2))
               a3 = sigma_3 *(1. + tau5/(eps_weno + S3))
               asum = a1 + a2 + a3

               w1 = a1 / asum
               w2 = a2 / asum
               w3 = a3 / asum

               C(i,j,k) = w1 * ( d1p3  * vm1   &
                               - d7p6  * vm2   &
                               + d11p6 * vm3 ) &
                        + w2 * ( -d1p6 * vm2   &
                               + d5p6  * vm3   &
                               + d1p3  * vm4 ) &
                        + w3 * ( d1p3  * vm3   &
                               + d5p6  * vm4   &
                               - d1p6  * vm5 )


               vp1 = ( phi(i,j+3,k) - phi(i,j+2,k) ) / dy
               vp2 = vm5
               vp3 = vm4
               vp4 = vm3
               vp5 = vm2

               S1 = d13p12 * (   vp1 - 2.*vp2 +    vp3)**2 &
                    + d1p4 * (   vp1 - 4.*vp2 + 3.*vp3)**2
                                                       
               S2 = d13p12 * (   vp2 - 2.*vp3 +    vp4)**2 &
                    + d1p4 * (   vp2 -             vp4)**2
                                                       
               S3 = d13p12 * (   vp3 - 2.*vp4 +    vp5)**2 &
                    + d1p4 * (3.*vp3 - 4.*vp4 +    vp5)**2

               tau5 = abs(S1-S3)

               a1 = sigma_1 *(1. + tau5/(eps_weno + S1))
               a2 = sigma_2 *(1. + tau5/(eps_weno + S2))
               a3 = sigma_3 *(1. + tau5/(eps_weno + S3))
               asum = a1 + a2 + a3

               w1 = a1 / asum
               w2 = a2 / asum
               w3 = a3 / asum

               D(i,j,k) = w1 * ( d1p3  * vp1   &
                               - d7p6  * vp2   &
                               + d11p6 * vp3 ) &
                        + w2 * ( -d1p6 * vp2   &
                               + d5p6  * vp3   &
                               + d1p3  * vp4 ) &
                        + w3 * ( d1p3  * vp3   &
                               + d5p6  * vp4   &
                               - d1p6  * vp5 )

               
               !----- In z-direction -----!
               
               if ( ( k .GE. 3 ) .and. ( k .LE. kmax-2 )) then

                  vm1 = ( phi(i,j,k-2) - phi(i,j,k-3) ) / dz
                  vm2 = ( phi(i,j,k-1) - phi(i,j,k-2) ) / dz
                  vm3 = ( phi(i,j,k)   - phi(i,j,k-1) ) / dz
                  vm4 = ( phi(i,j,k+1) - phi(i,j,k)   ) / dz
                  vm5 = ( phi(i,j,k+2) - phi(i,j,k+1) ) / dz

                  S1 = d13p12 * (   vm1 - 2.*vm2 +    vm3)**2 &
                       + d1p4 * (   vm1 - 4.*vm2 + 3.*vm3)**2
                                       
                  S2 = d13p12 * (   vm2 - 2.*vm3 +    vm4)**2 &
                       + d1p4 * (   vm2 -    vm4         )**2
                                                         
                  S3 = d13p12 * (   vm3 - 2.*vm4 +    vm5)**2 &
                       + d1p4 * (3.*vm3 - 4.*vm4 +    vm5)**2

                  tau5 = abs(S1-S3)
                  
                  a1 = sigma_1 *(1. + tau5/(eps_weno + S1))
                  a2 = sigma_2 *(1. + tau5/(eps_weno + S2))
                  a3 = sigma_3 *(1. + tau5/(eps_weno + S3))
                  asum = a1 + a2 + a3

                  w1 = a1 / asum
                  w2 = a2 / asum
                  w3 = a3 / asum

                  E(i,j,k) = w1 * ( d1p3  * vm1   &
                                  - d7p6  * vm2   &
                                  + d11p6 * vm3 ) &
                           + w2 * ( -d1p6 * vm2   &
                                  + d5p6  * vm3   &
                                  + d1p3  * vm4 ) &
                           + w3 * ( d1p3  * vm3   &
                                  + d5p6  * vm4   &
                                  - d1p6  * vm5 )


                  vp1 = ( phi(i,j,k+3) - phi(i,j,k+2) ) / dx
                  vp2 = vm5
                  vp3 = vm4
                  vp4 = vm3
                  vp5 = vm2

                  S1 = d13p12 * (   vp1 - 2.*vp2 +    vp3)**2 &
                       + d1p4 * (   vp1 - 4.*vp2 + 3.*vp3)**2
                                                          
                  S2 = d13p12 * (   vp2 - 2.*vp3 +    vp4)**2 &
                       + d1p4 * (   vp2 -             vp4)**2
                                                          
                  S3 = d13p12 * (   vp3 - 2.*vp4 +    vp5)**2 &
                       + d1p4 * (3.*vp3 - 4.*vp4 +    vp5)**2

                  tau5 = abs(S1-S3)

                  a1 = sigma_1 *(1. + tau5/(eps_weno + S1))
                  a2 = sigma_2 *(1. + tau5/(eps_weno + S2))
                  a3 = sigma_3 *(1. + tau5/(eps_weno + S3))
                  asum = a1 + a2 + a3

                  w1 = a1 / asum
                  w2 = a2 / asum
                  w3 = a3 / asum

                  F(i,j,k) = w1 * ( d1p3  * vp1   &
                                  - d7p6  * vp2   &
                                  + d11p6 * vp3 ) &
                           + w2 * ( -d1p6 * vp2   &
                                  + d5p6  * vp3   &
                                  + d1p3  * vp4 ) &
                           + w3 * ( d1p3  * vp3   &
                                  + d5p6  * vp4   &
                                  - d1p6  * vp5 )

               else  !# 1st order upwinding difference

                  E(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / dz
                  F(i,j,k) = ( phi(i,j,k+1) - phi(i,j,k) ) / dz

               endif

            endif ! end narrow band

         enddo
      enddo
   enddo


   return
 end subroutine WENO5z_deri



end module mod_hyperPDE
