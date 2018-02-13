module mod_ab

  use mod_mom
  use mod_common
  use mod_common_mpi
  use mod_param
  use mod_interface
  use mod_bound
  use mod_rib

  implicit none

  private
  public ab2, surf_ten

contains


  !-- RHS of NS
  !
  subroutine ab2(duold,dvold,dwold,pold)
    !                                                       
    ! Update velocity using 2nd order Adams-Bashforth and   
    ! a pressure-corretion method (Dodd & Ferrante JCP2014).
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    integer i,j,k,l

    real, dimension(1:imax,1:jmax,1:kmax) :: dunew,dvnew,dwnew
    real, dimension(0:i1,0:j1,0:k1) :: phat,phat_jx,phat_jy,phat_jz
    real, dimension(0:i1,0:j1,0:k1) :: curvature,curv_temp

    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) ::  y_temp

    real, dimension(1:,1:,1:), intent(inout) :: duold,dvold,dwold
    real, dimension(0:i1,0:j1,0:k1), intent(in) :: pold

  
    !- #1  Calculate u*
    ! 
    !      u*-unew         n         n-1            
    !     -------- = 1.5*RU  - 0.5*RU    
    !        dt
    ! 
    !     where RU  = - convection + viscous diffusion + buoyancy.
    
    if (Riblet) call mask_rib_vel(unew,vnew,wnew)  ! mask the previous velocity
    
    call mom_main(dunew,dvnew,dwnew)

    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             
             dudt(i,j,k) = unew(i,j,k) + dt*(1.5*dunew(i,j,k) -0.5*duold(i,j,k))
             dvdt(i,j,k) = vnew(i,j,k) + dt*(1.5*dvnew(i,j,k) -0.5*dvold(i,j,k))
             dwdt(i,j,k) = wnew(i,j,k) + dt*(1.5*dwnew(i,j,k) -0.5*dwold(i,j,k))
             
          enddo
       enddo
    enddo

    duold = dunew
    dvold = dvnew
    dwold = dwnew

  
    !- #2 Surface tension
    
    select case (surface_tension_method)
    case('CSF')

       call cont_surf_force  ! A source term in NS (u* -> u**)
       
    case('GFM')
       
       ! Jump to be removed from pressure gradient
       
       phat_jx = 2.*p_x - p_xold
       phat_jy = 2.*p_y - p_yold
       phat_jz = 2.*p_z - p_zold

    end select

  
    !- #3 Calculate u*** due to pressure correction (not AB2 but 2nd order in time)
    ! 
    !      u***-u**         1        1         ~     
    !     ---------- = - (----- - ------)*grad p 
    !         dt           rho     rho0 
    !                                           
    !     where rho is local density at time level n+1, 
    ! 
    !     ~      n    n-1 
    !     p = 2*p  - p    is the predicted pressure at time level n+1.
    
    phat = 2.*pnew - pold

    call momxp(dunew,phat,phat_jx,rho_u)  
    call momyp(dvnew,phat,phat_jy,rho_v)
    call momzp(dwnew,phat,phat_jz,rho_w)

    select case(BC_in_y)
    case('Periodic')
       
       if (constPoi) then
          forceytot = -12.*visc/lz**2  ! constant pressure gradient for single-phase channel flow
          !forceytot = -1.*(v_bulk_desired-v_bulk)/dt  ! brutally enforce v_bulk = 1
       else
          forceytot = 0.
       end if
       
    case('InOut','Walls','OutOut')
       
       forceytot = 0.
       
    end select
       
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             
             dudt(i,j,k) = dudt(i,j,k) + dt*dunew(i,j,k)
             dvdt(i,j,k) = dvdt(i,j,k) + dt*dvnew(i,j,k) - forceytot*dt/rho_v(i,j,k)
             dwdt(i,j,k) = dwdt(i,j,k) + dt*dwnew(i,j,k)
             
          enddo
       enddo
    enddo

    !- #  Boundary condition
    
    call bounduvw(dudt,dvdt,dwdt)
    if (Riblet) call mask_rib_vel(dudt,dvdt,dwdt)  ! ensure BC for the prediction velocity
    

    return
  end subroutine ab2


  

  !-- Continuum surface force (Brackbill 1992)
  !
  subroutine cont_surf_force

    integer :: i,j,k,l
    real, dimension(0:i1,0:j1,0:k1) :: curvature,curv_temp
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) ::  y_temp
    
    
    !-- Obtain the entire curvature field
    
    do l=1,lmax

       y_temp = lvset(:,:,:,l)
       call get_curvature(y_temp,curv_temp)
       call boundc(curv_temp)

       !# Caution that the interface cannot be closer than ~3dx with the current scheme. 
       if (l .eq. 1) then
          curvature = curv_temp
       else
          where (curvature .eq. 0.) curvature = curv_temp
       endif

    enddo

    !call exact_curvature(curvature)  !# parasitic currents test

    !---- Calculate u**
    ! 
    !      u**-u*         n         n-1             
    !     ------- = 1.5*RS  - 0.5*RS    
    !        dt
    ! 
    !     where RS = explicit surface tension force.

    call surf_ten(curvature,surfxnew,surfynew,surfznew)

    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             dudt(i,j,k) = dudt(i,j,k) + dt*(1.5*surfxnew(i,j,k) -0.5*surfxold(i,j,k))
             dvdt(i,j,k) = dvdt(i,j,k) + dt*(1.5*surfynew(i,j,k) -0.5*surfyold(i,j,k))
             dwdt(i,j,k) = dwdt(i,j,k) + dt*(1.5*surfznew(i,j,k) -0.5*surfzold(i,j,k))
          enddo
       enddo
    enddo

    surfxold = surfxnew
    surfyold = surfynew
    surfzold = surfznew

       
    return
  end subroutine cont_surf_force


  

  !-- Surface tension force
  !  
  subroutine surf_ten(curvature,surfx,surfy,surfz)
    !                                                   
    ! Compute surface tension force (Brackbill JCP1992):
    !                                                   
    !             1       curv * grad Heaviside         
    !     f = ---------- -----------------------        
    !          We or Ca      (rho1 + rho2)/2            
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    integer i,j,k
    real, dimension(0:i1,0:j1,0:k1), intent(in) :: curvature
    real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: surfx, surfy, surfz


    do k=1,kmax
       do j=1,jmax
          do i=1,imax

             surfx(i,j,k) = 1./rho_u(i,j,k) *(curvature(i+1,j,k)+curvature(i,j,k))/2. &
                  *(jiemian(i+1,j,k)-jiemian(i,j,k))/dx
             surfy(i,j,k) = 1./rho_v(i,j,k) *(curvature(i,j+1,k)+curvature(i,j,k))/2. &
                  *(jiemian(i,j+1,k)-jiemian(i,j,k))/dy
             surfz(i,j,k) = 1./rho_w(i,j,k) *(curvature(i,j,k+1)+curvature(i,j,k))/2. &
                  *(jiemian(i,j,k+1)-jiemian(i,j,k))/dz

          enddo
       enddo
    enddo
    
    surfx = surfx /We
    surfy = surfy /We
    surfz = surfz /We


    return
  end subroutine surf_ten




end module mod_ab
