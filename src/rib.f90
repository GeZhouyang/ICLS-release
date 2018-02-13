module mod_common_rib

  use mod_param

  implicit none

  !-- Rib parameters (in terms of # of the cells) --!
  !
  ! The rib bound has to be 2 cells away from the proc bound (due to mask_rib_vel)
  
  integer, parameter :: l_R = jtot/2  ! length
  integer, parameter :: d_R = jtot/2  ! distance (or gap)
  integer, parameter :: h_R = ktot/3  ! height
  integer, parameter :: n_R = 1       ! number

  integer, dimension(n_R) :: j_l,j_r  ! edge index of the rib
  real, dimension(n_R) :: rib_l,rib_r ! rib indicator

  !-- Streamwise index/mask arrays --!
  
  real, dimension(0:j1)    :: j_uind, j_vind  ! u(w), v grids
  real, dimension(-2:j1+2) :: j_sind          ! scalar grids

  real, dimension(1:jmax) :: s_left_0, s_right_0  ! boundary points
  real, dimension(1:jmax) :: s_left_1, s_right_1  ! first inner points
  real, dimension(1:jmax) :: s_left_2, s_right_2  ! second inner points

  !-- Rib layers --!

  real, dimension(1:imax,1:jmax,1:kmax) :: riblay_v, riblay_w
  
end module mod_common_rib


!------------------------------------------------------------
!------------------------------------------------------------


module mod_rib

  use mod_param
  use mod_param_ls
  use mod_common_rib

  implicit none

  private
  public init_rib, mask_rib_vel, ls_main_rib, HOUC5_nb_rib,RK_nb_rib

contains

  subroutine init_rib

    use mod_common_mpi
    use mod_bound

    integer :: i,j, j_left,j_right, jv,k, m, n
    real :: j_um,j_ui,j_up, j_vm,j_vi,j_vp, j_sl,j_sm,j_si,j_sp,j_sq
    
    !-- Streamwise index arrays of the rib --!
    
    j_uind = 0.  ! 0 denotes fluid
    j_vind = 0.  ! 0 denotes fluid
    j_sind = 0.  ! 0 denotes fluid

    do n = 0,n_R
       j_left =  n*(l_R+d_R)-l_R/2
       j_right = n*(l_R+d_R)+l_R/2
       do j = 0,j1        
          jv = mod(myid,dims(2))*jmax + j
          if (jv .ge. j_left .and. jv .le. j_right) then  ! global location determines local array
             j_vind(j) = 1. ! staggered from cell center
             if (jv .gt. j_left) j_uind(j) = 1.
          endif        
       enddo
    enddo

    do n = 0,n_R
       j_left =  n*(l_R+d_R)-l_R/2
       j_right = n*(l_R+d_R)+l_R/2
       do j = -2,j1+2  ! wider
          jv = mod(myid,dims(2))*jmax + j
          if (jv .gt. j_left .and. jv .le. j_right) then  ! global location determines local array
             j_sind(j) = 1.
          endif
       enddo
    enddo

    !-- Edge index of the rib --!  ## parallelization limit: max 1 edge per proc ##

    j_l(:) = 0  ! local index for the left edge of the rib
    j_r(:) = 0  !                ...  right ...
    rib_l = 0.  ! non-zero if the current proc has a left edge
    rib_r = 0.  !                               ...  right ...
    
    do n = 1,n_R
       do j = 1,jmax
          if (j_sind(j) .eq. 1. .and. j_sind(j-1) .eq. 0.) then
             j_l(n) = j
             rib_l(n) = 1.
          endif
          if (j_sind(j) .eq. 1. .and. j_sind(j+1) .eq. 0.) then
             j_r(n) = j
             rib_r(n) = 1.
          endif
       enddo
    enddo

    ! equalize all procs (local info contained in rib_l and rib_r)
    
    call mpi_allreduce(mpi_in_place, j_l, n_R, mpi_real8, mpi_max, comm_cart,error)
    call mpi_allreduce(mpi_in_place, j_r, n_R, mpi_real8, mpi_max, comm_cart,error)

    !-- Left/Right mask arrays --!

    s_left_0 = 0.
    s_left_1 = 0.
    s_left_2 = 0.

    s_right_0 = 0.
    s_right_1 = 0.
    s_right_2 = 0.

    do j = 1,jmax

       j_sl = j_sind(j-2)
       j_sm = j_sind(j-1)
       j_si = j_sind(j)
       j_sp = j_sind(j+1)
       j_sq = j_sind(j+2)

       if (j_si .eq. 1. .and. j_sp .eq. 0.) s_left_0(j) = 1.
       if (j_sm .eq. 1. .and. j_si .eq. 0.) s_left_1(j) = 1.
       if (j_sl .eq. 1. .and. j_sm .eq. 0.) s_left_2(j) = 1.
       
       if (j_si .eq. 1. .and. j_sm .eq. 0.) s_right_0(j) = 1.
       if (j_sp .eq. 1. .and. j_si .eq. 0.) s_right_1(j) = 1.
       if (j_sq .eq. 1. .and. j_sp .eq. 0.) s_right_2(j) = 1.
       
    enddo

    !-- Layers surounding the rib --!
    
    riblay_v = 0.
    riblay_w = 0.

    do k = 1,kmax
       do j = 1,jmax

          j_um = j_uind(j-1)
          j_ui = j_uind(j)
          j_up = j_uind(j+1)
          j_vm = j_vind(j-1)
          j_vi = j_vind(j)
          j_vp = j_vind(j+1)
          
          do i = 1,imax
             
             if (k .le. h_R) then

                if (k .eq. h_R .and. j_vm*j_vi*j_vp .eq. 1.)          riblay_v(i,j,k) = 1.
                if (j_ui .eq. 1. .and. j_up .eq. 0.)                  riblay_v(i,j,k) = 2.
                if (j_ui .eq. 1. .and. j_um .eq. 0.)                  riblay_v(i,j,k) = 3.
                if (k .lt. h_R .and. j_ui .eq. 0. .and. j_up .eq. 1.) riblay_v(i,j,k) = 4.
                if (k .lt. h_R .and. j_ui .eq. 1. .and. j_up .eq. 0.) riblay_v(i,j,k) = 5.

                if (k .eq. h_R .and. j_ui .eq. 1.)                    riblay_w(i,j,k) = 1.
                if (k .lt. h_R .and. j_ui .eq. 1. .and. j_up .eq. 0.) riblay_w(i,j,k) = 2.
                if (k .lt. h_R .and. j_ui .eq. 1. .and. j_um .eq. 0.) riblay_w(i,j,k) = 3.
                if (k .eq. h_R .and. j_um*j_ui*j_up .eq. 1.)          riblay_w(i,j,k) = 4.

             endif
          enddo
       enddo
    enddo
    
    return
  end subroutine init_rib
  
  
  !------------------------------------------------------------
  
  
  subroutine mask_rib_vel(u,v,w)

    ! Enforce the exact boundary condition for velocity on the rib surface. (Stress IBM.)
    ! The staggered velocity inside is also specified, such that the *approximate* 
    ! boundary condition is fulfilled after the projection.

    use mpi
    use mod_common_mpi
    use mod_bound

    integer :: i,j,k, jv
    real, dimension(0:,0:,0:) :: u,v,w
       
    do k = 0,k1
       do j = 1,jmax

          jv = mod(myid,dims(2))*jmax + j

          do i = 1,imax

             if (k .le. h_R) then  ! underneath
                
                if (riblay_v(i,j,k) .eq. 1.) v(i,j,k)   = -v(i,j,k+1)
                if (riblay_v(i,j,k) .eq. 2.) v(i,j,k)   = 0.
                if (riblay_v(i,j,k) .eq. 3.) v(i,j-1,k) = 0.
                if (riblay_v(i,j,k) .eq. 4.) v(i,j+1,k) = -v(i,j-1,k)
                if (riblay_v(i,j,k) .eq. 5.) v(i,j-1,k) = -v(i,j+1,k)

                if (riblay_w(i,j,k) .eq. 1.) w(i,j,k)   = 0.
                if (riblay_w(i,j,k) .eq. 2.) w(i,j,k)   = -w(i,j+1,k)
                if (riblay_w(i,j,k) .eq. 3.) w(i,j,k)   = -w(i,j-1,k)
                if (riblay_w(i,j,k) .eq. 4.) w(i,j,k-1) = -w(i,j,k+1)
                
             endif

          enddo
       enddo
    enddo

    ! pass x halo data

    call updthalos(u,1)
    call updthalos(v,1)
    call updthalos(w,1)

    ! pass y halo data

    call updthalos(u,2)
    call updthalos(v,2)
    call updthalos(w,2)


    return
  end subroutine mask_rib_vel
  
  
  !------------------------------------------------------------
  
  
  subroutine ls_main_rib(istep)

    ! Main level set subroutine tailored to riblet.

    use mod_common
    use mod_common_mpi
    use mod_misc
    use mod_contact_line_common

    integer, intent(in) :: istep
    integer :: l, iter

    real :: massloss, residual
    
    real, dimension(num_cp) :: app_cont_angle
    real, dimension(1:imax,1:jmax,1:kmax) :: dfdt
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phi,dphi, phi_1,phi_2,phi_3

    !-- Modify the cell-center velocity
    
    call mask_rib_vel(unew,vnew,wnew)
    call average_vel_at_center

    do l = 1,lmax

       if (Contact_line) then
          call contact_line_main_rib(istep, app_cont_angle) !## works only for one level set
       endif

       phi = lvset(:,:,:,l)
       
       !-- Construct the tube within the ribs

       call tube_construction_rib(phi)

       !-- Narrow-band Advection

       call HOUC5_nb_rib(uccnew,vccnew,wccnew, phi, dfdt)
       call RK_nb_rib(phi_1, phi, 1.,dfdt)

       call HOUC5_nb_rib(uccnew,vccnew,wccnew, phi_1, dfdt)
       call RK_nb_rib(phi_2, 3./4.*phi + 1./4.*phi_1, 1./4.,dfdt)

       call HOUC5_nb_rib(uccnew,vccnew,wccnew, phi_2, dfdt)
       call RK_nb_rib(phi_3, 1./3.*phi + 2./3.*phi_2, 2./3.,dfdt)

       lvset(:,:,:,l) = phi_3
       
       !-- Mass correction

       if (mass_correction) then
          maxchk_glb = 0.
          if (mod(istep,inflate_step) .eq. 0) then

             massloss = (volume(l) - init_volume(l))/init_volume(l)
             if (abs(massloss) .gt. 1e-5) then  ! a threshold value

                if (myid .eq. 0) write(*,'(A,I3)') '!-- Inflate interface',l
                call inflation_rib(l)

             endif

          endif
       endif       

       !-- Re-initialization

       if (re_initialization) then
          if (mod(istep,reinit_step) .eq. 0) then
             
             call ls_reinit_rib(l,1,iter,residual)
             call boundsls_rib(lvset(:,:,:,l))
             if (myid .eq. 0) then
                write(*,'(A,I3,A,I5,A,ES10.4)') '!-- Reinit LS',l, &
                     ' for',iter, ' iterations w/ max residual = ',residual
             endif
             
          endif
       endif

       !---- Contact angle ----!
       
       if (Contact_line) call enforce_contact_angle_rib(l,app_cont_angle)

    enddo

    !-- Update Viscosity, Density, and Volume

    call update_materials_rib
    

    return
  end subroutine ls_main_rib


  !------------------------------------------------------------
  

  subroutine contact_line_main_rib(istep, app_cont_angle)

    use mod_contact_line
    use mod_contact_line_common

    ! Compute the contact point velocity:  vel_c
    ! from the apparent contact angle:     app_cont_angle
    ! using a micro model (tabulated results of a separate phase-field simulation).
    ! Then spread the slip velocity around the contact point using a bump funtion.

    integer, intent(in) :: istep
    real, dimension(num_cp), intent(out) :: app_cont_angle

    real, dimension(num_cp) :: cont_loc, vel_c
    real, dimension(ibump,num_cp) :: vel_t,vel_n
    

    !-- Contact point locations & Apparent contact angles --!
    
    call xp_ls_comp_angle_rib(istep, cont_loc,app_cont_angle)

    !-- Contact point velocity --!
    
    call comp_cp_vel(app_cont_angle,vel_c)

    !-- Spreaded slip velocity --!
    
    call spread_slip_vel(istep, app_cont_angle,vel_c, vel_t,vel_n)
    call mask_adv_vel_rib(vel_t,vel_n)

    !-- Output --!

    if (mod(istep,100) .eq. 0) then
       call output_contact(istep,cont_loc,app_cont_angle,vel_c)
    endif

    return
  end subroutine contact_line_main_rib


  !------------------------------------------------------------
  

  subroutine xp_ls_comp_angle_rib(istep, cont_loc, cont_ang)
    
    ! Extrapolate the level set to the wall and
    ! compute the contact point locations and apparent contact angles
    ! ## only for y-normal ##

    use mod_common
    use mod_common_mpi
    use mod_contact_line_common

    integer, intent(in) :: istep
    real, dimension(num_cp), intent(out) :: cont_loc, cont_ang

    integer, parameter :: xp_method = 1
    integer :: i,j,k, jp, jl,jr
    real :: ribl,ribr
    real :: ph,dphdz,d2phdz2, dphdy,d2phdy2
    real :: lowerbound,coory,coorz, d1,d2,d3,d4, e0,e1,e2, theta
    real, dimension(10) :: phi_v
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phi

    i = imax/2  ! 2D
    lowerbound = coords(2)*ly/(dims(2)*1.)
    
    phi = lvset(:,:,:,1)

    d2c_xp1 = 0.
    d2c_xp2 = 0.
    
    cont_loc = 0.
    cont_ang = 0.
    
    !-- Extrapolate level set to contact surfaces

    select case (cont_surf_normal)  
    case('y')

       jr = j_r(1)      ! ## only for 1 rib now
       jl = j_l(1)      ! ## ...
       ribr = rib_r(1)  ! ## ...
       ribl = rib_l(1)  ! ## ...

       select case (xp_method)
       case(1)  ! linear upwind

          do k = -2,k1+2
             d2c_xp1(k) = ribr*( 1.5*phi(i,jr+1,k) - 0.5*phi(i,jr+2,k) )  ! right edge of the rib
             d2c_xp2(k) = ribl*( 1.5*phi(i,jl-1,k) - 0.5*phi(i,jl-2,k) )  ! left ...
          enddo

       case(2)  ! least squares upwind ## Need to modify the pseudo inverse later ##
          
       end select

       ! equalize info in all procs
       
       call mpi_allreduce(mpi_in_place, d2c_xp1,k1+1, mpi_real8, mpi_sum, comm_cart,error)
       call mpi_allreduce(mpi_in_place, d2c_xp2,k1+1, mpi_real8, mpi_sum, comm_cart,error)

       !-- Interpolate the interface location on contact surfaces
       !   (now works only if one contact point per surface)

       !-- Right edge
       
       e1 = 0.
       e2 = 0.

       do k = 1,kmax-1

          coorz = (k-0.5)*dz

          ! adjacent to the edge

          d1 = d2c_xp1(k)
          d2 = d2c_xp1(k+1)
          if (d1*d2 .lt. 0.) then
             e0 = dz/( 1. + abs(d2/d1) )
             e1 = coorz + e0
          elseif (d1 .eq. 0.) then
             e1 = coorz
          endif

          ! half layer into the fluid

          d3 = phi(i,jr+1,k)
          d4 = phi(i,jr+1,k+1)
          if (d3*d4 .lt. 0.) then
             e0 = dz/( 1. + abs(d4/d3) )
             e2 = coorz + e0
          elseif (d3 .eq. 0.) then
             e2 = coorz
          endif

       enddo

       !-- Compute the contact angle (assume flat interface)

       if (e1 .ne. e2) then
          theta = atan( dy/2./(e2-e1) )  ! measured from fluid 1 (phi > 0) 
       else
          theta = pi/2.
       endif
       if (theta .lt. 0.) theta = pi + theta

       !-- Contact point location & contact angle (local proc)

       cont_loc(1) = ribr*e1
       cont_ang(1) = ribr*theta

       !-- Left edge

       e1 = 0.
       e2 = 0.

       do k = 1,kmax-1

          coorz = (k-0.5)*dz

          ! adjacent to the edge

          d1 = d2c_xp2(k)
          d2 = d2c_xp2(k+1)
          if (d1*d2 .lt. 0.) then
             e0 = dz/( 1. + abs(d2/d1) )
             e1 = coorz + e0
          elseif (d1 .eq. 0.) then
             e1 = coorz
          endif

          ! half layer into the fluid

          d3 = phi(i,jl-1,k)
          d4 = phi(i,jl-1,k+1)
          if (d3*d4 .lt. 0.) then
             e0 = dz/( 1. + abs(d4/d3) )
             e2 = coorz + e0
          elseif (d3 .eq. 0.) then
             e2 = coorz
          endif

       enddo

       !-- Compute the contact angle (assume flat interface)

       if (e1 .ne. e2) then
          theta = atan( dy/2./(e2-e1) )  ! measured from fluid 1 (phi > 0) 
       else
          theta = pi/2.
       endif
       if (theta .lt. 0.) theta = pi + theta

       !-- Contact point location & contact angle (local proc)

       cont_loc(2) = ribl*e1
       cont_ang(2) = ribl*theta

       ! equalize info in all procs
       
       call mpi_allreduce(mpi_in_place, cont_loc,num_cp, mpi_real8, mpi_sum, comm_cart,error)
       call mpi_allreduce(mpi_in_place, cont_ang,num_cp, mpi_real8, mpi_sum, comm_cart,error)

       !write(*,*) '////////////', myid,cont_loc(1),cont_ang(1)
       !write(*,*) '\\\\\\\\\\\\', myid,cont_loc(2),cont_ang(2)

    end select
    

    return
  end subroutine xp_ls_comp_angle_rib


  !------------------------------------------------------------------------------------------


  subroutine mask_adv_vel_rib(vel_t,vel_n)

    use mod_common
    use mod_contact_line_common

    !  Update advection velocity within the narrow band

    real, dimension(ibump,num_cp), intent(in) :: vel_t,vel_n
    
    integer :: i,j,k
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phi

    
    select case (cont_surf_normal)
    case('y')
    
       phi = lvset(:,:,:,1)

       j = j_r(1)+1  ! adjacent to the right edge
       do k = 1,kmax
          do i = 1,imax
             if (abs(phi(i,j,k)) .le. NarrowBand_2) then  !# maybe I should do v=+v?

                wccnew(i,j,k) = vel_t(j,1)
                vccnew(i,j,k) = vel_n(j,1)

             endif
          enddo
       enddo
       
       j = j_l(1)-1  ! adjacent to the left edge
       do k = 1,kmax
          do i = 1,imax
             if (abs(phi(i,j,k)) .le. NarrowBand_2) then  !# maybe I should do v=+v?

                wccnew(i,j,k) = vel_t(j,2)
                vccnew(i,j,k) = -vel_n(j,2)  !## mirror

             endif
          enddo
       enddo

    end select
    
    
    return
  end subroutine mask_adv_vel_rib
  

  !------------------------------------------------------------------------------------------


  subroutine enforce_contact_angle_rib(label,app_cont_angle)

    use mod_common
    use mod_hyperPDE
    use mod_contact_line_common

    integer, intent(in) :: label
    real, dimension(num_cp), intent(in) :: app_cont_angle

    integer :: i,j,k, jr,jl,jm,jn
    real :: a,b,c, ribr,ribl, theta, d0,d1,d2
    real, dimension(-2:i1+2,-2:2,-2:k1+2) :: phi
    
    ! The code below is for the second simplest case where the walls are
    ! vertical, ie. nw = (0,+/-1,0), such that the problem reduces to 1D
    ! but parallelization needs to be considered.

    ! Coefficients for the second-order upwind scheme (positive speed)
    ! dfdy(j) = ( a*f(j)+b*f(jm)+c*f(jn) )/dz.
    ! Only one set of coefs is needed due to symmetry.

    select case (cont_surf_normal)
    case('y')
       
       a = 1.5
       b = -2.
       c = 0.5

       ribr = rib_r(1)
       ribl = rib_l(1)
       
       !-- Right edge --!

       if (ribr .eq. 1.) then

          theta = app_cont_angle(1)
          d0 = cos(theta)
       
          jr = j_r(1)
          phi = lvset(:,jr-2:jr+2,:,label)
       
          do k = -2,k1+2
             do j = 0,-2,-1
                jm = j + 1 
                jn = j + 2
                do i = -2,i1+2
       
                   if (abs(phi(i,jn,k)) .lt. NarrowBand_2) then
       
                      d1 = 1./(1./dy*a)
                      d2 = -1./dy*( b*phi(i,jm,k)+c*phi(i,jn,k) ) + d0
                      phi(i,j,k) = d1*d2
       
                   endif
       
                enddo
             enddo
          enddo
       
          lvset(:,jr-2:jr,:,label) = phi(:,-2:0,:)

       endif
       
       !-- Left edge --!

       if (ribl .eq. 1.) then

          theta = app_cont_angle(2)
          d0 = cos(theta)
       
          jl = j_l(1)
          phi = lvset(:,jl-2:jl+2,:,label)
       
          do k = -2,k1+2
             do j = 0,2
                jm = j - 1
                jn = j - 2
                do i = -2,i1+2
       
                   if (abs(phi(i,jn,k)) .lt. NarrowBand_2) then
       
                      d1 = 1./(1./dy*a)
                      d2 = -1./dy*( b*phi(i,jm,k)+c*phi(i,jn,k) ) + d0
                      phi(i,j,k) = d1*d2
       
                   endif
       
                enddo
             enddo
          enddo
       
          lvset(:,jl:jl+2,:,label) = phi(:,0:2,:)

       endif

    end select
    

    return
  end subroutine enforce_contact_angle_rib


  !------------------------------------------------------------------------------------------


  subroutine tube_construction_rib(phi)

    use mod_common

    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(inout) :: phi

    integer :: i,j,k, cnt
    integer, dimension(imax,jmax,kmax) :: mask

    real :: xsum,ysum,zsum
    real, dimension(-2:2) :: xscan,yscan,zscan    

    
    !-- Constant if outside of tube (regardless of riblet)

    do k = -2,k1+2
       do j = -2,j1+2
          do i = -2,i1+2
             
             if (phi(i,j,k) .gt. NarrowBand_2) then
                phi(i,j,k) = NarrowBand_2
             elseif (phi(i,j,k) .lt. -NarrowBand_2) then
                phi(i,j,k) = -NarrowBand_2
             endif

          enddo
       enddo
    enddo

    !-- Mask tube interior and skin (~2dx thick)
    
    cnt = 0

    do k = 1,kmax
       do j = 1,jmax
          do i = 1,imax

             if (j_sind(j) .eq. 0.) then  ! only in fluid domain
                
                if (abs(phi(i,j,k)) .lt. NarrowBand_2) then
                   mask(i,j,k) = 1
                   cnt = cnt + 1
                else
                   xscan(:) = abs(phi(i-2:i+2,j,k)) - NarrowBand_2
                   xsum = sum(xscan)

                   !yscan(:) = abs(phi(i,j-2:j+2,k)) - NarrowBand_2
                   !ysum = sum(yscan)
                   !-- rib tube only determined by zscan (for vertical walls)
                   yscan = 0.
                   ysum = 0.

                   zscan(:) = abs(phi(i,j,k-2:k+2)) - NarrowBand_2
                   zsum = sum(zscan)

                   if (xsum .eq. 0. .and. ysum .eq. 0. .and. zsum .eq. 0.) then
                      mask(i,j,k) = 0
                   else
                      mask(i,j,k) = 1
                      cnt = cnt + 1
                   endif
                endif

             else
                mask(i,j,k) = 0
             endif

          enddo
       enddo
    enddo

    nb_cnt = cnt
    
    !-- Create index arrays
    
    if (allocated(nb_i)) deallocate(nb_i)
    if (allocated(nb_j)) deallocate(nb_j)
    if (allocated(nb_k)) deallocate(nb_k)

    allocate (nb_i(nb_cnt))
    allocate (nb_j(nb_cnt))
    allocate (nb_k(nb_cnt))

    cnt = 0

    do k = 1,kmax
       do j = 1,jmax
          do i = 1,imax

             if (mask(i,j,k) .eq. 1) then

                cnt = cnt + 1

                nb_i(cnt) = i
                nb_j(cnt) = j
                nb_k(cnt) = k

             endif

          enddo
       enddo
    enddo

    !call boundsls(phi)

    return
  end subroutine tube_construction_rib

  
  !------------------------------------------------------------
  

  subroutine HOUC5_nb_rib(vel_u,vel_v,vel_w, phi, dfdt)

    use mod_common
    use mod_common_mpi

    !-- Advection velocities (cell-center)
    real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: vel_u,vel_v,vel_w

    !-- Scalar to be transported
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in) :: phi

    !-- Temporal derivative
    real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: dfdt

    integer :: i,j,k, cnt
    real :: udfdx,vdfdy,wdfdz, dfx,dfy,dfz

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

       if (s_left_1(j) .eq. 1. .or. s_right_1(j) .eq. 1.) then  ! 1st upwind for inner layer 1

          if ( vel_v(i,j,k) .GE. 0. ) then
             dfy = phi(i,j,k) - phi(i,j-1,k)
          else
             dfy = phi(i,j+1,k) - phi(i,j,k)                    
          endif

          vdfdy = -vel_v(i,j,k)*dfy/dy

       elseif (s_left_2(j) .eq. 1. .or. s_right_2(j) .eq. 1.) then  ! HOUC3 for inner layer 2

          if ( vel_v(i,j,k) .GE. 0. ) then
             dfy = (1./6.)*phi(i,j-2,k)-phi(i,j-1,k)+(1./2.)*phi(i,j,k)+(1./3.)*phi(i,j+1,k)  
          else
             dfy =-(1./6.)*phi(i,j+2,k)+phi(i,j+1,k)-(1./2.)*phi(i,j,k)-(1./3.)*phi(i,j-1,k)                          
          endif

          vdfdy = -vel_v(i,j,k)*dfy/dy
          
       else

          if ( vel_v(i,j,k) .GE. 0. ) then
             dfy = 1./60.*( -2.*phi(i,j-3,k)+15.*phi(i,j-2,k)-60.*phi(i,j-1,k)+ &
                  20.*phi(i,j,k)  +30.*phi(i,j+1,k)-3. *phi(i,j+2,k)  )        
          else
             dfy = 1./60.*(  2.*phi(i,j+3,k)-15.*phi(i,j+2,k)+60.*phi(i,j+1,k)- &
                  20.*phi(i,j,k)  -30.*phi(i,j-1,k)+3. *phi(i,j-2,k)  )          
          endif

          vdfdy = -vel_v(i,j,k)*dfy/dy

       endif

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
  end subroutine HOUC5_nb_rib

  
  !------------------------------------------------------------
  

  subroutine RK_nb_rib(phi_out, phi_in, factor,dfdt)

    use mod_common

    real, intent(in) :: factor
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in) :: phi_in
    real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: dfdt
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(out) :: phi_out

    integer :: i,j,k, cnt

    phi_out = phi_in

    do cnt = 1,nb_cnt

       i = nb_i(cnt)
       j = nb_j(cnt)
       k = nb_k(cnt)

       phi_out(i,j,k) = phi_in(i,j,k) + factor*dt*dfdt(i,j,k)

    enddo

    call boundsls_rib(phi_out)

    return
  end subroutine RK_nb_rib


  !------------------------------------------------------------
  
  
  subroutine boundsls_rib(phi)

    use mpi
    use mod_common_mpi
    use mod_bound

    real, dimension(-2:,-2:,-2:), intent(inout) :: phi

    integer :: i,j,k, jv
    real :: a1,a2
    real :: frontbound,backbound
    
    !-- Linear extrapolation at z-walls
    
    do j = -2,j1+2
       do i = -2,i1+2

          phi(i,j,0)  = 2.*phi(i,j,1)- phi(i,j,2)     
          phi(i,j,k1) = 2.*phi(i,j,kmax)- phi(i,j,kmax-1)  

       enddo
    enddo

    !##!-- Linear extrapolation at rib walls
    !##
    !##do k = 0,k1
    !##   do j = 1,jmax
    !##
    !##      a1 = s_left_0(j)
    !##      a2 = s_right_0(j)
    !##      
    !##      do i = 1,imax
    !##
    !##         if (a1 .eq. 1.) phi(i,j,k) = 2.*phi(i,j+1,k) - phi(i,j+2,k)
    !##         if (a2 .eq. 1.) phi(i,j,k) = 2.*phi(i,j-1,k) - phi(i,j-2,k)
    !##         
    !##      enddo
    !##   enddo
    !##enddo

    !-- Periodic in x and y (after riblet)

    call updthalos_ls(phi,1)
    call updthalos_ls(phi,2)

    !-- If not periodic in y  --!

    if (open_end_in_y) then
       
       frontbound  = (coords(2)  )*ly/(1.*dims(2)) ! front boundary of process myid
       backbound   = (coords(2)+1)*ly/(1.*dims(2)) ! back  boundary of process myid

       if (frontbound .eq. 0.) then
          do k=1,kmax
             do i=1,imax

                phi(i,0,k)  = 2.*phi(i,1,k) - phi(i,2,k)  ! linear extrapolation
                phi(i,-1,k) = 2.*phi(i,0,k) - phi(i,1,k)
                phi(i,-2,k) = 2.*phi(i,-1,k)- phi(i,0,k)

             enddo
          enddo
       endif

       if (backbound .eq. ly) then
          do k=1,kmax
             do i=1,imax

                phi(i,j1,k)   = 2.*phi(i,jmax,k)- phi(i,jmax-1,k)  ! linear extrapolation
                phi(i,j1+1,k) = 2.*phi(i,j1,k)  - phi(i,jmax,k)
                phi(i,j1+2,k) = 2.*phi(i,j1+1,k)- phi(i,j1,k)

             enddo
          enddo
       endif

    endif
    

    return
  end subroutine boundsls_rib


  !------------------------------------------------------------
  
  
  subroutine inflation_rib(l)

    use mod_common
    use mod_common_mpi
    use mod_misc
    use mod_interface
    use mod_bound

    integer, intent(in) :: l

    integer :: i,j,k
    real :: flag_s, phi_s, dd, temp, p_l, kappamax

    real, parameter :: s1 = 1.5*dz

    real, dimension(lmax) :: area, loss_rate
    real, dimension(1:imax,1:jmax,1:kmax) :: delta, heav_rib, u_lc,v_lc,w_lc, dfdt
    real, dimension(0:i1,0:j1,0:k1) :: dheav, dummy, kappa
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phi,dphi, phi_1,phi_2,phi_3


    dummy = 0.

    area(l) = 0.
    dheav = 1.
    phi = lvset(:,:,:,l)

    call get_curvature(phi,kappa)

    !##---- disregard rib boundary layers
    do k = 1,kmax
       do j = 1,jmax
          
          flag_s = s_left_0(j) + s_left_1(j) + s_right_0(j) + s_right_1(j)

          if (flag_s .ne. 0.) then
             do i = 1,imax
                kappa(i,j,k) = 0.
             enddo
          endif

       enddo
    enddo
    !##----
    
    call boundc(kappa)

    !## ----
    do k = 1,kmax
       do j = 1,jmax
          do i = 1,imax
             heav_rib(i,j,k) = 1. - j_sind(j)  ! 1 for fluid, 0 for rib
          enddo
       enddo
    enddo
    !## ----

    select case (ispeed_dependence)

    case('unfm')

       delta = Dirac('trig', l, s1)
       dheav = Heaviside('0','trig',l, s1)

       !-- Surface area --!

       !##temp = sum(delta*normal%m) *dx*dy*dz  !# i included abs grad phi here
       temp = sum(delta*heav_rib*normal%m) *dx*dy*dz  !######## 21 June 2017
       call MPI_allreduce(temp, area(l), 1,mpi_real8,mpi_sum,comm_cart,error)
       if (TwoD) area(l) = area(l)/lx

       !-- Mass loss --!

       temp = sum(dheav(1:imax,1:jmax,1:kmax))
       call MPI_allreduce(temp, volume(l), 1,mpi_real8,mpi_sum,comm_cart,error)
       volume(l) = 1. - volume(l)/(1.*itot*jtot*ktot)

       loss_rate(l) = lx*ly*lz*(volume(l)-init_volume(l))/(inflate_step*dt)
       if (TwoD) loss_rate(l) = loss_rate(l)/lx

       !-- Normal cell-center velocity --!

       temp = loss_rate(l)/area(l)

       !# not robust below
       ! u_lc(:,:,:) = - temp*delta(:,:,:)*normal(:,:,:)%x
       ! v_lc(:,:,:) = - temp*delta(:,:,:)*normal(:,:,:)%y
       ! w_lc(:,:,:) = - temp*delta(:,:,:)*normal(:,:,:)%z

       maxchk_glb = temp/s1  ! maximal magnitude of correction velocity


       !-- Normal velocity --!

       call boundp(dheav)

       do k=1,kmax
          do j=1,jmax
             do i=1,imax
                u_l(i,j,k) = - temp*(dheav(i+1,j,k)-dheav(i,j,k))/dx
                v_l(i,j,k) = - temp*(dheav(i,j+1,k)-dheav(i,j,k))/dy
                w_l(i,j,k) = - temp*(dheav(i,j,k+1)-dheav(i,j,k))/dz
             enddo
          enddo
       enddo

       call bounduvw(u_l,v_l,w_l)

       !-- Average at cell centers --!

       do k= 1,kmax
          do j= 1,jmax
             do i= 1,imax
                u_lc(i,j,k) = (u_l(i-1,j,k)+u_l(i,j,k))/2. *heav_rib(i,j,k)  !######## 21 June 2017
                v_lc(i,j,k) = (v_l(i,j-1,k)+v_l(i,j,k))/2. *heav_rib(i,j,k)  !##
                w_lc(i,j,k) = (w_l(i,j,k-1)+w_l(i,j,k))/2. *heav_rib(i,j,k)  !##
             enddo
          enddo
       enddo


    case('curv')

       delta = Dirac('trig', l,s1)
       dheav = Heaviside('0','trig',l,s1)

       !-- Surface area --!

       !temp = sum(delta*kappa(1:imax,1:jmax,1:kmax)) *dx*dy*dz
       temp = sum(delta*kappa(1:imax,1:jmax,1:kmax)*normal%m) *dx*dy*dz !# include abs grad phi here
       call MPI_allreduce(temp, area(l), 1,mpi_real8,mpi_sum,comm_cart,error)
       if (TwoD) area(l) = area(l)/lx

       !-- Mass loss --!

       temp = sum(dheav(1:imax,1:jmax,1:kmax))
       call MPI_allreduce(temp, volume(l), 1,mpi_real8,mpi_sum,comm_cart,error)
       volume(l) = 1. - volume(l)/(1.*itot*jtot*ktot)

       loss_rate(l) = lx*ly*lz*(volume(l)-init_volume(l))/(inflate_step*dt)
       if (TwoD) loss_rate(l) = loss_rate(l)/lx

       !-- Normal cell-center velocity --!

       temp = loss_rate(l)/area(l)

       ! not robust below
       !u_lc(:,:,:) = - temp*delta(:,:,:)*normal(:,:,:)%x*kappa(1:imax,1:jmax,1:kmax)
       !v_lc(:,:,:) = - temp*delta(:,:,:)*normal(:,:,:)%y*kappa(1:imax,1:jmax,1:kmax)
       !w_lc(:,:,:) = - temp*delta(:,:,:)*normal(:,:,:)%z*kappa(1:imax,1:jmax,1:kmax)

       maxchk_glb = temp/s1  ! maximal magnitude of correction velocity


       !-- Normal velocity --!  # a different way

       call boundp(dheav)

       do k=1,kmax
          do j=1,jmax
             do i=1,imax
                u_l(i,j,k) = - temp*(dheav(i+1,j,k)-dheav(i,j,k))/dx
                v_l(i,j,k) = - temp*(dheav(i,j+1,k)-dheav(i,j,k))/dy
                w_l(i,j,k) = - temp*(dheav(i,j,k+1)-dheav(i,j,k))/dz
             enddo
          enddo
       enddo

       call bounduvw(u_l,v_l,w_l)

       !-- Average at cell centers --!

       do k= 1,kmax
          do j= 1,jmax
             do i= 1,imax
                u_lc(i,j,k) = (u_l(i-1,j,k)+u_l(i,j,k))/2.*kappa(i,j,k)
                v_lc(i,j,k) = (v_l(i,j-1,k)+v_l(i,j,k))/2.*kappa(i,j,k)
                w_lc(i,j,k) = (w_l(i,j,k-1)+w_l(i,j,k))/2.*kappa(i,j,k)
             enddo
          enddo
       enddo

    end select


    !-- RK3
    
    select case (i_space)
    case('HOUC5')

       call HOUC5_nb_rib(u_lc,v_lc,w_lc, phi, dfdt)
       call RK_nb_rib(phi_1, phi, 1.,dfdt)

       call HOUC5_nb_rib(u_lc,v_lc,w_lc, phi_1, dfdt)
       call RK_nb_rib(phi_2, 3./4.*phi + 1./4.*phi_1, 1./4.,dfdt)

       call HOUC5_nb_rib(u_lc,v_lc,w_lc, phi_2, dfdt)
       call RK_nb_rib(phi_3, 1./3.*phi + 2./3.*phi_2, 2./3.,dfdt)

    case('WENO5')

       !call WENO5_nb_rib(u_lc,v_lc,w_lc, phi, dfdt)
       !call RK_nb_rib(phi_1, phi, 1.,dfdt)
       !
       !call WENO5_nb_rib(u_lc,v_lc,w_lc, phi_1, dfdt)
       !call RK_nb_rib(phi_2, 3./4.*phi + 1./4.*phi_1, 1./4.,dfdt)
       !
       !call WENO5_nb_rib(u_lc,v_lc,w_lc, phi_2, dfdt)
       !call RK_nb_rib(phi_3, 1./3.*phi + 2./3.*phi_2, 2./3.,dfdt)
    end select

    lvset(:,:,:,l) = phi_3


    return
  end subroutine inflation_rib


  !------------------------------------------------------------
  
  
  subroutine ls_reinit_rib(marker,maxiter,iter,residual)
    !                                                          
    ! Get signed distance function (phi_0 --> phi_new) by evolving
    !                                                          
    !   dø                                                     
    ! ------ = sign(ø0)*( 1 - |grad ø| )                       
    !  dtau                                                    
    !                                                          
    ! iteratively using EU1+RUSM1 or RK2/RK3+WENO5.            
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    use mod_common
    use mod_common_mpi
    use mod_hyperPDE

    integer i,j,k, iterations
    integer, intent(in) :: marker,maxiter
    integer, intent(out) :: iter
    real, intent(out) :: residual

    integer, dimension(2) :: ind_loc

    logical, parameter :: validation = .false.

    real :: dtau,smooth_width, grad_min,grad_max,c1

    real, dimension(1:imax,1:jmax,1:kmax) :: sgn
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: RHS, phi_0, phi_n1, phi_n2, phi_new


    phi_0(:,:,:) = lvset(:,:,:,marker)  ! time level n
    phi_new = phi_0                     ! time level n+1

    !---- Smoothed Sign Function ----!

    smooth_width = 1.5*dz

    do k=1,kmax
       do j=1,jmax
          do i=1,imax

             if (j_sind(j) .eq. 0.) then  ! only in fluid
                
                if (phi_0(i,j,k) .lt. -dz) then
                   sgn(i,j,k) = -1.
                elseif (phi_0(i,j,k) .gt. dz) then  
                   sgn(i,j,k) = 1.
                else
                   sgn(i,j,k) = phi_0(i,j,k)/sqrt(phi_0(i,j,k)**2 + dz**2)
                endif

                ! !# based on Heaviside of each level set
                !
                ! if (phi_0(i,j,k) .lt. -smooth_width ) then         
                !     sgn(i,j,k) = -1.
                ! elseif (phi_0(i,j,k) .gt. smooth_width ) then 
                !     sgn(i,j,k) =  1.
                ! else
                !     sgn(i,j,k) = phi_0(i,j,k)/smooth_width + &
                !                  (1./picon)*sin(phi_0(i,j,k)*picon/smooth_width)
                ! endif

             else

                sgn(i,j,k) = 0.

             endif

          enddo
       enddo
    enddo


    !---- Reinitialization Loop ----!

    dtau = 0.5*dz

    do iterations = 1,maxiter      

       select case (ls_reinit_time)
       case('RK2')  ! Mid-Point RK2

          call WENO5_reinit(sgn,phi_0,RHS)
          phi_n1 = phi_0 + 0.5*dtau*RHS
          phi_n2 = phi_0 + dtau*RHS
          call boundsls_rib(phi_n1)
          call boundsls_rib(phi_n2)

          call WENO5_reinit(sgn,phi_n2,RHS)
          phi_new = phi_n1 + 0.5*dtau*RHS
          call boundsls_rib(phi_new)

       case('RK3')  ! SSP RK3

          select case (ls_reinit_space)
          case('WENO5')

             call WENO5_reinit(sgn,phi_0,RHS)
             phi_n1 = phi_0 + dtau*RHS
             call boundsls_rib(phi_n1)

             call WENO5_reinit(sgn,phi_n1,RHS)
             phi_n2 = 3./4.*phi_0 + 1./4.*(phi_n1 + dtau*RHS)
             call boundsls_rib(phi_n2)

             call WENO5_reinit(sgn,phi_n2,RHS)
             phi_new = 1./3.*phi_0 + 2./3.*(phi_n2 + dtau*RHS)
             call boundsls_rib(phi_new)

          case('WENO5z')

             call WENO5z_reinit(sgn,phi_0,RHS)
             phi_n1 = phi_0 + dtau*RHS
             call boundsls_rib(phi_n1)

             call WENO5z_reinit(sgn,phi_n1,RHS)
             phi_n2 = 3./4.*phi_0 + 1./4.*(phi_n1 + dtau*RHS)
             call boundsls_rib(phi_n2)

             call WENO5z_reinit(sgn,phi_n2,RHS)
             phi_new = 1./3.*phi_0 + 2./3.*(phi_n2 + dtau*RHS)
             call boundsls_rib(phi_new)
          end select

       end select

       !-- check residual

       residual = maxval(abs(phi_new-phi_0))
       call mpi_allreduce(mpi_in_place,residual,1,mpi_real8,mpi_max,comm_cart,error)
       if (residual .lt. dz**3) then
          iter = iterations
          exit
       else
          iter = maxiter
       endif

       !-- check level set evolution when validating

       if (validation) then
          if (myid .eq. 0) write(*,*) 'iteration = ',iterations, residual
          !-- check min and max |grad phi|
          call get_curvature(phi_new,curv_cmn)
          grad_min = minval(normal%m)
          grad_max = maxval(normal%m)
          call mpi_allreduce(mpi_in_place,grad_min,1,mpi_real8,mpi_min,comm_cart,error)
          call mpi_allreduce(mpi_in_place,grad_max,1,mpi_real8,mpi_max,comm_cart,error)
          ind_loc = minloc(normal(itot/2,:,:)%m)
          c1 = phi_new(itot/2,ind_loc(1),ind_loc(2))
          if (abs(c1) .lt. NarrowBand_2) then
             write(6,'(A,1F8.3,A,F8.3,A,1F8.3,A,1F8.3)') 'At y =',(myid*jmax+ind_loc(1))*dy,' z =',ind_loc(2)*dz, &
                  ', got min |grad phi| =',grad_min, ' w/ phi =',c1/dz
          endif
          ind_loc = maxloc(normal(itot/2,:,:)%m)
          c1 = phi_new(itot/2,ind_loc(1),ind_loc(2))
          if (abs(c1) .lt. NarrowBand_2) then
             write(6,'(A,1F8.3,A,F8.3,A,1F8.3,A,1F8.3)') 'At y =',(myid*jmax+ind_loc(1))*dy,' z =',ind_loc(2)*dz, &
                  ', got max |grad phi| =',grad_max, ' w/ phi =',c1/dz
          endif
          if (myid .eq. 0) write(6,*)
          !-- output contours
          if (mod(iterations,5) .eq. 0) then
             lvset(:,:,:,marker) = phi_new(:,:,:)
             call post2d(iterations)
          endif
       endif

       !-- update

       phi_0 = phi_new

    enddo

    !---- Update Level Set ----!

    lvset(:,:,:,marker) = phi_new(:,:,:)


    return
  end subroutine ls_reinit_rib
  

  !------------------------------------------------------------


  subroutine update_materials_rib
    !                                                               
    ! Update material properties based on their Heaviside functions.           
    ! ------                                                        
    ! hnew_mls: Heaviside function of level set                     
    ! miu,rho : Viscosity, density at cell centers                  
    ! volume  : Volume fraction of phase 2 (eg. droplets or bubbles)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    use mod_common
    use mod_misc
    use mod_common_mpi

    integer :: i,j,k,l
    real :: half_width,volume_all
    real :: ph,php, theta
    real, dimension(0:i1,0:j1,0:k1) :: h_rho, buffer, h_temp

    integer, parameter :: average_method = 0
    character(len=4), parameter :: shape = 'trig'
    ! coresponding to the shape of the Heaviside function


    !-- Retrieve the single level set

    slset = lvset(0:i1,0:j1,0:k1,1)
    if (lmax .gt. 1) then
       do l = 1,lmax
          slset = min(slset, lvset(0:i1,0:j1,0:k1,l))
       enddo
    endif

    !-- Interface defined by Heaviside (1 in fluid 1, 0 in fluid 2)

    do l = 1,lmax
       hnew_mls(:,:,:,l) = Heaviside('0','trig', l, 1.5*dz)
    enddo

    jiemian = hnew_mls(:,:,:,1)
    if (lmax .gt. 1) then
       do l=1,lmax
          jiemian = min(jiemian, hnew_mls(:,:,:,l))
       enddo
    endif

    !-- Cell-center viscosity

    miu = miu1*jiemian + miu2*(1. - jiemian)

    !-- Volume fraction (of fluid 2) 

    do l = 1,lmax

       h_temp = Heaviside('0','trig', l, 1.5*dz)
       volume(l) = sum(h_temp(1:imax,1:jmax,1:kmax))

       call mpi_allreduce(volume(l),volume_all,1,mpi_real8,mpi_sum,comm_cart,error)
       volume(l) = 1. - volume_all/(1.*itot*jtot*ktot)

       if (l .eq. 1) then
          h_rho = h_temp
       else
          h_rho = min(h_rho,h_temp)
       endif

    enddo

    !-- Cell-center density and average density

    rho = rho1*h_rho + rho2*(1.-h_rho)
    rho_avr = rho1*(1. - sum(volume(:))) + rho2*sum(volume(:))

    !-- cell-face density

    if (rho1 .eq. rho2) then

       rho_u = rho1
       rho_v = rho1
       rho_w = rho1

    else

       select case (average_method)
       case(0)  ! arithmetic mean (Dodd & Ferrante)

          do k = 1,kmax
             do j = 1,jmax
                do i = 1,imax

                   rho_u(i,j,k) = (rho(i,j,k)+rho(i+1,j,k))/2.
                   rho_v(i,j,k) = (rho(i,j,k)+rho(i,j+1,k))/2.
                   rho_w(i,j,k) = (rho(i,j,k)+rho(i,j,k+1))/2.

                enddo
             enddo
          enddo

       case(1)  ! harmonic mean (Tanguy et al)
          ! (sharper than smearing with regularized Heaviside, adjust dt accordingly)

          do k = 1,kmax
             do j = 1,jmax
                do i = 1,imax

                   ph = slset(i,j,k)

                   if (ph .gt. 3.*dz) then  ! safely in fluid 1

                      rho_u(i,j,k) = rho1
                      rho_v(i,j,k) = rho1
                      rho_w(i,j,k) = rho1

                   elseif (ph .lt. -3.*dz) then  ! safely in fluid 2

                      rho_u(i,j,k) = rho2
                      rho_v(i,j,k) = rho2
                      rho_w(i,j,k) = rho2

                   else

                      php = slset(i+1,j,k)

                      if (ph*php .lt. 0.) then  ! cross interface in x-dir

                         theta = abs(php)/(abs(ph)+abs(php))
                         rho_u(i,j,k) = rho2*theta + rho1*(1.-theta)
                      else
                         if (ph .gt. 0.) then
                            rho_u(i,j,k) = rho1
                         else
                            rho_u(i,j,k) = rho2
                         endif
                      endif

                      php = slset(i,j+1,k)

                      if (ph*php .lt. 0.) then  ! cross interface in y-dir

                         theta = abs(php)/(abs(ph)+abs(php))
                         rho_v(i,j,k) = rho2*theta + rho1*(1.-theta)
                      else
                         if (ph .gt. 0.) then
                            rho_v(i,j,k) = rho1
                         else
                            rho_v(i,j,k) = rho2
                         endif
                      endif

                      php = slset(i,j,k+1)

                      if (ph*php .lt. 0.) then  ! cross interface in z-dir

                         theta = abs(php)/(abs(ph)+abs(php))
                         rho_w(i,j,k) = rho2*theta + rho1*(1.-theta)
                      else
                         if (ph .gt. 0.) then
                            rho_w(i,j,k) = rho1
                         else
                            rho_w(i,j,k) = rho2
                         endif
                      endif

                   endif

                enddo
             enddo
          enddo

       case(2)  ! based on the cell-face Heaviside

          half_width = 1.5*dz

          !-- u-grid --!

          h_temp = Heaviside('1',shape, 1 ,half_width)
          if (lmax .gt. 1) then
             do l = 1,lmax
                buffer = Heaviside('1',shape, l ,half_width)
                h_temp = min(h_temp,buffer)
             enddo
          endif

          rho_u = rho2 + (rho1-rho2)*h_temp(1:imax,1:jmax,1:kmax) ! note the size difference (same for v and w)

          !-- v-grid --!

          h_temp = Heaviside('2',shape, 1 ,half_width)
          if (lmax .gt. 1) then
             do l = 1,lmax
                buffer = Heaviside('2',shape, l ,half_width)
                h_temp = min(h_temp,buffer)
             enddo
          endif

          rho_v = rho2 + (rho1-rho2)*h_temp(1:imax,1:jmax,1:kmax)

          !-- w-grid --!

          h_temp = Heaviside('3',shape, 1 ,half_width)
          if (lmax .gt. 1) then
             do l = 1,lmax
                buffer = Heaviside('3',shape, l ,half_width)
                h_temp = min(h_temp,buffer)
             enddo
          endif

          rho_w = rho2 + (rho1-rho2)*h_temp(1:imax,1:jmax,1:kmax)

       end select

    endif

    return
  end subroutine update_materials_rib
  



end module mod_rib
