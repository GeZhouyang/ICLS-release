module mod_contact_line_common

  use mod_param
  
  implicit none

  integer, parameter :: num_cp = 2  ! number of contact points
  integer, parameter :: num_s = 2   ! now only works if = 2
  ! Number of surfaces containing at least one contact point.

  character(len=1), parameter :: cont_surf_normal = 'y'
  ! Orientation of the contact surface, specified by its normal direction.

  real, dimension(:), allocatable :: d2c_xp1,d2c_xp2
  ! Distance to the interface using extrapolated level set values on the contact surface.
  
  integer :: ibump  ! bump function array size
    
  real, parameter :: static_cont_angle = 140./180.*pi ! measured from fluid 1 (not used)
  real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: H_w  ! Heaviside for walls (not used)
  
  !-- Modeling parameters --!
  
  real, parameter :: h0 = 3.*dz  ! height of the bump
  real, parameter :: half_width = 6.*h0  ! half width of the bump

  
end module mod_contact_line_common


!-----------------------------------------------


module mod_contact_line

  use mod_param
  use mod_contact_line_common

  implicit none

  private
  public init_contact_line, contact_line_main, &
         comp_cp_vel, spread_slip_vel, enforce_contact_angle, output_contact
  
contains

  
  !------------------------------------------------------------------------------------------

  
  subroutine init_contact_line

    use decomp_2d
    use mpi
    use mod_common_mpi

    integer :: i,j,k
    real :: lowerbound, coory

    !-- Initialize the contact line, i.e. a 1D level set defined on the boundary

    select case(cont_surf_normal)
    case('z')
       
       ibump = jmax

       allocate (d2c_xp1(-2:j1+2))
       allocate (d2c_xp2(-2:j1+2))
       d2c_xp1 = 0.
       d2c_xp2 = 0.
       
    case('y')
       
       ibump = kmax

       allocate (d2c_xp1(-2:k1+2))
       allocate (d2c_xp2(-2:k1+2))
       d2c_xp1 = 0.
       d2c_xp2 = 0.
       
    end select

    !-- Identify the solid region

    H_w(:,:,-2:0) = 1.
    H_w(:,:,1:kmax) = 0.
    H_w(:,:,k1:k1+2) = 0.
    

    return
  end subroutine init_contact_line


  !------------------------------------------------------------------------------------------
  

  subroutine contact_line_main(istep, app_cont_angle)

    ! Compute the contact point velocity:  vel_c
    ! from the apparent contact angle:     app_cont_angle
    ! using a micro model (tabulated results of a separate phase-field simulation).
    ! Then spread the slip velocity around the contact point using a bump funtion.

    integer, intent(in) :: istep
    real, dimension(num_cp), intent(out) :: app_cont_angle

    real, dimension(num_cp) :: cont_loc, vel_c
    real, dimension(ibump,num_cp) :: vel_t,vel_n
    

    !-- Contact point locations & Apparent contact angles --!
    
    call xp_ls_comp_angle(istep, cont_loc,app_cont_angle)

    !-- Contact point velocity --!
    
    call comp_cp_vel(app_cont_angle,vel_c)

    !-- Spreaded slip velocity --!
    
    call spread_slip_vel(istep, app_cont_angle,vel_c, vel_t,vel_n)
    call mask_adv_vel(vel_t,vel_n)

    !-- Output --!

    if (mod(istep,100) .eq. 0) then
       call output_contact(istep,cont_loc,app_cont_angle,vel_c)
    endif

    return
  end subroutine contact_line_main


  !------------------------------------------------------------------------------------------
  

  subroutine xp_ls_comp_angle(istep, cont_loc, cont_ang)
    
    ! Extrapolate the level set to the wall and
    ! compute the contact point locations and apparent contact angles
    ! ## only for z-normal ##

    use mod_common
    use mod_common_mpi

    integer, intent(in) :: istep
    real, dimension(num_cp), intent(out) :: cont_loc, cont_ang

    integer, parameter :: xp_method = 2
    integer :: i,j,k, jp, jl,jr
    real :: ribl,ribr
    real :: ph,dphdz,d2phdz2
    real :: lowerbound,coory,coorz, d1,d2,d3,d4, e0,e1,e2, theta
    real, dimension(10) :: phi_v
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phi

    i = imax/2  ! 2D
    lowerbound = coords(2)*ly/(dims(2)*1.)
    
    phi = lvset(:,:,:,1)

    d2c_xp1 = 0.
    d2c_xp2 = 0.
    
    !-- Extrapolate level set to contact surfaces

    select case (cont_surf_normal)
    case('z')

       select case (xp_method)
       case(1)  ! linear upwind

          do j = 0,j1
             d2c_xp1(j) = 1.5*phi(i,j,1)    - 0.5*phi(i,j,2)       ! lower side
             d2c_xp2(j) = 1.5*phi(i,j,kmax) - 0.5*phi(i,j,kmax-1)  ! upper side
          enddo

       case(2)  ! least squares upwind

          do j = 0,j1
             phi_v = matmul(pseudo_inv_Ab_lower, reshape(phi(i-1:i+1,j-1:j+1,1:3), (/27/)))
             ph      = phi_v(1)
             dphdz   = phi_v(4)
             d2phdz2 = phi_v(7)
             d2c_xp1(j) = ph - dphdz*(dz/2.) + 0.5*d2phdz2*(dz/2.)**2  ! lower side

             phi_v = matmul(pseudo_inv_Ab_upper, reshape(phi(i-1:i+1,j-1:j+1,kmax-2:kmax), (/27/)))
             ph      = phi_v(1)
             dphdz   = phi_v(4)
             d2phdz2 = phi_v(7)
             d2c_xp2(j) = ph + dphdz*(dz/2.) + 0.5*d2phdz2*(dz/2.)**2  ! upper side
          enddo

       end select

       !-- Interpolate the interface location on contact surfaces
       !   (now works only if one contact point per surface)
       
       !## This chunk below shall be coded in a nicer way.
       
       !-- Lower surface
       
       e1 = 0.
       e2 = 0.

       do j = 1,jmax

          coory = lowerbound + (j-0.5)*dy

          ! lower side

          d1 = d2c_xp1(j)
          d2 = d2c_xp1(j+1)
          if (d1*d2 .lt. 0.) then
             e0 = dy/( 1. + abs(d2/d1) )
             e1 = coory + e0
          elseif (d1 .eq. 0.) then
             e1 = coory
          endif

          ! half layer up

          d3 = phi(i,j,1)
          d4 = phi(i,j+1,1)
          if (d3*d4 .lt. 0.) then
             e0 = dy/( 1. + abs(d4/d3) )
             e2 = coory + e0
          elseif (d3 .eq. 0.) then
             e2 = coory
          endif

       enddo

       call mpi_allreduce(mpi_in_place,e1,1,mpi_real8,mpi_max,comm_cart,error)
       call mpi_allreduce(mpi_in_place,e2,1,mpi_real8,mpi_max,comm_cart,error)

       !-- Compute the contact angle (assume flat interface)

       if (e1 .ne. e2) then
          theta = atan( dz/2./(e2-e1) )  ! measured from fluid 1 (phi > 0) 
       else
          theta = pi/2.
       endif
       if (theta .lt. 0.) theta = pi + theta

       !-- Contact point location & contact angle

       cont_loc(1) = e1
       cont_ang(1) = theta

       !-- Higher surface

       e1 = 0.
       e2 = 0.

       do j = 1,jmax

          coory = lowerbound + (j-0.5)*dy

          ! higher side

          d1 = d2c_xp2(j)
          d2 = d2c_xp2(j+1)
          if (d1*d2 .lt. 0.) then
             e0 = dy/( 1. + abs(d2/d1) )
             e1 = coory + e0
          elseif (d1 .eq. 0.) then
             e1 = coory
          endif

          ! half layer down

          d3 = phi(i,j,kmax)
          d4 = phi(i,j+1,kmax)
          if (d3*d4 .lt. 0.) then
             e0 = dy/( 1. + abs(d4/d3) )
             e2 = coory + e0
          elseif (d3 .eq. 0.) then
             e2 = coory
          endif

       enddo

       call mpi_allreduce(mpi_in_place, e1,1, mpi_real8, mpi_max, comm_cart,error)
       call mpi_allreduce(mpi_in_place, e2,1, mpi_real8, mpi_max, comm_cart,error)

       !-- Compute the contact angle (assume flat interface)

       if (e1 .ne. e2) then
          theta = atan( dz/2./(e2-e1) )  ! measured from fluid 1 (phi > 0) 
       else
          theta = pi/2.
       endif
       if (theta .lt. 0.) theta = pi + theta

       !-- Contact point location & contact angle

       cont_loc(2) = e1
       cont_ang(2) = theta
       
    end select
    

    return
  end subroutine xp_ls_comp_angle


  !------------------------------------------------------------------------------------------


  subroutine comp_cp_vel(theta,vel_c)

    ! Compute the single value contact point velocity as a function of the apparent contact angle,
    !
    !     U_c = f(theta),
    !
    ! from tabulated results of the micro model (Kronbichler & Kress),
    ! using a piece-wise cubic fit (coefs pre-determined in Matlab).
    
    real, dimension(num_cp), intent(in) :: theta
    real, dimension(num_cp), intent(out) :: vel_c

    integer :: i
    integer, parameter :: num_case = 2
    real :: p1,p2,p3,p4
    real :: x1,x2,x3, xc

    select case (num_case)
    case(0)  ! ex case in Martin & Gunilla (static angle 140 degree measured from fluid 1)
       
       do i = 1,num_cp

          x1 = theta(i)   ! in [rad]
          x2 = x1*x1
          x3 = x1*x2
          xc = 2.3419     ! the value separating two fits

          if (x1 .lt. xc) then

             p1 = -0.025174
             p2 =  0.046098
             p3 = -0.038078
             p4 =  0.15154

          else

             p1 = -0.41431
             p2 =  3.37
             p3 = -9.1416
             p4 =  8.2421

          endif

          vel_c(i) = p1*x3 + p2*x2 + p3*x1 + p4

       enddo

    case(1)  ! micro-cavity #1 (static angle 80 degree measure from fluid 1)
       
       do i = 1,num_cp

          x1 = theta(i)   ! in [rad]
          x2 = x1*x1
          x3 = x1*x2
          xc = 1.2566     ! the value separating two fits

          if (x1 .lt. xc) then

             p1 = -0.5095
             p2 = -2.125
             p3 =  0.51681  
             p4 =  4.5347 

          else

             p1 = -0.24572
             p2 =  3.9317 
             p3 = -15.384
             p4 =  14.457 

          endif

          vel_c(i) = p1*x3 + p2*x2 + p3*x1 + p4

       enddo

    case(2)  ! micro-cavity #2 (static angle 80 degree measure from fluid 1)
       
       do i = 1,num_cp

          x1 = theta(i)   ! in [rad]
          x2 = x1*x1
          x3 = x1*x2
          xc = 1.6755     ! the value separating two fits

          if (x1 .lt. xc) then

             p1 = -0.54828
             p2 =  0.93619
             p3 = -0.85725  
             p4 =  0.86382

          else

             p1 = -0.1209
             p2 = -2.674 
             p3 =  7.8861
             p4 = -5.6664 

          endif

          vel_c(i) = p1*x3 + p2*x2 + p3*x1 + p4

       enddo

    case(4)  ! micro-cavity #4 (static angle 80 degree measure from fluid 1)
       
       do i = 1,num_cp

          x1 = theta(i)   ! in [rad]
          x2 = x1*x1
          x3 = x1*x2
          xc = 2.1642     ! the value separating two fits

          if (x1 .lt. xc) then

             p1 = -0.8827
             p2 =  2.1089
             p3 = -1.9326  
             p4 =  1.0223 

          else

             p1 =  397.78
             p2 = -2798.6 
             p3 =  6350.7
             p4 = -5060.1 

          endif

          vel_c(i) = p1*x3 + p2*x2 + p3*x1 + p4

       enddo

    end select
    

    return
  end subroutine comp_cp_vel

  
  !------------------------------------------------------------------------------------------
  
  
  subroutine spread_slip_vel(istep, app_cont_angle,vel_c, vel_t,vel_n)

    integer, intent(in) :: istep
    real, dimension(num_cp), intent(in) :: app_cont_angle,vel_c
    real, dimension(ibump,num_cp), intent(out) :: vel_t,vel_n

    integer :: i,n
    integer, parameter :: ichk = iout2d

    real :: a1,b1,c1,d1, a2,b2,c2,d2 ! coeffs in velocity eqs
    real :: phi,theta, S,C, vel_r,vel_p, vel
    
    real, dimension(ibump) :: bump_angle


    do n = 1,num_cp

       phi = app_cont_angle(n)
       vel = vel_c(n)
       
       call comp_coef(phi, a1,b1,c1,d1,a2,b2,c2,d2)  ! coef as func of the apparent contact angle
    
       call bump_around_contact(phi, bump_angle)  ! bump around the contact point
    
       !-- Compute the slip velocity at those points
       
       do i = 1,ibump

          theta = bump_angle(i)
          S = sin(theta)
          C = cos(theta)

          ! polar coords (with the contact line fixed, the wall moving)

          if (theta .gt. phi) then

             vel_r = (b1 + d1*theta - c1)*S - (a1 + c1*theta + d1)*C
             vel_p = (a1 + c1*theta     )*S + (b1 + d1*theta     )*C

          else

             vel_r = (b2 + d2*theta - c2)*S - (a2 + c2*theta + d2)*C
             vel_p = (a2 + c2*theta     )*S + (b2 + d2*theta     )*C

          endif

          ! Cartesian coords (time -1 to switch to the lab reference frame)

          vel_t(i,n) = -1.*vel*(vel_r*C - vel_p*S - 1.)
          vel_n(i,n) = -1.*vel*(vel_r*S + vel_p*C)

       end do

    enddo
    
    
    return
  end subroutine spread_slip_vel

  
  !::::::::::::::::::::::::::::::::::::::::::::::::::::

  
  subroutine comp_coef(phi, a1,b1,c1,d1,a2,b2,c2,d2)

    ! Compute coefficients a,b,c,d

    real, intent(in) :: phi
    real, intent(out):: a1,b1,c1,d1, a2,b2,c2,d2
    
    real :: phi2, S,C,SC,S2, delta,Q, D, e1,e2
    

    phi2 = phi**2
    
    S = sin(phi)
    C = cos(phi)
    SC = S*C
    S2 = S**2
        
    delta = phi - pi
    !##Q = miu2  ! ambiguous (not sure if should be the inverse)
    Q = miu1/miu2  ! mu_A/mu_B, A is where the contact measured from
    
    D = (SC-phi)*(delta**2-S2) + Q*(delta-SC)*(phi2-S2)

    e1 = S2-delta*phi+Q*(phi2-S2)
    c1 = S2*e1/D
    d1 = SC*(e1-pi*tan(phi))/D
    a1 = -1.-pi*c1-d1
    b1 = -pi*d1

    e2 = S2-delta**2+Q*(delta*phi-S2)
    c2 = S2*e2/D
    d2 = SC*(e2-Q*pi*tan(phi))/D
    a2 = -1.-d2
    b2 = 0.
   

    return
  end subroutine comp_coef


  !::::::::::::::::::::::::::::::::::::::::::::::::::::
  

  subroutine bump_around_contact(phi, bump_angle)

    real, intent(in) :: phi
    real, dimension(1:ibump), intent(out) :: bump_angle

    integer :: i
    real :: x_val,y_val, c
    

    do i = 1,ibump

       x_val = d2c_xp1(i)/sin(phi)  ! distance to the contact point

       if (x_val .le. -half_width) then

          bump_angle(i) = pi

       elseif (x_val .ge. half_width) then

          bump_angle(i) = 0.

       else
          
          if (x_val .ne. 0.) then

             c = 1. - (x_val/half_width)**2
             y_val = h0*exp(1.-1./c)
       
             bump_angle(i) = atan(y_val/x_val)
             if (x_val .lt. 0.) bump_angle(i) = bump_angle(i) + pi
             if (bump_angle(i) .lt. 1e-15) bump_angle(i) = 0.  ! remove extremely small values
          
          else
             bump_angle(i) = pi/2.
          
          endif
          
       endif

    enddo

    return
  end subroutine bump_around_contact


  !------------------------------------------------------------------------------------------


  subroutine mask_adv_vel(vel_t,vel_n)

    use mod_common

    !  Update advection velocity within the narrow band

    real, dimension(ibump,num_cp), intent(in) :: vel_t,vel_n
    
    integer :: i,j,k
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phi

    
    select case (cont_surf_normal)
    case('z')
    
       phi = lvset(:,:,:,1)

       k = 1  ! bottom layer
       do j = 1,jmax
          do i = 1,imax
             if (abs(phi(i,j,k)) .le. NarrowBand_2) then  !# maybe I should do v=+v?

                vccnew(i,j,k) = vel_t(j,1)
                wccnew(i,j,k) = vel_n(j,1)

             endif
          enddo
       enddo
       
       k = kmax  ! top layer
       do j = 1,jmax
          do i = 1,imax
             if (abs(phi(i,j,k)) .le. NarrowBand_2) then  !# maybe I should do v=+v?

                vccnew(i,j,k) = vel_t(j,2)
                wccnew(i,j,k) = - vel_n(j,2)  !## mirror

             endif
          enddo
       enddo

    end select
    
    
    return
  end subroutine mask_adv_vel
  

  !------------------------------------------------------------------------------------------


  subroutine enforce_contact_angle(label,app_cont_angle)

    use mod_common
    use mod_hyperPDE

    integer, intent(in) :: label
    real, dimension(num_cp), intent(in) :: app_cont_angle

    integer :: i,j,k, km,kn
    real :: a,b,c, theta, d0,d1,d2
    real, dimension(-2:i1+2,-2:j1+2,-2:2) :: phi
    
    ! The code below is for the simplest case where the walls are
    ! horizontal, ie. nw = (0,0,1), such that the problem reduces
    ! to 1D and it is independent of parallelization.

    ! Coefficients for the second-order upwind scheme (positive speed)
    ! dfdz(k) = ( a*f(k)+b*f(km)+c*f(kn) )/dz.
    ! Only one set of coefs is needed due to symmetry.

    select case (cont_surf_normal)
    case('z')
       
       a = 1.5
       b = -2.
       c = 0.5

       !-- Lower wall --!

       theta = app_cont_angle(1)
       d0 = cos(theta)
       
       phi = lvset(:,:,-2:2,label)

       do k = 0,-2,-1
          km = k+1
          kn = k+2
          do j = -2,j1+2
             do i = -2,i1+2

                if (abs(phi(i,j,kn)) .lt. NarrowBand_2) then

                   d1 = 1./(1./dz*a)
                   d2 = -1./dz*( b*phi(i,j,km)+c*phi(i,j,kn) ) + d0
                   phi(i,j,k) = d1*d2

                endif

             enddo
          enddo
       enddo

       lvset(:,:,-2:0,label) = phi(:,:,-2:0)

       !-- Upper wall --!

       theta = app_cont_angle(2)
       d0 = cos(theta)

       phi = lvset(:,:,k1-2:k1+2,label)

       do k = 0,2
          km = k-1
          kn = k-2
          do j = -2,j1+2
             do i = -2,i1+2

                if (abs(phi(i,j,kn)) .lt. NarrowBand_2) then

                   d1 = 1./(1./dz*a)
                   d2 = -1./dz*( b*phi(i,j,km)+c*phi(i,j,kn) ) + d0
                   phi(i,j,k) = d1*d2

                endif

             enddo
          enddo
       enddo

       lvset(:,:,k1:k1+2,label) = phi(:,:,0:2)

    end select
    

    return
  end subroutine enforce_contact_angle


  !------------------------------------------------------------------------------------------


  subroutine output_contact(istep,cont_loc,app_cont_angle,vel_c)

    use mod_common_mpi

    integer, intent(in) :: istep
    real, dimension(num_cp), intent(in) :: cont_loc,app_cont_angle,vel_c

    integer :: fn
    real :: t

    fn = 57  ! file number
    t = istep*timestep

    if (myid .eq. 0) then

       ! Plot the global array

       open(fn,file=datadir//'contact_line.txt',position='append')
       write(fn,'(7ES16.8)') t, cont_loc(1),cont_loc(2), &
            app_cont_angle(1),app_cont_angle(2), &
            vel_c(1), vel_c(2)                                
       close(fn)

    endif

    
    return
  end subroutine output_contact


  !------------------------------------------------------------------------------------------


  subroutine plot_analytic_vel(istep, angle,vel_t,vel_n)

    use mod_common_mpi

    integer, intent(in) :: istep
    real, dimension(:), intent(in) :: angle, vel_t,vel_n

    character(len=7) :: filenumber
    integer :: i, l_y, i_g, ibump_tot
    integer :: npoints, sub_rank
    integer, dimension(2) :: l_coords  ! starts from 0
    
    real :: lowerbound,coory, theta
    real, dimension(:), allocatable :: angle_buf, vel_t_buf,vel_n_buf  ! buffer
    real, dimension(:), allocatable :: angle_tot, vel_t_tot,vel_n_tot  ! global
    

    lowerbound = coords(2)*ly/(dims(2)*1.)
    
    select case(cont_surf_normal)
    case('z')
       npoints = jmax+6 ! total # of (including halos) points in each processor
       ibump_tot = jtot
       allocate ( angle_buf(jtot) )
       allocate ( vel_t_buf(jtot) )
       allocate ( vel_n_buf(jtot) )
       allocate ( angle_tot(jtot) )
       allocate ( vel_t_tot(jtot) )
       allocate ( vel_n_tot(jtot) )
       
    end select
    
    l_coords(1) = 0
    
    do l_y = 0,dims(2)-1  ! Loop thru y ranks
       l_coords(2) = l_y
       
       if (l_y .ne. 0) then
          
          ! Determine the rank of current sub-domain

          call MPI_CART_RANK(comm_cart,l_coords,sub_rank,error)

          ! Pass data from sub_rank to rank 0
    
          if (myid .eq. sub_rank) then
    
            call MPI_SSEND(angle, npoints, MPI_REAL8, 0, 61,comm_cart,error)
            call MPI_SSEND(vel_t, npoints, MPI_REAL8, 0, 62,comm_cart,error)
            call MPI_SSEND(vel_n, npoints, MPI_REAL8, 0, 63,comm_cart,error)
    
          endif

          if (myid .eq. 0) then

            ! Processor 0 recieves data
    
            call MPI_RECV(angle_buf, npoints, MPI_REAL8, sub_rank, 61,comm_cart,status,error)
            call MPI_RECV(vel_t_buf, npoints, MPI_REAL8, sub_rank, 62,comm_cart,status,error)
            call MPI_RECV(vel_n_buf, npoints, MPI_REAL8, sub_rank, 63,comm_cart,status,error)
    
            ! Processor 0 distributes global values
            
            do i = 1,ibump
               i_g = i + l_coords(2)*jmax

               angle_tot(i_g) = angle_buf(i)
               vel_t_tot(i_g) = vel_t_buf(i)
               vel_n_tot(i_g) = vel_n_buf(i)
            enddo
    
          endif        
       endif
    enddo
    
    if (myid .eq. 0) then

       ! Finally processor 0 updates itself

       do i = 1,ibump
          angle_tot(i)   = angle(i)
          vel_t_tot(i)   = vel_t(i)
          vel_n_tot(i)   = vel_n(i)
       enddo

       ! Plot the global array

       write(filenumber,'(i7.7)') istep
       open(51,file=datadir//'ana_vel'//filenumber//'.txt')

       do i = 1,ibump_tot
          theta = angle_tot(i)
          coory = (i-0.5)*dy !#- contact_loc !## perhaps should substract cont_loc
          write(51,'(4ES16.8)') coory, theta/pi*180., vel_t_tot(i),vel_n_tot(i)
       enddo

       close(51)

    endif

    
    return
  end subroutine plot_analytic_vel

  
end module mod_contact_line
