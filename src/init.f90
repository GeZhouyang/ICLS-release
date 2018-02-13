module mod_init

  use mod_common
  use decomp_2d
  use mod_param
  use mpi
  use mod_common_mpi
  use mod_misc

  implicit none

  private
  public init_program

contains

  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  !

  subroutine init_program(begin, duold,dvold,dwold,pold,curv_restart)

    use mod_loadflds
    use mod_initsolver
    use mod_bound
    use mod_chkdiv
    use mod_param_ls
    use mod_interface
    use mod_ab
    use mod_mom
    use mod_rib
    use mod_contact_line

    integer, intent(in) :: begin

    real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: duold,dvold,dwold
    real, dimension(0:i1,0:j1,0:k1), intent(out) :: pold, curv_restart

    call write_log_header

    !---- Velocity, Pressure, Level Set

    if (begin .eq. 0) then
       time = 0. 
       call init_vel_p
       call init_levelset
    else
       call init_vel_p ! still need to re-init the vel in the case of partial inflow
       call loadflds(0,begin)
       if (myid .eq. 0) write(6,*) 'Restarting from step ',begin
    endif

    if (Riblet) then
       call init_static_ls
       call init_rib
    endif

    if (Contact_line) then
       call init_contact_line
    endif

    call init_matrix  ! least squares coef matrix
    call initsolver   ! pressure solver coef (real and Fourier spaces)

    !# before, pp2 => p

    dt = timestep  ! time step
    
    !---- Boundary Conditions (Physical and Computational)
    
    call bounduvw(unew,vnew,wnew)
    call bounduvw(dudt,dvdt,dwdt)
    call chkdiv(begin)
    call boundmls(lvset)
    call boundp(pnew)
    if (Riblet) call mask_rib_vel(unew,vnew,wnew)

    call average_vel_at_center  ! cell-center velocity

    !call ls_bound_error  ! check extrap. error

    !---- Phase, Density, Viscosity, and Volume

    call update_materials
    init_volume = volume
    
    !---- Obtain (n-1) values when (re)starting AB2 (effectively Euler's method)
    
    call mom_main(duold,dvold,dwold)

    ! only for semi-Lagrangian
    select case (ls_adv_space)
    case('semiLag')
       uccold = uccnew
       vccold = vccnew
       wccold = wccnew
    end select

    select case (surface_tension_method)
    case('CSF')
       call get_curvature(ynew,curv_restart)
       call boundc(curv_restart)
       call surf_ten(curv_restart,surfxold,surfyold,surfzold)

    case('GFM')
       pold = pnew
       p_xold = p_x
       p_yold = p_y
       p_zold = p_z
       dvndn = 0.  ! normal derivative of normal velocity
    end select


    
    return
  end subroutine init_program


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Initialize velocity and pressure
  !
  subroutine init_vel_p
    !                               
    ! unew,vnew,wnew : Velocity     
    ! dudt,dvdt,dwdt : Acceleration 
    ! pnew           : Pressure     
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    use mod_loadflds

    integer :: i,j,k, iglob,jglob,kglob

    real :: coory,coorz, lowerbound
    real :: v_upper,v_lower,shear_rate
    real :: num, rn1,rn2,rn3
    real :: umean,vmean,wmean
    real :: umean_all,vmean_all,wmean_all
    real :: alpha,beta,pC1,pD1,CsPgrad  ! two-layer flow parameters

    real, dimension(kmax) :: v
 
 
    lowerbound =  coords(2)*ly/(dims(2)*1.)

    select case(iniu)

    case('zer')

       forall(i=0:i1,j=0:j1,k=0:k1)
          unew(i,j,k) = 0.
          vnew(i,j,k) = 0.
          wnew(i,j,k) = 0.
          dudt(i,j,k) = 0.
          dvdt(i,j,k) = 0.
          dwdt(i,j,k) = 0.
          pnew(i,j,k) = 0.
       end forall
  
  
    case('uni')

       forall(i=0:i1,j=0:j1,k=0:k1)
          unew(i,j,k) = 0.
          vnew(i,j,k) = 1.
          wnew(i,j,k) = 0.  
          dudt(i,j,k) = 0.
          dvdt(i,j,k) = 0.
          dwdt(i,j,k) = 0.
          pnew(i,j,k) = 0.
       end forall


    case('poi')

       do k = 0,k1
          coorz = (k-0.5)*dz/lz ! [-dz/2 , 1+dz/2]
          do j = 0,j1
             do i = 0,i1
                unew(i,j,k) = 0.
                vnew(i,j,k) = 6.*v_bulk_desired*coorz*(1.-coorz)  ! max 1.5*bulk
                wnew(i,j,k) = 0.
                dudt(i,j,k) = 0.
                dvdt(i,j,k) = 0.
                dwdt(i,j,k) = 0.
                pnew(i,j,k) = 0.
             enddo
          enddo
       enddo
  
    case('cou')

       v_upper = 1.5  ! upper plate velocity
       v_lower =-1.5  ! lower plate velocity
       shear_rate = (v_upper - v_lower)/lz

       do k = 0,k1
          coorz = (k-0.5)*dz ! [-dz/2 , lz+dz/2]
          do j = 0,j1
             do i = 0,i1
                unew(i,j,k) = 0.
                vnew(i,j,k) = v_lower + shear_rate*coorz
                wnew(i,j,k) = 0.
                dudt(i,j,k) = 0.
                dvdt(i,j,k) = 0.
                dwdt(i,j,k) = 0.
                pnew(i,j,k) = 0.
             enddo
          enddo
       enddo       

    case('sep') !# 2D deformation vortex (Serpentine)

       do k = 0,k1
          coorz = (k-0.5)*dz/lz ! [-dz/2 , 1+dz/2]
          do j = 0,j1
             coory = lowerbound +(j)*dy/lz ! [0 , ly/lz+dy]
             do i = 0,i1
                unew(i,j,k) = 0.
                vnew(i,j,k) = -2. * (sin(pi*coory))**2 * sin(pi*coorz)*cos(pi*coorz)

                dudt(i,j,k) = 0.
                dvdt(i,j,k) = 0.
                dwdt(i,j,k) = 0.
                pnew(i,j,k) = 0.
             enddo
          enddo
       enddo
       do k = 0,k1
          coorz = (k)*dz/lz ! [0 , 1+dz]
          do j = 0,j1
             coory = lowerbound +(j-.5)*dy/lz ! [-dy/2 , ly/lz+dy/2]
             do i = 0,i1

                wnew(i,j,k) = 2. * sin(pi*coory)*cos(pi*coory) * (sin(pi*coorz))**2

             enddo
          enddo
       enddo


    case('rot') ! 2D rigid-body rotation

       do k = 0,k1
          coorz = (k-.5)*dz/lz ! [-dz/2 , 1+dz/2]
          do j = 0,j1
             coory = lowerbound + (j)*dy/lz ! [0 , ly/lz+dy]

             unew(:,j,k) = 0.
             vnew(:,j,k) =  2.*pi*(coorz-.5)
             dudt(:,j,k) = 0.
             dvdt(:,j,k) = 0.
             dwdt(:,j,k) = 0.
             pnew(:,j,k) = 0.
          enddo
       enddo
       do k = 0,k1
          coorz = (k)*dz/lz ! [0 , 1+dz]
          do j = 0,j1
             coory = lowerbound + (j-.5)*dy/lz ! [-dy/2 , ly/lz+dy/2]

             wnew(:,j,k) = -2.*pi*(coory-.5)        
          enddo
       enddo

  
    case('2lc')

       alpha = miu1/miu2
       beta = 0.2  !# interface height or lower layer thickness
       pC1 = - (1.-beta**2*(1.-alpha))/(1.-beta*(1.-alpha))
       pD1 = -1.-pC1
       CsPgrad = -10.*visc  !times visc to ensure velocities are of abuot order 1
       call random_seed ()  !# generating random numbers

       do k=0,k1
          coorz = (k-0.5)/dzi/lz ! normalised with channel height
          do j=0,j
             do i=0,i1
                !# exact solution for two-layer flow
                if (coorz .gt. beta) then
                   vnew(i,j,k) = CsPgrad/miu1/visc*(coorz*coorz+pC1*coorz+pD1)  !# upper layer
                else
                   vnew(i,j,k) = CsPgrad/miu2/visc*(coorz*coorz+pC1*coorz)      !# lower layer
                endif
                call random_number (num)  !#
                unew(i,j,k) = 0.
                wnew(i,j,k) = num*1.e-6  !# small perturbation
                dudt(i,j,k) = 0.
                dvdt(i,j,k) = 0.
                dwdt(i,j,k) = 0.
                pnew(i,j,k) = 0.
             enddo
          enddo
       enddo

       v_bulk = 0.
       do k=1,kmax
          do i=1,itot
             v_bulk = v_bulk+vnew(i,j,k)*dz*dx
          enddo
       enddo
       vnew(:,:,:) = lx*lz/v_bulk*vnew(:,:,:)  !#  make sure v_bulk = 1. initially

  
    case('log') ! Logarithmic profile
       
       !## ---- if(isfreeslip) then
       if (.not. ns_np) then
          do k=1,kmax
             v(k) = 2.5*log( (1./visc)*(k-0.5)/(2.*kmax) ) + 5.5
             if ((k-0.5)/kmax/visc .le. 11.6) v(k)=(k-0.5)/kmax/visc
          enddo
       else
          do k=1,kmax/2
             v(k) = 2.5*log( (1./visc)*(k-0.5)/(1.*kmax) ) + 5.5
             if ((k-0.5)/kmax/visc .le. 11.6) v(k)=(k-0.5)/kmax/visc
          enddo
          do k=1,kmax/2
             i = kmax+1-k
             v(kmax+1-k) = 2.5*log( (1./visc)* (k-0.5)/(1.*kmax) ) + 5.5
             if ((k-0.5)/kmax/visc .le. 11.6) v(kmax+1-k)=(k-0.5)/kmax/visc
          enddo
       endif
       !
       ! Below, random numbers are generated such that the initial field
       ! does not depend in the number of mpi tasks.
       !
       umean=0.
       vmean=0.
       wmean=0.
       do kglob=1,ktot
          do jglob=1,jtot
             do iglob=1,itot
                call random_number(rn1)
                call random_number(rn2)
                call random_number(rn3)
                i = iglob-imax*coords(1)
                j = jglob-jmax*coords(2)
                k = kglob
                if( &
                     (1.le.i.and.i.le.imax) .and. &
                     (1.le.j.and.j.le.jmax) &
                     ) then
                   unew(i,j,k)=15.*(rn1-0.5)
                   umean=umean+unew(i,j,k)
                   vnew(i,j,k)=v(k)*(1.+5.*(rn2-0.5))
                   vmean=vmean+vnew(i,j,k)
                   wnew(i,j,k)=15.*(rn3-0.5)
                   wmean=wmean+wnew(i,j,k)
                   dudt(i,j,k)=0.
                   dvdt(i,j,k)=0.
                   dwdt(i,j,k)=0.
                   pnew(i,j,k)=0.
                endif
             enddo
          enddo
       enddo

       call mpi_allreduce(umean,umean_all,1,mpi_real8,mpi_sum,comm_cart,error)
       call mpi_allreduce(vmean,vmean_all,1,mpi_real8,mpi_sum,comm_cart,error)
       call mpi_allreduce(wmean,wmean_all,1,mpi_real8,mpi_sum,comm_cart,error)
       umean_all = umean_all/(1.*itot*jtot*kmax)
       vmean_all = vmean_all/(1.*itot*jtot*kmax)
       wmean_all = wmean_all/(1.*itot*jtot*kmax)
       do k=1,kmax
          do j=1,jmax
             do i=1,imax
                unew(i,j,k) = unew(i,j,k)-umean_all !average = 0
                vnew(i,j,k) = vnew(i,j,k)/vmean_all !average = 1
                wnew(i,j,k) = wnew(i,j,k)-wmean_all !average = 0
             enddo
          enddo
       enddo

       
    case('fld')  ! load previous fld data

       forall(i=0:i1,j=0:j1,k=0:k1)
          dudt(i,j,k) = 0.
          dvdt(i,j,k) = 0.
          dwdt(i,j,k) = 0.
       end forall

       call init_from_fld  ! initialize unew,vnew,wnew,pnew (only 1:max)

    end select

    !-- Initialize inflow profile (if any) --!
    
    select case (BC_in_y)
    case('InOut')
       call inflow_profile
    end select
  
    !-- Initialize pressure jump --!
    
    p_xold = 0.
    p_yold = 0.
    p_zold = 0.
    p_x = 0.
    p_y = 0.
    p_z = 0.


    return
  end subroutine init_vel_p


  !------------------------------------
  ! Speicify an inflow velocity profile
  !
  subroutine inflow_profile

    integer :: i,j,k, i_glb, num
    integer, parameter :: num_chunk = 1  ! number of inflow chunks
    integer, dimension(num_chunk) :: i_start,i_end, k_start,k_end

    real :: ratio_i,ratio_k
    real :: v_mean, frontbound,coorz

    v_inflow = 0.

    !-- Starting and ending indices for i and k
    
    i_start(1) = itot/3
    i_end(1)   = itot/3*2

    k_start(1) = 0
    k_end(1)   = 9 !ktot

    !-- Ratios of the non-zero inflow in each Cartesian dir
    
    ratio_i = 1.0*(i_end(1)-i_start(1))/itot
    ratio_k = 1.0*(k_end(1)-k_start(1))/ktot
    
    v_mean = v_bulk_desired/(ratio_i*ratio_k)
    
    frontbound = coords(2)*ly/(1.*dims(2)) ! front boundary of process myid

    if (frontbound .eq. 0.) then
       do k = 1,kmax
          do num = 1,num_chunk
             if (k .gt. k_start(num) .and. k .le. k_end(num)) then ! k-dir is not parallelized
                coorz = (k - 0.5)*dz/(lz*ratio_k)
                do i = 1,imax
                   i_glb = i + coords(1)*imax  ! global i index
                   if (i_glb .gt. i_start(num) .and. i_glb .le. i_end(num)) then
                      v_inflow(i,k) = 6.*v_mean*coorz*(1.-coorz)  ! parabola of max = 1.5 (z from 0 to lz*ratio_k)
                   endif
                enddo ! i loop
             endif
          enddo ! num loop
       enddo! k loop
    endif

    return
  end subroutine inflow_profile

 
  !~~~~~~~~~~~~~~~~~~~~~
  ! Initialize level set
  !
  subroutine init_levelset
  !                                                           
  !     ynew    : Level set function                          
  !       + for phase 1                                       
  !       - for phase 2                                       
  !                                                           
  !     dydt    : No meaning (kept for historical reason)     
  !                                                           
  !     coor?   : Global coordinates in ? direction (from mpi)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    integer :: i,j,k,l
    real :: coorx,coory,coorz
    real :: lowerbound,upperbound,leftbound,rightbound

    real :: alpha,beta,gamma
    real :: phi_vc,phi_hc  ! Zalesak's disk
    real :: ypc,zpc,radius, gap,center_dist    ! pancake parameters
    real :: xdp,zdp,ydp, ydp1,ydp2,ydp3,zdp1   ! droplet parameters
    real :: phi_periodic1,phi_periodic2,phi_periodic3
    real :: xs,ys,zs,rs, s_theta,s_phi, sin_the,cos_the,sin_phi,cos_phi    ! Hadamard
    real :: coef_a,coef_b,coef_c,coef_d, vel_r,vel_t

    real, dimension(:), allocatable :: yz_angle, v_angle,h_angle        ! relative angles of droplets/bubbles
    real, dimension(lmax) :: xdp_mls,ydp_mls,zdp_mls,radius_mls  ! multiple level set
    
    ! sphere packings 
    real, dimension(6)  :: x_6pac,y_6pac,z_6pac
    real, dimension(7)  :: x_7pac,y_7pac,z_7pac
    real, dimension(8)  :: x_8pac,y_8pac,z_8pac
    real, dimension(9)  :: x_9pac,y_9pac,z_9pac
    real, dimension(10) :: x_10pac,y_10pac,z_10pac

    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phi_pancake, phi_slot
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phi_drp1, phi_drp2, phi_drp3

  
    ynew(:,:,:) = 0.
    dydt(:,:,:) = 0.
    volume = 0.
  
    !---- Bounds of Each Core ----!
    
    lowerbound =  coords(2)   *ly/(dims(2)*1.)
    upperbound = (coords(2)+1)*ly/(dims(2)*1.)
    leftbound  =  coords(1)   *lx/(dims(1)*1.)
    rightbound = (coords(1)+1)*lx/(dims(1)*1.)


    select case(inil)

    case('pan')

       ! Center (ydp,zdp) and Radius

       ypc = ly*0.5
       zpc = 0.5+0.25*sqrt(3.)
       radius = 0.5

       do k=-2,k1+2
          coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
          do j=-2,j1+2
             coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
             do i=-2,i1+2

                ynew(i,j,k) = sqrt((coory-ypc)**2+(coorz-zpc)**2) - radius
                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       lvset(:,:,:,1) = ynew(:,:,:)

  
    case('2pc')

       xdp = 0.5*lx
       zdp = 0.5*lz
       radius = 0.5
       ydp1= 0.5*ly -(radius + 3.*dy)
       ydp2= 0.5*ly +(radius + 3.*dy)

       do k=-2,k1+2
          coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
          do j=-2,j1+2
             coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
             do i=-2,i1+2

                !-- Define 2 pancakes separately

                phi_drp1(i,j,k) = sqrt( (coory-ydp1)**2+(coorz-zdp)**2 ) -radius
                phi_drp2(i,j,k) = sqrt( (coory-ydp2)**2+(coorz-zdp)**2 ) -radius

                !-- Consider periodic BC in y direction

                phi_periodic1 = sqrt( (coory-ydp1-ly)**2+(coorz-zdp)**2 )-radius
                phi_periodic2 = sqrt( (coory-ydp2-ly)**2+(coorz-zdp)**2 )-radius

                if (phi_drp1(i,j,k) .gt. 0.) phi_drp1(i,j,k) = min(phi_drp1(i,j,k), phi_periodic1)
                if (phi_drp2(i,j,k) .gt. 0.) phi_drp2(i,j,k) = min(phi_drp2(i,j,k), phi_periodic2)

                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       do k=-2,k1+2
          do j=-2,j1+2
             do i=-2,i1+2

                !-- Level set as min distance to interface

                if (phi_drp1(i,j,k) .le. 0.) ynew(i,j,k) = phi_drp1(i,j,k)
                if (phi_drp2(i,j,k) .le. 0.) ynew(i,j,k) = phi_drp2(i,j,k)

                if (phi_drp1(i,j,k) .gt. 0. .and. &
                     phi_drp2(i,j,k) .gt. 0.        ) then

                   ynew(i,j,k) = min( phi_drp1(i,j,k),phi_drp2(i,j,k) )

                endif

                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       lvset(:,:,:,1) = ynew(:,:,:)

       ! Initialize collision velocity

       where (phi_drp1(0:i1,0:j1,0:k1) .le. 1.5*dy) vnew = .5
       where (phi_drp2(0:i1,0:j1,0:k1) .le. 1.5*dy) vnew =-.5   


    case('3pc')

       xdp = .5*lx
       zdp = .5*lz
       radius = 1./6.*lz
       ydp1= .1*ly
       ydp2= ydp1 + 2.*radius+radius/5.
       ydp3= ydp2 + 2.*radius+radius/5.

       do k=-2,k1+2
          coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
          do j=-2,j1+2
             coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
             do i=-2,i1+2

                !-- Define 3 pancakes separately

                phi_drp1(i,j,k) = sqrt( (coory-ydp1)**2+(coorz-zdp)**2 ) -radius
                phi_drp2(i,j,k) = sqrt( (coory-ydp2)**2+(coorz-zdp)**2 ) -radius
                phi_drp3(i,j,k) = sqrt( (coory-ydp3)**2+(coorz-zdp)**2 ) -radius

                !-- Consider periodic BC in y direction

                phi_periodic1 = sqrt( (coory-ydp1-ly)**2+(coorz-zdp)**2 )-radius
                phi_periodic2 = sqrt( (coory-ydp2-ly)**2+(coorz-zdp)**2 )-radius
                phi_periodic3 = sqrt( (coory-ydp3-ly)**2+(coorz-zdp)**2 )-radius

                if (phi_drp1(i,j,k) .gt. 0.) phi_drp1(i,j,k) = min(phi_drp1(i,j,k), phi_periodic1)
                if (phi_drp2(i,j,k) .gt. 0.) phi_drp2(i,j,k) = min(phi_drp2(i,j,k), phi_periodic2)
                if (phi_drp3(i,j,k) .gt. 0.) phi_drp3(i,j,k) = min(phi_drp3(i,j,k), phi_periodic3)

                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       do k=-2,k1+2
          do j=-2,j1+2
             do i=-2,i1+2

                !-- Level set as min distance to interface

                if (phi_drp1(i,j,k) .le. 0.) ynew(i,j,k) = phi_drp1(i,j,k)
                if (phi_drp2(i,j,k) .le. 0.) ynew(i,j,k) = phi_drp2(i,j,k)
                if (phi_drp3(i,j,k) .le. 0.) ynew(i,j,k) = phi_drp3(i,j,k)

                if (phi_drp1(i,j,k) .gt. 0. .and. &
                     phi_drp2(i,j,k) .gt. 0. .and. &
                     phi_drp3(i,j,k) .gt. 0.        ) then

                   ynew(i,j,k) = min( phi_drp1(i,j,k),phi_drp2(i,j,k),phi_drp3(i,j,k) )

                endif

                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       lvset(:,:,:,1) = ynew(:,:,:)
  
  
    case('mpc')

       !-- Center of each droplet --!

       radius = 0.5
       xdp_mls = 0.5*lx
       ydp_mls(1) = 0.5*ly
       zdp_mls(1) = 0.25*lz

       gap = 3.*dy
       allocate ( yz_angle(3) )
       yz_angle = (/0., 45., 0./) /180.*pi

       do l = 2,lmax

          ydp_mls(l) = ydp_mls(l-1) + (2.*radius + gap)*cos(yz_angle(l))
          zdp_mls(l) = zdp_mls(l-1) + (2.*radius + gap)*sin(yz_angle(l))

       enddo

       do l = 1,lmax

          !-- Different level set for each pancake

          do k=-2,k1+2
             coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
             do j=-2,j1+2
                coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
                do i=-2,i1+2
                   coorx = leftbound + (i-0.5)/dxi ! [-dx*5/2 , lx+dx*5/2]

                   !-- Define 3 pancakes separately

                   lvset(i,j,k,l) = sqrt( (coory-ydp_mls(l))**2 + &
                        (coorz-zdp_mls(l))**2   ) -radius

                   !-- Consider periodic BC in y direction

                   phi_periodic1 = sqrt( (coory-ydp_mls(l)-ly)**2+(coorz-zdp_mls(l))**2 )-radius

                   if (lvset(i,j,k,l) .gt. 0.) lvset(i,j,k,l) = min(lvset(i,j,k,l),phi_periodic1)

                enddo
             enddo
          enddo

       enddo

       ynew(:,:,:) =  lvset(:,:,:,1)

  
    case('drp')

       ! Sphere centered at (xdp,ydp,zdp) with radius

       xdp = 0.5*lx
       ydp = 0.5*ly 
       zdp = 0.5*lz
       radius = 0.5      

       do k=-2,k1+2
          coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
          do j=-2,j1+2
             coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
             do i=-2,i1+2
                coorx = leftbound + (i-0.5)/dxi ! [-dx*5/2 , lx+dx*5/2]

                ynew(i,j,k) = sqrt((coorx-xdp)**2+(coory-ydp)**2+(coorz-zdp)**2) - radius

                !-- Consider periodic BC in y direction

                phi_periodic1 = sqrt((coorx-xdp)**2+(coory-ydp-ly)**2+(coorz-zdp)**2)-radius

                if (ynew(i,j,k) .gt. 0.) ynew(i,j,k) = min(ynew(i,j,k), phi_periodic1)

                !-- Consider periodic BC in x direction

                phi_periodic1 = sqrt((coorx-xdp-lx)**2+(coory-ydp)**2+(coorz-zdp)**2)-radius

                if (ynew(i,j,k) .gt. 0.) ynew(i,j,k) = min(ynew(i,j,k), phi_periodic1)

                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       lvset(:,:,:,1) = ynew(:,:,:)

  
    case('2dp')

       xdp = 0.5*lx
       zdp = 0.5*lz
       radius = 0.5
       ydp1= 0.5*ly -(radius + 1.*dy)
       ydp2= 0.5*ly +(radius + 1.*dy)

       do k=-2,k1+2
          coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
          do j=-2,j1+2
             coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
             do i=-2,i1+2
                coorx = leftbound + (i-0.5)/dxi ! [-dx*5/2 , lx+dx*5/2]

                !-- Define 2 spheres separately, differ only in center positions

                phi_drp1(i,j,k) = sqrt((coorx-xdp)**2+(coory-ydp1)**2+(coorz-zdp)**2) -radius
                phi_drp2(i,j,k) = sqrt((coorx-xdp)**2+(coory-ydp2)**2+(coorz-zdp)**2) -radius

                !-- Consider periodic BC in y direction

                phi_periodic1 = sqrt((coorx-xdp)**2+(coory-ydp1-ly)**2+(coorz-zdp)**2)-radius
                phi_periodic2 = sqrt((coorx-xdp)**2+(coory-ydp2-ly)**2+(coorz-zdp)**2)-radius

                if (phi_drp1(i,j,k) .gt. 0.) phi_drp1(i,j,k) = min(phi_drp1(i,j,k), phi_periodic1)
                if (phi_drp2(i,j,k) .gt. 0.) phi_drp2(i,j,k) = min(phi_drp2(i,j,k), phi_periodic2)

                !-- Consider periodic BC in x direction

                phi_periodic1 = sqrt((coorx-xdp-lx)**2+(coory-ydp1)**2+(coorz-zdp)**2)-radius
                phi_periodic2 = sqrt((coorx-xdp-lx)**2+(coory-ydp2)**2+(coorz-zdp)**2)-radius

                if (phi_drp1(i,j,k) .gt. 0.) phi_drp1(i,j,k) = min(phi_drp1(i,j,k), phi_periodic1)
                if (phi_drp2(i,j,k) .gt. 0.) phi_drp2(i,j,k) = min(phi_drp2(i,j,k), phi_periodic2)

                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       do k=-2,k1+2
          do j=-2,j1+2
             do i=-2,i1+2

                !-- Level set as min distance to interface

                if (phi_drp1(i,j,k) .le. 0.) ynew(i,j,k) = phi_drp1(i,j,k)
                if (phi_drp2(i,j,k) .le. 0.) ynew(i,j,k) = phi_drp2(i,j,k)

                if (phi_drp1(i,j,k) .gt. 0. .and. &
                     phi_drp2(i,j,k) .gt. 0.        ) then

                   ynew(i,j,k) = min( phi_drp1(i,j,k),phi_drp2(i,j,k) )

                endif

             enddo
          enddo
       enddo

       lvset(:,:,:,1) = ynew(:,:,:)

       ! Initialize collision velocity

       !     where (phi_drp1(0:i1,0:j1,0:k1) .le. 1.5*dy) vnew = .5
       !     where (phi_drp2(0:i1,0:j1,0:k1) .le. 1.5*dy) vnew =-.5     

  
    case('3dp')

       xdp = .5*lx
       zdp = .5*lz
       radius = 1./6.*lz
       ydp1= .1*ly
       ydp2= ydp1 + 2.*radius+radius/2.
       ydp3= ydp2 + 2.*radius+radius/2.

       do k=-2,k1+2
          coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
          do j=-2,j1+2
             coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
             do i=-2,i1+2
                coorx = leftbound + (i-0.5)/dxi ! [-dx*5/2 , lx+dx*5/2]

                !-- Define 3 spheres separately, differ only in center positions

                phi_drp1(i,j,k) = sqrt((coorx-xdp)**2+(coory-ydp1)**2+(coorz-zdp)**2) -radius
                phi_drp2(i,j,k) = sqrt((coorx-xdp)**2+(coory-ydp2)**2+(coorz-zdp)**2) -radius
                phi_drp3(i,j,k) = sqrt((coorx-xdp)**2+(coory-ydp3)**2+(coorz-zdp)**2) -radius

                !-- Consider periodic BC in y direction

                phi_periodic1 = sqrt((coorx-xdp)**2+(coory-ydp1-ly)**2+(coorz-zdp)**2)-radius
                phi_periodic2 = sqrt((coorx-xdp)**2+(coory-ydp2-ly)**2+(coorz-zdp)**2)-radius
                phi_periodic3 = sqrt((coorx-xdp)**2+(coory-ydp3-ly)**2+(coorz-zdp)**2)-radius

                if (phi_drp1(i,j,k) .gt. 0.) phi_drp1(i,j,k) = min(phi_drp1(i,j,k), phi_periodic1)
                if (phi_drp2(i,j,k) .gt. 0.) phi_drp2(i,j,k) = min(phi_drp2(i,j,k), phi_periodic2)
                if (phi_drp3(i,j,k) .gt. 0.) phi_drp3(i,j,k) = min(phi_drp3(i,j,k), phi_periodic3)

                !-- Consider periodic BC in x direction

                phi_periodic1 = sqrt((coorx-xdp-lx)**2+(coory-ydp1)**2+(coorz-zdp)**2)-radius
                phi_periodic2 = sqrt((coorx-xdp-lx)**2+(coory-ydp2)**2+(coorz-zdp)**2)-radius
                phi_periodic3 = sqrt((coorx-xdp-lx)**2+(coory-ydp3)**2+(coorz-zdp)**2)-radius

                if (phi_drp1(i,j,k) .gt. 0.) phi_drp1(i,j,k) = min(phi_drp1(i,j,k), phi_periodic1)
                if (phi_drp2(i,j,k) .gt. 0.) phi_drp2(i,j,k) = min(phi_drp2(i,j,k), phi_periodic2)
                if (phi_drp3(i,j,k) .gt. 0.) phi_drp3(i,j,k) = min(phi_drp3(i,j,k), phi_periodic3)

                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       do k=-2,k1+2
          do j=-2,j1+2
             do i=-2,i1+2

                !-- Level set as min distance to interface

                if (phi_drp1(i,j,k) .le. 0.) ynew(i,j,k) = phi_drp1(i,j,k)
                if (phi_drp2(i,j,k) .le. 0.) ynew(i,j,k) = phi_drp2(i,j,k)
                if (phi_drp3(i,j,k) .le. 0.) ynew(i,j,k) = phi_drp3(i,j,k)

                if (phi_drp1(i,j,k) .gt. 0. .and. &
                    phi_drp2(i,j,k) .gt. 0. .and. &
                    phi_drp3(i,j,k) .gt. 0.        ) then

                   ynew(i,j,k) = min( phi_drp1(i,j,k),phi_drp2(i,j,k),phi_drp3(i,j,k) )

                endif

                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       lvset(:,:,:,1) = ynew(:,:,:)

  
    case('mdp')

       !-- Center of each droplet --!

       radius_mls(1) = 0.5
       gap = 0. !2.0*dx
       
       xdp_mls(1) = 0.5*lx !-1.0*dx -0.5
       ydp_mls(1) = 1.5 !0.5*ly !-1.0*dy -0.5
       zdp_mls(1) = 0.3*lz !-1.0*dz -0.5

       allocate ( v_angle(3) )
       allocate ( h_angle(3) )
       v_angle = 80./180.*pi !pi/2.
       h_angle = (/0., 60., 75./) /180.*pi

       do l = 2,lmax

          radius_mls(l) = 0.5
          center_dist = radius_mls(l-1) + radius_mls(l) + gap
          
          xdp_mls(l) = xdp_mls(l-1) + center_dist*sin(v_angle(l))*cos(h_angle(l))
          ydp_mls(l) = ydp_mls(l-1) + center_dist*sin(v_angle(l))*sin(h_angle(l))
          zdp_mls(l) = zdp_mls(l-1) + center_dist*cos(v_angle(l))

       enddo

       !!## -- Up-Down
       !xdp_mls(2) = xdp_mls(1)
       !ydp_mls(2) = ydp_mls(1)
       !zdp_mls(2) = 0.5*lz +1.0*dz +0.5
       !!## --
       
       !!## -- tetrahedron
       !xdp_mls(4) = (xdp_mls(1)+xdp_mls(2)+xdp_mls(3))/3.
       !ydp_mls(4) = (ydp_mls(1)+ydp_mls(2)+ydp_mls(3))/3.
       !zdp_mls(4) = zdp_mls(3)+sqrt(2./3.)*center_dist       
       !!## -- polyhedron 5
       !xdp_mls(5) = (xdp_mls(1)+xdp_mls(2)+xdp_mls(3))/3.
       !ydp_mls(5) = (ydp_mls(1)+ydp_mls(2)+ydp_mls(3))/3.
       !zdp_mls(5) = zdp_mls(3)-sqrt(2./3.)*center_dist       
       !!## --
       

       do l = 1,lmax

          !-- Different level set for each drop

          do k=-2,k1+2
             coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
             do j=-2,j1+2
                coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
                do i=-2,i1+2
                   coorx = leftbound + (i-0.5)/dxi ! [-dx*5/2 , lx+dx*5/2]

                   !-- Define each drop separately

                   lvset(i,j,k,l) = sqrt( (coorx-xdp_mls(l))**2 + &
                        (coory-ydp_mls(l))**2 + &
                        (coorz-zdp_mls(l))**2   ) -radius_mls(l)

                   !-- Consider periodic BC in y direction

                   phi_periodic1 = sqrt((coorx-xdp_mls(l))**2    + &
                        (coory-ydp_mls(l)-ly)**2 + &
                        (coorz-zdp_mls(l))**2      ) -radius_mls(l)

                   if (lvset(i,j,k,l) .gt. 0.) lvset(i,j,k,l) = min(lvset(i,j,k,l),phi_periodic1)

                   !-- Consider periodic BC in x direction

                   phi_periodic2 = sqrt((coorx-xdp_mls(l)-lx)**2 + &
                        (coory-ydp_mls(l))**2    + &
                        (coorz-zdp_mls(l))**2      ) -radius_mls(l)

                   if (lvset(i,j,k,l) .gt. 0.) lvset(i,j,k,l) = min(lvset(i,j,k,l),phi_periodic2)

                enddo
             enddo
          enddo

       enddo

       ynew(:,:,:) =  lvset(:,:,:,1)

       
    case('pac') ! Rigid sphere packings of 6-10 drops (no gap)

       ! 6 drops
       
       xs = 0.6  ! x coordinate shift
       ys = 0.6  ! y coordinate shift
       zs = 0.6  ! z coordinate shift

       x_6pac(1) = xs + 0.0
       y_6pac(1) = ys + 0.0
       z_6pac(1) = zs + 0.0
       x_6pac(2) = xs + 0.5555555555555556
       y_6pac(2) = ys + 1.2830005981991683
       z_6pac(2) = zs + 0.9072184232530289
       x_6pac(3) = xs + 1.0
       y_6pac(3) = ys + 0.0
       z_6pac(3) = zs + 0.0
       x_6pac(4) = xs + 1.3333333333333333
       y_6pac(4) = ys + 0.7698003589195010
       z_6pac(4) = zs + 0.5443310539518174
       x_6pac(5) = xs + 0.5
       y_6pac(5) = ys + 0.8660254037844386
       z_6pac(5) = zs + 0.0
       x_6pac(6) = xs + 0.5
       y_6pac(6) = ys + 0.2886751345948129
       z_6pac(6) = zs + 0.8164965809277260
       
       ! 7 drops
       
       xs = 0.6  ! x coordinate shift
       ys = 0.6  ! y coordinate shift
       zs = 0.6  ! z coordinate shift

       x_7pac(1) = xs + 0.0
       y_7pac(1) = ys + 0.0
       z_7pac(1) = zs + 0.0
       x_7pac(2) = xs + 1.0925925925925926
       y_7pac(2) = ys + 0.6949586573578829
       z_7pac(2) = zs + 1.5120307054217148
       x_7pac(3) = xs + 1.0
       y_7pac(3) = ys + 0.0
       z_7pac(3) = zs + 0.0
       x_7pac(4) = xs + 0.5555555555555556
       y_7pac(4) = ys + 1.2830005981991683
       z_7pac(4) = zs + 0.9072184232530289
       x_7pac(5) = xs + 0.5
       y_7pac(5) = ys + 0.8660254037844386
       z_7pac(5) = zs + 0.0
       x_7pac(6) = xs + 1.3333333333333333
       y_7pac(6) = ys + 0.7698003589195010
       z_7pac(6) = zs + 0.5443310539518174
       x_7pac(7) = xs + 0.5
       y_7pac(7) = ys + 0.2886751345948129
       z_7pac(7) = zs + 0.8164965809277260

       ! 8 drops
       
       xs = 0.6  ! x coordinate shift
       ys = 0.9  ! y coordinate shift
       zs = 1.6  ! z coordinate shift

       x_8pac(1) = xs + 0.0
       y_8pac(1) = ys + 0.0
       z_8pac(1) = zs + 0.0
       x_8pac(2) = xs + 1.3888888888888888
       y_8pac(2) = ys + 0.8018753738744803
       z_8pac(2) = zs + 0.4536092116265145
       x_8pac(3) = xs + 1.3888888888888888
       y_8pac(3) = ys  -0.1603750747748960
       z_8pac(3) = zs  -0.9072184232530289
       x_8pac(4) = xs + 0.5555555555555556
       y_8pac(4) = ys + 1.2830005981991683
       z_8pac(4) = zs  -0.9072184232530289
       x_8pac(5) = xs + 1.0
       y_8pac(5) = ys + 0.0
       z_8pac(5) = zs + 0.0
       x_8pac(6) = xs + 0.5
       y_8pac(6) = ys + 0.8660254037844386
       z_8pac(6) = zs + 0.0
       x_8pac(7) = xs + 0.5
       y_8pac(7) = ys + 0.2886751345948129
       z_8pac(7) = zs  -0.8164965809277260
       x_8pac(8) = xs + 1.3333333333333333
       y_8pac(8) = ys + 0.7698003589195010
       z_8pac(8) = zs  -0.5443310539518174

       ! 9 drops
       
       xs = 1.1  ! x coordinate shift
       ys = 0.6  ! y coordinate shift
       zs = 0.6  ! z coordinate shift

       x_9pac(1) = xs + 0.0
       y_9pac(1) = ys + 0.0
       z_9pac(1) = zs + 0.0
       x_9pac(2) = xs + 1.0
       y_9pac(2) = ys + 0.0
       z_9pac(2) = zs + 0.0
       x_9pac(3) = xs  -0.5
       y_9pac(3) = ys + 0.8660254037844386
       z_9pac(3) = zs + 0.0
       x_9pac(4) = xs + 1.0
       y_9pac(4) = ys + 1.5396007178390021
       z_9pac(4) = zs + 0.5443310539518174
       x_9pac(5) = xs + 1.5
       y_9pac(5) = ys + 0.8660254037844386
       z_9pac(5) = zs + 0.0
       x_9pac(6) = xs + 0.0
       y_9pac(6) = ys + 1.5396007178390021
       z_9pac(6) = zs + 0.5443310539518174
       x_9pac(7) = xs + 0.0
       y_9pac(7) = ys + 0.5773502691896257
       z_9pac(7) = zs + 0.8164965809277260
       x_9pac(8) = xs + 1.0
       y_9pac(8) = ys + 0.5773502691896257
       z_9pac(8) = zs + 0.8164965809277260
       x_9pac(9) = xs + 0.5
       y_9pac(9) = ys + 0.8660254037844386
       z_9pac(9) = zs + 0.0

       ! 10 drops
       
       xs = 1.1  ! x coordinate shift
       ys = 0.6  ! y coordinate shift
       zs = 0.6  ! z coordinate shift

       x_10pac(1) = xs + 0.0
       y_10pac(1) = ys + 0.0
       z_10pac(1) = zs + 0.0
       x_10pac(2) = xs + 1.0
       y_10pac(2) = ys + 0.0
       z_10pac(2) = zs + 0.0
       x_10pac(3) = xs  -0.5
       y_10pac(3) = ys + 0.2886751345948129
       z_10pac(3) = zs + 0.8164965809277260
       x_10pac(4) = xs + 0.0
       y_10pac(4) = ys + 0.5773502691896257
       z_10pac(4) = zs + 1.6329931618554521
       x_10pac(5) = xs + 1.5
       y_10pac(5) = ys + 0.2886751345948129
       z_10pac(5) = zs + 0.8164965809277260
       x_10pac(6) = xs + 1.0
       y_10pac(6) = ys + 0.5773502691896257
       z_10pac(6) = zs + 1.6329931618554521
       x_10pac(7) = xs + 0.5
       y_10pac(7) = ys + 0.8660254037844386
       z_10pac(7) = zs + 0.0
       x_10pac(8) = xs + 0.0
       y_10pac(8) = ys + 1.1547005383792515
       z_10pac(8) = zs + 0.8164965809277260
       x_10pac(9) = xs + 1.0
       y_10pac(9) = ys + 1.1547005383792515
       z_10pac(9) = zs + 0.8164965809277260
       x_10pac(10)= xs + 0.5
       y_10pac(10)= ys + 0.2886751345948129
       z_10pac(10)= zs + 0.8164965809277260

       do l = 1,lmax

          if (lmax .eq. 6) then             
             xdp_mls(l) = x_6pac(l)
             ydp_mls(l) = y_6pac(l)
             zdp_mls(l) = z_6pac(l)

          elseif (lmax .eq. 7) then    
             xdp_mls(l) = x_7pac(l)
             ydp_mls(l) = y_7pac(l)
             zdp_mls(l) = z_7pac(l)

          elseif (lmax .eq. 8) then  
             xdp_mls(l) = x_8pac(l)
             ydp_mls(l) = y_8pac(l)
             zdp_mls(l) = z_8pac(l)

          elseif (lmax .eq. 9) then  
             xdp_mls(l) = x_9pac(l)
             ydp_mls(l) = y_9pac(l)
             zdp_mls(l) = z_9pac(l)

          elseif (lmax .eq. 10) then
             xdp_mls(l) = x_10pac(l)
             ydp_mls(l) = y_10pac(l)
             zdp_mls(l) = z_10pac(l)
          endif

          !-- Different level set for each drop

          do k=-2,k1+2
             coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
             do j=-2,j1+2
                coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
                do i=-2,i1+2
                   coorx = leftbound + (i-0.5)/dxi ! [-dx*5/2 , lx+dx*5/2]

                   !-- Define each drop separately

                   lvset(i,j,k,l) = sqrt( (coorx-xdp_mls(l))**2 + &
                        (coory-ydp_mls(l))**2 + &
                        (coorz-zdp_mls(l))**2   ) - 0.5

                   !-- Consider periodic BC in y direction

                   phi_periodic1 = sqrt((coorx-xdp_mls(l))**2    + &
                        (coory-ydp_mls(l)-ly)**2 + &
                        (coorz-zdp_mls(l))**2      ) - 0.5

                   if (lvset(i,j,k,l) .gt. 0.) lvset(i,j,k,l) = min(lvset(i,j,k,l),phi_periodic1)

                   !-- Consider periodic BC in x direction

                   phi_periodic2 = sqrt((coorx-xdp_mls(l)-lx)**2 + &
                        (coory-ydp_mls(l))**2    + &
                        (coorz-zdp_mls(l))**2      ) - 0.5

                   if (lvset(i,j,k,l) .gt. 0.) lvset(i,j,k,l) = min(lvset(i,j,k,l),phi_periodic2)

                enddo
             enddo
          enddo

       enddo

       ynew(:,:,:) =  lvset(:,:,:,1)


    case('elp')

       ! Ellipse (distorted)

       coef_a = 4.**2
       coef_b = 2.**2
       ypc = 3.5
       zpc = 2.
       coef_c = 0.1

       do k=-2,k1+2
          coorz = (k-0.5)/dzi - lz/2. ! [lz/2-dz*5/2 , lz/2+dz*5/2]
          do j=-2,j1+2
             coory = lowerbound + (j-0.5)/dyi - ly/2. ! [ly/2-dy*5/2 , ly/2+dy*5/2]
             do i=-2,i1+2

                !-- distortion factor
                coef_d = 1. + ((coory-3.5)**2 + (coorz-2.)**2)*0.05
                !coef_d = 1.
                ynew(i,j,k) = coef_d*(sqrt(coory**2/coef_a + coorz**2/coef_b) - 1.)

             enddo
          enddo
       enddo

       lvset(:,:,:,1) = ynew(:,:,:)
       reset = lvset


    case('zal') ! Zalesak's disk   

       ypc = ly*0.5
       zpc = lz*0.75
       radius = 0.15*ly

       do k=-2,k1+2
          coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
          do j=-2,j1+2
             coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
             do i=-2,i1+2

                phi_pancake(i,j,k) = sqrt((coory-ypc)**2+(coorz-zpc)**2) - radius  ! pancake

                phi_vc = min(abs(coory-0.475),abs(coory-0.525))  ! vertical cut
                phi_hc = min(abs(coorz-0.85),abs(coorz-0.6))   ! horizontal cut
                if (abs(coory-0.5) .le. 0.025 .and. abs(coorz-0.725) .le. .125) then  ! slot
                   phi_slot(i,j,k) = min(phi_vc,phi_hc)
                else
                   if (abs(coory-0.5) .le. 0.025) then
                      phi_slot(i,j,k) = -phi_hc
                   elseif (abs(coorz-0.725) .le. .125) then
                      phi_slot(i,j,k) = -phi_vc
                   else
                      phi_slot(i,j,k) = -sqrt(phi_vc**2+phi_hc**2)
                   endif
                endif

             enddo
          enddo
       enddo

       do k=-2,k1+2
          coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
          do j=-2,j1+2
             coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
             do i=-2,i1+2

                if (phi_pancake(i,j,k) .gt. 0.) ynew(i,j,k) = phi_pancake(i,j,k)
                if (phi_slot(i,j,k) .gt. 0.) ynew(i,j,k) = phi_slot(i,j,k)
                if (phi_pancake(i,j,k) .le. 0. .and. phi_slot(i,j,k) .le. 0.) then
                   ynew(i,j,k) = -min( abs(phi_pancake(i,j,k)),abs(phi_slot(i,j,k)) )
                endif

             enddo
          enddo
       enddo

       lvset(:,:,:,1) = ynew(:,:,:)
  
  
    case('2lc')

       beta = 0.2*lz  !# must be consistent with subroutine init
       do k=-2,k1+2
          coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
          do j=-2,j1+2
             coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
             do i=-2,i1+2

                !ynew(i,j,k) = coorz-beta ! Interface is a straight line at coorz = beta.
                ynew(i,j,k) = coorz - 0.47
                !ynew(i,j,k) = coory-25.
                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       lvset(:,:,:,1) = ynew(:,:,:)


    case('sin')

       alpha = 0.05  ! amplitude
       beta  = 1.00  ! period
       gamma = 0.50  ! mean value
       
       do k = -2,k1+2
          coorz = (k-0.5)/dzi
          do j = -2,j1+2
             coory = lowerbound + (j-0.5)/dyi
             do i = -2,i1+2

                ynew(i,j,k) = coorz - alpha*cos(coory*2.*pi/beta) - gamma
                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       lvset(:,:,:,1) = ynew(:,:,:)


    case('any')  ! anything (temporary testing)

       alpha = 0.05  ! amplitude
       beta  = 1.00  ! period
       gamma = 0.50  ! mean value
       
       do k = -2,k1+2
          coorz = (k-0.5)/dzi
          do j = -2,j1+2
             coory = lowerbound + (j-0.5)/dyi
             do i = -2,i1+2

                !ynew(i,j,k) = coorz - alpha*cos(coory*2.*pi/beta) - gamma  ! wave
                ynew(i,j,k) = 0.5-sqrt((coory-0.5)**2+(coorz-0.5-sqrt(3.)/4.)**2) ! circle
                if (coorz .gt. 1.) ynew(i,j,k) = NarrowBand_2
                phi_pancake(i,j,k) = coorz - gamma  ! straight line
                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       do k = -2,k1+2
          do j = -2,j1+2
             coory = lowerbound + (j-0.5)/dyi
             do i = -2,i1+2

                if (coory .lt. 0.25 .or. coory .gt. 0.75) then
                   ynew(i,j,k) = phi_pancake(i,j,k)
                endif

             enddo
          enddo
       enddo

       lvset(:,:,:,1) = ynew(:,:,:)


    case('hdm')  ! Hadamard's droplet

       ! Sphere centered at (xdp,ydp,zdp) with radius

       xdp = 0.5*lx
       ydp = 0.5*ly 
       zdp = 0.5*lz
       radius = 0.5       

       do k=-2,k1+2
          coorz = (k-0.5)*dz ! [-dz*5/2 , lz+dz*5/2]
          do j=-2,j1+2
             coory = lowerbound + (j-0.5)*dy ! [-dy*5/2 , ly+dy*5/2]
             do i=-2,i1+2
                coorx = leftbound + (i-0.5)*dx ! [-dx*5/2 , lx+dx*5/2]

                ynew(i,j,k) = sqrt((coorx-xdp)**2+(coory-ydp)**2+(coorz-zdp)**2) - radius

                !-- Consider periodic BC in y direction

                phi_periodic1 = sqrt((coorx-xdp)**2+(coory-ydp-ly)**2+(coorz-zdp)**2)-radius

                if (ynew(i,j,k) .gt. 0.) ynew(i,j,k) = min(ynew(i,j,k), phi_periodic1)

                !-- Consider periodic BC in x direction

                phi_periodic1 = sqrt((coorx-xdp-lx)**2+(coory-ydp)**2+(coorz-zdp)**2)-radius

                if (ynew(i,j,k) .gt. 0.) ynew(i,j,k) = min(ynew(i,j,k), phi_periodic1)

                dydt(i,j,k) = 0.

             enddo
          enddo
       enddo

       lvset(:,:,:,1) = ynew(:,:,:)

       ! Re-initialize the velocity field based on Hadamard's solution

       coef_d = -(rho2-rho1)*9.81/(12.*miu1+18.*miu2) ! original
       !coef_d = -1./(4.*(1.+3.*miu2/miu1)*radius**2) ! normalized s.t. O(vel)~1
       coef_c = -(3.+2.*miu2/miu1)*coef_d*radius**2
       coef_b = -(2.+3.*miu2/miu1)*coef_d*radius**3
       coef_a = miu2/miu1*coef_d*radius**5

       do k = 0,k1
          coorz = (k-0.5)*dz
          do j = 0,j1
             coory = lowerbound + (j-0.5)*dy
             do i = 0,i1
                coorx = leftbound + (i-0.5)*dx

                xs = coorx-xdp
                ys = coory-ydp
                zs = coorz-zdp
                rs = lvset(i,j,k,1) + radius

                cos_the = ys/rs
                sin_the = sqrt(1.-cos_the**2)  ! always positive
                sin_phi = xs/sqrt(xs**2+zs**2)
                cos_phi = zs/sqrt(xs**2+zs**2)
                
                if (rs .gt. radius) then  ! outside

                   vel_r = 2.*(coef_a/rs**3 + coef_b/rs)*cos_the
                   vel_t =    (coef_a/rs**3 - coef_b/rs)*sin_the

                else  ! inside

                   vel_r = 2.*(coef_c +    coef_d*rs**2)*cos_the
                   vel_t =-2.*(coef_c + 2.*coef_d*rs**2)*sin_the

                endif

                ! cell-center velocities in the droplet frame (not a bad init for cell-face)

                unew(i,j,k) = (vel_r*sin_the + vel_t*cos_the)*sin_phi
                vnew(i,j,k) =  vel_r*cos_the - vel_t*sin_the - 2.*(coef_c + coef_d*radius**2)
                wnew(i,j,k) = (vel_r*sin_the + vel_t*cos_the)*cos_phi

             enddo
          enddo
       enddo

    end select


    return
  end subroutine init_levelset

  
  !~~~~~~~~~~~~~~~~~~~~~~~~
  !
  subroutine init_static_ls
  !                                                             
  ! Initialize a static level set that is not advected later on.
  ! It might be used as a reference or for visualization.       
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    integer :: i,j,k,l
    integer :: l_pil
    parameter (l_pil = 2)
    
    real :: lowerbound,upperbound,leftbound,rightbound    
    real :: ypc,zpc, radius
    real :: coory,coorz
    real, dimension(l_pil) :: ypil,zpil
    real :: phi_vc,phi_hc, w
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2,l_pil) :: phi_pil
    real :: eps,phi
    
    character(len=3) :: static_shape
    parameter (static_shape = 'pil')


    !---- Bounds of Each Core ----!
    
    lowerbound =  coords(2)   *ly/(dims(2)*1.)
    upperbound = (coords(2)+1)*ly/(dims(2)*1.)
    leftbound  =  coords(1)   *lx/(dims(1)*1.)
    rightbound = (coords(1)+1)*lx/(dims(1)*1.)

    !---- Define the static level set ----!
    
    select case(static_shape)

    case('pan')

       ypc = ly*0.5
       zpc = lz*0.5
       radius = 0.5

       do k = -2,k1+2
          coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
          do j = -2,j1+2
             coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
             do i = -2,i1+2

                fixls(i,j,k) = sqrt((coory-ypc)**2+(coorz-zpc)**2) - radius

             enddo
          enddo
       enddo

    case('sin')

       do k = -2,k1+2
          coorz = (k-0.5)/dzi
          do j = -2,j1+2
             coory = lowerbound + (j-0.5)/dyi
             do i = -2,i1+2

                fixls(i,j,k) = coorz - 0.1*cos(coory*2.*pi)-lz/8.

             enddo
          enddo
       enddo       
       
    case('pil') ! pillars
       
       do l = 1,l_pil
          ypil(l) = 0. +(l-1)*ly  ! streamwise locations
          zpil(l) = 0.5  ! height of the pillar
       enddo

       w = 0.5 ! width of the pillar

       !---- define each pillar ----!
       
       do l = 1,l_pil
          do k = -2,k1+2
             coorz = (k-0.5)/dzi ! [-dz*5/2 , lz+dz*5/2]
             do j = -2,j1+2
                coory = lowerbound + (j-0.5)/dyi ! [-dy*5/2 , ly+dy*5/2]
                do i = -2,i1+2

                   phi_vc = min(abs(coory-(ypil(l)-w/2.)),abs(coory-(ypil(l)+w/2.)))  ! vertical cut
                   phi_hc = min(abs(coorz-zpil(l)),abs(coorz))   ! horizontal cut
                   if (abs(coory-ypil(l)) .le. w/2. .and. coorz .le. zpil(l)) then  ! pillar
                      phi_pil(i,j,k,l) = -min(phi_vc,phi_hc)
                   else
                      if (abs(coory-ypil(l)) .le. w/2.) then
                         phi_pil(i,j,k,l) = phi_hc
                      elseif (coorz .le. zpil(l)) then
                         phi_pil(i,j,k,l) = phi_vc
                      else
                         phi_pil(i,j,k,l) = sqrt(phi_vc**2+phi_hc**2)
                      endif
                   endif

                enddo
             enddo
          enddo
       enddo
       
       fixls = phi_pil(:,:,:,1)
       do l=2,l_pil
          fixls = min(fixls, phi_pil(:,:,:,l))
       enddo

    end select

    !---- Corresponding Heaviside for visualization ----!

    fixhv = Heaviside('4','trig', 1,1.5*dz)  ! the 3rd entry doesn't matter

    
    return
  end subroutine init_static_ls

  
  !~~~~~~~~~~~~~~~~~~~~~
  !
  subroutine init_matrix
  !                                                                          
  ! Pre-compute a least-squares matrix (pseudo_inv_A) for level set/curvature
  ! Ax = b --> x = (pseudo_inv_A)b                                           
  ! Can result in ~90% speedup in get/re_curvature.f90                       
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    integer :: dist_x, dist_y, dist_z, counter
    real, dimension(7,4) :: Al ! linear coef. matrix
    real, dimension(4,4) :: ATAl, ATAl_inv
    real, dimension(27,10) :: A, Ab ! quadratic coef. matrix
    real, dimension(10,10) :: ATA, ATA_inv, ATAb,ATAb_inv

 
    !-- Define matrix Al --!
    
    Al(:,:) = 0.
    Al(:,1) = 1.
    Al(3,2) = -dx
    Al(5,2) =  dx
    Al(2,3) = -dy
    Al(6,3) =  dy
    Al(1,4) = -dz
    Al(7,4) =  dz

    !-- Compute the pseudo inverse of Al --!

    ATAl = matmul(transpose(Al), Al)
    ATAl_inv = inv(ATAl)
    pseudo_inv_Al = matmul(ATAl_inv, transpose(Al))  ! stored in common.f90


    !-- Define matrix A --!

    A(:,1)  = 1.

    counter = 1  ! line counter of matrix A
    do dist_z = -1,1
       do dist_y = -1,1
          do dist_x = -1,1                 

             A(counter,2) = dx*dist_x
             A(counter,3) = dy*dist_y
             A(counter,4) = dz*dist_z

             counter = counter + 1

          enddo
       enddo
    enddo

    A(:,5) = 0.5 * A(:,2)**2
    A(:,6) = 0.5 * A(:,3)**2
    A(:,7) = 0.5 * A(:,4)**2

    A(:,8)  = A(:,2)*A(:,3)
    A(:,9)  = A(:,2)*A(:,4)
    A(:,10) = A(:,3)*A(:,4)

    !-- Compute the pseudo inverse of A --!

    ATA = matmul(transpose(A), A)
    ATA_inv= inv(ATA)
    pseudo_inv_A = matmul(ATA_inv, transpose(A))  ! stored in common.f90

    
    !-- Define matrix Ab (sided stencil) --!

    ! lower side
    
    Ab(:,1)  = 1.

    counter = 1  ! line counter of matrix Ab
    do dist_z = 0,2
       do dist_y = -1,1
          do dist_x = -1,1                 

             Ab(counter,2) = dx*dist_x
             Ab(counter,3) = dy*dist_y
             Ab(counter,4) = dz*dist_z

             counter = counter + 1

          enddo
       enddo
    enddo

    Ab(:,5) = 0.5 * Ab(:,2)**2
    Ab(:,6) = 0.5 * Ab(:,3)**2
    Ab(:,7) = 0.5 * Ab(:,4)**2

    Ab(:,8)  = Ab(:,2)*Ab(:,3)
    Ab(:,9)  = Ab(:,2)*Ab(:,4)
    Ab(:,10) = Ab(:,3)*Ab(:,4)

    !-- Compute the pseudo inverse of Ab --!

    ATAb = matmul(transpose(Ab), Ab)
    ATAb_inv= inv(ATAb)
    pseudo_inv_Ab_lower = matmul(ATAb_inv, transpose(Ab))  ! stored in common.f90

    ! upper side
    
    Ab(:,1)  = 1.

    counter = 1  ! line counter of matrix Ab
    do dist_z = -2,0  
       do dist_y = -1,1
          do dist_x = -1,1                 

             Ab(counter,2) = dx*dist_x
             Ab(counter,3) = dy*dist_y
             Ab(counter,4) = dz*dist_z

             counter = counter + 1

          enddo
       enddo
    enddo

    Ab(:,5) = 0.5 * Ab(:,2)**2
    Ab(:,6) = 0.5 * Ab(:,3)**2
    Ab(:,7) = 0.5 * Ab(:,4)**2

    Ab(:,8)  = Ab(:,2)*Ab(:,3)
    Ab(:,9)  = Ab(:,2)*Ab(:,4)
    Ab(:,10) = Ab(:,3)*Ab(:,4)

    !-- Compute the pseudo inverse of Ab --!

    ATAb = matmul(transpose(Ab), Ab)
    ATAb_inv= inv(ATAb)
    pseudo_inv_Ab_upper = matmul(ATAb_inv, transpose(Ab))  ! stored in common.f90


    return 
  end subroutine init_matrix


end module mod_init
