module mod_misc

  use mod_common
  use mod_param
  use mod_param_ls

  implicit none
  
  private

  public write_log_header, update_materials, Heaviside,Dirac, &
         inv, average_vel_at_center,max_vel_at_center, &
         average_stream_vel, check_steady_state, &
         serpentine, Prosperetti, spurious_currents, lennardjones, Grace, Taylor, &
         ls_bound_error

contains

  
  !-------------
  !-------------

  
  subroutine write_log_header

    use mod_common_mpi
    use mod_chkdt

    real :: dtmax

    call chkdt(dtmax)

    if (myid .eq. 0) then

       write(6,'(A)') '****************** Launching a new simulation ******************'
       write(6,'(A)') ' '
       write(6,'(A,1F5.1,A,1F5.1,A,1F5.1)') 'Domain size = ',lx,' x ',ly,' x ',lz
       write(6,'(A,1I5,A,1I5,A,1I5)')       'Resolution  = ',itot,' x ',jtot,' x ',ktot
       write(6,'(A,1I5)') 'MPI process(es) = ',dims(1)*dims(2)
       write(6,'(A)') ' '
       write(6,'(A,1ES10.4)') 'Initial time step computed by chkdt = ', dtmax
       write(6,'(A,1ES10.4)') 'Actual time step set in param = ', timestep
       write(6,'(A)') ' '
       write(6,'(A,1ES10.4,A,1ES10.4)') 'Reference phase viscosity = ',miu1,', density = ',rho1
       write(6,'(A,1ES10.4,A,1ES10.4)') 'Second phase viscosity    = ',miu2,', density = ',rho2
       write(6,'(A)') ' '
       write(6,'(A,1ES8.2,A,1ES8.2,A,1ES8.2)') 'Re = ',Reb/lz,', We = ',We,', Ca = ',We*visc
       write(6,'(A)') ' '
       write(6,'(A,1I4)') 'Number of level set(s) = ',lmax
       write(6,'(A)') ' '
       write(6,'(A)') 'Numerical schemes/methods:'
       write(6,'(A,A)') '  NS_conv = ',NS_conv
       write(6,'(A,A)') '  NS_diff_discretization = ',NS_diff_discretization
       write(6,'(A,A)') '  P_jump  = ',surface_tension_method
       write(6,'(A,A)') '  LS_time = ',ls_adv_time
       write(6,'(A,A)') '  LS_adv  = ',ls_adv_space
       if (re_initialization) write(6,'(A,1I6)') '  LS_reinit_step = ', reinit_step
       if (mass_correction) then
          write(6,'(A,1I6)') '  Inflate_step   = ',inflate_step
          write(6,'(A,A)') '  Speed dependence = ',ispeed_dependence
          write(6,'(A,A)') '  Inflate_adv      = ',i_space
       end if
       write(6,'(A)') ' '
       write(6,'(A)') '******************* Running the main program *******************'
       write(6,'(A)') ' '
       write(6,'(A)') '! Initializing ... '
       write(6,'(A)') ' '

    endif
    
    
    return
  end subroutine write_log_header
   

  !-------------
  !-------------

  subroutine update_materials

    use mod_common_mpi
                                                                   
    ! Update material properties based on their Heaviside functions.           
    ! ------                                                        
    ! hnew_mls: Heaviside function of level set                     
    ! miu,rho : Viscosity, density at cell centers                  
    ! volume  : Volume fraction of phase 2 (eg. droplets or bubbles)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  end subroutine update_materials
   

  !-------------
  !-------------

  
  function Heaviside(var,shape, l,eps)  result(H)

    ! Compute discrete Heaviside function

    character(len=1), intent(in) :: var   ! independent variable
    character(len=4), intent(in) :: shape ! shape of the function
    integer, intent(in) :: l              ! multiple level set marker
    real, intent(in) :: eps               ! half spreading width

    integer :: i,j,k
    real :: phi, dys
    real, dimension(0:i1,0:j1,0:k1) :: buffer, H
    logical :: dynamic_support
    parameter (dynamic_support = .true.)


    select case (var)
       
    case('0') ! cell-center
       
       buffer(:,:,:) = lvset(0:i1,0:j1,0:k1,l)  

    case('1') ! u-grid
       
       do k = 0,k1
          do j = 0,j1
             do i = 0,i1

                buffer(i,j,k) = (lvset(i,j,k,l) +lvset(i+1,j,k,l))/2.

             enddo
          enddo
       enddo

    case('2') ! v-grid
       
       do k = 0,k1
          do j = 0,j1
             do i = 0,i1

                buffer(i,j,k) = (lvset(i,j,k,l) +lvset(i,j+1,k,l))/2.

             enddo
          enddo
       enddo
       
    case('3') ! w-grid
       
       do k = 0,k1
          do j = 0,j1
             do i = 0,i1

                buffer(i,j,k) = (lvset(i,j,k,l) +lvset(i,j,k+1,l))/2.

             enddo
          enddo
       enddo
       
    case('4') ! cell-center using fixls
       
       buffer(:,:,:) = fixls(0:i1,0:j1,0:k1)  
       
    end select
    
    select case (shape)
    case('jump') ! discontinuous
       
       do k = 0,k1
          do j = 0,j1
             do i = 0,i1

                phi = buffer(i,j,k)  

                if (phi .lt. 0.) then
                   H(i,j,k) = 0.
                else
                   H(i,j,k) = 1.
                endif

             enddo
          enddo
       enddo
       
    case('line')  ! linear

       dys = eps

       do k = 0,k1
          do j = 0,j1
             do i = 0,i1

                phi = buffer(i,j,k)  

                if (dynamic_support) dys = dys*(abs(normal(i,j,k)%x) &
                                              + abs(normal(i,j,k)%y) &
                                              + abs(normal(i,j,k)%z))

                if (phi .lt. -dys) then
                   H(i,j,k) = 0.

                elseif (phi .lt. dys) then
                   H(i,j,k) = .5 + .5*phi/dys

                else
                   H(i,j,k) = 1.
                endif

             enddo
          enddo
       enddo

    case('para')  ! parabolic

       dys = eps

       do k = 0,k1
          do j = 0,j1
             do i = 0,i1

                phi = buffer(i,j,k)  

                if (dynamic_support) dys = dys*(abs(normal(i,j,k)%x) &
                     + abs(normal(i,j,k)%y) &
                     + abs(normal(i,j,k)%z))

                if (phi .lt. -dys) then
                   H(i,j,k) = 0.

                elseif (phi .lt. 0.) then
                   H(i,j,k) = .5 + phi/dys + .5*phi**2/dys**2

                elseif (phi .lt. dys) then
                   H(i,j,k) = .5 + phi/dys - .5*phi**2/dys**2

                else
                   H(i,j,k) = 1.
                endif

             enddo
          enddo
       enddo

    case('trig')  ! trigonometric

       do k = 0,k1
          do j = 0,j1
             do i = 0,i1

                phi = buffer(i,j,k)  

                if (phi .lt. -eps) then
                   H(i,j,k) = 0.

                elseif (phi .lt. eps) then
                   H(i,j,k) = .5 + .5*phi/eps + .5/pi*sin(pi*phi/eps)

                else
                   H(i,j,k) = 1.
                endif

             enddo
          enddo
       enddo

    end select

    return
  end function Heaviside


  !-------------
  !-------------

  function Dirac(flag, l,eps)  result(del)

    ! Compute discrete Heaviside function

    character(len=4), intent(in) :: flag  ! different discretizations
    integer, intent(in) :: l              ! multiple level set marker
    real, intent(in) :: eps               ! half spreading width

    integer :: i,j,k
    real :: phi, dys
    real, dimension(1:imax,1:jmax,1:kmax) :: del
    logical :: dynamic_support
    parameter (dynamic_support = .false.)


    select case (flag)

    case('line')  ! linear

       dys = eps

       do k = 1,kmax
          do j = 1,jmax
             do i = 1,imax

                phi = lvset(i,j,k,l)  

                if (dynamic_support) dys = dys*( abs(normal(i,j,k)%x) &
                     + abs(normal(i,j,k)%y) &
                     + abs(normal(i,j,k)%z) )

                if (phi .lt. -dys) then
                   del(i,j,k) = 0.

                elseif (phi .lt. 0.) then
                   del(i,j,k) = 1./dys + phi/dys**2

                elseif (phi .lt. dys) then
                   del(i,j,k) = 1./dys - phi/dys**2

                else
                   del(i,j,k) = 0.
                endif

             enddo
          enddo
       enddo

    case('trig')  ! trig

       do k = 1,kmax
          do j = 1,jmax
             do i = 1,imax

                phi = lvset(i,j,k,l)  

                if (phi .lt. -eps) then
                   del(i,j,k) = 0.

                elseif (phi .lt. eps) then
                   del(i,j,k) = .5/eps + .5/eps*cos(pi*phi/eps)

                else
                   del(i,j,k) = 0.
                endif

             enddo
          enddo
       enddo

    end select

    return
  end function Dirac


  !-------------
  !-------------

  function inv(A)  result(Ainv)

    ! Compute matrix inverse

    real, dimension(:,:), intent(in) :: A
    real, dimension(size(A,1),size(A,2)) :: Ainv

    real, dimension(size(A,1)) :: work
    integer, dimension(size(A,1)) :: ipiv
    integer :: n, info

    external dgetrf
    external dgetri

    Ainv = A
    n = size(A,1)

    call dgetrf(n, n, Ainv, n, ipiv, info)

    if (info/=0) then
       write(6,*) 'Lapack matrix is singular.'
       stop
    endif

    call dgetri(n, Ainv, n, ipiv, work, n, info)

    if (info/=0) then
       write(6,*) 'Lapack failed to inverse the matrix.'
       stop
    endif

    return
  end function inv


  !-------------
  !-------------

  subroutine average_vel_at_center

    integer :: i,j,k

    do k= 1,kmax
       do j= 1,jmax
          do i= 1,imax

             uccnew(i,j,k) = (unew(i-1,j,k)+unew(i,j,k))/2.
             vccnew(i,j,k) = (vnew(i,j-1,k)+vnew(i,j,k))/2.
             wccnew(i,j,k) = (wnew(i,j,k-1)+wnew(i,j,k))/2.

          enddo
       enddo
    enddo

    return
  end subroutine average_vel_at_center


  !-------------
  !-------------

  subroutine max_vel_at_center(vel_mag_max_all)

    use mpi
    use mod_common_mpi

    real, intent(out) :: vel_mag_max_all

    integer :: i,j,k
    real :: vel_mag_max
    real,  dimension(1:imax,1:jmax,1:kmax) :: vel_mag

    do k= 1,kmax
       do j= 1,jmax
          do i= 1,imax

             vel_mag(i,j,k) = sqrt(uccnew(i,j,k)**2 + &
                                   vccnew(i,j,k)**2 + & 
                                   wccnew(i,j,k)**2 )
          enddo
       enddo
    enddo

    vel_mag_max = maxval(vel_mag)
    call mpi_allreduce(vel_mag_max,vel_mag_max_all,1,mpi_real8,mpi_max,comm_cart,error)

    return
  end subroutine max_vel_at_center

  
  !-------------
  !-------------

  subroutine average_stream_vel(v_avr)

    use mpi
    use mod_common_mpi

    real, intent(inout) :: v_avr

    
    v_avr = sum(vnew(1:imax,1:jmax,1:kmax))
    call mpi_allreduce(mpi_in_place, v_avr ,1,mpi_real8,mpi_sum,comm_cart,error) 
    v_avr = v_avr/(itot*jtot*kmax)  ! average streamwise velocity

    return
  end subroutine average_stream_vel

  
  !-------------
  !-------------

  subroutine check_steady_state(istep)
    ! Determine if reaching the steady state by checking the pressure gradient correction
    ! (in the correc.f90) which should vanish (ie. u,v,w,p all become steady).
    !------------------------------------------------------------------------------------

    use mpi
    use mod_common_mpi
    use mod_loadflds
    use mod_write_vtk

    integer, intent(in) :: istep
    integer :: i,j,k
    integer :: flag, flag_u,flag_v,flag_w
    
    real :: L2_u,L2_v,L2_w
    real, parameter :: tolerance = dz**4  ! second-order accuracy in L2

    flag_u = 0
    flag_v = 0
    flag_w = 0
    
    L2_u = 0.
    L2_v = 0.
    L2_w = 0.
    
    do k = 1,kmax
       do j = 1,jmax
          do i = 1,imax

             L2_u = L2_u + (unew(i,j,k) - dudt(i,j,k))**2
             L2_v = L2_v + (vnew(i,j,k) - dvdt(i,j,k))**2
             L2_w = L2_w + (wnew(i,j,k) - dwdt(i,j,k))**2

          enddo
       enddo
    enddo

    call mpi_allreduce(mpi_in_place, L2_u ,1,mpi_real8,mpi_sum,comm_cart,error)
    call mpi_allreduce(mpi_in_place, L2_v ,1,mpi_real8,mpi_sum,comm_cart,error)
    call mpi_allreduce(mpi_in_place, L2_w ,1,mpi_real8,mpi_sum,comm_cart,error)

    L2_u = L2_u/(itot*jtot*kmax)
    L2_v = L2_v/(itot*jtot*kmax)
    L2_w = L2_w/(itot*jtot*kmax)

    if (L2_u .lt. tolerance) flag_u = 1
    if (L2_v .lt. tolerance) flag_v = 1
    if (L2_w .lt. tolerance) flag_w = 1

    flag = flag_u*flag_v*flag_w

    if (flag .eq. 0) then
       
       if (mod(istep,ioutchk) .eq. 0) then
          if (myid .eq. 0) then
             write(6,'(A)') '! Velocity increment in the L2 norm:'
             write(6,'(A,1ES10.4,A,1ES10.4,A,1ES10.4)') '! L2_u = ',L2_u, ', L2_v = ',L2_v,' , L2_w = ',L2_w
          endif
       endif
       
    else
       
       if(myid .eq. 0) then
          write(6,'(A)') '!---------------------------------------------'
          write(6,'(A)') '! Reached steady state. Saving output files...'
       endif
       call loadflds(1,istep)
       call write_vtk(istep)
       if(myid .eq. 0) then
          write(6,*)
          write(6,'(A)') '************** Simulation completed. Congrats! \o/ *************'
       endif
       call decomp_2d_finalize
       call MPI_FINALIZE(error)
       stop
       
    endif
    
    return
  end subroutine check_steady_state

  
  !-------------
  !-------------

  subroutine serpentine(time)

    use mpi
    use mod_common_mpi
    
    integer i,j,k
    real period, coorz, coory
    real lowerbound
    real, intent(in) :: time

    lowerbound =  coords(2)*ly/(dims(2)*1.)

    !---- Periodic Vortex test case ----!

    period = 8.

    do k = 1,kmax
       coorz = (k-0.5)*dz/lz !# normalised channel height [0,1]
       do j = 1,jmax
          coory =  lowerbound + (j-0.5)*dy/lz !# normalised channel height [0,1]
          do i = 1,imax
             uccnew(i,j,k) = 0.
             vccnew(i,j,k) = -2. * (sin(pi*coory))**2 * sin(pi*coorz)*cos(pi*coorz)*cos(pi*time/period)
             wccnew(i,j,k) =  2. * sin(pi*coory)*cos(pi*coory) * (sin(pi*coorz))**2*cos(pi*time/period)     
          enddo
       enddo
    enddo


    return
  end subroutine serpentine


  !-------------
  !-------------

  subroutine Prosperetti(wave_amp)

    use mod_common_mpi

    !-- Compute the amplitude of a cosine wave at y=0
    !-- to verify w/ analytical solution in Prosperetti 1981 Phys. Fluids.
    !-- Note that only parallelization in y direction is allowed.

    real, intent(out) :: wave_amp

    if (myid .eq. 0) wave_amp = -(lvset(2,0,kmax/2,1)+lvset(2,1,kmax/2,1)+dz)/2.


    return
  end subroutine Prosperetti


  !-------------
  !-------------

  subroutine spurious_currents(Ca)

    use mod_common_mpi

    !-- Compute the spurious currents in a static droplet test
    !-- Desjardins et al JCP2008, Herrmann JCP2008

    real, intent(out) :: Ca
    integer :: i,j,k
    real :: mu,sigma
    real :: u,v,w, vel_max,vel_max_all
    real, dimension(1:imax,1:jmax,1:kmax) :: vel

    mu = 0.1
    sigma = 1.

    do k = 1,kmax
       do j = 1,jmax
          do i = 1,imax

             !-- Velocity magnitude at cell center --!

             u = (unew(i,j,k) + unew(i-1,j,k))/2.
             v = (vnew(i,j,k) + vnew(i,j-1,k))/2.
             w = (wnew(i,j,k) + wnew(i,j,k-1))/2.

             vel(i,j,k) = sqrt(u**2 + v**2 + w**2)

          enddo
       enddo
    enddo

    vel_max = maxval(vel)
    call mpi_allreduce(vel_max,vel_max_all,1,mpi_real8,mpi_max,comm_cart,error)

    !-- Capillary number as a measure of the spurious currents --!

    Ca = mu*vel_max_all/sigma


    return
  end subroutine spurious_currents



  !-------------
  !-------------

  subroutine lennardjones(r,pot)

    !---- Compute LJ potential:
    ! 
    !                   r(eq)^p     p   r(eq)^q
    !     pot = epl ( ---------- - --- ---------- )
    !                  r(m,n)^p     q   r(m,n)^q
    !     where
    !          epl    : strength of LJ potential
    !          r_eq   : equilibrium center-to-center distance
    !          r(m,n) : center distance between m and n
    !          r_cof  : cut-off center-to-center distance

    integer p,q
    real epl, pq, radius,r_eq,r_mn,r_cof,pot_cof
    real, intent(in) :: r
    real, intent(out) :: pot


    p = 18
    q = 6
    epl = 3.
    pq = p*1./(q*1.)

    radius = lz/6.
    r_eq = 2.*radius + 3.5*dz
    r_mn = 2.*radius + r
    r_cof = 2.*radius + 6.*dz

    pot_cof = epl*( r_eq**p/r_cof**p - pq* r_eq**q/r_cof**q) !potential at cut-off
    pot = epl*( r_eq**p/r_mn**p - pq* r_eq**q/r_mn**q) - pot_cof !make potential zero at cut-off
    if (r_mn .gt. r_cof) pot = 0. !make potential zero after cut-off


    return
  end subroutine lennardjones


  !-------------
  !-------------


  subroutine Grace(Re_comp)

    !-- Compute the terminal Reynolds number of a rising bubble
    !-- Original experiment by Grace.

    real, intent(out) :: Re_comp

    Re_comp = Reb/lz*v_disp(1)
    ! Reb/lz gives bubble Reynolds
    ! times a non-dimensional average velocity in y-dir (the rising direction)


    return
  end subroutine Grace


  !-------------
  !-------------


  subroutine Taylor(D)

    use mod_common_mpi

    !-- Compute the deformation parameter defined by G. I. Taylor.
    !-- Applied to a single droplet in a simple shear flow.

    real, intent(out) :: D

    real :: cutoff, d1,d1_all, d2,d2_all
    real :: length,breadth
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phim,phin,super_phi


    cutoff = NarrowBand_2  ! max distance considered

    phim =  fixls  ! static level set
    phin = -lvset(:,:,:,1)  ! droplet (negative since they overlap)

    super_phi = 0.

    where( abs(phim) .lt.  1.*dz .and. &
           abs(phin) .lt. cutoff)

       super_phi = phim + phin

    end where

    d1 = minval(super_phi)
    call mpi_allreduce(d1,d1_all,1,mpi_real8,mpi_min,comm_cart,error)
    breadth = 0.5 + d1_all  ! 0.5 is the radius of the reference circle

    d2 = maxval(super_phi)
    call mpi_allreduce(d2,d2_all,1,mpi_real8,mpi_max,comm_cart,error)
    length = 0.5 + d2_all  ! 0.5 is the radius of the reference circle

    D = (length-breadth)/(length+breadth)


    return
  end subroutine Taylor


  !-------------
  !-------------
  

  subroutine ls_bound_error

    use mod_common_mpi

    integer :: i,j,k, km
    real :: phi_bound,phi_exact, Linf_max
    real :: lowerbound, y_contact, phi,phip,phim
    real :: n_wall_y,n_wall_z,n_contact_y,n_contact_z, contact_angle
    real, dimension(0:j1) :: phi_contact
    real, dimension(imax,jmax) :: Linf

    !-- Compare ghost values with the exact
    
    Linf = 0.

    k = 0  ! 0 or k1
    do j = 1,jmax
       do i = 1,imax

          phi_bound = lvset(i,j,k,1)
          phi_exact = ynew(i,j,k)

          if (abs(phi_bound) .lt. NarrowBand_2-dx) then
             Linf(i,j) = abs(phi_bound - phi_exact)
          endif
          
       enddo
    enddo

    Linf_max = maxval(Linf)
    call mpi_allreduce(mpi_in_place,Linf_max,1,mpi_real8,mpi_max,comm_cart,error)

    if (myid .eq. 0) then
       write(6,*) '~~~~~~~~'
       write(6,'(A,1ES16.8)') 'L infinite of boundary level set = ',Linf_max
    endif

    !-- Compare the contact point position with the exact

    lowerbound = coords(2)*ly/(dims(2)*1.)
    
    k = 1
    km = k-1
    i = itot/2
    do j = 0,j1
       phi_contact(j) = (lvset(i,j,k,1)+lvset(i,j,km,1))/2.
    enddo
    y_contact = 0.
    do j = 1,jmax
       phi = phi_contact(j)
       phip= phi_contact(j+1)
       if (phi*phip .lt. 0.) then
          y_contact = lowerbound+(j-0.5)*dy + dy*abs(phi)/(abs(phi)+abs(phip))  ! global position
       endif
    enddo
    call mpi_allreduce(mpi_in_place,y_contact,1,mpi_real8,mpi_max,comm_cart,error)  ! update all procs
    if (myid .eq. 0) then
       write(6,'(A,1ES16.8)') 'Contact point error = ',y_contact-0.5
    endif

    !-- Calculate the contact angle

    n_wall_y = 0.
    n_wall_z = -1.
    n_contact_y = -1e8 ! initialize with a small value
    n_contact_z = -1e8
    
    k = 1
    km = k-1
    i = itot/2
    do j = 1,jmax
       phi = phi_contact(j)
       phip= phi_contact(j+1)
       if (phi*phip .lt. 0.) then
          phim = phi_contact(j-1)
          n_contact_y = (phip-phim)/(2.*dy)  ! doesn't matter for a horizontal wall
          n_contact_z = (lvset(i,j,k,1)-lvset(i,j,km,1))/dz
       endif
    enddo
    call mpi_allreduce(mpi_in_place,n_contact_y,1,mpi_real8,mpi_max,comm_cart,error)  ! update all procs
    call mpi_allreduce(mpi_in_place,n_contact_z,1,mpi_real8,mpi_max,comm_cart,error)
    contact_angle = acos(n_wall_y*n_contact_y + n_wall_z*n_contact_z)
    if (myid .eq. 0) then
       write(6,'(A,1ES16.8)') 'Contact angle error = ',contact_angle-pi/4.*3
    endif
    
    
    return
  end subroutine ls_bound_error


end module mod_misc
