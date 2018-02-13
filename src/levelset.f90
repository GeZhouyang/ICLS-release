module mod_levelset

  use mod_common
  use mod_common_mpi
  use mod_param
  use mod_param_ls
  use mod_bound
  use mod_output
  use mod_init
  use mod_interface
  use mod_hyperPDE
  use mod_misc
  use mod_contact_line

  implicit none

  private
  public ls_main

contains


  !~~~~~~~~~~~~~~~~~~~~~~~! 
  !                       !
  !      Main script      !
  !                       !

  subroutine ls_main(istep)

    use mod_rib
    use mod_contact_line_common

    integer, intent(in) :: istep

    integer :: l
    real, dimension(num_cp) :: app_cont_angle

    !---- Riblet Level Set ----!
    
    if (Riblet) then
       call ls_main_rib(istep)
       return
    endif

    !---- Multiple Level Set ----!
    
    do l = 1,lmax

       if (Contact_line) call contact_line_main(istep, app_cont_angle) !## works only for one level set ##!

       !-- Narrow-band Advection (#0 #1)
       
       call tube_construction(l)
       call ls_advection(l)     
       
       !-- Mass Correction (#2)
       
       if (mass_correction) call ls_mass_correction(istep,l)

       !-- Re-initialization (#3)
       
       if (re_initialization) call ls_reinit_main(istep,l)

       !-- Contact angle
       
       if (Contact_line) call enforce_contact_angle(l,app_cont_angle)

    enddo

    !---- Update Viscosity, Density, and Volume ----! 

    call update_materials
    

    return
  end subroutine ls_main



  

 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !     
  ! #0

  subroutine tube_construction(l)

    ! Construct the narrow band tube of mls #l

    integer, intent(in) :: l

    integer :: i,j,k, cnt
    integer, dimension(imax,jmax,kmax) :: mask

    real :: xsum,ysum,zsum
    real, dimension(-2:2) :: xscan,yscan,zscan

    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phi

    !-- Previous level set
    
    phi = lvset(:,:,:,l)

    !-- Constant if outside of tube

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

             if (abs(phi(i,j,k)) .lt. NarrowBand_2) then

                mask(i,j,k) = 1
                cnt = cnt + 1

             else

                xscan(:) = abs(phi(i-2:i+2,j,k)) - NarrowBand_2
                xsum = sum(xscan)

                yscan(:) = abs(phi(i,j-2:j+2,k)) - NarrowBand_2
                ysum = sum(yscan)

                zscan(:) = abs(phi(i,j,k-2:k+2)) - NarrowBand_2
                zsum = sum(zscan)

                if (xsum .eq. 0. .and. ysum .eq. 0. .and. zsum .eq. 0.) then
                   mask(i,j,k) = 0
                else
                   mask(i,j,k) = 1
                   cnt = cnt + 1
                endif

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

    call boundsls(phi)  ! can comment, since the first triple loop is -2~+2 ?

    !-- Narrow banded level set
    
    lvset(:,:,:,l) = phi

    return
  end subroutine tube_construction




  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~
  ! 
  ! #1

  subroutine ls_advection(l)

    use mod_rib

    integer, intent(in) :: l

    real, dimension(1:imax,1:jmax,1:kmax) :: dfdt
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phi_0,phi_1,phi_2,phi_3

    !-- Unadvected level set

    phi_0 = lvset(:,:,:,l)

    !-- Advection schemes
    
    select case (ls_adv_space)
    case('HOUC5')

       select case (ls_adv_time)
       case('RK2')
          call HOUC5_nb(uccnew,vccnew,wccnew, phi_0, dfdt)
          call RK_nb(phi_1, phi_0, 1./2.,dfdt)
          call RK_nb(phi_2, phi_0, 1.,   dfdt)

          call HOUC5_nb(uccnew,vccnew,wccnew, phi_2, dfdt)
          call RK_nb(phi_3, phi_1, 1./2.,dfdt)

       case('RK3')
          call HOUC5_nb(uccnew,vccnew,wccnew, phi_0, dfdt)
          call RK_nb(phi_1, phi_0, 1.,dfdt)

          call HOUC5_nb(uccnew,vccnew,wccnew, phi_1, dfdt)
          call RK_nb(phi_2, 3./4.*phi_0 + 1./4.*phi_1, 1./4.,dfdt)

          call HOUC5_nb(uccnew,vccnew,wccnew, phi_2, dfdt)
          call RK_nb(phi_3, 1./3.*phi_0 + 2./3.*phi_2, 2./3.,dfdt)
       end select

    case('HOUC5rib')

       call HOUC5_nb_rib(uccnew,vccnew,wccnew, phi_0, dfdt)
       call RK_nb_rib(phi_1, phi_0, 1.,dfdt)

       call HOUC5_nb_rib(uccnew,vccnew,wccnew, phi_1, dfdt)
       call RK_nb_rib(phi_2, 3./4.*phi_0 + 1./4.*phi_1, 1./4.,dfdt)

       call HOUC5_nb_rib(uccnew,vccnew,wccnew, phi_2, dfdt)
       call RK_nb_rib(phi_3, 1./3.*phi_0 + 2./3.*phi_2, 2./3.,dfdt)
       
    case('WENO5')

       select case (ls_adv_time)
       case('RK2')
          call WENO5_nb(uccnew,vccnew,wccnew, phi_0, dfdt)
          call RK_nb(phi_1, phi_0, 1./2.,dfdt)
          call RK_nb(phi_2, phi_0, 1.,   dfdt)

          call WENO5_nb(uccnew,vccnew,wccnew, phi_2, dfdt)
          call RK_nb(phi_3, phi_1, 1./2.,dfdt)

       case('RK3')
          call WENO5_nb(uccnew,vccnew,wccnew, phi_0, dfdt)
          call RK_nb(phi_1, phi_0, 1.,dfdt)

          call WENO5_nb(uccnew,vccnew,wccnew, phi_1, dfdt)
          call RK_nb(phi_2, 3./4.*phi_0 + 1./4.*phi_1, 1./4.,dfdt)

          call WENO5_nb(uccnew,vccnew,wccnew, phi_2, dfdt)
          call RK_nb(phi_3, 1./3.*phi_0 + 2./3.*phi_2, 2./3.,dfdt)
       end select

    case('mWENO5c')

       select case (ls_adv_time)
       case('RK2')
          call mWENO5c_nb(unew,vnew,wnew, phi_0, dfdt)
          call RK_nb(phi_1, phi_0, 1./2.,dfdt)
          call RK_nb(phi_2, phi_0, 1.,dfdt)

          call mWENO5c_nb(unew,vnew,wnew, phi_2, dfdt)
          call RK_nb(phi_3, phi_1, 1./2.,dfdt)

       case('RK3')
          call mWENO5c_nb(unew,vnew,wnew, phi_0, dfdt)
          call RK_nb(phi_1, phi_0, 1.,dfdt)

          call mWENO5c_nb(unew,vnew,wnew, phi_1, dfdt)
          call RK_nb(phi_2, 3./4.*phi_0 + 1./4.*phi_1, 1./4.,dfdt)

          call mWENO5c_nb(unew,vnew,wnew, phi_2, dfdt)
          call RK_nb(phi_3, 1./3.*phi_0 + 2./3.*phi_2, 2./3.,dfdt)
       end select
       
    case('mWENO5nc')

       select case (ls_adv_time)
       case('RK2')
          call mWENO5nc_nb(uccnew,vccnew,wccnew, phi_0, dfdt)
          call RK_nb(phi_1, phi_0, 1./2.,dfdt)
          call RK_nb(phi_2, phi_0, 1.,dfdt)

          call mWENO5nc_nb(uccnew,vccnew,wccnew, phi_2, dfdt)
          call RK_nb(phi_3, phi_1, 1./2.,dfdt)

       case('RK3')
          call mWENO5nc_nb(uccnew,vccnew,wccnew, phi_0, dfdt)
          call RK_nb(phi_1, phi_0, 1.,dfdt)

          call mWENO5nc_nb(uccnew,vccnew,wccnew, phi_1, dfdt)
          call RK_nb(phi_2, 3./4.*phi_0 + 1./4.*phi_1, 1./4.,dfdt)

          call mWENO5nc_nb(uccnew,vccnew,wccnew, phi_2, dfdt)
          call RK_nb(phi_3, 1./3.*phi_0 + 2./3.*phi_2, 2./3.,dfdt)
       end select

    case('semiLag') ! Not accurate. Perhaps convenient with AMR. NOT recommended.

       call semiLagrangian(l,phi_0, phi_3) ! build-in RK2
       call boundsls(phi_3)

    end select

    !-- Advected level set

    lvset(:,:,:,l) = phi_3
    
    return
  end subroutine ls_advection
    

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! 
  ! #1.1
  
  subroutine RK_nb(phi_out, phi_in, factor,dfdt)

    real, intent(in) :: factor
    real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: dfdt
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in) :: phi_in
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(out) :: phi_out

    integer :: i,j,k, cnt

    phi_out = phi_in

    do cnt = 1,nb_cnt

       i = nb_i(cnt)
       j = nb_j(cnt)
       k = nb_k(cnt)

       phi_out(i,j,k) = phi_in(i,j,k) + factor*dt*dfdt(i,j,k)

    enddo

    call boundsls(phi_out)


    return
  end subroutine RK_nb





  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! #2

  subroutine ls_mass_correction(istep,l)

    integer, intent(in) :: istep,l

    real :: massloss
    real, parameter :: eps = 1e-5  ! mass correction threshold
    
    maxchk_glb = 0.
    
    if (mod(istep,inflate_step) .eq. 0) then

       massloss = (volume(l) - init_volume(l))/init_volume(l)
       
       if (abs(massloss) .gt. eps) then

          if (myid .eq. 0) write(*,'(A,I3)') '!-- Inflate interface',l
          call inflation(l)

       endif

    endif

    return
  end subroutine ls_mass_correction


  !~~~~~~~~~~~~~~~~~~~~~~
  !                   
  ! #2.1
  
  subroutine inflation(l)

    integer, intent(in) :: l

    integer :: i,j,k
    real :: phi_s, dd, temp, p_l, kappamax

    real, parameter :: s1 = 1.5*dz

    real, dimension(lmax) :: area, loss_rate
    real, dimension(1:imax,1:jmax,1:kmax) :: delta, u_lc,v_lc,w_lc, dfdt
    real, dimension(0:i1,0:j1,0:k1) :: dheav, dummy, kappa
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: phi,dphi, phi_1,phi_2,phi_3

    dummy = 0.

    area(l) = 0.
    dheav = 1.
    phi = lvset(:,:,:,l)

    call get_curvature(phi,kappa)
    call boundc(kappa)

    select case (ispeed_dependence)

    case('unfm')

       delta = Dirac('trig', l, s1)
       dheav = Heaviside('0','trig',l, s1)

       !-- Surface area --!

       temp = sum(delta*normal%m) *dx*dy*dz  !# i included abs grad phi here
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
                u_lc(i,j,k) = (u_l(i-1,j,k)+u_l(i,j,k))/2.
                v_lc(i,j,k) = (v_l(i,j-1,k)+v_l(i,j,k))/2.
                w_lc(i,j,k) = (w_l(i,j,k-1)+w_l(i,j,k))/2.
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

       call HOUC5_nb(u_lc,v_lc,w_lc, phi, dfdt)
       call RK_nb(phi_1, phi, 1.,dfdt)

       call HOUC5_nb(u_lc,v_lc,w_lc, phi_1, dfdt)
       call RK_nb(phi_2, 3./4.*phi + 1./4.*phi_1, 1./4.,dfdt)

       call HOUC5_nb(u_lc,v_lc,w_lc, phi_2, dfdt)
       call RK_nb(phi_3, 1./3.*phi + 2./3.*phi_2, 2./3.,dfdt)

    case('WENO5')

       call WENO5_nb(u_lc,v_lc,w_lc, phi, dfdt)
       call RK_nb(phi_1, phi, 1.,dfdt)

       call WENO5_nb(u_lc,v_lc,w_lc, phi_1, dfdt)
       call RK_nb(phi_2, 3./4.*phi + 1./4.*phi_1, 1./4.,dfdt)

       call WENO5_nb(u_lc,v_lc,w_lc, phi_2, dfdt)
       call RK_nb(phi_3, 1./3.*phi + 2./3.*phi_2, 2./3.,dfdt)
    end select

    lvset(:,:,:,l) = phi_3

    return
  end subroutine inflation





  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! #3

  subroutine ls_reinit_main(istep,l)

    integer, intent(in) :: istep,l

    integer :: iter
    real :: residual

    if (mod(istep,reinit_step) .eq. 0) then

       call ls_reinit(l,reinit_iter_max, iter,residual)
       call boundsls(lvset(:,:,:,l))
       
       if (myid .eq. 0) then
          write(*,'(A,I3,A,I5,A,ES10.4)') '!-- Reinit LS',l, &
               ' for',iter, ' iterations w/ max residual = ',residual
       endif

    endif
          
    return
  end subroutine ls_reinit_main

  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !                                  
  ! #3.1
  
  subroutine ls_reinit(marker,maxiter,iter,residual)
                                                           
    ! Get signed distance function (phi_0 --> phi_new) by evolving
    !                                                          
    !   dø                                                     
    ! ------ = sign(ø0)*( 1 - |grad ø| )                       
    !  dtau                                                    
    !                                                          
    ! iteratively using combinations of EU1/RK2/RK3 and RUSM1/RUSM2/RUSM4/WENO5.
    ! (Notes: RUSM4 is not working. WENO5 is recommeded. See Ge et al, JCP2018.)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    

    integer, intent(in) :: marker,maxiter
    integer, intent(out) :: iter
    real, intent(out) :: residual

    integer :: i,j,k, iterations
    integer, dimension(2) :: ind_loc
    integer, dimension(1:imax,1:jmax,1:kmax) :: Indx_x,Indx_y,Indx_z

    logical, parameter :: validation = .false.

    real :: dtau,smooth_width, grad_min,grad_max,c1

    real, dimension(1:imax,1:jmax,1:kmax) :: sgn,region,Dist
    real, dimension(1:imax,1:jmax,1:kmax) :: flag_x,flag_y,flag_z
    real, dimension(1:imax,1:jmax,1:kmax) :: dist_x,dist_y,dist_z
    real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: RHS, phi_0, phi_n1, phi_n2, phi_new


    phi_0(:,:,:) = lvset(:,:,:,marker)  ! time level n
    phi_new = phi_0                     ! time level n+1

    select case (ls_reinit_space)
    case('RUSM1','RUSM2')

       call RussoSmereka_init(phi_0, sgn,region,Dist)

    case('RUSM4')

       call RussoSmereka_init_4(phi_0, sgn,region,Dist, flag_x,flag_y,flag_z,dist_x,dist_y,dist_z)

    case('WENO5','WENO5z')
       
       !---- Smoothed Sign Function ----!
       
       smooth_width = 1.5*dz

       do k=1,kmax
          do j=1,jmax
             do i=1,imax   

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

             enddo
          enddo
       enddo
       
    end select

    !---- Reinitialization Loop ----!

    dtau = 0.5*dz

    do iterations = 1,maxiter      

       select case (ls_reinit_time)
       case('EU1')  ! forward Euler

          select case (ls_reinit_space)
          case('RUSM1')

             call RussoSmereka_iterate_1(sgn,region,Dist,phi_0, RHS)
             phi_new = phi_0 +dtau*RHS
             call boundsls(phi_new)

          case('RUSM2')
          
             call RussoSmereka_iterate_2(sgn,region,Dist,phi_0, RHS)
             phi_new = phi_0 +dtau*RHS
             call boundsls(phi_new)

          case('RUSM4')
          
             call RussoSmereka_iterate_4(sgn,region,Dist, flag_x,flag_y,flag_z,dist_x,dist_y,dist_z,phi_0, RHS)
             phi_new = phi_0 +dtau*RHS
             call boundsls(phi_new)
          end select
          
       case('RK2')  ! Mid-Point RK2

          select case (ls_reinit_space)
          case('RUSM1')
          
             call RussoSmereka_iterate_1(sgn,region,Dist,phi_0, RHS)
             phi_n1 = phi_0 + dtau*RHS
             call boundsls(phi_n1)
          
             call RussoSmereka_iterate_1(sgn,region,Dist,phi_n1, RHS)
             phi_n2 = phi_n1 + dtau*RHS
             call boundsls(phi_n2)
             
             phi_new = (phi_n1 + phi_n2)/2.
          
          case('RUSM2')
          
             call RussoSmereka_iterate_2(sgn,region,Dist,phi_0, RHS)
             phi_n1 = phi_0 + 0.5*dtau*RHS
             phi_n2 = phi_0 + dtau*RHS
             call boundsls(phi_n1)
             call boundsls(phi_n2)
          
             call RussoSmereka_iterate_2(sgn,region,Dist,phi_n2, RHS)
             phi_new = phi_n1 + 0.5*dtau*RHS
             call boundsls(phi_new)

          case('WENO5')

             call WENO5_reinit(sgn,phi_0,RHS)
             phi_n1 = phi_0 + 0.5*dtau*RHS
             phi_n2 = phi_0 + dtau*RHS
             call boundsls(phi_n1)
             call boundsls(phi_n2)

             call WENO5_reinit(sgn,phi_n2,RHS)
             phi_new = phi_n1 + 0.5*dtau*RHS
             call boundsls(phi_new)

          end select

       case('RK3')  ! SSP RK3

          select case (ls_reinit_space)             
          case('RUSM1')
          
             call RussoSmereka_iterate_1(sgn,region,Dist,phi_0, RHS)
             phi_n1 = phi_0 + dtau*RHS
             call boundsls(phi_n1)
          
             call RussoSmereka_iterate_1(sgn,region,Dist,phi_n1, RHS)
             phi_n2 = 3./4.*phi_0 + 1./4.*(phi_n1 + dtau*RHS)
             call boundsls(phi_n2)
          
             call RussoSmereka_iterate_1(sgn,region,Dist,phi_n2, RHS)
             phi_new = 1./3.*phi_0 + 2./3.*(phi_n2 + dtau*RHS)
             call boundsls(phi_new)

          case('RUSM4')
          
             call RussoSmereka_iterate_4(sgn,region,Dist, flag_x,flag_y,flag_z,dist_x,dist_y,dist_z,phi_0, RHS)
             phi_n1 = phi_0 + dtau*RHS
             call boundsls(phi_n1)

             call RussoSmereka_iterate_4(sgn,region,Dist,flag_x,flag_y,flag_z,dist_x,dist_y,dist_z,phi_n1, RHS)
             phi_n2 = 3./4.*phi_0 + 1./4.*(phi_n1 + dtau*RHS)
             call boundsls(phi_n2)

             call RussoSmereka_iterate_4(sgn,region,Dist,flag_x,flag_y,flag_z,dist_x,dist_y,dist_z,phi_n2, RHS)
             phi_new = 1./3.*phi_0 + 2./3.*(phi_n2 + dtau*RHS)
             call boundsls(phi_new)

          case('WENO5')

             call WENO5_reinit(sgn,phi_0,RHS)
             phi_n1 = phi_0 + dtau*RHS
             call boundsls(phi_n1)

             call WENO5_reinit(sgn,phi_n1,RHS)
             phi_n2 = 3./4.*phi_0 + 1./4.*(phi_n1 + dtau*RHS)
             call boundsls(phi_n2)

             call WENO5_reinit(sgn,phi_n2,RHS)
             phi_new = 1./3.*phi_0 + 2./3.*(phi_n2 + dtau*RHS)
             call boundsls(phi_new)

          case('WENO5z')

             call WENO5z_reinit(sgn,phi_0,RHS)
             phi_n1 = phi_0 + dtau*RHS
             call boundsls(phi_n1)

             call WENO5z_reinit(sgn,phi_n1,RHS)
             phi_n2 = 3./4.*phi_0 + 1./4.*(phi_n1 + dtau*RHS)
             call boundsls(phi_n2)

             call WENO5z_reinit(sgn,phi_n2,RHS)
             phi_new = 1./3.*phi_0 + 2./3.*(phi_n2 + dtau*RHS)
             call boundsls(phi_new)
          end select

       end select

       !-- check residual

       residual = maxval(abs(phi_new-phi_0))
       call mpi_allreduce(mpi_in_place,residual,1,mpi_real8,mpi_max,comm_cart,error)
       if (residual .lt. dz**5) then
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
  end subroutine ls_reinit



  
end module mod_levelset
