module mod_gfm

  use mod_param

  implicit none

  private
  public gfm4p, gfm4diff

contains


  !------------------------------------------------------
  !-- The ghost fluid method (GFM) for the computation of
  !-- surface tenseion across a sharp interface.
  !
  subroutine gfm4p(p_jump)

    use mod_common
    use mod_interface
    use mod_bound
    use mod_common_mpi
    
    integer :: i,j,k,l, nx,ny,nz,flag, cnt

    real :: nondi, factor,idvndn, t1,tmax
    real :: icurv, jump_x, jump_y, jump_z
    real :: cmax,cmin, cmax_all,cmin_all !#

    real, dimension(-1:1) :: ph,cv,ps, jp_x, jp_y, jp_z

    real, dimension(0:i1,0:j1,0:k1) :: curvature, p_exjump, cux,cuy,cuz, re_curv !#
    real, dimension(0:i1,0:j1,0:k1,1:lmax) :: px_sep,py_sep,pz_sep, pj_sep

    real, dimension(0:i1,0:j1,0:k1), intent(out) :: p_jump


    nondi = 1./We  ! non-dimensional surface tension coef

    select case (NS_diff_discretization)
    case('del')
       factor = 0.

    case('GFM')
       factor = 2.*visc*(miu2-miu1)
       if (myid .eq. 0) then
          write(*,*) 'This method seems unrobust. Use at your own risk.'
          write(*,*) 'Program aborted. (Comment these if you insist.)'
       endif
       call mpi_finalize(error)
       stop
    end select


    select case (surface_tension_method)
    case('CSF') ! add nothing
       p_jump = 0.

    case('GFM')

       call ipot  ! get minimal distances between interfaces

       !-- Update jump at time level (n-1) --!

       p_xold = p_x
       p_yold = p_y
       p_zold = p_z

       !-- Compute jump at time level (n) --!

       px_sep = 0.
       py_sep = 0.
       pz_sep = 0.
       pj_sep = 0.
       p_jump = 0.
       p_exjump = 0.
       p_x = 0. 
       p_y = 0.
       p_z = 0.

       do l=1,lmax

          !-- Obtain cell-center curvature of each drop --!

          call get_curvature(lvset(:,:,:,l),curv_cmn)
          call boundc(curv_cmn)
          !#
          if (.not. TwoD) then
             call re_curvature(lvset(:,:,:,l),curv_cmn,re_curv)
             curv_cmn = re_curv
             call boundc(curv_cmn)
          endif
          !#

          do k=1,kmax
             do j=1,jmax
                do i=1,imax

                   if ( abs(lvset(i,j,k,l)) .le. 2.*dz) then ! jump only occurs across cut cells

                      jump_x = 0.
                      jump_y = 0.
                      jump_z = 0.
                      jp_x = 0.
                      jp_y = 0.
                      jp_z = 0.

                      !-- x component --! 
                      
                      do nx = -1,1,2

                         ph(-1:1) =    lvset(i-1:i+1,j,k,l) 
                         cv(-1:1) = curv_cmn(i-1:i+1,j,k)

                         if ( ph(0)*ph(nx) .lt. 0.) then

                            call icurvature(nx, ph,cv, icurv)

                            !-- jump due to surface tension
                            
                            jp_x(nx) = nondi*(-icurv)

                            !-- jump due to discontinuous viscosity

                            !#idvndn = ( dvndn(i,j,k)*abs(ph(nx)) + dvndn(i+nx,j,k)*abs(ph(0)) ) &
                            !#     /( abs(ph(0)) + abs(ph(nx)) )
                            !#
                            !#jump_x = jump_x + factor*idvndn

                            if (nx .eq. 1) px_sep(i,j,k,l) = jp_x(nx)

                         endif
                      enddo
                      
                      jump_x = sum(jp_x)
                      
                      !-- y component --! 

                      do ny = -1,1,2

                         ph(-1:1) =    lvset(i,j-1:j+1,k,l) 
                         cv(-1:1) = curv_cmn(i,j-1:j+1,k)

                         if ( ph(0)*ph(ny) .lt. 0.) then

                            call icurvature(ny, ph,cv, icurv)

                            !-- jump due to surface tension

                            jp_y(ny) = nondi*(-icurv)

                            !-- jump due to discontinuous viscosity

                            !#idvndn = ( dvndn(i,j,k)*abs(ph(ny)) + dvndn(i,j+ny,k)*abs(ph(0)) ) &
                            !#     /( abs(ph(0)) + abs(ph(ny)) )
                            !#
                            !#jump_y = jump_y + factor*idvndn

                            if (ny .eq. 1) py_sep(i,j,k,l) = jp_y(ny)       

                         endif
                      enddo
                      
                      jump_y = sum(jp_y)

                      !-- z component --! 

                      do nz = -1,1,2

                         ph(-1:1) =    lvset(i,j,k-1:k+1,l) 
                         cv(-1:1) = curv_cmn(i,j,k-1:k+1)

                         if ( ph(0)*ph(nz) .lt. 0.) then

                            call icurvature(nz, ph,cv, icurv)

                            !-- jump due to surface tension

                            jp_z(nz) = nondi*(-icurv)

                            !-- jump due to discontinuous viscosity

                            !#idvndn = ( dvndn(i,j,k)*abs(ph(nz)) + dvndn(i,j,k+nz)*abs(ph(0)) ) &
                            !#     /( abs(ph(0)) + abs(ph(nz)) )
                            !#
                            !#jump_z = jump_z + factor*idvndn                          

                            if (nz .eq. 1) pz_sep(i,j,k,l) = jp_z(nz)

                         endif
                      enddo
                      
                      jump_z = sum(jp_z)
                      

                      !-- "Laplacian" of pressure jump (separately) --! 

                      if (lvset(i,j,k,l) .gt. 0.) then  ! jump higher

                         pj_sep(i,j,k,l) =   (jump_x/dx**2 + jump_y/dy**2 + jump_z/dz**2)

                      elseif (lvset(i,j,k,l) .lt. 0.) then  ! drop lower

                         pj_sep(i,j,k,l) = - (jump_x/dx**2 + jump_y/dy**2 + jump_z/dz**2)
                         px_sep(i,j,k,l) = - px_sep(i,j,k,l)
                         py_sep(i,j,k,l) = - py_sep(i,j,k,l)
                         pz_sep(i,j,k,l) = - pz_sep(i,j,k,l)

                      else
                         pj_sep(i,j,k,l) = 0.
                         px_sep(i,j,k,l) = 0.
                         py_sep(i,j,k,l) = 0.
                         pz_sep(i,j,k,l) = 0.
                         write(*,*) 'Happens to be at the interface (highly unlikely)!'
                      endif

                      !-- Account for close neighbors --! 

                      if (l .eq. 1) then  ! initialize

                         p_jump(i,j,k) = pj_sep(i,j,k,l)
                         p_x(i,j,k) = px_sep(i,j,k,l)
                         p_y(i,j,k) = py_sep(i,j,k,l)
                         p_z(i,j,k) = pz_sep(i,j,k,l)

                      else ! superimpose each separate jump

                         p_jump(i,j,k) = p_jump(i,j,k) + pj_sep(i,j,k,l)
                         p_x(i,j,k) = p_x(i,j,k) + px_sep(i,j,k,l)
                         p_y(i,j,k) = p_y(i,j,k) + py_sep(i,j,k,l)
                         p_z(i,j,k) = p_z(i,j,k) + pz_sep(i,j,k,l)

                      endif

                   endif ! jump shell
                enddo ! i
             enddo ! j
          enddo ! k
       enddo ! l

       !!###########
       !
       !      cmax = maxval(cuz(:,:,:))
       !      call mpi_allreduce(cmax,cmax_all,1,mpi_real8,mpi_max,comm_cart,error)
       !      cmin = minval(cuz(:,:,:))
       !      call mpi_allreduce(cmin,cmin_all,1,mpi_real8,mpi_min,comm_cart,error)
       !
       !
       !      if (myid .eq. 2) then
       !        open(60,file=datadir//'curv.txt',position='append')
       !        write(60,'(2ES12.4)') cmax_all,cmin_all
       !        close(60)
       !      endif
       !
       !!###########

       !-- Extral jump due to depletion --!

       if (attractive) then 
          call exjump(p_exjump)
          p_jump = p_jump + p_exjump
       endif
    end select


    return
  end subroutine gfm4p



  !-----------------------------------------------
  !-- A hydrodynamic model for the depletion force
  !-- using MLS and GFM
  !
  subroutine exjump(p_jump)

    use mod_common
    use mod_common_mpi
    use mod_bound

    integer :: i,j,k,l, m,n, nx,ny,nz
    real :: phi, phx,phy,phz
    real :: jump_x, jump_y, jump_z

    real :: del_r, dvol,vol_ex, pot_tot, cp, c0,c1, osp
    
    real, dimension(-1:1) :: fo, pf, jp_x, jp_y, jp_z
    real, dimension(0:i1,0:j1,0:k1) :: phim,phin, flago, pfactor, px_osp,py_osp,pz_osp

    !-- Model Parameters --!

    real, parameter :: rsm = 2.0*dz    ! radius of surfactant micelles
    real, parameter :: rs1 = 0.1*rsm   ! reduce osp when less than this distance
    real, parameter :: rs2 = 0.0*rsm   ! negative osp when less than this distance
    real, parameter :: osp0 = -40.     ! osmotic pressure due to surfactant micelles

    real, dimension(0:i1,0:j1,0:k1), intent(out) :: p_jump

    px_osp = 0.
    py_osp = 0.
    pz_osp = 0.
    p_jump = 0.
    
    dvol = dx*dy*dz

    do m = 1,lmax-1
       phim = lvset(0:i1,0:j1,0:k1,m)  ! droplet m

       do n = m+1,lmax
          phin = lvset(0:i1,0:j1,0:k1,n)  ! droplet n

          !-- Reduce magnitude if too close --!
          
          if (dmin(m,n) .gt. rs1) then
             osp = osp0
          else
             del_r = (rs1-dmin(m,n))/(rs1-rs2)
             osp = osp0*(1.-del_r)
          endif

          !-- Identify overlap of surfactant shell --!

          where (phim .le. rsm .and. phin .le. rsm)
             flago = 1.  ! inside overlap
          elsewhere
             flago = 0.  ! outside overlap
          end where
          
          call boundc(flago)

          !-- Compute overlap volume --!

          vol_ex = sum( flago(1:imax,1:jmax,1:kmax) )*dvol  ! a rough count
          call mpi_allreduce(mpi_in_place,vol_ex,1,mpi_real8,mpi_sum,comm_cart,error)

          pot_tot = osp*vol_ex  ! total potential energy

          !-- Compute the spatially varying osmotic pressure --!

          pfactor = 0.
          
          do k=1,kmax
             do j=1,jmax
                do i=1,imax
                   if (flago(i,j,k) .eq. 1.) then
                      
                      pfactor(i,j,k) = (phim(i,j,k) + phin(i,j,k))/(2.*rsm) -1.

                   endif
                enddo
             enddo
          enddo

          call boundc(pfactor)

          c0 = sum( pfactor(1:imax,1:jmax,1:kmax) )*dvol
          call mpi_allreduce(mpi_in_place,c0,1,mpi_real8,mpi_sum,comm_cart,error)

          cp = pot_tot/c0  ! a constant prefactor

          !-- Impose the osmotic pressure jump --!

          px_osp = 0.
          py_osp = 0.
          pz_osp = 0.
                      
          do k=1,kmax
             do j=1,jmax
                do i=1,imax 

                   if (abs(slset(i,j,k)) .le. 2.*rsm) then ! jump within this shell

                      jump_x = 0.
                      jump_y = 0.
                      jump_z = 0.
                      jp_x = 0.
                      jp_y = 0.
                      jp_z = 0.

                      !-- x component
                      
                      do nx = -1,1,2
                        
                         fo(-1:1) =   flago(i-1:i+1,j,k)
                         pf(-1:1) = pfactor(i-1:i+1,j,k)
                         
                         if ( fo(0) .ne. fo(nx) ) then

                            if ( fo(0) .eq. 1.) then  ! (i,j,k) is inside overlap
                               c1 = pf(0)
                            else
                               c1 = pf(nx)
                            endif
                            
                            jp_x(nx) = cp*c1
                            if (nx .eq. 1) px_osp(i,j,k) = jp_x(nx)
                            
                         endif
                      enddo
                      
                      jump_x = sum(jp_x)

                      !-- y component
                      
                      do ny = -1,1,2
                         
                         fo(-1:1) =   flago(i,j-1:j+1,k)
                         pf(-1:1) = pfactor(i,j-1:j+1,k)
                         
                         if ( fo(0) .ne. fo(ny) ) then

                            if ( fo(0) .eq. 1.) then  ! (i,j,k) is inside overlap
                               c1 = pf(0)
                            else
                               c1 = pf(ny)
                            endif
                            
                            jp_y(ny) = cp*c1
                            if (ny .eq. 1) py_osp(i,j,k) = jp_y(ny)
                            
                         endif
                      enddo
                      
                      jump_y = sum(jp_y)

                      !-- z component
                                            
                      do nz = -1,1,2
                         
                         fo(-1:1) =   flago(i,j,k-1:k+1)
                         pf(-1:1) = pfactor(i,j,k-1:k+1)
                         
                         if ( fo(0) .ne. fo(nz) ) then

                            if ( fo(0) .eq. 1.) then  ! (i,j,k) is inside overlap
                               c1 = pf(0)
                            else
                               c1 = pf(nz)
                            endif
                            
                            jp_z(nz) = cp*c1
                            if (nz .eq. 1) pz_osp(i,j,k) = jp_z(nz)
                            
                         endif
                      enddo
                      
                      jump_z = sum(jp_z)
                      
                      !-- "Laplacian" of osmotic pressure jump (superimposed) --!

                      if (flago(i,j,k) .gt. 0.) then  ! jump from overlap to outside

                         p_jump(i,j,k) = p_jump(i,j,k) + (jump_x/dx**2 + jump_y/dy**2 + jump_z/dz**2)

                      else  ! jump from outside to overlap

                         p_jump(i,j,k) = p_jump(i,j,k) - (jump_x/dx**2 + jump_y/dy**2 + jump_z/dz**2)
                         px_osp(i,j,k) = - px_osp(i,j,k)
                         py_osp(i,j,k) = - py_osp(i,j,k)
                         pz_osp(i,j,k) = - pz_osp(i,j,k)

                      endif                      
                      
                      !-- Correction to the pressure gradient --!
                      
                      p_x(i,j,k) = p_x(i,j,k) + px_osp(i,j,k)
                      p_y(i,j,k) = p_y(i,j,k) + py_osp(i,j,k)
                      p_z(i,j,k) = p_z(i,j,k) + pz_osp(i,j,k)
                      
                   endif ! jump shell
                enddo ! i
             enddo ! j
          enddo ! k
       enddo  ! n
    enddo ! m


    return
  end subroutine exjump

  
  
  
  !---------------------------------------------------
  !-- The lengthy code below is for the computation of
  !-- a discontinuous viscosity using GFM
  !
  subroutine ls_c2f(dir,phi)

    ! Average level set values from cell centers to faces

    use mod_common
    use mod_bound

    character(len=6), intent(in) :: dir   
    real, dimension(0:i1,0:j1,0:k1), intent(out) :: phi

    integer :: i,j,k, cnt

    phi = lvset(0:i1,0:j1,0:k1, 1)  !## now only consider the first level set ##temp

    select case (dir)

    case('u-grid')

       do cnt = 1,nb_cnt
          i = nb_i(cnt)
          j = nb_j(cnt)
          k = nb_k(cnt)
          phi(i,j,k) = ( lvset(i,j,k,1) +lvset(i+1,j,k,1) )/2.  !##temp
       enddo

    case('v-grid')

       do cnt = 1,nb_cnt
          i = nb_i(cnt)
          j = nb_j(cnt)
          k = nb_k(cnt)
          phi(i,j,k) = ( lvset(i,j,k,1) +lvset(i,j+1,k,1) )/2.  !##temp
       enddo

    case('w-grid')

       do cnt = 1,nb_cnt
          i = nb_i(cnt)
          j = nb_j(cnt)
          k = nb_k(cnt)
          phi(i,j,k) = ( lvset(i,j,k,1) +lvset(i,j,k+1,1) )/2.  !##temp
       enddo

    end select

    call boundc(phi) !# wall halos untouched

    return
  end subroutine ls_c2f


  !------------------------
  !
  !
  subroutine gradjump(jump)

    ! Compute the jumpy velocity-gradient matrix

    use mod_common
    use mod_bound
    use mod_interface
    use mod_newtypes

    type (real_3by3_matrix), dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(out) :: jump !temp bigger than needed

    integer :: i,j,k, cnt, m

    real :: mag

    real, dimension(3) :: dummy, tan_vec1,tan_vec2
    real, dimension(1,3) :: nor_vec, temp
    real, dimension(3,3) :: grad_vel, nor_sq,tan_sq, subt, term1,term2,term3
    real, dimension(0:i1,0:j1,0:k1) :: u,v,w


    u = 0.
    v = 0.
    w = 0.

    subt = 0.

    dvndn = 0.

    ! initialize gradient jump matrix

    forall(i=-2:i1+2,j=-2:j1+2,k=-2:k1+2)
       jump(i,j,k)%grad = subt
    end forall

    if (miu1 .eq. miu2) return

    ! obtain cell-center velocities

    do k = 1,kmax
       do j = 1,jmax
          do i = 1,imax   
             u(i,j,k) = (unew(i,j,k) + unew(i-1,j,k))/2.
             v(i,j,k) = (vnew(i,j,k) + vnew(i,j-1,k))/2.
             w(i,j,k) = (wnew(i,j,k) + wnew(i,j,k-1))/2.
          enddo
       enddo
    enddo

    call bounduvw(u,v,w)

    ! obtain normals (also curvature but not used now, can improve efficiency later)

    call get_curvature(lvset(:,:,:,1),curv_cmn)  !## now only the first level set ##temp
    call boundc(curv_cmn)

    ! construct matrices within a narrow band

    do cnt = 1,nb_cnt

       i = nb_i(cnt)
       j = nb_j(cnt)
       k = nb_k(cnt)

       ! velocity gradient matrix

       grad_vel(1,1) = (u(i+1,j,k) - u(i-1,j,k))/2./dx
       grad_vel(1,2) = (u(i,j+1,k) - u(i,j-1,k))/2./dy
       grad_vel(1,3) = (u(i,j,k+1) - u(i,j,k-1))/2./dz
       grad_vel(2,1) = (v(i+1,j,k) - v(i-1,j,k))/2./dx
       grad_vel(2,2) = (v(i,j+1,k) - v(i,j-1,k))/2./dy
       grad_vel(2,3) = (v(i,j,k+1) - v(i,j,k-1))/2./dz
       grad_vel(3,1) = (w(i+1,j,k) - w(i-1,j,k))/2./dx
       grad_vel(3,2) = (w(i,j+1,k) - w(i,j-1,k))/2./dy
       grad_vel(3,3) = (w(i,j,k+1) - w(i,j,k-1))/2./dz

       ! normal and tangential vectors

       nor_vec(1,1) = normal(i,j,k)%x
       nor_vec(1,2) = normal(i,j,k)%y
       nor_vec(1,3) = normal(i,j,k)%z

       dummy(:) = abs(nor_vec(1,:))
       m = minloc(dummy, dim=1)
       dummy(:) = 0.
       dummy(m) = 1.

       tan_vec1(1) =   nor_vec(1,2)*dummy(3) - nor_vec(1,3)*dummy(2)
       tan_vec1(2) = - nor_vec(1,1)*dummy(3) + nor_vec(1,3)*dummy(1)
       tan_vec1(3) =   nor_vec(1,1)*dummy(2) - nor_vec(1,2)*dummy(1)
       mag = dot_product(tan_vec1, tan_vec1)
       tan_vec1 = tan_vec1/sqrt(mag)

       tan_vec2(1) =   nor_vec(1,2)*tan_vec1(3) - nor_vec(1,3)*tan_vec1(2)
       tan_vec2(2) = - nor_vec(1,1)*tan_vec1(3) + nor_vec(1,3)*tan_vec1(1)
       tan_vec2(3) =   nor_vec(1,1)*tan_vec1(2) - nor_vec(1,2)*tan_vec1(1)

       ! corresponding matrices

       nor_sq = matmul(transpose(nor_vec),nor_vec)

       tan_sq(1,:) = 0.
       tan_sq(2,:) = tan_vec1(:)
       tan_sq(3,:) = tan_vec2(:)

       ! first-derivative jump matrix (w/o viscosity term)

       subt = matmul(grad_vel, transpose(tan_sq))
       term1 = matmul(subt, tan_sq)

       subt = matmul(nor_sq, grad_vel)
       term2 = matmul(subt, nor_sq)

       term3 = matmul(transpose(tan_sq), tan_sq)
       subt = matmul(term3, transpose(grad_vel))
       term3 = matmul(subt, nor_sq)

       jump(i,j,k)%grad = term1 + term2 - term3

       ! normal derivative of normal velocity (part of pressure jump)

       temp = matmul(nor_vec,grad_vel)
       dvndn(i,j,k) = dot_product( temp(1,:) ,nor_vec(1,:) )

    enddo

    call boundj(jump)
    call boundc(dvndn)


    return
  end subroutine gradjump


  !---------------------------------------------------------
  !
  !
  subroutine gfm4diff(RDu,RDv,RDw, u,v,w, rho_u,rho_v,rho_w)

    use mod_newtypes

    real, dimension(0:,0:,0:), intent(in) :: u,v,w
    real, dimension(1:imax,1:jmax,1:kmax), intent(in) :: rho_u, rho_v, rho_w
    real, dimension(1:imax,1:jmax,1:kmax), intent(out) :: RDu,RDv,RDw

    integer :: i,j,k, im,ip,jm,jp,km,kp, nb

    real :: s1, theta, miu_pls,miu_mns,miuh
    real :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
    real :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
    real :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real :: RD1,RD2,RD3,RD4,RD5
    real :: J_l,J_m,J_r, Ji

    real, dimension(-1:1) :: ph, jmean
    real, dimension(-1:1) :: dudx,dudy,dudz, dvdx,dvdy,dvdz, dwdx,dwdy,dwdz
    real, dimension(0:i1,0:j1,0:k1) :: phi

    type (real_3by3_matrix), dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: jump


    !== Obtain velocity-gradient-jump matrix at cell centers ==!

    call gradjump(jump)  ! visc jump is NOT multiplied


    !-- v-component --! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    call ls_c2f('v-grid',phi)  ! get level set on the v grid
    ph = 0.                    ! level set buffer
    s1 = 3.*dz                 ! tube half width

    dvdx = 0.
    dudy = 0.
    dvdy = 0.
    dvdz = 0.
    dwdy = 0.

    do k = 1,kmax
       do j = 1,jmax
          do i = 1,imax

             ip = i + 1
             jp = j + 1
             kp = k + 1
             im = i - 1
             jm = j - 1
             km = k - 1

             ! continuous first derivatives

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

             RD1 = (dvdxp - dvdxm)/dx     ! d/dx(dv/dx)
             RD3 = (dvdyp - dvdym)/dy     ! d/dy(dv/dy)
             RD4 = (dvdzp - dvdzm)/dz     ! d/dz(dv/dz)        

             ph(0) = phi(i,j,k)

             if (ph(0) .gt. s1) then       ! safely in fluid 1

                RDv(i,j,k) = miu1*(RD1+RD3+RD4)/rho1

             elseif (ph(0) .lt. -s1) then  ! safely in fluid 2

                RDv(i,j,k) = miu2*(RD1+RD3+RD4)/rho2

             else

                !- 2nd x-derivative -!

                if (TwoD) then

                   RD1 = 0.
                   RD2 = 0.

                else

                   ph(-1) = phi(i-1,j,k)
                   ph( 1) = phi(i+1,j,k)

                   if (minval(ph) .gt. 0.) then      ! in fluid 1

                      RD1 = miu1*RD1
                      RD2 = 0.

                   elseif (maxval(ph) .lt. 0.) then  ! in fluid 2

                      RD1 = miu2*RD1
                      RD2 = 0.

                   else

                      if (ph(0) .gt. 0.) then
                         miu_pls = miu1
                         miu_mns = miu2
                      else
                         miu_pls = miu2
                         miu_mns = miu1
                      endif

                      if ( ph(-1)*ph(0) .le. 0.) then
                         nb = -1
                         jmean(-1) = (jump(i-1,j,k)%grad(2,1) + jump(i-1,j+1,k)%grad(2,1))/2.!#
                      else
                         nb = 1
                         jmean(1) = (jump(i+1,j,k)%grad(2,1) + jump(i+1,j+1,k)%grad(2,1))/2.!#
                      endif

                      theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                      jmean(0) = (jump(i,j,k)%grad(2,1) + jump(i,j+1,k)%grad(2,1))/2.
                      Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                      Ji = Ji*(miu_pls-miu_mns)

                      miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                      dvdx(-1) = dvdxm
                      dvdx( 1) = dvdxp

                      dvdx( nb) = miuh*dvdx(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                      dvdx(-nb) = miu_pls*dvdx(-nb)

                      RD1 = (dvdx(1)-dvdx(-1))/dx  ! d/dx(miu*dv/dx)

                      ! Below L and R locates on the vertices

                      ph(-1) = (phi(i,j,k)+phi(i-1,j,k))/2.
                      ph( 1) = (phi(i,j,k)+phi(i+1,j,k))/2.

                      if (minval(ph) .gt. 0.) then  ! in fluid 1

                         RD2 = miu1*(dudyp-dudym)/dx

                      elseif ( maxval(ph) .lt. 0. ) then  ! in fluid 2

                         RD2 = miu2*(dudyp-dudym)/dx

                      else  ! between L and R

                         if ( ph(-1)*ph(0) .le. 0.) then
                            nb = -1
                            jmean(-1) = (jump(i-1,j,k)%grad(1,2) + jump(i-1,j+1,k)%grad(1,2) &
                                 +jump(i,j  ,k)%grad(1,2) + jump(i,j+1  ,k)%grad(1,2))/4.
                         elseif ( ph(1)*ph(0) .lt. 0.) then
                            nb = 1
                            jmean(1) = (jump(i+1,j,k)%grad(1,2) + jump(i+1,j+1,k)%grad(1,2) &
                                 +jump(i,j  ,k)%grad(1,2) + jump(i,j+1  ,k)%grad(1,2))/4.
                         endif

                         theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                         jmean(0) = (jump(i,j,k)%grad(1,2) + jump(i,j+1,k)%grad(1,2))/2.
                         Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                         Ji = Ji*(miu_pls-miu_mns)

                         miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                         dudy(-1) = dudym
                         dudy( 1) = dudyp

                         dudy( nb) = miuh*dudy(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                         dudy(-nb) = miu_pls*dudy(-nb)

                         RD2 = (dudy(1)-dudy(-1))/dx  ! d/dx(miu*du/dy)

                      endif

                   endif
                endif


                !- 2nd y-derivative -!

                ph(-1) = phi(i,j-1,k)
                ph( 1) = phi(i,j+1,k)

                if (minval(ph) .gt. 0.) then      ! in fluid 1

                   RD3 = miu1*RD3

                elseif (maxval(ph) .lt. 0.) then  ! in fluid 2

                   RD3 = miu2*RD3

                else

                   if (ph(0) .gt. 0.) then
                      miu_pls = miu1
                      miu_mns = miu2
                   else
                      miu_pls = miu2
                      miu_mns = miu1
                   endif

                   if ( ph(-1)*ph(0) .le. 0.) then
                      nb = -1
                      jmean(-1) = (jump(i,j,k)%grad(2,2) + jump(i,j-1,k)%grad(2,2))/2.
                   else
                      nb = 1
                      jmean(1) = (jump(i,j+1,k)%grad(2,2) + jump(i,j+2,k)%grad(2,2))/2.
                   endif

                   theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                   jmean(0) = (jump(i,j,k)%grad(2,2) + jump(i,j+1,k)%grad(2,2))/2.
                   Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                   Ji = Ji*(miu_pls-miu_mns)

                   miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                   dvdy(-1) = dvdym
                   dvdy( 1) = dvdyp

                   dvdy( nb) = miuh*dvdy(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                   dvdy(-nb) = miu_pls*dvdy(-nb)

                   RD3 = 2.*(dvdy(1)-dvdy(-1))/dy  ! d/dy(2.*miu*dv/dy)

                endif

                !- 2nd z-derivative -!

                ph(-1) = phi(i,j,k-1)
                ph( 1) = phi(i,j,k+1)

                if (minval(ph) .gt. 0.) then      ! in fluid 1

                   RD4 = miu1*RD4
                   RD5 = 0.

                elseif (maxval(ph) .lt. 0.) then  ! in fluid 2

                   RD4 = miu2*RD4
                   RD5 = 0.

                else

                   if (ph(0) .gt. 0.) then
                      miu_pls = miu1
                      miu_mns = miu2
                   else
                      miu_pls = miu2
                      miu_mns = miu1
                   endif

                   if ( ph(-1)*ph(0) .le. 0.) then
                      nb = -1
                      jmean(-1) = (jump(i,j,k-1)%grad(2,3) + jump(i,j+1,k-1)%grad(2,3))/2.
                   else
                      nb = 1
                      jmean(1) = (jump(i,j,k+1)%grad(2,3) + jump(i,j+1,k+1)%grad(2,3))/2.
                   endif

                   theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                   jmean(0) = (jump(i,j,k)%grad(2,3) + jump(i,j+1,k)%grad(2,3))/2.
                   Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                   Ji = Ji*(miu_pls-miu_mns)

                   miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                   dvdz(-1) = dvdzm
                   dvdz( 1) = dvdzp

                   dvdz( nb) = miuh*dvdz(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                   dvdz(-nb) = miu_pls*dvdz(-nb)

                   RD4 = (dvdz(1)-dvdz(-1))/dz  ! d/dz(miu*dv/dz)

                   ! Below L and R locates on the vertices, separated by dz (not 2dz)

                   ph(-1) = (phi(i,j,k)+phi(i,j,k-1))/2.
                   ph( 1) = (phi(i,j,k)+phi(i,j,k+1))/2.

                   if (minval(ph) .gt. 0.) then  ! in fluid 1

                      RD5 = miu1*(dwdyp-dwdym)/dz

                   elseif ( maxval(ph) .lt. 0. ) then  ! in fluid 2

                      RD5 = miu2*(dwdyp-dwdym)/dz

                   else  ! between L and R

                      if ( ph(-1)*ph(0) .le. 0.) then
                         nb = -1
                         jmean(-1) = (jump(i,j,k-1)%grad(2,3) + jump(i,j+1,k-1)%grad(2,3) &
                              +jump(i,j  ,k)%grad(2,3) + jump(i,j+1  ,k)%grad(2,3))/4.
                      elseif ( ph(1)*ph(0) .lt. 0.) then
                         nb = 1
                         jmean(1) = (jump(i,j,k+1)%grad(2,3) + jump(i,j+1,k+1)%grad(2,3) &
                              +jump(i,j  ,k)%grad(2,3) + jump(i,j+1  ,k)%grad(2,3))/4.
                      endif

                      theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                      jmean(0) = (jump(i,j,k)%grad(2,3) + jump(i,j+1,k)%grad(2,3))/2.
                      Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                      Ji = Ji*(miu_pls-miu_mns)

                      miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                      dwdy(-1) = dwdym
                      dwdy( 1) = dwdyp

                      dwdy( nb) = miuh*dwdy(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                      dwdy(-nb) = miu_pls*dwdy(-nb)

                      RD5 = (dwdy(1)-dwdy(-1))/dz  ! d/dz(miu*dw/dy)

                   endif
                endif

                RDv(i,j,k) = (RD1+RD2+RD3+RD4+RD5)/rho_v(i,j,k)

             endif

          enddo
       enddo
    enddo


    !-- w-component --! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    call ls_c2f('w-grid',phi)  ! get level set on the w grid
    ph = 0.                    ! level set buffer
    s1 = 3.*dz                 ! tube half width

    dwdx = 0.
    dudz = 0.
    dwdy = 0.
    dvdz = 0.
    dwdz = 0.

    do k = 1,kmax
       do j = 1,jmax
          do i = 1,imax

             ip = i + 1
             jp = j + 1
             kp = k + 1
             im = i - 1
             jm = j - 1
             km = k - 1

             ! continuous first derivatives

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

             RD1 = (dwdxp - dwdxm)/dx  ! d/dx(dw/dx)
             RD3 = (dwdyp - dwdym)/dy  ! d/dy(dw/dy)
             RD5 = (dwdzp - dwdzm)/dz  ! d/dz(dw/dz)         

             ph(0) = phi(i,j,k)

             if (ph(0) .gt. s1) then       ! (most likely) safely in fluid 1

                RDw(i,j,k) = miu1*(RD1+RD3+RD5)/rho1

             elseif (ph(0) .lt. -s1) then  ! (most likely) safely in fluid 2

                RDw(i,j,k) = miu2*(RD1+RD3+RD5)/rho2

             else              

                !- 2nd x-derivative -!

                if (TwoD) then

                   RD1 = 0.
                   RD2 = 0.

                else

                   ph(-1) = phi(i-1,j,k)
                   ph( 1) = phi(i+1,j,k)

                   if (minval(ph) .gt. 0.) then      ! in fluid 1

                      RD1 = miu1*RD1
                      RD2 = 0.

                   elseif (maxval(ph) .lt. 0.) then  ! in fluid 2

                      RD1 = miu2*RD1
                      RD2 = 0.

                   else

                      if (ph(0) .gt. 0.) then
                         miu_pls = miu1
                         miu_mns = miu2
                      else
                         miu_pls = miu2
                         miu_mns = miu1
                      endif

                      if ( ph(-1)*ph(0) .le. 0.) then
                         nb = -1
                         jmean(-1) = (jump(i-1,j,k)%grad(3,1) + jump(i-1,j,k+1)%grad(3,1))/2.!#
                      else
                         nb = 1
                         jmean(1) = (jump(i+1,j,k)%grad(3,1) + jump(i+1,j,k+1)%grad(3,1))/2.!#
                      endif

                      theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                      jmean(0) = (jump(i,j,k)%grad(3,1) + jump(i,j,k+1)%grad(3,1))/2.
                      Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                      Ji = Ji*(miu_pls-miu_mns)

                      miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                      dwdx(-1) = dwdxm
                      dwdx( 1) = dwdxp

                      dwdx( nb) = miuh*dwdx(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                      dwdx(-nb) = miu_pls*dwdx(-nb)

                      RD1 = (dwdx(1)-dwdx(-1))/dx  ! d/dx(miu*dw/dx)

                      ! Below L and R locates on the vertices

                      ph(-1) = (phi(i,j,k)+phi(i-1,j,k))/2.
                      ph( 1) = (phi(i,j,k)+phi(i+1,j,k))/2.

                      if (minval(ph) .gt. 0.) then  ! in fluid 1

                         RD2 = miu1*(dudzp-dudzm)/dx

                      elseif ( maxval(ph) .lt. 0. ) then  ! in fluid 2

                         RD2 = miu2*(dudzp-dudzm)/dx

                      else  ! between L and R

                         if ( ph(-1)*ph(0) .le. 0.) then
                            nb = -1
                            jmean(-1) = (jump(i-1,j,k)%grad(1,3) + jump(i-1,j,k+1)%grad(1,3) &
                                 +jump(i,j  ,k)%grad(1,3) + jump(i,j  ,k+1)%grad(1,3))/4.
                         elseif ( ph(1)*ph(0) .lt. 0.) then
                            nb = 1
                            jmean(1) = (jump(i+1,j,k)%grad(1,3) + jump(i+1,j,k+1)%grad(1,3) &
                                 +jump(i,j  ,k)%grad(1,3) + jump(i,j  ,k+1)%grad(1,3))/4.
                         endif

                         theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                         jmean(0) = (jump(i,j,k)%grad(1,3) + jump(i,j,k+1)%grad(1,3))/2.
                         Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                         Ji = Ji*(miu_pls-miu_mns)

                         miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                         dudz(-1) = dudzm
                         dudz( 1) = dudzp

                         dudz( nb) = miuh*dudz(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                         dudz(-nb) = miu_pls*dudz(-nb)

                         RD2 = (dudz(1)-dudz(-1))/dx  ! d/dx(miu*du/dz)

                      endif

                   endif
                endif


                !- 2nd y-derivative -!

                ph(-1) = phi(i,j-1,k)
                ph( 1) = phi(i,j+1,k)

                if (minval(ph) .gt. 0.) then      ! in fluid 1

                   RD3 = miu1*RD3
                   RD4 = 0.

                elseif (maxval(ph) .lt. 0.) then  ! in fluid 2

                   RD3 = miu2*RD3
                   RD4 = 0.

                else

                   if (ph(0) .gt. 0.) then
                      miu_pls = miu1
                      miu_mns = miu2
                   else
                      miu_pls = miu2
                      miu_mns = miu1
                   endif

                   if ( ph(-1)*ph(0) .le. 0.) then
                      nb = -1
                      jmean(-1) = (jump(i,j-1,k)%grad(3,2) + jump(i,j-1,k+1)%grad(3,2))/2.
                   else
                      nb = 1
                      jmean(1) = (jump(i,j+1,k)%grad(3,2) + jump(i,j+1,k+1)%grad(3,2))/2.
                   endif

                   theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                   jmean(0) = (jump(i,j,k)%grad(3,2) + jump(i,j,k+1)%grad(3,2))/2.
                   Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                   Ji = Ji*(miu_pls-miu_mns)

                   miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                   dwdy(-1) = dwdym
                   dwdy( 1) = dwdyp

                   dwdy( nb) = miuh*dwdy(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                   dwdy(-nb) = miu_pls*dwdy(-nb)

                   RD3 = (dwdy(1)-dwdy(-1))/dy  ! d/dy(miu*dw/dy)

                   ! Below L and R locates on the vertices

                   ph(-1) = (phi(i,j,k)+phi(i,j-1,k))/2.
                   ph( 1) = (phi(i,j,k)+phi(i,j+1,k))/2.

                   if (minval(ph) .gt. 0.) then  ! in fluid 1

                      RD4 = miu1*(dvdzp-dvdzm)/dy

                   elseif ( maxval(ph) .lt. 0. ) then  ! in fluid 2

                      RD4 = miu2*(dvdzp-dvdzm)/dy

                   else  ! between L and R

                      if ( ph(-1)*ph(0) .le. 0.) then
                         nb = -1
                         jmean(-1) = (jump(i,j-1,k)%grad(2,3) + jump(i,j-1,k+1)%grad(2,3) &
                              +jump(i,j  ,k)%grad(2,3) + jump(i,j  ,k+1)%grad(2,3))/4.
                      elseif ( ph(1)*ph(0) .lt. 0.) then
                         nb = 1
                         jmean(1) = (jump(i,j+1,k)%grad(2,3) + jump(i,j+1,k+1)%grad(2,3) &
                              +jump(i,j  ,k)%grad(2,3) + jump(i,j  ,k+1)%grad(2,3))/4.
                      endif

                      theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                      jmean(0) = (jump(i,j,k)%grad(2,3) + jump(i,j,k+1)%grad(2,3))/2.
                      Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                      Ji = Ji*(miu_pls-miu_mns)

                      miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                      dvdz(-1) = dvdzm
                      dvdz( 1) = dvdzp

                      dvdz( nb) = miuh*dvdz(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                      dvdz(-nb) = miu_pls*dvdz(-nb)

                      RD4 = (dvdz(1)-dvdz(-1))/dy  ! d/dy(miu*dv/dz)

                   endif

                endif


                !- 2nd z-derivative -!

                ph(-1) = phi(i,j,k-1)
                ph( 1) = phi(i,j,k+1)

                if (minval(ph) .gt. 0.) then      ! in fluid 1

                   RD5 = miu1*RD5

                elseif (maxval(ph) .lt. 0.) then  ! in fluid 2

                   RD5 = miu2*RD5

                else

                   if (ph(0) .gt. 0.) then
                      miu_pls = miu1
                      miu_mns = miu2
                   else
                      miu_pls = miu2
                      miu_mns = miu1
                   endif

                   if ( ph(-1)*ph(0) .le. 0.) then
                      nb = -1
                      jmean(-1) = (jump(i,j,k)%grad(3,3) + jump(i,j,k-1)%grad(3,3))/2.
                   else
                      nb = 1
                      jmean(1) = (jump(i,j,k+1)%grad(3,3) + jump(i,j,k+2)%grad(3,3))/2.
                   endif

                   theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                   jmean(0) = (jump(i,j,k)%grad(3,3) + jump(i,j,k+1)%grad(3,3))/2.
                   Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                   Ji = Ji*(miu_pls-miu_mns)

                   miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                   dwdz(-1) = dwdzm
                   dwdz( 1) = dwdzp

                   dwdz( nb) = miuh*dwdz(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                   dwdz(-nb) = miu_pls*dwdz(-nb)

                   RD5 = 2.*(dwdz(1)-dwdz(-1))/dz  ! d/dz(2.*miu*dw/dz)

                endif

                RDw(i,j,k) = (RD1+RD2+RD3+RD4+RD5)/rho_w(i,j,k)

             endif

          enddo
       enddo
    enddo


    !-- u-component --! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    if (TwoD) then
       RDu = 0.
       return
    endif

    call ls_c2f('u-grid',phi)  ! get level set on the u grid
    ph = 0.                    ! level set buffer
    s1 = 3.*dz                 ! tube half width

    dudx = 0.
    dudy = 0.
    dvdx = 0.
    dudz = 0.
    dwdx = 0.


    do k = 1,kmax
       do j = 1,jmax
          do i = 1,imax

             ip = i + 1
             jp = j + 1
             kp = k + 1
             im = i - 1
             jm = j - 1
             km = k - 1

             ! continuous first derivatives

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

             RD1 = (dudxp -dudxm)/dx  ! d/dx(du/dx)
             RD2 = (dudyp -dudym)/dy  ! d/dy(du/dy)
             RD4 = (dudzp -dudzm)/dz  ! d/dz(du/dz)

             ph(0) = phi(i,j,k)

             if (ph(0) .gt. s1) then       ! (most likely) safely in fluid 1

                RDu(i,j,k) = miu1*(RD1+RD2+RD4)/rho1

             elseif (ph(0) .lt. -s1) then  ! (most likely) safely in fluid 2

                RDu(i,j,k) = miu2*(RD1+RD2+RD4)/rho2

             else

                !- 2nd x-derivative -!

                ph(-1) = phi(i-1,j,k)
                ph( 1) = phi(i+1,j,k)

                if (minval(ph) .gt. 0.) then      ! in fluid 1

                   RD1 = miu1*RD1

                elseif (maxval(ph) .lt. 0.) then  ! in fluid 2

                   RD1 = miu2*RD1

                else

                   if (ph(0) .gt. 0.) then
                      miu_pls = miu1
                      miu_mns = miu2
                   else
                      miu_pls = miu2
                      miu_mns = miu1
                   endif

                   if ( ph(-1)*ph(0) .le. 0.) then
                      nb = -1
                      jmean(-1) = (jump(i,j,k)%grad(1,1) + jump(i-1,j,k)%grad(1,1))/2.
                   else
                      nb = 1
                      jmean(1) = (jump(i+1,j,k)%grad(1,1) + jump(i+2,j,k)%grad(1,1))/2.
                   endif

                   theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                   jmean(0) = (jump(i,j,k)%grad(1,1) + jump(i+1,j,k)%grad(1,1))/2.
                   Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                   Ji = Ji*(miu_pls-miu_mns)

                   miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                   dudx(-1) = dudxm
                   dudx( 1) = dudxp

                   dudx( nb) = miuh*dudx(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                   dudx(-nb) = miu_pls*dudx(-nb)

                   RD1 = 2.*(dudx(1)-dudx(-1))/dx  ! d/dx(2.*miu*du/dx)

                endif


                !- 2nd y-derivative -!

                ph(-1) = phi(i,j-1,k)
                ph( 1) = phi(i,j+1,k)

                if (minval(ph) .gt. 0.) then      ! in fluid 1

                   RD2 = miu1*RD2
                   RD3 = 0.

                elseif (maxval(ph) .lt. 0.) then  ! in fluid 2

                   RD2 = miu2*RD2
                   RD3 = 0.

                else

                   if (ph(0) .gt. 0.) then
                      miu_pls = miu1
                      miu_mns = miu2
                   else
                      miu_pls = miu2
                      miu_mns = miu1
                   endif

                   if ( ph(-1)*ph(0) .le. 0.) then
                      nb = -1
                      jmean(-1) = (jump(i,j-1,k)%grad(1,2) + jump(i+1,j-1,k)%grad(1,2))/2.
                   else
                      nb = 1
                      jmean(1) = (jump(i,j+1,k)%grad(1,2) + jump(i+1,j+1,k)%grad(1,2))/2.
                   endif

                   theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                   jmean(0) = (jump(i,j,k)%grad(1,2) + jump(i+1,j,k)%grad(1,2))/2.
                   Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                   Ji = Ji*(miu_pls-miu_mns)

                   miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                   dudy(-1) = dudym
                   dudy( 1) = dudyp

                   dudy( nb) = miuh*dudy(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                   dudy(-nb) = miu_pls*dudy(-nb)

                   RD2 = (dudy(1)-dudy(-1))/dy  ! d/dy(miu*du/dy)

                   ! Below L and R locates on the vertices

                   ph(-1) = (phi(i,j,k)+phi(i,j-1,k))/2.
                   ph( 1) = (phi(i,j,k)+phi(i,j+1,k))/2.

                   if (minval(ph) .gt. 0.) then  ! in fluid 1

                      RD3 = miu1*(dvdxp-dvdxm)/dy

                   elseif ( maxval(ph) .lt. 0. ) then  ! in fluid 2

                      RD3 = miu2*(dvdxp-dvdxm)/dy

                   else  ! between L and R

                      if ( ph(-1)*ph(0) .le. 0.) then
                         nb = -1
                         jmean(-1) = (jump(i,j-1,k)%grad(2,1) + jump(i+1,j-1,k)%grad(2,1) &
                              +jump(i,j  ,k)%grad(2,1) + jump(i+1,j  ,k)%grad(2,1))/4.
                      elseif ( ph(1)*ph(0) .lt. 0.) then
                         nb = 1
                         jmean(1) = (jump(i,j+1,k)%grad(2,1) + jump(i+1,j+1,k)%grad(2,1) &
                              +jump(i,j  ,k)%grad(2,1) + jump(i+1,j  ,k)%grad(2,1))/4.
                      endif

                      theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                      jmean(0) = (jump(i,j,k)%grad(2,1) + jump(i+1,j,k)%grad(2,1))/2.
                      Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                      Ji = Ji*(miu_pls-miu_mns)

                      miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                      dvdx(-1) = dvdxm
                      dvdx( 1) = dvdxp

                      dvdx( nb) = miuh*dvdx(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                      dvdx(-nb) = miu_pls*dvdx(-nb)

                      RD3 = (dvdx(1)-dvdx(-1))/dy  ! d/dy(miu*dv/dx)

                   endif
                endif


                !- 2nd z-derivative -!

                ph(-1) = phi(i,j,k-1)
                ph( 1) = phi(i,j,k+1)

                if (minval(ph) .gt. 0.) then      ! in fluid 1

                   RD4 = miu1*RD4
                   RD5 = 0.

                elseif (maxval(ph) .lt. 0.) then  ! in fluid 2

                   RD4 = miu2*RD4
                   RD5 = 0.

                else

                   if (ph(0) .gt. 0.) then
                      miu_pls = miu1
                      miu_mns = miu2
                   else
                      miu_pls = miu2
                      miu_mns = miu1
                   endif

                   if ( ph(-1)*ph(0) .le. 0.) then
                      nb = -1
                      jmean(-1) = (jump(i,j,k-1)%grad(1,3) + jump(i+1,j,k-1)%grad(1,3))/2.
                   else
                      nb = 1
                      jmean(1) = (jump(i,j,k+1)%grad(1,3) + jump(i+1,j,k+1)%grad(1,3))/2.
                   endif

                   theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                   jmean(0) = (jump(i,j,k)%grad(1,3) + jump(i+1,j,k)%grad(1,3))/2.
                   Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                   Ji = Ji*(miu_pls-miu_mns)

                   miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                   dudz(-1) = dudzm
                   dudz( 1) = dudzp

                   dudz( nb) = miuh*dudz(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                   dudz(-nb) = miu_pls*dudz(-nb)

                   RD4 = (dudz(1)-dudz(-1))/dz  ! d/dz(miu*du/dz)

                   ! Below L and R locates on the vertices

                   ph(-1) = (phi(i,j,k)+phi(i,j,k-1))/2.
                   ph( 1) = (phi(i,j,k)+phi(i,j,k+1))/2.

                   if (minval(ph) .gt. 0.) then  ! in fluid 1

                      RD5 = miu1*(dwdxp-dwdxm)/dz

                   elseif ( maxval(ph) .lt. 0. ) then  ! in fluid 2

                      RD5 = miu2*(dwdxp-dwdxm)/dz

                   else  ! between L and R

                      if ( ph(-1)*ph(0) .le. 0.) then
                         nb = -1
                         jmean(-1) = (jump(i,j,k-1)%grad(3,1) + jump(i+1,j,k-1)%grad(3,1) &
                              +jump(i,j  ,k)%grad(3,1) + jump(i+1,j  ,k)%grad(3,1))/4.
                      elseif ( ph(1)*ph(0) .lt. 0.) then
                         nb = 1
                         jmean(1) = (jump(i,j,k+1)%grad(3,1) + jump(i+1,j,k+1)%grad(3,1) &
                              +jump(i,j  ,k)%grad(3,1) + jump(i+1,j  ,k)%grad(3,1))/4.
                      endif

                      theta = abs(ph(nb))/( abs(ph(nb)) +abs(ph(0)) )  ! linear interpln.

                      jmean(0) = (jump(i,j,k)%grad(3,1) + jump(i+1,j,k)%grad(3,1))/2.
                      Ji = theta*jmean(0) +(1.-theta)*jmean(nb)  ! interpolated jump at the interface
                      Ji = Ji*(miu_pls-miu_mns)

                      miuh = miu_pls*miu_mns/(miu_pls*theta+miu_mns*(1.-theta))  ! effective viscosity

                      dwdx(-1) = dwdxm
                      dwdx( 1) = dwdxp

                      dwdx( nb) = miuh*dwdx(nb) + miuh/miu_mns*theta*Ji  ! interface btw 0 and nb
                      dwdx(-nb) = miu_pls*dwdx(-nb)

                      RD5 = (dwdx(1)-dwdx(-1))/dz  ! d/dz(miu*dw/dx)

                   endif
                endif

                RDu(i,j,k) = (RD1+RD2+RD3+RD4+RD5)/rho_u(i,j,k)

             endif

          enddo
       enddo
    enddo



    return
  end subroutine gfm4diff



end module mod_gfm
