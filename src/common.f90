module mod_common

 use mod_param

 implicit none

 real ,dimension(0:i1,0:j1,0:k1) :: unew,vnew,wnew,pnew,dudt,dvdt,dwdt
 real ,dimension(0:i1,0:j1,0:k1) :: forcex,forcey,forcez,forces
 real ,dimension(0:i1,0:j1,0:k1) :: jiemian, miu,rho, u_l,v_l,w_l  !#

 real, dimension(1:imax,1:jmax,1:kmax) :: uccnew,vccnew,wccnew  ! cell-center velocities
 real, dimension(1:imax,1:jmax,1:kmax) :: uccold,vccold,wccold  ! cell-center velocities
 real, dimension(1:imax,1:jmax,1:kmax) :: rho_u,rho_v,rho_w     ! cell-face densities

 !---- Level Set ----!
 
 real ,dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: ynew,dydt, fixls    ! 3 grid cells overlaping in all directions
 real ,dimension(-2:i1+2,-2:j1+2,-2:k1+2,1:lmax) :: lvset,reset  ! multiple level set
 real ,dimension(0:i1,0:j1,0:k1) :: slset, curv_cmn              ! single level set and curvature

 type nb_vector
   real :: m, x,y,z
 end type nb_vector

 type (nb_vector), dimension(1:imax,1:jmax,1:kmax) :: normal !# can include curv

 real ,dimension(0:i1,0:j1,0:k1) :: fixhv                  ! static Heaviside (c. static level set)
 real ,dimension(0:i1,0:j1,0:k1,1:lmax) :: hnew_mls, mvof  ! multiple Heaviside or VOF for phase indication
 real, dimension(1:imax,1:jmax,1:kmax) :: cutcell

 integer :: nb_cnt
 integer, dimension(:), allocatable :: nb_i,nb_j,nb_k  ! narrow band index array
 
 !---- Navier-Stokes ----!

 real ,dimension(0:i1,0:k1) :: v_inflow
 
 real, dimension(1:imax,1:jmax,1:kmax) :: surfxold,surfyold,surfzold, surfxnew,surfynew,surfznew  ! surface tension force (CSF)

 real, dimension(0:i1,0:j1,0:k1) :: dvndn  ! normal derivative of normal velocity
 
 real, dimension(0:i1,0:j1,0:k1) :: p_x, p_y, p_z, &       ! pressure jump in x,y,z (GFM)
                                    p_xold, p_yold, p_zold
 
 
 !---- Dispersed Entity ----!
 
 real :: maxchk_glb, &  !# curvature error, temporary
         x_disp_avr, y_disp_avr, z_disp_avr, &  !# average of all droplets
         u_disp_avr, v_disp_avr, w_disp_avr  !# average of all droplets
 
 real, dimension(1:lmax) :: volume, init_volume, &
                            tiji, init_tiji, tijiloss, &     ! mass remedy stuff
                            x_disp, y_disp, z_disp, &        ! Center of Gravity
                            u_disp, v_disp, w_disp           ! Average velocity
 
 real, dimension(lmax,lmax) :: dmin  ! minimal distances between interfaces
 
 
 !---- Miscellaneous ----!
 
 real :: total_time  ! timing
 
 real, dimension(4,7) :: pseudo_inv_Al  ! linear least squares matrix
 real, dimension(10,27) :: pseudo_inv_A, pseudo_inv_Ab_lower,pseudo_inv_Ab_upper ! quadratic least squares matrix

 
 !--------------  My additions end  --------------!
 
 
 real(mytype) :: time,dt

 !#real(mytype) :: wi(itot+15), wj(jtot+15), wk(2*kmax+15)
 !#real, dimension(imax,jmax) :: xyrt
 !#real, dimension(kmax) :: a,b,c
 !#real, dimension(imax,jmax,kmax) :: xzrt
 
 real :: forcextot,forceytot,forceztot
 real :: u_bulk,v_bulk,w_bulk
 real :: wallshearold
 real :: dpdx_sumrk
 real :: rkparalpha
 integer :: rkiter
 


 
end module mod_common



module mod_common_mpi

 use mpi
 use decomp_2d

 implicit none

 integer :: myid,xhalo,yhalo,restartp,rankcw, xhalo_ls, yhalo_ls !# Jan 2016
 integer :: xhalo_jump, yhalo_jump  !# Oct 2016
 integer :: xhalo_1d, yhalo_1d
 integer :: comm_cart
 
 logical periods(3),reorder
 integer error,request,status(MPI_STATUS_SIZE)
 integer right,rightfront,front,leftfront,left,leftback,back,rightback
 integer, dimension(0:8) :: neighbor
 integer, dimension(1:2) :: coords
 real :: boundleftmyid,boundfrontmyid
 
 real :: ORMx, ORMy, ORMz


end module mod_common_mpi



module mod_newtypes

! Definition of new data types
  
  implicit none

  
  type real_3by3_matrix
     real, dimension(3,3) :: grad
  end type real_3by3_matrix

  

end module mod_newtypes
