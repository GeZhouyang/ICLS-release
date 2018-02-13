module mod_param

 use decomp_2d

 implicit none


 !---- Parallelization ----!
 !                         !
 ! 2D (xy) decomposition   !
 
 integer, parameter :: ndims = 2
 integer, dimension(ndims), parameter :: dims = (/1,8/)


 !---- Discretization ----!
 !                        !
 ! wall normal: z or k    !
 ! streamwise:  y or j    !
 ! spanwise:    x or i    !
 
 !integer, parameter :: cellpD = 32  ! cells per diameter (droplet)
 !integer, parameter :: ktot = cellpD*3  ! total number of cells
 !integer, parameter :: jtot = cellpD*4
 !integer, parameter :: itot = cellpD*4

 integer, parameter :: ktot = 96, jtot = ktot, itot = 4, cellpD = 32  ! for some 2D tests

 integer, parameter :: imax = itot/dims(1), jmax = jtot/dims(2), kmax = ktot   ! process internal cells
 integer, parameter :: i1   = imax+1,       j1   = jmax+1,       k1   = kmax+1
 integer, parameter :: it1  = itot+1,       jt1  = jtot+1,       kt1  = ktot+1

 real, parameter :: lz = ktot/(cellpD*1.)  ! domain size
 real, parameter :: ly = jtot/(cellpD*1.)
 real, parameter :: lx = itot/(cellpD*1.)
 
 real, parameter :: dxi = itot/lx, dyi = jtot/ly, dzi = ktot/lz
 real, parameter :: dx  = 1./dxi,  dy  = 1./dyi,  dz  = 1./dzi  ! grid spacing
 
 real, parameter :: NarrowBand_1 = 3.*dz      ! cut-off band (not used)
 real, parameter :: NarrowBand_2 = 6.*dz      ! advection band (level set)


 !---- Constants ----!
 !                   !
 
 real, parameter :: timestep = 5e-4, CFL = timestep/dz
 
 real, parameter :: v_bulk_desired = 1.  ! desired streamwise bulk (mean) velocity
 real, parameter :: pi = acos(-1.), picon = acos(-1.)
 

 !---- Material Properties ----!
 !                             !
 ! non-dimensional             !
 ! reference phase: 1          !
   
 real, parameter :: miu1 = 1., rho1 = 1. ! DON'T change (by definition)
 real, parameter :: miu2 = miu1*1.e-1
 real, parameter :: rho2 = rho1*10.
 
 real :: rho_avr
 real, parameter :: rho0 = min(rho1,rho2)  ! constant density in FFT solver (Dodd & Ferrante)
 

 !---- Non-Dimensionalization ----!
 !                                !
  
 real, parameter :: Re  = 1e1	        ! Reynolds number = rho1*U*L/mu1
 real, parameter :: We  = 1e0      	! Weber number = rho1*U^2*L/sigma
 real, parameter :: iFr = 0.		! Froud number = U^2/gL, iFr is its inverse
 
 real, parameter :: visc = 1./Re	! 1/Re in the Navier-Stokes
 real, parameter :: Ca = We*visc	! Capillary number = sigma/(miu*U)
 real, parameter :: Reb = lz*Re	        ! rescale Re by lz (useless, kept for historical reasons)

 real, dimension(3), parameter :: g_unit = (/0.,0.,0./)
 ! unit vector of gravatational acceleration eg. /0.,0.,-1./
 

 !-----  Output Settings -----!
 !                            !
 ! output ? every iout? steps !
 
 integer, parameter :: ioutchk = 20, iout1d = 1/1
 integer, parameter :: iout2d  = 4000, iout3d = 10000
 integer, parameter :: ioutfld = 50000
 
 
 !---- Flow Field Settings ----!
 !                             !
 
 logical, parameter :: solveNS = .true.
 logical, parameter :: TwoD = .true.  ! 2D in y-z plane (dummy periodicity in x)
 logical, parameter :: ns_np = .true.  ! no-slip/no-penetration in z-dir
 logical, parameter :: constPoi = .true.  ! constant bulk flow (effective only if BC_in_y = Periodic)
 logical, parameter :: Riblet = .false.  ! internal Cartesian walls (stress IBM)
 logical, parameter :: stop_at_steady = .false.  ! check steady state (u,v,w,p) and stop if so
  
 character(len=3), parameter :: iniu = 'poi'
 !                                                          
 ! Initial velocity field:                                  
 !                                                          
 ! 'zer' --> zero velocity everywhere                       
 ! 'uni' --> uniform streamwise flow                        
 ! 'poi' --> channel Poiseuille flow              
 ! 'cou' --> plane Couette flow                    
 ! 'sep' --> Serpentine 2D deformation vortex               
 ! 'rot' --> 2D rigid-body rotation, not working  
 ! '2lc' --> 2-parallel-mixing-layer (w/ wall normal noise) 
 ! 'log' --> logarithmic profile + random noise
 ! 'fld' --> load and initialize from fld_init in datadir
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
 character(len=8), parameter :: BC_in_y = 'Periodic'
 !
 ! Boundary condition in the y direction,        
 ! which determines the pressure solver.         
 !                                               
 ! 'Periodic' --> periodic                       
 ! 'InOut'    --> inflow/outflow                         
 ! 'OutOut'   --> outflow/outflow (open boundaries)                
 ! 'Walls'    --> two walls                         
 ! (By default, x is periodic, z is wall-bounded)
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 character(len=8), parameter :: BC_in_z = 'Dirichl'
 !                                                          
 ! Boundary condition in the z direction:                                  
 !                      
 ! 'Dirichl' --> Dirichlet BC (prescribed velocity)                
 ! 'shr-drv' --> shear-driven (mixed Dirichlet & Neumann)
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 character(len=5), parameter :: NS_conv = 'CNTR2'
 !                                                    
 ! Spatial discretization for the NS convective terms:
 !                                                    
 ! 'CNTR2' --> 2nd order central difference           
 ! 'WENO5' --> 5th order Weighted ENO                 
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
 character(len=3), parameter :: surface_tension_method = 'GFM'
 !
 ! Computation of surface tension:  
 !                                  
 ! 'CSF' --> Continuum Surface Force
 ! 'GFM' --> Ghost Fluid Method
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 character(len=3), parameter :: NS_diff_discretization = 'del'
 !                                                   
 ! Spatial discretization for the NS diffusion terms:
 !                                                   
 ! 'del' --> Dirac delta formulation (some smearing)             
 ! 'GFM' --> Ghost fluid method (inaccurate and unrobust)            
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 
 !---- Level Set Settings ----!
 !                            !
 
 logical, parameter :: solveLS = .true.
 logical, parameter :: re_initialization = .true.
 logical, parameter :: mass_correction = .true.
 logical, parameter :: attractive = .false.  ! MLS-GFM depletion attraction
 logical, parameter :: Contact_line = .false.
 logical, parameter :: open_end_in_y = .false.  ! if not enclosed
 
 integer, parameter :: lmax = 2  ! number of multiple level sets
 
 character(len=3), parameter :: inil = 'mpc'
 !                                                                
 ! Initial level set field:                                       
 !                                                                
 ! 'zal' --> Zalesak's disk
 ! 'elp' --> Distorted ellipse
 ! '2lc' --> 2-parallel-mixing-layer                              
 ! 'sin' --> small amplitute sin-wave                             
 ! 'mpc' --> multiple pancakes (MLS) in y-z plane                 
 ! 'mdp' --> multiple spherical droplets (MLS)                    
 ! 'hdm' --> Hadamard's solution for a spherical droplet          
 ! 'pan','2pc','3pc' --> 1 or 2 or 3 pancake(s) (SLS) in y-z plane
 ! 'drp','2dp','3dp' --> 1 or 2 or 3 spherical droplet(s) (SLS)   
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 !----  Miscenalleous  ----!
 !                         !
 
 character(len=5), parameter :: datadir = 'data/'
 character(len=5), parameter :: partdir = 'part/'


end module mod_param



module mod_param_ls
  !
  ! More settings for level set
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  implicit none

  !-- Advection --!
  !               !
  character(len=3), parameter :: ls_adv_time = 'RK3'
  !                                             
  ! Temporal integration of level set advection:
  !                                             
  ! 'RK2' --> 2nd order mid-point Runge-Kutta   
  ! 'RK3' --> 3rd order SSP Runge-Kutta         
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  character(len=8), parameter :: ls_adv_space = 'HOUC5'
  !                                                  
  ! Spatial discretization of level set gradient:    
  !                                                  
  ! 'HOUC5'    --> 5th order High-Order Upstream Central
  ! 'HOUC5rib' --> HOUC5 for Riblet                  
  ! 'WENO5'    --> 5th order Weighted ENO               
  ! 'mWENO5c'  --> Modified WENO5 (conservative form) 
  ! 'mWENO5nc' --> Modified WENO5 (non-conservative) 
  ! 'semiLag'  --> Semi-Lagrangian (NOT recommended)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !-- Reinitialization --!
  !                      !
  integer, parameter :: reinit_step = 20
  ! Frequency of reinitialization (depends on CFL, marching speed = 1)

  integer, parameter :: reinit_iter_max = 1
  ! Maximum number of reinitialization iterations (pseudo time steps)

  character(len=3), parameter :: ls_reinit_time = 'RK3'
  !                                                    
  ! Temporal integration of level set reinitialization:
  !                                                                                    
  ! 'EU1' --> 1st order forward Euler                   
  ! 'RK2' --> 2nd order mid-point Runge-Kutta          
  ! 'RK3' --> 3rd order SSP Runge-Kutta                
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  character(len=6), parameter :: ls_reinit_space = 'WENO5'
  !                                                  
  ! Spatial discretization of level set gradient:    
  !                 
  ! 'RUSM1'  --> 1st order subcell fix by Russo & Smereka
  ! 'RUSM2'  --> 2nd order subcell fix by Russo & Smereka
  ! 'RUSM4'  --> 4th order subcell fix by Russo & Smereka (NOT working)
  ! 'WENO5'  --> 5th order Weighted ENO (classical)
  ! 'WENO5z' --> WENO5-Z (Brazil) 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !-- Mass correction --!
  !                     !
  integer, parameter :: inflate_step = 1*reinit_step
  character(len=4), parameter :: ispeed_dependence = 'curv'
  !                                         !
  ! Dependence of the correction velocity:  !
  !                                         !
  ! 'unfm' --> uniform                      !
  ! 'curv' --> proportional to curvature    !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  character(len=6), parameter :: i_space = 'HOUC5'
  !                                                  
  ! Spatial discretization of level set gradient:    
  !                                                  
  ! 'HOUC5' --> 5th order High-Order Upstream Central
  ! 'WENO5' --> 5th order Weighted ENO               
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

end module mod_param_ls
