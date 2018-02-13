!                  ___                         ___     
!    ___          /  /\                       /  /\    
!   /  /\        /  /:/                      /  /:/_   
!  /  /:/       /  /:/       ___     ___    /  /:/ /\  
! /__/::\      /  /:/  ___  /__/\   /  /\  /  /:/ /::\ 
! \__\/\:\__  /__/:/  /  /\ \  \:\ /  /:/ /__/:/ /:/\:\
!    \  \:\/\ \  \:\ /  /:/  \  \:\  /:/  \  \:\/:/~/:/
!     \__\::/  \  \:\  /:/    \  \:\/:/    \  \::/ /:/ 
!     /__/:/    \  \:\/:/      \  \::/      \__\/ /:/  
!     \__\/      \  \::/        \__\/         /__/:/   
!                 \__\/                       \__\/    
!                                                          
!                                                       
! ICLS - Interface-Correction Level Set                 
!------------------------------------------------------

program ICLS

  use mod_param
  use mod_param_ls
  use mod_common
  use mod_common_mpi
  use mod_initmpi
  use mod_init
  use mod_bound
  use mod_chkdiv
  use mod_chkdt
  use mod_loadflds
  use mod_mom
  use mod_ab
  use mod_levelset
  use mod_interface
  use mod_gfm
  use mod_solver
  use mod_fillps
  use mod_correc
  use mod_output
  use mod_write_vtk
  use mod_stopwatch
  use mod_misc
  
  implicit none

  integer :: i,j,k,l
  integer :: iopt, begin,nstep,istep

  real :: dtmax
  real :: ls_time, ls_time_max, ls_stopwatch_max
  real :: ab2_time, ab2_time_max, ab2_stopwatch_max
  real :: pj_time, pj_time_max, pj_stopwatch_max
  real :: comp_time, comp_time_min,comp_time_max
  real :: stopwatch_min,stopwatch_max
  real :: norm,norm_all,norm_old
  real :: v_bulk_all, massloss, quck, vel_max

  real, dimension(1:imax,1:jmax,1:kmax) :: duold,dvold,dwold
  real, dimension(0:i1,0:j1,0:k1) :: pold,p_jump, curv_restart

  real, pointer, dimension(:,:,:) :: pp2 ! uncomment iff solver2d is used
  real, target,  dimension(0:i1,0:j1,0:k1) :: p
 
 
  !--------------------------------------------!
  !              INITIALIZATION                !
  !--------------------------------------------!
  
  iopt = 0
  if(iopt .eq. 1) then
     call MPI_INIT(error)
     call decomp_2d_init(itot, jtot, ktot, 0,0)
     call MPI_FINALIZE(error)
     stop
  else
     call initmpi
  endif
  
  begin = 0 
  nstep = 4000

  call init_program(begin, duold,dvold,dwold,pold,curv_restart)
  
  pp2 => p(1:imax,1:jmax,1:kmax) ! uncomment iff solver2d is used

  !---- Output initial conditions ----!
  
  call post2d(begin)
  !call write_vtk(begin)  
  
  !---- MPI timer ----!
  
  ls_stopwatch_max = 0.
  ab2_stopwatch_max = 0.
  pj_stopwatch_max = 0.
  stopwatch_min = 0.
  stopwatch_max = 0.
  total_time = 0.


  !--------------------------------------------!
  !                 MAIN LOOP                  !
  !--------------------------------------------!
  
  do istep = begin+1,nstep

     if (myid .eq. 0) write(6,*) 'istep =', istep
     comp_time = mpi_wtime()
     time = time + dt

     !---- Level Set ----!

     if (solveLS) then
        
        call ls_main(istep)
        
        ls_time = mpi_wtime() - comp_time
        call mpi_allreduce(ls_time,ls_time_max,1,mpi_real8,mpi_max,comm_cart,error)
        
     endif

     !---- Navier-Stokes ----!
     
     if (solveNS) then
        
        call average_stream_vel(v_bulk)  ! calculate average streamwise velocity

        !-- Intermediate Velocity --!
        
        ab2_time = mpi_wtime()

        call ab2(duold,dvold,dwold,pold)

        ab2_time = mpi_wtime() - ab2_time
        call mpi_allreduce(ab2_time,ab2_time_max,1,mpi_real8,mpi_max,comm_cart,error)

        !-- RHS of Poisson --!

        pj_time = mpi_wtime()
        
        call fillps(p)
        if (solveLS) then
           call gfm4p(p_jump)
           p = p + p_jump
        endif

        !-- Pressure --!

        select case(BC_in_y)
        case ('Periodic')
           call solver2d(pp2)  ! pp2 => p
           !      call solver1d(p)
        case ('InOut','Walls','OutOut')
           call solver2d_cos(pp2)
        end select
        
        call boundp(p)

        !-- New Velocity --!

        call correc(p)  
        call bounduvw(unew,vnew,wnew)

        !-- Update & Normalize Pressure --!
        
        pold(:,:,:) = pnew(:,:,:) 
        pnew(:,:,:) = p(:,:,:)

        call boundp(pnew)

        select case(BC_in_y)
        case('Periodic','Walls')
           
           k=1 ! near bottom wall
           norm = sum(sum(pnew(:,:,k),2),1)
           call mpi_allreduce(norm,norm_all,1,mpi_real8,mpi_sum,comm_cart,error)
           norm = norm_all/(1.*itot*jtot)
           pnew(:,:,:) = pnew(:,:,:) - norm

           if (istep .eq. begin+1) then
              pold(:,:,:) = pold(:,:,:) - norm  
           else
              pold(:,:,:) = pold(:,:,:) + norm_old - norm  ! re-normalize pold with the new norm
           endif
           norm_old = norm  ! update norm_old after correcting pold

        end select

        pj_time = mpi_wtime() - pj_time
        call mpi_allreduce(pj_time,pj_time_max,1,mpi_real8,mpi_max,comm_cart,error)

        call average_vel_at_center ! update center velocities
        if (stop_at_steady) call check_steady_state(istep)
        
     else
        
        select case (iniu)
        case('sep')
           call serpentine(time)  ! 2D periodic deformation (purely level set)
        end select

     endif


     !--------------------------------------------!
     !                   OUTPUT                   !
     !--------------------------------------------!
     
     if (isnan(vnew(1,1,1))) then
        if (myid .eq. 0) write(6,*) 'Got NaN... Terminate brutally. :('
        stop
     endif

     if (mod(istep,ioutchk) .eq. 0) then
        call chkdiv(istep)
        call chkdt(dtmax)
        if (myid .eq. 0) write(6,'(A,1ES12.4)') '! dt calculated by chkdt =', dtmax
        dt = timestep

        !-- Dispersed file --!

        if (.false.) then
           call disp_cg     ! center of gravity
           call disp_speed  ! average velocity
           
           if (myid .eq. 0) then
             open(61,file=datadir//'dispersed.txt',position='append')
             write(61,'(7ES16.8)') 1.*istep,x_disp_avr,y_disp_avr,z_disp_avr,&
                                            u_disp_avr,v_disp_avr,w_disp_avr
             close(61)
           endif
        endif
        
        !-- Timestep file (quick check) --!
        
        quck = 0.
        !call Prosperetti(quck)
        !call spurious_currents(quck)
        !call Grace(quck)
        !quck = v_bulk
        !quck = maxval(dmin)/dz  ! max distance between droplets
        quck = dmin(1,2)/dz
        !call Taylor(quck)

        call max_vel_at_center(vel_max)

        if (myid .eq. 0) then
           write(6,'(A,1ES12.4,A,1ES12.4)') '! Actual dt              =',dt,'  max CFL =', CFL*vel_max
           write(6,'(A,1ES12.4)') '! Time =', time

           massloss = 0.   
           do l = 1,lmax
              massloss = massloss + (1. - volume(l)/init_volume(l))
           enddo

           open(29,file=datadir//'timestep.txt',position='append')
           write(29,'(6ES16.8)') 1.*istep,time,dt, massloss, quck, maxchk_glb
           close(29)
        endif
     endif
     
     !-- 2D file (mid-yz plane) --!
     
     if (mod((istep),iout2d) .eq. 0) then
        call post2d(istep) 
     endif

     !-- 3D file --!

     if (mod(istep,iout3d) .eq. 0) then
        !      call post3d(istep)
        call write_vtk(istep)   ! read in paraview
     endif

     !-- Restarting file --!
     
     if (mod(istep,ioutfld) .eq. 0) then
        call loadflds(1,istep)  ! flow field and level set field
     endif

     !-- Timing --!
     
     comp_time = MPI_WTIME()-comp_time
     call mpi_allreduce(comp_time,comp_time_min,1,mpi_real8,mpi_min,comm_cart,error)
     call mpi_allreduce(comp_time,comp_time_max,1,mpi_real8,mpi_max,comm_cart,error)

     call stopwatch(istep,'Total', comp_time_min,comp_time_max, stopwatch_min,stopwatch_max)
     call time_breakdown(istep,'LvSet', ls_time_max, ls_stopwatch_max, stopwatch_max)
     call time_breakdown(istep,'exAB2',ab2_time_max,ab2_stopwatch_max, stopwatch_max)
     call time_breakdown(istep,'Projt', pj_time_max, pj_stopwatch_max, stopwatch_max)


  enddo  ! end of main loop


  if(myid .eq. 0) then
     write(6,'(A)') ' '
     write(6,'(A)') '************** Simulation completed. Congrats! \o/ *************'
  end if

  call decomp_2d_finalize
  call MPI_FINALIZE(error)


  stop
end program
