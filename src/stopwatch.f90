module mod_stopwatch

 implicit none

 private
 public stopwatch, time_breakdown
 contains

 subroutine stopwatch(istep,module_name, &
                      comp_time_min,comp_time_max, stpw_min,stpw_max)

   use mod_param
   use mod_common
   use mod_common_mpi
   use mod_initmpi
   
   integer, intent(in) :: istep
   character(len=5), intent(in) :: module_name
   real, intent(in) :: comp_time_min,comp_time_max
   real, intent(inout) :: stpw_min,stpw_max


   stpw_min = stpw_min + comp_time_min
   stpw_max = stpw_max + comp_time_max

   if(myid .eq. 0) then

     if (mod(istep,ioutchk) .eq. 0) then
  
       if (stpw_min .lt. 60.) then
  
         write(6,'(A,A,A,1ES10.4,A)') '! ',module_name,' computing time (minimum) = ', stpw_min,' s'
         write(6,'(A,1ES10.4,A)')               '!                      (maximum) = ', stpw_max,' s'
  
       elseif (stpw_min .lt. 3600.) then
  
         write(6,'(A,A,A,1ES10.4,A)') '! ',module_name,' computing time (minimum) = ', stpw_min/60.,' min'
         write(6,'(A,1ES10.4,A)')               '!                      (maximum) = ', stpw_max/60.,' min'
  
       elseif (stpw_min .lt. 24.*3600.) then
  
         write(6,'(A,A,A,1ES10.4,A)') '! ',module_name,' computing time (minimum) = ', stpw_min/3600.,' hours'
         write(6,'(A,1ES10.4,A)')               '!                      (maximum) = ', stpw_max/3600.,' hours'
  
       else
  
         write(6,'(A,A,A,1ES10.4,A)') '! ',module_name,' computing time (minimum) = ', stpw_min/3600./24.,' days'
         write(6,'(A,1ES10.4,A)')               '!                      (maximum) = ', stpw_max/3600./24.,' days'
  
       endif
  
     endif

   endif

  return
 end subroutine stopwatch



 !~~~~~~~~~~~~~~~~~~~~~~

 subroutine time_breakdown(istep,module_name, &
                             comp_time_max,stpw_max, total_time_max)

   use mod_param
   use mod_common
   use mod_common_mpi
   use mod_initmpi
   
   integer, intent(in) :: istep
   character(len=5), intent(in) :: module_name
   real, intent(in) :: comp_time_max, total_time_max
   real, intent(inout) :: stpw_max
   real :: t_ratio


   stpw_max = stpw_max + comp_time_max
   t_ratio = stpw_max/total_time_max*100.

   if(myid .eq. 0) then

     if (mod(istep,ioutchk) .eq. 0) then
  
       write(6,'(A,A,A,1F5.2,A)') '! ',module_name,' contribution (maximum): ', t_ratio,' %'
  
     endif

   endif

  return
 end subroutine time_breakdown


end module mod_stopwatch
