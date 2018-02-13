module mod_chkdiv
  
  implicit none
  
  private
  public chkdiv
  
contains
  
  subroutine chkdiv(istep)
    
    use mpi
    use mod_param
    use mod_common
    use mod_common_mpi

    integer, intent(in) :: istep
    integer :: i,j,k,im,jm,km
    real :: div,divtot,divtot_all,divmax(2),divmax_all(2)
    
    divmax(1) = 0.
    divmax(2) = 1.*myid
    divtot = 0.
    im = 0
    jm = 0
    km = 0
    
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             
             div =( wnew(i,j,k)-wnew(i,j,k-1) )*dzi + &
                  ( vnew(i,j,k)-vnew(i,j-1,k) )*dyi + &
                  ( unew(i,j,k)-unew(i-1,j,k) )*dxi

             divtot = divtot + div
             
             div = abs(div)
             if(div.gt.divmax(1)) then
                divmax(1) = div
                im = i
                jm = j
                km = k
             endif

          enddo
       enddo
    enddo

    call mpi_allreduce(divtot,divtot_all,1,mpi_real8,mpi_sum,MPI_COMM_WORLD,error)
    call mpi_allreduce(divmax,divmax_all,1,mpi_2double_precision,mpi_maxloc,MPI_COMM_WORLD,error)

    if (myid .eq. int(divmax_all(2))) then
       write(6,111) zstart(1)-1+im, zstart(2)-1+jm,km
       write(6,222) divtot_all,divmax_all(1),int(divmax_all(2))
111    format('! Maximal divergence at i = ',I5,' j = ', I5,' k = ',I5)
222    format('! Divergence: Tot = ',ES13.6,' L_infty = ',ES13.6,' Rank = ',I5)
    endif
    
    if (divmax_all(1) .gt. 1e-12) then
       if (istep .ne. 0) then
          if (myid .eq. divmax_all(2)) then
             write(6,*) 'Error: Magnitude of the maximal divergence > 1e-12.'
             write(6,*) 'Carefully check the following values, or try a smaller timestep.'
             write(6,*) '****************************************************************'
             write(6,'(A,3ES16.8)') 'u = ',unew(im-1,jm,km),unew(im,jm,km),unew(im,jm,km)-unew(im-1,jm,km)
             write(6,'(A,3ES16.8)') 'v = ',vnew(im,jm-1,km),vnew(im,jm,km),vnew(im,jm,km)-vnew(im,jm-1,km)
             write(6,'(A,3ES16.8)') 'w = ',wnew(im,jm,km-1),wnew(im,jm,km),wnew(im,jm,km)-wnew(im,jm,km-1)
             write(6,'(A,1ES16.8)') 'phi(i,j,k) = ',slset(im,jm,km)/dz
             write(6,'(A,2ES16.8)') 'phi(i-/+1) = ',slset(im-1,jm,km)/dz,slset(im+1,jm,km)/dz
             write(6,'(A,2ES16.8)') 'phi(j-/+1) = ',slset(im,jm-1,km)/dz,slset(im,jm+1,km)/dz
             write(6,'(A,2ES16.8)') 'phi(k-/+1) = ',slset(im,jm,km-1)/dz,slset(im,jm,km+1)/dz
             write(6,'(A,3ES16.8)') 'p_jump x, y, z = ',p_x(im,jm,km),p_y(im,jm,km),p_z(im,jm,km)
             write(6,'(A,3ES16.8)') 'p_jump x -1:+1 = ',p_x(im-1,jm,km),p_y(im,jm,km),p_y(im+1,jm,km)
             write(6,'(A,3ES16.8)') 'p_jump y -1:+1 = ',p_y(im,jm-1,km),p_y(im,jm,km),p_y(im,jm+1,km)
             write(6,'(A,3ES16.8)') 'p_jump z -1:+1 = ',p_z(im,jm,km-1),p_z(im,jm,km),p_z(im,jm,km+1)
             write(6,'(A,1ES16.8)') 'curv(i,j,k) = ',curv_cmn(im,jm,km)
             write(6,'(A,2ES16.8)') 'curv(i-/+1) = ',curv_cmn(im-1,jm,km),curv_cmn(im+1,jm,km)
             write(6,'(A,2ES16.8)') 'curv(j-/+1) = ',curv_cmn(im,jm-1,km),curv_cmn(im,jm+1,km)
             write(6,'(A,2ES16.8)') 'curv(k-/+1) = ',curv_cmn(im,jm,km-1),curv_cmn(im,jm,km+1)
             write(6,*) '***************************************************************'
             write(6,*) 'Program aborted.'
          endif
          call mpi_finalize(error)
          stop
       endif
    endif
    
    return
  end subroutine chkdiv
  
end module mod_chkdiv
