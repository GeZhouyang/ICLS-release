module mod_output

use mod_common
use mod_common_mpi
use decomp_2d
use decomp_2d_io
!use mod_post
use mod_zredistribute
implicit none
private
public post1d,post2d,post3d,writeallpart
contains



subroutine post1d(time,istep,nstep)
implicit none

integer, intent(in) :: istep,nstep
real, intent(in) :: time
real, dimension(0:k1) :: um,vm,wm,pm, &
                         ur,vr,wr,pr, &
                         uw
real, dimension(0:k1) :: um_all,vm_all,wm_all,pm_all, &
                         ur_all,vr_all,wr_all,pr_all, &
                         uw_all,vw,vw_all
integer :: i,j,k
real :: dpdx, dpdx_all
!
do k=1,kmax
  um(k) = 0.
  vm(k) = 0.
  wm(k) = 0.
  pm(k) = 0.
  do j=1,jmax
    do i=1,imax
      um(k) = um(k) + unew(i,j,k)
      vm(k) = vm(k) + vnew(i,j,k)
      wm(k) = wm(k) + wnew(i,j,k)
      pm(k) = pm(k) + pnew(i,j,k)
    enddo
  enddo
enddo
!
! no-slip and no-penetration b.c.'s
!
um(0)    = -um(1)
vm(0)    = -vm(1)
wm(0)    =  0.
pm(0)    =  pm(1)      !Neumann
um(k1)   = -um(kmax)
vm(k1)   = -vm(kmax)
wm(kmax) =  0.
wm(k1)   =  wm(kmax-1) !dw/dz=0 (not used in code)
pm(k1)   =  pm(kmax)   !Neumann

call mpi_allreduce(um,um_all,k1+1,mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(vm,vm_all,k1+1,mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(wm,wm_all,k1+1,mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(pm,pm_all,k1+1,mpi_real8,mpi_sum,comm_cart,error)

um(:) = um_all(:)/(1.*itot*jtot)
vm(:) = vm_all(:)/(1.*itot*jtot)
wm(:) = wm_all(:)/(1.*itot*jtot)
pm(:) = pm_all(:)/(1.*itot*jtot)
!
! rms of velocities
!
do k=1,kmax
  ur(k) = 0.
  vr(k) = 0.
  wr(k) = 0.
  pr(k) = 0.
  do j=1,jmax
    do i=1,imax
      ur(k) = ur(k) + (unew(i,j,k)-um(k))**2
      vr(k) = vr(k) + (vnew(i,j,k)-vm(k))**2
      wr(k) = wr(k) + (wnew(i,j,k)-wm(k))**2
      pr(k) = pr(k) + (pnew(i,j,k)-pm(k))**2
    enddo
  enddo
enddo
!No-slip and no-penetration b.c.'s
ur(0)    = -ur(1)
vr(0)    = -vr(1)
wr(0)    =  0.
pr(0)    =  pr(1)        ! Neumann
ur(k1)   = -ur(kmax)
vr(k1)   = -vr(kmax)
wr(kmax) =  0.
wr(k1)   =  wr(kmax-1)   ! dw/dz=0
pr(k1)   =  pr(kmax)     ! Neumann
call mpi_reduce(ur,ur_all,k1+1,mpi_real8,mpi_sum,0,comm_cart,error)
call mpi_reduce(vr,vr_all,k1+1,mpi_real8,mpi_sum,0,comm_cart,error)
call mpi_reduce(wr,wr_all,k1+1,mpi_real8,mpi_sum,0,comm_cart,error)
call mpi_reduce(pr,pr_all,k1+1,mpi_real8,mpi_sum,0,comm_cart,error)
ur(:) = ur_all(:)/(1.*itot*jtot)
vr(:) = vr_all(:)/(1.*itot*jtot)
wr(:) = wr_all(:)/(1.*itot*jtot)
pr(:) = pr_all(:)/(1.*itot*jtot)
!
do k=1,kmax
  uw(k)=0.
  do j=1,jmax
    do i=1,imax
      uw(k) = uw(k) + 0.25*(unew(i,j,k+1)-um(k+1) + unew(i,j,k)-um(k))* &
                           (wnew(i+1,j,k)-wm(k  ) + wnew(i,j,k)-wm(k))
    enddo
  enddo
enddo

call mpi_allreduce(uw(0),uw_all(0),k1+1,mpi_real8,mpi_sum, &
                   comm_cart,error)
uw(1:kmax) = -uw_all(1:kmax)/(1.*itot*jtot)
!
! BC
!
uw(0)    = 0.
uw(kmax) = 0.
uw(k1)   = uw(kmax-1) !d(uw)/dz = 0 at solid wall
!
do k=1,kmax
  vw(k)=0.
  do j=1,jmax
    do i=1,imax
      vw(k) = vw(k) + 0.25*(vnew(i,j,k+1)-vm(k+1) + vnew(i,j,k)-vm(k))* &
                           (wnew(i,j+1,k)-wm(k  ) + wnew(i,j,k)-wm(k))
    enddo
  enddo
enddo

call mpi_allreduce(vw(0),vw_all(0),k1+1,mpi_real8,mpi_sum, &
                   comm_cart,error)
vw(1:kmax) = -vw_all(1:kmax)/(1.*itot*jtot)
!
! BC
!
vw(0)    = 0.
vw(kmax) = 0.
vw(k1)   = vw(kmax-1) !d(uw)/dz = 0 at solid wall
!
if (myid .eq. 0) then
  open(unit=30,file=datadir//'progress.out')
  write(30,'(A6,E16.8,A8,I7,A5,E16.8)') 'time =', time, &
   ' istep =', istep, ' dt =',dt
  write(30,'(A10,E16.8,A2)') 'progress =',100.*istep/(1.*nstep),' %'
  close(30)
!  !velocity profiles !CREATE ONE FILE PER SAVE
  open(40,file=datadir//'vel_stats.out')
  write(40,*)
  do k=1,kmax
    write(40,'(10E15.7)') (k-0.5)/dzi, (k-0.5)/dzi/visc, um(k), &
                                vm(k),(wm(k)+wm(k-1))*0.5,pm(k), &
                                ur(k),vr(k),0.5*(wr(k)+wr(k-1)),pr(k)
  enddo
  close(40)
!
! various stress profiles
!
  open(42,file=datadir//'str_stats.out')
  write(42,*)
  do k=0,kmax
    write(42,'(5E15.7)') k/dzi, k/dzi/visc, &
      visc * ( vm(k+1)-vm(k) )*dzi,  vw(k), &
      visc * ( vm(k+1)-vm(k) )*dzi + vw(k)
  enddo
  close(42)
endif


!
return
end subroutine post1d
!




subroutine post2d(istep)
implicit none

integer, intent(in) :: istep
real, dimension(imax,jmax,kmax) :: aux
real    :: u4(1:jmax,1:kmax)
real    :: v4(1:jmax,1:kmax)
real    :: w4(1:jmax,1:kmax)
real    :: p4(1:jmax,1:kmax)
real    :: y4(1:jmax,1:kmax)
real    :: h4(1:jmax,1:kmax) !#
real    :: cc4(1:jmax,1:kmax) !#
real    :: indie4(1:jmax,1:kmax) !#

real    :: u4r(0:(jtot+1),0:(kmax+1))
real    :: v4r(0:(jtot+1),0:(kmax+1))
real    :: w4r(0:(jtot+1),0:(kmax+1))
real    :: p4r(0:(jtot+1),0:(kmax+1))
real    :: y4r(0:(jtot+1),0:(kmax+1))
real    :: h4r(0:(jtot+1),0:(kmax+1)) !#
real    :: cc4r(0:(jtot+1),0:(kmax+1)) !#
real    :: indie4r(0:(jtot+1),0:(kmax+1)) !#

integer, parameter :: k_shift = 0  !##

integer :: ksol,nprocs
parameter(nprocs=dims(1)*dims(2),ksol=kmax/nprocs)
real ::   w2(0:i1,0:j1,0:k1)
real ::   unew2(itot,jtot,ksol)
real ::   vnew2(itot,jtot,ksol)
real ::   wnew2(itot,jtot,ksol)
real ::   pnew2(itot,jtot,ksol)

integer :: l,q,rank,sub_rank,sub_coords(ndims)
integer :: npoints
integer :: i,j,k,im,jm,km,it,jt,iloc,jloc
 character*3 number
 character*7 filenumber

!# plot in paraview 
!aux(:,:,:) = ynew(1:imax,1:jmax,1:kmax)
!call write2dplane(aux,1,itot/2,'yzy',istep)

!aux(:,:,:) = 0.5*(vnew(1:imax,1+1:jmax+1,1:kmax) + vnew(1:imax,1:jmax,1:kmax))
!call write2dplane(aux,1,itot/2,'yzv',istep)
 
      i = itot/2
      q=-1
  114 q=q+1
      if ( (1+q*imax) .lt. i ) go to 114
      q=q-1
      iloc = i - q*imax
      sub_coords(1) = q
      sub_coords(2) = 0
      call MPI_CART_RANK(comm_cart,sub_coords,rank,error) !rank of process to which data has to be sent
      do l=1,dims(2)-1
        sub_coords(1) = q
        sub_coords(2) = l
        call MPI_CART_RANK(comm_cart,sub_coords,sub_rank,error)
        if (myid .eq. sub_rank) then
          do k=1,kmax
            do j=1,jmax

              u4(j,k)= 0.5*(unew(iloc,j,k)+unew(-1+iloc,j,k))
              v4(j,k)= vnew(iloc,j,k) 
!              v4(j,k)= v_l(iloc,j,k)  !#
              w4(j,k)= wnew(iloc,j,k) !w_l(iloc,j,k)
!              w4(j,k)= w_l(iloc,j,k)  !#
              p4(j,k)= pnew(iloc,j,k)
              y4(j,k)= lvset(iloc,j,k +k_shift ,1) !##
              h4(j,k)= fixhv(iloc,j,k) !#
              cc4(j,k)= cutcell(iloc,j,k) !#
              indie4(j,k)= jiemian(iloc,j,k)  !#

            enddo
          enddo
          call MPI_SSEND(u4,jmax*kmax,MPI_REAL8,rank,7,comm_cart,error)
          call MPI_SSEND(v4,jmax*kmax,MPI_REAL8,rank,8,comm_cart,error)
          call MPI_SSEND(w4,jmax*kmax,MPI_REAL8,rank,9,comm_cart,error)
          call MPI_SSEND(p4,jmax*kmax,MPI_REAL8,rank,10,comm_cart,error)
          call MPI_SSEND(y4,jmax*kmax,MPI_REAL8,rank,11,comm_cart,error)
          call MPI_SSEND(h4,jmax*kmax,MPI_REAL8,rank,12,comm_cart,error) !#
          call MPI_SSEND(cc4,jmax*kmax,MPI_REAL8,rank,13,comm_cart,error) !#
          call MPI_SSEND(indie4,jmax*kmax,MPI_REAL8,rank,14,comm_cart,error) !#

        endif
        if (myid .eq. rank) then
          call MPI_RECV(u4,jmax*kmax,MPI_REAL8,sub_rank,7,  &
                        comm_cart,status,error)
          call MPI_RECV(v4,jmax*kmax,MPI_REAL8,sub_rank,8,  &
                        comm_cart,status,error)
          call MPI_RECV(w4,jmax*kmax,MPI_REAL8,sub_rank,9,  &
                        comm_cart,status,error)
          call MPI_RECV(p4,jmax*kmax,MPI_REAL8,sub_rank,10, &
                        comm_cart,status,error)
          call MPI_RECV(y4,jmax*kmax,MPI_REAL8,sub_rank,11, &
                        comm_cart,status,error)
          call MPI_RECV(h4,jmax*kmax,MPI_REAL8,sub_rank,12, &
                        comm_cart,status,error) !#
          call MPI_RECV(cc4,jmax*kmax,MPI_REAL8,sub_rank,13, &
                        comm_cart,status,error) !#
          call MPI_RECV(indie4,jmax*kmax,MPI_REAL8,sub_rank,14, &
                        comm_cart,status,error) !#

          do k=1,kmax
            do j=1,jmax

              u4r(j+l*jmax,k)=u4(j,k)
              v4r(j+l*jmax,k)=v4(j,k)
              w4r(j+l*jmax,k)=w4(j,k)
              p4r(j+l*jmax,k)=p4(j,k)
              y4r(j+l*jmax,k)=y4(j,k)
              h4r(j+l*jmax,k)=h4(j,k) !#
              cc4r(j+l*jmax,k)=cc4(j,k) !#
              indie4r(j+l*jmax,k)=indie4(j,k) !#

            enddo
          enddo
        endif
      enddo
      if (myid .eq. rank) then
        do k=1,kmax
          do j=1,jmax
            u4r(j,k)= 0.5*(unew(iloc,j,k)+unew(-1+iloc,j,k))
            v4r(j,k)= vnew(iloc,j,k)
!            v4r(j,k)= v_l(iloc,j,k)  !#
            w4r(j,k)= wnew(iloc,j,k) !w_l(iloc,j,k)
!            w4r(j,k)= w_l(iloc,j,k)  !#
            p4r(j,k)= pnew(iloc,j,k)
            y4r(j,k)= lvset(iloc,j,k +k_shift,1) !##
            h4r(j,k)=fixhv(iloc,j,k) !#
            cc4r(j,k)= cutcell(iloc,j,k) !#
            indie4r(j,k)= jiemian(iloc,j,k)  !#

          enddo
          !periodic b.c.
          u4r(0,k)     =u4r(jtot,k)
        enddo
        !b.c.
        do j=1,jtot
          w4r(j,0)=0.
        enddo

!!#
!        do k=1,kmax
!          do j=1,jmax
!            h4r(j,k)=0.
!            do i=1,imax

!              h4r(j,k)=h4r(j,k)+hnew(i,j,k)
!            enddo
!            h4r(j,k)=h4r(j,k)/imax
!          enddo
!        enddo
!!#

        if(mod(istep,1).eq.0) then 

          write(filenumber,'(i7.7)') istep
!#          write(filenumber,'(i7.7)') iterations  !# comment if not validating level set re-init
          open(43,file=datadir//'yz_mid'//filenumber//'.dat')
          write(43,*) 'VARIABLES = "y","z","u","v","w","p","s","indie","fixhv"'  !#
          write(43,*) 'ZONE T="Zone1"',' I=',jtot+2*npoints,' J=',kmax,', F=POINT'
!'
          do k=1,kmax
            km=k-1
            do jt=-npoints+1,jtot+npoints
              j=jt
              if (j .lt. 1)    j = j+jtot
              if (j .gt. jtot) j = j-jtot
              jm=jt-1
              if (jm .lt. 1)    jm = jm+jtot
              if (jm .gt. jtot) jm = jm-jtot
              write(43,'(9E15.7)') (jt-0.5)/dyi,(k-0.5)/dzi,u4r(j,k), &
                                   0.5*(v4r(j,k)+v4r(jm,k)),                    &
                                   0.5*(w4r(j,k)+w4r(j,km)),p4r(j,k),&
                                   y4r(j,k), indie4r(j,k), h4r(j,k)  !cc4r(j,k)
            enddo
          enddo
          close(43)


        endif

      endif

!!!!!!!XY
!      do k=1,k1
!        do j=0,j1
!          do i=0,i1
!            w2(i,j,k) = 0.5*(wnew(i,j,k)+wnew(i,j,k-1))
!          enddo
!        enddo
!      enddo
!!     B.c.'s: zero velocity at walls
!      do j=0,j1
!        do i=0,i1
!          w2(i,j,0)    = 0.
!          w2(i,j,kmax) = 0.
!        enddo
!      enddo
!      call zredistribute(unew,unew2,0)
!      call zredistribute(vnew,vnew2,0)
!      call zredistribute(w2  ,wnew2,0)
!      call zredistribute(pnew,pnew2,0)
!      k = kmax/2 
!      l = -1
!  199 l = l+1
!      if ( (1+l*ksol) .lt. k ) go to 199
!      l = l-1
!      rank = l
!      if (myid .eq. rank) then
!        write(number,'(i3.3)') k
!        k    = k - rank*ksol

!        if(mod(istep,100).eq.0) then 
!        open(42,file=datadir//'xy'//number//'time'//filenumber//'.dat')
!        write(filenumber,'(i7.7)') istep
!        write(42,*) 'VARIABLES = "x","y","u","v","w","p"'
!        write(42,*) 'ZONE T="Zone1"',' I=',itot+2*npoints,' J=',jtot+2*npoints,', F=POINT'
!        do jt=-npoints+1,jtot+npoints
!          j=jt
!          if (j .lt. 1)    j = j+jtot
!          if (j .gt. jtot) j = j-jtot
!          jm=jt-1
!          if (jm .lt. 1)    jm = jm+jtot
!          if (jm .gt. jtot) jm = jm-jtot
!          do it=-npoints+1,itot+npoints
!            i=it
!            if (i .lt. 1)    i = i+itot
!            if (i .gt. itot) i = i-itot
!            im=it-1
!            if (im .lt. 1)    im = im+itot
!            if (im .gt. itot) im = im-itot
!            write(42,'(6E15.7)') (it-0.5)/dxi/L_ref,(jt-0.5)/dyi/L_ref,0.5*(unew2(i,j,k)+unew2(im,j,k))/U_ref, &
!                                0.5*(vnew2(i,j,k)+vnew2(i,jm,k))/U_ref,                                      &
!                                wnew2(i,j,k)/U_ref,pnew2(i,j,k)/(U_ref**2.)                                   
!          enddo
!        enddo
!        close(42)
!        endif
!      endif



return
end subroutine post2d





subroutine write2dplane(var,norm,islice,name,istep)
implicit none
real, intent(in), dimension(imax,jmax,kmax) :: var
integer, intent(in) :: norm,islice,istep
character(len=3) :: name
character(len=4) :: slicechar
character(len=7) :: fldnum

write(fldnum,'(i7.7)') istep
write(slicechar,'(i4.4)') islice
select case(norm)
case(1) !normal to x --> yz plane
  call decomp_2d_write_plane(3,var,norm,islice,datadir//name//'_xeq_'//slicechar//'_fld_'//fldnum//'.out')
case(2) !normal to y --> zx plane
  call decomp_2d_write_plane(3,var,norm,islice,datadir//name//'_yeq_'//slicechar//'_fld_'//fldnum//'.out')
case(3) !normal to z --> xy plane
  call decomp_2d_write_plane(3,var,norm,islice,datadir//name//'_zeq_'//slicechar//'_fld_'//fldnum//'.out')
end select

return
end subroutine write2dplane





subroutine post3d(istep)
implicit none
integer, intent(in) :: istep
real, allocatable, dimension(:,:,:) :: var1,var2,var3,var4
integer, parameter :: nprocs=dims(1)*dims(2),ksol=kmax/nprocs
integer :: p
!!
!! pressure
!!
!allocate(var1(1:imax,1:jmax,1:kmax))
!var1(:,:,:) = pnew(1:imax,1:jmax,1:kmax)
!call write3dscal(istep,imax,jmax,kmax,var1,'pre')
!!
!! u velocity
!!
!var1(:,:,:) = unew(1:imax,1:jmax,1:kmax)
!call write3dscal(istep,imax,jmax,kmax,var1,'vex')
!!
!! v velocity
!!
!var1(:,:,:) = vnew(1:imax,1:jmax,1:kmax)
!call write3dscal(istep,imax,jmax,kmax,var1,'vey')
!!
!! w velocity
!!
!var1(:,:,:) = wnew(1:imax,1:jmax,1:kmax)
!call write3dscal(istep,imax,jmax,kmax,var1,'vez')
!!
!! dissipation (rate-of-strain)
!!
!call strain_rate(unew,vnew,wnew,var1)
!call write3dscal(istep,imax,jmax,kmax,var1,'dis')
!!
!! enstrophy
!!
!allocate(var2(1:imax,1:jmax,1:kmax))
!call enstrophy(unew,vnew,wnew,var2)
!!
!! Q-criterion
!!
!allocate(var3(1:imax,1:jmax,1:kmax))
!call q_criterion(var2,var1,var3)
!call write3dscal(istep,imax,jmax,kmax,var3,'qcr')
!!
!! R (third invariant of the velocity-gradient tensor
!!
!allocate(var4(1:imax,1:jmax,1:kmax))
!call compute_r(unew,vnew,wnew,var4)
!!
!! swirling strength (lambda_{ci})
!!
!call swirl(var3,var4,var1)
!call write3dscal(istep,imax,jmax,kmax,var3,'swr')
!!
return
end subroutine post3d



subroutine write3dscal(istep,n1,n2,n3,var,name)
implicit none
integer, intent(in) :: istep,n1,n2,n3
real, intent(in), dimension(n1,n2,n3) :: var
character(len=3), intent(in) :: name
integer :: fh
integer(kind=MPI_OFFSET_KIND) :: filesize,disp
character :: istepchar*7
!
write(istepchar,'(i7.7)') istep
call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//name//istepchar, &
     MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
filesize = 0_MPI_OFFSET_KIND
call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
disp = 0_MPI_OFFSET_KIND
call decomp_2d_write_var(fh,disp,3,var)
call MPI_FILE_CLOSE(fh,error)
!
return
end subroutine write3dscal









subroutine writeallpart(istep,irank,nprank)
implicit none
integer,intent(in) :: istep,irank,nprank
character(len=1), parameter :: lf = char(10)
character(len=400) :: buffer
character(len=50) :: filename,filename2
integer :: ivtk
integer :: p
integer :: e_io,indent
real :: dummy
integer(kind=4) :: idummy
integer(kind=1) :: idummy1
integer :: offst
integer :: counter
offst = 0
!

!ivtk = 99+irank
!write(filename,fmt='(A,(i6.6),A,(i7.7),A)') datadir//'allpart_rk',irank,'_fld',istep,'.vtu'
!
!open(unit = ivtk, file = trim(filename), form = 'binary',access='sequential',  &
!     action = 'write', convert ='big_endian', recordtype='stream', buffered = 'yes', &
!     iostat = e_io)
!write(unit=ivtk,iostat=e_io)'<?xml version="1.0"?>'//lf
!write(unit=ivtk,iostat=e_io)'<VTKFile type="'//'UnstructuredGrid'//'" version="0.1" byte_order="BigEndian">'//lf
!indent = 2
!write(buffer,fmt='(A)',iostat=e_io)repeat(' ',indent)//'<'//'UnstructuredGrid'//'>'
!write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!  indent = indent + 2
!  write(buffer,fmt='(A,(I6),A,(I6),A)',iostat=e_io)repeat(' ',indent)//'<'//'Piece'//' NumberOfPoints="',nprank,'"'// &
!                                                                              ' NumberOfCells="',nprank,'">'
!  write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!  indent = indent + 2
!      write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'<PointData Vectors="vectors">'//lf
!      indent = indent + 2
!      write(buffer,fmt='(A,(I12),A)')repeat(' ',indent)//'<DataArray Name="velocity" NumberOfComponents="3" type="Float64" format="appended" offset="',offst,'">'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      offst = offst + nprank*3*sizeof(dummy)+sizeof(idummy) 
!      write(buffer,fmt='(A)')repeat(' ',indent)//'</DataArray'//'>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      write(buffer,fmt='(A,(I12),A)')repeat(' ',indent)//'<DataArray Name="omegaxyz" NumberOfComponents="3" type="Float64" format="appended" offset="',offst,'">'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      offst = offst + nprank*3*sizeof(dummy)+sizeof(idummy) 
!      write(buffer,fmt='(A)')repeat(' ',indent)//'</DataArray'//'>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      write(buffer,fmt='(A,(I12),A)')repeat(' ',indent)//'<DataArray Name="thetaxyz" NumberOfComponents="3" type="Float64" format="appended" offset="',offst,'">'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      offst = offst + nprank*3*sizeof(dummy)+sizeof(idummy) 
!      write(buffer,fmt='(A)')repeat(' ',indent)//'</DataArray'//'>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!    indent = indent - 2
!    write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'</PointData>'//lf
!    write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'<Points>'//lf
!      indent = indent + 2
!      write(buffer,fmt='(A,(I12),A)')repeat(' ',indent)//'<DataArray Name="Point" NumberOfComponents="3" type="Float64" format="appended" offset="',offst,'">'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      offst = offst + nprank*3*sizeof(dummy)+sizeof(idummy) 
!      write(buffer,fmt='(A)')repeat(' ',indent)//'</DataArray'//'>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      indent = indent - 2
!      write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'</Points>'//lf
!  write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'<Cells>'//lf
!  indent = indent + 2
!  write(buffer,fmt='(A,(I8),A)')repeat(' ',indent)//'<DataArray type="Int32" Name="connectivity" format="appended" offset="',offst,'"></DataArray>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!  offst = offst + sizeof(idummy) + nprank*sizeof(idummy)
!  write(buffer,fmt='(A,(I8),A)')repeat(' ',indent)//'<DataArray type="Int32" Name="offsets" format="appended" offset="',offst,'"></DataArray>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!  offst = offst + sizeof(idummy) + nprank*sizeof(idummy)
!  write(buffer,fmt='(A,(I8),A)')repeat(' ',indent)//'<DataArray type="Int8" Name="types" format="appended" offset="',offst,'"></DataArray>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!  offst = offst + sizeof(idummy) + nprank*sizeof(idummy1)
!  indent = indent - 2

!  write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'</Cells>'//lf
!    indent = indent - 2
!    write(buffer,fmt='(A)',iostat=e_io)repeat(' ',indent)//'</'//'Piece'//'>'
!    write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!    indent = indent - 2
!  write(buffer,fmt='(A)',iostat=e_io)repeat(' ',indent)//'</'//'UnstructuredGrid'//'>'
!  write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!  write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'<AppendedData encoding="raw">'//lf
!  write(unit=ivtk,iostat=e_io)'_'
!  write(unit=ivtk,iostat=e_io) 3*nprank*8!sizeof(dummy) !position and 3 vectors
!  do p=1,pmax
!    if(ap(p)%mslv.gt.0) then
!      write(ivtk,iostat=e_io) ap(p)%u,ap(p)%v,ap(p)%w
!    endif
!  enddo
!  write(unit=ivtk,iostat=e_io) 3*nprank*8!*sizeof(dummy) !position and 3 vectors
!  do p=1,pmax
!    if(ap(p)%mslv.gt.0) then
!      write(ivtk,iostat=e_io) ap(p)%omx,ap(p)%omy,ap(p)%omz
!    endif
!  enddo
!  write(unit=ivtk,iostat=e_io) 3*nprank*8!*sizeof(dummy) !position and 3 vectors
!  do p=1,pmax
!    if(ap(p)%mslv.gt.0) then
!      write(ivtk,iostat=e_io) ap(p)%omx,ap(p)%omy,ap(p)%omz
!    endif
!  enddo
!  write(unit=ivtk,iostat=e_io) 3*nprank*8!sizeof(dummy) !position
!  do p=1,pmax
!    if(ap(p)%mslv.gt.0) then
!      write(ivtk,iostat=e_io) ap(p)%x,ap(p)%y,ap(p)%z
!    endif
!  enddo
!  write(unit=ivtk,iostat=e_io) nprank*4!sizeof(idummy)
!  write(unit=ivtk,iostat=e_io) (p-1,p=1,nprank)
!  write(unit=ivtk,iostat=e_io) nprank*4!sizeof(idummy)
!  write(unit=ivtk,iostat=e_io) (p,p=1,nprank)
!  write(unit=ivtk,iostat=e_io) nprank*1!sizeof(idummy1)
!  idummy1 = 1
!  write(unit=ivtk,iostat=e_io) (idummy1,p=1,nprank)
!  write(unit=ivtk,iostat=e_io) lf
!  write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'</AppendedData>'//lf
!  indent = indent - 2
!write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'</VTKFile>'//lf
!close(ivtk)
!
!if(myid.eq.0) then
!  offst = 0
!  write(filename,fmt='(A,(i6.6),A,(i7.7),A)') datadir//'allpart_fld',istep,'.pvtu'
!  open(unit = ivtk, file = trim(filename), form = 'binary',access='sequential',  &
!       action = 'write', convert ='big_endian', recordtype='stream', buffered = 'yes', &
!       iostat = e_io)
!  write(unit=ivtk,iostat=e_io)'<?xml version="1.0"?>'//lf
!  write(unit=ivtk,iostat=e_io)'<VTKFile type="'//'PUnstructuredGrid'//'" version="0.1" byte_order="BigEndian">'//lf
!  indent = 2
!  write(buffer,fmt='(A)',iostat=e_io)repeat(' ',indent)//'<'//'PUnstructuredGrid GhostLevel="0"'//'>'
!  write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!  indent = indent + 2
!      write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'<PPointData Vectors="vectors">'//lf
!      indent = indent + 2
!      write(buffer,fmt='(A,(I12),A)')repeat(' ',indent)//'<PDataArray Name="velocity" NumberOfComponents="3" type="Float64" format="appended" offset="',offst,'">'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      offst = offst + nprank*3*sizeof(dummy)+sizeof(idummy) 
!      write(buffer,fmt='(A)')repeat(' ',indent)//'</PDataArray'//'>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      write(buffer,fmt='(A,(I12),A)')repeat(' ',indent)//'<PDataArray Name="omegaxyz" NumberOfComponents="3" type="Float64" format="appended" offset="',offst,'">'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      offst = offst + nprank*3*sizeof(dummy)+sizeof(idummy) 
!      write(buffer,fmt='(A)')repeat(' ',indent)//'</PDataArray'//'>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      write(buffer,fmt='(A,(I12),A)')repeat(' ',indent)//'<PDataArray Name="thetaxyz" NumberOfComponents="3" type="Float64" format="appended" offset="',offst,'">'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      offst = offst + nprank*3*sizeof(dummy)+sizeof(idummy) 
!      write(buffer,fmt='(A)')repeat(' ',indent)//'</PDataArray'//'>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!    indent = indent - 2
!    write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'</PPointData>'//lf
!    write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'<PPoints>'//lf
!      indent = indent + 2
!      write(buffer,fmt='(A,(I12),A)')repeat(' ',indent)//'<PDataArray Name="Point" NumberOfComponents="3" type="Float64" format="appended" offset="',offst,'">'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      offst = offst + nprank*3*sizeof(dummy)+sizeof(idummy) 
!      write(buffer,fmt='(A)')repeat(' ',indent)//'</PDataArray'//'>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!      indent = indent - 2
!      write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'</PPoints>'//lf
!  write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'<PCells>'//lf
!  indent = indent + 2
!  write(buffer,fmt='(A,(I8),A)')repeat(' ',indent)//'<PDataArray type="Int32" Name="connectivity" format="appended" offset="',offst,'"></PDataArray>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!  offst = offst + sizeof(idummy) + nprank*sizeof(idummy)
!  write(buffer,fmt='(A,(I8),A)')repeat(' ',indent)//'<PDataArray type="Int32" Name="offsets" format="appended" offset="',offst,'"></PDataArray>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!  offst = offst + sizeof(idummy) + nprank*sizeof(idummy)
!  write(buffer,fmt='(A,(I8),A)')repeat(' ',indent)//'<PDataArray type="Int8" Name="types" format="appended" offset="',offst,'"></PDataArray>'
!      write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!  offst = offst + sizeof(idummy) + nprank*sizeof(idummy1)
!  indent = indent - 2
!  write(unit=ivtk,iostat=e_io)repeat(' ',indent)//'</PCells>'//lf
!  do p=0,product(dims)-1
!    write(filename2,fmt='(A,(i6.6),A,(i7.7),A)') 'allpart_rk',p,'_fld',istep,'.vtu'
!    write(unit=ivtk,iostat=e_io) repeat(' ',indent)//'<Piece Source="'//trim(filename2)//'"/>'//lf
!  enddo
!  indent = indent - 2
!  write(buffer,fmt='(A)',iostat=e_io) repeat(' ',indent)//'<'//'/PUnstructuredGrid'//'>'
!  write(unit=ivtk,iostat=e_io)trim(buffer)//lf
!  indent = indent - 2
!  write(buffer,fmt='(A)',iostat=e_io) repeat(' ',indent)//'</VTKFile>'
!  write(unit=ivtk,iostat=e_io)trim(buffer)
!  close(ivtk)
!endif
!
return
end subroutine writeallpart
!
end module mod_output

