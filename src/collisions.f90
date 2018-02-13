module mod_collisions !VER O QUE ACONTECE SE TIVER NDT=40 E ITERMAX = 1
use mpi
use mod_common
use mod_common_mpi
implicit none
private
public collisions, lubrication
contains
!
subroutine collisions(p,q,qq,idp,idq,dist,deltax,deltay,deltaz, &
                      x_nb,y_nb,z_nb,                           &
                      u_nb,v_nb,w_nb,omx_nb,omy_nb,omz_nb,dtsub,kkk, &
                      crx_nb,cry_nb,crz_nb,radius_nb,touchx_nb, &
                      touchy_nb,touchz_nb )
                    
implicit none
integer, intent(in) :: p,q,qq,idp,idq,kkk
real, intent(in) :: dist,deltax,deltay,deltaz, &
                    x_nb,y_nb,z_nb,u_nb,v_nb,w_nb,omx_nb,omy_nb,omz_nb, &
                    dtsub,crx_nb,cry_nb,crz_nb,radius_nb,touchx_nb, &
                    touchy_nb,touchz_nb

integer :: i
real :: distold,dtabs,radius1,radius2
real :: nx,ny,nz,nxold,nyold,nzold,tx,ty,tz
real :: dxn,dyn,dzn,dxt,dyt,dzt
real :: dun,dvn,dwn,dut,dvt,dwt,vabs,vtabs
real :: fxn,fyn,fzn,fnabs,fxt,fyt,fzt,ftabs
real :: psi,psi_crit
real :: deltan,kn,etan,kt,etat,muc
real :: torqx,torqy,torqz,torqxn,torqyn,torqzn
real :: hx,hy,hz,aa,bb,cc,nx1,ny1,nz1,nx2,ny2,nz2
real :: h11,h12,h13, &
        h21,h22,h23, &
        h31,h32,h33
real :: meffn,mefft
real :: dxt_aux,dyt_aux,dzt_aux

!
! vector (nx,ny,nz) is normal unit vector, pointing from particle p towards particle q
!
 nx = deltax/dist
 ny = deltay/dist
 nz = deltaz/dist


if (idq.le.np) then
  deltan = ap(p)%Relative_Centers(kkk,4) + radius_nb - dist
  meffn = meffn_ss
  mefft = mefft_ss
  kn = kn_ss
  kt = kt_ss
  etan = etan_ss
  etat = etat_ss
  muc = muc_ss
  psi_crit = psi_crit_ss
else
  deltan = ap(p)%Relative_Centers(kkk,4) - dist
  meffn = meffn_sw
  mefft = mefft_sw
  kn = kn_sw
  kt = kt_sw
  etan = etan_sw
  etat = etat_sw
  muc = muc_sw
  psi_crit = psi_crit_sw
endif
!
! relative velocity between particle p and q (Nb: p - q !)
!
radius1 = (ap(p)%x - ap(p)%Relative_Centers(kkk,5))**2.
radius1 = radius1 + (ap(p)%y - ap(p)%Relative_Centers(kkk,6))**2.
radius1 = radius1 + (ap(p)%z - ap(p)%Relative_Centers(kkk,7))**2.
radius1 = sqrt(radius1)

nx1 = -(ap(p)%x - ap(p)%Relative_Centers(kkk,5)) / radius1
ny1 = -(ap(p)%y - ap(p)%Relative_Centers(kkk,6)) / radius1
nz1 = -(ap(p)%z - ap(p)%Relative_Centers(kkk,7)) / radius1


if (idq.le.np) then

  radius2 = (x_nb - touchx_nb)**2.
  radius2 = radius2 + (y_nb - touchy_nb)**2.
  radius2 = radius2 + (z_nb - touchz_nb)**2.
  radius2 = sqrt(radius2)

  nx2 = (x_nb - touchx_nb) / radius2
  ny2 = (y_nb - touchy_nb) / radius2
  nz2 = (z_nb - touchz_nb) / radius2

else
 
  radius2 = 0.

  nx2 = 0.
  ny2 = 0.
  nz2 = 0.

endif

dun = (ap(p)%u+(ap(p)%omy*radius1*nz1) - (ap(p)%omz*radius1*ny1)) - &
       (u_nb + (omy_nb*radius2*(-nz2)) - (omz_nb*radius2*(-ny2)))
dvn = (ap(p)%v-(ap(p)%omx*radius1*nz1) + (ap(p)%omz*radius1*nx1)) - &
       (v_nb - (omx_nb*radius2*(-nz2)) + (omz_nb*radius2*(-nx2)))
dwn = (ap(p)%w+(ap(p)%omx*radius1*ny1) - (ap(p)%omy*radius1*nx1)) - &
       (w_nb + (omx_nb*radius2*(-ny2)) - (omy_nb*radius2*(-nx2)))

vabs = dun*nx + dvn*ny + dwn*nz
 
!
! computation of contact forces based on soft-sphere model
!
fxn = - kn*(deltan*nx) - etan*(vabs*nx)
fyn = - kn*(deltan*ny) - etan*(vabs*ny)
fzn = - kn*(deltan*nz) - etan*(vabs*nz)

!write(6,*) 'deltan',deltan*nx,deltan*ny,deltan*nz

if (deltan .lt. 0.) then
  fxn = 0.
  fyn = 0.
  fzn = 0.
endif

torqxn = radius1*(ny1*fzn-nz1*fyn)
torqyn = radius1*(nz1*fxt-nx1*fzt)
torqzn = radius1*(nx1*fyt-ny1*fxt)


!
! relative tangential velocity
!
dut = dun - vabs*nx
dvt = dvn - vabs*ny
dwt = dwn - vabs*nz
!
vtabs = sqrt(dut**2.+dvt**2.+dwt**2.)
!
! prediction of contact instant skipped for now for simplicity's sake
!
distold = sqrt(op(p)%dx(qq)**2.+op(p)%dy(qq)**2.+op(p)%dz(qq)**2.)
if(distold.eq.0) then
  nxold = nx
  nyold = ny
  nzold = nz
else
  nxold = op(p)%dx(qq)/distold
  nyold = op(p)%dy(qq)/distold
  nzold = op(p)%dz(qq)/distold
endif
!
! computation of rotation matrix hij
!
if(ap(p)%firstc(qq).eq.idq) then
  psi = ap(p)%psi(qq)
  hx = ny*nzold-nyold*nz
  hy = nz*nxold-nzold*nx
  hz = nx*nyold-nxold*ny
  aa = sqrt(hx**2. + hy**2. + hz **2.)
  cc = cos(asin(aa))
  bb = 1. - cc
  if(aa.eq.0.) then
    hx = 0.
    hy = 0.
    hz = 0. 
  else
    hx = hx/aa
    hy = hy/aa
    hz = hz/aa
  endif
  h11 = bb*hx**2. + cc
  h12 = bb*hx*hy - aa*hz
  h13 = bb*hx*hz + aa*hy
  h21 = bb*hx*hy + aa*hz
  h22 = bb*hy**2. + cc
  h23 = bb*hy*hz - aa*hx
  h31 = bb*hx*hz - aa*hy
  h32 = bb*hy*hz + aa*hx
  h33 = bb*hz**2. + cc
!  
! tangential displacement evolved in time
!
  dxt = op(p)%dxt(qq)*h11 + op(p)%dyt(qq)*h21 + op(p)%dzt(qq)*h31 + &
        0.5*dtsub*rkcoeffab(rkiter)*(dut + &
        op(p)%dut(qq)*h11 + op(p)%dvt(qq)*h21 + op(p)%dwt(qq)*h31)
  dyt = op(p)%dxt(qq)*h12 + op(p)%dyt(qq)*h22 + op(p)%dzt(qq)*h32 + &
        0.5*dtsub*rkcoeffab(rkiter)*(dvt + &
        op(p)%dut(qq)*h12 + op(p)%dvt(qq)*h22 + op(p)%dwt(qq)*h32)
  dzt = op(p)%dxt(qq)*h13 + op(p)%dyt(qq)*h23 + op(p)%dzt(qq)*h33 + &
        0.5*dtsub*rkcoeffab(rkiter)*(dwt + &
        op(p)%dut(qq)*h13 + op(p)%dvt(qq)*h23 + op(p)%dwt(qq)*h33)
else
!  ap(p)%firstc(qq) = idq
  psi = vtabs/abs(vabs)
  ap(p)%psi(qq) = psi
  dxt = 0.
  dyt = 0.
  dzt = 0.
endif
  dtabs = sqrt(dxt**2.+dyt**2.+dzt**2.)
!
! computation of tangential force
!
!if(psi.lt.psi_crit) then ! stick
!  fxt = - (kt*dxt) - (etat*dut)
!  fyt = - (kt*dyt) - (etat*dvt)
!  fzt = - (kt*dzt) - (etat*dwt)
!  print*,'stick'
!else ! slip
!  fnabs = sqrt(fxn**2.+fyn**2.+fzn**2.)
!  fxt = - muc*fnabs*tx
!  fyt = - muc*fnabs*ty
!  fzt = - muc*fnabs*tz
!  print*,'slip'
!endif
!
fxt = - (kt*dxt) - (etat*dut)
fyt = - (kt*dyt) - (etat*dvt)
fzt = - (kt*dzt) - (etat*dwt)
fnabs = sqrt(fxn**2.+fyn**2.+fzn**2.)
ftabs = sqrt(fxt**2.+fyt**2.+fzt**2.)
tx = - fxt/ftabs
ty = - fyt/ftabs
tz = - fzt/ftabs
if(ftabs.le.muc*fnabs) then
!  write(*,*) 'Sticking'
else
dxt_aux = muc*fnabs*tx/kt !!!!!
dyt_aux = muc*fnabs*ty/kt !!!!!
dzt_aux = muc*fnabs*tz/kt !!!!!
if(sqrt(dxt_aux**2.+dyt_aux**2.+dzt_aux**2.).lt.dtabs) then !!!!!!
  dxt = dxt_aux  !!!!!
  dyt = dyt_aux  !!!!!
  dzt = dzt_aux  !!!!!
endif !!!!!
fxt = - (kt*dxt) - (etat*dut)
fyt = - (kt*dyt) - (etat*dvt)
fzt = - (kt*dzt) - (etat*dwt)
ftabs = sqrt(fxt**2.+fyt**2.+fzt**2.)
  if(ftabs.le.muc*fnabs) then
    !we were lucky
!    write(*,*) 'frontier between stick and slip'
  else
    fxt = -muc*fnabs*tx
    fyt = -muc*fnabs*ty
    fzt = -muc*fnabs*tz
!    write(*,*) 'Slipping'
  endif
endif
!
! computation of collision torque
!
torqx = radius1*(ny1*fzt-nz1*fyt)
torqy = radius1*(nz1*fxt-nx1*fzt)
torqz = radius1*(nx1*fyt-ny1*fxt)
!
! add contribution from this contact to the total contact forces/torques
!
ap(p)%colfx = ap(p)%colfx + fxn + fxt
ap(p)%colfy = ap(p)%colfy + fyn + fyt
ap(p)%colfz = ap(p)%colfz + fzn + fzt
!
ap(p)%coltx = ap(p)%coltx + torqxn + torqx 
ap(p)%colty = ap(p)%colty + torqyn + torqy  
ap(p)%coltz = ap(p)%coltz + torqzn + torqz  

!write(6,*) 'fxn...,fyn...,fzn...',fxn,fyn,fzn,idp
!write(6,*) 'fxt...,fyt...,fzt...',fxt,fyt,fzt,idp
!write(6,*) 'torqxn,torqyn,torqzn',torqxn,torqyn,torqzn,idp
!write(6,*) 'torqx.,torqy.,torqz.',torqx,torqy,torqz,idp
!
! update variables needed for integrating the tangential displacement
!
ap(p)%dx(qq) = deltax
ap(p)%dy(qq) = deltay
ap(p)%dz(qq) = deltaz
ap(p)%dxt(qq) = dxt
ap(p)%dyt(qq) = dyt
ap(p)%dzt(qq) = dzt
ap(p)%dut(qq) = dut
ap(p)%dvt(qq) = dvt
ap(p)%dwt(qq) = dwt
!
return
end subroutine collisions
!
subroutine lubrication(p,q,idp,idq,nx,ny,nz,eps,kappa, &
                       x_nb,y_nb,z_nb,u_nb,v_nb,w_nb,  &
                       omx_nb,omy_nb,omz_nb,kkk,       &  
                       crx_nb,cry_nb,crz_nb,radius_nb, &
                       touchx_nb,touchy_nb,touchz_nb )

implicit none
integer, intent(in) :: p,q,idp,idq,kkk
real, intent(in) :: nx,ny,nz,eps,kappa,   &
                    x_nb,y_nb,z_nb,       &                   
                    u_nb,v_nb,w_nb,       &
                    omx_nb,omy_nb,omz_nb, &
                    crx_nb,cry_nb,crz_nb, &
                    radius_nb,touchx_nb,  &
                    touchy_nb,touchz_nb
                    
!
! 1-> squeezing direction; 2 and 3 -> shearing directions
! see paper by dance and Maxey
!
real :: a11,a11a,a11b          
real :: u1a,u2a,u3a,omg1a,omg2a,omg3a 
real :: u1b,u2b,u3b,omg1b,omg2b,omg3b
real :: t1x,t1y,t1z,t2x,t2y,t2z
real :: f1,f2,f3,t1,t2,t3
real :: dun,dvn,dwn,dut,dvt,dwt,vabs,vtabs
real :: radius1,nx1,ny1,nz1
real :: radius2,nx2,ny2,nz2
real :: Vel_Ctr_p(1:3),Vel_Ctr_q(1:3),coeff_f
!
! determine the unit vectors t1x and t2x from nx,ny,nz
! and the velocity differences
!

radius1 = (ap(p)%x - ap(p)%Relative_Centers(kkk,5))**2.
radius1 = radius1 + (ap(p)%y - ap(p)%Relative_Centers(kkk,6))**2.
radius1 = radius1 + (ap(p)%z - ap(p)%Relative_Centers(kkk,7))**2.
radius1 = sqrt(radius1)

nx1 = -(ap(p)%x - ap(p)%Relative_Centers(kkk,5)) / radius1
ny1 = -(ap(p)%y - ap(p)%Relative_Centers(kkk,6)) / radius1
nz1 = -(ap(p)%z - ap(p)%Relative_Centers(kkk,7)) / radius1



if (idq.le.np) then

  radius2 = (x_nb - touchx_nb)**2.
  radius2 = radius2 + (y_nb - touchy_nb)**2.
  radius2 = radius2 + (z_nb - touchz_nb)**2.
  radius2 = sqrt(radius2)

  nx2 = (x_nb - touchx_nb) / radius2
  ny2 = (y_nb - touchy_nb) / radius2
  nz2 = (z_nb - touchz_nb) / radius2

else
 
  radius2 = 0.

  nx2 = 0.
  ny2 = 0.
  nz2 = 0.

endif


dun = (ap(p)%u+(ap(p)%omy*radius1*nz1) - (ap(p)%omz*radius1*ny1)) - &
       (u_nb + (omy_nb*radius2*(-nz2)) - (omz_nb*radius2*(-ny2)))
dvn = (ap(p)%v-(ap(p)%omx*radius1*nz1) + (ap(p)%omz*radius1*nx1)) - &
       (v_nb - (omx_nb*radius2*(-nz2)) + (omz_nb*radius2*(-nx2)))
dwn = (ap(p)%w+(ap(p)%omx*radius1*ny1) - (ap(p)%omy*radius1*nx1)) - &
       (w_nb + (omx_nb*radius2*(-ny2)) - (omy_nb*radius2*(-nx2)))

vabs = dun*nx + dvn*ny + dwn*nz


dut = dun - vabs*nx
dvt = dvn - vabs*ny
dwt = dwn - vabs*nz
!
vtabs = sqrt(dut**2.+dvt**2.+dwt**2.)
!
! t1i has the direction of the tangential velocity at contact point
!
t1x = dut/vtabs
t1y = dvt/vtabs
t1z = dwt/vtabs
!
! t2i results from the cross product of ni and t1i
!
t2x = ny*t1z-nz*t1y
t2y = nz*t1x-nx*t1z
t2z = nx*t1y-ny*t1x



Vel_Ctr_p(1) = ap(p)%u + ap(p)%omy*(ap(p)%Relative_Centers(kkk,3)-ap(p)%z) &
                       - ap(p)%omz*(ap(p)%Relative_Centers(kkk,2)-ap(p)%y)
Vel_Ctr_p(2) = ap(p)%v + ap(p)%omz*(ap(p)%Relative_Centers(kkk,1)-ap(p)%x) &
                       - ap(p)%omx*(ap(p)%Relative_Centers(kkk,3)-ap(p)%z)
Vel_Ctr_p(3) = ap(p)%w + ap(p)%omx*(ap(p)%Relative_Centers(kkk,2)-ap(p)%y) &
                       - ap(p)%omy*(ap(p)%Relative_Centers(kkk,1)-ap(p)%x)

Vel_Ctr_q(1) = u_nb + omy_nb*(crz_nb-z_nb) - omz_nb*(cry_nb-y_nb)                    
Vel_Ctr_q(2) = v_nb + omz_nb*(crx_nb-x_nb) - omx_nb*(crz_nb-z_nb)                    
Vel_Ctr_q(3) = w_nb + omx_nb*(cry_nb-y_nb) - omy_nb*(crx_nb-x_nb)
                    

!
! compute velocity differences in local reference frame
!
u1a = Vel_Ctr_p(1)*nx  + Vel_Ctr_p(2)*ny  + Vel_Ctr_p(3)*nz
u1b = Vel_Ctr_q(1)*nx  + Vel_Ctr_q(2)*ny  + Vel_Ctr_q(3)*nz
u2a = Vel_Ctr_p(1)*t1x + Vel_Ctr_p(2)*t1y + Vel_Ctr_p(3)*t1z
u2b = Vel_Ctr_q(1)*t1x + Vel_Ctr_q(2)*t1y + Vel_Ctr_q(3)*t1z
u3a = Vel_Ctr_p(1)*t2x + Vel_Ctr_p(2)*t2y + Vel_Ctr_p(3)*t2z
u3b = Vel_Ctr_q(1)*t2x + Vel_Ctr_q(2)*t2y + Vel_Ctr_q(3)*t2z
omg1a = ap(p)%omx*nx + ap(p)%omy*ny + ap(p)%omz*nz
omg1b = omx_nb*nx + omy_nb*ny + omz_nb*nz
omg2a = ap(p)%omx*t1x + ap(p)%omy*t1y + ap(p)%omz*t1z
omg2b = omx_nb*t1x + omy_nb*t1y + omz_nb*t1z
omg3a = ap(p)%omx*t2x + ap(p)%omy*t2y + ap(p)%omz*t2z
omg3b = omx_nb*t2x + omy_nb*t2y + omz_nb*t2z
!
! compute lubrication forces
!
if(idq.gt.np) then ! particle-wall collision

  if(eps.lt.eps_sat_pw) then

    a11 = a11_sat_pw

  else

    a11 = -1./eps+1./5.*log(eps)+1./21.*eps*log(eps)-.9713

  endif

  f1 = (a11*u1a) ! squeezing

else

  if(eps.lt.eps_sat_pp) then

    a11  = -1./( ((1.-kappa)**2.)*0.001 ) + &            
           ((1. - (7.*kappa) + (kappa**2.))/(5.*((1.-kappa)**3.)))*log(0.001) + &
           ((1.-18.*kappa-29.*(kappa**2.)-18.*(kappa**3.)+(kappa**4.))/(21.*((1.-kappa)**4.)))*0.001*log(0.001) &
           - 0.712484968798678 !- 0.130097750043269*0.001
    a11a = a11

    a11b = -a11

  else

    a11  = -1./(((1.-kappa)**2.)*eps) + &
           ((1.-(7.*kappa)+(kappa**2.))/(5.*((1.-kappa)**3.)))*log(eps) + &
           ((1.-18.*kappa-29.*(kappa**2.)-18.*(kappa**3.)+(kappa**4.))/(21.*((1.-kappa)**4.)))*eps*log(eps) &
           - 0.712484968798678 !- 0.130097750043269*eps
    a11a = a11

    a11b = -a11

  endif
  f1 = (a11a*u1a+a11b*u1b) ! squeezing 
endif
!

coeff_f = 6.*pi*visc*ap(p)%Relative_Centers(kkk,4)

if(idq.gt.np) then ! particle-wall collision

  a11 = a11_ini_pw
  f1 = coeff_f*(f1-(a11*u1a))

else

  a11  = -1./( ((1.-kappa)**2.)*0.025 ) + &            
              ((1. - (7.*kappa) + (kappa**2.))/(5.*((1.-kappa)**3.)))*log(0.025) + &
              ((1.-18.*kappa-29.*(kappa**2.)-18.*(kappa**3.)+(kappa**4.))/(21.*((1.-kappa)**4.)))*0.025*log(0.025) &
             - 0.712484968798678 !- 0.130097750043269*eps2
              !- 0.638779693707812 - (6.282706246033640E-002)*0.025
  a11a = a11
  a11b = -a11
!
endif

f1 = coeff_f*(f1-(a11a*u1a+a11b*u1b))

!
! project forces in global reference frame
! non-normal interactions neglected for now
!
ap(p)%colfx = ap(p)%colfx + f1*nx! + f2*t1x + f3*t2x
ap(p)%colfy = ap(p)%colfy + f1*ny! + f2*t1y + f3*t2y
ap(p)%colfz = ap(p)%colfz + f1*nz! + f2*t1z + f3*t2z
!ap(p)%coltx = ap(p)%colfx + t1*nx + t2*t1x + t3*t2x
!ap(p)%colty = ap(p)%colfx + t1*ny + t2*t1y + t3*t2y
!ap(p)%coltz = ap(p)%colfz + t1*nz + t2*t1z + t3*t2z
!
!write(6,*) 'ap(p)%colfz',ap(p)%colfz

return
end subroutine lubrication
!
end module mod_collisions
