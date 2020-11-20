! THERMODRIFT VELOCITIES following Geppert & Wiebecke 1991
! vr  = -c*Qhall/B*(dT/dr - omega*tau bphi 1/r(dT/d theta) )
! vth = -c*Qhall/B*(1/r(dT/d theta) + omega*tau bphi (dT/dr) )
! vphi= -c*Qhall/B*(omega*tau*br(1/r)(dT/d theta) - omega*tau btheta dT/dr )
! where br, bth, bphi are the normalized components of MF

      SUBROUTINE THERMOELECTRIC(pc)
      implicit none
      include '../decl/dim2.h'
      include '../decl/grid.h'
      include '../decl/varsB.h'
      include '../decl/thermo.h'

      integer i,j,k,l
      real*8 pc(kd,ld)
      real*8 gradtri(kd,ld),gradtthi(kd,ld)
      real*8 omegatau2(0:nang+1,2*lcore:np+2)
      real*8 tem(4,kd,ld)

      do i=1,nang
      do j=2*lcore,np
	if (etab(i,j) .ne. 0.) omegatau2(i,j)=fh(j)*bm(i,j)/etab(i,j)
      enddo
      enddo

      do k=1,kmax
      do l=lcore-2,lmax
	gradtri(k,l) = ( pc(k,l+1) - pc(k,l) )/lr(j)
	gradtthi(k,l)= ( pc(k+1,l) - pc(k,l) )/lth(j)
      enddo
      enddo
      gradtri (kd,:)= gradtri(kmax,:)
      gradtthi(kd,:)=-gradtthi(kmax,:)
      gradtri (:,ld)= gradtri(:,lmax)
      gradtthi(:,ld)= gradtthi(:,lmax)

      gradtr =0.
      gradtth=0.
      vrte   =0.
      vthte  =0.
      vphite =0.

      do k=1,kmax
      do l=lcore-1,lmax-2
	i=2*k-1
	j=2*l-1
	gradtr(i,j) = 0.25d0*(gradtri(k,l)+gradtri(k+1,l)+
     &		gradtri(k,l-1)+gradtri(k+1,l-1))
	gradtth(i,j)= gradtthi(k,l)

	i=2*k-1
	j=2*l
	gradtr(i,j) = 0.5d0*(gradtri(k,l)+gradtri(k+1,l))
	gradtth(i,j)= 0.5d0*(gradtthi(k,l)+gradtthi(k,l+1))

       if (k .ne. 1) then
	i=2*k-2
	j=2*l
	gradtr(i,j) = gradtri(k,l)
	gradtth(i,j)= 0.25d0*(gradtthi(k,l)+gradtthi(k,l+1)+
     &		gradtthi(k-1,l)+gradtthi(k-1,l+1))

	i=2*k-2
	j=2*l-1
	gradtr(i,j) = 0.5d0*(gradtri(k,l)+gradtri(k,l-1))
	gradtth(i,j)= 0.5d0*(gradtthi(k,l)+gradtthi(k-1,l))
       endif

      enddo
      enddo

      do i=1,nang
      do j=2*lcore+1,np-2
	if (bm(i,j) .ne. 0.) then
	 vrte(i,j)  =-qfac*qhall(i,j)/bm(i,j)*
     &	(gradtr(i,j)-omegatau2(i,j)*bphi(i,j)*gradtth(i,j)/bm(i,j))
	 vthte(i,j) =-qfac*qhall(i,j)/bm(i,j)*
     &	(gradtth(i,j)+omegatau2(i,j)*bphi(i,j)*gradtr(i,j)/bm(i,j))
	 vphite(i,j)=-qfac*qhall(i,j)*fh(j)/etab(i,j)*
     &	(br(i,j)*gradtth(i,j)-bth(i,j)*gradtr(i,j))/bm(i,j) 
	endif

! 	if (vrte(i,j) .ne. vrte(i,j) .or. 
!      &		vthte(i,j) .ne. vthte(i,j) .or.
!      &		vphite(i,j) .ne. vphite(i,j)) then
! 	print*,i,j,vrte(i,j),vthte(i,j),vphite(i,j)
! 	endif
      enddo
      enddo

!        thloss=0.
!        do k=1,kmax
!        do l=lcore+1,lmax-1
! 	   i=2*k-2
! 	   j=2*l-1
! 	   thloss(k,l)=qpar(i,j)*6.86e-13* ( 
!      &	(pc(k+1,l)-pc(k-1,l))/(r(l)*(z(k+1)-z(k-1)))*jth(i,j) + 
!      &  (pc(k,l+1)-pc(k,l-1))/(r(l+1)-r(l-1))*jr(i,j)   )
!        enddo
!        enddo

! 	  if (it .eq. 1) then
! 	  open(unit=11,file='out/thermopower.d')
! 	  write(11,*) tyear
! 	  do k=1,kmax
! 	  do l=lc,lmax
! 	    write(11,*) z(k),r(l),qcj(k,l),emnu(k,l),thloss(k,l)
! 	  enddo
! 	  enddo
! 	  close(11)
! 	  endif


       END SUBROUTINE
