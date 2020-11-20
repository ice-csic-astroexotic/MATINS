!-----------------------------------------------------------------------
! Additional subprograms related to the boundary conditions.
! Some of there were previously included in bevol.f and binit.f.
!
! Changes:
! - New subroutine VACUUMBC_CENTER for vacuum conditions in the center.
!
! Contents:
! A) Subroutines:
! GAULEG
! GETBL
! GET_GREEN
! INVERT_BC
! LEGPOL
! VACUUMBC
! VACUUMBC_CENTER
! B) Functions:
! FINT
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
        subroutine GET_GREEN(m)

	implicit none
	integer i,j,k,m,mi,resi
	parameter(resi=10)
	real*8 pi,fint,thm(m)
	real*8 thi(resi*m),dthi

! Common block shared with INVERT_BC.
	real*8 fw(1000,1000),dthm
	common /fw_green/ fw,dthm

	mi=resi*m
	if(mi.gt.1000)then
	   write(*,*)"GET_GREEN: Parameter too high!"
	   stop
	endif
	pi=dacos(-1d0)
	dthm=pi/dble(m)
	dthi=pi/dble(mi)
	do j=1,m
	   thm(j)=dble(j-0.5d0)*dthm
	enddo
	do k=1,mi
	   thi(k)=dble(k-0.5d0)*dthi
	enddo

	fw=0d0
	do i=1,m
	  do j=1,m
        do k=1,mi
		if(thi(k).ge.thm(j)-0.5d0*dthm.and.thi(k).le.thm(j)+0.5d0*dthm)then
		    fw(i,j)=fw(i,j)+fint(thm(i),thi(k),pi)*dthi
		 endif
	    enddo
	  enddo
	enddo

	return
	end
!-----------------------------------------------------------------------
! Green's external vacuum boundary conditions.
!-----------------------------------------------------------------------
	subroutine INVERT_BC(m,brad,btheta)

	use math, only: ludcmp, lubksb

	implicit none
	integer i,j,m
	integer*4 indx(m)
	real*8 pi,d
	real*8 b(m),brad(m),btheta(m+1)
	real*8 matpsi(m,m)

! Common block shared with GET_GREEN.
	real*8 fw(1000,1000),dthm
	common /fw_green/ fw,dthm

	pi=dacos(-1d0)

	do i=1,m
	   do j=1,m
	      matpsi(i,j)=0.5d0*fw(i,j)
	      if(j.eq.i) matpsi(i,j)=matpsi(i,j)+2d0*pi
	   enddo
	enddo

	b=0d0
	do i=1,m
	   do j=1,m
	      b(i)=b(i)-brad(j)*fw(i,j)
	   enddo
	enddo

	call ludcmp(m, matpsi, indx, d)
	call lubksb(m, matpsi, indx, b)

	do i=2,m
	   btheta(i)=(b(i)-b(i-1))/dthm
	enddo

	btheta(1)=0d0
	btheta(m+1)=0d0

	return
	end
	
!
!-----------------------------------------------------------------------
	real*8 function fint(th,thp,pi)

	implicit none
	integer i,nphi
	parameter(nphi=1000)
	real*8 th,thp,sphi(nphi),dphi,pi
	real*8 snt,sntp
	real*8 g1,g2,intk

	snt=dsin(th)
	sntp=dsin(thp)
	g1=1d0-dcos(th-thp)
	g2=2d0*snt*sntp
	dphi=0.5d0*pi/nphi
	do i=1,nphi
	   sphi(i)=dsin((i-0.5d0)*dphi)
	enddo

	intk=0d0
	do i=1,nphi
	   intk=intk+dphi/dsqrt(g1+g2*sphi(i)**2)
	enddo

	fint=dsqrt(8d0)*sntp*intk

	return
	end
